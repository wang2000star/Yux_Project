#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>
#include <omp.h>
#include <immintrin.h> // 包含 SIMD 指令集
#include <array>
#include <vector>
#include <algorithm> // 用于 std::copy
#include <iomanip>   // 包含此头文件以使用 std::setw 和 std::setfill#include <iomanip> // 包含此头文件以使用 std::setw 和 std::setfill
using namespace std;

const double clocks2seconds = 1. / CLOCKS_PER_SEC;

// key-switching parameters
using iksP = TFHEpp::lvl10param;
// bootstrapping parameters
using bkP = TFHEpp::lvl02param;
// private key-switching parameters
using privksP = TFHEpp::lvl21param;
using TLWE_0 = TFHEpp::TLWE<typename bkP::domainP>;
using TLWE_1 = TFHEpp::TLWE<typename privksP::targetP>;   // level 1
using TRLWE_1 = TFHEpp::TRLWE<typename privksP::targetP>; // level 1
const uint32_t byte_mul2[8] = {0, 0, 0, 0, 0, 0, 0, 0};

/*
16进制转2进制，8bits，
binstr[0]是最低位，binstr[7]是最高位，
有一个性质：00~7F的最高位是0，80~FF的最高位是1
所以sboxtable前128个元素的最高位是0，后128个元素的最高位是1
*/
void HexToBinStr(int hex, int *bin_str, int n)
{
    for (int i = 0; i < n; ++i)
    {
        bin_str[i] = (hex >> i) & 1; // 位移操作替代 % 2
    }
}
// 2进制转16进制
void BinStrToHex(const int *bin_str, int &dec_hex, int n)
{
    for (int i = 0; i < n; ++i)
    {
        dec_hex |= (bin_str[i] << i); // 使用位操作重构
    }
}
// 两个二维数组a,b相加，结果保存在数组result中，数组大小为8*(n+1)
template <class P>
void XOR_Two(P &result, P &a, P &b, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int num = 0; num < bkP::domainP::n + 1; num++)
        {
            result[i][num] = a[i][num] + b[i][num];
        }
    }
}
template <class P>
void TRlWE1_XOR_Two(P &result, P &a, P &b)
{
        for (int k=0; k< 2; k++){
        for (int num = 0; num < 1024; num++)
        {
            result[k][num] = a[k][num] + b[k][num];
        }
        }
}
void TRLWE1_XOR_Four(TRLWE_1 &result, TRLWE_1 &a, TRLWE_1 &b, TRLWE_1 &c, TRLWE_1 &d)
{
    TRlWE1_XOR_Two(result, a, b);
    TRlWE1_XOR_Two(result, result, c);
    TRlWE1_XOR_Two(result, result, d);
}
// 四个一维数组a,b,c,d相加，结果保存在数组result中，数组大小为8*(n+1)
template <class P>
void XOR_Four(P &result, P &a, P &b, P &c, P &d)
{
    XOR_Two<P>(result, a, b, 8);
    XOR_Two<P>(result, result, c, 8);
    XOR_Two<P>(result, result, d, 8);
}

// 定义4个s盒，分别用于固定低四位为0或高四位为0的8位多项式相乘
static const unsigned char AmBm[16][16] = {
    // 0    1     2     3     4     5     6     7     8     9     A     B     C     D     E     F
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // 0
    0x00, 0x1b, 0x36, 0x2d, 0x6c, 0x77, 0x5a, 0x41, 0xd8, 0xc3, 0xee, 0xf5, 0xb4, 0xaf, 0x82, 0x99, // 1
    0x00, 0x36, 0x6c, 0x5a, 0xd8, 0xee, 0xb4, 0x82, 0xab, 0x9d, 0xc7, 0xf1, 0x73, 0x45, 0x1f, 0x29, // 2
    0x00, 0x2d, 0x5a, 0x77, 0xb4, 0x99, 0xee, 0xc3, 0x73, 0x5e, 0x29, 0x04, 0xc7, 0xea, 0x9d, 0xb0, // 3
    0x00, 0x6c, 0xd8, 0xb4, 0xab, 0xc7, 0x73, 0x1f, 0x4d, 0x21, 0x95, 0xf9, 0xe6, 0x8a, 0x3e, 0x52, // 4
    0x00, 0x77, 0xee, 0x99, 0xc7, 0xb0, 0x29, 0x5e, 0x95, 0xe2, 0x7b, 0x0c, 0x52, 0x25, 0xbc, 0xcb, // 5
    0x00, 0x5a, 0xb4, 0xee, 0x73, 0x29, 0xc7, 0x9d, 0xe6, 0xbc, 0x52, 0x08, 0x95, 0xcf, 0x21, 0x7b, // 6
    0x00, 0x41, 0x82, 0xc3, 0x1f, 0x5e, 0x9d, 0xdc, 0x3e, 0x7f, 0xbc, 0xfd, 0x21, 0x60, 0xa3, 0xe2, // 7
    0x00, 0xd8, 0xab, 0x73, 0x4d, 0x95, 0xe6, 0x3e, 0x9a, 0x42, 0x31, 0xe9, 0xd7, 0x0f, 0x7c, 0xa4, // 8
    0x00, 0xc3, 0x9d, 0x5e, 0x21, 0xe2, 0xbc, 0x7f, 0x42, 0x81, 0xdf, 0x1c, 0x63, 0xa0, 0xfe, 0x3d, // 9
    0x00, 0xee, 0xc7, 0x29, 0x95, 0x7b, 0x52, 0xbc, 0x31, 0xdf, 0xf6, 0x18, 0xa4, 0x4a, 0x63, 0x8d, // A
    0x00, 0xf5, 0xf1, 0x04, 0xf9, 0x0c, 0x08, 0xfd, 0xe9, 0x1c, 0x18, 0xed, 0x10, 0xe5, 0xe1, 0x14, // B
    0x00, 0xb4, 0x73, 0xc7, 0xe6, 0x52, 0x95, 0x21, 0xd7, 0x63, 0xa4, 0x10, 0x31, 0x85, 0x42, 0xf6, // C
    0x00, 0xaf, 0x45, 0xea, 0x8a, 0x25, 0xcf, 0x60, 0x0f, 0xa0, 0x4a, 0xe5, 0x85, 0x2a, 0xc0, 0x6f, // D
    0x00, 0x82, 0x1f, 0x9d, 0x3e, 0xbc, 0x21, 0xa3, 0x7c, 0xfe, 0x63, 0xe1, 0x42, 0xc0, 0x5d, 0xdf, // E
    0x00, 0x99, 0x29, 0xb0, 0x52, 0xcb, 0x7b, 0xe2, 0xa4, 0x3d, 0x8d, 0x14, 0xf6, 0x6f, 0xdf, 0x46, // F
};
static const unsigned char AmBl[16][16] = {
    // 0    1     2     3     4     5     6     7     8     9     A     B     C     D     E     F
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // 0
    0x00, 0x10, 0x20, 0x30, 0x40, 0x50, 0x60, 0x70, 0x80, 0x90, 0xa0, 0xb0, 0xc0, 0xd0, 0xe0, 0xf0, // 1
    0x00, 0x20, 0x40, 0x60, 0x80, 0xa0, 0xc0, 0xe0, 0x1b, 0x3b, 0x5b, 0x7b, 0x9b, 0xbb, 0xdb, 0xfb, // 2
    0x00, 0x30, 0x60, 0x50, 0xc0, 0xf0, 0xa0, 0x90, 0x9b, 0xab, 0xfb, 0xcb, 0x5b, 0x6b, 0x3b, 0x0b, // 3
    0x00, 0x40, 0x80, 0xc0, 0x1b, 0x5b, 0x9b, 0xdb, 0x36, 0x76, 0xb6, 0xf6, 0x2d, 0x6d, 0xad, 0xed, // 4
    0x00, 0x50, 0xa0, 0xf0, 0x5b, 0x0b, 0xfb, 0xab, 0xb6, 0xe6, 0x16, 0x46, 0xed, 0xbd, 0x4d, 0x1d, // 5
    0x00, 0x60, 0xc0, 0xa0, 0x9b, 0xfb, 0x5b, 0x3b, 0x2d, 0x4d, 0xed, 0x8d, 0xb6, 0xd6, 0x76, 0x16, // 6
    0x00, 0x70, 0xe0, 0x90, 0xdb, 0xab, 0x3b, 0x4b, 0xad, 0xdd, 0x4d, 0x3d, 0x76, 0x06, 0x96, 0xe6, // 7
    0x00, 0x80, 0x1b, 0x9b, 0x36, 0xb6, 0x2d, 0xad, 0x6c, 0xec, 0x77, 0xf7, 0x5a, 0xda, 0x41, 0xc1, // 8
    0x00, 0x90, 0x3b, 0xab, 0x76, 0xe6, 0x4d, 0xdd, 0xec, 0x7c, 0xd7, 0x47, 0x9a, 0x0a, 0xa1, 0x31, // 9
    0x00, 0xa0, 0x5b, 0xfb, 0xb6, 0x16, 0xed, 0x4d, 0x77, 0xd7, 0x2c, 0x8c, 0xc1, 0x61, 0x9a, 0x3a, // A
    0x00, 0xb0, 0x7b, 0xcb, 0xf6, 0x46, 0x8d, 0x3d, 0xf7, 0x47, 0x8c, 0x3c, 0x01, 0xb1, 0x7a, 0xca, // B
    0x00, 0xc0, 0x9b, 0x5b, 0x2d, 0xed, 0xb6, 0x76, 0x5a, 0x9a, 0xc1, 0x01, 0x77, 0xb7, 0xec, 0x2c, // C
    0x00, 0xd0, 0xbb, 0x6b, 0x6d, 0xbd, 0xd6, 0x06, 0xda, 0x0a, 0x61, 0xb1, 0xb7, 0x67, 0x0c, 0xdc, // D
    0x00, 0xe0, 0xdb, 0x3b, 0xad, 0x4d, 0x76, 0x96, 0x41, 0xa1, 0x9a, 0x7a, 0xec, 0x0c, 0x37, 0xd7, // E
    0x00, 0xf0, 0xfb, 0x0b, 0xed, 0x1d, 0x16, 0xe6, 0xc1, 0x31, 0x3a, 0xca, 0x2c, 0xdc, 0xd7, 0x27, // F
};
static const unsigned char AlBm[16][16] = {
    // 0    1     2     3     4     5     6     7     8     9     A     B     C     D     E     F
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // 0
    0x00, 0x10, 0x20, 0x30, 0x40, 0x50, 0x60, 0x70, 0x80, 0x90, 0xa0, 0xb0, 0xc0, 0xd0, 0xe0, 0xf0, // 1
    0x00, 0x20, 0x40, 0x60, 0x80, 0xa0, 0xc0, 0xe0, 0x1b, 0x3b, 0x5b, 0x7b, 0x9b, 0xbb, 0xdb, 0xfb, // 2
    0x00, 0x30, 0x60, 0x50, 0xc0, 0xf0, 0xa0, 0x90, 0x9b, 0xab, 0xfb, 0xcb, 0x5b, 0x6b, 0x3b, 0x0b, // 3
    0x00, 0x40, 0x80, 0xc0, 0x1b, 0x5b, 0x9b, 0xdb, 0x36, 0x76, 0xb6, 0xf6, 0x2d, 0x6d, 0xad, 0xed, // 4
    0x00, 0x50, 0xa0, 0xf0, 0x5b, 0x0b, 0xfb, 0xab, 0xb6, 0xe6, 0x16, 0x46, 0xed, 0xbd, 0x4d, 0x1d, // 5
    0x00, 0x60, 0xc0, 0xa0, 0x9b, 0xfb, 0x5b, 0x3b, 0x2d, 0x4d, 0xed, 0x8d, 0xb6, 0xd6, 0x76, 0x16, // 6
    0x00, 0x70, 0xe0, 0x90, 0xdb, 0xab, 0x3b, 0x4b, 0xad, 0xdd, 0x4d, 0x3d, 0x76, 0x06, 0x96, 0xe6, // 7
    0x00, 0x80, 0x1b, 0x9b, 0x36, 0xb6, 0x2d, 0xad, 0x6c, 0xec, 0x77, 0xf7, 0x5a, 0xda, 0x41, 0xc1, // 8
    0x00, 0x90, 0x3b, 0xab, 0x76, 0xe6, 0x4d, 0xdd, 0xec, 0x7c, 0xd7, 0x47, 0x9a, 0x0a, 0xa1, 0x31, // 9
    0x00, 0xa0, 0x5b, 0xfb, 0xb6, 0x16, 0xed, 0x4d, 0x77, 0xd7, 0x2c, 0x8c, 0xc1, 0x61, 0x9a, 0x3a, // A
    0x00, 0xb0, 0x7b, 0xcb, 0xf6, 0x46, 0x8d, 0x3d, 0xf7, 0x47, 0x8c, 0x3c, 0x01, 0xb1, 0x7a, 0xca, // B
    0x00, 0xc0, 0x9b, 0x5b, 0x2d, 0xed, 0xb6, 0x76, 0x5a, 0x9a, 0xc1, 0x01, 0x77, 0xb7, 0xec, 0x2c, // C
    0x00, 0xd0, 0xbb, 0x6b, 0x6d, 0xbd, 0xd6, 0x06, 0xda, 0x0a, 0x61, 0xb1, 0xb7, 0x67, 0x0c, 0xdc, // D
    0x00, 0xe0, 0xdb, 0x3b, 0xad, 0x4d, 0x76, 0x96, 0x41, 0xa1, 0x9a, 0x7a, 0xec, 0x0c, 0x37, 0xd7, // E
    0x00, 0xf0, 0xfb, 0x0b, 0xed, 0x1d, 0x16, 0xe6, 0xc1, 0x31, 0x3a, 0xca, 0x2c, 0xdc, 0xd7, 0x27, // F
};
static const unsigned char AlBl[16][16] = {
    // 0    1     2     3     4     5     6     7     8     9     A     B     C     D     E     F
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // 0
    0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, // 1
    0x00, 0x02, 0x04, 0x06, 0x08, 0x0a, 0x0c, 0x0e, 0x10, 0x12, 0x14, 0x16, 0x18, 0x1a, 0x1c, 0x1e, // 2
    0x00, 0x03, 0x06, 0x05, 0x0c, 0x0f, 0x0a, 0x09, 0x18, 0x1b, 0x1e, 0x1d, 0x14, 0x17, 0x12, 0x11, // 3
    0x00, 0x04, 0x08, 0x0c, 0x10, 0x14, 0x18, 0x1c, 0x20, 0x24, 0x28, 0x2c, 0x30, 0x34, 0x38, 0x3c, // 4
    0x00, 0x05, 0x0a, 0x0f, 0x14, 0x11, 0x1e, 0x1b, 0x28, 0x2d, 0x22, 0x27, 0x3c, 0x39, 0x36, 0x33, // 5
    0x00, 0x06, 0x0c, 0x0a, 0x18, 0x1e, 0x14, 0x12, 0x30, 0x36, 0x3c, 0x3a, 0x28, 0x2e, 0x24, 0x22, // 6
    0x00, 0x07, 0x0e, 0x09, 0x1c, 0x1b, 0x12, 0x15, 0x38, 0x3f, 0x36, 0x31, 0x24, 0x23, 0x2a, 0x2d, // 7
    0x00, 0x08, 0x10, 0x18, 0x20, 0x28, 0x30, 0x38, 0x40, 0x48, 0x50, 0x58, 0x60, 0x68, 0x70, 0x78, // 8
    0x00, 0x09, 0x12, 0x1b, 0x24, 0x2d, 0x36, 0x3f, 0x48, 0x41, 0x5a, 0x53, 0x6c, 0x65, 0x7e, 0x77, // 9
    0x00, 0x0a, 0x14, 0x1e, 0x28, 0x22, 0x3c, 0x36, 0x50, 0x5a, 0x44, 0x4e, 0x78, 0x72, 0x6c, 0x66, // A
    0x00, 0x0b, 0x16, 0x1d, 0x2c, 0x27, 0x3a, 0x31, 0x58, 0x53, 0x4e, 0x45, 0x74, 0x7f, 0x62, 0x69, // B
    0x00, 0x0c, 0x18, 0x14, 0x30, 0x3c, 0x28, 0x24, 0x60, 0x6c, 0x78, 0x74, 0x50, 0x5c, 0x48, 0x44, // C
    0x00, 0x0d, 0x1a, 0x17, 0x34, 0x39, 0x2e, 0x23, 0x68, 0x65, 0x72, 0x7f, 0x5c, 0x51, 0x46, 0x4b, // D
    0x00, 0x0e, 0x1c, 0x12, 0x38, 0x36, 0x24, 0x2a, 0x70, 0x7e, 0x6c, 0x62, 0x48, 0x46, 0x54, 0x5a, // E
    0x00, 0x0f, 0x1e, 0x11, 0x3c, 0x33, 0x22, 0x2d, 0x78, 0x77, 0x66, 0x69, 0x44, 0x4b, 0x5a, 0x55, // F
};
// 生成S盒的tlwe加密
// A*B = (Am+Al) * (Bm+Bl) = AmBm + AmBl + AlBm + AlBl
// 所以tlwe(A*B) = tlwe(AmBm) + tlwe(AmBl) + tlwe(AlBm) + tlwe(AlBl)
void MakeAmBmTable(std::vector<TRLWE_1> &Table, const TFHEpp::Key<privksP::targetP> &key)
{
    // dec-> bit
    int Sbox_binary[256][8];
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            int bin_str[8];
            HexToBinStr(AmBm[i][j], bin_str, 8);
            for (int k = 0; k < 8; k++)
            {
                Sbox_binary[i * 16 + j][k] = bin_str[k];
            }
        }
    }
    for (int k = 0; k < 2; k++)
    {
        TFHEpp::Polynomial<typename privksP::targetP> poly;
        for (int i = 0; i < 128; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                poly[i * 8 + j] = static_cast<typename privksP::targetP::T>(Sbox_binary[k * 128 + i][j]);
            }
        }
        Table[k] = TFHEpp::trlweSymIntEncrypt<privksP::targetP>(poly, privksP::targetP::alpha, key);
    }
}
//
void MakeAmBlTable(std::vector<TRLWE_1> &Table, const TFHEpp::Key<privksP::targetP> &key)
{
    // dec-> bit
    int Sbox_binary[256][8];
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            int bin_str[8];
            HexToBinStr(AmBl[i][j], bin_str, 8);
            for (int k = 0; k < 8; k++)
            {
                Sbox_binary[i * 16 + j][k] = bin_str[k];
            }
        }
    }
    for (int k = 0; k < 2; k++)
    {
        TFHEpp::Polynomial<typename privksP::targetP> poly;
        for (int i = 0; i < 128; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                poly[i * 8 + j] = static_cast<typename privksP::targetP::T>(Sbox_binary[k * 128 + i][j]);
            }
        }
        Table[k] = TFHEpp::trlweSymIntEncrypt<privksP::targetP>(poly, privksP::targetP::alpha, key);
    }
}
//
void MakeAlBmTable(std::vector<TRLWE_1> &Table, const TFHEpp::Key<privksP::targetP> &key)
{
    // dec-> bit
    int Sbox_binary[256][8];
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            int bin_str[8];
            HexToBinStr(AlBm[i][j], bin_str, 8);
            for (int k = 0; k < 8; k++)
            {
                Sbox_binary[i * 16 + j][k] = bin_str[k];
            }
        }
    }
    for (int k = 0; k < 2; k++)
    {
        TFHEpp::Polynomial<typename privksP::targetP> poly;
        for (int i = 0; i < 128; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                poly[i * 8 + j] = static_cast<typename privksP::targetP::T>(Sbox_binary[k * 128 + i][j]);
            }
        }
        Table[k] = TFHEpp::trlweSymIntEncrypt<privksP::targetP>(poly, privksP::targetP::alpha, key);
    }
}
//
void MakeAlBlTable(std::vector<TRLWE_1> &Table, const TFHEpp::Key<privksP::targetP> &key)
{
    // dec-> bit
    int Sbox_binary[256][8];
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            int bin_str[8];
            HexToBinStr(AlBl[i][j], bin_str, 8);
            for (int k = 0; k < 8; k++)
            {
                Sbox_binary[i * 16 + j][k] = bin_str[k];
            }
        }
    }
    for (int k = 0; k < 2; k++)
    {
        TFHEpp::Polynomial<typename privksP::targetP> poly;
        for (int i = 0; i < 128; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                poly[i * 8 + j] = static_cast<typename privksP::targetP::T>(Sbox_binary[k * 128 + i][j]);
            }
        }
        Table[k] = TFHEpp::trlweSymIntEncrypt<privksP::targetP>(poly, privksP::targetP::alpha, key);
    }
}
// 字节替换混合打包实现
// 计算 NX2
const privksP::targetP::T NX2 = 2 * privksP::targetP::n;
// 定义全局常量 bara
const privksP::targetP::T bara[7] = {NX2 - 8 * (1 << 0), NX2 - 8 * (1 << 1), NX2 - 8 * (1 << 2), NX2 - 8 * (1 << 3), NX2 - 8 * (1 << 4), NX2 - 8 * (1 << 5), NX2 - 8 * (1 << 6)};
void CipherSubBytesMixedPacking(TRLWE_1 &result, std::vector<TRLWE_1> &Table, std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> &select)
{
    TFHEpp::CMUXFFT<typename privksP::targetP>(result, select[7], Table[1], Table[0]);
    TFHEpp::BlindRotate_LUT<privksP>(result, bara, select, 7); //, resultOfCMUX);
}
// 8比特乘法的查表实现
void tfhe_mul(TRLWE_1 &result1, TRLWE_1 &result2, TRLWE_1 &result3, TRLWE_1 &result4,
              std::vector<TRLWE_1> &Table1, std::vector<TRLWE_1> &Table2, std::vector<TRLWE_1> &Table3, std::vector<TRLWE_1> &Table4,
              std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> &select1, std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> &select2)
{
    //
    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> ambm(8);
    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> ambl(8);
    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> albm(8);
    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> albl(8);
    // 将select1的后四位和select2的后四位值赋给ambm，即ambm={select2[4],select2[5],select2[6],select2[7],select1[4],select1[5],select1[6],select1[7]}
    for (int i = 0; i < 4; i++)
    {
        ambm[i] = select2[i + 4];
        ambm[i + 4] = select1[i + 4];
    }

    // 将select1的前四位和select2的后四位值赋给ambl，即ambl={select2[4],select2[5],select2[6],select2[7],select1[0],select1[1],select1[2],select1[3]}
    for (int i = 0; i < 4; i++)
    {
        ambl[i] = select2[i + 4];
        ambl[i + 4] = select1[i];
    }
    // 将select1的后四位和select2的前四位值赋给albm，即albm={select2[0],select2[1],select2[2],select2[3],select1[4],select1[5],select1[6],select1[7]}
    for (int i = 0; i < 4; i++)
    {
        albm[i] = select2[i];
        albm[i + 4] = select1[i + 4];
    }
    // 将select1的前四位和select2的前四位值赋给albl，即albl={select2[0],select2[1],select2[2],select2[3],select1[0],select1[1],select1[2],select1[3]}
    for (int i = 0; i < 4; i++)
    {
        albl[i] = select2[i];
        albl[i + 4] = select1[i];
    }
    CipherSubBytesMixedPacking(result1, Table1, ambm);
    CipherSubBytesMixedPacking(result2, Table2, ambl);
    CipherSubBytesMixedPacking(result3, Table3, albm);
    CipherSubBytesMixedPacking(result4, Table4, albl);
}
// 解密 TLWE 密文的函数
void decryptTLWE(const std::vector<std::vector<TLWE_0>> &cipher, const TFHEpp::SecretKey &sk, std::array<unsigned char, 16> &decrypt)
{
    for (int i = 0; i < 16; i++)
    {
        int dec_bin[8]; // 存储解密后的二进制数组

        // 解密每个字节的 8 个比特位
        for (int j = 0; j < 8; j++)
        {
            dec_bin[j] = TFHEpp::tlweSymIntDecrypt<typename bkP::domainP>(cipher[i][j], sk.key.get<typename bkP::domainP>());
        }

        // 将二进制数组转换为对应的字符（unsigned char）
        int temp_int = 0;                                  // 用来存储解密后的整数
        BinStrToHex(dec_bin, temp_int, 8);                 // 将二进制转换为整数
        decrypt[i] = static_cast<unsigned char>(temp_int); // 存储为解密后的 unsigned char
    }
}
// 检查tlwe解密结果是否正确
int check_decrypt(std::array<unsigned char, 16> &decrypt, std::array<unsigned char, 16> &plain)
{
    int flag = 1;
    for (int i = 0; i < 16; i++)
    {
        if (decrypt[i] != plain[i])
        {
            flag = 0;
            break;
        }
    }
    return flag;
}
// 明文的有限域乘法
unsigned char mul(unsigned char a, unsigned char b) {
    unsigned char p = 0;
    unsigned char counter;
    unsigned char hi_bit_set;
    for(counter = 0; counter < 8; counter++) {
          if((b & 1) == 1) 
                p ^= a;
          hi_bit_set = (a & 0x80);
          a <<= 1;
          if(hi_bit_set == 0x80) 
                a ^= 0x1b;          
          b >>= 1;
      }
      return p;
}
int main()
{
// 生成tlwe密钥，自举密钥，key-switching密钥等等
#if 1
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    // Generate key
    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey; // 用于给s盒进行TRLWE加密
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    ek.emplacebkfft<TFHEpp::lvl01param>(*sk); // used for identitybootstrapping
#endif

// 对四个S盒进行tlwe加密
#if 1
    std::vector<TRLWE_1> Table1(2);
    std::vector<TRLWE_1> Table2(2);
    std::vector<TRLWE_1> Table3(2);
    std::vector<TRLWE_1> Table4(2);
    MakeAmBmTable(Table1, sk->key.get<privksP::targetP>());
    MakeAmBlTable(Table2, sk->key.get<privksP::targetP>());
    MakeAlBmTable(Table3, sk->key.get<privksP::targetP>());
    MakeAlBlTable(Table4, sk->key.get<privksP::targetP>());
#endif

// 生成consByte
#if 1
    std::vector<TLWE_0> consByte(8);
    for (int i = 0; i < 8; i++)
    {
        // encrypt 0
        consByte[i] = TFHEpp::tlweSymIntEncrypt<typename bkP::domainP>((typename bkP::domainP::T)byte_mul2[i], bkP::domainP::alpha,
                                                                       sk->key.get<typename bkP::domainP>());
    }
#endif

    // 明文
    std::array<unsigned char, 16> plain = {0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34};
    // 明文运算
    std::array<unsigned char, 16> plain2;
    for (int i = 0; i < 16; i++)
    {
        plain2[i] = mul(plain[i], plain[(i+1)%16]);
    }
// 对密文进行tlwe加密
#if 1
    std::cout << "明文tlwe加密开始" << endl;
    std::vector<std::vector<TLWE_0>> cipher;
    cipher.resize(16);
    for (int i = 0; i < 16; i++)
    {
        cipher[i].resize(8);
    }
    // #pragma omp parallel for
    for (int i = 0; i < 16; i++)
    {
        int bin_str[8];
        HexToBinStr(plain[i], bin_str, 8);
        // #pragma omp parallel for num_threads(8)
        for (int j = 0; j < 8; j++)
        {
            cipher[i][j] = TFHEpp::tlweSymIntEncrypt<typename bkP::domainP>((typename bkP::domainP::T)bin_str[j], bkP::domainP::alpha,
                                                                            sk->key.get<typename bkP::domainP>());
        }
    }
    std::cout << "明文tlwe加密完成" << endl;
#endif

// 密文的tlwe解密验证
// 对cipher进行解密，与state进行比较
#if 1
    std::array<unsigned char, 16> decrypt;
    std::cout << "密文tlwe解密开始" << endl;
    decryptTLWE(cipher, *sk, decrypt);
    std::cout << "密文tlwe解密完成" << endl;
    int flag = check_decrypt(decrypt, plain);

    if (flag)
    {
        std::cout << "密文tlwe解密成功!!!" << endl;
    }
    else
    {
        std::cout << "密文tlwe解密失败" << endl;
        return 0;
    }
#endif

    // 前面已经用tlwe对plain进行加密得到了cipher，现在求每个元素的平方

    std::vector<std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>>> bootedTRGSW;
    bootedTRGSW.resize(16);
    for (int i = 0; i < 16; i++)
    {
        bootedTRGSW[i].resize(8);
    }

    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 8; j++)
            TFHEpp::SM4_CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTRGSW[i][j],
                                                                    cipher[i][j], ek);
    }
    std::vector<TRLWE_1> lut_result1(16); //
    std::vector<TRLWE_1> lut_result2(16); //
    std::vector<TRLWE_1> lut_result3(16); //
    std::vector<TRLWE_1> lut_result4(16); //
    std::cout << "密文平方tlwe加密开始" << endl;
    for (int i = 0; i < 16; i++)
    {
        tfhe_mul(lut_result1[i], lut_result2[i], lut_result3[i], lut_result4[i], Table1, Table2, Table3, Table4, bootedTRGSW[i], bootedTRGSW[(i+1)%16]);
    }
    // SampleExtract level 1
    std::vector<std::vector<TLWE_1>> Sbox_value1;
    std::vector<std::vector<TLWE_1>> Sbox_value2;
    std::vector<std::vector<TLWE_1>> Sbox_value3;
    std::vector<std::vector<TLWE_1>> Sbox_value4;
    Sbox_value1.resize(16);
    Sbox_value2.resize(16);
    Sbox_value3.resize(16);
    Sbox_value4.resize(16);
    for (int i = 0; i < Sbox_value1.size(); i++)
    {
        Sbox_value1[i].resize(8);
        Sbox_value2[i].resize(8);
        Sbox_value3[i].resize(8);
        Sbox_value4[i].resize(8);
    }
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            TFHEpp::SampleExtractIndex<typename privksP::targetP>(Sbox_value1[i][j], lut_result1[i], j);
            TFHEpp::SampleExtractIndex<typename privksP::targetP>(Sbox_value2[i][j], lut_result2[i], j);
            TFHEpp::SampleExtractIndex<typename privksP::targetP>(Sbox_value3[i][j], lut_result3[i], j);
            TFHEpp::SampleExtractIndex<typename privksP::targetP>(Sbox_value4[i][j], lut_result4[i], j);
        }
    }
    // Key Switch to LWE B  on level 0
    std::vector<std::vector<TLWE_0>> cipher1;
    std::vector<std::vector<TLWE_0>> cipher2;
    std::vector<std::vector<TLWE_0>> cipher3;
    std::vector<std::vector<TLWE_0>> cipher4;
    cipher1.resize(16);
    cipher2.resize(16);
    cipher3.resize(16);
    cipher4.resize(16);
    for (int i = 0; i < 16; i++)
    {
        cipher1[i].resize(8);
        cipher2[i].resize(8);
        cipher3[i].resize(8);
        cipher4[i].resize(8);
    }

    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            // level 1 -> level 0
            TFHEpp::IdentityKeySwitch<iksP>(cipher1[i][j], Sbox_value1[i][j], ek.getiksk<iksP>());
            TFHEpp::IdentityKeySwitch<iksP>(cipher2[i][j], Sbox_value2[i][j], ek.getiksk<iksP>());
            TFHEpp::IdentityKeySwitch<iksP>(cipher3[i][j], Sbox_value3[i][j], ek.getiksk<iksP>());
            TFHEpp::IdentityKeySwitch<iksP>(cipher4[i][j], Sbox_value4[i][j], ek.getiksk<iksP>());
        }
    }
    for(int i=0;i<16;i++){
        XOR_Four(cipher[i],cipher1[i],cipher2[i],cipher3[i],cipher4[i]);
    }
    // tfhe平方完成
    std::cout << "密文平方tlwe加密完成" << endl;
    // 调用解密函数
    decryptTLWE(cipher, *sk, decrypt);
    // 检查解密结果是否正确
    flag = check_decrypt(decrypt, plain2);
    if (flag)
    {
        std::cout << "密文平方tlwe解密成功!!!" << endl;
    }
    else
    {
        std::cout << "密文平方tlwe解密失败" << endl;
        return 0;
    }
    return 0;
}
//8bit乘法那里应该可以优化，目前是先用分成4个S盒，然后分别对四个值密钥转换，最后进行tlwe加法
// 这样效率有点低，可以考虑s盒查表后对四个结果进行trgsw加法，最后进行密钥转换