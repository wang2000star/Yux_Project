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
// 声明并初始化全局变量
std::vector<std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>>> bootedTRGSW(16, std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>>(8));
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
    for (int k = 0; k < 2; k++)
    {
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
void XOR_Four(P &result, P &a, P &b, P &c, P &d, int N)
{
    XOR_Two<P>(result, a, b, N);
    XOR_Two<P>(result, result, c, N);
    XOR_Two<P>(result, result, d, N);
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
// 恒等变换表
void MakeIdentTable(std::vector<TRLWE_1> &Table, const TFHEpp::Key<privksP::targetP> &key)
{
    // dec-> bit
    unsigned char Sbox_Ident[16][16];
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            Sbox_Ident[i][j] = i * 16 + j;
        }
    }
    int Sbox_binary[256][8];
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            int bin_str[8];
            HexToBinStr(Sbox_Ident[i][j], bin_str, 8);
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
void tfhe_mul(TRLWE_1 &result,
              std::vector<TRLWE_1> &Table1, std::vector<TRLWE_1> &Table2, std::vector<TRLWE_1> &Table3, std::vector<TRLWE_1> &Table4,
              std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> &select1, std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> &select2)
{
    // 分割select1和select2
    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> ambm(8);
    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> ambl(8);
    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> albm(8);
    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> albl(8);
    // 将select1的后四位和select2的后四位值赋给ambm，即ambm={select2[4],select2[5],select2[6],select2[7],select1[4],select1[5],select1[6],select1[7]}
    // 将select1的前四位和select2的后四位值赋给ambl，即ambl={select2[4],select2[5],select2[6],select2[7],select1[0],select1[1],select1[2],select1[3]}
    // 将select1的后四位和select2的前四位值赋给albm，即albm={select2[0],select2[1],select2[2],select2[3],select1[4],select1[5],select1[6],select1[7]}
    // 将select1的前四位和select2的前四位值赋给albl，即albl={select2[0],select2[1],select2[2],select2[3],select1[0],select1[1],select1[2],select1[3]}
    for (int i = 0; i < 4; i++)
    {
        ambm[i] = select2[i + 4];
        ambm[i + 4] = select1[i + 4];
        ambl[i] = select2[i + 4];
        ambl[i + 4] = select1[i];
        albm[i] = select2[i];
        albm[i + 4] = select1[i + 4];
        albl[i] = select2[i];
        albl[i + 4] = select1[i];
    }
    // 计算四个表的乘法
    TRLWE_1 result1;
    TRLWE_1 result2;
    TRLWE_1 result3;
    TRLWE_1 result4;
    CipherSubBytesMixedPacking(result1, Table1, ambm);
    CipherSubBytesMixedPacking(result2, Table2, ambl);
    CipherSubBytesMixedPacking(result3, Table3, albm);
    CipherSubBytesMixedPacking(result4, Table4, albl);
    TRLWE1_XOR_Four(result, result1, result2, result3, result4);
}
// 实现s盒
void Sbox(std::vector<std::vector<TLWE_0>> &cipher,
          TRLWE_1 &RoundConstantTRLWE,
          std::vector<TRLWE_1> &Table1, std::vector<TRLWE_1> &Table2, std::vector<TRLWE_1> &Table3, std::vector<TRLWE_1> &Table4, std::vector<TRLWE_1> &Table5,
          TFHEpp::EvalKey &ek)
{

    std::vector<std::vector<TLWE_1>> Sbox_value;
    Sbox_value.resize(16);
    for (int i = 0; i < 16; i++)
    {
        Sbox_value[i].resize(8);
    }
    std::vector<std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>>> bootedTRGSW;
    bootedTRGSW.resize(16);
    for (int i = 0; i < 16; i++)
    {
        bootedTRGSW[i].resize(8);
    }
    std::vector<TRLWE_1> lut_result(16);
    int index[16] = {2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13};
    for (int r = 0; r < 2; r++)
    {
#pragma omp parallel for
        for (int i = 0; i < 16; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                TFHEpp::SM4_CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTRGSW[i][j], cipher[i][j], ek);
            }
        }

        for (int i = 0; i < 16; i++)
        {
            // 0123变为1230，4567变为5674，89ab变为ab89，cdef变为efcd
            CipherSubBytesMixedPacking(lut_result[i], Table5, bootedTRGSW[index[i]]);
        }
        std::vector<TRLWE_1> lut_result0(8);
        tfhe_mul(lut_result0[0], Table1, Table2, Table3, Table4, bootedTRGSW[1], bootedTRGSW[2]);
        tfhe_mul(lut_result0[1], Table1, Table2, Table3, Table4, bootedTRGSW[2], bootedTRGSW[3]);
        tfhe_mul(lut_result0[2], Table1, Table2, Table3, Table4, bootedTRGSW[5], bootedTRGSW[6]);
        tfhe_mul(lut_result0[3], Table1, Table2, Table3, Table4, bootedTRGSW[6], bootedTRGSW[7]);
        tfhe_mul(lut_result0[4], Table1, Table2, Table3, Table4, bootedTRGSW[9], bootedTRGSW[10]);
        tfhe_mul(lut_result0[5], Table1, Table2, Table3, Table4, bootedTRGSW[10], bootedTRGSW[11]);
        tfhe_mul(lut_result0[6], Table1, Table2, Table3, Table4, bootedTRGSW[13], bootedTRGSW[14]);
        tfhe_mul(lut_result0[7], Table1, Table2, Table3, Table4, bootedTRGSW[14], bootedTRGSW[15]);

#if 1

        TRlWE1_XOR_Two(lut_result[2], lut_result[2], lut_result[1]);
        TRlWE1_XOR_Two(lut_result[2], lut_result[2], lut_result0[0]);
        TRlWE1_XOR_Two(lut_result[6], lut_result[6], lut_result[5]);
        TRlWE1_XOR_Two(lut_result[6], lut_result[6], lut_result0[2]);
        TRlWE1_XOR_Two(lut_result[10], lut_result[10], lut_result[9]);
        TRlWE1_XOR_Two(lut_result[10], lut_result[10], lut_result0[4]);
        TRlWE1_XOR_Two(lut_result[14], lut_result[14], lut_result[13]);
        TRlWE1_XOR_Two(lut_result[14], lut_result[14], lut_result0[6]);     

        TRlWE1_XOR_Two(lut_result[3], lut_result[3], lut_result[2]);
        TRlWE1_XOR_Two(lut_result[3], lut_result[3], lut_result0[1]);
        TRlWE1_XOR_Two(lut_result[7], lut_result[7], lut_result[6]);
        TRlWE1_XOR_Two(lut_result[7], lut_result[7], lut_result0[3]);
        TRlWE1_XOR_Two(lut_result[11], lut_result[11], lut_result[10]);
        TRlWE1_XOR_Two(lut_result[11], lut_result[11], lut_result0[5]);
        TRlWE1_XOR_Two(lut_result[15], lut_result[15], lut_result[14]);
        TRlWE1_XOR_Two(lut_result[15], lut_result[15], lut_result0[7]);

        TRlWE1_XOR_Two(lut_result[2], lut_result[2], RoundConstantTRLWE);
        TRlWE1_XOR_Two(lut_result[6], lut_result[6], RoundConstantTRLWE);
        TRlWE1_XOR_Two(lut_result[10], lut_result[10], RoundConstantTRLWE);
        TRlWE1_XOR_Two(lut_result[14], lut_result[14], RoundConstantTRLWE);

        
#endif
        // SampleExtract level 1
        for (int i = 0; i < 16; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                // SampleExtractIndex
                TFHEpp::SampleExtractIndex<typename privksP::targetP>(Sbox_value[i][j], lut_result[i], j);
                // Key Switch to LWE B on level 0
                TFHEpp::IdentityKeySwitch<iksP>(cipher[i][j], Sbox_value[i][j], ek.getiksk<iksP>());
            }
        }
    }  
}

// decLinearLayer
//  输入x[16]，输出x[16]
//  temp0=x<<<0,temp3=x<<<3,temp4=x<<<4,temp8=x<<<8,temp9=x<<<9,temp12=x<<<12,temp14=x<<<14
//  x=temp0^temp3^temp4^temp8^temp9^temp12^temp14
void cipherinvLinearLayer(std::vector<std::vector<TLWE_0>> &cipher, std::vector<TRLWE_1> &Table5, TFHEpp::EvalKey &ek)
{
    // 开始自举
#if 1
    std::vector<std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>>> bootedTRGSW;
    bootedTRGSW.resize(16);
    for (int i = 0; i < 16; i++)
    {
        bootedTRGSW[i].resize(8);
    }
#pragma omp parallel for
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            TFHEpp::SM4_CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTRGSW[i][j], cipher[i][j], ek);
        }
    }
    std::vector<TRLWE_1> lut_result(16); //
    for (int i = 0; i < 16; i++)
    {
        CipherSubBytesMixedPacking(lut_result[i], Table5, bootedTRGSW[i]);
    }
    // SampleExtract level 1
    std::vector<std::vector<TLWE_1>> Sbox_value;
    Sbox_value.resize(16);
    for (int i = 0; i < Sbox_value.size(); i++)
    {
        Sbox_value[i].resize(8);
    }
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            TFHEpp::SampleExtractIndex<typename privksP::targetP>(Sbox_value[i][j], lut_result[i], j);
        }
    }
    // Key Switch to LWE B  on level 0
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            // level 1 -> level 0
            TFHEpp::IdentityKeySwitch<iksP>(cipher[i][j], Sbox_value[i][j], ek.getiksk<iksP>());
        }
    }
#endif
    //
    std::vector<std::vector<TLWE_0>> temp;
    temp.resize(16);
    for (int i = 0; i < 16; i++)
    {
        temp[i].resize(8);
    }
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            TFHEpp::HomCOPY<typename bkP::domainP>(temp[i][j], cipher[i][j]);
        }
    }
    std::vector<TLWE_0> t(8);
    for (int i = 0; i < 16; i++)
    {
        XOR_Four(t, temp[i], temp[(i + 3) % 16], temp[(i + 4) % 16], temp[(i + 8) % 16], 8);
        XOR_Four(cipher[i], t, temp[(i + 9) % 16], temp[(i + 12) % 16], temp[(i + 14) % 16], 8);
    }
}
// 密钥扩展的tlwe加密
template <typename bkP>
void encryptRoundKeys(const unsigned char *RoundKey, std::vector<std::vector<TFHEpp::TLWE<typename bkP::domainP>>> &rk, const TFHEpp::SecretKey &sk, int nRoundKeys)
{
    for (int i = 0; i < 16 * nRoundKeys; i++)
    {
        int bin_str[8];
        rk[i].resize(8);
        HexToBinStr(RoundKey[i], bin_str, 8);
        for (int k = 0; k < 8; k++)
        {
            // encrypt TLWE in level 0
            rk[i][k] = TFHEpp::tlweSymIntEncrypt<typename bkP::domainP>(
                static_cast<typename bkP::domainP::T>(bin_str[k]),
                bkP::domainP::alpha,
                sk.key.get<typename bkP::domainP>());
        }
    }
}
// 解密 TLWE 密文的函数
void decryptTLWE(const std::vector<std::vector<TLWE_0>> &cipher, const TFHEpp::SecretKey &sk, unsigned char decrypt[16])
{
#pragma omp parallel for
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
int check_decrypt(unsigned char decrypt[16], unsigned char plain[16])
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
// 函数用于从二维数组中快速切片
std::vector<std::vector<TLWE_0>> slice_2d(const std::vector<std::vector<TLWE_0>> &matrix, int start_row)
{
    std::vector<std::vector<TLWE_0>> slice; // 创建切片的目标容器
    slice.resize(16);
    for (int i = 0; i < 16; i++)
    {
        slice[i].resize(8);
    }
#pragma omp parallel for
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            // 用copy
            TFHEpp::HomCOPY<typename bkP::domainP>(slice[i][j], matrix[start_row + i][j]);
        }
    }
    return slice; // 返回切片
}

extern int KeyExpansion(unsigned char RoundKey[], int ROUND, int blockByte, unsigned char Key[]);
extern void addRoundKey(unsigned char state[], unsigned char RoundKey[], int round);
extern void encSboxFi(unsigned char state[], int begin);
extern void encLinearLayer(unsigned char state[]);
static const unsigned char roundConstant = 0xCD;
extern unsigned char mul(unsigned char a, unsigned char b);
int main()
{

    std::cout << "Program start" << std::endl;
    // omp_set_num_threads(16); // 设置线程数为16
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
// 对5个S盒进行tlwe加密
#if 1
    std::vector<TRLWE_1> Table1(2);
    std::vector<TRLWE_1> Table2(2);
    std::vector<TRLWE_1> Table3(2);
    std::vector<TRLWE_1> Table4(2);
    std::vector<TRLWE_1> Table5(2);
    MakeAmBmTable(Table1, sk->key.get<privksP::targetP>());
    MakeAmBlTable(Table2, sk->key.get<privksP::targetP>());
    MakeAlBmTable(Table3, sk->key.get<privksP::targetP>());
    MakeAlBlTable(Table4, sk->key.get<privksP::targetP>());
    MakeIdentTable(Table5, sk->key.get<privksP::targetP>());
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
    unsigned char plain[16] = {0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34};
    // 初始密钥
    unsigned char key[16] = {0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};
// 密钥扩展
#if 1
    unsigned char RoundKey[208]; // 16*13
    int ROUND = 12;              // Yux加密轮数，含最后一轮，不含白化轮
    int blockByte = 16;
    int nRoundKeys = KeyExpansion(RoundKey, ROUND, blockByte, key);
    // 将RoundKey转为二维数组RoundKey2D
    unsigned char RoundKey2D[13][16];

    for (int i = 0; i < 13; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            RoundKey2D[i][j] = RoundKey[i * 16 + j];
        }
    }

#endif
// 明文Yux加密
#if 1

    // 用轮密钥加密明文
    unsigned char state[16];
    // 令state=plain
    for (int i = 0; i < 16; i++)
    {
        state[i] = plain[i];
    }

    // 定义一个二维数组，来存储Yux每一轮加密后的state,12行16列。不是tlwe密文
    unsigned char Yux_state[3 * ROUND - 1][16];
    // std::cout << "Yux encryption start" << std::endl;

    // 白化轮
    addRoundKey(state, RoundKey, 0); // 0~12
    // 将state赋值给Yux_state
    for (int i = 0; i < 16; i++)
    {
        Yux_state[0][i] = state[i];
    }
    // 11 轮常规 Yux 轮
    for (int r = 1; r < ROUND; r++)
    {
        // S Layer -- 4 sbox
        for (int i = 0; i < 4; i++)
        {
            encSboxFi(state, i * 4);
            encSboxFi(state, i * 4);
            encSboxFi(state, i * 4);
            encSboxFi(state, i * 4);
        }
        for (int i = 0; i < 16; i++)
        {
            Yux_state[3 * r - 2][i] = state[i];
        }
        // linear Layer
        encLinearLayer(state);
        for (int i = 0; i < 16; i++)
        {
            Yux_state[3 * r - 1][i] = state[i];
        }
        addRoundKey(state, RoundKey, r);
        for (int i = 0; i < 16; i++)
        {
            Yux_state[3 * r][i] = state[i];
        }
    }
    // 最后一轮
    for (int i = 0; i < 4; i++)
    {
        encSboxFi(state, i * 4);
        encSboxFi(state, i * 4);
        encSboxFi(state, i * 4);
        encSboxFi(state, i * 4);
    }
    for (int i = 0; i < 16; i++)
    {
        Yux_state[34][i] = state[i];
    }
    addRoundKey(state, RoundKey, 12);
    std::cout << "Yux加密完成" << std::endl;
#endif
// 对密文进行tlwe加密
#if 1
    std::cout << "Yux密文tlwe加密开始" << std::endl;
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
        HexToBinStr(state[i], bin_str, 8);
        for (int j = 0; j < 8; j++)
        {
            cipher[i][j] = TFHEpp::tlweSymIntEncrypt<typename bkP::domainP>((typename bkP::domainP::T)bin_str[j], bkP::domainP::alpha,
                                                                            sk->key.get<typename bkP::domainP>());
        }
    }
    std::cout << "Yux密文tlwe加密完成" << std::endl;
#endif
// 密文的tlwe解密验证
// 对cipher进行解密，与state进行比较
#if 1
    unsigned char decrypt[16];
    std::cout << "密文tlwe解密开始" << std::endl;
    decryptTLWE(cipher, *sk, decrypt);
    std::cout << "密文tlwe解密完成" << std::endl;
    int flag = check_decrypt(decrypt, state);
    if (flag)
    {
        std::cout << "密文tlwe解密成功!!!" << std::endl;
    }
    else
    {
        std::cout << "密文tlwe解密失败" << std::endl;
        return 0;
    }
#endif
// Yux密钥拓展的tlwe加解密
#if 1
    std::vector<std::vector<TLWE_0>> rk;
    rk.resize(16 * nRoundKeys);
    encryptRoundKeys<bkP>(RoundKey, rk, *sk, nRoundKeys);
#endif
    std::cout << "Yux密钥扩展的tlwe加密完成" << std::endl;
/* 前面以得到密钥拓展tlwe加密后的轮密钥rk[208][8]
   下面解密前面的轮密钥rk，得到decrypt_rk[208][8]，
   将其转为每行转为4为hex整数，再变为13行16列，最后与未加密的拓展密钥RoundKey2D进行逐元素比较
*/
#if 1
    std::array<std::array<int, 16>, 13> decrypt_rk;
    std::cout << "Yuxkey的tlwe解密开始" << std::endl;
    // #pragma omp parallel for
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < nRoundKeys; j++)
        {
            int dec_hex = 0;
            int dec_bin[8];
            for (int k = 0; k < 8; k++)
            {
                // typename P::T a = TFHEpp::tlweSymIntDecrypt<typename bkP::domainP>();
                dec_bin[k] = TFHEpp::tlweSymIntDecrypt<typename bkP::domainP>(rk[16 * j + i][k], sk->key.get<typename bkP::domainP>());
                // bootsSymDecrypt(&rk[0][i][j], key);
            }
            // 将二进制转为十六进制
            BinStrToHex(dec_bin, dec_hex, 8);
            decrypt_rk[j][i] = dec_hex;
        }
    }
    std::cout << "Yuxkey的tlwe解密完成" << std::endl;
    // 将解密后的轮密钥与原始轮密钥进行比较
    flag = 1;
    for (int i = 0; i < nRoundKeys; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            if (RoundKey2D[i][j] != decrypt_rk[i][j])
            {
                std::cout << "Yuxkey的tlwe解密失败" << std::endl;
                flag = 0;
                return 0;
            }
        }
    }
    if (flag)
    {
        std::cout << "Yuxkey的tlwe解密成功!!!" << std::endl;
    }
#endif
// 前面已经得到了AES密文且用tlwe对其进行加密得到了cipher，现在使用tlwe下的AES密钥扩展对cipher进行AES解密，这样就间接得到了用tlwe直接加密明文结果
#if 1
    // 将常数roundConstant加密成tlwe;
    std::vector<TLWE_0> RoundConstantTLWE;
    RoundConstantTLWE.resize(8);
    int bin_str[8];
    HexToBinStr(roundConstant, bin_str, 8);
    // #pragma omp parallel for num_threads(8)
    for (int j = 0; j < 8; j++)
    {
        RoundConstantTLWE[j] = TFHEpp::tlweSymIntEncrypt<typename bkP::domainP>((typename bkP::domainP::T)bin_str[j], bkP::domainP::alpha,
                                                                                sk->key.get<typename bkP::domainP>());
    }
    // 将RoundConstantTLWE电路自举为TRLWE
    TRLWE_1 RoundConstantTRLWE;
    std::vector<std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>>> bootedTRGSW_constant;
    bootedTRGSW_constant.resize(1);
    for (int i = 0; i < 1; i++)
    {
        bootedTRGSW_constant[i].resize(8);
    }
    for (int j = 0; j < 8; j++)
    {
        TFHEpp::SM4_CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTRGSW_constant[0][j], RoundConstantTLWE[j], ek);
    }
    CipherSubBytesMixedPacking(RoundConstantTRLWE, Table5, bootedTRGSW_constant[0]);
#endif
// 最后一轮第12轮解密
#if 1
    // 定义总时间
    double total_time = 0;
    double total_time_Sbox = 0;
    double total_time_XOR = 0;
    std::vector<std::vector<TLWE_0>> rk_slice;
    rk_slice.resize(16);
    for (int i = 0; i < 16; i++)
    {
        rk_slice[i].resize(8);
    }
    std::cout << "第" << ROUND << "轮解密开始" << std::endl;
    // 记录开始时间
    auto start = std::chrono::system_clock::now();
    rk_slice = slice_2d(rk, 16 * nRoundKeys - 16);
    //  使用切片进行 XOR_Two 操作
    // #pragma omp parallel for
    for (int i = 0; i < 16; i++)
    {
        XOR_Two(cipher[i], cipher[i], rk_slice[i], 8);
    }
    // 记录解密验证时间
    auto start1 = std::chrono::system_clock::now();
    // 解密第12轮XOR
    decryptTLWE(cipher, *sk, decrypt);
    flag = check_decrypt(decrypt, Yux_state[34]);
    if (flag)
    {
        std::cout << "第" << ROUND << "轮XOR解密成功!!!" << std::endl;
    }
    else
    {
        std::cout << "第" << ROUND << "轮XOR解密失败" << std::endl;
        return 0;
    }
    // 记录解密验证时间
    auto end1 = std::chrono::system_clock::now();
    // 逆S盒操作
    Sbox(cipher, RoundConstantTRLWE, Table1, Table2, Table3, Table4, Table5, ek);
    // 记录解密验证时间
    auto start2 = std::chrono::system_clock::now();
    // 调用解密函数
    decryptTLWE(cipher, *sk, decrypt);
    //输出解密结果
    for (int i = 0; i < 16; i++)
    {
        std::cout << hex << (int)decrypt[i] << " ";
    }
    std::cout << std::endl;
    // 输出Yux_state[33]
    for (int i = 0; i < 16; i++)
    {
        std::cout << hex << (int)Yux_state[33][i] << " ";
    }
    std::cout << std::endl;
    // 检查解密结果是否正确
    flag = check_decrypt(decrypt, Yux_state[33]);
    if (flag)
    {
        std::cout << "第" << ROUND << "轮Sbox解密成功!!!" << std::endl;
    }
    else
    {
        std::cout << "第" << ROUND << "轮Sbox解密失败" << std::endl;
        // return 0;
    }
    // 记录解密验证时间
    auto end2 = std::chrono::system_clock::now();
    // 计算XOR时间
    std::chrono::duration<double> elapsed_seconds_XOR = start1 - start;
    // 计算Sbox时间
    std::chrono::duration<double> elapsed_seconds_Sbox = start2 - end1;
    // 计算解密时间
    std::chrono::duration<double> elapsed_seconds = elapsed_seconds_Sbox + elapsed_seconds_XOR;
    std::cout << "第" << ROUND << "轮解密完成" << std::endl;
    std::cout << "第" << ROUND << "轮XOR用时：" << elapsed_seconds_XOR.count() << "s" << std::endl;
    std::cout << "第" << ROUND << "轮Sbox用时：" << elapsed_seconds_Sbox.count() << "s" << std::endl;
    std::cout << "第" << ROUND << "轮解密用时：" << elapsed_seconds.count() << "s" << std::endl;
    total_time += elapsed_seconds.count();
    total_time_Sbox += elapsed_seconds_Sbox.count();
    total_time_XOR += elapsed_seconds_XOR.count();
#endif
    auto start3 = std::chrono::system_clock::now();
    auto end3 = std::chrono::system_clock::now();
    auto start4 = std::chrono::system_clock::now();
    auto end4 = std::chrono::system_clock::now();

    double total_time_Linear = 0;
    std::chrono::duration<double> elapsed_seconds_Linear;
    double total_time_LinearBootstrap = 0;
    std::chrono::duration<double> elapsed_seconds_LinearBootstrap;

    // 常规轮1~11解密
    std::vector<TRLWE_1> lut_result(16); //
#if 1
    for (int r = ROUND - 1; r > 0; r--) // 11 rounds
    {
        std::cout << "第" << r << "轮解密开始" << std::endl;
        // 记录开始时间
        start1 = std::chrono::system_clock::now();
        rk_slice = slice_2d(rk, 16 * r);
        // 使用切片进行 XOR_Two 操作
        // #pragma omp parallel for
        for (int j = 0; j < 16; j++)
        {
            XOR_Two(cipher[j], cipher[j], rk_slice[j], 8);
        }
        // 记录解密验证时间
        end1 = std::chrono::system_clock::now();
        // 解密验证
        decryptTLWE(cipher, *sk, decrypt);
        flag = check_decrypt(decrypt, Yux_state[3 * r - 1]);
        if (flag)
        {
            std::cout << "第" << r << "轮XOR解密成功!!!!!!" << std::endl;
        }
        else
        {
            std::cout << "第" << r << "轮XOR解密失败" << std::endl;
            // return 0;
        }
        // 记录解密验证时间
        start2 = std::chrono::system_clock::now();
        if (0)
        {
            for (int i = 0; i < 16; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    TFHEpp::SM4_CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTRGSW[i][j], cipher[i][j], ek);
                }
            }
            for (int i = 0; i < 16; i++)
            {
                CipherSubBytesMixedPacking(lut_result[i], Table5, bootedTRGSW[i]);
            }
            // SampleExtract level 1
            std::vector<std::vector<TLWE_1>> Sbox_value(16, std::vector<TLWE_1>(8));
            for (int i = 0; i < 16; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    TFHEpp::SampleExtractIndex<typename privksP::targetP>(Sbox_value[i][j], lut_result[i], j);
                }
            }
            // Key Switch to LWE B on level 0
            for (int i = 0; i < 16; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    // level 1 -> level 0
                    TFHEpp::IdentityKeySwitch<iksP>(cipher[i][j], Sbox_value[i][j], ek.getiksk<iksP>());
                }
            }
        }
        // 记录时间
        end2 = std::chrono::system_clock::now();
        // 进行线性层逆变换
        cipherinvLinearLayer(cipher, Table5, ek);
        // 记录解密验证时间
        start3 = std::chrono::system_clock::now();
        // 解密验证
        decryptTLWE(cipher, *sk, decrypt);
        flag = check_decrypt(decrypt, Yux_state[3 * r - 2]);
        if (flag)
        {
            std::cout << "第" << r << "轮线性层解密成功!!!!!!" << std::endl;
        }
        else
        {
            std::cout << "第" << r << "轮线性层解密失败" << std::endl;
            // return 0;
        }
        // 记录解密验证时间
        end3 = std::chrono::system_clock::now();
        // 进行s盒逆变换
        Sbox(cipher, RoundConstantTRLWE, Table1, Table2, Table3, Table4, Table5, ek);
        // 记录解密验证时间
        start4 = std::chrono::system_clock::now();
        // 调用解密函数
        decryptTLWE(cipher, *sk, decrypt);
        // 检查解密结果是否正确
        flag = check_decrypt(decrypt, Yux_state[3 * r - 3]);
        if (flag)
        {
            std::cout << "第" << r << "轮Sbox解密成功!!!" << std::endl;
        }
        else
        {
            std::cout << "第" << r << "轮Sbox解密失败" << std::endl;
            // return 0;
        }
        // 记录解密验证时间也是结束时间
        end4 = std::chrono::system_clock::now();

        // 计算XOR时间
        elapsed_seconds_XOR = end1 - start1;
        // 计算LinearBootstrap时间
        elapsed_seconds_LinearBootstrap = end2 - start2;
        // 计算Linear时间
        elapsed_seconds_Linear = start3 - end2;
        // 计算Sbox时间
        elapsed_seconds_Sbox = start4 - end3;
        // 计算解密时间，四个时间相加
        elapsed_seconds = elapsed_seconds_XOR + elapsed_seconds_LinearBootstrap + elapsed_seconds_Linear + elapsed_seconds_Sbox;
        std::cout << "第" << r << "轮解密完成" << std::endl;
        std::cout << "第" << r << "轮XOR用时：" << elapsed_seconds_XOR.count() << "s" << std::endl;
        std::cout << "第" << r << "轮LinearBootstrap用时：" << elapsed_seconds_LinearBootstrap.count() << "s" << std::endl;
        std::cout << "第" << r << "轮Linear用时：" << elapsed_seconds_Linear.count() << "s" << std::endl;
        std::cout << "第" << r << "轮Sbox用时：" << elapsed_seconds_Sbox.count() << "s" << std::endl;
        std::cout << "第" << r << "轮解密用时：" << elapsed_seconds.count() << "s" << std::endl;
        total_time += elapsed_seconds.count();
        total_time_XOR += elapsed_seconds_XOR.count();
        total_time_LinearBootstrap += elapsed_seconds_LinearBootstrap.count();
        total_time_Linear += elapsed_seconds_Linear.count();
        total_time_Sbox += elapsed_seconds_Sbox.count();
    }
#endif
/*
//输出decrypt
    for (int i = 0; i < 16; i++)
    {
        printf("%02x ", decrypt[i]);
    }
    std::cout << std::endl;

    // 对decrypt重新tlwe加密
    std::vector<std::vector<TLWE_0>> cipher2;
    cipher2.resize(16);
    for (int j = 0; j < 16; j++)
    {
        cipher2[j].resize(8);
    }
    for (int j = 0; j < 16; j++)
    {
        int bin_str[8];
        HexToBinStr(decrypt[j], bin_str, 8);
        for (int k = 0; k < 8; k++)
        {
            cipher2[j][k] = TFHEpp::tlweSymIntEncrypt<typename bkP::domainP>((typename bkP::domainP::T)bin_str[k], bkP::domainP::alpha,
                                                                             sk->key.get<typename bkP::domainP>());
        }
    }
    // 对cipher2进行XOR
    rk_slice = slice_2d(rk, 0);

#pragma omp parallel for
    for (int i = 0; i < 16; i++)
    {
        XOR_Two(cipher2[i], cipher2[i], rk_slice[i], 8);
    }
    // 解密验证
    decryptTLWE(cipher2, *sk, decrypt);
    // 输出解密结果
    for (int j = 0; j < 16; j++)
    {
        printf("%02x ", decrypt[j]);
    }
    std::cout << std::endl;
*/
// 白化轮即第0轮解密
#if 1
    std::cout << "第0轮解密开始" << std::endl;
    // 记录开始时间
    start = std::chrono::system_clock::now();
    // 使用切片进行 XOR_Two 操作
    rk_slice = slice_2d(rk, 0);
    // #pragma omp parallel for
    for (int i = 0; i < 16; i++)
    {
        XOR_Two(cipher[i], cipher[i], rk_slice[i], 8);
    }
    // 记录解密验证时间
    start1 = std::chrono::system_clock::now();
    // 调用解密函数
    decryptTLWE(cipher, *sk, decrypt);
    // 检查解密结果是否正确
    flag = check_decrypt(decrypt, plain);
    if (flag)
    {
        std::cout << "第0轮XOR解密成功!!!" << std::endl;
    }
    else
    {
        std::cout << "第0轮XOR解密失败" << std::endl;
    }
    // 记录解密验证时间也是结束时间
    end1 = std::chrono::system_clock::now();
    elapsed_seconds = end1 - start - (end1 - start1);
    std::cout << "第0轮解密完成" << std::endl;
    std::cout << "第0轮解密(XOR)用时：" << elapsed_seconds.count() << "s" << std::endl;
    total_time += elapsed_seconds.count();
    total_time_XOR += elapsed_seconds.count();
    // 输出各项总时间和百分比
    std::cout << "XOR总时间：" << total_time_XOR << "s, " << std::fixed << std::setprecision(6) << total_time_XOR / total_time * 100 << "%" << std::endl;
    std::cout << "LinearBootstrap总时间：" << total_time_LinearBootstrap << "s, " << std::fixed << std::setprecision(6) << total_time_LinearBootstrap / total_time * 100 << "%" << std::endl;
    std::cout << "Linear总时间：" << total_time_Linear << "s, " << std::fixed << std::setprecision(6) << total_time_Linear / total_time * 100 << "%" << std::endl;
    std::cout << "Sbox总时间：" << total_time_Sbox << "s, " << std::fixed << std::setprecision(6) << total_time_Sbox / total_time * 100 << "%" << std::endl;
    // 输出总时间
    std::cout << "解密总时间：" << total_time << "s" << std::endl;
#endif
    return 0;
}
// 总时间 = XOR时间 + Sbox时间 + Linear自举时间 + Linear，其中XOR时间和Linear时间忽略不计。
// Sbox时间占80%，Linear自举时间占20%，Sbox相当于4*12次自举，Linear相当于12次自举，总自举次数为60次
// 优化方向：1.减少Sbox时间，2.减少Linear自举时间
// 设想：Sbox减少为12次自举，Linear前面不需要自举。总自举次数减少到12次
// 总时间减少到原来的1/5，如果原时间为160s，则减少到32s