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

// AES 加密算法的实现，查大表法Tbox
#if 1
// AES密钥扩展
extern long AESKeyExpansion(unsigned char roundKeySchedule[],
                            unsigned char key[], int keyBits);
// AES S 盒（TBox01），和优化后的查表（TBox02 和 TBox03）
const int TBox01[256] = {0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76, 0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0, 0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15, 0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75, 0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84, 0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF, 0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8, 0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2, 0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73, 0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB, 0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C, 0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79, 0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08, 0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A, 0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E, 0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF, 0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16};
const int TBox02[256] = {0xC6, 0xF8, 0xEE, 0xF6, 0xFF, 0xD6, 0xDE, 0x91, 0x60, 0x02, 0xCE, 0x56, 0xE7, 0xB5, 0x4D, 0xEC, 0x8F, 0x1F, 0x89, 0xFA, 0xEF, 0xB2, 0x8E, 0xFB, 0x41, 0xB3, 0x5F, 0x45, 0x23, 0x53, 0xE4, 0x9B, 0x75, 0xE1, 0x3D, 0x4C, 0x6C, 0x7E, 0xF5, 0x83, 0x68, 0x51, 0xD1, 0xF9, 0xE2, 0xAB, 0x62, 0x2A, 0x08, 0x95, 0x46, 0x9D, 0x30, 0x37, 0x0A, 0x2F, 0x0E, 0x24, 0x1B, 0xDF, 0xCD, 0x4E, 0x7F, 0xEA, 0x12, 0x1D, 0x58, 0x34, 0x36, 0xDC, 0xB4, 0x5B, 0xA4, 0x76, 0xB7, 0x7D, 0x52, 0xDD, 0x5E, 0x13, 0xA6, 0xB9, 0x00, 0xC1, 0x40, 0xE3, 0x79, 0xB6, 0xD4, 0x8D, 0x67, 0x72, 0x94, 0x98, 0xB0, 0x85, 0xBB, 0xC5, 0x4F, 0xED, 0x86, 0x9A, 0x66, 0x11, 0x8A, 0xE9, 0x04, 0xFE, 0xA0, 0x78, 0x25, 0x4B, 0xA2, 0x5D, 0x80, 0x05, 0x3F, 0x21, 0x70, 0xF1, 0x63, 0x77, 0xAF, 0x42, 0x20, 0xE5, 0xFD, 0xBF, 0x81, 0x18, 0x26, 0xC3, 0xBE, 0x35, 0x88, 0x2E, 0x93, 0x55, 0xFC, 0x7A, 0xC8, 0xBA, 0x32, 0xE6, 0xC0, 0x19, 0x9E, 0xA3, 0x44, 0x54, 0x3B, 0x0B, 0x8C, 0xC7, 0x6B, 0x28, 0xA7, 0xBC, 0x16, 0xAD, 0xDB, 0x64, 0x74, 0x14, 0x92, 0x0C, 0x48, 0xB8, 0x9F, 0xBD, 0x43, 0xC4, 0x39, 0x31, 0xD3, 0xF2, 0xD5, 0x8B, 0x6E, 0xDA, 0x01, 0xB1, 0x9C, 0x49, 0xD8, 0xAC, 0xF3, 0xCF, 0xCA, 0xF4, 0x47, 0x10, 0x6F, 0xF0, 0x4A, 0x5C, 0x38, 0x57, 0x73, 0x97, 0xCB, 0xA1, 0xE8, 0x3E, 0x96, 0x61, 0x0D, 0x0F, 0xE0, 0x7C, 0x71, 0xCC, 0x90, 0x06, 0xF7, 0x1C, 0xC2, 0x6A, 0xAE, 0x69, 0x17, 0x99, 0x3A, 0x27, 0xD9, 0xEB, 0x2B, 0x22, 0xD2, 0xA9, 0x07, 0x33, 0x2D, 0x3C, 0x15, 0xC9, 0x87, 0xAA, 0x50, 0xA5, 0x03, 0x59, 0x09, 0x1A, 0x65, 0xD7, 0x84, 0xD0, 0x82, 0x29, 0x5A, 0x1E, 0x7B, 0xA8, 0x6D, 0x2C};
const int TBox03[256] = {0xA5, 0x84, 0x99, 0x8D, 0x0D, 0xBD, 0xB1, 0x54, 0x50, 0x03, 0xA9, 0x7D, 0x19, 0x62, 0xE6, 0x9A, 0x45, 0x9D, 0x40, 0x87, 0x15, 0xEB, 0xC9, 0x0B, 0xEC, 0x67, 0xFD, 0xEA, 0xBF, 0xF7, 0x96, 0x5B, 0xC2, 0x1C, 0xAE, 0x6A, 0x5A, 0x41, 0x02, 0x4F, 0x5C, 0xF4, 0x34, 0x08, 0x93, 0x73, 0x53, 0x3F, 0x0C, 0x52, 0x65, 0x5E, 0x28, 0xA1, 0x0F, 0xB5, 0x09, 0x36, 0x9B, 0x3D, 0x26, 0x69, 0xCD, 0x9F, 0x1B, 0x9E, 0x74, 0x2E, 0x2D, 0xB2, 0xEE, 0xFB, 0xF6, 0x4D, 0x61, 0xCE, 0x7B, 0x3E, 0x71, 0x97, 0xF5, 0x68, 0x00, 0x2C, 0x60, 0x1F, 0xC8, 0xED, 0xBE, 0x46, 0xD9, 0x4B, 0xDE, 0xD4, 0xE8, 0x4A, 0x6B, 0x2A, 0xE5, 0x16, 0xC5, 0xD7, 0x55, 0x94, 0xCF, 0x10, 0x06, 0x81, 0xF0, 0x44, 0xBA, 0xE3, 0xF3, 0xFE, 0xC0, 0x8A, 0xAD, 0xBC, 0x48, 0x04, 0xDF, 0xC1, 0x75, 0x63, 0x30, 0x1A, 0x0E, 0x6D, 0x4C, 0x14, 0x35, 0x2F, 0xE1, 0xA2, 0xCC, 0x39, 0x57, 0xF2, 0x82, 0x47, 0xAC, 0xE7, 0x2B, 0x95, 0xA0, 0x98, 0xD1, 0x7F, 0x66, 0x7E, 0xAB, 0x83, 0xCA, 0x29, 0xD3, 0x3C, 0x79, 0xE2, 0x1D, 0x76, 0x3B, 0x56, 0x4E, 0x1E, 0xDB, 0x0A, 0x6C, 0xE4, 0x5D, 0x6E, 0xEF, 0xA6, 0xA8, 0xA4, 0x37, 0x8B, 0x32, 0x43, 0x59, 0xB7, 0x8C, 0x64, 0xD2, 0xE0, 0xB4, 0xFA, 0x07, 0x25, 0xAF, 0x8E, 0xE9, 0x18, 0xD5, 0x88, 0x6F, 0x72, 0x24, 0xF1, 0xC7, 0x51, 0x23, 0x7C, 0x9C, 0x21, 0xDD, 0xDC, 0x86, 0x85, 0x90, 0x42, 0xC4, 0xAA, 0xD8, 0x05, 0x01, 0x12, 0xA3, 0x5F, 0xF9, 0xD0, 0x91, 0x58, 0x27, 0xB9, 0x38, 0x13, 0xB3, 0x33, 0xBB, 0x70, 0x89, 0xA7, 0xB6, 0x22, 0x92, 0x20, 0x49, 0xFF, 0x78, 0x7A, 0x8F, 0xF8, 0x80, 0x17, 0xDA, 0x31, 0xC6, 0xB8, 0xC3, 0xB0, 0x77, 0x11, 0xCB, 0xFC, 0xD6, 0x3A};

using TBoxArray = std::array<int, 4>;

// TBox 函数：查表优化的轮函数
TBoxArray TBox_0(unsigned char x)
{
    return {TBox02[x], TBox01[x], TBox01[x], TBox03[x]};
}

TBoxArray TBox_1(unsigned char x)
{
    return {TBox03[x], TBox02[x], TBox01[x], TBox01[x]};
}

TBoxArray TBox_2(unsigned char x)
{
    return {TBox01[x], TBox03[x], TBox02[x], TBox01[x]};
}

TBoxArray TBox_3(unsigned char x)
{
    return {TBox01[x], TBox01[x], TBox03[x], TBox02[x]};
}

// XOR 五个 TBoxArray，优化的查表法
TBoxArray XOR_five(const TBoxArray &a, const TBoxArray &b, const TBoxArray &c, const TBoxArray &d, const TBoxArray &e)
{
    TBoxArray result;
    for (int i = 0; i < 4; i++)
    {
        result[i] = a[i] ^ b[i] ^ c[i] ^ d[i] ^ e[i];
    }
    return result;
}

// 白化轮：初始密钥加操作
// 白化轮函数
void AES_WhiteRound(const std::array<unsigned char, 16> &plain, std::array<unsigned char, 16> &state, const std::array<unsigned char, 16> &round_key)
{
    for (int i = 0; i < 16; ++i)
    {
        state[i] = plain[i] ^ round_key[i];
    }
}

// AES 常规轮，使用查表优化法
void AES_round(std::array<unsigned char, 16> &y, const std::array<unsigned char, 16> &key)
{
    std::array<TBoxArray, 4> keys = {
        TBoxArray{key[0], key[1], key[2], key[3]},
        TBoxArray{key[4], key[5], key[6], key[7]},
        TBoxArray{key[8], key[9], key[10], key[11]},
        TBoxArray{key[12], key[13], key[14], key[15]}};

    std::array<TBoxArray, 4> temps;
    temps[0] = XOR_five(TBox_0(y[0]), TBox_1(y[5]), TBox_2(y[10]), TBox_3(y[15]), keys[0]);
    temps[1] = XOR_five(TBox_0(y[4]), TBox_1(y[9]), TBox_2(y[14]), TBox_3(y[3]), keys[1]);
    temps[2] = XOR_five(TBox_0(y[8]), TBox_1(y[13]), TBox_2(y[2]), TBox_3(y[7]), keys[2]);
    temps[3] = XOR_five(TBox_0(y[12]), TBox_1(y[1]), TBox_2(y[6]), TBox_3(y[11]), keys[3]);

    for (int i = 0; i < 4; ++i)
    {
        y[i] = temps[0][i];
        y[i + 4] = temps[1][i];
        y[i + 8] = temps[2][i];
        y[i + 12] = temps[3][i];
    }
}

// 最后一轮：SubBytes、ShiftRows 和 AddRoundKey
void AES_LastRound(std::array<unsigned char, 16> &y, const std::array<unsigned char, 16> &key)
{
    // SubBytes：查表法替换
    for (int i = 0; i < 16; i++)
    {
        y[i] = TBox01[y[i]];
    }

    // ShiftRows + AddRoundKey
    std::array<unsigned char, 16> temp = y;

    // 第 0 行左移 0 位（保持原样）
    y[0] = y[0] ^ key[0];
    y[4] = y[4] ^ key[4];
    y[8] = y[8] ^ key[8];
    y[12] = y[12] ^ key[12];

    // 第 1 行左移 1 位
    y[1] = temp[5] ^ key[1];
    y[5] = temp[9] ^ key[5];
    y[9] = temp[13] ^ key[9];
    y[13] = temp[1] ^ key[13];

    // 第 2 行左移 2 位
    y[2] = temp[10] ^ key[2];
    y[6] = temp[14] ^ key[6];
    y[10] = temp[2] ^ key[10];
    y[14] = temp[6] ^ key[14];

    // 第 3 行左移 3 位
    y[3] = temp[15] ^ key[3];
    y[7] = temp[3] ^ key[7];
    y[11] = temp[7] ^ key[11];
    y[15] = temp[11] ^ key[15];
}
#endif

// key-switching parameters
using iksP = TFHEpp::lvl10param;
// bootstrapping parameters
using bkP = TFHEpp::lvl02param;
// private key-switching parameters
using privksP = TFHEpp::lvl21param;

using TLWE_0 = TFHEpp::TLWE<typename bkP::domainP>;
using TLWE_1 = TFHEpp::TLWE<typename privksP::targetP>; // level 1

using TRLWE_1 = TFHEpp::TRLWE<typename privksP::targetP>; // level 1

const uint32_t byte_mul2[8] = {0, 0, 0, 0, 0, 0, 0, 0};

// 计算 NX2
const privksP::targetP::T NX2 = 2 * privksP::targetP::n;

// 定义全局常量 bara
const privksP::targetP::T bara[7] = {
    NX2 - 8 * (1 << 0),
    NX2 - 8 * (1 << 1),
    NX2 - 8 * (1 << 2),
    NX2 - 8 * (1 << 3),
    NX2 - 8 * (1 << 4),
    NX2 - 8 * (1 << 5),
    NX2 - 8 * (1 << 6)};

// AES逆S盒
static const unsigned char InvSboxTable[16][16] = {
    //  0     1     2     3     4     5     6     7     8     9     A     B     C     D     E     F
    0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb, // 0
    0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb, // 1
    0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e, // 2
    0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25, // 3
    0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92, // 4
    0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84, // 5
    0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06, // 6
    0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b, // 7
    0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73, // 8
    0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e, // 9
    0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b, // A
    0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4, // B
    0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f, // C
    0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef, // D
    0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61, // E
    0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d  // F
};

const uint32_t Rcon[10] = {0x01000000, 0x02000000, 0x04000000, 0x08000000, 0x10000000, 0x20000000, 0x40000000, 0x80000000, 0x1B000000, 0x36000000};

// 16进制转2进制，8bits，
// binstr[0]是最低位，binstr[7]是最高位，
// 有一个性质：00~7F的最高位是0，80~FF的最高位是1
// 所以sboxtable前128个元素的最高位是0，后128个元素的最高位是1
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
// 四个一维数组a,b,c,d相加，结果保存在数组result中，数组大小为8*(n+1)
template <class P>
void XOR_Four(P &result, P &a, P &b, P &c, P &d)
{
    XOR_Two<P>(result, a, b, 8);
    XOR_Two<P>(result, result, c, 8);
    XOR_Two<P>(result, result, d, 8);
}
// 生成逆S盒的tlwe加密
void MakeInvSBoxTable(std::vector<TRLWE_1> &Table, const TFHEpp::Key<privksP::targetP> &key)
{
    // dec-> bit
    int Sbox_binary[256][8];
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            int bin_str[8];
            HexToBinStr(InvSboxTable[i][j], bin_str, 8);
            for (int k = 0; k < 8; k++)
            {
                Sbox_binary[i * 16 + j][k] = bin_str[k];
                // cout<< Sbox_binary[i*16+j][k] << " ";
            }
            // cout << endl;
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
void CipherSubBytesMixedPacking(TRLWE_1 &result, std::vector<TRLWE_1> &Table, std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> &select)
{
    TFHEpp::CMUXFFT<typename privksP::targetP>(result, select[7], Table[1], Table[0]);
    TFHEpp::BlindRotate_LUT<privksP>(result, bara, select, 7); //, resultOfCMUX);
}
// 逆行移位
void CipherShiftRowsInverse(std::vector<std::vector<TLWE_0>> &cipher)
{
    // 临时矩阵存储当前的 cipher 状态
    std::vector<std::vector<TLWE_0>> temp = cipher;

    for (int j = 0; j < 8; j++)
    {
        // 第 0 行保持不变
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[0][j], temp[0][j]);
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[4][j], temp[4][j]);
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[8][j], temp[8][j]);
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[12][j], temp[12][j]);

        // 第 1 行向右移 1 个字节
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[1][j], temp[13][j]);
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[5][j], temp[1][j]);
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[9][j], temp[5][j]);
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[13][j], temp[9][j]);

        // 第 2 行向右移 2 个字节
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[2][j], temp[10][j]);
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[6][j], temp[14][j]);
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[10][j], temp[2][j]);
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[14][j], temp[6][j]);

        // 第 3 行向右移 3 个字节
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[3][j], temp[7][j]);
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[7][j], temp[11][j]);
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[11][j], temp[15][j]);
        TFHEpp::HomCOPY<typename bkP::domainP>(cipher[15][j], temp[3][j]);
    }
}
// Helper functions to perform multiplication by 02 in GF(2^8)
void CipherMul2(std::vector<TLWE_0> &byte, std::vector<TLWE_0> &consByte)
{
    // byte((b7 b6 b5 b4 b3 b2 b1)*x = b6 b5 b4 b3 b2 b1 b0 b7 + 000b7 b70 b7 0
    // consByte = 0000 0000
    std::vector<TLWE_0> temp(8);

    for (int i = 1; i < 8; i++)
    {
        TFHEpp::HomCOPY<typename bkP::domainP>(temp[i], byte[i - 1]);
    }
    TFHEpp::HomCOPY<typename bkP::domainP>(temp[0], byte[7]);

    for (int i = 0; i < 8; i++)
    {
        if (i == 1 || i == 3 || i == 4)
        {
            // 000b7 b70 b7 0
            TFHEpp::HomCOPY<typename bkP::domainP>(consByte[i], byte[7]);
        }
    }
    XOR_Two(byte, temp, consByte, 8);
}
// 逆列混淆
void CipherInvMixColumns(std::vector<std::vector<TLWE_0>> &cipher, std::vector<TLWE_0> &consByte)
{
    /*
    [0e  0b  0d  09] [ s00  s01  s02  s03]
    |09  0e  0b  0d| | s10  s11  s12  s13|
    |0d  09  0e  0b| | s20  s21  s22  s23|
    [0b  0d  09  0e] [ s30  s31  s32  s33]
    e=14=x^3+x^2+x; b=11=x^3+x+1; d=13=x^3+x^2+1; 9=x^3+1
    Tmp = s[0][i] ^ s[1][i] ^ s[2][i] ^ s[3][i];
    y[i][0]= Tmp ^ xtime(s[0][i] ^ s[1][i]) ^ x^2time(s[0][i]^s[2][i]) ^ x^3time(Tmp) ^ s[0][i];
    y[i][1]= Tmp ^ xtime(s[1][i] ^ s[2][i]) ^ x^2time(s[1][i]^s[3][i]) ^ x^3time(Tmp) ^ s[1][i];
    y[i][2]= Tmp ^ xtime(s[2][i] ^ s[3][i]) ^ x^2time(s[2][i]^s[0][i]) ^ x^3time(Tmp) ^ s[2][i];
    y[i][3]= Tmp ^ xtime(s[3][i] ^ s[0][i]) ^ x^2time(s[3][i]^s[1][i]) ^ x^3time(Tmp) ^ s[3][i];
    */
    // 提前初始化临时向量以避免在循环中重复创建
    std::vector<TLWE_0> Tmp(8);
    std::vector<TLWE_0> Tmp3(8);
    std::vector<TLWE_0> Tm(8);
    std::vector<TLWE_0> Tm1(8);
    std::vector<TLWE_0> Tm2(8);
    std::vector<std::vector<TLWE_0>> T(4);
    for (int i = 0; i < 4; i++)
    {
        // 将cipher[4*i],cipher[4*i+1],cipher[4*i+2],cipher[4*i+3]拷贝到T中
        for (int j = 0; j < 4; j++)
        {
            T[j].resize(8);
            for (int k = 0; k < 8; k++)
            {
                TFHEpp::HomCOPY<typename bkP::domainP>(T[j][k], cipher[4 * i + j][k]);
            }
        }
        // Tmp = cipher[0][i] ^ cipher[1][i] ^ cipher[2][i] ^ cipher[3][i];
        XOR_Four(Tmp, T[0], T[1], T[2], T[3]);

        // Tmp3 = x^3 * Tmp
        Tmp3 = Tmp;
        for (int j = 0; j < 3; j++)
        {
            CipherMul2(Tmp3, consByte); // 执行三次 CipherMul2 来实现 x^3 乘法
        }

        for (int row = 0; row < 4; row++)
        {
            // 将当前行的数据拷贝到 Tm
            for (int j = 0; j < 8; j++)
            {
                TFHEpp::HomCOPY<typename bkP::domainP>(Tm[j], T[row][j]);
            }

            // 计算 Tm1 = Tm ^ cipher[row+1][i]
            XOR_Two(Tm1, Tm, T[(row + 1) % 4], 8);
            CipherMul2(Tm1, consByte); // Tm1 = x * Tm1

            // 计算 Tm2 = Tm ^ cipher[row+2][i]
            XOR_Two(Tm2, Tm, T[(row + 2) % 4], 8);
            CipherMul2(Tm2, consByte); // Tm2 = x * Tm2
            CipherMul2(Tm2, consByte); // 再次乘以 x 得到 x^2 * Tm2

            // 最终计算 cipher[4*i+row] = cipher[4*i+row] ^ Tmp ^ Tm1 ^ Tm2 ^ Tmp3;
            XOR_Four(cipher[4 * i + row], T[row], Tmp, Tm1, Tm2);
            XOR_Two(cipher[4 * i + row], cipher[4 * i + row], Tmp3, 8);
        }
    }
}
// tlwe初始密钥下的拓展
template <class bkp>
void RotWord_SubBytes(std::vector<std::vector<TLWE_0>> &cipher, std::vector<TRLWE_1> &Table, const TFHEpp::EvalKey &ek, int N)
{
    // 循环左移8位
    std::vector<std::vector<TLWE_0>> temp = cipher;
    // 行移位，让第0行变为第3行，第1行变为第0行，第2行变为第1行，第3行变为第2行
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            TFHEpp::HomCOPY<typename bkP::domainP>(cipher[i][j], temp[(i + 1) % 4][j]);
        }
    }
    // N=4
    std::vector<std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>>> bootedTRGSW;
    bootedTRGSW.resize(N);
    for (int i = 0; i < N; i++)
    {
        bootedTRGSW[i].resize(8);
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 8; j++)

            TFHEpp::SM4_CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTRGSW[i][j],
                                                                    cipher[i][j], ek);
    }

    std::vector<TRLWE_1> lut_result(N); //
    for (int i = 0; i < N; i++)
    {
        CipherSubBytesMixedPacking(lut_result[i], Table, bootedTRGSW[i]);
    }
    std::vector<std::vector<TLWE_1>> Sbox_value;
    Sbox_value.resize(N);
    for (int i = 0; i < Sbox_value.size(); i++)
    {
        Sbox_value[i].resize(8);
    }
    // SampleExtract level 1
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            TFHEpp::SampleExtractIndex<typename privksP::targetP>(Sbox_value[i][j], lut_result[i], j);
        }
    }
    std::vector<std::vector<TLWE_0>> B;
    B.resize(N);
    for (int i = 0; i < N; i++)
    {
        B[i].resize(8);
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            // level 1 -> level 0
            TFHEpp::IdentityKeySwitch<iksP>(B[i][j], Sbox_value[i][j], ek.getiksk<iksP>());
        }
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            TFHEpp::HomCOPY<typename bkP::domainP>(cipher[i][j], B[i][j]);
        }
    }
}
// 密钥扩展的tlwe加密
template <typename bkP>
void encryptRoundKeys(const unsigned char *RoundKey, std::vector<std::vector<TFHEpp::TLWE<typename bkP::domainP>>> &rk, const TFHEpp::SecretKey &sk)
{
    for (int i = 0; i < 176; i++)
    {
        int bin_str[8];
        rk[i].resize(8);
        HexToBinStr(RoundKey[i], bin_str, 8);

        for (int j = 0; j < 8; j++)
        {
            // encrypt TLWE in level 0
            rk[i][j] = TFHEpp::tlweSymIntEncrypt<typename bkP::domainP>(
                static_cast<typename bkP::domainP::T>(bin_str[j]),
                bkP::domainP::alpha,
                sk.key.get<typename bkP::domainP>());
        }
    }
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
// 检查解密结果是否正确
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
// 函数用于从二维数组中快速切片
std::vector<std::vector<TLWE_0>> slice_2d(const std::vector<std::vector<TLWE_0>> &matrix, int start_row)
{
    std::vector<std::vector<TLWE_0>> slice; // 创建切片的目标容器
    slice.resize(16);
    for (int i = 0; i < 16; i++)
    {
        slice[i].resize(8);
    }
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
// 对逆s盒进行TRLWE加密
#if 1
    std::vector<TRLWE_1> Table(2);
    MakeInvSBoxTable(Table, sk->key.get<privksP::targetP>()); // 利用level1的trlwe key
#endif
// 输出S盒Table,4行，1024列，每个元素是个整数
#if 0
    for (const auto &trlwe : Table)
    {
        for (const auto &poly : trlwe)
        {
            for (const auto &coeff : poly)
            {
                std::cout << coeff << " ";
            }
            std::cout << std::endl;
        }
    }
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
// 定义AES明文和初始密钥
#if 1
    // AES明文
    std::array<unsigned char, 16> plain = {0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34};
    // AES初始密钥
    std::array<unsigned char, 16> aeskey = {0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};
#endif
// AES轮密钥生成
#if 1
    unsigned char RoundKey[176];
    long nRoundKeys = AESKeyExpansion(RoundKey, aeskey.data(), 128);
    std::cout << "rounds of key: " << nRoundKeys << endl;
    // 将RoundKey转为二维数组，11行，16列
    std::array<std::array<unsigned char, 16>, 11> RoundKey2D;
    for (int i = 0; i < 11; ++i)
    {
        for (int j = 0; j < 16; ++j)
        {
            RoundKey2D[i][j] = RoundKey[i * 16 + j];
        }
    }
    int n_round = nRoundKeys - 1; // 10，加密轮数=轮密钥数-1
#endif
// AES加密明文
#if 1
    // 用轮密钥加密明文
    std::array<unsigned char, 16> state;
    // 定义一个二维数组，来存储AES每一轮加密后的state,10行16列。不是tlwe密文
    std::array<std::array<unsigned char, 16>, 10> AES_state;
    std::cout << "AES加密开始" << endl;
    // 白化轮
    AES_WhiteRound(plain, state, RoundKey2D[0]);
    // 将state赋值给AES_state
    for (int i = 0; i < 16; i++)
    {
        AES_state[0][i] = state[i];
    }
    // 9 轮常规 AES 轮
    for (int i = 0; i < 9; i++)
    {
        AES_round(state, RoundKey2D[i + 1]); // 每一轮都用相同的 key（在实际中应该扩展轮密钥）
                                             // 将state赋值给AES_state
        for (int j = 0; j < 16; j++)
        {
            AES_state[i + 1][j] = state[j];
        }
    }
    // 最后一轮
    AES_LastRound(state, RoundKey2D[10]);
    std::cout << "AES加密完成" << endl;
    // 输出AES加密后的密文
    std::cout << "AES加密后的密文：";
    for (int i = 0; i < 16; i++)
    {
        std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)state[i] << " ";
    }
    std::cout << std::endl;
#endif
// 对密文进行tlwe加密
#if 1
    std::cout << "AES密文tlwe加密开始" << endl;
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
        // #pragma omp parallel for num_threads(8)
        for (int j = 0; j < 8; j++)
        {
            cipher[i][j] = TFHEpp::tlweSymIntEncrypt<typename bkP::domainP>((typename bkP::domainP::T)bin_str[j], bkP::domainP::alpha,
                                                                            sk->key.get<typename bkP::domainP>());
        }
    }
    std::cout << "AES密文tlwe加密完成" << endl;
#endif
// AES密文的tlwe加解密验证
// 对cipher进行解密，与state进行比较
#if 1
    std::array<unsigned char, 16> decrypt;
    std::cout << "AES密文tlwe解密开始" << endl;
    decryptTLWE(cipher, *sk, decrypt);
    std::cout << "AES密文tlwe解密完成" << endl;
    int flag = check_decrypt(decrypt, state);

    if (flag)
    {
        std::cout << "AES密文tlwe解密成功!!!" << endl;
    }
    else
    {
        std::cout << "AES密文tlwe解密失败" << endl;
        return 0;
    }
#endif
// 1、初始密钥tlwe加密，并进行tlwe下的拓展
#if 0
    std::vector<TLWE_0> AES_Key_tlwed;
    AES_Key_tlwed.resize(128);
    // tlwe加密AES初始密钥
    std::cout << "密钥扩展开始" << endl;
    for (int i = 0; i < 16; i++)
    {
        int bin_str[8];
        HexToBinStr(aeskey[i], bin_str, 8);
        for (int j = 0; j < 8; j++)
        {
            AES_Key_tlwed[8 * i + j] = TFHEpp::tlweSymIntEncrypt<typename bkP::domainP>(
                static_cast<typename bkP::domainP::T>(bin_str[j]),
                bkP::domainP::alpha,
                sk->key.get<typename bkP::domainP>());
        }
    }
    // 对AES_Key_tlwed进行AES密钥扩展
    std::vector<std::vector<TLWE_0>> AES_RoundKey_tlwed;
    AES_RoundKey_tlwed.resize(44);
    for (int i = 0; i < 44; i++)
    {
        AES_RoundKey_tlwed[i].resize(32);
    }

    // AES_RoundKey_tlwed[0] = AES_Key_tlwed;copy
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 32; j++)
        {
            AES_RoundKey_tlwed[i][j] = AES_Key_tlwed[i * 32 + j];
        }
    }

    for (int i = 4; i < 44; i++)
    {
        std::vector<std::vector<TLWE_0>> temp;
        temp.resize(4);
        for (int j = 0; j < 4; j++)
        {
            temp[j].resize(8);
        }
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 8; k++)
            {
                temp[j][k] = AES_RoundKey_tlwed[i - 1][j * 8 + k];
            }
        }
        if (i % 4 == 0)
        {
            // 将Rcon[i / 4 - 1]转为32位2进制再进行tlwe加密
            int bin_str[32];
            HexToBinStr(Rcon[i / 4 - 1], bin_str, 32);

            std::vector<std::vector<TLWE_0>> Rcon_tlwed;
            Rcon_tlwed.resize(4);
            for (int j = 0; j < 4; j++)
            {
                Rcon_tlwed[j].resize(8);
            }
            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < 8; k++)
                {
                    Rcon_tlwed[j][k] = TFHEpp::tlweSymIntEncrypt<typename bkP::domainP>(
                        static_cast<typename bkP::domainP::T>(bin_str[j * 8 + k]),
                        bkP::domainP::alpha,
                        sk->key.get<typename bkP::domainP>());
                }
            }
            // 将 temp 转换为 std::vector<std::vector<TLWE_0>>
            std::vector<std::vector<TLWE_0>> SR_temp;
            SR_temp.resize(4);
            for (int j = 0; j < 4; j++)
            {
                SR_temp[j].resize(8);
            }
            RotWord_SubBytes<bkP>(SR_temp, Table, ek, 4); // 显式指定模板参数 bkP，4=32/8
            for (int j = 0; j < 4; j++)
            {
                XOR_Two(temp[j], SR_temp[j], Rcon_tlwed[j], 8);
            }
        }

        std::vector<TLWE_0> AES_temp;
        AES_temp.resize(32);
        // 令temp转为AES_temp
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 8; k++)
            {
                // 用homecopy
                TFHEpp::HomCOPY<typename bkP::domainP>(AES_temp[j * 8 + k], temp[j][k]);
            }
        }
        XOR_Two(AES_RoundKey_tlwed[i], AES_RoundKey_tlwed[i - 4], AES_temp, 32);
    }
    // 将 rk 定义为 176 行 8 列的二维数组
    std::vector<std::vector<TLWE_0>> rk;
    rk.resize(176);
    for (int i = 0; i < 176; i++)
    {
        rk[i].resize(8);
    }

    for (int i = 0; i < 176; i++) // 176 行
    {
        for (int j = 0; j < 8; j++) // 8 列
        {
            int src_row = (i * 8 + j) / 32;                  // 将目标位置映射到源数组行索引
            int src_col = (i * 8 + j) % 32;                  // 将目标位置映射到源数组列索引
            rk[i][j] = AES_RoundKey_tlwed[src_row][src_col]; // 从源数组获取值
        }
    }
#endif
// 2、AES密钥拓展的tlwe加解密
// 这里与上面的区别在于上面是对初始密钥进行tlwe加密再进行扩展，这里是对扩展后的密钥进行tlwe加解密。
// 显然后面的实现更简单，但是实际应用中如果只发送了初始密钥的加密结果，那么需要对加密后的初始密钥进行扩展，这里的实现就是对加密后的初始密钥进行扩展，则只能用上面的方法
#if 1
    std::vector<std::vector<TLWE_0>> rk;
    rk.resize(176);
    encryptRoundKeys<bkP>(RoundKey, rk, *sk);
#endif
    std::cout << "AES密钥扩展的tlwe加密完成" << endl;
/* 前面两种方法均可以得到密钥拓展tlwe加密后的轮密钥rk[176][8]
   下面解密前面的轮密钥rk，得到decrypt_rk[176][8]，
   将其转为每行转为4为hex整数，再变为11行16列，最后与未加密的拓展密钥RoundKey2D进行逐元素比较
*/
#if 1
    std::array<std::array<int, 16>, 11> decrypt_rk;
    std::cout << "AESkey的tlwe解密开始" << std::endl;
    for (int i = 0; i < 11; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            int dec_hex = 0;
            int dec_bin[8];
            for (int k = 0; k < 8; k++)
            {
                // typename P::T a = TFHEpp::tlweSymIntDecrypt<typename bkP::domainP>();
                dec_bin[k] = TFHEpp::tlweSymIntDecrypt<typename bkP::domainP>(rk[16 * i + j][k], sk->key.get<typename bkP::domainP>());
                // bootsSymDecrypt(&rk[0][i][j], key);
            }
            // 将二进制转为十六进制
            BinStrToHex(dec_bin, dec_hex, 8);
            decrypt_rk[i][j] = dec_hex;
        }
    }
    std::cout << "AESkey的tlwe解密完成" << endl;
    // 将解密后的轮密钥与原始轮密钥进行比较
    // 按行输出RounKey2D
    // 输出加密前和解密后的轮密钥
    flag = 1;
    for (int i = 0; i < 11; i++)
    {
        for (int j = 0; j < 16; j++)
        {

            if (RoundKey2D[i][j] != decrypt_rk[i][j])
            {
                std::cout << "AESkey的tlwe解密失败" << std::endl;
                // 输出i,j
                std::cout << "i:" << i << " j:" << j << std::endl;
                flag = 0;
                return 0;
            }
        }
    }
    if (flag)
    {
        std::cout << "AESkey的tlwe解密成功!!!" << std::endl;
    }
#endif
// 前面已经得到了AES密文且用tlwe对其进行加密得到了cipher，现在使用tlwe下的AES密钥扩展对cipher进行AES解密，这样就间接得到了用tlwe直接加密明文结果
#if 1
    std::vector<std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>>> bootedTRGSW;
    bootedTRGSW.resize(16);
    for (int i = 0; i < 16; i++)
    {
        bootedTRGSW[i].resize(8);
    }
#endif
// 最后一轮第10轮解密
#if 1

    std::vector<std::vector<TLWE_0>> rk_slice;
    rk_slice.resize(16);
    for (int i = 0; i < 16; i++)
    {
        rk_slice[i].resize(8);
    }
    std::cout << "第10轮解密开始" << endl;
    rk_slice = slice_2d(rk, 160);
    // 使用切片进行 XOR_Two 操作
    for (int i = 0; i < 16; i++)
    {
        XOR_Two(cipher[i], cipher[i], rk_slice[i], 8);
    }
    // 让state和 RoundKey2D[10]进行异或操作，tlwe解密后的结果和AES_state[9]进行比较
    /*
    std::array<unsigned char, 16> t;
    for (int i = 0; i < 16; i++)
    {
        t[i] = RoundKey2D[10][i] ^ state[i];
    }
    // 解密第10轮XOR
    decryptTLWE(cipher, *sk, decrypt);
    flag = check_decrypt(decrypt, t);
    if (flag)
    {
        std::cout << "第10轮tlweXOR解密成功" << endl;
    }
    else
    {
        std::cout << "第10轮tlweXOR解密失败" << endl;
        return 0;
    }*/
    // CipherShiftRowsInverse 操作
    CipherShiftRowsInverse(cipher);
    // CipherSubBytes();
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 8; j++)
            TFHEpp::SM4_CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTRGSW[i][j],
                                                                    cipher[i][j], ek);
    }
    std::vector<TRLWE_1> lut_result(16); //
    for (int i = 0; i < 16; i++)
    {
        CipherSubBytesMixedPacking(lut_result[i], Table, bootedTRGSW[i]);
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
    std::cout << "第10轮解密完成" << endl;
    // 第9轮解密
    std::cout << "第9轮tlwe解密开始" << endl;
    // 调用解密函数
    decryptTLWE(cipher, *sk, decrypt);
    std::cout << "第9轮tlwe解密完成" << endl;
    // 检查解密结果是否正确
    flag = check_decrypt(decrypt, AES_state[9]);
    if (flag)
    {
        std::cout << "第9轮tlwe解密成功!!!" << endl;
    }
    else
    {
        std::cout << "第9轮tlwe解密失败" << endl;
        return 0;
    }
#endif
// 常规轮1~9解密
#if 1
    for (int i = n_round - 1; i > 0; i--) // 9 rounds
    {
        std::cout << "第" << i << "轮解密开始" << endl;
        rk_slice = slice_2d(rk, 16 * i);
        // 使用切片进行 XOR_Two 操作
        for (int j = 0; j < 16; j++)
        {
            XOR_Two(cipher[j], cipher[j], rk_slice[j], 8);
        }
        // 进行列混淆逆变换
        CipherInvMixColumns(cipher, consByte); // 使用正确定义的函数
        // CipherShiftRowsInverse 操作
        CipherShiftRowsInverse(cipher);
        // CipherSubBytes();
        for (int i = 0; i < 16; i++)
        {
            for (int j = 0; j < 8; j++)
                TFHEpp::SM4_CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTRGSW[i][j],
                                                                        cipher[i][j], ek);
        }
        std::vector<TRLWE_1> lut_result(16); //
        for (int i = 0; i < 16; i++)
        {
            CipherSubBytesMixedPacking(lut_result[i], Table, bootedTRGSW[i]);
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
        std::cout << "第" << i << "轮解密完成" << endl;
        std::cout << "第" << i - 1 << "轮tlwe解密开始" << endl;
        // 调用解密函数
        decryptTLWE(cipher, *sk, decrypt);
        std::cout << "第" << i - 1 << "轮tlwe解密完成" << endl;
        // 检查解密结果是否正确
        flag = check_decrypt(decrypt, AES_state[i - 1]);
        if (flag)
        {
            std::cout << "第" << i - 1 << "轮tlwe解密成功!!!" << endl;
        }
        else
        {
            std::cout << "第" << i - 1 << "轮tlwe解密失败" << endl;
            return 0;
        }
    }
#endif
// 白化轮即第0轮解密
#if 1
    std::cout << "第0轮解密开始" << endl;
    // 使用切片进行 XOR_Two 操作
    rk_slice = slice_2d(rk, 0);
    for (int i = 0; i < 16; i++)
    {
        XOR_Two(cipher[i], cipher[i], rk_slice[i], 8);
    }
    std::cout << "第0轮解密完成" << endl;
#endif
// 至此，我们得到了用tlwe加密的AES明文的结果，下面将其用tlwe私钥解密，检验解密结果是否正确
// tlwe明文解密
#if 1
    std::cout << "tlwe明文解密开始" << std::endl;
    // 调用解密函数
    decryptTLWE(cipher, *sk, decrypt);
    std::cout << "tlwe明文解密完成" << endl;
    // 检查解密结果是否正确
    flag = check_decrypt(decrypt, plain);
    if (flag)
    {
        std::cout << "tlwe明文解密成功!!!" << endl;
    }
    else
    {
        std::cout << "tlwe明文解密失败" << endl;
    }
#endif
    return 0;
}