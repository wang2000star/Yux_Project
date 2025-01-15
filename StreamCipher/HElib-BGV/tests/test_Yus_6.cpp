#include <iostream>
#include <filesystem>
#include <vector>
#include <array>
#include <atomic>
#include <cmath>
#include <chrono>
#include <fstream>
#include <memory>
#include <random>
#include <climits>
#include <omp.h>
#include <cassert>

#include <NTL/ZZX.h>
// #include <NTL/GF2X.h>
#include <helib/helib.h>
#include <helib/ArgMap.h>
#include <helib/DoubleCRT.h>

#include "helib/CModulus.h"
#include "helib/powerful.h"
#include <immintrin.h> // 包含 SIMD 指令集支持

#include "../Symmetric/Yus_p.hpp"
#include "../utils/tool.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

namespace fs = std::filesystem;
/************************************************************************
  long p;          // plaintext primeplain_mod;
  long m;          // m-th cyclotomic polynomial
  long r;          // Lifting [defualt = 1]
  long bits;       // bits in the ciphertext modulus chain
  long c;          // columns in the key-switching matrix [default=2]
  long d;          // Degree of the field extension [default=1]
  long k;          // Security parameter [default=80]
  long s;          // Minimum number of slots [default=0]
************************************************************************/
// p^d = 1 mod m,d=1,slots=\phi(m)/d=\phi(m);m=65536=2^16,\phi(m)=2^15=32768
// 更一般的，应该有d|ord_p(m)，slots=\phi(m)/ord_p(m)
//!!!!!!!!!!!!!!!!
constexpr unsigned BlockWords = 36;      // 分组密钥字长度=KeyWords
constexpr unsigned BlockPlainWords = 24; // 明文分组字长度
constexpr double TruncRate = BlockPlainWords / (double)BlockWords;
// ===============模式设置================
constexpr bool Rkflag = 1;        // true/1表示乘法，false/0表示加法，指的是随机向量和密钥间的操作
constexpr bool deflag = 0;        // true/1表示进行每一步解密验证，false/0表示不进行每一步解密验证
constexpr bool ompflag = 0;       // true/1表示使用OpenMP并行编码，false/0表示不使用OpenMP并行编码
constexpr bool symkeyflag = 0;    // true/1表示对称密钥同态解密验证加密，false/0表示不验证
constexpr bool KeyStreamflag = 0; // true/1表示密钥流同态解密验证，false/0表示不验证
constexpr bool plainflag = 0;     // true/1表示明文同态解密验证，false/0表示不验证
// 参数设置，paramMap[Nr-4][idx]
constexpr unsigned Nr = 6; // 轮数
constexpr long idx = 0;
constexpr unsigned Sbox_depth = 1 * Nr; // S盒深度
// 当c=2时，Qbits=1.5*bits,当c=3时，Qbits=1.5*bits - 100
// 硬编码参数值
constexpr tuple<long, long, long, long> paramMap[1][8] = {
    {// Nr = 4
     // {p, log2(m), bits, c}
     {65537, 16, 330, 3},  // 0 *
     {163841, 15, 240, 2}, // 1
     {65537, 14, 220, 2},  // 2 *
     {163841, 14, 230, 2}, // 3
     {65537, 15, 350, 2},  // 4
     {163841, 15, 350, 2}, // 5
     {65537, 15, 400, 2},  // 6
     {163841, 15, 400, 2}} // 7
};
// p=k*m+1
//  2^10=1024,2^11=2048,2^12=4096,2^13=8192,2^14=16384,2^15=32768,2^16=65536,
// 当电脑内存有限时，log2Para_m太大，会导致内存不足，出现terminate called recursively错误，从而终止程序.
// 此外，即使正常运行，由于内存一直临界，会导致程序运行速度变慢，时间测量不准确。

constexpr long log2Para_m = get<1>(paramMap[0][idx]) - 0;
constexpr long Para_p = get<0>(paramMap[0][idx]);    // plaintext prime
constexpr long Para_m = 1 << log2Para_m;             // cyclotomic polynomial
constexpr long phi_m = Para_m >> 1;                  // phi(m)
constexpr long Para_bits = get<2>(paramMap[0][idx]); // bits in the ciphertext modulus chain
constexpr long Para_c = get<3>(paramMap[0][idx]);    // columns in the key-switching matrix
constexpr long Para_r = 1;                           // Lifting [defualt = 1]
//!!!!!!!!!!!!!!!
constexpr long nslots = phi_m;              // 槽数
constexpr unsigned PlainBlock = nslots - 0; // 明文分组数,应该PlainBlock<=nslots
constexpr unsigned len3 = BlockWords / 3;
// 计算 log2 的 constexpr 函数
constexpr unsigned int log2_constexpr(unsigned long long n, unsigned int p = 0)
{
    return (n <= 1) ? p : log2_constexpr(n / 2, p + 1);
}
constexpr long PlainMod = Para_p;                               // 明文模数
constexpr unsigned Wordbits = log2_constexpr(PlainMod - 1) + 1; // 字节比特长度=ceil(log2(PlainMod-1))
constexpr unsigned randbits = Wordbits - 1;
constexpr unsigned BlockSize = Wordbits * BlockWords;    // 分组比特长度=BlockWords*Wordbits
constexpr unsigned NrBlockWords = BlockWords * (Nr + 1); // Nr轮分组密钥字节长度

constexpr long KeyStreamWords = BlockWords * PlainBlock;  // 密钥流字节长度
constexpr long PlainWords = BlockPlainWords * PlainBlock; // 明文字节长度

constexpr long Plainbits = Wordbits * PlainWords; // 明文比特长度

constexpr long max_prime_size = (1ULL << (Wordbits - 1)) - 1;

constexpr unsigned NonceSize = 32;                           // Nonce比特长度
constexpr long counter_begin = 0;                            // 计数器起始值
constexpr long counter_end = PlainBlock + counter_begin - 1; // 计数器结束值

YusP yusP(PlainMod); // 构建明文对称加密实例

// Linear transformation
void HE_M(vector<Ctxt> &eData)
{
    vector<Ctxt> temp = eData;
    Ctxt temp0_1 = temp[0];
    temp0_1 += temp[1];
    Ctxt temp0_2 = temp[0];
    temp0_2 += temp[2];
    Ctxt temp0_3 = temp[0];
    temp0_3 += temp[3];
    Ctxt temp1_2 = temp[1];
    temp1_2 += temp[2];
    Ctxt temp1_3 = temp[1];
    temp1_3 += temp[3];
    Ctxt temp2_3 = temp[2];
    temp2_3 += temp[3];
    Ctxt temp0_1_2 = temp0_1;
    temp0_1_2 += temp[2];
    Ctxt temp0_1_3 = temp0_1;
    temp0_1_3 += temp[3];
    Ctxt temp0_2_3 = temp0_2;
    temp0_2_3 += temp[3];
    Ctxt temp1_2_3 = temp1_2;
    temp1_2_3 += temp[3];
    Ctxt temp0_1_2_3 = temp0_1_2;
    temp0_1_2_3 += temp[3];
    Ctxt temp4_5 = temp[4];
    temp4_5 += temp[5];
    Ctxt temp4_6 = temp[4];
    temp4_6 += temp[6];
    Ctxt temp4_7 = temp[4];
    temp4_7 += temp[7];
    Ctxt temp5_6 = temp[5];
    temp5_6 += temp[6];
    Ctxt temp5_7 = temp[5];
    temp5_7 += temp[7];
    Ctxt temp6_7 = temp[6];
    temp6_7 += temp[7];
    Ctxt temp4_5_6 = temp4_5;
    temp4_5_6 += temp[6];
    Ctxt temp4_5_7 = temp4_5;
    temp4_5_7 += temp[7];
    Ctxt temp4_6_7 = temp4_6;
    temp4_6_7 += temp[7];
    Ctxt temp5_6_7 = temp5_6;
    temp5_6_7 += temp[7];
    Ctxt temp4_5_6_7 = temp4_5_6;
    temp4_5_6_7 += temp[7];
    Ctxt temp8_9 = temp[8];
    temp8_9 += temp[9];
    Ctxt temp8_10 = temp[8];
    temp8_10 += temp[10];
    Ctxt temp8_11 = temp[8];
    temp8_11 += temp[11];
    Ctxt temp9_10 = temp[9];
    temp9_10 += temp[10];
    Ctxt temp9_11 = temp[9];
    temp9_11 += temp[11];
    Ctxt temp10_11 = temp[10];
    temp10_11 += temp[11];
    Ctxt temp8_9_10 = temp8_9;
    temp8_9_10 += temp[10];
    Ctxt temp8_9_11 = temp8_9;
    temp8_9_11 += temp[11];
    Ctxt temp8_10_11 = temp8_10;
    temp8_10_11 += temp[11];
    Ctxt temp9_10_11 = temp9_10;
    temp9_10_11 += temp[11];
    Ctxt temp8_9_10_11 = temp8_9_10;
    temp8_9_10_11 += temp[11];
    Ctxt temp12_13 = temp[12];
    temp12_13 += temp[13];
    Ctxt temp12_14 = temp[12];
    temp12_14 += temp[14];
    Ctxt temp12_15 = temp[12];
    temp12_15 += temp[15];
    Ctxt temp13_14 = temp[13];
    temp13_14 += temp[14];
    Ctxt temp13_15 = temp[13];
    temp13_15 += temp[15];
    Ctxt temp14_15 = temp[14];
    temp14_15 += temp[15];
    Ctxt temp12_13_14 = temp12_13;
    temp12_13_14 += temp[14];
    Ctxt temp12_13_15 = temp12_13;
    temp12_13_15 += temp[15];
    Ctxt temp12_14_15 = temp12_14;
    temp12_14_15 += temp[15];
    Ctxt temp13_14_15 = temp13_14;
    temp13_14_15 += temp[15];
    Ctxt temp12_13_14_15 = temp12_13_14;
    temp12_13_14_15 += temp[15];
    Ctxt temp16_17 = temp[16];
    temp16_17 += temp[17];
    Ctxt temp16_18 = temp[16];
    temp16_18 += temp[18];
    Ctxt temp16_19 = temp[16];
    temp16_19 += temp[19];
    Ctxt temp17_18 = temp[17];
    temp17_18 += temp[18];
    Ctxt temp17_19 = temp[17];
    temp17_19 += temp[19];
    Ctxt temp18_19 = temp[18];
    temp18_19 += temp[19];
    Ctxt temp16_17_18 = temp16_17;
    temp16_17_18 += temp[18];
    Ctxt temp16_17_19 = temp16_17;
    temp16_17_19 += temp[19];
    Ctxt temp16_18_19 = temp16_18;
    temp16_18_19 += temp[19];
    Ctxt temp17_18_19 = temp17_18;
    temp17_18_19 += temp[19];
    Ctxt temp16_17_18_19 = temp16_17_18;
    temp16_17_18_19 += temp[19];
    Ctxt temp20_21 = temp[20];
    temp20_21 += temp[21];
    Ctxt temp20_22 = temp[20];
    temp20_22 += temp[22];
    Ctxt temp20_23 = temp[20];
    temp20_23 += temp[23];
    Ctxt temp21_22 = temp[21];
    temp21_22 += temp[22];
    Ctxt temp21_23 = temp[21];
    temp21_23 += temp[23];
    Ctxt temp22_23 = temp[22];
    temp22_23 += temp[23];
    Ctxt temp20_21_22 = temp20_21;
    temp20_21_22 += temp[22];
    Ctxt temp20_21_23 = temp20_21;
    temp20_21_23 += temp[23];
    Ctxt temp20_22_23 = temp20_22;
    temp20_22_23 += temp[23];
    Ctxt temp21_22_23 = temp21_22;
    temp21_22_23 += temp[23];
    Ctxt temp20_21_22_23 = temp20_21_22;
    temp20_21_22_23 += temp[23];
    Ctxt temp24_25 = temp[24];
    temp24_25 += temp[25];
    Ctxt temp24_26 = temp[24];
    temp24_26 += temp[26];
    Ctxt temp24_27 = temp[24];
    temp24_27 += temp[27];
    Ctxt temp25_26 = temp[25];
    temp25_26 += temp[26];
    Ctxt temp25_27 = temp[25];
    temp25_27 += temp[27];
    Ctxt temp26_27 = temp[26];
    temp26_27 += temp[27];
    Ctxt temp24_25_26 = temp24_25;
    temp24_25_26 += temp[26];
    Ctxt temp24_25_27 = temp24_25;
    temp24_25_27 += temp[27];
    Ctxt temp24_26_27 = temp24_26;
    temp24_26_27 += temp[27];
    Ctxt temp25_26_27 = temp25_26;
    temp25_26_27 += temp[27];
    Ctxt temp24_25_26_27 = temp24_25_26;
    temp24_25_26_27 += temp[27];
    Ctxt temp28_29 = temp[28];
    temp28_29 += temp[29];
    Ctxt temp28_30 = temp[28];
    temp28_30 += temp[30];
    Ctxt temp28_31 = temp[28];
    temp28_31 += temp[31];
    Ctxt temp29_30 = temp[29];
    temp29_30 += temp[30];
    Ctxt temp29_31 = temp[29];
    temp29_31 += temp[31];
    Ctxt temp30_31 = temp[30];
    temp30_31 += temp[31];
    Ctxt temp28_29_30 = temp28_29;
    temp28_29_30 += temp[30];
    Ctxt temp28_29_31 = temp28_29;
    temp28_29_31 += temp[31];
    Ctxt temp28_30_31 = temp28_30;
    temp28_30_31 += temp[31];
    Ctxt temp29_30_31 = temp29_30;
    temp29_30_31 += temp[31];
    Ctxt temp28_29_30_31 = temp28_29_30;
    temp28_29_30_31 += temp[31];
    Ctxt temp32_33 = temp[32];
    temp32_33 += temp[33];
    Ctxt temp32_34 = temp[32];
    temp32_34 += temp[34];
    Ctxt temp32_35 = temp[32];
    temp32_35 += temp[35];
    Ctxt temp33_34 = temp[33];
    temp33_34 += temp[34];
    Ctxt temp33_35 = temp[33];
    temp33_35 += temp[35];
    Ctxt temp34_35 = temp[34];
    temp34_35 += temp[35];
    Ctxt temp32_33_34 = temp32_33;
    temp32_33_34 += temp[34];
    Ctxt temp32_33_35 = temp32_33;
    temp32_33_35 += temp[35];
    Ctxt temp32_34_35 = temp32_34;
    temp32_34_35 += temp[35];
    Ctxt temp33_34_35 = temp33_34;
    temp33_34_35 += temp[35];
    Ctxt temp32_33_34_35 = temp32_33_34;
    temp32_33_34_35 += temp[35];

    eData[0] += temp1_3;
    eData[0] += temp4_5_6_7;
    eData[0] += temp8_11;
    eData[0] += temp14_15;
    eData[0] += temp16_17_19;
    eData[0] += temp20_21_22;
    eData[0] += temp24_25;
    eData[0] += temp29_30_31;
    eData[0] += temp33_34_35;
    eData[1] += temp0_2_3;
    eData[1] += temp4_6;
    eData[1] += temp8_10;
    eData[1] += temp12_13_15;
    eData[1] += temp17_18;
    eData[1] += temp20_21_22_23;
    eData[1] += temp24_25_26;
    eData[1] += temp28_31;
    eData[1] += temp32_33_34;
    eData[2] = temp[1];
    eData[2] += temp4_5_7;
    eData[2] += temp8_9_10;
    eData[2] += temp12_14;
    eData[2] += temp16_17_18_19;
    eData[2] += temp20_21_23;
    eData[2] += temp25_26_27;
    eData[2] += temp28_29_30_31;
    eData[2] += temp32_33_35;
    eData[3] += temp0_1_2;
    eData[3] += temp4_6_7;
    eData[3] += temp8_9_10_11;
    eData[3] += temp[14];
    eData[3] += temp17_18_19;
    eData[3] += temp20_22_23;
    eData[3] += temp24_25_27;
    eData[3] += temp[28];
    eData[3] += temp32_33_34;
    eData[4] += temp0_1_3;
    eData[4] += temp5_6_7;
    eData[4] += temp9_11;
    eData[4] += temp13_15;
    eData[4] += temp16_18;
    eData[4] += temp20_21_23;
    eData[4] += temp24_25_26_27;
    eData[4] += temp28_29_31;
    eData[4] += temp34_35;
    eData[5] = temp0_2;
    eData[5] += temp4_7;
    eData[5] += temp8_10_11;
    eData[5] += temp12_13_15;
    eData[5] += temp17_19;
    eData[5] += temp20_21_22_23;
    eData[5] += temp24_26;
    eData[5] += temp28_29_30_31;
    eData[5] += temp32_33_34_35;
    eData[6] += temp0_1_3;
    eData[6] += temp4_5_7;
    eData[6] += temp9_10_11;
    eData[6] += temp12_13_14;
    eData[6] += temp[17];
    eData[6] += temp20_21_22_23;
    eData[6] += temp25_26_27;
    eData[6] += temp28_30_31;
    eData[6] += temp[35];
    eData[7] += temp1_2_3;
    eData[7] += temp4_6;
    eData[7] += temp8_9_10;
    eData[7] += temp12_14;
    eData[7] += temp16_18_19;
    eData[7] += temp21_23;
    eData[7] += temp24_26_27;
    eData[7] += temp28_29_30_31;
    eData[7] += temp32_34;
    eData[8] = temp0_1_2_3;
    eData[8] += temp5_7;
    eData[8] += temp10_11;
    eData[8] += temp13_14_15;
    eData[8] += temp16_18;
    eData[8] += temp20_22_23;
    eData[8] += temp24_25_26_27;
    eData[8] += temp29_31;
    eData[8] += temp32_33_34_35;
    eData[9] += temp2_3;
    eData[9] += temp4_6_7;
    eData[9] += temp8_10;
    eData[9] += temp12_13_14_15;
    eData[9] += temp16_17;
    eData[9] += temp20_23;
    eData[9] += temp24_25_26;
    eData[9] += temp28_29_30_31;
    eData[9] += temp33_34;
    eData[10] += temp[1];
    eData[10] += temp4_5_6_7;
    eData[10] += temp9_11;
    eData[10] += temp12_13_15;
    eData[10] += temp17_19;
    eData[10] += temp21_22;
    eData[10] += temp24_26_27;
    eData[10] += temp29_30_31;
    eData[10] += temp32_33_34_35;
    eData[11] = temp0_1_2_3;
    eData[11] += temp4_5_6;
    eData[11] += temp8_10;
    eData[11] += temp13_14;
    eData[11] += temp16_17_18_19;
    eData[11] += temp21_23;
    eData[11] += temp25_26_27;
    eData[11] += temp28_29_30;
    eData[11] += temp32_34_35;
    eData[12] += temp0_1;
    eData[12] += temp5_6_7;
    eData[12] += temp9_10_11;
    eData[12] += temp13_15;
    eData[12] += temp16_17_18_19;
    eData[12] += temp20_23;
    eData[12] += temp26_27;
    eData[12] += temp28_29_31;
    eData[12] += temp32_33_34;
    eData[13] += temp0_1_2;
    eData[13] += temp4_7;
    eData[13] += temp8_9_10;
    eData[13] += temp12_14_15;
    eData[13] += temp16_18;
    eData[13] += temp20_22;
    eData[13] += temp24_25_27;
    eData[13] += temp29_30;
    eData[13] += temp32_33_34_35;
    eData[14] = temp1_2_3;
    eData[14] += temp4_5_6_7;
    eData[14] += temp8_9_11;
    eData[14] += temp[13];
    eData[14] += temp16_17_19;
    eData[14] += temp20_21_22;
    eData[14] += temp24_26;
    eData[14] += temp28_29_30_31;
    eData[14] += temp32_33_35;
    eData[15] += temp0_1_3;
    eData[15] += temp[4];
    eData[15] += temp8_9_10;
    eData[15] += temp12_13_14;
    eData[15] += temp16_18_19;
    eData[15] += temp20_21_22_23;
    eData[15] += temp[26];
    eData[15] += temp29_30_31;
    eData[15] += temp32_34_35;
    eData[16] += temp0_1_2_3;
    eData[16] += temp4_5_7;
    eData[16] += temp10_11;
    eData[16] += temp12_13_15;
    eData[16] += temp17_18_19;
    eData[16] += temp21_23;
    eData[16] += temp25_27;
    eData[16] += temp28_30;
    eData[16] += temp32_33_35;
    eData[17] = temp0_2;
    eData[17] += temp4_5_6_7;
    eData[17] += temp8_9_10_11;
    eData[17] += temp12_14;
    eData[17] += temp16_19;
    eData[17] += temp20_22_23;
    eData[17] += temp24_25_27;
    eData[17] += temp29_31;
    eData[17] += temp32_33_34_35;
    eData[18] += temp1_2_3;
    eData[18] += temp4_6_7;
    eData[18] += temp[11];
    eData[18] += temp12_13_15;
    eData[18] += temp16_17_19;
    eData[18] += temp21_22_23;
    eData[18] += temp24_25_26;
    eData[18] += temp[29];
    eData[18] += temp32_33_34_35;
    eData[19] += temp0_2_3;
    eData[19] += temp4_5_6_7;
    eData[19] += temp8_10;
    eData[19] += temp13_14_15;
    eData[19] += temp16_18;
    eData[19] += temp20_21_22;
    eData[19] += temp24_26;
    eData[19] += temp28_30_31;
    eData[19] += temp33_35;
    eData[20] = temp0_1_2_3;
    eData[20] += temp5_7;
    eData[20] += temp8_9_10_11;
    eData[20] += temp12_13_14_15;
    eData[20] += temp17_19;
    eData[20] += temp22_23;
    eData[20] += temp25_26_27;
    eData[20] += temp28_30;
    eData[20] += temp32_34_35;
    eData[21] += temp0_1_2;
    eData[21] += temp4_5_6_7;
    eData[21] += temp9_10;
    eData[21] += temp14_15;
    eData[21] += temp16_18_19;
    eData[21] += temp20_22;
    eData[21] += temp24_25_26_27;
    eData[21] += temp28_29;
    eData[21] += temp32_35;
    eData[22] += temp0_2_3;
    eData[22] += temp5_6_7;
    eData[22] += temp8_9_10_11;
    eData[22] += temp[13];
    eData[22] += temp16_17_18_19;
    eData[22] += temp21_23;
    eData[22] += temp24_25_27;
    eData[22] += temp29_31;
    eData[22] += temp33_34;
    eData[23] = temp1_2_3;
    eData[23] += temp4_5_6;
    eData[23] += temp8_10_11;
    eData[23] += temp12_13_14_15;
    eData[23] += temp16_17_18;
    eData[23] += temp20_22;
    eData[23] += temp25_26;
    eData[23] += temp28_29_30_31;
    eData[23] += temp33_35;
    eData[24] += temp2_3;
    eData[24] += temp4_5_7;
    eData[24] += temp8_9_10;
    eData[24] += temp12_13;
    eData[24] += temp17_18_19;
    eData[24] += temp21_22_23;
    eData[24] += temp25_27;
    eData[24] += temp28_29_30_31;
    eData[24] += temp32_35;
    eData[25] += temp0_1_3;
    eData[25] += temp5_6;
    eData[25] += temp8_9_10_11;
    eData[25] += temp12_13_14;
    eData[25] += temp16_19;
    eData[25] += temp20_21_22;
    eData[25] += temp24_26_27;
    eData[25] += temp28_30;
    eData[25] += temp32_34;
    eData[26] = temp0_2;
    eData[26] += temp4_5_6_7;
    eData[26] += temp8_9_11;
    eData[26] += temp13_14_15;
    eData[26] += temp16_17_18_19;
    eData[26] += temp20_21_23;
    eData[26] += temp[25];
    eData[26] += temp28_29_31;
    eData[26] += temp32_33_34;
    eData[27] += temp[2];
    eData[27] += temp5_6_7;
    eData[27] += temp8_10_11;
    eData[27] += temp12_13_15;
    eData[27] += temp[16];
    eData[27] += temp20_21_22;
    eData[27] += temp24_25_26;
    eData[27] += temp28_30_31;
    eData[27] += temp32_33_34_35;
    eData[28] += temp1_3;
    eData[28] += temp4_6;
    eData[28] += temp8_9_11;
    eData[28] += temp12_13_14_15;
    eData[28] += temp16_17_19;
    eData[28] += temp22_23;
    eData[28] += temp24_25_27;
    eData[28] += temp29_30_31;
    eData[28] += temp33_35;
    eData[29] = temp0_1_3;
    eData[29] += temp5_7;
    eData[29] += temp8_9_10_11;
    eData[29] += temp12_14;
    eData[29] += temp16_17_18_19;
    eData[29] += temp20_21_22_23;
    eData[29] += temp24_26;
    eData[29] += temp28_31;
    eData[29] += temp32_34_35;
    eData[30] += temp0_1_2;
    eData[30] += temp[5];
    eData[30] += temp8_9_10_11;
    eData[30] += temp13_14_15;
    eData[30] += temp16_18_19;
    eData[30] += temp[23];
    eData[30] += temp24_25_27;
    eData[30] += temp28_29_31;
    eData[30] += temp33_34_35;
    eData[31] += temp0_2;
    eData[31] += temp4_6_7;
    eData[31] += temp9_11;
    eData[31] += temp12_14_15;
    eData[31] += temp16_17_18_19;
    eData[31] += temp20_22;
    eData[31] += temp25_26_27;
    eData[31] += temp28_30;
    eData[31] += temp32_33_34;
    eData[32] = temp1_2_3;
    eData[32] += temp4_6;
    eData[32] += temp8_10_11;
    eData[32] += temp12_13_14_15;
    eData[32] += temp17_19;
    eData[32] += temp20_21_22_23;
    eData[32] += temp24_25_26_27;
    eData[32] += temp29_31;
    eData[32] += temp34_35;
    eData[33] += temp0_1_2_3;
    eData[33] += temp4_5;
    eData[33] += temp8_11;
    eData[33] += temp12_13_14;
    eData[33] += temp16_17_18_19;
    eData[33] += temp21_22;
    eData[33] += temp26_27;
    eData[33] += temp28_30_31;
    eData[33] += temp32_34;
    eData[34] += temp0_1_3;
    eData[34] += temp5_7;
    eData[34] += temp9_10;
    eData[34] += temp12_14_15;
    eData[34] += temp17_18_19;
    eData[34] += temp20_21_22_23;
    eData[34] += temp[25];
    eData[34] += temp28_29_30_31;
    eData[34] += temp33_35;
    eData[35] = temp1_2;
    eData[35] += temp4_5_6_7;
    eData[35] += temp9_11;
    eData[35] += temp13_14_15;
    eData[35] += temp16_17_18;
    eData[35] += temp20_22_23;
    eData[35] += temp24_25_26_27;
    eData[35] += temp28_29_30;
    eData[35] += temp32_34;
}
// Compute the constants for Sbox,(x0,x1,x2)——>(x0,x0*x2+x1,-x0*x1+x0*x2+x2)
void HE_Sbox(vector<Ctxt> &eData)
{
    // (x0,x1,x2)——> (x0,x0*x2+x1,-x0*x1+x0*x2+x2)
    Ctxt c01 = eData[1];
    Ctxt c02 = eData[2];
    c01.multiplyBy(eData[0]);
    // c01*=eData[0];
    c02.multiplyBy(eData[0]);
    // c02*=eData[0];

    eData[1] += c02;
    eData[2] -= c01;
    eData[2] += c02;
    // omp_set_num_threads(12); // 设置线程数为12
    // #pragma omp parallel for
    for (long j = 3; j < BlockWords; j += 3)
    {
        c01 = eData[j + 1];
        c01.multiplyBy(eData[j]);

        c02 = eData[j + 2];
        c02.multiplyBy(eData[j]);

        eData[j + 1] += c02;
        eData[j + 2] -= c01;
        eData[j + 2] += c02;
    }
    // c01.cleanUp();
    // c02.cleanUp();
}
// Compute the constants for the last Sbox,(x0,x1,x2)——>(x0^3,x0*x2+x1,-x0*x1+x0*x2+x2)
void HE_Last_Sbox(vector<Ctxt> &eData)
{
    // (x0,x1,x2)——> (x0^3,x0*x2+x1,-x0*x1+x0*x2+x2)
    Ctxt c01 = eData[1];
    Ctxt c02 = eData[2];
    c01.multiplyBy(eData[0]);
    // c01*=eData[0];
    c02.multiplyBy(eData[0]);
    // c02*=eData[0];

    eData[1] += c02;
    eData[2] -= c01;
    eData[2] += c02;
    eData[0].cube();
    // omp_set_num_threads(12); // 设置线程数为12
    // #pragma omp parallel for
    for (long j = 3; j < BlockWords; j += 3)
    {
        c01 = eData[j + 1];
        c01.multiplyBy(eData[j]);

        c02 = eData[j + 2];
        c02.multiplyBy(eData[j]);

        eData[j + 1] += c02;
        eData[j + 2] -= c01;
        eData[j + 2] += c02;
        eData[j].cube();
    }
    // c01.cleanUp();
    // c02.cleanUp();
}

int main()
{
    std::cout << "Nr: " << Nr << std::endl;
    std::cout << "BlockWord: " << BlockWords << std::endl;
    std::cout << "BlockPlainWord: " << BlockPlainWords << std::endl;
    std::cout << "PlainMod: " << PlainMod << std::endl;
    std::cout << "PlainBlock: " << PlainBlock << std::endl;
    std::cout << "nslots: " << nslots << std::endl;
    std::cout << "Para_bits: " << Para_bits << std::endl;
    std::cout << "Para_c: " << Para_c << std::endl;
    //=============客户端offline阶段================
    // 定义初始向量
    vector<long> IV(BlockWords);
    for (unsigned i = 0; i < BlockWords; i++)
    {
        IV[i] = i + 1;
    }
    // 生成随机对称密钥
    random_device rd;
    vector<long> SymKey(BlockWords);
    for (unsigned i = 0; i < BlockWords; i++)
    {
        SymKey[i] = rd() % PlainMod;
    }
    std::cout << "SymKey generated." << std::endl;

    // Generating Public Key and encrypting the symmetric key
    auto start = std::chrono::high_resolution_clock::now();
    shared_ptr<Context> context(ContextBuilder<BGV>()
                                    .m(Para_m)
                                    .p(Para_p)
                                    .r(Para_r)
                                    .bits(Para_bits)
                                    .c(Para_c)
                                    .buildPtr());
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds_context = end - start;
    std::cout << "Context generation time: " << elapsed_seconds_context.count() << "s\n";

    auto start_PubKey = std::chrono::high_resolution_clock::now();
    SecKey secretKey(*context);
    secretKey.GenSecKey();
    unique_ptr<PubKey> publicKey = std::make_unique<helib::PubKey>(secretKey);
    // Compute key-switching matrices that we need
    // helib::addSome1DMatrices(secretKey);
    auto end_PubKey = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds_PubKey = end_PubKey - start_PubKey;
    std::cout << "PublicKey generation time: " << elapsed_seconds_PubKey.count() << "s\n";

    // 创建PAlgebra对象
    const helib::PAlgebra &zMStar = context->getZMStar();
    // long minv = -4;
    long root = 3; // m-th root of unity modulo p
    if (Para_m == 65536)
    {
        root = 3;
    }
    if (Para_m == 32768)
    {
        root = 9;
    }
    if (Para_m == 16384)
    {
        root = 81;
    }
    // 初始化Cmodulus对象
    helib::Cmodulus cmodulus(zMStar, Para_p, root);

    // 输出 context
    context->printout();
    std::cout << std::endl;

    long Qbits = context->bitSizeOfQ();
    double SecLevel = context->securityLevel();

    auto start_keyEncryption = std::chrono::high_resolution_clock::now();
    vector<Ctxt> encryptedSymKey;
    encryptSymKey(encryptedSymKey, SymKey, publicKey, cmodulus, nslots);
    auto end_keyEncryption = std::chrono::high_resolution_clock::now();
    double keyEncryption = std::chrono::duration<double>(end_keyEncryption - start_keyEncryption).count();
    std::cout << "SymKey FHE time: " << keyEncryption << "s\n";
    //  解密验证
    if (symkeyflag)
    {
        if (!verify_encryptSymKey(encryptedSymKey, SymKey, secretKey, cmodulus))
        {
            std::cerr << "Symmetric key encryption failed!" << std::endl;
            return 0;
        }
        std::cout << "Symmetric key encryption succeeded!" << std::endl;
    }
            // 获取当前源文件名
    std::string cppFileName = __FILE__;
    std::string baseFileName = cppFileName.substr(cppFileName.find_last_of("/\\") + 1);
    std::string outputFileName = baseFileName + ".txt";
    // 设置输出文件路径
    std::string dirPath = "../tests";
    std::string filePath;
    if (!fs::exists(dirPath))
    {
        filePath = outputFileName;
    }
    else
    {
        filePath = dirPath + "/" + outputFileName;
    }
    std::ofstream outfile(filePath, std::ios::app);
    if (!outfile)
    {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return 0;
    }
    // 获取当前时间
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    // 格式化时间并写入文件
    outfile << "Current time: " << std::put_time(std::localtime(&now_time), "%Y-%m-%d %H:%M:%S") << std::endl;
    // 获取当前 CPU 型号并写入文件
    std::string cpu_model = get_cpu_model();
    outfile << "CPU model: " << cpu_model << std::endl;
    // 获取内存信息并写入文件
    std::string memory_info = get_memory_info();
    outfile << "Memory info: " << std::endl
            << memory_info;
    // 获取系统版本信息并写入文件
    std::string os_version = get_os_version();
    outfile << "OS version: " << os_version;
    // 获取环境信息
    std::string environment_info = get_environment_info();
    outfile << "Environment info: " << environment_info << std::endl; 
    for (int test = 0; test < 10; test++)
    {
        std::cout << "--------------- Test = " << test << "---------------"<< std::endl;
        // Generating key stream
        std::vector<long> NonceSet(PlainBlock);
        vector<std::vector<long>> RoundKeySet(Nr + 1, std::vector<long>(KeyStreamWords, 0));
        std::vector<long> KeyStream(PlainWords);
        std::vector<long> RoundKey(BlockWords);
        long nonce;
        long block_num;
        std::vector<long> state(BlockWords); // 初始化 state
        Keccak_HashInstance shake128_2;
        std::cout << "Generating KeyStream..." << std::endl;
        auto start_keyStream = std::chrono::high_resolution_clock::now();
        uint64_t start_cycle1 = rdtsc();
        for (long counter = counter_begin; counter <= counter_end; counter++)
        {
            nonce = generate_secure_random_int(NonceSize);
            random_init_shake(nonce, counter, shake128_2);
            block_num = counter - counter_begin;
            NonceSet[block_num] = nonce;
            // 逐轮进行加密
            for (unsigned r = 0; r <= Nr; r++)
            {
                // 计算RoundKey
                for (unsigned i = 0; i < BlockWords; ++i)
                {
                    RoundKey[i] = (SymKey[i] * generate_random_field_element(shake128_2, false, max_prime_size, PlainMod)) % PlainMod;
                }
                // 将RoundKey 复制到RoundKeySet
                // 测试使用
                if (deflag)
                {
                    memcpy(&RoundKeySet[r][BlockWords * (block_num)], RoundKey.data(), BlockWords * sizeof(long));
                }
                if (r == 0)
                { // 初始轮
                    for (unsigned i = 0; i < BlockWords; i++)
                    {
                        state[i] = (RoundKey[i] + IV[i]) % PlainMod;
                    }
                }
                else if (r < Nr)
                { // 常规轮

                    yusP.M36(state);  // MixColumns
                    yusP.Sbox(state); // S盒
                    for (unsigned i = 0; i < BlockWords; i++)
                    {
                        state[i] = (state[i] + RoundKey[i]) % PlainMod;
                    }
                }
                else
                {                     // 最后一轮
                    yusP.M36(state);  // Liner
                    yusP.Sbox(state); // S盒
                    for (unsigned i = 0; i < BlockWords; i++)
                    {
                        state[i] = (state[i] + RoundKey[i]) % PlainMod;
                    }
                    yusP.M36(state); // Liner
                    memcpy(&KeyStream[(block_num)*BlockPlainWords], state.data(), BlockPlainWords * sizeof(long));
                }
            }
        }
        uint64_t end_cycle1 = rdtsc();
        auto end_keyStream = std::chrono::high_resolution_clock::now();
        double Client_offtime = std::chrono::duration_cast<std::chrono::duration<double>>(end_keyStream - start_keyStream).count();
        std::cout << "Encryption offline total time: " << Client_offtime << "s\n";
        uint64_t cycle_count1 = end_cycle1 - start_cycle1;
        std::cout << "Encryption offline total cycles: " << cycle_count1 << std::endl;
        // 将KeyStream同态密文写入文件
        // if (!writeEncryptedSymKey(encryptedSymKey, "Client_encryptedSymKey.bin"))
        // {
        //     return false;
        // }
        // std::cout << "Client_encryptedSymKey.bin has been written to file." << std::endl;

        //=============服务端offline阶段================
        std::cout << "Generating XOF stream..." << std::endl;
        std::vector<vec_long> Xset(NrBlockWords);
        for (int i = 0; i < NrBlockWords; i++)
        {
            Xset[i].SetLength(nslots);
            for (int j = 0; j < nslots; j++)
            {
                Xset[i][j] = 0;
            }
        }
        long rB;
        Keccak_HashInstance shake128_3;
        auto start_XOF = std::chrono::high_resolution_clock::now();
        for (long counter = counter_begin; counter <= counter_end; counter++)
        {
            block_num = counter - counter_begin;
            nonce = NonceSet[block_num];
            random_init_shake(nonce, counter, shake128_3);
            // 逐轮进行加密
            for (unsigned r = 0; r <= Nr; r++)
            {
                rB = r * BlockWords;
                // 计算 Xc
                for (unsigned i = 0; i < BlockWords; ++i)
                {
                    Xset[rB + i][block_num] = generate_random_field_element(shake128_3, false, max_prime_size, PlainMod);
                }
            }
        }
        if (PlainBlock == 1)
        {
            for (int i = 0; i < NrBlockWords; i++)
            {
                for (int j = 1; j < nslots; j++)
                {
                    Xset[i][j] = Xset[i][0];
                }
            }
        }
        auto end_XOF = std::chrono::high_resolution_clock::now();
        double XOF_time = std::chrono::duration<double>(end_XOF - start_XOF).count();
        std::cout << "XOF stream Generation time: " << XOF_time << "s\n";

        zz_pX encodedtemp;
        int noise_budget = min_noise_budget(encryptedSymKey);
        std::cout << "noise budget initially: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            continue;
        }
        // 生成 encryptedKeyStream
        // 定义Add_time、Sbox_time、Linear_time
        double Sbox_time = 0, Linear_time = 0, Add_time = 0;
        auto start_sbox = std::chrono::high_resolution_clock::now();
        auto end_sbox = std::chrono::high_resolution_clock::now();
        auto start_linear = std::chrono::high_resolution_clock::now();
        auto end_linear = std::chrono::high_resolution_clock::now();
        auto start_add = std::chrono::high_resolution_clock::now();
        auto end_add = std::chrono::high_resolution_clock::now();
        vector<double> sbox_set(Nr);
        vector<double> linear_set(Nr + 1);
        vector<Ctxt> encryptedKeyStream = encryptedSymKey;
        std::cout << "whiteround start" << std::endl;
        zzX encodedXseti;
        start_add = std::chrono::high_resolution_clock::now();
        for (long i = 0; i < BlockWords; i++)
        { // encrypt the encoded key
            cmodulus.iFFT(encodedtemp, Xset[i]);
            convert(encodedXseti, encodedtemp);
            encryptedKeyStream[i].multByConstant(encodedXseti);
            encryptedKeyStream[i].addConstant(IV[i]);
        }
        end_add = std::chrono::high_resolution_clock::now();
        Add_time += std::chrono::duration<double>(end_add - start_add).count();
        // 输出 Add_time
        std::cout << "whiteround time: " << Add_time << "s\n";
        noise_budget = min_noise_budget(encryptedKeyStream);
        std::cout << "noise budget after Whiteround: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            continue;
        }
        // 明文密钥流
        vector<long> KeyStream2(KeyStreamWords);
        if (deflag)
        {
            // 对IV和RoundKeySet进行异或
            for (long i = 0; i < KeyStreamWords; i++)
            {
                KeyStream2[i] = (IV[i % BlockWords] + RoundKeySet[0][i]) % PlainMod;
            }
            // 使用 verifyDecryption 函数解密并验证 KeyStream2
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords, PlainBlock, nslots, Para_p))
            {
                std::cerr << "Decryption verification failed for KeyStream2." << std::endl;
                continue;
            }
            std::cout << "Decryption verification succeeded for whiteround." << std::endl;
        }
        Ctxt tmpCtxt(*publicKey);
        for (long r = 1; r < Nr; r++)
        {
            std::cout << "Round " << r << " start" << std::endl;
            start_linear = std::chrono::high_resolution_clock::now();
            // Linear Layer
            HE_M(encryptedKeyStream);
            end_linear = std::chrono::high_resolution_clock::now();
            Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
            linear_set[r - 1] = std::chrono::duration<double>(end_linear - start_linear).count();
            noise_budget = min_noise_budget(encryptedKeyStream);
            std::cout << "noise budget after linear: " << noise_budget << std::endl;
            if (noise_budget <= 0)
            {
                std::cerr << "noise budget is not enough!!!" << std::endl;
                continue;
            }
            if (deflag)
            {
                for (int i = 0; i < PlainBlock; i++)
                {
                    vector<long> tmp(BlockWords);
                    for (int j = 0; j < BlockWords; j++)
                    {
                        tmp[j] = KeyStream2[i * BlockWords + j];
                    }
                    yusP.M36(tmp);
                    for (int j = 0; j < BlockWords; j++)
                    {
                        KeyStream2[i * BlockWords + j] = tmp[j];
                    }
                }
                if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords, PlainBlock, nslots, Para_p))
                {
                    std::cerr << "Decryption verification failed for KeyStream2 Linear Layer." << std::endl;
                    continue;
                }
                std::cout << "Decryption verification succeeded for KeyStream2 Linear Layer." << std::endl;
            }
            start_sbox = std::chrono::high_resolution_clock::now();
            // S Layer
            HE_Sbox(encryptedKeyStream);
            end_sbox = std::chrono::high_resolution_clock::now();
            Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
            sbox_set[r - 1] = std::chrono::duration<double>(end_sbox - start_sbox).count();
            noise_budget = min_noise_budget(encryptedKeyStream);
            std::cout << "noise budget after sbox: " << noise_budget << std::endl;
            if (noise_budget <= 0)
            {
                std::cerr << "noise budget is not enough!!!" << std::endl;
                continue;
            }
            if (deflag)
            {
                yusP.Sbox(KeyStream2);
                if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords, PlainBlock, nslots, Para_p))
                {
                    std::cerr << "Decryption verification failed for KeyStream2 Sbox." << std::endl;
                    continue;
                }
                std::cout << "Decryption verification succeeded for KeyStream2 Sbox." << std::endl;
            }
            // Round Key Addition
            start_add = std::chrono::high_resolution_clock::now();
            // omp_set_num_threads(12); // 设置线程数为12
            // #pragma omp parallel for
            for (long j = 0; j < BlockWords; j++)
            {
                cmodulus.iFFT(encodedtemp, Xset[r * BlockWords + j]);
                convert(encodedXseti, encodedtemp);
                tmpCtxt = encryptedSymKey[j];
                tmpCtxt.multByConstant(encodedXseti);
                encryptedKeyStream[j] += tmpCtxt;
            }
            end_add = std::chrono::high_resolution_clock::now();
            Add_time += std::chrono::duration<double>(end_add - start_add).count();
            noise_budget = min_noise_budget(encryptedKeyStream);
            std::cout << "noise budget after Add: " << noise_budget << std::endl;
            if (noise_budget <= 0)
            {
                std::cerr << "noise budget is not enough!!!" << std::endl;
                continue;
            }
            if (deflag)
            {
                for (long i = 0; i < KeyStreamWords; i++)
                {
                    KeyStream2[i] = (KeyStream2[i] + RoundKeySet[r][i]) % PlainMod;
                }
                if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords, PlainBlock, nslots, Para_p))
                {
                    std::cerr << "Decryption verification failed for KeyStream2 Round Key Addition." << std::endl;
                    continue;
                }
                std::cout << "Decryption verification succeeded for KeyStream2 Round Key Addition." << std::endl;
            }
        }
        // 最后一轮
        std::cout << "the last Round " << Nr << " start" << std::endl;
        start_linear = std::chrono::high_resolution_clock::now();
        // Linear Layer
        HE_M(encryptedKeyStream);
        end_linear = std::chrono::high_resolution_clock::now();
        Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
        linear_set[Nr - 1] = std::chrono::duration<double>(end_linear - start_linear).count();
        noise_budget = min_noise_budget(encryptedKeyStream);
        std::cout << "noise budget after linear: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            continue;
        }
        if (deflag)
        {
            for (int i = 0; i < PlainBlock; i++)
            {
                vector<long> tmp(BlockWords);
                for (int j = 0; j < BlockWords; j++)
                {
                    tmp[j] = KeyStream2[i * BlockWords + j];
                }
                yusP.M36(tmp);
                for (int j = 0; j < BlockWords; j++)
                {
                    KeyStream2[i * BlockWords + j] = tmp[j];
                }
            }
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords, PlainBlock, nslots, Para_p))
            {
                std::cerr << "Decryption verification failed for KeyStream2 Linear Layer." << std::endl;
                continue;
            }
            std::cout << "Decryption verification succeeded for KeyStream2 Linear Layer." << std::endl;
        }
        start_sbox = std::chrono::high_resolution_clock::now();
        // S Layer
        HE_Sbox(encryptedKeyStream);
        // HE_Last_Sbox(encryptedKeyStream);
        end_sbox = std::chrono::high_resolution_clock::now();
        Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
        sbox_set[Nr - 1] = std::chrono::duration<double>(end_sbox - start_sbox).count();
        noise_budget = min_noise_budget(encryptedKeyStream);
        std::cout << "noise budget after sbox: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            continue;
        }
        if (deflag)
        {
            yusP.Sbox(KeyStream2);
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords, PlainBlock, nslots, Para_p))
            {
                std::cerr << "Decryption verification failed for KeyStream2 Sbox." << std::endl;
                continue;
            }
            std::cout << "Decryption verification succeeded for KeyStream2 Sbox." << std::endl;
        }
        // add
        start_add = std::chrono::high_resolution_clock::now();
        for (long j = 0; j < BlockWords; j++)
        {
            cmodulus.iFFT(encodedtemp, Xset[Nr * BlockWords + j]);
            convert(encodedXseti, encodedtemp);
            tmpCtxt = encryptedSymKey[j];
            tmpCtxt.multByConstant(encodedXseti);
            encryptedKeyStream[j] += tmpCtxt;
        }
        end_add = std::chrono::high_resolution_clock::now();
        Add_time += std::chrono::duration<double>(end_add - start_add).count();
        noise_budget = min_noise_budget(encryptedKeyStream);
        std::cout << "noise budget after Add: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            continue;
        }
        if (deflag)
        {
            for (long i = 0; i < KeyStreamWords; i++)
            {
                KeyStream2[i] = (KeyStream2[i] + RoundKeySet[Nr][i]) % PlainMod;
            }
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords, PlainBlock, nslots, Para_p))
            {
                std::cerr << "Decryption verification failed for KeyStream2 Round Key Addition." << std::endl;
                continue;
            }
            std::cout << "Decryption verification succeeded for KeyStream2 Round Key Addition." << std::endl;
        }
        start_linear = std::chrono::high_resolution_clock::now();
        // Linear Layer
        HE_M(encryptedKeyStream);
        end_linear = std::chrono::high_resolution_clock::now();
        Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
        linear_set[Nr] = std::chrono::duration<double>(end_linear - start_linear).count();
        noise_budget = min_noise_budget(encryptedKeyStream);
        std::cout << "noise budget after linear: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            continue;
        }
        if (deflag)
        {
            for (int i = 0; i < PlainBlock; i++)
            {
                vector<long> tmp(BlockWords);
                for (int j = 0; j < BlockWords; j++)
                {
                    tmp[j] = KeyStream2[i * BlockWords + j];
                }
                yusP.M36(tmp);
                for (int j = 0; j < BlockWords; j++)
                {
                    KeyStream2[i * BlockWords + j] = tmp[j];
                }
            }
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords, PlainBlock, nslots, Para_p))
            {
                std::cerr << "Decryption verification failed for KeyStream2 Linear Layer." << std::endl;
                continue;
            }
            std::cout << "Decryption verification succeeded for KeyStream2 Linear Layer." << std::endl;
        }
        // 截断密钥流
        encryptedKeyStream.erase(encryptedKeyStream.begin() + BlockPlainWords, encryptedKeyStream.end());
        // 输出 XOF_time,Add_time、Sbox_time、Linear_time
        std::cout << "XOF time: " << XOF_time << "s\n";
        std::cout << "Add time: " << Add_time << "s\n";
        std::cout << "Sbox time: " << Sbox_time << "s\n";
        std::cout << "Linear time: " << Linear_time << "s\n";
        // 计算总时间
        double Server_offtime = XOF_time + Add_time + Sbox_time + Linear_time;
        std::cout << "Server offline total time: " << Server_offtime << "s\n";
        std::cout << "sbox_timeset: " << sbox_set << endl;
        std::cout << "linear_timeset: " << linear_set << endl;
        if (KeyStreamflag)
        {
            if (!verifyDecryption(encryptedKeyStream, KeyStream, secretKey, cmodulus, BlockPlainWords, PlainBlock, nslots, Para_p))
            {
                std::cerr << "Decryption verification failed for KeyStream." << std::endl;
                continue;
            }
            std::cout << "Decryption verification succeeded for KeyStream." << std::endl;
        }
        // 客户端在线
        // 生成随机对称明文流，只用于测试
        random_device rd;
        vector<long> PlainStream(PlainWords);
        for (int i = 0; i < PlainWords; i++)
        {
            PlainStream[i] = rd() % PlainMod;
        }
        // 加密
        vector<vec_long> CipherStream(BlockPlainWords);
        for (int i = 0; i < BlockPlainWords; i++)
        {
            CipherStream[i].SetLength(nslots);
            for (int j = 0; j < nslots; j++)
            {
                CipherStream[i][j] = 0;
            }
        }
        long jBP;
        auto start_ClientOnline = std::chrono::high_resolution_clock::now();
        uint64_t start_cycle2 = rdtsc();
        for (int j = 0; j < PlainBlock; j++)
        {
            jBP = j * BlockPlainWords;
            for (int i = 0; i < BlockPlainWords; i++)
            {
                CipherStream[i][j] = (PlainStream[jBP + i] - KeyStream[jBP + i]) % PlainMod;
            }
        }
        if (PlainBlock == 1)
        {

            for (int i = 0; i < BlockPlainWords; i++)
            {
                for (int j = 1; j < nslots; j++)
                {
                    CipherStream[i][j] = CipherStream[i][0];
                }
            }
        }
        uint64_t end_cycle2 = rdtsc();
        auto end_ClientOnline = std::chrono::high_resolution_clock::now();
        uint64_t cycle_count2 = end_cycle2 - start_cycle2;
        uint64_t cycle_count = cycle_count1 + cycle_count2;
        double Client_ontime = std::chrono::duration<double>(end_ClientOnline - start_ClientOnline).count();
        std::cout << "Client onine total time:" << std::chrono::duration<double>(end_ClientOnline - start_ClientOnline).count() << "s\n";
        double Client_totaltime = Client_offtime + Client_ontime;
        std::cout << "Client total time: " << Client_totaltime << "s\n";
        // 输出吞吐量
        double Cli_throughput = (8 * cycle_count) / Plainbits; // 单位：Cycle/Byte
        std::cout << "Client Throughput: " << Cli_throughput << "Cycle/Byte\n";
        // 服务端在线
        // 同态加密
        vector<Ctxt> encrypedPlainStream = encryptedKeyStream;
        // 对CipherStream进行编码
        ZZX encodedCipherStreami;
        auto start_ServerOnline = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < BlockPlainWords; i++)
        {
            cmodulus.iFFT(encodedtemp, CipherStream[i]);
            convert(encodedCipherStreami, encodedtemp);
            encrypedPlainStream[i].addConstant(encodedCipherStreami);
        }
        auto end_ServerOnline = std::chrono::high_resolution_clock::now();
        double server_ontime = std::chrono::duration<double>(end_ServerOnline - start_ServerOnline).count();
        std::cout << "Server onine total time:" << server_ontime << "s\n";
        noise_budget = min_noise_budget(encryptedKeyStream);
        std::cout << "noise budget after online: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            continue;
        }
        // 同态解密验证
        // for (int i = 0; i < encryptedKeyStream.size(); i++)
        // {
        //     encryptedKeyStream[i].bringToSet(encryptedKeyStream[i].naturalPrimeSet());
        // }
        if (plainflag)
        {
            if (!verifyDecryption(encrypedPlainStream, PlainStream, secretKey, cmodulus, BlockPlainWords, PlainBlock, nslots, Para_p))
            {
                std::cerr << "Decryption verification failed for encrypedPlainStream." << std::endl;
                continue;
            }
            std::cout << "Decryption verification succeeded for encrypedPlainStream." << std::endl;
        }
        // 计算吞吐量,KiB/min
        double Server_totaltime = Server_offtime + server_ontime;
        double Ser_throughput = (Plainbits * 60) / (pow(2, 13) * Server_totaltime);
        std::cout << "Server total time: " << Server_totaltime << "s\n";
        std::cout << "Server Throughput: " << Ser_throughput << "KiB/min\n";

        outfile << std::left << std::setw(3) << Nr
                << std::left << std::setw(12) << Para_p
                << std::left << std::setw(10) << nslots
                << std::left << std::setw(10) << PlainBlock
                << std::left << std::setw(10) << Qbits
                << std::left << std::setw(10) << SecLevel
                << std::left << std::setw(15) << cycle_count1
                << std::left << std::setw(15) << cycle_count2
                << std::left << std::setw(15) << cycle_count
                << std::fixed << std::setprecision(3)
                << std::left << std::setw(15) << Cli_throughput
                << std::left << std::setw(15) << Server_offtime
                << std::left << std::setw(15) << server_ontime
                << std::left << std::setw(15) << Server_totaltime
                << std::left << std::setw(15) << Ser_throughput
                << std::left << std::setw(10) << noise_budget
                << std::endl;
        
        std::cout << "Test " << test << " finished." << std::endl;
    }outfile.close();
    return 0;
}