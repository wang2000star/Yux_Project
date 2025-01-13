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


#include <immintrin.h> // 包含 SIMD 指令集支持

#include <SEAL-4.1/seal/seal.h>

#include "../utils/random_bit.hpp"
#include "../Symmetric/Yus_p.hpp"
#include "../utils/tool.hpp"

using namespace std;
using namespace seal;

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
constexpr bool Rkflag = 1;     // true/1表示乘法，false/0表示加法，指的是随机向量和密钥间的操作
constexpr bool deflag = 0;     // true/1表示进行每一步解密验证，false/0表示不进行每一步解密验证
constexpr bool ompflag = 0;    // true/1表示使用OpenMP并行编码，false/0表示不使用OpenMP并行编码
constexpr bool symkeyflag = 0; // true/1表示对称密钥同态解密验证加密，false/0表示不验证
constexpr bool KeyStreamflag = 0;  // true/1表示密钥流同态解密验证，false/0表示不验证
constexpr bool plainflag = 0;      // true/1表示明文同态解密验证，false/0表示不验证
// 参数设置，paramMap[Nr-4][idx]
constexpr unsigned Nr = 4; // 轮数
constexpr long idx = 0;
constexpr unsigned Sbox_depth = 1 * Nr; // S盒深度
// 当c=2时，Qbits=1.5*bits,当c=3时，Qbits=1.5*bits - 100
// 硬编码参数值
constexpr tuple<long, long> paramMap[1][8] = {
    {// Nr = 4
     // {p, log2(m)}
     {65537, 15},   // 0 *
     {163841, 15},  // 1
     {65537, 14},   // 2 *
     {163841, 14},  // 3
     {65537, 15},   // 4
     {163841, 15},  // 5
     {65537, 15},   // 6
     {163841, 15}
    } // 7
    };
// p=k*m+1
//  2^10=1024,2^11=2048,2^12=4096,2^13=8192,2^14=16384,2^15=32768,2^16=65536,
// 当电脑内存有限时，log2Para_m太大，会导致内存不足，出现terminate called recursively错误，从而终止程序.
// 此外，即使正常运行，由于内存一直临界，会导致程序运行速度变慢，时间测量不准确。

constexpr long log2Para_m = get<1>(paramMap[0][idx]) - 0;
constexpr long Para_p = get<0>(paramMap[0][idx]);    // plaintext prime
constexpr long Para_m = 1 << log2Para_m;                  // cyclotomic polynomial
constexpr long phi_m = Para_m >> 1;                       // phi(m)
constexpr long Para_r = 1;                                // Lifting [defualt = 1]
//!!!!!!!!!!!!!!!
constexpr long nslots = phi_m;             // 槽数
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
constexpr unsigned NrBlockWords = BlockWords * (Nr + 1);  // Nr轮分组密钥字节长度

constexpr long KeyStreamWords = BlockWords * PlainBlock; // 密钥流字节长度
constexpr long PlainWords = BlockPlainWords * PlainBlock;      // 明文字节长度

constexpr long Plainbits = Wordbits * PlainWords;        // 明文比特长度

constexpr long max_prime_size = (1ULL << (Wordbits-1)) - 1;

constexpr unsigned NonceSize = 32;                           // Nonce比特长度
constexpr long counter_begin = 0;                            // 计数器起始值
constexpr long counter_end = PlainBlock + counter_begin - 1; // 计数器结束值

YusP yusP(PlainMod); // 构建明文对称加密实例

// Linear transformation
void HE_M(vector<Ciphertext> &eData,Evaluator &evaluator,RelinKeys &relin_keys)
{
    vector<Ciphertext> temp = eData;
    Ciphertext temp0_1 = temp[0];
    evaluator.add_inplace(temp0_1, temp[1]);
    Ciphertext temp0_2 = temp[0];
    evaluator.add_inplace(temp0_2, temp[2]);
    Ciphertext temp0_3 = temp[0];
    evaluator.add_inplace(temp0_3, temp[3]);
    Ciphertext temp1_2 = temp[1];
    evaluator.add_inplace(temp1_2, temp[2]);
    Ciphertext temp1_3 = temp[1];
    evaluator.add_inplace(temp1_3, temp[3]);
    Ciphertext temp2_3 = temp[2];
    evaluator.add_inplace(temp2_3, temp[3]);
    Ciphertext temp0_1_2 = temp0_1;
    evaluator.add_inplace(temp0_1_2, temp[2]);
    Ciphertext temp0_1_3 = temp0_1;
    evaluator.add_inplace(temp0_1_3, temp[3]);
    Ciphertext temp0_2_3 = temp0_2;
    evaluator.add_inplace(temp0_2_3, temp[3]);
    Ciphertext temp1_2_3 = temp1_2;
    evaluator.add_inplace(temp1_2_3, temp[3]);
    Ciphertext temp0_1_2_3 = temp0_1_2;
    evaluator.add_inplace(temp0_1_2_3, temp[3]);
    Ciphertext temp4_5 = temp[4];
    evaluator.add_inplace(temp4_5, temp[5]);
    Ciphertext temp4_6 = temp[4];
    evaluator.add_inplace(temp4_6, temp[6]);
    Ciphertext temp4_7 = temp[4];
    evaluator.add_inplace(temp4_7, temp[7]);
    Ciphertext temp5_6 = temp[5];
    evaluator.add_inplace(temp5_6, temp[6]);
    Ciphertext temp5_7 = temp[5];
    evaluator.add_inplace(temp5_7, temp[7]);
    Ciphertext temp6_7 = temp[6];
    evaluator.add_inplace(temp6_7, temp[7]);
    Ciphertext temp4_5_6 = temp4_5;
    evaluator.add_inplace(temp4_5_6, temp[6]);
    Ciphertext temp4_5_7 = temp4_5;
    evaluator.add_inplace(temp4_5_7, temp[7]);
    Ciphertext temp4_6_7 = temp4_6;
    evaluator.add_inplace(temp4_6_7, temp[7]);
    Ciphertext temp5_6_7 = temp5_6;
    evaluator.add_inplace(temp5_6_7, temp[7]);
    Ciphertext temp4_5_6_7 = temp4_5_6;
    evaluator.add_inplace(temp4_5_6_7, temp[7]);
    Ciphertext temp8_9 = temp[8];
    evaluator.add_inplace(temp8_9, temp[9]);
    Ciphertext temp8_10 = temp[8];
    evaluator.add_inplace(temp8_10, temp[10]);
    Ciphertext temp8_11 = temp[8];
    evaluator.add_inplace(temp8_11, temp[11]);
    Ciphertext temp9_10 = temp[9];
    evaluator.add_inplace(temp9_10, temp[10]);
    Ciphertext temp9_11 = temp[9];
    evaluator.add_inplace(temp9_11, temp[11]);
    Ciphertext temp10_11 = temp[10];
    evaluator.add_inplace(temp10_11, temp[11]);
    Ciphertext temp8_9_10 = temp8_9;
    evaluator.add_inplace(temp8_9_10, temp[10]);
    Ciphertext temp8_9_11 = temp8_9;
    evaluator.add_inplace(temp8_9_11, temp[11]);
    Ciphertext temp8_10_11 = temp8_10;
    evaluator.add_inplace(temp8_10_11, temp[11]);
    Ciphertext temp9_10_11 = temp9_10;
    evaluator.add_inplace(temp9_10_11, temp[11]);
    Ciphertext temp8_9_10_11 = temp8_9_10;
    evaluator.add_inplace(temp8_9_10_11, temp[11]);
    Ciphertext temp12_13 = temp[12];
    evaluator.add_inplace(temp12_13, temp[13]);
    Ciphertext temp12_14 = temp[12];
    evaluator.add_inplace(temp12_14, temp[14]);
    Ciphertext temp12_15 = temp[12];
    evaluator.add_inplace(temp12_15, temp[15]);
    Ciphertext temp13_14 = temp[13];
    evaluator.add_inplace(temp13_14, temp[14]);
    Ciphertext temp13_15 = temp[13];
    evaluator.add_inplace(temp13_15, temp[15]);
    Ciphertext temp14_15 = temp[14];
    evaluator.add_inplace(temp14_15, temp[15]);
    Ciphertext temp12_13_14 = temp12_13;
    evaluator.add_inplace(temp12_13_14, temp[14]);
    Ciphertext temp12_13_15 = temp12_13;
    evaluator.add_inplace(temp12_13_15, temp[15]);
    Ciphertext temp12_14_15 = temp12_14;
    evaluator.add_inplace(temp12_14_15, temp[15]);
    Ciphertext temp13_14_15 = temp13_14;
    evaluator.add_inplace(temp13_14_15, temp[15]);
    Ciphertext temp12_13_14_15 = temp12_13_14;
    evaluator.add_inplace(temp12_13_14_15, temp[15]);
    Ciphertext temp16_17 = temp[16];
    evaluator.add_inplace(temp16_17, temp[17]);
    Ciphertext temp16_18 = temp[16];
    evaluator.add_inplace(temp16_18, temp[18]);
    Ciphertext temp16_19 = temp[16];
    evaluator.add_inplace(temp16_19, temp[19]);
    Ciphertext temp17_18 = temp[17];
    evaluator.add_inplace(temp17_18, temp[18]);
    Ciphertext temp17_19 = temp[17];
    evaluator.add_inplace(temp17_19, temp[19]);
    Ciphertext temp18_19 = temp[18];
    evaluator.add_inplace(temp18_19, temp[19]);
    Ciphertext temp16_17_18 = temp16_17;
    evaluator.add_inplace(temp16_17_18, temp[18]);
    Ciphertext temp16_17_19 = temp16_17;
    evaluator.add_inplace(temp16_17_19, temp[19]);
    Ciphertext temp16_18_19 = temp16_18;
    evaluator.add_inplace(temp16_18_19, temp[19]);
    Ciphertext temp17_18_19 = temp17_18;
    evaluator.add_inplace(temp17_18_19, temp[19]);
    Ciphertext temp16_17_18_19 = temp16_17_18;
    evaluator.add_inplace(temp16_17_18_19, temp[19]);
    Ciphertext temp20_21 = temp[20];
    evaluator.add_inplace(temp20_21, temp[21]);
    Ciphertext temp20_22 = temp[20];
    evaluator.add_inplace(temp20_22, temp[22]);
    Ciphertext temp20_23 = temp[20];
    evaluator.add_inplace(temp20_23, temp[23]);
    Ciphertext temp21_22 = temp[21];
    evaluator.add_inplace(temp21_22, temp[22]);
    Ciphertext temp21_23 = temp[21];
    evaluator.add_inplace(temp21_23, temp[23]);
    Ciphertext temp22_23 = temp[22];
    evaluator.add_inplace(temp22_23, temp[23]);
    Ciphertext temp20_21_22 = temp20_21;
    evaluator.add_inplace(temp20_21_22, temp[22]);
    Ciphertext temp20_21_23 = temp20_21;
    evaluator.add_inplace(temp20_21_23, temp[23]);
    Ciphertext temp20_22_23 = temp20_22;
    evaluator.add_inplace(temp20_22_23, temp[23]);
    Ciphertext temp21_22_23 = temp21_22;
    evaluator.add_inplace(temp21_22_23, temp[23]);
    Ciphertext temp20_21_22_23 = temp20_21_22;
    evaluator.add_inplace(temp20_21_22_23, temp[23]);
    Ciphertext temp24_25 = temp[24];
    evaluator.add_inplace(temp24_25, temp[25]);
    Ciphertext temp24_26 = temp[24];
    evaluator.add_inplace(temp24_26, temp[26]);
    Ciphertext temp24_27 = temp[24];
    evaluator.add_inplace(temp24_27, temp[27]);
    Ciphertext temp25_26 = temp[25];
    evaluator.add_inplace(temp25_26, temp[26]);
    Ciphertext temp25_27 = temp[25];
    evaluator.add_inplace(temp25_27, temp[27]);
    Ciphertext temp26_27 = temp[26];
    evaluator.add_inplace(temp26_27, temp[27]);
    Ciphertext temp24_25_26 = temp24_25;
    evaluator.add_inplace(temp24_25_26, temp[26]);
    Ciphertext temp24_25_27 = temp24_25;
    evaluator.add_inplace(temp24_25_27, temp[27]);
    Ciphertext temp24_26_27 = temp24_26;
    evaluator.add_inplace(temp24_26_27, temp[27]);
    Ciphertext temp25_26_27 = temp25_26;
    evaluator.add_inplace(temp25_26_27, temp[27]);
    Ciphertext temp24_25_26_27 = temp24_25_26;
    evaluator.add_inplace(temp24_25_26_27, temp[27]);
    Ciphertext temp28_29 = temp[28];
    evaluator.add_inplace(temp28_29, temp[29]);
    Ciphertext temp28_30 = temp[28];
    evaluator.add_inplace(temp28_30, temp[30]);
    Ciphertext temp28_31 = temp[28];
    evaluator.add_inplace(temp28_31, temp[31]);
    Ciphertext temp29_30 = temp[29];
    evaluator.add_inplace(temp29_30, temp[30]);
    Ciphertext temp29_31 = temp[29];
    evaluator.add_inplace(temp29_31, temp[31]);
    Ciphertext temp30_31 = temp[30];
    evaluator.add_inplace(temp30_31, temp[31]);
    Ciphertext temp28_29_30 = temp28_29;
    evaluator.add_inplace(temp28_29_30, temp[30]);
    Ciphertext temp28_29_31 = temp28_29;
    evaluator.add_inplace(temp28_29_31, temp[31]);
    Ciphertext temp28_30_31 = temp28_30;
    evaluator.add_inplace(temp28_30_31, temp[31]);
    Ciphertext temp29_30_31 = temp29_30;
    evaluator.add_inplace(temp29_30_31, temp[31]);
    Ciphertext temp28_29_30_31 = temp28_29_30;
    evaluator.add_inplace(temp28_29_30_31, temp[31]);
    Ciphertext temp32_33 = temp[32];
    evaluator.add_inplace(temp32_33, temp[33]);
    Ciphertext temp32_34 = temp[32];
    evaluator.add_inplace(temp32_34, temp[34]);
    Ciphertext temp32_35 = temp[32];
    evaluator.add_inplace(temp32_35, temp[35]);
    Ciphertext temp33_34 = temp[33];
    evaluator.add_inplace(temp33_34, temp[34]);
    Ciphertext temp33_35 = temp[33];
    evaluator.add_inplace(temp33_35, temp[35]);
    Ciphertext temp34_35 = temp[34];
    evaluator.add_inplace(temp34_35, temp[35]);
    Ciphertext temp32_33_34 = temp32_33;
    evaluator.add_inplace(temp32_33_34, temp[34]);
    Ciphertext temp32_33_35 = temp32_33;
    evaluator.add_inplace(temp32_33_35, temp[35]);
    Ciphertext temp32_34_35 = temp32_34;
    evaluator.add_inplace(temp32_34_35, temp[35]);
    Ciphertext temp33_34_35 = temp33_34;
    evaluator.add_inplace(temp33_34_35, temp[35]);
    Ciphertext temp32_33_34_35 = temp32_33_34;
    evaluator.add_inplace(temp32_33_34_35, temp[35]);

    evaluator.add_inplace(eData[0], temp1_3);
    evaluator.add_inplace(eData[0], temp4_5_6_7);
    evaluator.add_inplace(eData[0], temp8_11);
    evaluator.add_inplace(eData[0], temp14_15);
    evaluator.add_inplace(eData[0], temp16_17_19);
    evaluator.add_inplace(eData[0], temp20_21_22);
    evaluator.add_inplace(eData[0], temp24_25);
    evaluator.add_inplace(eData[0], temp29_30_31);
    evaluator.add_inplace(eData[0], temp33_34_35);
    evaluator.add_inplace(eData[1], temp0_2_3);
    evaluator.add_inplace(eData[1], temp4_6);
    evaluator.add_inplace(eData[1], temp8_10);
    evaluator.add_inplace(eData[1], temp12_13_15);
    evaluator.add_inplace(eData[1], temp17_18);
    evaluator.add_inplace(eData[1], temp20_21_22_23);
    evaluator.add_inplace(eData[1], temp24_25_26);
    evaluator.add_inplace(eData[1], temp28_31);
    evaluator.add_inplace(eData[1], temp32_33_34);
    eData[2] = temp[1];
    evaluator.add_inplace(eData[2], temp4_5_7);
    evaluator.add_inplace(eData[2], temp8_9_10);
    evaluator.add_inplace(eData[2], temp12_14);
    evaluator.add_inplace(eData[2], temp16_17_18_19);
    evaluator.add_inplace(eData[2], temp20_21_23);
    evaluator.add_inplace(eData[2], temp25_26_27);
    evaluator.add_inplace(eData[2], temp28_29_30_31);
    evaluator.add_inplace(eData[2], temp32_33_35);
    evaluator.add_inplace(eData[3], temp0_1_2);
    evaluator.add_inplace(eData[3], temp4_6_7);
    evaluator.add_inplace(eData[3], temp8_9_10_11);
    evaluator.add_inplace(eData[3], temp[14]);
    evaluator.add_inplace(eData[3], temp17_18_19);
    evaluator.add_inplace(eData[3], temp20_22_23);
    evaluator.add_inplace(eData[3], temp24_25_27);
    evaluator.add_inplace(eData[3], temp[28]);
    evaluator.add_inplace(eData[3], temp32_33_34);
    evaluator.add_inplace(eData[4], temp0_1_3);
    evaluator.add_inplace(eData[4], temp5_6_7);
    evaluator.add_inplace(eData[4], temp9_11);
    evaluator.add_inplace(eData[4], temp13_15);
    evaluator.add_inplace(eData[4], temp16_18);
    evaluator.add_inplace(eData[4], temp20_21_23);
    evaluator.add_inplace(eData[4], temp24_25_26_27);
    evaluator.add_inplace(eData[4], temp28_29_31);
    evaluator.add_inplace(eData[4], temp34_35);
    eData[5] = temp0_2;
    evaluator.add_inplace(eData[5], temp4_7);
    evaluator.add_inplace(eData[5], temp8_10_11);
    evaluator.add_inplace(eData[5], temp12_13_15);
    evaluator.add_inplace(eData[5], temp17_19);
    evaluator.add_inplace(eData[5], temp20_21_22_23);
    evaluator.add_inplace(eData[5], temp24_26);
    evaluator.add_inplace(eData[5], temp28_29_30_31);
    evaluator.add_inplace(eData[5], temp32_33_34_35);
    evaluator.add_inplace(eData[6], temp0_1_3);
    evaluator.add_inplace(eData[6], temp4_5_7);
    evaluator.add_inplace(eData[6], temp9_10_11);
    evaluator.add_inplace(eData[6], temp12_13_14);
    evaluator.add_inplace(eData[6], temp[17]);
    evaluator.add_inplace(eData[6], temp20_21_22_23);
    evaluator.add_inplace(eData[6], temp25_26_27);
    evaluator.add_inplace(eData[6], temp28_30_31);
    evaluator.add_inplace(eData[6], temp[35]);
    evaluator.add_inplace(eData[7], temp1_2_3);
    evaluator.add_inplace(eData[7], temp4_6);
    evaluator.add_inplace(eData[7], temp8_9_10);
    evaluator.add_inplace(eData[7], temp12_14);
    evaluator.add_inplace(eData[7], temp16_18_19);
    evaluator.add_inplace(eData[7], temp21_23);
    evaluator.add_inplace(eData[7], temp24_26_27);
    evaluator.add_inplace(eData[7], temp28_29_30_31);
    evaluator.add_inplace(eData[7], temp32_34);
    eData[8] = temp0_1_2_3;
    evaluator.add_inplace(eData[8], temp5_7);
    evaluator.add_inplace(eData[8], temp10_11);
    evaluator.add_inplace(eData[8], temp13_14_15);
    evaluator.add_inplace(eData[8], temp16_18);
    evaluator.add_inplace(eData[8], temp20_22_23);
    evaluator.add_inplace(eData[8], temp24_25_26_27);
    evaluator.add_inplace(eData[8], temp29_31);
    evaluator.add_inplace(eData[8], temp32_33_34_35);
    evaluator.add_inplace(eData[9], temp2_3);
    evaluator.add_inplace(eData[9], temp4_6_7);
    evaluator.add_inplace(eData[9], temp8_10);
    evaluator.add_inplace(eData[9], temp12_13_14_15);
    evaluator.add_inplace(eData[9], temp16_17);
    evaluator.add_inplace(eData[9], temp20_23);
    evaluator.add_inplace(eData[9], temp24_25_26);
    evaluator.add_inplace(eData[9], temp28_29_30_31);
    evaluator.add_inplace(eData[9], temp33_34);
    evaluator.add_inplace(eData[10], temp[1]);
    evaluator.add_inplace(eData[10], temp4_5_6_7);
    evaluator.add_inplace(eData[10], temp9_11);
    evaluator.add_inplace(eData[10], temp12_13_15);
    evaluator.add_inplace(eData[10], temp17_19);
    evaluator.add_inplace(eData[10], temp21_22);
    evaluator.add_inplace(eData[10], temp24_26_27);
    evaluator.add_inplace(eData[10], temp29_30_31);
    evaluator.add_inplace(eData[10], temp32_33_34_35);
    eData[11] = temp0_1_2_3;
    evaluator.add_inplace(eData[11], temp4_5_6);
    evaluator.add_inplace(eData[11], temp8_10);
    evaluator.add_inplace(eData[11], temp13_14);
    evaluator.add_inplace(eData[11], temp16_17_18_19);
    evaluator.add_inplace(eData[11], temp21_23);
    evaluator.add_inplace(eData[11], temp25_26_27);
    evaluator.add_inplace(eData[11], temp28_29_30);
    evaluator.add_inplace(eData[11], temp32_34_35);
    evaluator.add_inplace(eData[12], temp0_1);
    evaluator.add_inplace(eData[12], temp5_6_7);
    evaluator.add_inplace(eData[12], temp9_10_11);
    evaluator.add_inplace(eData[12], temp13_15);
    evaluator.add_inplace(eData[12], temp16_17_18_19);
    evaluator.add_inplace(eData[12], temp20_23);
    evaluator.add_inplace(eData[12], temp26_27);
    evaluator.add_inplace(eData[12], temp28_29_31);
    evaluator.add_inplace(eData[12], temp32_33_34);
    evaluator.add_inplace(eData[13], temp0_1_2);
    evaluator.add_inplace(eData[13], temp4_7);
    evaluator.add_inplace(eData[13], temp8_9_10);
    evaluator.add_inplace(eData[13], temp12_14_15);
    evaluator.add_inplace(eData[13], temp16_18);
    evaluator.add_inplace(eData[13], temp20_22);
    evaluator.add_inplace(eData[13], temp24_25_27);
    evaluator.add_inplace(eData[13], temp29_30);
    evaluator.add_inplace(eData[13], temp32_33_34_35);
    eData[14] = temp1_2_3;
    evaluator.add_inplace(eData[14], temp4_5_6_7);
    evaluator.add_inplace(eData[14], temp8_9_11);
    evaluator.add_inplace(eData[14], temp[13]);
    evaluator.add_inplace(eData[14], temp16_17_19);
    evaluator.add_inplace(eData[14], temp20_21_22);
    evaluator.add_inplace(eData[14], temp24_26);
    evaluator.add_inplace(eData[14], temp28_29_30_31);
    evaluator.add_inplace(eData[14], temp32_33_35);
    evaluator.add_inplace(eData[15], temp0_1_3);
    evaluator.add_inplace(eData[15], temp[4]);
    evaluator.add_inplace(eData[15], temp8_9_10);
    evaluator.add_inplace(eData[15], temp12_13_14);
    evaluator.add_inplace(eData[15], temp16_18_19);
    evaluator.add_inplace(eData[15], temp20_21_22_23);
    evaluator.add_inplace(eData[15], temp[26]);
    evaluator.add_inplace(eData[15], temp29_30_31);
    evaluator.add_inplace(eData[15], temp32_34_35);
    evaluator.add_inplace(eData[16], temp0_1_2_3);
    evaluator.add_inplace(eData[16], temp4_5_7);
    evaluator.add_inplace(eData[16], temp10_11);
    evaluator.add_inplace(eData[16], temp12_13_15);
    evaluator.add_inplace(eData[16], temp17_18_19);
    evaluator.add_inplace(eData[16], temp21_23);
    evaluator.add_inplace(eData[16], temp25_27);
    evaluator.add_inplace(eData[16], temp28_30);
    evaluator.add_inplace(eData[16], temp32_33_35);
    eData[17] = temp0_2;
    evaluator.add_inplace(eData[17], temp4_5_6_7);
    evaluator.add_inplace(eData[17], temp8_9_10_11);
    evaluator.add_inplace(eData[17], temp12_14);
    evaluator.add_inplace(eData[17], temp16_19);
    evaluator.add_inplace(eData[17], temp20_22_23);
    evaluator.add_inplace(eData[17], temp24_25_27);
    evaluator.add_inplace(eData[17], temp29_31);
    evaluator.add_inplace(eData[17], temp32_33_34_35);
    evaluator.add_inplace(eData[18], temp1_2_3);
    evaluator.add_inplace(eData[18], temp4_6_7);
    evaluator.add_inplace(eData[18], temp[11]);
    evaluator.add_inplace(eData[18], temp12_13_15);
    evaluator.add_inplace(eData[18], temp16_17_19);
    evaluator.add_inplace(eData[18], temp21_22_23);
    evaluator.add_inplace(eData[18], temp24_25_26);
    evaluator.add_inplace(eData[18], temp[29]);
    evaluator.add_inplace(eData[18], temp32_33_34_35);
    evaluator.add_inplace(eData[19], temp0_2_3);
    evaluator.add_inplace(eData[19], temp4_5_6_7);
    evaluator.add_inplace(eData[19], temp8_10);
    evaluator.add_inplace(eData[19], temp13_14_15);
    evaluator.add_inplace(eData[19], temp16_18);
    evaluator.add_inplace(eData[19], temp20_21_22);
    evaluator.add_inplace(eData[19], temp24_26);
    evaluator.add_inplace(eData[19], temp28_30_31);
    evaluator.add_inplace(eData[19], temp33_35);
    eData[20] = temp0_1_2_3;
    evaluator.add_inplace(eData[20], temp5_7);
    evaluator.add_inplace(eData[20], temp8_9_10_11);
    evaluator.add_inplace(eData[20], temp12_13_14_15);
    evaluator.add_inplace(eData[20], temp17_19);
    evaluator.add_inplace(eData[20], temp22_23);
    evaluator.add_inplace(eData[20], temp25_26_27);
    evaluator.add_inplace(eData[20], temp28_30);
    evaluator.add_inplace(eData[20], temp32_34_35);
    evaluator.add_inplace(eData[21], temp0_1_2);
    evaluator.add_inplace(eData[21], temp4_5_6_7);
    evaluator.add_inplace(eData[21], temp9_10);
    evaluator.add_inplace(eData[21], temp14_15);
    evaluator.add_inplace(eData[21], temp16_18_19);
    evaluator.add_inplace(eData[21], temp20_22);
    evaluator.add_inplace(eData[21], temp24_25_26_27);
    evaluator.add_inplace(eData[21], temp28_29);
    evaluator.add_inplace(eData[21], temp32_35);
    evaluator.add_inplace(eData[22], temp0_2_3);
    evaluator.add_inplace(eData[22], temp5_6_7);
    evaluator.add_inplace(eData[22], temp8_9_10_11);
    evaluator.add_inplace(eData[22], temp[13]);
    evaluator.add_inplace(eData[22], temp16_17_18_19);
    evaluator.add_inplace(eData[22], temp21_23);
    evaluator.add_inplace(eData[22], temp24_25_27);
    evaluator.add_inplace(eData[22], temp29_31);
    evaluator.add_inplace(eData[22], temp33_34);
    eData[23] = temp1_2_3;
    evaluator.add_inplace(eData[23], temp4_5_6);
    evaluator.add_inplace(eData[23], temp8_10_11);
    evaluator.add_inplace(eData[23], temp12_13_14_15);
    evaluator.add_inplace(eData[23], temp16_17_18);
    evaluator.add_inplace(eData[23], temp20_22);
    evaluator.add_inplace(eData[23], temp25_26);
    evaluator.add_inplace(eData[23], temp28_29_30_31);
    evaluator.add_inplace(eData[23], temp33_35);
    evaluator.add_inplace(eData[24], temp2_3);
    evaluator.add_inplace(eData[24], temp4_5_7);
    evaluator.add_inplace(eData[24], temp8_9_10);
    evaluator.add_inplace(eData[24], temp12_13);
    evaluator.add_inplace(eData[24], temp17_18_19);
    evaluator.add_inplace(eData[24], temp21_22_23);
    evaluator.add_inplace(eData[24], temp25_27);
    evaluator.add_inplace(eData[24], temp28_29_30_31);
    evaluator.add_inplace(eData[24], temp32_35);
    evaluator.add_inplace(eData[25], temp0_1_3);
    evaluator.add_inplace(eData[25], temp5_6);
    evaluator.add_inplace(eData[25], temp8_9_10_11);
    evaluator.add_inplace(eData[25], temp12_13_14);
    evaluator.add_inplace(eData[25], temp16_19);
    evaluator.add_inplace(eData[25], temp20_21_22);
    evaluator.add_inplace(eData[25], temp24_26_27);
    evaluator.add_inplace(eData[25], temp28_30);
    evaluator.add_inplace(eData[25], temp32_34);
    eData[26] = temp0_2;
    evaluator.add_inplace(eData[26], temp4_5_6_7);
    evaluator.add_inplace(eData[26], temp8_9_11);
    evaluator.add_inplace(eData[26], temp13_14_15);
    evaluator.add_inplace(eData[26], temp16_17_18_19);
    evaluator.add_inplace(eData[26], temp20_21_23);
    evaluator.add_inplace(eData[26], temp[25]);
    evaluator.add_inplace(eData[26], temp28_29_31);
    evaluator.add_inplace(eData[26], temp32_33_34);
    evaluator.add_inplace(eData[27], temp[2]);
    evaluator.add_inplace(eData[27], temp5_6_7);
    evaluator.add_inplace(eData[27], temp8_10_11);
    evaluator.add_inplace(eData[27], temp12_13_15);
    evaluator.add_inplace(eData[27], temp[16]);
    evaluator.add_inplace(eData[27], temp20_21_22);
    evaluator.add_inplace(eData[27], temp24_25_26);
    evaluator.add_inplace(eData[27], temp28_30_31);
    evaluator.add_inplace(eData[27], temp32_33_34_35);
    evaluator.add_inplace(eData[28], temp1_3);
    evaluator.add_inplace(eData[28], temp4_6);
    evaluator.add_inplace(eData[28], temp8_9_11);
    evaluator.add_inplace(eData[28], temp12_13_14_15);
    evaluator.add_inplace(eData[28], temp16_17_19);
    evaluator.add_inplace(eData[28], temp22_23);
    evaluator.add_inplace(eData[28], temp24_25_27);
    evaluator.add_inplace(eData[28], temp29_30_31);
    evaluator.add_inplace(eData[28], temp33_35);
    eData[29] = temp0_1_3;
    evaluator.add_inplace(eData[29], temp5_7);
    evaluator.add_inplace(eData[29], temp8_9_10_11);
    evaluator.add_inplace(eData[29], temp12_14);
    evaluator.add_inplace(eData[29], temp16_17_18_19);
    evaluator.add_inplace(eData[29], temp20_21_22_23);
    evaluator.add_inplace(eData[29], temp24_26);
    evaluator.add_inplace(eData[29], temp28_31);
    evaluator.add_inplace(eData[29], temp32_34_35);
    evaluator.add_inplace(eData[30], temp0_1_2);
    evaluator.add_inplace(eData[30], temp[5]);
    evaluator.add_inplace(eData[30], temp8_9_10_11);
    evaluator.add_inplace(eData[30], temp13_14_15);
    evaluator.add_inplace(eData[30], temp16_18_19);
    evaluator.add_inplace(eData[30], temp[23]);
    evaluator.add_inplace(eData[30], temp24_25_27);
    evaluator.add_inplace(eData[30], temp28_29_31);
    evaluator.add_inplace(eData[30], temp33_34_35);
    evaluator.add_inplace(eData[31], temp0_2);
    evaluator.add_inplace(eData[31], temp4_6_7);
    evaluator.add_inplace(eData[31], temp9_11);
    evaluator.add_inplace(eData[31], temp12_14_15);
    evaluator.add_inplace(eData[31], temp16_17_18_19);
    evaluator.add_inplace(eData[31], temp20_22);
    evaluator.add_inplace(eData[31], temp25_26_27);
    evaluator.add_inplace(eData[31], temp28_30);
    evaluator.add_inplace(eData[31], temp32_33_34);
    eData[32] = temp1_2_3;
    evaluator.add_inplace(eData[32], temp4_6);
    evaluator.add_inplace(eData[32], temp8_10_11);
    evaluator.add_inplace(eData[32], temp12_13_14_15);
    evaluator.add_inplace(eData[32], temp17_19);
    evaluator.add_inplace(eData[32], temp20_21_22_23);
    evaluator.add_inplace(eData[32], temp24_25_26_27);
    evaluator.add_inplace(eData[32], temp29_31);
    evaluator.add_inplace(eData[32], temp34_35);
    evaluator.add_inplace(eData[33], temp0_1_2_3);
    evaluator.add_inplace(eData[33], temp4_5);
    evaluator.add_inplace(eData[33], temp8_11);
    evaluator.add_inplace(eData[33], temp12_13_14);
    evaluator.add_inplace(eData[33], temp16_17_18_19);
    evaluator.add_inplace(eData[33], temp21_22);
    evaluator.add_inplace(eData[33], temp26_27);
    evaluator.add_inplace(eData[33], temp28_30_31);
    evaluator.add_inplace(eData[33], temp32_34);
    evaluator.add_inplace(eData[34], temp0_1_3);
    evaluator.add_inplace(eData[34], temp5_7);
    evaluator.add_inplace(eData[34], temp9_10);
    evaluator.add_inplace(eData[34], temp12_14_15);
    evaluator.add_inplace(eData[34], temp17_18_19);
    evaluator.add_inplace(eData[34], temp20_21_22_23);
    evaluator.add_inplace(eData[34], temp[25]);
    evaluator.add_inplace(eData[34], temp28_29_30_31);
    evaluator.add_inplace(eData[34], temp33_35);
    eData[35] = temp1_2;
    evaluator.add_inplace(eData[35], temp4_5_6_7);
    evaluator.add_inplace(eData[35], temp9_11);
    evaluator.add_inplace(eData[35], temp13_14_15);
    evaluator.add_inplace(eData[35], temp16_17_18);
    evaluator.add_inplace(eData[35], temp20_22_23);
    evaluator.add_inplace(eData[35], temp24_25_26_27);
    evaluator.add_inplace(eData[35], temp28_29_30);
    evaluator.add_inplace(eData[35], temp32_34);
}
// Compute the constants for Sbox,(x0,x1,x2)——>(x0,x0*x2+x1,-x0*x1+x0*x2+x2)
void HE_Sbox(vector<Ciphertext> &eData, Evaluator &evaluator,RelinKeys &relin_keys)
{
    // (x0,x1,x2)——> (x0,x0*x2+x1,-x0*x1+x0*x2+x2)
    Ciphertext c01 = eData[1];
    Ciphertext c02 = eData[2];
    evaluator.multiply_inplace(c01, eData[0]);
    evaluator.relinearize_inplace(c01, relin_keys);
    // c01*=eData[0];
    evaluator.multiply_inplace(c02, eData[0]);
    evaluator.relinearize_inplace(c02, relin_keys);
    // c02*=eData[0];

    evaluator.add_inplace(eData[1], c02);
    evaluator.sub_inplace(eData[2], c01);
    evaluator.add_inplace(eData[2], c02);
    // omp_set_num_threads(12); // 设置线程数为12
    // #pragma omp parallel for
    for (long j = 3; j < BlockWords; j += 3)
    {
        c01 = eData[j + 1];
        evaluator.multiply_inplace(c01, eData[j]);
        evaluator.relinearize_inplace(c01, relin_keys);

        c02 = eData[j + 2];
        evaluator.multiply_inplace(c02, eData[j]);
        evaluator.relinearize_inplace(c02, relin_keys);

        evaluator.add_inplace(eData[j + 1], c02);
        evaluator.sub_inplace(eData[j + 2], c01);
        evaluator.add_inplace(eData[j + 2], c02);
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
    //=============客户端offline阶段================
    // 定义初始向量
    vector<long> IV(BlockWords);
    for (unsigned i = 0; i < BlockWords; i++)
    {
        IV[i] = i + 1;
    }


    // Generating Public Key and encrypting the symmetric key
    auto start = std::chrono::high_resolution_clock::now();
    EncryptionParameters parms(scheme_type::bfv);
    long poly_modulus_degree = phi_m;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    vector<int> bit_sizes = {40,40, 40 ,40, 40, 40}; // 你可以根据需要调整比特大小
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, bit_sizes));
    //parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(PlainMod);
    //SEALContext context(parms, true, sec_level_type::tc128);
    SEALContext context(parms, true, sec_level_type::none);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds_context = end - start;
    std::cout << "Context generation time: " << elapsed_seconds_context.count() << "s\n";
    print_parameters(context);

    auto start_PubKey = std::chrono::high_resolution_clock::now();
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    auto end_PubKey = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds_PubKey = end_PubKey - start_PubKey;
    std::cout << "PublicKey generation time: " << elapsed_seconds_PubKey.count() << "s\n";

    BatchEncoder batch_encoder(context);

    long nslots = batch_encoder.slot_count();

    auto &plain_modulus = parms.plain_modulus();

    long Qbits = context.key_context_data()->total_coeff_modulus_bit_count();

    long SecLevel = 80;

        // 生成随机对称密钥
    random_device rd;
    vector<long> SymKey(BlockWords);
    for (unsigned i = 0; i < BlockWords; i++)
    {
            SymKey[i] = rd() % PlainMod;
    }
    std::cout << "SymKey generated." << std::endl;

    auto start_keyEncryption = std::chrono::high_resolution_clock::now();
    vector<Ciphertext> encryptedSymKey(BlockWords);
    encryptSymKey(encryptedSymKey, SymKey, batch_encoder, encryptor, nslots);
    auto end_keyEncryption = std::chrono::high_resolution_clock::now();
    double keyEncryption = std::chrono::duration<double>(end_keyEncryption - start_keyEncryption).count();
    std::cout << "SymKey FHE time: " << keyEncryption << "s\n";
    //  解密验证
    if (symkeyflag)
    {
        if (!verify_encryptSymKey(encryptedSymKey, SymKey, batch_encoder, decryptor))
        {
            return 0;
        }
        std::cout << "Symmetric key encryption succeeded!" << std::endl;
    }
    // Generating key stream
    auto start_keyStream = std::chrono::high_resolution_clock::now();

    std::vector<long> NonceSet(PlainBlock);
    std::vector<long> RoundKeySet(KeyStreamWords * (Nr + 1));
    std::vector<long> KeyStream(PlainWords);
    RandomBit<BlockWords * randbits> randomBit(Nr);
    long X;
    std::vector<long> RoundKey(BlockWords);
    bool bit_array[randbits];
    long nonce;
    long block_num;
    long ir;
    std::vector<long> state(BlockWords); // 初始化 state
    Keccak_HashInstance shake128_2;
    std::cout << "Generating KeyStream..." << std::endl;
    uint64_t start_cycle = rdtsc();
    for (long counter = counter_begin; counter <= counter_end; counter++)
    {
        nonce = generate_secure_random_int(NonceSize);
        //randomBit.generate_Instance_all_new(nonce, counter);
        //auto &RanVecs = randomBit.roundconstants;
        random_init_shake(nonce, counter, shake128_2);
        block_num = counter - counter_begin;
        NonceSet[block_num] = nonce;
        // 逐轮进行加密
        for (unsigned r = 0; r <= Nr; r++)
        {
            // 计算RoundKey
            for (unsigned i = 0; i < BlockWords; ++i)
            {
                RoundKey[i] = (SymKey[i] * generate_random_field_element(shake128_2, false,max_prime_size, PlainMod))% PlainMod;
            }
            // 将RoundKey 复制到RoundKeySet
            // 测试使用
            if (deflag){
            memcpy(&RoundKeySet[KeyStreamWords * r + BlockWords * (block_num)], RoundKey.data(), BlockWords * sizeof(long));
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

                yusP.M36(state);  // 线性变换
                yusP.Sbox(state); // S盒
                for (unsigned i = 0; i < BlockWords; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
            }
            else
            {                      // 最后一轮
                yusP.M36(state); // 线性变换
                // yusP.Sbox_last(state); // S盒
                yusP.Sbox(state); // S盒
                for (unsigned i = 0; i < BlockWords; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
                yusP.M36(state); // 线性变换
                memcpy(&KeyStream[(block_num)*BlockPlainWords], state.data(), BlockPlainWords * sizeof(long));
            }
        }
    }
    uint64_t end_cycle = rdtsc();
    auto end_keyStream = std::chrono::high_resolution_clock::now();
    double Client_offtime = std::chrono::duration_cast<std::chrono::duration<double>>(end_keyStream - start_keyStream).count();
    std::cout << "Encryption offline total time: " << Client_offtime << "s\n";
    uint64_t cycle_count = end_cycle - start_cycle;
    std::cout << "Encryption offline total cycles: " << cycle_count << std::endl;
    // 将KeyStream同态密文写入文件
    // if (!writeEncryptedSymKey(encryptedSymKey, "Client_encryptedSymKey.bin"))
    // {
    //     return false;
    // }
    // std::cout << "Client_encryptedSymKey.bin has been written to file." << std::endl;

    //=============服务端offline阶段================
        for (int test = 0; test < 10; test++)
    {
    std::cout << "Generating XOF stream..." << std::endl;
    std::vector<vector<long>> Xset(NrBlockWords, vector<long>(nslots,0));
    long rB;
    // vec_long Xc;
    // Xc.SetLength(BlockWords);
    Keccak_HashInstance shake128_3;
    auto start_XOF = std::chrono::high_resolution_clock::now();
    for (long counter = counter_begin; counter <= counter_end; counter++)
    {
        block_num = counter - counter_begin;
        nonce = NonceSet[block_num];
        //randomBit.generate_Instance_all_new(nonce, counter);
        //auto &RanVecs = randomBit.roundconstants;
        random_init_shake(nonce, counter, shake128_3);
        // 逐轮进行加密
        for (unsigned r = 0; r <= Nr; r++)
        {
                        rB = r * BlockWords;
            // 计算 Xc
            for (unsigned i = 0; i < BlockWords; ++i)
            {
                Xset[rB + i][block_num] = generate_random_field_element(shake128_3, false,max_prime_size, PlainMod);
            }
        }
    }
    if (PlainBlock == 1){
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

    vector<Plaintext> encodedXset(NrBlockWords);

    auto start_Xset = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < NrBlockWords; i++)
        {
           batch_encoder.encode(Xset[i], encodedXset[i]);
        }
    auto end_Xset = std::chrono::high_resolution_clock::now();
    double Encode_time = std::chrono::duration<double>(end_Xset - start_Xset).count();
    std::cout << "encode time: " << Encode_time << "s\n";
    std::cout << "encode time/slots:" << Encode_time / (NrBlockWords + 2 * len3) * pow(10, 6) << "us\n";

    Ciphertext tmpCiphertext;
    int noise_budget = min_noise_budget(encryptedSymKey,decryptor);
    std::cout << "noise budget initially: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    // 计算 encryptedRoundKeySet
    long eRk_len = BlockWords * (Nr + 1);
    vector<Ciphertext> encryptedRoundKeySet(eRk_len, tmpCiphertext);
    // for (int i = 0; i < eRk_len; i++)
    // {
    //     encryptedRoundKeySet[i] = encryptedSymKey[i % BlockWords];
    // }
    int index = 0;
    for (int i = 0; i <= Nr; ++i)
    {
        for (int j = 0; j < BlockWords; ++j)
        {
            encryptedRoundKeySet[index++] = encryptedSymKey[j];
        }
    }
    auto start_RoundKeySet_FHE = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < eRk_len; i++)
    {
        evaluator.multiply_plain_inplace(encryptedRoundKeySet[i], encodedXset[i]);
    }
    auto end_RoundKeySet_FHE = std::chrono::high_resolution_clock::now();
    double RoundKey_time = std::chrono::duration<double>(end_RoundKeySet_FHE - start_RoundKeySet_FHE).count();
    std::cout << "RoundKeySet FHE succeeded! Time: " << RoundKey_time << "s\n";
    noise_budget = min_noise_budget(encryptedRoundKeySet,decryptor);
    std::cout << "noise budget after RoundKeySet FHE: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    // 使用 verifyDecryption 函数解密并验证 RoundKeySet
    // if (deflag)
    // {
    //     if (!verifyDecryption(encryptedRoundKeySet, RoundKeySet,batch_encoder, decryptor,BlockWords,PlainBlock, nslots, PlainMod));
    //     {
    //         std::cerr << "Decryption verification failed for RoundKeySet." << std::endl;
    //         return 0;
    //     }
    //     std::cout << "Decryption verification succeeded for RoundKeySet." << std::endl;
    // }

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
    vector<Ciphertext> encryptedKeyStream(encryptedRoundKeySet.begin(), encryptedRoundKeySet.begin() + BlockWords);
    std::cout << "whiteround start" << std::endl;
    start_add = std::chrono::high_resolution_clock::now();
    for (long i = 0; i < BlockWords; i++)
    { // encrypt the encoded key
        evaluator.add_plain_inplace(encryptedKeyStream[i], encodedXset[i]);
    }
    end_add = std::chrono::high_resolution_clock::now();
    Add_time += std::chrono::duration<double>(end_add - start_add).count();
    // 输出 Add_time
    std::cout << "whiteround time: " << Add_time << "s\n";
    noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
    std::cout << "noise budget after Whiteround: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    // 明文密钥流
    vector<long> KeyStream2(KeyStreamWords);
    if (deflag)
    {
        // 对IV和RoundKeySet进行异或
        for (long i = 0; i < KeyStreamWords; i++)
        {
            KeyStream2[i] = (IV[i % BlockWords] + RoundKeySet[i]) % PlainMod;
        }
        // 使用 verifyDecryption 函数解密并验证 KeyStream2
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor,BlockWords,PlainBlock, nslots, PlainMod))
        {
            std::cerr << "Decryption verification failed for KeyStream2." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for whiteround." << std::endl;
    }
    for (long r = 1; r < Nr; r++)
    {
        std::cout << "Round " << r << " start" << std::endl;
        start_linear = std::chrono::high_resolution_clock::now();
        // Linear Layer
        HE_M(encryptedKeyStream ,evaluator,relin_keys);
        end_linear = std::chrono::high_resolution_clock::now();
        Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
        linear_set[r - 1] = std::chrono::duration<double>(end_linear - start_linear).count();
        noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
        std::cout << "noise budget after linear: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            return 0;
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
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor,BlockWords,PlainBlock, nslots, PlainMod))
            {
                std::cerr << "Decryption verification failed for KeyStream2 Linear Layer." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream2 Linear Layer." << std::endl;
        }
        start_sbox = std::chrono::high_resolution_clock::now();
        // S Layer
        HE_Sbox(encryptedKeyStream,evaluator,relin_keys);
        end_sbox = std::chrono::high_resolution_clock::now();
        Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
        sbox_set[r - 1] = std::chrono::duration<double>(end_sbox - start_sbox).count();
        noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
        std::cout << "noise budget after sbox: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            return 0;
        }
        if (deflag)
        {
            yusP.Sbox(KeyStream2);
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor,BlockWords,PlainBlock, nslots, PlainMod))
            {
                std::cerr << "Decryption verification failed for KeyStream2 Sbox." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream2 Sbox." << std::endl;
        }
        // Round Key Addition
        start_add = std::chrono::high_resolution_clock::now();
        // omp_set_num_threads(12); // 设置线程数为12
        // #pragma omp parallel for
        for (long j = 0; j < BlockWords; j++)
        {
            evaluator.add_plain_inplace(encryptedKeyStream[j], encodedXset[r * BlockWords + j]);
        }
        end_add = std::chrono::high_resolution_clock::now();
        Add_time += std::chrono::duration<double>(end_add - start_add).count();
        noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
        std::cout << "noise budget after Add: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            return 0;
        }
        if (deflag)
        {
            for (long i = 0; i < KeyStreamWords; i++)
            {
                KeyStream2[i] = (KeyStream2[i] + RoundKeySet[r * KeyStreamWords + i]) % PlainMod;
            }
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor,BlockWords,PlainBlock, nslots, PlainMod))
            {
                std::cerr << "Decryption verification failed for KeyStream2 Round Key Addition." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream2 Round Key Addition." << std::endl;
        }
    }
// 最后一轮
#if (1)
    std::cout << "the last Round " << Nr << " start" << std::endl;
    start_linear = std::chrono::high_resolution_clock::now();
    // Linear Layer
    HE_M(encryptedKeyStream,evaluator,relin_keys);
    end_linear = std::chrono::high_resolution_clock::now();
    Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
    linear_set[Nr - 1] = std::chrono::duration<double>(end_linear - start_linear).count();
    noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
    std::cout << "noise budget after linear: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
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
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor,BlockWords,PlainBlock, nslots, PlainMod))
        {
            std::cerr << "Decryption verification failed for KeyStream2 Linear Layer." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream2 Linear Layer." << std::endl;
    }
    start_sbox = std::chrono::high_resolution_clock::now();
    // S Layer
    HE_Sbox(encryptedKeyStream,evaluator,relin_keys);
    // HE_Last_Sbox(encryptedKeyStream);
    end_sbox = std::chrono::high_resolution_clock::now();
    Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
    sbox_set[Nr - 1] = std::chrono::duration<double>(end_sbox - start_sbox).count();
    noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
    std::cout << "noise budget after sbox: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    if (deflag)
    {
        // yusP.Sbox_last(KeyStream2);
        yusP.Sbox(KeyStream2);
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor,BlockWords,PlainBlock, nslots, PlainMod))
        {
            std::cerr << "Decryption verification failed for KeyStream2 Sbox." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream2 Sbox." << std::endl;
    }
    // add
    start_add = std::chrono::high_resolution_clock::now();
    for (long j = 0; j < BlockWords; j++)
    {
        evaluator.add_plain_inplace(encryptedKeyStream[j], encodedXset[Nr * BlockWords + j]);
    }
    end_add = std::chrono::high_resolution_clock::now();
    Add_time += std::chrono::duration<double>(end_add - start_add).count();
    noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
    std::cout << "noise budget after Add: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    if (deflag)
    {
        for (long i = 0; i < KeyStreamWords; i++)
        {
            KeyStream2[i] = (KeyStream2[i] + RoundKeySet[Nr * KeyStreamWords + i]) % PlainMod;
        }
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor,BlockWords,PlainBlock, nslots, PlainMod))
        {
            std::cerr << "Decryption verification failed for KeyStream2 Round Key Addition." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream2 Round Key Addition." << std::endl;
    }
    start_linear = std::chrono::high_resolution_clock::now();
    // Linear Layer
    HE_M(encryptedKeyStream,evaluator,relin_keys);
    end_linear = std::chrono::high_resolution_clock::now();
    Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
    linear_set[Nr] = std::chrono::duration<double>(end_linear - start_linear).count();
    noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
    std::cout << "noise budget after linear: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
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
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor,BlockWords,PlainBlock, nslots, PlainMod))
        {
            std::cerr << "Decryption verification failed for KeyStream2 Linear Layer." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream2 Linear Layer." << std::endl;
    }
    //截断密钥流
     encryptedKeyStream.erase(encryptedKeyStream.begin() + BlockPlainWords, encryptedKeyStream.end());

#endif

    // 输出 XOF_time,Encode_time,Add_time、Sbox_time、Linear_time
    std::cout << "XOF time: " << XOF_time << "s\n";
    std::cout << "Encode time: " << Encode_time << "s\n";
    std::cout << "RoundKey time: " << RoundKey_time << "s\n";
    std::cout << "Add time: " << Add_time << "s\n";
    std::cout << "Sbox time: " << Sbox_time << "s\n";
    std::cout << "Linear time: " << Linear_time << "s\n";
    // 计算总时间
    double Server_offtime = XOF_time + Encode_time + RoundKey_time + Add_time + Sbox_time + Linear_time;
    std::cout << "Server offline total time: " << Server_offtime << "s\n";
    std::cout << "sbox_timeset: ";
    for (const auto &time : sbox_set) {
        std::cout << time << " ";
    }
    std::cout << endl;
    std::cout << "linear_timeset: ";
    for (const auto &time : linear_set) {
        std::cout << time << " ";
    }
    std::cout << endl;
    if (KeyStreamflag)
    {
        if (!verifyDecryption(encryptedKeyStream, KeyStream, batch_encoder, decryptor,BlockPlainWords,PlainBlock, nslots, PlainMod))
        {
            std::cerr << "Decryption verification failed for KeyStream." << std::endl;
            return 0;
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
    vector<vector<long>> CipherStream(BlockPlainWords, vector<long>(nslots,0));
    auto start_ClientOnline = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < PlainBlock; j++)
    {
        for (int i = 0; i < BlockPlainWords; i++)
        {
            CipherStream[i][j] = ((PlainStream[j * BlockPlainWords + i] - KeyStream[j * BlockPlainWords + i]) % PlainMod + PlainMod) % PlainMod;
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
    auto end_ClientOnline = std::chrono::high_resolution_clock::now();
    double Client_ontime = std::chrono::duration<double>(end_ClientOnline - start_ClientOnline).count();
    std::cout << "Client onine total time:" << std::chrono::duration<double>(end_ClientOnline - start_ClientOnline).count() << "s\n";
    double Client_totaltime = Client_offtime + Client_ontime;
    std::cout << "Client total time: " << Client_totaltime << "s\n";
    // 服务端在线
    // 同态加密
    vector<Ciphertext> encrypedPlainStream = encryptedKeyStream;
    // 对CipherStream进行编码
    vector<Plaintext> encodedCipherStream(BlockPlainWords);
    auto start_ServerOnline = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < BlockPlainWords; i++)
    {
        batch_encoder.encode(CipherStream[i], encodedCipherStream[i]);
    }
    for (int i = 0; i < BlockPlainWords; i++)
    {
        evaluator.add_plain_inplace(encrypedPlainStream[i],encodedCipherStream[i]);
    }
    auto end_ServerOnline = std::chrono::high_resolution_clock::now();
    double server_ontime = std::chrono::duration<double>(end_ServerOnline - start_ServerOnline).count();
    std::cout << "Server onine total time:" << server_ontime << "s\n";
    noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
    std::cout << "noise budget after online: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    // 同态解密验证
    // for (int i = 0; i < encryptedKeyStream.size(); i++)
    // {
    //     encryptedKeyStream[i].bringToSet(encryptedKeyStream[i].naturalPrimeSet());
    // }
    if (plainflag)
    {
        if (!verifyDecryption(encrypedPlainStream, PlainStream,batch_encoder, decryptor,BlockPlainWords,PlainBlock, nslots, PlainMod))
        {
            std::cerr << "Decryption verification failed for encrypedPlainStream." << std::endl;
            return 0;
        }
    }
    std::cout << "Decryption verification succeeded for encrypedPlainStream." << std::endl;
    // 计算吞吐量,KiB/min
    double Server_totaltime = Server_offtime + server_ontime;
    double throughput = (Plainbits * 60) / (pow(2, 13) * Server_totaltime);
    std::cout << "Server total time: " << Server_totaltime << "s\n";
    std::cout << "Throughput: " << throughput << "KiB/min\n";
    std::string dirPath = "../tests";
    std::string filePath;
    if (!fs::exists(dirPath))
    {
        filePath = "test_Yus_4.txt";
    }
    else
    {
        filePath = "../tests/test_Yus_4.txt";
    }
    std::ofstream outfile(filePath, std::ios::app);
    if (!outfile)
    {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return 0;
    }
    outfile << std::left << std::setw(3) << Nr
            << std::left << std::setw(12) << Para_p
            << std::left << std::setw(10) << nslots
            << std::left << std::setw(10) << PlainBlock
            << std::left << std::setw(6) << Qbits
            << std::fixed << std::setprecision(3)
            << std::left << std::setw(14) << SecLevel
            << std::left << std::setw(7) << XOF_time
            << std::left << std::setw(10) << Encode_time
            << std::left << std::setw(12) << RoundKey_time
            << std::left << std::setw(7) << Add_time
            << std::left << std::setw(8) << Sbox_time
            << std::left << std::setw(10) << Linear_time
            << std::left << std::setw(10) << Server_offtime
            << std::left << std::setw(10) << server_ontime
            << std::left << std::setw(10) << Server_totaltime
            << std::left << std::setw(20) << throughput
            << std::left << std::setw(10) << noise_budget
            << std::endl;
    outfile.close();
    std::cout << "test_Yus_4.txt updated." << std::endl;
    }
    return 0;
}