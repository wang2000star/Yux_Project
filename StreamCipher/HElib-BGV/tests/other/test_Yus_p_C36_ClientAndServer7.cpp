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
#include <NTL/GF2X.h>
#include <helib/helib.h>
#include <helib/ArgMap.h>
#include <helib/DoubleCRT.h>

#include "random_bit.hpp"
#include "Yus_p.hpp"
#include "tool.hpp"

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
constexpr long BlockByte = 36; // 分组字节长度
// ===============模式设置================
static bool Rkflag = 1;     // true/1表示乘法，false/0表示加法，指的是随机向量和密钥间的操作
static bool deflag = 0;     // true/1表示进行每一步解密验证，false/0表示不进行每一步解密验证
static bool ompflag = 0;    // true/1表示使用OpenMP并行编码，false/0表示不使用OpenMP并行编码
static bool symkeyflag = 0; // true/1表示对称密钥同态解密验证加密，false/0表示不验证
static bool plainflag = 0;  // true/1表示对称密文同态解密验证，false/0表示不验证
// 参数设置，paramMap[Nr-4][idx]
static constexpr unsigned Nr = 4; // 轮数
constexpr long idx = 2;
constexpr unsigned Sbox_depth = 1 * Nr; // S盒深度
// 当c=2时，Qbits=1.5*bits,当c=3时，Qbits=1.5*bits - 100
// 硬编码参数值
constexpr tuple<long, long, long, long> paramMap[5][8] = {
    {// Nr = 4
     // {p, log2(m), bits, c}
     {65537, 15, 230, 2},   // 0 *
     {163841, 15, 240, 2},  // 1 
     {65537, 14, 220, 2},   // 2 *
     {163841, 14, 230, 2},  // 3 
     {65537, 15, 350, 2},   // 4
     {163841, 15, 350, 2},  // 5
     {65537, 15, 400, 2},   // 6
     {163841, 15, 400, 2}}, // 7
    {
        // Nr = 5
        // {p, log2(m), bits, c}
        {65537, 15, 240, 2},  // 0 *
        {163841, 15, 280, 2}, // 1 
        {65537, 16, 300, 2},  // 2
        {65537, 16, 350, 2},  // 3
        {0, 0, 0, 0},         // 填充空位
        {0, 0, 0, 0},         // 填充空位
        {0, 0, 0, 0},         // 填充空位
        {0, 0, 0, 0}          // 填充空位
    },
    {
        // Nr = 6
        // {p, log2(m), bits, c}
        {65537, 16, 340, 2}, // 0 *
        {65537, 16, 360, 2}, // 
        {65537, 16, 370, 2},
        {0, 0, 0, 0}, // 填充空位
        {0, 0, 0, 0}, // 填充空位
        {0, 0, 0, 0}, // 填充空位
        {0, 0, 0, 0}, // 填充空位
        {0, 0, 0, 0}  // 填充空位
    },
    {
        // Nr = 7
        // {p, log2(m), bits, c}
        {65537, 16, 400, 2}, // 0 *
        {0, 0, 0, 0},        // 填充空位
        {0, 0, 0, 0},        // 填充空位
        {0, 0, 0, 0},        // 填充空位
        {0, 0, 0, 0},        // 填充空位
        {0, 0, 0, 0},        // 填充空位
        {0, 0, 0, 0},        // 填充空位
        {0, 0, 0, 0}         // 填充空位
    },
    {
        // Nr = 8
        // {p, log2(m), bits, c}
        {65537, 16, 430, 2}, // 0 *
        {786433, 17, 450, 2}, // slots太大，不建议使用
        {0, 0, 0, 0},         // 填充空位
        {0, 0, 0, 0},         // 填充空位
        {0, 0, 0, 0},         // 填充空位
        {0, 0, 0, 0},         // 填充空位
        {0, 0, 0, 0},         // 填充空位
        {0, 0, 0, 0}          // 填充空位
    }};
// p=k*m+1
//  2^10=1024,2^11=2048,2^12=4096,2^13=8192,2^14=16384,2^15=32768,2^16=65536,
// 当电脑内存有限时，log2Para_m太大，会导致内存不足，出现terminate called recursively错误，从而终止程序.
// 此外，即使正常运行，由于内存一直临界，会导致程序运行速度变慢，时间测量不准确。

constexpr long log2Para_m = get<1>(paramMap[Nr - 4][idx]) - 0;
constexpr long Para_p = get<0>(paramMap[Nr - 4][idx]);    // plaintext prime
constexpr long Para_m = 1 << log2Para_m;                  // cyclotomic polynomial
constexpr long phi_m = Para_m >> 1;                       // phi(m)=nlsots
constexpr long Para_bits = get<2>(paramMap[Nr - 4][idx]); // bits in the ciphertext modulus chain
constexpr long Para_c = get<3>(paramMap[Nr - 4][idx]);    // columns in the key-switching matrix

//!!!!!!!!!!!!!!!
constexpr unsigned PlainBlock = phi_m - 0; // 明文分组数,应该PlainBlock<=phi_m

// 计算 log2 的 constexpr 函数
constexpr unsigned int log2_constexpr(unsigned long long n, unsigned int p = 0)
{
    return (n <= 1) ? p : log2_constexpr(n / 2, p + 1);
}
constexpr long PlainMod = Para_p;                               // 明文模数
constexpr unsigned Bytebits = log2_constexpr(PlainMod - 1) + 1; // 字节比特长度=ceil(log2(PlainMod-1))

constexpr unsigned BlockSize = Bytebits * BlockByte; // 分组比特长度=BlockByte*Bytebits

static const long PlainByte = BlockByte * PlainBlock; // 明文字节长度
static const long Plainbits = Bytebits * PlainByte;   // 明文比特长度
static const long PlainSize = BlockSize * PlainBlock; // 明文比特长度

static const unsigned NonceSize = 32;                           // Nonce比特长度
static const long counter_begin = 0;                            // 计数器起始值
static const long counter_end = PlainBlock + counter_begin - 1; // 计数器结束值

YusP yusP(PlainMod); // 构建明文对称加密实例


void encodeTo36Ctxt(vector<ZZX> &encData, const vector<long> &data, const EncryptedArray &ea)
{
    long R = data.size() / PlainByte;
    long nCtxt = BlockByte * R;
    long data_size = data.size();
    long ea_size = ea.size();
    encData.resize(nCtxt);
#if (opmflag)
    omp_set_num_threads(12); // 设置线程数为12
#pragma omp parallel for
#endif
    for (long i = 0; i < BlockByte; i++)
    {
        vector<long> slots(ea_size, 0);
        for (long r = 0; r < R; r++)
        {
            for (long j = 0; j < PlainBlock; j++)
            {
                long byteIdx = j * BlockByte + i + r * PlainByte;
                slots[j] = data[byteIdx];
            }
            ea.encode(encData[r * BlockByte + i], slots);
        }
    }
}
// encodeTo36Ctxt对应的解码
void decodeTo36Ctxt(vector<long> &data, const vector<vector<long>> &encData,
                    const EncryptedArray &ea)
{
    long R = encData.size() / BlockByte;
    long data_size = R * PlainByte;
    data.resize(data_size);
    omp_set_num_threads(12); // 设置线程数为12
#pragma omp parallel for
    for (long j = 0; j < PlainBlock; j++)
    {
        for (long r = 0; r < R; r++)
        {
            for (long i = 0; i < BlockByte; i++)
            { // i is the ciphertext number
                // j is the block number in this ctxt
                long byteIdx = j * BlockByte + i + r * PlainByte;
                data[byteIdx] = encData[r * BlockByte + i][j];
            }
        }
    }
}
// 函数：解密并验证密文是否正确，需要解码
// 函数：解密并验证密文是否正确
bool verifyDecryption36(const std::vector<Ctxt> &encryptedVec, const vector<long> &originalVec, const SecKey &secretKey,
                        const EncryptedArray &ea)
{
    auto start_decrypt = std::chrono::steady_clock::now();
    vector<long> decryptedVec;
    std::vector<std::vector<long>> decryptedPolys(encryptedVec.size());
    omp_set_num_threads(12); // 设置线程数为12
#pragma omp parallel for
    for (std::size_t i = 0; i < encryptedVec.size(); ++i)
    {
        ea.decrypt(encryptedVec[i], secretKey, decryptedPolys[i]);
    }
    // 解码
    decodeTo36Ctxt(decryptedVec, decryptedPolys, ea);
    // 验证解密结果
    bool isDecryptedVecCorrect = std::equal(decryptedVec.begin(), decryptedVec.end(), originalVec.begin());
    auto end_decrypt = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification finished! Time: " << elapsed_seconds.count() << "s\n";
    // 如果解密结果不正确，输出第一个错误的位置
    if (!isDecryptedVecCorrect)
    {
        for (size_t i = 0; i < BlockByte; i++)
        {
            if (decryptedVec[i] != originalVec[i])
            {
                std::cout << "Error at position " << i << ": " << decryptedVec[i] << " != " << originalVec[i] << std::endl;
                // break;
            }
        }
    }
    return isDecryptedVecCorrect;
}
void encryptSymKey(vector<Ctxt> &encryptedSymKey, const vector<long> &SymKey, unique_ptr<PubKey> &pk, EncryptedArray &ea)
{
    long nslots = ea.size();
    // 加密
    encryptedSymKey.resize(BlockByte, Ctxt(*pk));
    for (long i = 0; i < BlockByte; i++)
    { // encrypt the encoded key
        vector<long> slotsData(nslots, SymKey[i]);
        ea.encrypt(encryptedSymKey[i], *pk, slotsData);
    }
}
bool verify_encryptSymKey(vector<Ctxt> &encryptedSymKey, const vector<long> &SymKey, const SecKey &secretKey, EncryptedArray &ea)
{
    auto start_decrypt = std::chrono::steady_clock::now();
    vector<long> decryptedSymKey(BlockByte);
    omp_set_num_threads(12); // 设置线程数为12
#pragma omp parallel for
    for (long i = 0; i < BlockByte; i++)
    { // encrypt the encoded key
        vector<long> slotsData;
        ea.decrypt(encryptedSymKey[i], secretKey, slotsData);
        decryptedSymKey[i] = slotsData[0];
    }
    bool isDecryptedSymKeyCorrect = std::equal(SymKey.begin(), SymKey.end(), decryptedSymKey.begin());
    auto end_decrypt = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification finished! Time: " << elapsed_seconds.count() << "s\n";
    return isDecryptedSymKeyCorrect;
}
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

    eData[0] += temp1_2_3;
    eData[0] += temp4_5;
    eData[0] += temp8_9_10_11;
    eData[0] += temp14_15;
    eData[0] += temp16_17_18;
    eData[0] += temp[22];
    eData[0] += temp24_26_27;
    eData[0] += temp28_29_30_31;
    eData[0] += temp[33];
    eData[1] = temp0_2_3;
    eData[1] += temp4_6_7;
    eData[1] += temp9_10_11;
    eData[1] += temp13_14_15;
    eData[1] += temp[17];
    eData[1] += temp20_21_22_23;
    eData[1] += temp24_25_26_27;
    eData[1] += temp29_30;
    eData[1] += temp32_33_35;
    eData[2] += temp1_3;
    eData[2] += temp[6];
    eData[2] += temp8_10_11;
    eData[2] += temp[12];
    eData[2] += temp16_17;
    eData[2] += temp20_23;
    eData[2] += temp26_27;
    eData[2] += temp29_31;
    eData[2] += temp[32];
    eData[3] = temp[1];
    eData[3] += temp5_6_7;
    eData[3] += temp9_10;
    eData[3] += temp12_14_15;
    eData[3] += temp17_19;
    eData[3] += temp20_21_23;
    eData[3] += temp25_26;
    eData[3] += temp28_29_30_31;
    eData[3] += temp34_35;
    eData[4] += temp[1];
    eData[4] += temp5_6_7;
    eData[4] += temp8_9;
    eData[4] += temp12_13_14_15;
    eData[4] += temp18_19;
    eData[4] += temp20_21_22;
    eData[4] += temp[26];
    eData[4] += temp28_30_31;
    eData[4] += temp32_33_34_35;
    eData[5] = temp0_1_3;
    eData[5] += temp4_6_7;
    eData[5] += temp8_10_11;
    eData[5] += temp13_14_15;
    eData[5] += temp17_18_19;
    eData[5] += temp[21];
    eData[5] += temp24_25_26_27;
    eData[5] += temp28_29_30_31;
    eData[5] += temp33_34;
    eData[6] += temp[0];
    eData[6] += temp5_7;
    eData[6] += temp[10];
    eData[6] += temp12_14_15;
    eData[6] += temp[16];
    eData[6] += temp20_21;
    eData[6] += temp24_27;
    eData[6] += temp30_31;
    eData[6] += temp33_35;
    eData[7] = temp2_3;
    eData[7] += temp[5];
    eData[7] += temp9_10_11;
    eData[7] += temp13_14;
    eData[7] += temp16_18_19;
    eData[7] += temp21_23;
    eData[7] += temp24_25_27;
    eData[7] += temp29_30;
    eData[7] += temp32_33_34_35;
    eData[8] += temp0_1_2_3;
    eData[8] += temp[5];
    eData[8] += temp9_10_11;
    eData[8] += temp12_13;
    eData[8] += temp16_17_18_19;
    eData[8] += temp22_23;
    eData[8] += temp24_25_26;
    eData[8] += temp[30];
    eData[8] += temp32_34_35;
    eData[9] = temp1_2;
    eData[9] += temp4_5_7;
    eData[9] += temp8_10_11;
    eData[9] += temp12_14_15;
    eData[9] += temp17_18_19;
    eData[9] += temp21_22_23;
    eData[9] += temp[25];
    eData[9] += temp28_29_30_31;
    eData[9] += temp32_33_34_35;
    eData[10] += temp1_3;
    eData[10] += temp[4];
    eData[10] += temp9_11;
    eData[10] += temp[14];
    eData[10] += temp16_18_19;
    eData[10] += temp[20];
    eData[10] += temp24_25;
    eData[10] += temp28_31;
    eData[10] += temp34_35;
    eData[11] = temp0_1_2_3;
    eData[11] += temp6_7;
    eData[11] += temp[9];
    eData[11] += temp13_14_15;
    eData[11] += temp17_18;
    eData[11] += temp20_22_23;
    eData[11] += temp25_27;
    eData[11] += temp28_29_31;
    eData[11] += temp33_34;
    eData[12] += temp0_2_3;
    eData[12] += temp4_5_6_7;
    eData[12] += temp[9];
    eData[12] += temp13_14_15;
    eData[12] += temp16_17;
    eData[12] += temp20_21_22_23;
    eData[12] += temp26_27;
    eData[12] += temp28_29_30;
    eData[12] += temp[34];
    eData[13] = temp0_1_2_3;
    eData[13] += temp5_6;
    eData[13] += temp8_9_11;
    eData[13] += temp12_14_15;
    eData[13] += temp16_18_19;
    eData[13] += temp21_22_23;
    eData[13] += temp25_26_27;
    eData[13] += temp[29];
    eData[13] += temp32_33_34_35;
    eData[14] += temp2_3;
    eData[14] += temp5_7;
    eData[14] += temp[8];
    eData[14] += temp13_15;
    eData[14] += temp[18];
    eData[14] += temp20_22_23;
    eData[14] += temp[24];
    eData[14] += temp28_29;
    eData[14] += temp32_35;
    eData[15] = temp1_2;
    eData[15] += temp4_5_6_7;
    eData[15] += temp10_11;
    eData[15] += temp[13];
    eData[15] += temp17_18_19;
    eData[15] += temp21_22;
    eData[15] += temp24_26_27;
    eData[15] += temp29_31;
    eData[15] += temp32_33_35;
    eData[16] += temp[2];
    eData[16] += temp4_6_7;
    eData[16] += temp8_9_10_11;
    eData[16] += temp[13];
    eData[16] += temp17_18_19;
    eData[16] += temp20_21;
    eData[16] += temp24_25_26_27;
    eData[16] += temp30_31;
    eData[16] += temp32_33_34;
    eData[17] = temp0_1_2_3;
    eData[17] += temp4_5_6_7;
    eData[17] += temp9_10;
    eData[17] += temp12_13_15;
    eData[17] += temp16_18_19;
    eData[17] += temp20_22_23;
    eData[17] += temp25_26_27;
    eData[17] += temp29_30_31;
    eData[17] += temp[33];
    eData[18] += temp0_3;
    eData[18] += temp6_7;
    eData[18] += temp9_11;
    eData[18] += temp[12];
    eData[18] += temp17_19;
    eData[18] += temp[22];
    eData[18] += temp24_26_27;
    eData[18] += temp[28];
    eData[18] += temp32_33;
    eData[19] = temp0_1_3;
    eData[19] += temp5_6;
    eData[19] += temp8_9_10_11;
    eData[19] += temp14_15;
    eData[19] += temp[17];
    eData[19] += temp21_22_23;
    eData[19] += temp25_26;
    eData[19] += temp28_30_31;
    eData[19] += temp33_35;
    eData[20] += temp0_1_2;
    eData[20] += temp[6];
    eData[20] += temp8_10_11;
    eData[20] += temp12_13_14_15;
    eData[20] += temp[17];
    eData[20] += temp21_22_23;
    eData[20] += temp24_25;
    eData[20] += temp28_29_30_31;
    eData[20] += temp34_35;
    eData[21] = temp[1];
    eData[21] += temp4_5_6_7;
    eData[21] += temp8_9_10_11;
    eData[21] += temp13_14;
    eData[21] += temp16_17_19;
    eData[21] += temp20_22_23;
    eData[21] += temp24_26_27;
    eData[21] += temp29_30_31;
    eData[21] += temp33_34_35;
    eData[22] += temp0_1;
    eData[22] += temp4_7;
    eData[22] += temp10_11;
    eData[22] += temp13_15;
    eData[22] += temp[16];
    eData[22] += temp21_23;
    eData[22] += temp[26];
    eData[22] += temp28_30_31;
    eData[22] += temp[32];
    eData[23] = temp1_3;
    eData[23] += temp4_5_7;
    eData[23] += temp9_10;
    eData[23] += temp12_13_14_15;
    eData[23] += temp18_19;
    eData[23] += temp[21];
    eData[23] += temp25_26_27;
    eData[23] += temp29_30;
    eData[23] += temp32_34_35;
    eData[24] += temp2_3;
    eData[24] += temp4_5_6;
    eData[24] += temp[10];
    eData[24] += temp12_14_15;
    eData[24] += temp16_17_18_19;
    eData[24] += temp[21];
    eData[24] += temp25_26_27;
    eData[24] += temp28_29;
    eData[24] += temp32_33_34_35;
    eData[25] = temp1_2_3;
    eData[25] += temp[5];
    eData[25] += temp8_9_10_11;
    eData[25] += temp12_13_14_15;
    eData[25] += temp17_18;
    eData[25] += temp20_21_23;
    eData[25] += temp24_26_27;
    eData[25] += temp28_30_31;
    eData[25] += temp33_34_35;
    eData[26] += temp[0];
    eData[26] += temp4_5;
    eData[26] += temp8_11;
    eData[26] += temp14_15;
    eData[26] += temp17_19;
    eData[26] += temp[20];
    eData[26] += temp25_27;
    eData[26] += temp[30];
    eData[26] += temp32_34_35;
    eData[27] = temp0_2_3;
    eData[27] += temp5_7;
    eData[27] += temp8_9_11;
    eData[27] += temp13_14;
    eData[27] += temp16_17_18_19;
    eData[27] += temp22_23;
    eData[27] += temp[25];
    eData[27] += temp29_30_31;
    eData[27] += temp33_34;
    eData[28] += temp0_1_2_3;
    eData[28] += temp6_7;
    eData[28] += temp8_9_10;
    eData[28] += temp[14];
    eData[28] += temp16_18_19;
    eData[28] += temp20_21_22_23;
    eData[28] += temp[25];
    eData[28] += temp29_30_31;
    eData[28] += temp32_33;
    eData[29] = temp1_2_3;
    eData[29] += temp5_6_7;
    eData[29] += temp[9];
    eData[29] += temp12_13_14_15;
    eData[29] += temp16_17_18_19;
    eData[29] += temp21_22;
    eData[29] += temp24_25_27;
    eData[29] += temp28_30_31;
    eData[29] += temp32_34_35;
    eData[30] += temp0_2_3;
    eData[30] += temp[4];
    eData[30] += temp8_9;
    eData[30] += temp12_15;
    eData[30] += temp18_19;
    eData[30] += temp21_23;
    eData[30] += temp[24];
    eData[30] += temp29_31;
    eData[30] += temp[34];
    eData[31] = temp1_2;
    eData[31] += temp4_6_7;
    eData[31] += temp9_11;
    eData[31] += temp12_13_15;
    eData[31] += temp17_18;
    eData[31] += temp20_21_22_23;
    eData[31] += temp26_27;
    eData[31] += temp[29];
    eData[31] += temp33_34_35;
    eData[32] += temp0_1;
    eData[32] += temp4_5_6_7;
    eData[32] += temp10_11;
    eData[32] += temp12_13_14;
    eData[32] += temp[18];
    eData[32] += temp20_22_23;
    eData[32] += temp24_25_26_27;
    eData[32] += temp[29];
    eData[32] += temp33_34_35;
    eData[33] = temp0_2_3;
    eData[33] += temp5_6_7;
    eData[33] += temp9_10_11;
    eData[33] += temp[13];
    eData[33] += temp16_17_18_19;
    eData[33] += temp20_21_22_23;
    eData[33] += temp25_26;
    eData[33] += temp28_29_31;
    eData[33] += temp32_34_35;
    eData[34] += temp[2];
    eData[34] += temp4_6_7;
    eData[34] += temp[8];
    eData[34] += temp12_13;
    eData[34] += temp16_19;
    eData[34] += temp22_23;
    eData[34] += temp25_27;
    eData[34] += temp[28];
    eData[34] += temp33_35;
    eData[35] = temp1_2_3;
    eData[35] += temp5_6;
    eData[35] += temp8_10_11;
    eData[35] += temp13_15;
    eData[35] += temp16_17_19;
    eData[35] += temp21_22;
    eData[35] += temp24_25_26_27;
    eData[35] += temp30_31;
    eData[35] += temp[33];

    eData[0] += temp0_1_3;
    eData[0] += temp[15];
    eData[0] += temp17_18;
    eData[0] += temp[22];
    eData[0] += temp24_27;
    eData[1] += temp6_7;
    eData[1] += temp[9];
    eData[1] += temp[13];
    eData[1] += temp20_21;
    eData[1] += temp25_27;
    eData[1] += temp[30];
    eData[1] += temp[33];
    eData[2] += temp8_10_11;
    eData[2] += temp[20];
    eData[2] += temp[32];
    eData[3] += temp[12];
    eData[3] += temp17_19;
    eData[3] += temp25_26;
    eData[3] += temp[31];
    eData[3] += temp[35];
    eData[4] += temp4_5_7;
    eData[4] += temp[19];
    eData[4] += temp21_22;
    eData[4] += temp[26];
    eData[4] += temp28_31;
    eData[5] += temp[1];
    eData[5] += temp10_11;
    eData[5] += temp[13];
    eData[5] += temp[17];
    eData[5] += temp24_25;
    eData[5] += temp29_31;
    eData[5] += temp[34];
    eData[6] += temp[0];
    eData[6] += temp12_14_15;
    eData[6] += temp[24];
    eData[7] += temp[3];
    eData[7] += temp[16];
    eData[7] += temp21_23;
    eData[7] += temp29_30;
    eData[7] += temp[35];
    eData[8] += temp8_9_11;
    eData[8] += temp[23];
    eData[8] += temp25_26;
    eData[8] += temp[30];
    eData[8] += temp32_35;
    eData[9] += temp[2];
    eData[9] += temp[5];
    eData[9] += temp14_15;
    eData[9] += temp[17];
    eData[9] += temp[21];
    eData[9] += temp28_29;
    eData[9] += temp33_35;
    eData[10] += temp[4];
    eData[10] += temp16_18_19;
    eData[10] += temp[28];
    eData[11] += temp[3];
    eData[11] += temp[7];
    eData[11] += temp[20];
    eData[11] += temp25_27;
    eData[11] += temp33_34;
    eData[12] += temp0_3;
    eData[12] += temp12_13_15;
    eData[12] += temp[27];
    eData[12] += temp29_30;
    eData[12] += temp[34];
    eData[13] += temp1_3;
    eData[13] += temp[6];
    eData[13] += temp[9];
    eData[13] += temp18_19;
    eData[13] += temp[21];
    eData[13] += temp[25];
    eData[13] += temp32_33;
    eData[14] += temp[8];
    eData[14] += temp20_22_23;
    eData[14] += temp[32];
    eData[15] += temp1_2;
    eData[15] += temp[7];
    eData[15] += temp[11];
    eData[15] += temp[24];
    eData[15] += temp29_31;
    eData[16] += temp[2];
    eData[16] += temp4_7;
    eData[16] += temp16_17_19;
    eData[16] += temp[31];
    eData[16] += temp33_34;
    eData[17] += temp0_1;
    eData[17] += temp5_7;
    eData[17] += temp[10];
    eData[17] += temp[13];
    eData[17] += temp22_23;
    eData[17] += temp[25];
    eData[17] += temp[29];
    eData[18] += temp[0];
    eData[18] += temp[12];
    eData[18] += temp24_26_27;
    eData[19] += temp5_6;
    eData[19] += temp[11];
    eData[19] += temp[15];
    eData[19] += temp[28];
    eData[19] += temp33_35;
    eData[20] += temp1_2;
    eData[20] += temp[6];
    eData[20] += temp8_11;
    eData[20] += temp20_21_23;
    eData[20] += temp[35];
    eData[21] += temp4_5;
    eData[21] += temp9_11;
    eData[21] += temp[14];
    eData[21] += temp[17];
    eData[21] += temp26_27;
    eData[21] += temp[29];
    eData[21] += temp[33];
    eData[22] += temp[4];
    eData[22] += temp[16];
    eData[22] += temp28_30_31;
    eData[23] += temp1_3;
    eData[23] += temp9_10;
    eData[23] += temp[15];
    eData[23] += temp[19];
    eData[23] += temp[32];
    eData[24] += temp[3];
    eData[24] += temp5_6;
    eData[24] += temp[10];
    eData[24] += temp12_15;
    eData[24] += temp24_25_27;
    eData[25] += temp[1];
    eData[25] += temp8_9;
    eData[25] += temp13_15;
    eData[25] += temp[18];
    eData[25] += temp[21];
    eData[25] += temp30_31;
    eData[25] += temp[33];
    eData[26] += temp[8];
    eData[26] += temp[20];
    eData[26] += temp32_34_35;
    eData[27] += temp[0];
    eData[27] += temp5_7;
    eData[27] += temp13_14;
    eData[27] += temp[19];
    eData[27] += temp[23];
    eData[28] += temp[7];
    eData[28] += temp9_10;
    eData[28] += temp[14];
    eData[28] += temp16_19;
    eData[28] += temp28_29_31;
    eData[29] += temp[1];
    eData[29] += temp[5];
    eData[29] += temp12_13;
    eData[29] += temp17_19;
    eData[29] += temp[22];
    eData[29] += temp[25];
    eData[29] += temp34_35;
    eData[30] += temp0_2_3;
    eData[30] += temp[12];
    eData[30] += temp[24];
    eData[31] += temp[4];
    eData[31] += temp9_11;
    eData[31] += temp17_18;
    eData[31] += temp[23];
    eData[31] += temp[27];
    eData[32] += temp[11];
    eData[32] += temp13_14;
    eData[32] += temp[18];
    eData[32] += temp20_23;
    eData[32] += temp32_33_35;
    eData[33] += temp2_3;
    eData[33] += temp[5];
    eData[33] += temp[9];
    eData[33] += temp16_17;
    eData[33] += temp21_23;
    eData[33] += temp[26];
    eData[33] += temp[29];
    eData[34] += temp4_6_7;
    eData[34] += temp[16];
    eData[34] += temp[28];
    eData[35] += temp[8];
    eData[35] += temp13_15;
    eData[35] += temp21_22;
    eData[35] += temp[27];
    eData[35] += temp[31];
}
// Compute Sbox
void HE_Sbox(vector<Ctxt> &eData)
{
    // for (int i = 0; i < eData.size(); i++)
    // {
    //     eData[i].cleanUp();
    // }
    // #pragma omp parallel for
    vector<Ctxt> temp = eData;
    // 一次迭代，(x0,x1,x2,x3)——> (x0, -x0*x1+x1+x0*x2+x0*x3, x1+x0*x2-x2+x0*x3, -x0*x1-x2-x0*x3+x3)
    Ctxt t01 = temp[0];
    t01.multiplyBy(temp[1]);
    Ctxt t02 = temp[0];
    t02.multiplyBy(temp[2]);
    Ctxt t03 = temp[0];
    t03.multiplyBy(temp[3]);
    t02 += t03;
    eData[1] -= t01;
    eData[1] += t02;
    t02 += temp[1];
    t02 -= temp[2];
    eData[2] = t02;
    eData[3] -= t01;
    eData[3] -= temp[2];
    eData[3] -= t03;
    for (int i = 4; i < BlockByte; i += 4)
    {
        t01 = temp[i];
        t01.multiplyBy(temp[i + 1]);
        t02 = temp[i];
        t02.multiplyBy(temp[i + 2]);
        t03 = temp[i];
        t03.multiplyBy(temp[i + 3]);
        t02 += t03;
        eData[i + 1] -= t01;
        eData[i + 1] += t02;
        t02 += temp[i + 1];
        t02 -= temp[i + 2];
        eData[i + 2] = t02;
        eData[i + 3] -= t01;
        eData[i + 3] -= temp[i + 2];
        eData[i + 3] -= t03;
    }
}
int main()
{
    std::cout << "Nr: " << Nr << std::endl;
    //=============客户端offline阶段================
    // 定义初始向量
    vector<long> IV(BlockByte);
    for (unsigned i = 0; i < BlockByte; i++)
    {
        IV[i] = i + 1;
    }
    // 生成随机对称密钥
    GF2X rnd;
    int Bytebitsdiv8 = ceil(Bytebits / 8);
    vector<uint8_t> SymKey0(Bytebitsdiv8 * BlockByte);
    random(rnd, 8 * SymKey0.size());
    BytesFromGF2X(SymKey0.data(), rnd, SymKey0.size());
    vector<long> SymKey(BlockByte);
    for (unsigned i = 0; i < BlockByte; i++)
    {
        SymKey[i] = 0;
        for (unsigned j = 0; j < Bytebitsdiv8; j++)
        {
            SymKey[i] += (SymKey0[Bytebitsdiv8 * i + j] << (8 * j));
        }
        SymKey[i] %= PlainMod;
    }

    std::cout << "SymKey generated." << std::endl;
    //========
    // Generating symmetric key and key stream
    auto start_keyStream = std::chrono::steady_clock::now();

    std::vector<long> NonceSet(PlainBlock);
    std::vector<long> Xset(PlainByte * (Nr + 1));
    std::vector<long> RoundKeySet(PlainByte * (Nr + 1));
    std::vector<long> KeyStream(PlainByte);

    long total_tasks = counter_end - counter_begin + 1;
    std::atomic<long> completed_tasks(0);
    // 定义进度打印的粒度，例如每完成 1% 的任务打印一次
    long progress_step = total_tasks / 100;
    if (progress_step == 0)
        progress_step = 1; // 防止除零
    RandomBit<BlockSize> randomBit(Nr);
    std::cout << "Generating KeyStream..." << std::endl;
    omp_set_num_threads(12); // 设置线程数为12
#pragma omp parallel for firstprivate(randomBit)
    for (long counter = counter_begin; counter <= counter_end; counter++)
    {
        long nonce = generate_secure_random_int(NonceSize);
        randomBit.generate_Instance_all_new(nonce, counter);
        auto &RanVecs = randomBit.roundconstants;
        // long nonce = counter;
        // std::vector<std::bitset<544>> RanVecs(Nr + 1);

        NonceSet[counter - counter_begin] = nonce;
        // 使用 std::array 代替 vector 并固定大小
        std::vector<long> state(BlockByte); // 初始化 state

        // 逐轮进行加密
        for (unsigned r = 0; r <= Nr; r++)
        {
            std::vector<long> X(BlockByte);
            std::vector<long> RoundKey(BlockByte);
            uint64_t temp;
            // 计算 X 和 RoundKey
            for (unsigned i = 0; i < BlockByte; ++i)
            {
                bool bit_array[Bytebits];
                for (unsigned j = 0; j < Bytebits; ++j)
                {
                    bit_array[j] = RanVecs[r][i * Bytebits + j];
                }
                BinStrToHex(bit_array, temp, Bytebits);
                // 强制转换为 long 类型
                X[i] = static_cast<long>(temp % PlainMod);
                if (Rkflag)
                {
                    RoundKey[i] = (SymKey[i] * X[i]) % PlainMod;
                }
                else
                {
                    RoundKey[i] = (SymKey[i] + X[i]) % PlainMod;
                }
            }

            // 将 X 和 RoundKey 复制到 Xset 和 RoundKeySet
            memcpy(&Xset[PlainByte * r + BlockByte * (counter - counter_begin)], X.data(), BlockByte * sizeof(long));
            memcpy(&RoundKeySet[PlainByte * r + BlockByte * (counter - counter_begin)], RoundKey.data(), BlockByte * sizeof(long));
            if (r == 0)
            { // 初始轮
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (RoundKey[i] + IV[i]) % PlainMod;
                }
            }
            else if (r < Nr)
            {                       // 常规轮
                yusP.M36_7(state);  // 线性变换
                yusP.Sbox_7(state); // S盒
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
            }
            else
            {                       // 最后一轮
                yusP.M36_7(state);  // 线性变换
                yusP.Sbox_7(state); // S盒
                yusP.M36_7(state);  // 线性变换
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
                memcpy(&KeyStream[(counter - counter_begin) * BlockByte], state.data(), BlockByte * sizeof(long));
            }
        }

        // 更新已完成的任务数
        long local_completed = ++completed_tasks;

        // 定期打印进度
        if (local_completed % progress_step == 0 || local_completed == total_tasks)
        {
#pragma omp critical
            {
                std::cout << "Progress: " << (local_completed * 100) / total_tasks << "% completed.\r" << std::flush;
            }
        }
    }
    std::cout << std::endl; // 完成后换行

    auto end_keyStream = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_keyStream = end_keyStream - start_keyStream;
    std::cout << "KeyStream Generation time: " << elapsed_seconds_keyStream.count() << "s\n";

    // Generating Public Key and encrypting the symmetric key

    long p = Para_p;
    long m = Para_m;
    long r = 1;
    long bits = Para_bits;
    long c = Para_c;
    long d = 1; // slots = phi(m)/d = phi(m) = 32768 = PlainBlock
    long k = 128;
    long s = 1;

    if (!m)
        m = FindM(k, bits, c, p, d, s, 0);
    auto start = std::chrono::steady_clock::now();

    // Context context = ContextBuilder<BGV>()
    //                                 .m(m)
    //                                 .p(p)
    //                                 .r(r)
    //                                 .bits(bits)
    //                                 .c(c)
    //                                 .buildPtr();
    shared_ptr<Context> context(ContextBuilder<BGV>()
                                    .m(m)
                                    .p(p)
                                    .r(r)
                                    .bits(bits)
                                    .c(c)
                                    .buildPtr());
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_context = end - start;
    std::cout << "Context generation time: " << elapsed_seconds_context.count() << "s\n";

    auto start_PubKey = std::chrono::steady_clock::now();
    SecKey secretKey(*context);
    secretKey.GenSecKey();
    unique_ptr<PubKey> publicKey = std::make_unique<helib::PubKey>(secretKey);
    helib::EncryptedArray ea(context->getEA());
    long nslots = ea.size();
    // if (nslots > PlainBlock)
    // {
    //     std::cerr << "nslots > PlainBlock" << std::endl;
    //     return false;
    // }
    std::cout << "p=" << p << std::endl;
    std::cout << "m=" << m << std::endl;
    std::cout << "nslots=" << nslots << std::endl;
    std::cout << "bits=" << bits << std::endl;
    std::cout << "c=" << c << std::endl;
    auto end_PubKey = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_PubKey = end_PubKey - start_PubKey;
    std::cout << "PublicKey generation time: " << elapsed_seconds_PubKey.count() << "s\n";
    // return 0;
    // 输出 context
    context->printout();
    std::cout << std::endl;
    long Qbits = context->bitSizeOfQ();
    double SecurityLevel = context->securityLevel();
    // 输出 securityLevel
    std::cout << "Security: " << context->securityLevel() << std::endl;

    auto start_keyEncryption = std::chrono::steady_clock::now();
    vector<Ctxt> encryptedSymKey;
    encryptSymKey(encryptedSymKey, SymKey, publicKey, ea);
    auto end_keyEncryption = std::chrono::steady_clock::now();
    double keyEncryption = std::chrono::duration<double>(end_keyEncryption - start_keyEncryption).count();
    std::cout << "SymKey FHE time: " << keyEncryption << "s\n";
    //  解密验证
    if (symkeyflag)
    {
        if (!verify_encryptSymKey(encryptedSymKey, SymKey, secretKey, ea))
        {
            return 0;
        }
        std::cout << "Symmetric key encryption succeeded!" << std::endl;
    }
    // 离线客户端时间=KeyStream Generation time+PublicKey generation and SymKey FHE time
    double total_time_off = elapsed_seconds_keyStream.count() + elapsed_seconds_PubKey.count() + elapsed_seconds_PubKey.count() + keyEncryption;
    std::cout << "Encryption offline total time: " << total_time_off << "s\n";

    //=============服务端offline阶段================
    // 计算 encryptedRoundKeySet
    vector<Ctxt> encryptedRoundKeySet;
    Ctxt tmpCtxt(*publicKey);
    long eRk_len = BlockByte * (Nr + 1);
    encryptedRoundKeySet.resize(eRk_len, tmpCtxt);
    for (int i = 0; i < eRk_len; i++)
    {
        encryptedRoundKeySet[i] = encryptedSymKey[i % BlockByte];
    }
    vector<ZZX> encodedXset;
    auto m1 = std::chrono::steady_clock::now();
    encodeTo36Ctxt(encodedXset, Xset, ea); // encode as HE plaintext
    auto m2 = std::chrono::steady_clock::now();
    double Encode_time = std::chrono::duration<double>(m2 - m1).count();
    std::cout << "encodeTo36Ctxt time: " << std::chrono::duration<double>(m2 - m1).count() << "s\n";

    auto start_RoundKeySet_FHE = std::chrono::steady_clock::now();
    if (Rkflag)
    {
        for (int i = 0; i < eRk_len; i++)
        {
            encryptedRoundKeySet[i].multByConstant(encodedXset[i]);
        }
    }
    else
    {
        for (int i = 0; i < eRk_len; i++)
        {
            encryptedRoundKeySet[i].addConstant(encodedXset[i]);
        }
    }
    auto end_RoundKeySet_FHE = std::chrono::steady_clock::now();
    double RoundKey_time = std::chrono::duration<double>(end_RoundKeySet_FHE - start_RoundKeySet_FHE).count();
    std::cout << "RoundKeySet FHE succeeded! Time: " << RoundKey_time << "s\n";
    // // 使用 verifyDecryption 函数解密并验证 RoundKeySet
    // if (deflag)
    // {
    //     if (!verifyDecryption36(encryptedRoundKeySet, RoundKeySet, secretKey, ea))
    //     {
    //         std::cerr << "Decryption verification failed for RoundKeySet." << std::endl;
    //         return 0;
    //     }
    //     std::cout << "Decryption verification succeeded for RoundKeySet." << std::endl;
    // }
    // 生成 encryptedKeyStream
    // 定义Add_time、Sbox_time、Linear_time
    double Sbox_time = 0, Linear_time = 0, Add_time = 0;

    vector<Ctxt> encryptedKeyStream;
    encryptedKeyStream.resize(BlockByte, tmpCtxt);
    std::copy(encryptedRoundKeySet.begin(), encryptedRoundKeySet.begin() + BlockByte, encryptedKeyStream.begin());

    vector<long> expandedIV(BlockByte * PlainBlock);
    for (long j = 0; j < PlainBlock; j++)
    {
        memcpy(&expandedIV[BlockByte * j], IV.data(), BlockByte * sizeof(long));
    }
    // 对expanded进行simd编码，这样会返回nRoundKeys个多项式数组即encoded，nRoundKeys=encoded.length()
    vector<ZZX> encoded_expandedIV;
    auto m3 = std::chrono::steady_clock::now();
    encodeTo36Ctxt(encoded_expandedIV, expandedIV, ea); // encode as HE plaintext
    auto m4 = std::chrono::steady_clock::now();
    std::cout << "encodeTo36Ctxt time: " << std::chrono::duration<double>(m4 - m3).count() << "s\n";

    std::cout << "whiteround start" << std::endl;
    std::cout << "noise budget: " << encryptedKeyStream[0].bitCapacity() << std::endl;
    auto start_roundkey = std::chrono::high_resolution_clock::now();
    for (long i = 0; i < BlockByte; i++)
    { // encrypt the encoded key
        encryptedKeyStream[i].addConstant(encoded_expandedIV[i]);
    }
    auto end_roundkey = std::chrono::high_resolution_clock::now();
    Add_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
    // 输出 Add_time
    std::cout << "whiteround time: " << Add_time << "s\n";
    std::cout << "noise budget: " << encryptedKeyStream[0].bitCapacity() << std::endl;
    // 测试
    // vector<Ctxt> test1 = encryptedKeyStream;
    // auto start_test1 = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < 1000; i++)
    // {
    //     test1[i % 36] = encryptedKeyStream[i % 36];
    // }
    // auto end_test1 = std::chrono::high_resolution_clock::now();
    // double test_time1 = std::chrono::duration<double>(end_test1 - start_test1).count();
    // vector<Ctxt> test2 = encryptedKeyStream;
    // auto start_test2 = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < 1000; i++)
    // {
    //     test2[i % 36] += encryptedKeyStream[i % 36];
    // }
    // auto end_test2 = std::chrono::high_resolution_clock::now();
    // double test_time2 = std::chrono::duration<double>(end_test2 - start_test2).count();
    // std::cout << "test1 time: " << test_time1 << "s\n";
    // std::cout << "test2 time: " << test_time2 << "s\n";

    // 明文密钥流
    vector<long> KeyStream2(PlainByte);
    if (deflag)
    {
        // 对IV和RoundKeySet进行异或
        for (long i = 0; i < PlainByte; i++)
        {
            KeyStream2[i] = (expandedIV[i] + RoundKeySet[i]) % PlainMod;
        }
        // 使用 verifyDecryption 函数解密并验证 KeyStream
        if (!verifyDecryption36(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for whiteround." << std::endl;
    }
    auto start_sbox = std::chrono::high_resolution_clock::now();
    auto start_linear = std::chrono::high_resolution_clock::now();
    auto end_sbox = std::chrono::high_resolution_clock::now();
    auto end_linear = std::chrono::high_resolution_clock::now();

    for (long r = 1; r < Nr; r++)
    {
        std::cout << "Round " << r << " start" << std::endl;
        start_linear = std::chrono::high_resolution_clock::now();
        // Linear Layer
        HE_M(encryptedKeyStream);
        end_linear = std::chrono::high_resolution_clock::now();
        Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
        if (deflag)
        {
            for (int i = 0; i < PlainBlock; i++)
            {
                vector<long> tmp(BlockByte);
                for (int j = 0; j < BlockByte; j++)
                {
                    tmp[j] = KeyStream2[i * BlockByte + j];
                }
                yusP.M36_7(tmp);
                for (int j = 0; j < BlockByte; j++)
                {
                    KeyStream2[i * BlockByte + j] = tmp[j];
                }
            }
            if (!verifyDecryption36(encryptedKeyStream, KeyStream2, secretKey, ea))
            {
                std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
        }
        start_sbox = std::chrono::high_resolution_clock::now();
        // S Layer
        HE_Sbox(encryptedKeyStream);
        end_sbox = std::chrono::high_resolution_clock::now();
        Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
        if (deflag)
        {
            yusP.Sbox_7(KeyStream2);
            if (!verifyDecryption36(encryptedKeyStream, KeyStream2, secretKey, ea))
            {
                std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
        }
        start_roundkey = std::chrono::high_resolution_clock::now();
        // omp_set_num_threads(12); // 设置线程数为12
        // #pragma omp parallel for
        for (long j = 0; j < BlockByte; j++)
        {
            encryptedKeyStream[j] += encryptedRoundKeySet[r * BlockByte + j];
        }
        end_roundkey = std::chrono::high_resolution_clock::now();
        Add_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
        if (deflag)
        {
            for (long i = 0; i < PlainByte; i++)
            {
                KeyStream2[i] = (KeyStream2[i] + RoundKeySet[r * PlainByte + i]) % PlainMod;
            }
            if (!verifyDecryption36(encryptedKeyStream, KeyStream2, secretKey, ea))
            {
                std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;
        }
        std::cout << "noise budget: " << encryptedKeyStream[0].bitCapacity() << std::endl;
    }
// 最后一轮
#if (1)
    std::cout << "Round " << Nr << " start" << std::endl;
    start_linear = std::chrono::high_resolution_clock::now();
    // Linear Layer
    HE_M(encryptedKeyStream);
    end_linear = std::chrono::high_resolution_clock::now();
    Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
    if (deflag)
    {
        for (int i = 0; i < PlainBlock; i++)
        {
            vector<long> tmp(BlockByte);
            for (int j = 0; j < BlockByte; j++)
            {
                tmp[j] = KeyStream2[i * BlockByte + j];
            }
            yusP.M36_7(tmp);
            for (int j = 0; j < BlockByte; j++)
            {
                KeyStream2[i * BlockByte + j] = tmp[j];
            }
        }
        if (!verifyDecryption36(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
    }
    start_sbox = std::chrono::high_resolution_clock::now();
    // S Layer
    HE_Sbox(encryptedKeyStream);
    end_sbox = std::chrono::high_resolution_clock::now();
    Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
    if (deflag)
    {
        yusP.Sbox_7(KeyStream2);
        if (!verifyDecryption36(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
    }
    start_linear = std::chrono::high_resolution_clock::now();
    // Linear Layer
    HE_M(encryptedKeyStream);
    end_linear = std::chrono::high_resolution_clock::now();
    Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
    if (deflag)
    {
        for (int i = 0; i < PlainBlock; i++)
        {
            vector<long> tmp(BlockByte);
            for (int j = 0; j < BlockByte; j++)
            {
                tmp[j] = KeyStream2[i * BlockByte + j];
            }
            yusP.M36_7(tmp);
            for (int j = 0; j < BlockByte; j++)
            {
                KeyStream2[i * BlockByte + j] = tmp[j];
            }
        }
        if (!verifyDecryption36(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
    }
    // add
    start_roundkey = std::chrono::high_resolution_clock::now();
    for (long j = 0; j < BlockByte; j++)
    {
        encryptedKeyStream[j] += encryptedRoundKeySet[Nr * BlockByte + j];
    }
    end_roundkey = std::chrono::high_resolution_clock::now();
    Add_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
    if (deflag)
    {
        for (long i = 0; i < PlainByte; i++)
        {
            KeyStream2[i] = (KeyStream2[i] + RoundKeySet[Nr * PlainByte + i]) % PlainMod;
        }
        if (!verifyDecryption36(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;
    }
#endif
    std::cout << "noise budget: " << encryptedKeyStream[0].bitCapacity() << std::endl;
    if (encryptedKeyStream[0].bitCapacity() <= 0)
    {
        std::cerr << "noise budget is 0" << std::endl;
        return 0;
    }
    // 输出 Add_time、Sbox_time、Linear_time
    std::cout << "RoundKey time: " << Add_time << "s\n";
    std::cout << "Sbox time: " << Sbox_time << "s\n";
    std::cout << "Linear Layer time: " << Linear_time << "s\n";
    // 计算总时间
    double total_time = Encode_time + RoundKey_time + Add_time + Sbox_time + Linear_time;
    std::cout << "Server offline total time: " << total_time << "s\n";
    // 计算吞吐量,KiB/min
    double throughput = (Plainbits * 60) / (pow(2, 13) * total_time);
    std::cout << "Throughput: " << throughput << "KiB/min\n";

    if (plainflag)
    {
        for (int i = 0; i < encryptedKeyStream.size(); i++)
        {
            encryptedKeyStream[i].bringToSet(encryptedKeyStream[i].naturalPrimeSet());
        }
        if (!verifyDecryption36(encryptedKeyStream, KeyStream, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream." << std::endl;
    }
    // 将total_time, throughput, Nr, p, nslots, bits, c, Add_time, Sbox_time, Linear_time, RoundKey_time写入文件test_Yus_p_C32_ClientAndServer2.txt,如果已存在则追加
    // 检查路径是否存在
    std::string dirPath = "../tests";
    std::string filePath;
    if (!fs::exists(dirPath))
    {
        filePath = "test_Yus_p_C36_ClientAndServer7.txt";
    }
    else
    {
        filePath = "../tests/test_Yus_p_C36_ClientAndServer7.txt";
    }
    std::ofstream outfile(filePath, std::ios::app);
    if (!outfile)
    {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return 0;
    }
    outfile << std::left << std::setw(3) << Nr
            << std::left << std::setw(10) << p
            << std::left << std::setw(10) << nslots
            << std::left << std::setw(5) << bits
            << std::left << std::setw(4) << c
            << std::left << std::setw(6) << Qbits
            << std::fixed << std::setprecision(3)
            << std::left << std::setw(14) << SecurityLevel
            << std::left << std::setw(10) << Encode_time
            << std::left << std::setw(12) << RoundKey_time
            << std::left << std::setw(7) << Add_time
            << std::left << std::setw(8) << Sbox_time
            << std::left << std::setw(10) << Linear_time
            << std::left << std::setw(9) << total_time
            << std::left << std::setw(10) << throughput
            << std::endl;
    outfile.close();
    std::cout << "test_Yus_p_C36_ClientAndServer7.txt updated." << std::endl;
    return 0;
}