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
constexpr unsigned BlockByte = 36; // 分组字节长度
constexpr unsigned TruncByte = 24; // 截断字节长度
constexpr double TruncRate = TruncByte / (double)BlockByte;
// ===============模式设置================
constexpr bool Rkflag = 1;     // true/1表示乘法，false/0表示加法，指的是随机向量和密钥间的操作
constexpr bool deflag = 0;     // true/1表示进行每一步解密验证，false/0表示不进行每一步解密验证
constexpr bool ompflag = 0;    // true/1表示使用OpenMP并行编码，false/0表示不使用OpenMP并行编码
constexpr bool symkeyflag = 0; // true/1表示对称密钥同态解密验证加密，false/0表示不验证
constexpr bool plainflag = 1;  // true/1表示对称密文同态解密验证，false/0表示不验证
// 参数设置，paramMap[Nr-4][idx]
constexpr unsigned Nr = 4; // 轮数
constexpr long idx = 0;
constexpr unsigned Sbox_depth = 1 * Nr; // S盒深度
// 当c=2时，Qbits=1.5*bits,当c=3时，Qbits=1.5*bits - 100
// 硬编码参数值
constexpr tuple<long, long, long, long> paramMap[5][8] = {
    {// Nr = 4
     // {p, log2(m), bits, c}
     {65537, 15, 200, 3},   // 0 *
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
        {65537, 15, 236, 2},  // 0 *
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
        {65537, 15, 300, 3}, // 0 *
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
        {65537, 16, 550, 2}, // 0 *
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
        {65537, 16, 430, 2},  // 0 *
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
constexpr long phi_m = Para_m >> 1;                       // phi(m)
constexpr long Para_bits = get<2>(paramMap[Nr - 4][idx]); // bits in the ciphertext modulus chain
constexpr long Para_c = get<3>(paramMap[Nr - 4][idx]);    // columns in the key-switching matrix
constexpr long Para_r = 1;                                // Lifting [defualt = 1]
//!!!!!!!!!!!!!!!
constexpr long nslots = phi_m;             // 槽数
constexpr unsigned PlainBlock = phi_m - 0; // 明文分组数,应该PlainBlock<=phi_m
constexpr unsigned len3 = BlockByte / 3;
// 计算 log2 的 constexpr 函数
constexpr unsigned int log2_constexpr(unsigned long long n, unsigned int p = 0)
{
    return (n <= 1) ? p : log2_constexpr(n / 2, p + 1);
}
constexpr long PlainMod = Para_p;                               // 明文模数
constexpr unsigned Bytebits = log2_constexpr(PlainMod - 1) + 1; // 字节比特长度=ceil(log2(PlainMod-1))
constexpr unsigned randbits = Bytebits - 1;
constexpr unsigned BlockSize = Bytebits * BlockByte;    // 分组比特长度=BlockByte*Bytebits
constexpr unsigned NrBlockByte = BlockByte * (Nr + 1);  // Nr轮分组密钥字节长度
constexpr long PlainByte = BlockByte * PlainBlock;      // 明文字节长度
constexpr long TruncPlainByte = TruncByte * PlainBlock; // 截断后的明文字节长度
constexpr long Plainbits = Bytebits * PlainByte;        // 明文比特长度

constexpr unsigned NonceSize = 32;                           // Nonce比特长度
constexpr long counter_begin = 0;                            // 计数器起始值
constexpr long counter_end = PlainBlock + counter_begin - 1; // 计数器结束值

YusP yusP(PlainMod); // 构建明文对称加密实例

// 函数：对多项式的每个系数乘以整数 a 并取模 c
zzX multiplyAndMod(const zzX &a, long b)
{
    zzX res;
    res.SetLength(lsize(a));
    int len = lsize(a);
    for (long i = 0; i < len; ++i)
    {
        res[i] = (a[i] * b) % Para_p;
    }
    return res;
}
// void vector_addition_avx2(const long* __restrict a, const long* __restrict b, long* __restrict result, size_t size){
//     // 检查size是否为8的倍数，确保可以正确处理AVX2的256位寄存器
//     assert(size % 8 == 0);
//     __m256i va, vb, vr;
//     for (size_t i = 0; i < size; i += 8)
//     {
//         // 加载8个long整数到AVX寄存器
//         va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + i));
//         vb = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b + i));
//         // 使用AVX2指令进行向量加法
//         vr = _mm256_add_epi32(va, vb);
//         // 存储结果回内存
//         _mm256_storeu_si256(reinterpret_cast<__m256i*>(result + i), vr);
//     }
// }
int min_noise_budget(vector<Ctxt> &eData)
{
    int min_noise = 10000;
    int noise;
    for (int i = 0; i < eData.size(); i++)
    {
        noise = eData[i].bitCapacity();
        if (noise < min_noise)
        {
            min_noise = noise;
        }
    }
    return min_noise;
}
bool writeEncryptedSymKey(const vector<Ctxt> &encryptedSymKey, const std::string &filename)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open())
    {
        std::cerr << "Failed to open " << filename << " for writing" << std::endl;
        return false;
    }

    for (const auto &ctxt : encryptedSymKey)
    {
        ctxt.writeTo(out);
    }

    out.close();

    return true;
}
void decodeToCtxt(vector<long> &data, const vector<vec_long> &encData, const long CtxtByte)
{
    long R = encData.size() / CtxtByte;
    long AllByte = CtxtByte * PlainBlock;
    long data_size = R * AllByte;
    data.resize(data_size);
    long byteIdx;
    long rAB;
    long jCB;
    long rCB;
    //    omp_set_num_threads(16); // 设置线程数为16
    // #pragma omp parallel for
    for (long j = 0; j < nslots; j++)
    {
        jCB = j * CtxtByte;
        for (long r = 0; r < R; r++)
        {
            rAB = r * AllByte;
            rCB = r * CtxtByte;
            for (long i = 0; i < CtxtByte; i++)
            { // i is the ciphertext number
                // j is the block number in this ctxt
                byteIdx = jCB + i + rAB;
                data[byteIdx] = encData[rCB + i][j];
            }
        }
    }
}
// 函数：解密并验证密文是否正确，需要解码
bool verifyDecryption(const std::vector<Ctxt> &encryptedVec, const vector<long> &originalVec, const SecKey &secretKey,
                      const helib::Cmodulus &cmodulus, const long CtxtByte)
{
    vector tempVec = originalVec;
    for (int i = 0; i < originalVec.size(); i++)
    {
        tempVec[i] = (tempVec[i] + Para_p) % Para_p;
    }
    int size = encryptedVec.size();
    auto start_decrypt = std::chrono::high_resolution_clock::now();
    vector<long> decryptedVec = originalVec;
    std::vector<vec_long> decryptedPolys(size);
    vector<ZZX> Polys(size);
    //     omp_set_num_threads(16); // 设置线程数为16
    // #pragma omp parallel for
    for (std::size_t i = 0; i < size; ++i)
    {
        secretKey.Decrypt(Polys[i], encryptedVec[i]);
        cmodulus.FFT(decryptedPolys[i], Polys[i]);
    }
    decodeToCtxt(decryptedVec, decryptedPolys, CtxtByte);
    // 验证解密结果
    bool isDecryptedVecCorrect = std::equal(decryptedVec.begin(), decryptedVec.end(), tempVec.begin());
    auto end_decrypt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification finished! Time: " << elapsed_seconds.count() << "s\n";
    // 如果解密结果不正确，输出错误的位置
    if (!isDecryptedVecCorrect)
    {
        for (size_t i = 0; i < CtxtByte; i++)
        {
            if (decryptedVec[i] != tempVec[i])
            {
                std::cout << "Error at position " << i << ": " << decryptedVec[i] << " != " << originalVec[i] << std::endl;
                // break;
            }
        }
    }
    return isDecryptedVecCorrect;
}
void encryptSymKey(vector<Ctxt> &encryptedSymKey, vector<Ctxt> &encryptedSymKey01, vector<Ctxt> &encryptedSymKey02,
                   const vector<long> &SymKey, unique_ptr<PubKey> &pk, const helib::Cmodulus &cmodulus)
{
    vec_long slotsData;
    slotsData.SetLength(nslots);
    encryptedSymKey.resize(BlockByte, Ctxt(*pk));
    zz_pX encodedData;
    ZZX encodedData2;
    for (long i = 0; i < BlockByte; i++)
    { // encrypt the encoded key
        for (long j = 0; j < nslots; j++)
        {
            slotsData[j] = SymKey[i];
        }
        cmodulus.iFFT(encodedData, slotsData);
        conv(encodedData2, encodedData);
        pk->Encrypt(encryptedSymKey[i], encodedData2);
    }
    encryptedSymKey01.resize(len3, Ctxt(*pk));
    encryptedSymKey02.resize(len3, Ctxt(*pk));
    vec_long slotsData01;
    slotsData01.SetLength(nslots);
    vec_long slotsData02;
    slotsData02.SetLength(nslots);
    long k01, k02;
    zz_pX encodedData01;
    ZZX encodedData01_2;
    zz_pX encodedData02;
    ZZX encodedData02_2;
    for (long i = 0; i < BlockByte; i += 3)
    { // encrypt the encoded key01,k02
        k01 = (SymKey[i] * SymKey[i + 1]) % PlainMod;
        k02 = (SymKey[i] * SymKey[i + 2]) % PlainMod;
        for (long j = 0; j < nslots; j++)
        {
            slotsData01[j] = k01;
            slotsData02[j] = k02;
        }
        cmodulus.iFFT(encodedData01, slotsData01);
        conv(encodedData01_2, encodedData01);
        pk->Encrypt(encryptedSymKey01[i / 3], encodedData01_2);
        cmodulus.iFFT(encodedData02, slotsData02);
        conv(encodedData02_2, encodedData02);
        pk->Encrypt(encryptedSymKey02[i / 3], encodedData02_2);
    }
}
bool verify_encryptSymKey(vector<Ctxt> &encryptedSymKey, const vector<long> &SymKey, const SecKey &secretKey,
                          const helib::Cmodulus &cmodulus)
{
    auto start_decrypt = std::chrono::high_resolution_clock::now();
    vector<long> decryptedSymKey(BlockByte);
    ZZX encodedSymKey;
    vec_long slotsData;
    //     omp_set_num_threads(16); // 设置线程数为16
    // #pragma omp parallel for
    for (long i = 0; i < BlockByte; i++)
    { // encrypt the encoded key
        secretKey.Decrypt(encodedSymKey, encryptedSymKey[i]);
        cmodulus.FFT(slotsData, encodedSymKey);
        decryptedSymKey[i] = slotsData[0];
    }
    bool isDecryptedSymKeyCorrect = std::equal(SymKey.begin(), SymKey.end(), decryptedSymKey.begin());
    auto end_decrypt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification finished! Time: " << elapsed_seconds.count() << "s\n";
    return isDecryptedSymKeyCorrect;
}
// Linear transformation
void HE_M2(vector<Ctxt> &eData)
{
    vector<Ctxt> temp = eData;
    // 0,1,2,3
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
    // 4,5,6,7
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
    // 8,9,10,11
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
    // 12,13,14,15
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
    // 16,17,18,19
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
    // 20,21,22,23
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
    // 24,25,26,27
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
    // 28,29,30,31
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
    // 32,33,34,35
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

    eData[0] = temp1_2_3;
    eData[0] += temp4_5_6_7;
    eData[0] += temp9_11;
    eData[0] += temp[13];
    eData[0] += temp16_17_19;
    eData[0] += temp20_22;
    eData[0] += temp25_26_27;
    eData[0] += temp29_30;
    eData[0] += temp33_35;

    eData[1] += temp[3];
    eData[1] += temp4_6_7;
    eData[1] += temp8_9_10;
    eData[1] += temp[12];
    eData[1] += temp17_18_19;
    eData[1] += temp21_22;
    eData[1] += temp24_27;
    eData[1] += temp[31];
    eData[1] += temp33_34_35;

    eData[2] = temp0_3;
    eData[2] += temp[7];
    eData[2] += temp9_10_11;
    eData[2] += temp13_14_15;
    eData[2] += temp17_18;
    eData[2] += temp21_22_23;
    eData[2] += temp[26];
    eData[2] += temp[28];
    eData[2] += temp32_33_34;

    eData[3] = temp0_2;
    eData[3] += temp4_5_6_7;
    eData[3] += temp8_9_10;
    eData[3] += temp12_14;
    eData[3] += temp16_19;
    eData[3] += temp20_22_23;
    eData[3] += temp[25];
    eData[3] += temp28_29_30;
    eData[3] += temp32_33;

    eData[4] += temp0_1_2;
    eData[4] += temp6_7;
    eData[4] += temp9_10_11;
    eData[4] += temp12_13_15;
    eData[4] += temp20_21_22;
    eData[4] += temp24_25_27;
    eData[4] += temp[30];
    eData[4] += temp[34];

    eData[5] = temp0_1_3;
    eData[5] += temp[6];
    eData[5] += temp[10];
    eData[5] += temp12_13_14;
    eData[5] += temp16_17_18;
    eData[5] += temp20_21;
    eData[5] += temp24_25_26;
    eData[5] += temp29_31;
    eData[5] += temp[35];

    eData[6] = temp0_3;
    eData[6] += temp5_7;
    eData[6] += temp8_9_10_11;
    eData[6] += temp12_13_15;
    eData[6] += temp17_19;
    eData[6] += temp22_23;
    eData[6] += temp25_26;
    eData[6] += temp28_31;
    eData[6] += temp32_33_35;

    eData[7] += temp1_3;
    eData[7] += temp4_5;
    eData[7] += temp9_10;
    eData[7] += temp12_13_14_15;
    eData[7] += temp16_18;
    eData[7] += temp[23];
    eData[7] += temp24_25_27;
    eData[7] += temp28_30;
    eData[7] += temp[33];

    eData[8] = temp2_3;
    eData[8] += temp4_6;
    eData[8] += temp[9];
    eData[8] += temp13_15;
    eData[8] += temp16_17_19;
    eData[8] += temp20_21_23;
    eData[8] += temp24_27;
    eData[8] += temp28_29;
    eData[8] += temp32_34;

    eData[9] = temp0_2_3;
    eData[9] += temp[6];
    eData[9] += temp8_10_11;
    eData[9] += temp12_13_14_15;
    eData[9] += temp16_18;
    eData[9] += temp20_22;
    eData[9] += temp25_26;
    eData[9] += temp28_29_31;
    eData[9] += temp34_35;

    eData[10] += temp[0];
    eData[10] += temp4_6_7;
    eData[10] += temp[8];
    eData[10] += temp12_13_15;
    eData[10] += temp16_17_18_19;
    eData[10] += temp[21];
    eData[10] += temp26_27;
    eData[10] += temp28_30_31;
    eData[10] += temp[33];

    eData[11] = temp[1];
    eData[11] += temp5_6_7;
    eData[11] += temp[9];
    eData[11] += temp[12];
    eData[11] += temp16_18_19;
    eData[11] += temp20_22_23;
    eData[11] += temp24_26_27;
    eData[11] += temp30_31;
    eData[11] += temp32_35;

    eData[12] = temp1_2_3;
    eData[12] += temp5_6;
    eData[12] += temp9_11;
    eData[12] += temp13_14_15;
    eData[12] += temp16_17_18_19;
    eData[12] += temp21_23;
    eData[12] += temp[25];
    eData[12] += temp28_29_31;
    eData[12] += temp32_34;

    eData[13] += temp0_3;
    eData[13] += temp[7];
    eData[13] += temp9_10_11;
    eData[13] += temp[15];
    eData[13] += temp16_18_19;
    eData[13] += temp20_21_22;
    eData[13] += temp[24];
    eData[13] += temp29_30_31;
    eData[13] += temp33_34;

    eData[14] = temp[2];
    eData[14] += temp[4];
    eData[14] += temp8_9_10;
    eData[14] += temp12_15;
    eData[14] += temp[19];
    eData[14] += temp21_22_23;
    eData[14] += temp25_26_27;
    eData[14] += temp29_30;
    eData[14] += temp33_34_35;

    eData[15] = temp[1];
    eData[15] += temp4_5_6;
    eData[15] += temp8_9;
    eData[15] += temp12_14;
    eData[15] += temp16_17_18_19;
    eData[15] += temp20_21_22;
    eData[15] += temp24_26;
    eData[15] += temp28_31;
    eData[15] += temp32_34_35;

    eData[16] += temp0_1_3;
    eData[16] += temp[6];
    eData[16] += temp[10];
    eData[16] += temp12_13_14;
    eData[16] += temp18_19;
    eData[16] += temp21_22_23;
    eData[16] += temp24_25_27;
    eData[16] += temp32_33_34;

    eData[17] = temp0_1_2;
    eData[17] += temp5_7;
    eData[17] += temp[11];
    eData[17] += temp12_13_15;
    eData[17] += temp[18];
    eData[17] += temp[22];
    eData[17] += temp24_25_26;
    eData[17] += temp28_29_30;
    eData[17] += temp32_33;

    eData[18] = temp1_2;
    eData[18] += temp4_7;
    eData[18] += temp8_9_11;
    eData[18] += temp12_15;
    eData[18] += temp17_19;
    eData[18] += temp20_21_22_23;
    eData[18] += temp24_25_27;
    eData[18] += temp29_31;
    eData[18] += temp34_35;

    eData[19] += temp0_1_3;
    eData[19] += temp4_6;
    eData[19] += temp[9];
    eData[19] += temp13_15;
    eData[19] += temp16_17;
    eData[19] += temp21_22;
    eData[19] += temp24_25_26_27;
    eData[19] += temp28_30;
    eData[19] += temp[35];

    eData[20] = temp0_3;
    eData[20] += temp4_5;
    eData[20] += temp8_10;
    eData[20] += temp14_15;
    eData[20] += temp16_18;
    eData[20] += temp[21];
    eData[20] += temp25_27;
    eData[20] += temp28_29_31;
    eData[20] += temp32_33_35;

    eData[21] = temp1_2;
    eData[21] += temp4_5_7;
    eData[21] += temp10_11;
    eData[21] += temp12_14_15;
    eData[21] += temp[18];
    eData[21] += temp20_22_23;
    eData[21] += temp24_25_26_27;
    eData[21] += temp28_30;
    eData[21] += temp32_34;

    eData[22] += temp2_3;
    eData[22] += temp4_6_7;
    eData[22] += temp[9];
    eData[22] += temp[12];
    eData[22] += temp16_18_19;
    eData[22] += temp[20];
    eData[22] += temp24_25_27;
    eData[22] += temp28_29_30_31;
    eData[22] += temp[33];

    eData[23] = temp0_2_3;
    eData[23] += temp6_7;
    eData[23] += temp8_11;
    eData[23] += temp[13];
    eData[23] += temp17_18_19;
    eData[23] += temp[21];
    eData[23] += temp[24];
    eData[23] += temp28_30_31;
    eData[23] += temp32_34_35;

    eData[24] = temp[1];
    eData[24] += temp4_5_7;
    eData[24] += temp8_10;
    eData[24] += temp13_14_15;
    eData[24] += temp17_18;
    eData[24] += temp21_23;
    eData[24] += temp25_26_27;
    eData[24] += temp28_29_30_31;
    eData[24] += temp33_35;

    eData[25] += temp[0];
    eData[25] += temp5_6_7;
    eData[25] += temp9_10;
    eData[25] += temp12_15;
    eData[25] += temp[19];
    eData[25] += temp21_22_23;
    eData[25] += temp[27];
    eData[25] += temp28_30_31;
    eData[25] += temp32_33_34;

    eData[26] = temp1_2_3;
    eData[26] += temp5_6;
    eData[26] += temp9_10_11;
    eData[26] += temp[14];
    eData[26] += temp[16];
    eData[26] += temp20_21_22;
    eData[26] += temp24_27;
    eData[26] += temp[31];
    eData[26] += temp33_34_35;

    eData[27] = temp0_2;
    eData[27] += temp4_7;
    eData[27] += temp8_10_11;
    eData[27] += temp[13];
    eData[27] += temp16_17_18;
    eData[27] += temp20_21;
    eData[27] += temp24_26;
    eData[27] += temp28_29_30_31;
    eData[27] += temp32_33_34;

    eData[28] += temp0_1_3;
    eData[28] += temp8_9_10;
    eData[28] += temp12_13_15;
    eData[28] += temp[18];
    eData[28] += temp[22];
    eData[28] += temp24_25_26;
    eData[28] += temp30_31;
    eData[28] += temp33_34_35;

    eData[29] = temp0_1_2;
    eData[29] += temp4_5_6;
    eData[29] += temp8_9;
    eData[29] += temp12_13_14;
    eData[29] += temp17_19;
    eData[29] += temp[23];
    eData[29] += temp24_25_27;
    eData[29] += temp[30];
    eData[29] += temp[34];

    eData[30] = temp0_1_3;
    eData[30] += temp5_7;
    eData[30] += temp10_11;
    eData[30] += temp13_14;
    eData[30] += temp16_19;
    eData[30] += temp20_21_23;
    eData[30] += temp24_27;
    eData[30] += temp29_31;
    eData[30] += temp32_33_34_35;

    eData[31] += temp0_1_2_3;
    eData[31] += temp4_6;
    eData[31] += temp[11];
    eData[31] += temp12_13_15;
    eData[31] += temp16_18;
    eData[31] += temp[21];
    eData[31] += temp25_27;
    eData[31] += temp28_29;
    eData[31] += temp33_34;

    eData[32] = temp1_3;
    eData[32] += temp4_5_7;
    eData[32] += temp8_9_11;
    eData[32] += temp12_15;
    eData[32] += temp16_17;
    eData[32] += temp20_22;
    eData[32] += temp26_27;
    eData[32] += temp28_30;
    eData[32] += temp[33];

    eData[33] = temp0_1_2_3;
    eData[33] += temp4_6;
    eData[33] += temp8_10;
    eData[33] += temp13_14;
    eData[33] += temp16_17_19;
    eData[33] += temp22_23;
    eData[33] += temp24_26_27;
    eData[33] += temp[30];
    eData[33] += temp32_34_35;

    eData[34] += temp0_1_3;
    eData[34] += temp4_5_6_7;
    eData[34] += temp[9];
    eData[34] += temp14_15;
    eData[34] += temp16_18_19;
    eData[34] += temp[21];
    eData[34] += temp[24];
    eData[34] += temp28_30_31;
    eData[34] += temp[32];

    eData[35] = temp[0];
    eData[35] += temp4_6_7;
    eData[35] += temp8_10_11;
    eData[35] += temp12_14_15;
    eData[35] += temp18_19;
    eData[35] += temp20_23;
    eData[35] += temp[25];
    eData[35] += temp29_30_31;
    eData[35] += temp[33];

    if (0)
    {
        temp = eData;
        // 0,1,2,3
        temp0_1 = temp[0];
        temp0_1 += temp[1];
        temp0_2 = temp[0];
        temp0_2 += temp[2];
        temp0_3 = temp[0];
        temp0_3 += temp[3];
        temp1_2 = temp[1];
        temp1_2 += temp[2];
        temp1_3 = temp[1];
        temp1_3 += temp[3];
        temp2_3 = temp[2];
        temp2_3 += temp[3];
        temp0_1_2 = temp0_1;
        temp0_1_2 += temp[2];
        temp0_1_3 = temp0_1;
        temp0_1_3 += temp[3];
        temp0_2_3 = temp0_2;
        temp0_2_3 += temp[3];
        temp1_2_3 = temp1_2;
        temp1_2_3 += temp[3];
        temp0_1_2_3 = temp0_1_2;
        temp0_1_2_3 += temp[3];
        // 4,5,6,7
        temp4_5 = temp[4];
        temp4_5 += temp[5];
        temp4_6 = temp[4];
        temp4_6 += temp[6];
        temp4_7 = temp[4];
        temp4_7 += temp[7];
        temp5_6 = temp[5];
        temp5_6 += temp[6];
        temp5_7 = temp[5];
        temp5_7 += temp[7];
        temp6_7 = temp[6];
        temp6_7 += temp[7];
        temp4_5_6 = temp4_5;
        temp4_5_6 += temp[6];
        temp4_5_7 = temp4_5;
        temp4_5_7 += temp[7];
        temp4_6_7 = temp4_6;
        temp4_6_7 += temp[7];
        temp5_6_7 = temp5_6;
        temp5_6_7 += temp[7];
        temp4_5_6_7 = temp4_5_6;
        temp4_5_6_7 += temp[7];
        // 8,9,10,11
        temp8_9 = temp[8];
        temp8_9 += temp[9];
        temp8_10 = temp[8];
        temp8_10 += temp[10];
        temp8_11 = temp[8];
        temp8_11 += temp[11];
        temp9_10 = temp[9];
        temp9_10 += temp[10];
        temp9_11 = temp[9];
        temp9_11 += temp[11];
        temp10_11 = temp[10];
        temp10_11 += temp[11];
        temp8_9_10 = temp8_9;
        temp8_9_10 += temp[10];
        temp8_9_11 = temp8_9;
        temp8_9_11 += temp[11];
        temp8_10_11 = temp8_10;
        temp8_10_11 += temp[11];
        temp9_10_11 = temp9_10;
        temp9_10_11 += temp[11];
        temp8_9_10_11 = temp8_9_10;
        temp8_9_10_11 += temp[11];
        // 12,13,14,15
        temp12_13 = temp[12];
        temp12_13 += temp[13];
        temp12_14 = temp[12];
        temp12_14 += temp[14];
        temp12_15 = temp[12];
        temp12_15 += temp[15];
        temp13_14 = temp[13];
        temp13_14 += temp[14];
        temp13_15 = temp[13];
        temp13_15 += temp[15];
        temp14_15 = temp[14];
        temp14_15 += temp[15];
        temp12_13_14 = temp12_13;
        temp12_13_14 += temp[14];
        temp12_13_15 = temp12_13;
        temp12_13_15 += temp[15];
        temp12_14_15 = temp12_14;
        temp12_14_15 += temp[15];
        temp13_14_15 = temp13_14;
        temp13_14_15 += temp[15];
        temp12_13_14_15 = temp12_13_14;
        temp12_13_14_15 += temp[15];
        // 16,17,18,19
        temp16_17 = temp[16];
        temp16_17 += temp[17];
        temp16_18 = temp[16];
        temp16_18 += temp[18];
        temp16_19 = temp[16];
        temp16_19 += temp[19];
        temp17_18 = temp[17];
        temp17_18 += temp[18];
        temp17_19 = temp[17];
        temp17_19 += temp[19];
        temp18_19 = temp[18];
        temp18_19 += temp[19];
        temp16_17_18 = temp16_17;
        temp16_17_18 += temp[18];
        temp16_17_19 = temp16_17;
        temp16_17_19 += temp[19];
        temp16_18_19 = temp16_18;
        temp16_18_19 += temp[19];
        temp17_18_19 = temp17_18;
        temp17_18_19 += temp[19];
        temp16_17_18_19 = temp16_17_18;
        temp16_17_18_19 += temp[19];
        // 20,21,22,23
        temp20_21 = temp[20];
        temp20_21 += temp[21];
        temp20_22 = temp[20];
        temp20_22 += temp[22];
        temp20_23 = temp[20];
        temp20_23 += temp[23];
        temp21_22 = temp[21];
        temp21_22 += temp[22];
        temp21_23 = temp[21];
        temp21_23 += temp[23];
        temp22_23 = temp[22];
        temp22_23 += temp[23];
        temp20_21_22 = temp20_21;
        temp20_21_22 += temp[22];
        temp20_21_23 = temp20_21;
        temp20_21_23 += temp[23];
        temp20_22_23 = temp20_22;
        temp20_22_23 += temp[23];
        temp21_22_23 = temp21_22;
        temp21_22_23 += temp[23];
        temp20_21_22_23 = temp20_21_22;
        temp20_21_22_23 += temp[23];
        // 24,25,26,27
        temp24_25 = temp[24];
        temp24_25 += temp[25];
        temp24_26 = temp[24];
        temp24_26 += temp[26];
        temp24_27 = temp[24];
        temp24_27 += temp[27];
        temp25_26 = temp[25];
        temp25_26 += temp[26];
        temp25_27 = temp[25];
        temp25_27 += temp[27];
        temp26_27 = temp[26];
        temp26_27 += temp[27];
        temp24_25_26 = temp24_25;
        temp24_25_26 += temp[26];
        temp24_25_27 = temp24_25;
        temp24_25_27 += temp[27];
        temp24_26_27 = temp24_26;
        temp24_26_27 += temp[27];
        temp25_26_27 = temp25_26;
        temp25_26_27 += temp[27];
        temp24_25_26_27 = temp24_25_26;
        temp24_25_26_27 += temp[27];
        // 28,29,30,31
        temp28_29 = temp[28];
        temp28_29 += temp[29];
        temp28_30 = temp[28];
        temp28_30 += temp[30];
        temp28_31 = temp[28];
        temp28_31 += temp[31];
        temp29_30 = temp[29];
        temp29_30 += temp[30];
        temp29_31 = temp[29];
        temp29_31 += temp[31];
        temp30_31 = temp[30];
        temp30_31 += temp[31];
        temp28_29_30 = temp28_29;
        temp28_29_30 += temp[30];
        temp28_29_31 = temp28_29;
        temp28_29_31 += temp[31];
        temp28_30_31 = temp28_30;
        temp28_30_31 += temp[31];
        temp29_30_31 = temp29_30;
        temp29_30_31 += temp[31];
        temp28_29_30_31 = temp28_29_30;
        temp28_29_30_31 += temp[31];
        // 32,33,34,35
        temp32_33 = temp[32];
        temp32_33 += temp[33];
        temp32_34 = temp[32];
        temp32_34 += temp[34];
        temp32_35 = temp[32];
        temp32_35 += temp[35];
        temp33_34 = temp[33];
        temp33_34 += temp[34];
        temp33_35 = temp[33];
        temp33_35 += temp[35];
        temp34_35 = temp[34];
        temp34_35 += temp[35];
        temp32_33_34 = temp32_33;
        temp32_33_34 += temp[34];
        temp32_33_35 = temp32_33;
        temp32_33_35 += temp[35];
        temp32_34_35 = temp32_34;
        temp32_34_35 += temp[35];
        temp33_34_35 = temp33_34;
        temp33_34_35 += temp[35];
        temp32_33_34_35 = temp32_33_34;
        temp32_33_34_35 += temp[35];
        eData[0] = temp1_2_3;
        eData[0] += temp4_5_6_7;
        eData[0] += temp9_11;
        eData[0] += temp[13];
        eData[0] += temp16_17_19;
        eData[0] += temp20_22;
        eData[0] += temp25_26_27;
        eData[0] += temp29_30;
        eData[0] += temp33_35;

        eData[1] += temp[3];
        eData[1] += temp4_6_7;
        eData[1] += temp8_9_10;
        eData[1] += temp[12];
        eData[1] += temp17_18_19;
        eData[1] += temp21_22;
        eData[1] += temp24_27;
        eData[1] += temp[31];
        eData[1] += temp33_34_35;

        eData[2] = temp0_3;
        eData[2] += temp[7];
        eData[2] += temp9_10_11;
        eData[2] += temp13_14_15;
        eData[2] += temp17_18;
        eData[2] += temp21_22_23;
        eData[2] += temp[26];
        eData[2] += temp[28];
        eData[2] += temp32_33_34;

        eData[3] = temp0_2;
        eData[3] += temp4_5_6_7;
        eData[3] += temp8_9_10;
        eData[3] += temp12_14;
        eData[3] += temp16_19;
        eData[3] += temp20_22_23;
        eData[3] += temp[25];
        eData[3] += temp28_29_30;
        eData[3] += temp32_33;

        eData[4] += temp0_1_2;
        eData[4] += temp6_7;
        eData[4] += temp9_10_11;
        eData[4] += temp12_13_15;
        eData[4] += temp20_21_22;
        eData[4] += temp24_25_27;
        eData[4] += temp[30];
        eData[4] += temp[34];

        eData[5] = temp0_1_3;
        eData[5] += temp[6];
        eData[5] += temp[10];
        eData[5] += temp12_13_14;
        eData[5] += temp16_17_18;
        eData[5] += temp20_21;
        eData[5] += temp24_25_26;
        eData[5] += temp29_31;
        eData[5] += temp[35];

        eData[6] = temp0_3;
        eData[6] += temp5_7;
        eData[6] += temp8_9_10_11;
        eData[6] += temp12_13_15;
        eData[6] += temp17_19;
        eData[6] += temp22_23;
        eData[6] += temp25_26;
        eData[6] += temp28_31;
        eData[6] += temp32_33_35;

        eData[7] += temp1_3;
        eData[7] += temp4_5;
        eData[7] += temp9_10;
        eData[7] += temp12_13_14_15;
        eData[7] += temp16_18;
        eData[7] += temp[23];
        eData[7] += temp24_25_27;
        eData[7] += temp28_30;
        eData[7] += temp[33];

        eData[8] = temp2_3;
        eData[8] += temp4_6;
        eData[8] += temp[9];
        eData[8] += temp13_15;
        eData[8] += temp16_17_19;
        eData[8] += temp20_21_23;
        eData[8] += temp24_27;
        eData[8] += temp28_29;
        eData[8] += temp32_34;

        eData[9] = temp0_2_3;
        eData[9] += temp[6];
        eData[9] += temp8_10_11;
        eData[9] += temp12_13_14_15;
        eData[9] += temp16_18;
        eData[9] += temp20_22;
        eData[9] += temp25_26;
        eData[9] += temp28_29_31;
        eData[9] += temp34_35;

        eData[10] += temp[0];
        eData[10] += temp4_6_7;
        eData[10] += temp[8];
        eData[10] += temp12_13_15;
        eData[10] += temp16_17_18_19;
        eData[10] += temp[21];
        eData[10] += temp26_27;
        eData[10] += temp28_30_31;
        eData[10] += temp[33];

        eData[11] = temp[1];
        eData[11] += temp5_6_7;
        eData[11] += temp[9];
        eData[11] += temp[12];
        eData[11] += temp16_18_19;
        eData[11] += temp20_22_23;
        eData[11] += temp24_26_27;
        eData[11] += temp30_31;
        eData[11] += temp32_35;

        eData[12] = temp1_2_3;
        eData[12] += temp5_6;
        eData[12] += temp9_11;
        eData[12] += temp13_14_15;
        eData[12] += temp16_17_18_19;
        eData[12] += temp21_23;
        eData[12] += temp[25];
        eData[12] += temp28_29_31;
        eData[12] += temp32_34;

        eData[13] += temp0_3;
        eData[13] += temp[7];
        eData[13] += temp9_10_11;
        eData[13] += temp[15];
        eData[13] += temp16_18_19;
        eData[13] += temp20_21_22;
        eData[13] += temp[24];
        eData[13] += temp29_30_31;
        eData[13] += temp33_34;

        eData[14] = temp[2];
        eData[14] += temp[4];
        eData[14] += temp8_9_10;
        eData[14] += temp12_15;
        eData[14] += temp[19];
        eData[14] += temp21_22_23;
        eData[14] += temp25_26_27;
        eData[14] += temp29_30;
        eData[14] += temp33_34_35;

        eData[15] = temp[1];
        eData[15] += temp4_5_6;
        eData[15] += temp8_9;
        eData[15] += temp12_14;
        eData[15] += temp16_17_18_19;
        eData[15] += temp20_21_22;
        eData[15] += temp24_26;
        eData[15] += temp28_31;
        eData[15] += temp32_34_35;

        eData[16] += temp0_1_3;
        eData[16] += temp[6];
        eData[16] += temp[10];
        eData[16] += temp12_13_14;
        eData[16] += temp18_19;
        eData[16] += temp21_22_23;
        eData[16] += temp24_25_27;
        eData[16] += temp32_33_34;

        eData[17] = temp0_1_2;
        eData[17] += temp5_7;
        eData[17] += temp[11];
        eData[17] += temp12_13_15;
        eData[17] += temp[18];
        eData[17] += temp[22];
        eData[17] += temp24_25_26;
        eData[17] += temp28_29_30;
        eData[17] += temp32_33;

        eData[18] = temp1_2;
        eData[18] += temp4_7;
        eData[18] += temp8_9_11;
        eData[18] += temp12_15;
        eData[18] += temp17_19;
        eData[18] += temp20_21_22_23;
        eData[18] += temp24_25_27;
        eData[18] += temp29_31;
        eData[18] += temp34_35;

        eData[19] += temp0_1_3;
        eData[19] += temp4_6;
        eData[19] += temp[9];
        eData[19] += temp13_15;
        eData[19] += temp16_17;
        eData[19] += temp21_22;
        eData[19] += temp24_25_26_27;
        eData[19] += temp28_30;
        eData[19] += temp[35];

        eData[20] = temp0_3;
        eData[20] += temp4_5;
        eData[20] += temp8_10;
        eData[20] += temp14_15;
        eData[20] += temp16_18;
        eData[20] += temp[21];
        eData[20] += temp25_27;
        eData[20] += temp28_29_31;
        eData[20] += temp32_33_35;

        eData[21] = temp1_2;
        eData[21] += temp4_5_7;
        eData[21] += temp10_11;
        eData[21] += temp12_14_15;
        eData[21] += temp[18];
        eData[21] += temp20_22_23;
        eData[21] += temp24_25_26_27;
        eData[21] += temp28_30;
        eData[21] += temp32_34;

        eData[22] += temp2_3;
        eData[22] += temp4_6_7;
        eData[22] += temp[9];
        eData[22] += temp[12];
        eData[22] += temp16_18_19;
        eData[22] += temp[20];
        eData[22] += temp24_25_27;
        eData[22] += temp28_29_30_31;
        eData[22] += temp[33];

        eData[23] = temp0_2_3;
        eData[23] += temp6_7;
        eData[23] += temp8_11;
        eData[23] += temp[13];
        eData[23] += temp17_18_19;
        eData[23] += temp[21];
        eData[23] += temp[24];
        eData[23] += temp28_30_31;
        eData[23] += temp32_34_35;

        eData[24] = temp[1];
        eData[24] += temp4_5_7;
        eData[24] += temp8_10;
        eData[24] += temp13_14_15;
        eData[24] += temp17_18;
        eData[24] += temp21_23;
        eData[24] += temp25_26_27;
        eData[24] += temp28_29_30_31;
        eData[24] += temp33_35;

        eData[25] += temp[0];
        eData[25] += temp5_6_7;
        eData[25] += temp9_10;
        eData[25] += temp12_15;
        eData[25] += temp[19];
        eData[25] += temp21_22_23;
        eData[25] += temp[27];
        eData[25] += temp28_30_31;
        eData[25] += temp32_33_34;

        eData[26] = temp1_2_3;
        eData[26] += temp5_6;
        eData[26] += temp9_10_11;
        eData[26] += temp[14];
        eData[26] += temp[16];
        eData[26] += temp20_21_22;
        eData[26] += temp24_27;
        eData[26] += temp[31];
        eData[26] += temp33_34_35;

        eData[27] = temp0_2;
        eData[27] += temp4_7;
        eData[27] += temp8_10_11;
        eData[27] += temp[13];
        eData[27] += temp16_17_18;
        eData[27] += temp20_21;
        eData[27] += temp24_26;
        eData[27] += temp28_29_30_31;
        eData[27] += temp32_33_34;

        eData[28] += temp0_1_3;
        eData[28] += temp8_9_10;
        eData[28] += temp12_13_15;
        eData[28] += temp[18];
        eData[28] += temp[22];
        eData[28] += temp24_25_26;
        eData[28] += temp30_31;
        eData[28] += temp33_34_35;

        eData[29] = temp0_1_2;
        eData[29] += temp4_5_6;
        eData[29] += temp8_9;
        eData[29] += temp12_13_14;
        eData[29] += temp17_19;
        eData[29] += temp[23];
        eData[29] += temp24_25_27;
        eData[29] += temp[30];
        eData[29] += temp[34];

        eData[30] = temp0_1_3;
        eData[30] += temp5_7;
        eData[30] += temp10_11;
        eData[30] += temp13_14;
        eData[30] += temp16_19;
        eData[30] += temp20_21_23;
        eData[30] += temp24_27;
        eData[30] += temp29_31;
        eData[30] += temp32_33_34_35;

        eData[31] += temp0_1_2_3;
        eData[31] += temp4_6;
        eData[31] += temp[11];
        eData[31] += temp12_13_15;
        eData[31] += temp16_18;
        eData[31] += temp[21];
        eData[31] += temp25_27;
        eData[31] += temp28_29;
        eData[31] += temp33_34;

        eData[32] = temp1_3;
        eData[32] += temp4_5_7;
        eData[32] += temp8_9_11;
        eData[32] += temp12_15;
        eData[32] += temp16_17;
        eData[32] += temp20_22;
        eData[32] += temp26_27;
        eData[32] += temp28_30;
        eData[32] += temp[33];

        eData[33] = temp0_1_2_3;
        eData[33] += temp4_6;
        eData[33] += temp8_10;
        eData[33] += temp13_14;
        eData[33] += temp16_17_19;
        eData[33] += temp22_23;
        eData[33] += temp24_26_27;
        eData[33] += temp[30];
        eData[33] += temp32_34_35;

        eData[34] += temp0_1_3;
        eData[34] += temp4_5_6_7;
        eData[34] += temp[9];
        eData[34] += temp14_15;
        eData[34] += temp16_18_19;
        eData[34] += temp[21];
        eData[34] += temp[24];
        eData[34] += temp28_30_31;
        eData[34] += temp[32];

        eData[35] = temp[0];
        eData[35] += temp4_6_7;
        eData[35] += temp8_10_11;
        eData[35] += temp12_14_15;
        eData[35] += temp18_19;
        eData[35] += temp20_23;
        eData[35] += temp[25];
        eData[35] += temp29_30_31;
        eData[35] += temp[33];
    }
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
    for (long j = 3; j < BlockByte; j += 3)
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
    for (long j = 3; j < BlockByte; j += 3)
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
    std::cout << "BlockWord: " << BlockByte << std::endl;
    std::cout << "PlainMod: " << PlainMod << std::endl;
    std::cout << "PlainBlock: " << PlainBlock << std::endl;
    std::cout << "nslots: " << nslots << std::endl;
    std::cout << "Para_bits: " << Para_bits << std::endl;
    std::cout << "Para_c: " << Para_c << std::endl;
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

    // if (nslots > PlainBlock)
    // {
    //     std::cerr << "nslots > PlainBlock" << std::endl;
    //     return false;
    // }

    // 创建PAlgebra对象
    const helib::PAlgebra &zMStar = context->getZMStar();
    // long minv = -4;
    long root; // m-th root of unity modulo p
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
    double SecurityLevel = context->securityLevel();

    auto start_keyEncryption = std::chrono::high_resolution_clock::now();
    vector<Ctxt> encryptedSymKey;
    vector<Ctxt> encryptedSymKey01;
    vector<Ctxt> encryptedSymKey02;
    encryptSymKey(encryptedSymKey, encryptedSymKey01, encryptedSymKey02, SymKey, publicKey, cmodulus);
    auto end_keyEncryption = std::chrono::high_resolution_clock::now();
    double keyEncryption = std::chrono::duration<double>(end_keyEncryption - start_keyEncryption).count();
    std::cout << "SymKey FHE time: " << keyEncryption << "s\n";
    //  解密验证
    if (symkeyflag)
    {
        if (!verify_encryptSymKey(encryptedSymKey, SymKey, secretKey, cmodulus))
        {
            return 0;
        }
        std::cout << "Symmetric key encryption succeeded!" << std::endl;
    }
    // Generating key stream
    auto start_keyStream = std::chrono::high_resolution_clock::now();

    std::vector<long> NonceSet(PlainBlock);
    std::vector<long> RoundKeySet(PlainByte * (Nr + 1));
    std::vector<long> KeyStream(PlainByte);
    RandomBit<BlockByte * randbits> randomBit(Nr);
    long X;
    std::vector<long> RoundKey(BlockByte);
    bool bit_array[randbits];
    long nonce;
    long block_num;
    long ir;
    std::vector<long> state(BlockByte); // 初始化 state
    std::cout << "Generating KeyStream..." << std::endl;
    uint64_t start_cycle = rdtsc();
    for (long counter = counter_begin; counter <= counter_end; counter++)
    {
        nonce = generate_secure_random_int(NonceSize);
        randomBit.generate_Instance_all_new(nonce, counter);
        auto &RanVecs = randomBit.roundconstants;
        block_num = counter - counter_begin;
        NonceSet[block_num] = nonce;
        // 逐轮进行加密
        for (unsigned r = 0; r <= Nr; r++)
        {
            // 计算 Xc 和 RoundKey
            for (unsigned i = 0; i < BlockByte; ++i)
            {
                ir = i * randbits;
                for (unsigned j = 0; j < randbits; ++j)
                {
                    bit_array[j] = RanVecs[r][ir + j];
                }
                BinStrToHex(bit_array, X, randbits);
                // 强制转换为 long 类型
                // Xc[i] = static_cast<long>(temp);
                // if (Rkflag)
                // {
                //     RoundKey[i] = (SymKey[i] * Xc[i]) % PlainMod;
                // }
                // else
                // {
                //     RoundKey[i] = (SymKey[i] + Xc[i]) % PlainMod;
                // }
                RoundKey[i] = (SymKey[i] * X) % PlainMod;
            }
            // 将RoundKey 复制到RoundKeySet
            // 测试使用
            memcpy(&RoundKeySet[PlainByte * r + BlockByte * (block_num)], RoundKey.data(), BlockByte * sizeof(long));
            if (r == 0)
            { // 初始轮
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (RoundKey[i] + IV[i]) % PlainMod;
                }
            }
            else if (r < Nr)
            {                       // 常规轮
                yusP.Sbox_5(state); // S盒
                yusP.M36_5(state);  // 线性变换
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
            }
            else
            {                      // 最后一轮
                yusP.M36_5(state); // 线性变换
                // yusP.Sbox_5_last(state); // S盒
                yusP.Sbox_5(state); // S盒
                yusP.M36_5(state);  // 线性变换
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
                memcpy(&KeyStream[(block_num)*BlockByte], state.data(), BlockByte * sizeof(long));
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
    std::cout << "Generating XOF stream..." << std::endl;
    std::vector<vec_long> Xset01(len3);
    for (int i = 0; i < len3; i++)
    {
        Xset01[i].SetLength(PlainBlock);
    }
    std::vector<vec_long> Xset02 = Xset01;
    std::vector<vec_long> Xset(NrBlockByte);
    for (int i = 0; i < NrBlockByte; i++)
    {
        Xset[i].SetLength(PlainBlock);
    }
    long rB;
    vec_long Xc;
    Xc.SetLength(BlockByte);
    auto start_XOF = std::chrono::high_resolution_clock::now();
    for (long counter = counter_begin; counter <= counter_end; counter++)
    {
        block_num = counter - counter_begin;
        nonce = NonceSet[block_num];
        randomBit.generate_Instance_all_new(nonce, counter);
        auto &RanVecs = randomBit.roundconstants;
        // 逐轮进行加密
        for (unsigned r = 0; r <= Nr; r++)
        {
            // 计算 Xc
            for (unsigned i = 0; i < BlockByte; ++i)
            {
                ir = i * randbits;
                for (unsigned j = 0; j < randbits; ++j)
                {
                    bit_array[j] = RanVecs[r][ir + j];
                }
                BinStrToHex(bit_array, Xc[i], randbits);
            }
            // 将 Xc 复制到 Xset
            rB = r * BlockByte;
            for (int i = 0; i < BlockByte; i++)
            {
                Xset[rB + i][block_num] = Xc[i];
            }
            // if (r == 0)
            // { // 初始轮
            //     for (int i = 0; i < BlockByte; i += 3)
            //     {
            //         Xset01[i/3][block_num] = (Xc[i] * Xc[i + 1]) % PlainMod;
            //         Xset02[i/3][block_num] = (Xc[i] * Xc[i + 2]) % PlainMod;
            //     }
            // }

            if (r == 0)
            { // 初始轮
                // for (int i = 0; i < BlockByte; i += 3)
                // {
                //     Xset01[i/3][block_num] = (Xc[i] * Xc[i + 1]) % PlainMod;
                //     Xset02[i/3][block_num] = (Xc[i] * Xc[i + 2]) % PlainMod;
                // }
                Xset01[0][block_num] = (Xc[0] * Xc[1]) % PlainMod;
                Xset02[0][block_num] = (Xc[0] * Xc[2]) % PlainMod;
                Xset01[1][block_num] = (Xc[3] * Xc[4]) % PlainMod;
                Xset02[1][block_num] = (Xc[3] * Xc[5]) % PlainMod;
                Xset01[2][block_num] = (Xc[6] * Xc[7]) % PlainMod;
                Xset02[2][block_num] = (Xc[6] * Xc[8]) % PlainMod;
                Xset01[3][block_num] = (Xc[9] * Xc[10]) % PlainMod;
                Xset02[3][block_num] = (Xc[9] * Xc[11]) % PlainMod;
                Xset01[4][block_num] = (Xc[12] * Xc[13]) % PlainMod;
                Xset02[4][block_num] = (Xc[12] * Xc[14]) % PlainMod;
                Xset01[5][block_num] = (Xc[15] * Xc[16]) % PlainMod;
                Xset02[5][block_num] = (Xc[15] * Xc[17]) % PlainMod;
                Xset01[6][block_num] = (Xc[18] * Xc[19]) % PlainMod;
                Xset02[6][block_num] = (Xc[18] * Xc[20]) % PlainMod;
                Xset01[7][block_num] = (Xc[21] * Xc[22]) % PlainMod;
                Xset02[7][block_num] = (Xc[21] * Xc[23]) % PlainMod;
                Xset01[8][block_num] = (Xc[24] * Xc[25]) % PlainMod;
                Xset02[8][block_num] = (Xc[24] * Xc[26]) % PlainMod;
                Xset01[9][block_num] = (Xc[27] * Xc[28]) % PlainMod;
                Xset02[9][block_num] = (Xc[27] * Xc[29]) % PlainMod;
                Xset01[10][block_num] = (Xc[30] * Xc[31]) % PlainMod;
                Xset02[10][block_num] = (Xc[30] * Xc[32]) % PlainMod;
                Xset01[11][block_num] = (Xc[33] * Xc[34]) % PlainMod;
                Xset02[11][block_num] = (Xc[33] * Xc[35]) % PlainMod;
            }
        }
    }
    //         for (int j = 0; j < PlainBlock; j++)
    //     {
    //         Xset01[0][j] = (Xset[0][j] * Xset[1][j]) % PlainMod;
    //         Xset02[0][j] = (Xset[0][j] * Xset[2][j]) % PlainMod;
    //         Xset01[1][j] = (Xset[3][j] * Xset[4][j]) % PlainMod;
    //         Xset02[1][j] = (Xset[3][j] * Xset[5][j]) % PlainMod;
    //         Xset01[2][j] = (Xset[6][j] * Xset[7][j]) % PlainMod;
    //         Xset02[2][j] = (Xset[6][j] * Xset[8][j]) % PlainMod;
    //         Xset01[3][j] = (Xset[9][j] * Xset[10][j]) % PlainMod;
    //         Xset02[3][j] = (Xset[9][j] * Xset[11][j]) % PlainMod;
    //         Xset01[4][j] = (Xset[12][j] * Xset[13][j]) % PlainMod;
    //         Xset02[4][j] = (Xset[12][j] * Xset[14][j]) % PlainMod;
    //         Xset01[5][j] = (Xset[15][j] * Xset[16][j]) % PlainMod;
    //         Xset02[5][j] = (Xset[15][j] * Xset[17][j]) % PlainMod;
    //         Xset01[6][j] = (Xset[18][j] * Xset[19][j]) % PlainMod;
    //         Xset02[6][j] = (Xset[18][j] * Xset[20][j]) % PlainMod;
    //         Xset01[7][j] = (Xset[21][j] * Xset[22][j]) % PlainMod;
    //         Xset02[7][j] = (Xset[21][j] * Xset[23][j]) % PlainMod;
    //         Xset01[8][j] = (Xset[24][j] * Xset[25][j]) % PlainMod;
    //         Xset02[8][j] = (Xset[24][j] * Xset[26][j]) % PlainMod;
    //         Xset01[9][j] = (Xset[27][j] * Xset[28][j]) % PlainMod;
    //         Xset02[9][j] = (Xset[27][j] * Xset[29][j]) % PlainMod;
    //         Xset01[10][j] = (Xset[30][j] * Xset[31][j]) % PlainMod;
    //         Xset02[10][j] = (Xset[30][j] * Xset[32][j]) % PlainMod;
    //         Xset01[11][j] = (Xset[33][j] * Xset[34][j]) % PlainMod;
    //         Xset02[11][j] = (Xset[33][j] * Xset[35][j]) % PlainMod;
    // }
    auto end_XOF = std::chrono::high_resolution_clock::now();
    double XOF_time = std::chrono::duration<double>(end_XOF - start_XOF).count();
    std::cout << "XOF stream Generation time: " << XOF_time << "s\n";

    vector<zzX> encodedXset(NrBlockByte);
    vector<zzX> encodedXset01(len3);
    vector<zzX> encodedXset02(len3);
    zz_pX encodedtemp;
    auto start_Xset = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NrBlockByte; i++)
    {
        cmodulus.iFFT(encodedtemp, Xset[i]);
        convert(encodedXset[i], encodedtemp);
    }
    for (int i = 0; i < len3; i++)
    {
        cmodulus.iFFT(encodedtemp, Xset01[i]);
        convert(encodedXset01[i], encodedtemp);
        cmodulus.iFFT(encodedtemp, Xset02[i]);
        convert(encodedXset02[i], encodedtemp);
    }
    auto end_Xset = std::chrono::high_resolution_clock::now();
    double Encode_time = std::chrono::duration<double>(end_Xset - start_Xset).count();
    std::cout << "encode time: " << Encode_time << "s\n";
    std::cout << "encode time/slots:" << Encode_time / (NrBlockByte + 2 * len3) * pow(10, 6) << "us\n";

    Ctxt tmpCtxt(*publicKey);
    int noise_budget = min_noise_budget(encryptedSymKey);
    std::cout << "noise budget initially: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    // 计算 encryptedRoundKeySet
    long eRk_len = BlockByte * (Nr + 1);
    vector<Ctxt> encryptedRoundKeySet(eRk_len, tmpCtxt);
    // for (int i = 0; i < eRk_len; i++)
    // {
    //     encryptedRoundKeySet[i] = encryptedSymKey[i % BlockByte];
    // }
    int index = 0;
    for (int i = 0; i <= Nr; ++i)
    {
        for (int j = 0; j < BlockByte; ++j)
        {
            encryptedRoundKeySet[index++] = encryptedSymKey[j];
        }
    }
    auto start_RoundKeySet_FHE = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < eRk_len; i++)
    {
        encryptedRoundKeySet[i].multByConstant(encodedXset[i]);
    }
    // if (Rkflag)
    // {
    //     for (int i = 0; i < eRk_len; i++)
    //     {
    //         encryptedRoundKeySet[i].multByConstant(encodedXset[i]);
    //     }
    // }
    // else
    // {
    //     ZZX temp;
    //     for (int i = 0; i < eRk_len; i++)
    //     {
    //             for (long j = 0; j < nslots; ++j) {
    //                 SetCoeff(temp, j, encodedXset[i][j]);
    //             }
    //         encryptedRoundKeySet[i].addConstant(temp);
    //     }
    // }
    auto end_RoundKeySet_FHE = std::chrono::high_resolution_clock::now();
    double RoundKey_time = std::chrono::duration<double>(end_RoundKeySet_FHE - start_RoundKeySet_FHE).count();
    std::cout << "RoundKeySet FHE succeeded! Time: " << RoundKey_time << "s\n";
    noise_budget = min_noise_budget(encryptedRoundKeySet);
    std::cout << "noise budget after RoundKeySet FHE: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    // 使用 verifyDecryption 函数解密并验证 RoundKeySet
    // if (deflag)
    // {
    //     if (!verifyDecryption18(encryptedRoundKeySet, RoundKeySet, secretKey, cmodulus,BlockByte))
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
    vector<Ctxt> encryptedKeyStream(encryptedRoundKeySet.begin(), encryptedRoundKeySet.begin() + BlockByte);
    std::cout << "whiteround + sbox start" << std::endl;
    Ctxt K0Iv1R0(encryptedSymKey[0]);
    Ctxt K1Iv0R1(encryptedSymKey[0]);
    Ctxt K0Iv2R0(encryptedSymKey[0]);
    Ctxt K2Iv0R2(encryptedSymKey[0]);
    vector<Ctxt> tmp01 = encryptedSymKey01;
    vector<Ctxt> tmp02 = encryptedSymKey02;
    start_sbox = std::chrono::high_resolution_clock::now();
    for (long i = 0; i < BlockByte; i += 3) // BlockByte
    {
        index = i / 3;
        K1Iv0R1 = encryptedSymKey[i + 1];
        K1Iv0R1.multByConstant(multiplyAndMod(encodedXset[i + 1], IV[i]));
        K0Iv1R0 = encryptedSymKey[i + 0];
        K0Iv1R0.multByConstant(multiplyAndMod(encodedXset[i], IV[i + 1]));
        K2Iv0R2 = encryptedSymKey[i + 2];
        K2Iv0R2.multByConstant(multiplyAndMod(encodedXset[i + 2], IV[i]));
        K0Iv2R0 = encryptedSymKey[i + 0];
        K0Iv2R0.multByConstant(multiplyAndMod(encodedXset[i], IV[i + 2]));

        // 计算Sbox 0
        encryptedKeyStream[i].addConstant(IV[i]);

        // 计算Sbox 1
        tmp02[index].multByConstant(encodedXset02[index]);
        tmp02[index] += K2Iv0R2;
        tmp02[index] += K0Iv2R0;
        tmp02[index].addConstant((IV[i] * IV[i + 2]) % PlainMod);
        encryptedKeyStream[i + 1].addConstant(IV[i + 1]);
        encryptedKeyStream[i + 1] += tmp02[index];

        // 计算Sbox 2
        tmp01[index].multByConstant(encodedXset01[index]);
        tmp01[index] += K1Iv0R1;
        tmp01[index] += K0Iv1R0;
        tmp01[index].addConstant((IV[i] * IV[i + 1] - IV[i + 2]) % PlainMod);
        encryptedKeyStream[i + 2] -= tmp01[index];
        encryptedKeyStream[i + 2] += tmp02[index];
    }
    end_sbox = std::chrono::high_resolution_clock::now();
    Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
    sbox_set[0] = Sbox_time;
    // 输出 Sbox_time
    std::cout << "SboxAndWhiteround time: " << Sbox_time << "s\n";
    noise_budget = min_noise_budget(encryptedKeyStream);
    std::cout << "noise budget after SboxAndWhiteround: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    // 明文密钥流
    vector<long> KeyStream2(PlainByte);
    if (deflag)
    {
        // 对IV和RoundKeySet进行异或
        for (long i = 0; i < PlainByte; i++)
        {
            KeyStream2[i] = (IV[i % BlockByte] + RoundKeySet[i]) % PlainMod;
        }
        // sbox
        yusP.Sbox_5(KeyStream2);
        // 使用 verifyDecryption 函数解密并验证 KeyStream
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockByte))
        {
            std::cerr << "Decryption verification failed for KeyStream." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for whiteround." << std::endl;
    }
    for (long r = 1; r < Nr; r++)
    {
        std::cout << "Round " << r << " start" << std::endl;
        if (r > 1)
        {
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
                return 0;
            }
            if (deflag)
            {
                yusP.Sbox_5(KeyStream2);
                if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockByte))
                {
                    std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
                    return 0;
                }
                std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
            }
        }
        start_linear = std::chrono::high_resolution_clock::now();
        // Linear Layer
        HE_M2(encryptedKeyStream);
        end_linear = std::chrono::high_resolution_clock::now();
        Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
        linear_set[r - 1] = std::chrono::duration<double>(end_linear - start_linear).count();
        noise_budget = min_noise_budget(encryptedKeyStream);
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
                vector<long> tmp(BlockByte);
                for (int j = 0; j < BlockByte; j++)
                {
                    tmp[j] = KeyStream2[i * BlockByte + j];
                }
                yusP.M36_5(tmp);
                for (int j = 0; j < BlockByte; j++)
                {
                    KeyStream2[i * BlockByte + j] = tmp[j];
                }
            }
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockByte))
            {
                std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
        }
        // Round Key Addition
        start_add = std::chrono::high_resolution_clock::now();
        // omp_set_num_threads(12); // 设置线程数为12
        // #pragma omp parallel for
        for (long j = 0; j < BlockByte; j++)
        {
            encryptedKeyStream[j] += encryptedRoundKeySet[r * BlockByte + j];
        }
        end_add = std::chrono::high_resolution_clock::now();
        Add_time += std::chrono::duration<double>(end_add - start_add).count();
        noise_budget = min_noise_budget(encryptedKeyStream);
        std::cout << "noise budget after Add: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            return 0;
        }
        if (deflag)
        {
            for (long i = 0; i < PlainByte; i++)
            {
                KeyStream2[i] = (KeyStream2[i] + RoundKeySet[r * PlainByte + i]) % PlainMod;
            }
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockByte))
            {
                std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;
        }
    }
// 最后一轮
#if (1)
    std::cout << "the last Round " << Nr << " start" << std::endl;
    start_linear = std::chrono::high_resolution_clock::now();
    // Linear Layer
    HE_M2(encryptedKeyStream);
    end_linear = std::chrono::high_resolution_clock::now();
    Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
    linear_set[Nr - 1] = std::chrono::duration<double>(end_linear - start_linear).count();
    noise_budget = min_noise_budget(encryptedKeyStream);
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
            vector<long> tmp(BlockByte);
            for (int j = 0; j < BlockByte; j++)
            {
                tmp[j] = KeyStream2[i * BlockByte + j];
            }
            yusP.M36_5(tmp);
            for (int j = 0; j < BlockByte; j++)
            {
                KeyStream2[i * BlockByte + j] = tmp[j];
            }
        }
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockByte))
        {
            std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
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
        return 0;
    }
    if (deflag)
    {
        // yusP.Sbox_5_last(KeyStream2);
        yusP.Sbox_5(KeyStream2);
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockByte))
        {
            std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
    }
    start_linear = std::chrono::high_resolution_clock::now();
    // Linear Layer
    HE_M2(encryptedKeyStream);
    end_linear = std::chrono::high_resolution_clock::now();
    Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
    linear_set[Nr] = std::chrono::duration<double>(end_linear - start_linear).count();
    noise_budget = min_noise_budget(encryptedKeyStream);
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
            vector<long> tmp(BlockByte);
            for (int j = 0; j < BlockByte; j++)
            {
                tmp[j] = KeyStream2[i * BlockByte + j];
            }
            yusP.M36_5(tmp);
            for (int j = 0; j < BlockByte; j++)
            {
                KeyStream2[i * BlockByte + j] = tmp[j];
            }
        }
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockByte))
        {
            std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
    }
    // add
    start_add = std::chrono::high_resolution_clock::now();
    for (long j = 0; j < BlockByte; j++)
    {
        encryptedKeyStream[j] += encryptedRoundKeySet[Nr * BlockByte + j];
    }
    end_add = std::chrono::high_resolution_clock::now();
    Add_time += std::chrono::duration<double>(end_add - start_add).count();
    noise_budget = min_noise_budget(encryptedKeyStream);
    std::cout << "noise budget after Add: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    if (deflag)
    {
        for (long i = 0; i < PlainByte; i++)
        {
            KeyStream2[i] = (KeyStream2[i] + RoundKeySet[Nr * PlainByte + i]) % PlainMod;
        }
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockByte))
        {
            std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;
    }
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
    std::cout << "sbox_timeset: " << sbox_set << endl;
    std::cout << "linear_timeset: " << linear_set << endl;
    if (plainflag)
    {
        for (int i = 0; i < encryptedKeyStream.size(); i++)
        {
            encryptedKeyStream[i].bringToSet(encryptedKeyStream[i].naturalPrimeSet());
        }
        if (!verifyDecryption(encryptedKeyStream, KeyStream, secretKey, cmodulus, BlockByte))
        {
            std::cerr << "Decryption verification failed for KeyStream." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream." << std::endl;
    }
    // 客户端在线
    // 生成随机对称明文流，只用于测试
    random_device rd;
    vector<long> PlainStream(TruncPlainByte);
    for (int i = 0; i < TruncPlainByte; i++)
    {
        PlainStream[i] = rd() % PlainMod;
    }
    // 加密
    vector<vec_long> CipherStream(TruncByte);
    for (int i = 0; i < TruncByte; i++)
    {
        CipherStream[i].SetLength(PlainBlock);
    }
    auto start_ClientOnline = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < PlainBlock; j++)
    {
        for (int i = 0; i < TruncByte; i++)
        {
            CipherStream[i][j] = ((PlainStream[j * TruncByte + i] - KeyStream[j * BlockByte + i]) % PlainMod + PlainMod) % PlainMod;
        }
    }
    auto end_ClientOnline = std::chrono::high_resolution_clock::now();
    double Client_ontime = std::chrono::duration<double>(end_ClientOnline - start_ClientOnline).count();
    std::cout << "Client onine total time:" << std::chrono::duration<double>(end_ClientOnline - start_ClientOnline).count() << "s\n";
    double Client_totaltime = Client_offtime + Client_ontime;
    std::cout << "Client total time: " << Client_totaltime << "s\n";
    // 服务端在线
    // 同态加密
    vector<Ctxt> TruncencryptedKeyStream(encryptedKeyStream.begin(), encryptedKeyStream.begin() + TruncByte);
    vector<Ctxt> encrypedPlainStream = TruncencryptedKeyStream;
    // 对CipherStream进行编码
    vector<ZZX> encodedCipherStream(TruncByte);
    auto start_ServerOnline = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < TruncByte; i++)
    {
        cmodulus.iFFT(encodedtemp, CipherStream[i]);
        convert(encodedCipherStream[i], encodedtemp);
    }
    for (int i = 0; i < TruncByte; i++)
    {
        encrypedPlainStream[i].addConstant(encodedCipherStream[i]);
    }
    auto end_ServerOnline = std::chrono::high_resolution_clock::now();
    double server_ontime = std::chrono::duration<double>(end_ServerOnline - start_ServerOnline).count();
    std::cout << "Server onine total time:" << server_ontime << "s\n";
    noise_budget = min_noise_budget(encryptedKeyStream);
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
    if (!verifyDecryption(encrypedPlainStream, PlainStream, secretKey, cmodulus, TruncByte))
    {
        std::cerr << "Decryption verification failed for encrypedPlainStream." << std::endl;
        return 0;
    }
    std::cout << "Decryption verification succeeded for encrypedPlainStream." << std::endl;
    // 计算吞吐量,KiB/min
    double Server_totaltime = Server_offtime + server_ontime;
    double throughput = (Plainbits * 60) / (pow(2, 13) * Server_totaltime);
    throughput = throughput * TruncRate;
    std::cout << "Server total time: " << Server_totaltime << "s\n";
    std::cout << "Throughput: " << throughput << "KiB/min\n";
    std::string dirPath = "../tests";
    std::string filePath;
    if (!fs::exists(dirPath))
    {
        filePath = "test_Yus_p_C36_ClientAndServer5_2.txt";
    }
    else
    {
        filePath = "../tests/test_Yus_p_C36_ClientAndServer5_2.txt";
    }
    std::ofstream outfile(filePath, std::ios::app);
    if (!outfile)
    {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return 0;
    }
    outfile << std::left << std::setw(3) << Nr
            << std::left << std::setw(10) << Para_p
            << std::left << std::setw(10) << nslots
            << std::left << std::setw(5) << Para_bits
            << std::left << std::setw(4) << Para_c
            << std::left << std::setw(6) << Qbits
            << std::fixed << std::setprecision(3)
            << std::left << std::setw(14) << SecurityLevel
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
    std::cout << "test_Yus_p_C36_ClientAndServer5_2.txt updated." << std::endl;
    return 0;
}