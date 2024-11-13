#include <iostream>
#include <cstring>
#include <stdint.h>
#include <chrono>
#include <vector>
#include <omp.h>
#include <array>
#include <cstdint>
#include <cstring>
#include <atomic>
#include <cmath>
#include <chrono>
#include <fstream>
#include <memory>
#include <random>
#include <climits>
#include <thread>

#include <helib/helib.h>
#include <helib/ArgMap.h>
#include <helib/DoubleCRT.h>

#include "random_bit.hpp"
#include "Yus_p.hpp"
#include "tool.hpp"

using namespace std;
using namespace helib;
using namespace NTL;
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

static const unsigned Nr = 4;         // 轮数
static const unsigned Sbox_depth = 2; // S盒深度

constexpr long mValues[][4] = {
    // {p, log2(m)Max, bits, c}
    {65537, 16, 500, 2},    // 0
    {16776961, 8, 500, 2},  // 1
    {33553153, 8, 500, 2},  // 2
    {16777729, 9, 500, 2},  // 3
    {33551873, 9, 500, 2},  // 4
    {16780289, 10, 500, 2}, // 5
    {33510401, 10, 500, 2}, // 6
    {16807937, 11, 500, 2}, // 7
    {33540097, 11, 500, 2}, // 8
    {16863233, 12, 500, 2}, // 9
    {33550337, 12, 500, 2}, // 10
    {16801793, 13, 500, 2}, // 11
    {33349633, 13, 500, 2}, // 12
    {16760833, 14, 500, 2}, // 13
    {33538049, 14, 500, 2}, // 14
    {17006593, 15, 500, 2}, // 15
    {33128449, 15, 500, 2}, // 16，建议实际值<=14
    {17367041, 16, 500, 2}, // 17,
    {32440321, 16, 500, 2}, // 18, 建议实际值<=14
}; // p=k*m+1

// 当电脑内存有限时，log2Para_m太大，会导致内存不足，出现terminate called recursively错误，从而终止程序.
// 此外，即使正常运行，由于内存一直临界，会导致程序运行速度变慢，时间测量不准确。
//!!!!!!!!!!!!!!!
constexpr long idx = 0;
constexpr long log2Para_m = mValues[idx][1] - 1;

constexpr long Para_p = mValues[idx][0];    // plaintext prime
constexpr long Para_m = 1 << log2Para_m;    // cyclotomic polynomial
constexpr long phi_m = Para_m >> 1;         // phi(m)=nlsots
constexpr long Para_bits = mValues[idx][2]; // bits in the ciphertext modulus chain
constexpr long Para_c = mValues[idx][3];    // columns in the key-switching matrix

//!!!!!!!!!!!!!!!
constexpr unsigned PlainBlock = phi_m - 0; // 明文分组数,应该PlainBlock<=phi_m

// 计算 log2 的 constexpr 函数
constexpr unsigned int log2_constexpr(unsigned long long n, unsigned int p = 0)
{
    return (n <= 1) ? p : log2_constexpr(n / 2, p + 1);
}
constexpr long PlainMod = Para_p;
constexpr unsigned Bytebits = log2_constexpr(PlainMod - 1) + 1;
; // 字节比特长度=ceil(log2(PlainMod-1))

//!!!!!!!!!!!!!!!!
constexpr long BlockByte = 32; // 分组字节长度

constexpr unsigned BlockSize = Bytebits * BlockByte;  // 分组比特长度=BlockByte*Bytebits
static const long PlainByte = BlockByte * PlainBlock; // 明文字节长度
static const long Plainbits = Bytebits * PlainByte;   // 明文比特长度
static const long PlainSize = BlockSize * PlainBlock; // 明文比特长度

static const unsigned NonceSize = 32;                           // Nonce比特长度
static const long counter_begin = 0;                            // 计数器起始值
static const long counter_end = PlainBlock + counter_begin - 1; // 计数器结束值
static bool Rkflag = 1;                                         // true/1表示乘法，false/0表示加法

YusP yusP(PlainMod);

void encodeTo32Ctxt(vector<ZZX> &encData, const vector<long> &data, const EncryptedArray &ea)
{
    long R = data.size() / PlainByte;
    long nCtxt = BlockByte * R;
    long data_length = data.size();
    encData.resize(nCtxt);

    for (long r = 0; r < R; r++)
    {
        for (long i = 0; i < BlockByte; i++)
        {
            vector<long> slots(ea.size(), 0);
            for (long j = 0; j < PlainBlock; j++)
            {
                long byteIdx = j * BlockByte + i + r * PlainByte;
                if (byteIdx < data_length)
                {
                    slots[j] = data[byteIdx];
                }
            }
            ea.encode(encData[r * BlockByte + i], slots);
        }
    }
}
// encodeTo32Ctxt对应的解码
void decodeTo32Ctxt(vector<long> &data, const vector<vector<long>> &encData,
                    const EncryptedArray &ea)
{
    long R = encData.size() / BlockByte;
    long data_length = R * PlainByte;

    data.resize(data_length);

    for (long r = 0; r < R; r++)
    {
        for (long j = 0; j < PlainBlock; j++)
        {
            for (long i = 0; i < BlockByte; i++)
            { // i is the ciphertext number

                // j is the block number in this ctxt
                long byteIdx = j * BlockByte + i + r * PlainByte;
                if (byteIdx < data_length)
                {
                    data[byteIdx] = encData[r * BlockByte + i][j];
                }
            }
        }
    }
}
// 函数：解密并验证密文是否正确，需要解码
// 函数：解密并验证密文是否正确
bool verifyDecryption32(const std::vector<Ctxt> &encryptedVec, const vector<long> &originalVec, const SecKey &secretKey,
                        const EncryptedArray &ea)
{
    auto start_decrypt = std::chrono::steady_clock::now();

    vector<long> decryptedVec;
    std::vector<std::vector<long>> decryptedPolys(encryptedVec.size());
#pragma omp parallel for
    for (std::size_t i = 0; i < encryptedVec.size(); ++i)
    {
        ea.decrypt(encryptedVec[i], secretKey, decryptedPolys[i]);
    }
    // 解码
    decodeTo32Ctxt(decryptedVec, decryptedPolys, ea);
    // 验证解密结果
    bool isDecryptedVecCorrect = std::equal(decryptedVec.begin(), decryptedVec.end(), originalVec.begin());
    auto end_decrypt = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification finished! Time: " << elapsed_seconds.count() << "s\n";
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
// A[id0] + A[id1] + t[2]   + A[id3] + t[4]   + A[id5] + t[6]   + t[7];
// t[0]   + A[id1]          + t[3]   + A[id4] + t[5]   + A[id6] + t[7];
// t[0]   + t[1]   + A[id2] + A[id3] + t[4]   + A[id5] + t[6]   + A[id7];
// A[id0] + t[1]   + t[2]   + A[id3]          + t[5]   + A[id6] + t[7];
// t[0]   + A[id1] + t[2]   + t[3]   + A[id4] + A[id5] + t[6]   + A[id7];
// A[id0] + t[1]   + A[id2] + t[3]   + t[4]   + A[id5]          + t[7];
// t[0]   + A[id1] + t[2]   + A[id3] + t[4]   + t[5]   + A[id6] + A[id7];
//          t[1]   + A[id2] + t[3]   + A[id4] + t[5]   + t[6]   + A[id7];
void HE_MC_MR(vector<Ctxt> &eData)
{
    vector<Ctxt> temp = eData;
    vector<int> index = {0, 1, 2, 3, 4, 5, 6, 7};
    vector<Ctxt> t = eData;
    for (int i = 0; i < 4; i++)
    {
        int s = 8 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s;
        // 下面是60次加法
        t[id0] += t[id0];
        t[id1] += t[id1];
        t[id2] += t[id2];
        t[id3] += t[id3];
        t[id4] += t[id4];
        t[id5] += t[id5];
        t[id6] += t[id6];
        t[id7] += t[id7];
        eData[id0] += temp[id1];
        eData[id0] += t[id2];
        eData[id0] += temp[id3];
        eData[id0] += t[id4];
        eData[id0] += temp[id5];
        eData[id0] += t[id6];
        eData[id0] += t[id7];
        eData[id1] += t[id0];
        eData[id1] += t[id3];
        eData[id1] += temp[id4];
        eData[id1] += t[id5];
        eData[id1] += temp[id6];
        eData[id1] += t[id7];
        eData[id2] += t[id0];
        eData[id2] += t[id1];
        eData[id2] += temp[id3];
        eData[id2] += t[id4];
        eData[id2] += temp[id5];
        eData[id2] += t[id6];
        eData[id2] += temp[id7];
        eData[id3] += temp[id0];
        eData[id3] += t[id1];
        eData[id3] += t[id2];
        eData[id3] += t[id5];
        eData[id3] += temp[id6];
        eData[id3] += t[id7];
        eData[id4] += t[id0];
        eData[id4] += temp[id1];
        eData[id4] += t[id2];
        eData[id4] += t[id3];
        eData[id4] += temp[id5];
        eData[id4] += t[id6];
        eData[id4] += temp[id7];
        eData[id5] += temp[id0];
        eData[id5] += t[id1];
        eData[id5] += temp[id2];
        eData[id5] += t[id3];
        eData[id5] += t[id4];
        eData[id5] += t[id7];
        eData[id6] += t[id0];
        eData[id6] += temp[id1];
        eData[id6] += t[id2];
        eData[id6] += temp[id3];
        eData[id6] += t[id4];
        eData[id6] += t[id5];
        eData[id6] += temp[id7];
        eData[id7] += t[id1];
        eData[id7] += temp[id2];
        eData[id7] += t[id3];
        eData[id7] += temp[id4];
        eData[id7] += t[id5];
        eData[id7] += t[id6];
    }
    temp = eData;
    index = {0, 1, 8, 9, 16, 17, 24, 25};
    t = eData;
    for (int i = 0; i < 4; i++)
    {
        int s = 2 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s;
        // 下面是60次加法
        t[id0] += t[id0];
        t[id1] += t[id1];
        t[id2] += t[id2];
        t[id3] += t[id3];
        t[id4] += t[id4];
        t[id5] += t[id5];
        t[id6] += t[id6];
        t[id7] += t[id7];
        eData[id0] += temp[id1];
        eData[id0] += t[id2];
        eData[id0] += temp[id3];
        eData[id0] += t[id4];
        eData[id0] += temp[id5];
        eData[id0] += t[id6];
        eData[id0] += t[id7];
        eData[id1] += t[id0];
        eData[id1] += t[id3];
        eData[id1] += temp[id4];
        eData[id1] += t[id5];
        eData[id1] += temp[id6];
        eData[id1] += t[id7];
        eData[id2] += t[id0];
        eData[id2] += t[id1];
        eData[id2] += temp[id3];
        eData[id2] += t[id4];
        eData[id2] += temp[id5];
        eData[id2] += t[id6];
        eData[id2] += temp[id7];
        eData[id3] += temp[id0];
        eData[id3] += t[id1];
        eData[id3] += t[id2];
        eData[id3] += t[id5];
        eData[id3] += temp[id6];
        eData[id3] += t[id7];
        eData[id4] += t[id0];
        eData[id4] += temp[id1];
        eData[id4] += t[id2];
        eData[id4] += t[id3];
        eData[id4] += temp[id5];
        eData[id4] += t[id6];
        eData[id4] += temp[id7];
        eData[id5] += temp[id0];
        eData[id5] += t[id1];
        eData[id5] += temp[id2];
        eData[id5] += t[id3];
        eData[id5] += t[id4];
        eData[id5] += t[id7];
        eData[id6] += t[id0];
        eData[id6] += temp[id1];
        eData[id6] += t[id2];
        eData[id6] += temp[id3];
        eData[id6] += t[id4];
        eData[id6] += t[id5];
        eData[id6] += temp[id7];
        eData[id7] += t[id1];
        eData[id7] += temp[id2];
        eData[id7] += t[id3];
        eData[id7] += temp[id4];
        eData[id7] += t[id5];
        eData[id7] += t[id6];
    }
}
// HE_MC_MR2有问题，需要修改
void HE_MC_MR2(vector<Ctxt> &eData)
{
    // auto start = std::chrono::steady_clock::now();
    // for (int i = 0; i < eData.size(); i++)
    // {
    //     eData[i].cleanUp();
    // }
    // auto end = std::chrono::steady_clock::now();
    // std::cout << "cleanUp finished! Time: " << std::chrono::duration<double>(end - start).count() << "s\n";
    vector<Ctxt> temp = eData;
    vector<int> index = {0, 1, 2, 3, 4, 5, 6, 7};
    Ctxt D = eData[0];
    Ctxt t06 = eData[0];
    Ctxt t15 = eData[0];
    Ctxt t24 = eData[0];
    Ctxt t14 = eData[0];
    Ctxt t05 = eData[0];
    Ctxt D7 = eData[0];
    Ctxt D3 = eData[0];
    Ctxt D73 = eData[0];
    for (int i = 0; i < 4; i++)
    {
        int s = 8 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s;
        // 下面是33次加/减法, 17次赋值
        D = temp[id0];           // 0 次
        D += temp[id1];          // 1 次
        D += temp[id2];          // 2 次
        D += temp[id3];          // 3 次
        D += temp[id4];          // 4 次
        D += temp[id5];          // 5 次
        D += temp[id6];          // 6 次
        D += temp[id7];          // 7 次
        D7 = D;                  // 7 次
        D7 += temp[id7];         // 8 次
        D73 = D7;                // 8 次
        D73 += temp[id3];        // 9 次
        D3 = D;                  // 9 次
        D3 += temp[id3];         // 10 次
        t06 = temp[id0];         // 10 次
        t06 += temp[id6];        // 11 次
        t15 = temp[id1];         // 11 次
        t15 += temp[id5];        // 12 次
        t24 = temp[id2];         // 12 次
        t24 += temp[id4];        // 13 次
        t14 = temp[id1];         // 13 次
        t14 += temp[id4];        // 14 次
        t05 = temp[id0];         // 14 次
        t05 += temp[id5];        // 15 次
        eData[id0] = D7;         // 15 次
        eData[id0] += t24;       // 16 次
        eData[id0] += temp[id6]; // 17 次
        eData[id1] = D73;        // 17 次
        eData[id1] += t05;       // 18 次
        eData[id1] -= temp[id2]; // 19 次
        eData[id2] = t06;        // 19 次
        eData[id2] += D;         // 20 次
        eData[id2] += t14;       // 21 次
        eData[id3] = D7;         // 21 次
        eData[id3] += t15;       // 22 次
        eData[id3] += temp[id2]; // 23 次
        eData[id3] -= temp[id4]; // 24 次
        eData[id4] = D3;         // 24 次
        eData[id4] += t06;       // 25 次
        eData[id4] += temp[id2]; // 26 次
        eData[id5] = D73;        // 26 次
        eData[id5] += t14;       // 27 次
        eData[id5] -= temp[id5]; // 28 次
        eData[id6] = D;          // 28 次
        eData[id6] += t24;       // 29 次
        eData[id6] += t05;       // 30 次
        eData[id7] = D3;         // 30 次
        eData[id7] += t15;       // 31 次
        eData[id7] -= temp[id0]; // 32 次
        eData[id7] += temp[id6]; // 33 次
    }
    temp = eData;
    index = {0, 1, 8, 9, 16, 17, 24, 25};
    D = eData[0];
    t06 = eData[0];
    t15 = eData[0];
    t24 = eData[0];
    t14 = eData[0];
    t05 = eData[0];
    D7 = eData[0];
    D3 = eData[0];
    D73 = eData[0];
    for (int i = 0; i < 4; i++)
    {
        int s = 2 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s;
        // 下面是33次加/减法, 17次赋值
        D = temp[id0];           // 0 次
        D += temp[id1];          // 1 次
        D += temp[id2];          // 2 次
        D += temp[id3];          // 3 次
        D += temp[id4];          // 4 次
        D += temp[id5];          // 5 次
        D += temp[id6];          // 6 次
        D += temp[id7];          // 7 次
        D7 = D;                  // 7 次
        D7 += temp[id7];         // 8 次
        D73 = D7;                // 8 次
        D73 += temp[id3];        // 9 次
        D3 = D;                  // 9 次
        D3 += temp[id3];         // 10 次
        t06 = temp[id0];         // 10 次
        t06 += temp[id6];        // 11 次
        t15 = temp[id1];         // 11 次
        t15 += temp[id5];        // 12 次
        t24 = temp[id2];         // 12 次
        t24 += temp[id4];        // 13 次
        t14 = temp[id1];         // 13 次
        t14 += temp[id4];        // 14 次
        t05 = temp[id0];         // 14 次
        t05 += temp[id5];        // 15 次
        eData[id0] = D7;         // 15 次
        eData[id0] += t24;       // 16 次
        eData[id0] += temp[id6]; // 17 次
        eData[id1] = D73;        // 17 次
        eData[id1] += t05;       // 18 次
        eData[id1] -= temp[id2]; // 19 次
        eData[id2] = t06;        // 19 次
        eData[id2] += D;         // 20 次
        eData[id2] += t14;       // 21 次
        eData[id3] = D7;         // 21 次
        eData[id3] += t15;       // 22 次
        eData[id3] += temp[id2]; // 23 次
        eData[id3] -= temp[id4]; // 24 次
        eData[id4] = D3;         // 24 次
        eData[id4] += t06;       // 25 次
        eData[id4] += temp[id2]; // 26 次
        eData[id5] = D73;        // 26 次
        eData[id5] += t14;       // 27 次
        eData[id5] -= temp[id5]; // 28 次
        eData[id6] = D;          // 28 次
        eData[id6] += t24;       // 29 次
        eData[id6] += t05;       // 30 次
        eData[id7] = D3;         // 30 次
        eData[id7] += t15;       // 31 次
        eData[id7] -= temp[id0]; // 32 次
        eData[id7] += temp[id6]; // 33 次
    }
}
// Compute the constants for Sbox
void HE_Sbox(vector<Ctxt> &eData)
{
    for (int i = 0; i < eData.size(); i++)
    {
        eData[i].cleanUp();
    }
    vector<Ctxt> c = eData;
    // #pragma omp parallel for
    for (long j = 0; j < BlockByte; j += 2)
    {
        for (int i = 0; i < 2; i++)
        {
            Ctxt temp = c[j];
            // c0.frobeniusAutomorph(1);
            c[j].square();
            c[j] += c[j + 1];
            c[j + 1] = temp;
        }
        // #pragma omp critical
        {
            eData[j] = c[j];
            eData[j + 1] = c[j + 1];
        }
    }
}
int main()
{
    std::cout << "Nr=" << Nr << std::endl;
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
                yusP.MC32_2(state); // 列混淆
                yusP.MR32_2(state); // 行移位
                yusP.Sbox(state);   // S盒
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
            }
            else
            {                       // 最后一轮
                yusP.MC32_2(state); // 列混淆
                yusP.MR32_2(state); // 行移位
                yusP.Sbox(state);   // S盒
                yusP.MC32_2(state); // 列混淆
                yusP.MR32_2(state); // 行移位
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

    auto start_keyEncryption = std::chrono::steady_clock::now();
    vector<Ctxt> encryptedSymKey;
    encryptSymKey(encryptedSymKey, SymKey, publicKey, ea);
    auto end_keyEncryption = std::chrono::steady_clock::now();
    double keyEncryption = std::chrono::duration<double>(end_keyEncryption - start_keyEncryption).count();
    std::cout << "SymKey FHE time: " << keyEncryption << "s\n";

    // 解密验证
    if (!verify_encryptSymKey(encryptedSymKey, SymKey, secretKey, ea))
    {
        return 0;
    }
    std::cout << "Symmetric key encryption succeeded!" << std::endl;

    // 离线客户端时间=KeyStream Generation time+PublicKey generation and SymKey FHE time
    double total_time_off = elapsed_seconds_keyStream.count() + elapsed_seconds_PubKey.count() + elapsed_seconds_PubKey.count() + keyEncryption;
    std::cout << "Encryption offline total time: " << total_time_off << "s\n";
    //=============服务端offline阶段================
    // 计算 encryptedRoundKeySet
    auto start_RoundKeySet_FHE = std::chrono::steady_clock::now();
    vector<Ctxt> encryptedRoundKeySet;
    Ctxt tmpCtxt(*publicKey);
    long eRk_len = BlockByte * (Nr + 1);
    encryptedRoundKeySet.resize(eRk_len, tmpCtxt);

    for (int i = 0; i < eRk_len; i++)
    {
        encryptedRoundKeySet[i] = encryptedSymKey[i % BlockByte];
    }

    vector<ZZX> encodedXset;
    encodeTo32Ctxt(encodedXset, Xset, ea); // encode as HE plaintext
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
    double RoundKeySet_FHE_time = std::chrono::duration<double>(end_RoundKeySet_FHE - start_RoundKeySet_FHE).count();
    std::cout << "RoundKeySet FHE succeeded! Time: " << RoundKeySet_FHE_time << "s\n";
    // // 使用 verifyDecryption 函数解密并验证 RoundKeySet
    // if (!verifyDecryption32(encryptedRoundKeySet, RoundKeySet, secretKey, ea))
    // {
    //     std::cerr << "Decryption verification failed for RoundKeySet." << std::endl;
    //     return 0;
    // }
    // std::cout << "Decryption verification succeeded for RoundKeySet." << std::endl;

    // 生成 encryptedKeyStream
    // 定义roundkey_time、sbox_time、linear_layer_time
    double sbox_time = 0, linear_layer_time = 0, roundkey_time = 0;

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
    encodeTo32Ctxt(encoded_expandedIV, expandedIV, ea); // encode as HE plaintext

    std::cout << "whiteround start" << std::endl;
    auto start_roundkey = std::chrono::high_resolution_clock::now();
    for (long i = 0; i < BlockByte; i++)
    { // encrypt the encoded key
        encryptedKeyStream[i].addConstant(encoded_expandedIV[i]);
    }
    auto end_roundkey = std::chrono::high_resolution_clock::now();
    roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
    // 输出 roundkey_time
    std::cout << "whiteround time: " << roundkey_time << "s\n";
    // 明文密钥流
    vector<long> KeyStream2(PlainByte);
    // 对IV和RoundKeySet进行异或
    for (long i = 0; i < PlainByte; i++)
    {
        KeyStream2[i] = (expandedIV[i] + RoundKeySet[i]) % PlainMod;
    }
    // 使用 verifyDecryption 函数解密并验证 KeyStream
    if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for KeyStream." << std::endl;
        return 0;
    }
    std::cout << "Decryption verification succeeded for whiteround." << std::endl;

    auto start_sbox = std::chrono::high_resolution_clock::now();
    auto start_linear = std::chrono::high_resolution_clock::now();
    auto end_sbox = std::chrono::high_resolution_clock::now();
    auto end_linear = std::chrono::high_resolution_clock::now();

    for (long r = 1; r < Nr; r++)
    {
        std::cout << "Round " << r << " start" << std::endl;
        start_linear = std::chrono::high_resolution_clock::now();
        // #pragma omp parallel for
        // MC Layer + MR Layer
        HE_MC_MR(encryptedKeyStream);
        end_linear = std::chrono::high_resolution_clock::now();
        linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
        for (int i = 0; i < PlainBlock; i++)
        {
            vector<long> tmp(BlockByte);
            for (int j = 0; j < BlockByte; j++)
            {
                tmp[j] = KeyStream2[i * BlockByte + j];
            }
            yusP.MC32_2(tmp);
            yusP.MR32_2(tmp);
            for (int j = 0; j < BlockByte; j++)
            {
                KeyStream2[i * BlockByte + j] = tmp[j];
            }
        }
        if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
        start_sbox = std::chrono::high_resolution_clock::now();
        // S Layer
        HE_Sbox(encryptedKeyStream);
        end_sbox = std::chrono::high_resolution_clock::now();
        sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();

        yusP.Sbox(KeyStream2);
        if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
        start_roundkey = std::chrono::high_resolution_clock::now();
        // #pragma omp parallel for
        for (long j = 0; j < BlockByte; j++)
        {
            encryptedKeyStream[j] += encryptedRoundKeySet[r * BlockByte + j];
        }
        end_roundkey = std::chrono::high_resolution_clock::now();
        roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
        for (long i = 0; i < PlainByte; i++)
        {
            KeyStream2[i] = (KeyStream2[i] + RoundKeySet[r * PlainByte + i]) % PlainMod;
        }
        if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;
    }

    // 最后一轮
    std::cout << "Round " << Nr << " start" << std::endl;
    start_linear = std::chrono::high_resolution_clock::now();
    // MC Layer + MR Layer
    HE_MC_MR(encryptedKeyStream);
    end_linear = std::chrono::high_resolution_clock::now();
    linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
    // for (int i = 0; i < PlainBlock; i++)
    // {
    //     vector<long> tmp(BlockByte);
    //     for (int j = 0; j < BlockByte; j++)
    //     {
    //         tmp[j] = KeyStream2[i * BlockByte + j];
    //     }
    //     yusP.MC32_2(tmp);
    //     yusP.MR32_2(tmp);
    //     for (int j = 0; j < BlockByte; j++)
    //     {
    //         KeyStream2[i * BlockByte + j] = tmp[j];
    //     }
    // }
    // if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
    // {
    //     std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
    //     return 0;
    // }
    // std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
    start_sbox = std::chrono::high_resolution_clock::now();
    // S Layer
    HE_Sbox(encryptedKeyStream);
    end_sbox = std::chrono::high_resolution_clock::now();
    sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
    // yusP.Sbox(KeyStream2);
    // if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
    // {
    //     std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
    //     return 0;
    // }
    // std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
    start_linear = std::chrono::high_resolution_clock::now();
    // MC Layer + MR Layer
    HE_MC_MR(encryptedKeyStream);
    end_linear = std::chrono::high_resolution_clock::now();
    linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
    // for (int i = 0; i < PlainBlock; i++)
    // {
    //     vector<long> tmp(BlockByte);
    //     for (int j = 0; j < BlockByte; j++)
    //     {
    //         tmp[j] = KeyStream2[i * BlockByte + j];
    //     }
    //     yusP.MC32_2(tmp);
    //     yusP.MR32_2(tmp);
    //     for (int j = 0; j < BlockByte; j++)
    //     {
    //         KeyStream2[i * BlockByte + j] = tmp[j];
    //     }
    // }
    // if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
    // {
    //     std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
    //     return 0;
    // }
    // std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
    // add
    start_roundkey = std::chrono::high_resolution_clock::now();
    for (long j = 0; j < BlockByte; j++)
    {
        encryptedKeyStream[j] += encryptedRoundKeySet[Nr * BlockByte + j];
    }
    end_roundkey = std::chrono::high_resolution_clock::now();
    roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
    // for (long i = 0; i < PlainByte; i++)
    // {
    //     KeyStream2[i] = (KeyStream2[i] + RoundKeySet[Nr * PlainByte + i]) % PlainMod;
    // }
    // if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
    // {
    //     std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
    //     return 0;
    // }
    // std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;

    // 输出 roundkey_time、sbox_time、linear_layer_time
    std::cout << "RoundKey time: " << roundkey_time << "s\n";
    std::cout << "Sbox time: " << sbox_time << "s\n";
    std::cout << "Linear Layer time: " << linear_layer_time << "s\n";
    // 计算总时间
    double total_time = roundkey_time + sbox_time + linear_layer_time + RoundKeySet_FHE_time;
    std::cout << "Server offline total time: " << total_time << "s\n";
    // for (int i = 0; i < encryptedKeyStream.size(); i++)
    // {
    //     encryptedKeyStream[i].bringToSet(encryptedKeyStream[i].naturalPrimeSet());
    // }
    // 计算吞吐量,KB/min
    double throughput = (Plainbits * 60) / (pow(2, 13) * total_time);
    std::cout << "Throughput: " << throughput << "KB/min\n";
    if (!verifyDecryption32(encryptedKeyStream, KeyStream, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for KeyStream." << std::endl;
        return 0;
    }
    std::cout << "Decryption verification succeeded for KeyStream." << std::endl;
    // 将total_time, throughput, Nr, p, nslots, bits, c, roundkey_time, sbox_time, linear_layer_time, RoundKeySet_FHE_time写入文件test_Yus_p_C32_ClientAndServer2.txt,如果已存在则追加
    // 指定文件路径
    std::string filePath = "../tests/test_Yus_p_C32_ClientAndServer2.txt";
    std::ofstream outfile(filePath, std::ios::app);
    if (!outfile)
    {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return 0;
    }
    outfile << std::left << std::setw(20) << total_time
            << std::left << std::setw(20) << throughput
            << std::left << std::setw(10) << Nr
            << std::left << std::setw(20) << p
            << std::left << std::setw(20) << nslots
            << std::left << std::setw(15) << bits
            << std::left << std::setw(10) << c
            << std::left << std::setw(20) << roundkey_time
            << std::left << std::setw(20) << sbox_time
            << std::left << std::setw(20) << linear_layer_time
            << std::left << std::setw(20) << RoundKeySet_FHE_time
            << std::endl;
    outfile.close();
    std::cout << "test_Yus_p_C32_ClientAndServer2.txt updated." << std::endl;
    return 0;
}