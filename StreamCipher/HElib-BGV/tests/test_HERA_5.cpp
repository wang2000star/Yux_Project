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

#include "../utils/random_bit.hpp"
#include "../Symmetric/HERA.hpp"
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
constexpr unsigned BlockWords = 16;      // 分组密钥字长度=KeyWords
constexpr unsigned BlockPlainWords = 16; // 明文分组字长度
constexpr double TruncRate = BlockPlainWords / (double)BlockWords;
// ===============模式设置================
constexpr bool Rkflag = 1;     // true/1表示乘法，false/0表示加法，指的是随机向量和密钥间的操作
constexpr bool deflag = 0;     // true/1表示进行每一步解密验证，false/0表示不进行每一步解密验证
constexpr bool ompflag = 0;    // true/1表示使用OpenMP并行编码，false/0表示不使用OpenMP并行编码
constexpr bool symkeyflag = 0; // true/1表示对称密钥同态解密验证加密，false/0表示不验证
constexpr bool KeyStreamflag = 0;  // true/1表示密钥流同态解密验证，false/0表示不验证
constexpr bool plainflag = 0;      // true/1表示明文同态解密验证，false/0表示不验证
// 参数设置，paramMap[Nr-4][idx]
constexpr unsigned Nr = 5; // 轮数
constexpr long idx = 0;
constexpr unsigned Sbox_depth = 1 * Nr; // S盒深度
// 当c=2时，Qbits=1.5*bits,当c=3时，Qbits=1.5*bits - 100
// 硬编码参数值
constexpr tuple<long, long, long, long> paramMap[1][8] = {
    {
        // Nr = 5
        // {p, log2(m), bits, c}
        {65537, 16, 490, 2},  // 0 *
        {163841, 15, 280, 2}, // 1
        {65537, 16, 300, 2},  // 2
        {65537, 16, 350, 2},  // 3
        {0, 0, 0, 0},         // 填充空位
        {0, 0, 0, 0},         // 填充空位
        {0, 0, 0, 0},         // 填充空位
        {0, 0, 0, 0}          // 填充空位
    }};
// p=k*m+1
//  2^10=1024,2^11=2048,2^12=4096,2^13=8192,2^14=16384,2^15=32768,2^16=65536,
// 当电脑内存有限时，log2Para_m太大，会导致内存不足，出现terminate called recursively错误，从而终止程序.
// 此外，即使正常运行，由于内存一直临界，会导致程序运行速度变慢，时间测量不准确。

constexpr long log2Para_m = get<1>(paramMap[0][idx]) - 0;
constexpr long Para_p = get<0>(paramMap[0][idx]);    // plaintext prime
constexpr long Para_m = 1 << log2Para_m;                  // cyclotomic polynomial
constexpr long phi_m = Para_m >> 1;                       // phi(m)
constexpr long Para_bits = get<2>(paramMap[0][idx]); // bits in the ciphertext modulus chain
constexpr long Para_c = get<3>(paramMap[0][idx]);    // columns in the key-switching matrix
constexpr long Para_r = 1;                                // Lifting [defualt = 1]
//!!!!!!!!!!!!!!!
constexpr long nslots = phi_m;             // 槽数
constexpr unsigned PlainBlock = nslots - 0; // 明文分组数,应该PlainBlock<=phi_m
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

HERA hera(PlainMod); // 构建明文对称加密实例

void HE_MC(vector<Ctxt> &eData)
{
    vector<Ctxt> temp = eData;
    array<int, 8> index = {0, 1, 2, 3};
    Ctxt T2(eData[0]);
    Ctxt T3(eData[0]);
    for (int i = 0; i < 4; i++)
    {
        int s = 4 * i;
        T2 = temp[index[0] + s];
        T2 += temp[index[0] + s];
        T3 = temp[index[1] + s];
        T3 += temp[index[1] + s];
        T3 += temp[index[1] + s];
        eData[index[0] + s] = T2;
        eData[index[0] + s] += T3;
        eData[index[0] + s] += temp[index[2] + s];
        eData[index[0] + s] += temp[index[3] + s];

        T2 = temp[index[1] + s];
        T2 += temp[index[1] + s];
        T3 = temp[index[2] + s];
        T3 += temp[index[2] + s];
        T3 += temp[index[2] + s];
        eData[index[1] + s] = T2;
        eData[index[1] + s] += T3;
        eData[index[1] + s] += temp[index[3] + s];
        eData[index[1] + s] += temp[index[0] + s];

        T2 = temp[index[2] + s];
        T2 += temp[index[2] + s];
        T3 = temp[index[3] + s];
        T3 += temp[index[3] + s];
        T3 += temp[index[3] + s];
        eData[index[2] + s] = T2;
        eData[index[2] + s] += T3;
        eData[index[2] + s] += temp[index[0] + s];
        eData[index[2] + s] += temp[index[1] + s];

        T2 = temp[index[3] + s];
        T2 += temp[index[3] + s];
        T3 = temp[index[0] + s];
        T3 += temp[index[0] + s];
        T3 += temp[index[0] + s];
        eData[index[3] + s] = T2;
        eData[index[3] + s] += T3;
        eData[index[3] + s] += temp[index[1] + s];
        eData[index[3] + s] += temp[index[2] + s];
    }
}
void HE_MR(vector<Ctxt> &eData)
{
    vector<Ctxt> temp = eData;
    vector<int> index = {0, 4, 8, 12};
    Ctxt T2(eData[0]);
    Ctxt T3(eData[0]);
    for (int i = 0; i < 4; i++)
    {
        int s = i;
        T2 = temp[index[0] + s];
        T2 += temp[index[0] + s];
        T3 = temp[index[1] + s];
        T3 += temp[index[1] + s];
        T3 += temp[index[1] + s];
        eData[index[0] + s] = T2;
        eData[index[0] + s] += T3;
        eData[index[0] + s] += temp[index[2] + s];
        eData[index[0] + s] += temp[index[3] + s];

        T2 = temp[index[1] + s];
        T2 += temp[index[1] + s];
        T3 = temp[index[2] + s];
        T3 += temp[index[2] + s];
        T3 += temp[index[2] + s];
        eData[index[1] + s] = T2;
        eData[index[1] + s] += T3;
        eData[index[1] + s] += temp[index[3] + s];
        eData[index[1] + s] += temp[index[0] + s];

        T2 = temp[index[2] + s];
        T2 += temp[index[2] + s];
        T3 = temp[index[3] + s];
        T3 += temp[index[3] + s];
        T3 += temp[index[3] + s];
        eData[index[2] + s] = T2;
        eData[index[2] + s] += T3;
        eData[index[2] + s] += temp[index[0] + s];
        eData[index[2] + s] += temp[index[1] + s];

        T2 = temp[index[3] + s];
        T2 += temp[index[3] + s];
        T3 = temp[index[0] + s];
        T3 += temp[index[0] + s];
        T3 += temp[index[0] + s];
        eData[index[3] + s] = T2;
        eData[index[3] + s] += T3;
        eData[index[3] + s] += temp[index[1] + s];
        eData[index[3] + s] += temp[index[2] + s];
    }
}
// Compute the constants for Sbox
void HE_Sbox(vector<Ctxt> &eData)
{
    // #pragma omp parallel for
    for (long j = 0; j < BlockWords; j++)
    {
        eData[j].cube();
    }
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
    GF2X rnd;
    int Bytebitsdiv8 = ceil(Wordbits / 8);
    vector<uint8_t> SymKey0(Bytebitsdiv8 * BlockWords);
    random(rnd, 8 * SymKey0.size());
    BytesFromGF2X(SymKey0.data(), rnd, SymKey0.size());
    vector<long> SymKey(BlockWords);
    for (unsigned i = 0; i < BlockWords; i++)
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
    long root = 3; // m-th root of unity modulo p
    if (Para_m == 131072)
    {
        root = 3;
    }
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
    encryptSymKey(encryptedSymKey, SymKey, publicKey, cmodulus,nslots);
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
            // memcpy(&RoundKeySet[KeyStreamWords * r + BlockWords * (block_num)], RoundKey.data(), BlockWords * sizeof(long));
            if (r == 0)
            { // 初始轮
                for (unsigned i = 0; i < BlockWords; i++)
                {
                    state[i] = (RoundKey[i] + IV[i]) % PlainMod;
                }
            }
            else if (r < Nr)
            { // 常规轮

                hera.MC(state);  // MixColumns
                hera.MR(state);  // MixRows
                hera.Sbox(state); // S盒
                for (unsigned i = 0; i < BlockWords; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
            }
            else
            {                      // 最后一轮
                hera.MC(state); // MixColumns
                hera.MR(state); // MixRows
                hera.Sbox(state); // S盒
                hera.MC(state);  // MixColumns
                hera.MR(state);  // MixRows
                for (unsigned i = 0; i < BlockWords; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
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

    vector<zzX> encodedXset(NrBlockWords);
    zz_pX encodedtemp;
    auto start_Xset = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NrBlockWords; i++)
    {
        cmodulus.iFFT(encodedtemp, Xset[i]);
        convert(encodedXset[i], encodedtemp);
    }
    auto end_Xset = std::chrono::high_resolution_clock::now();
    double Encode_time = std::chrono::duration<double>(end_Xset - start_Xset).count();
    std::cout << "encode time: " << Encode_time << "s\n";
    std::cout << "encode time/slots:" << Encode_time / (NrBlockWords + 2 * len3) * pow(10, 6) << "us\n";

    Ctxt tmpCtxt(*publicKey);
    int noise_budget = min_noise_budget(encryptedSymKey);
    std::cout << "noise budget initially: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    // 计算 encryptedRoundKeySet
    long eRk_len = BlockWords * (Nr + 1);
    vector<Ctxt> encryptedRoundKeySet(eRk_len, tmpCtxt);
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
    //     if (!verifyDecryption(encryptedRoundKeySet, RoundKeySet, secretKey, cmodulus,BlockWords,PlainBlock,nlsots))
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
    vector<Ctxt> encryptedKeyStream(encryptedRoundKeySet.begin(), encryptedRoundKeySet.begin() + BlockWords);
    std::cout << "whiteround start" << std::endl;
    start_add = std::chrono::high_resolution_clock::now();
    for (long i = 0; i < BlockWords; i++)
    { // encrypt the encoded key
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
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords,PlainBlock,nslots,Para_p))
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
        HE_MC(encryptedKeyStream);
        HE_MR(encryptedKeyStream);
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
                vector<long> tmp(BlockWords);
                for (int j = 0; j < BlockWords; j++)
                {
                    tmp[j] = KeyStream2[i * BlockWords + j];
                }
                hera.MC(tmp);
                hera.MR(tmp);
                for (int j = 0; j < BlockWords; j++)
                {
                    KeyStream2[i * BlockWords + j] = tmp[j];
                }
            }
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords,PlainBlock,nslots,Para_p))
            {
                std::cerr << "Decryption verification failed for KeyStream2 Linear Layer." << std::endl;
                return 0;
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
            return 0;
        }
        if (deflag)
        {
            hera.Sbox(KeyStream2);
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords,PlainBlock,nslots,Para_p))
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
            encryptedKeyStream[j] += encryptedRoundKeySet[r * BlockWords + j];
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
            for (long i = 0; i < KeyStreamWords; i++)
            {
                KeyStream2[i] = (KeyStream2[i] + RoundKeySet[r * KeyStreamWords + i]) % PlainMod;
            }
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords,PlainBlock,nslots,Para_p))
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
    HE_MC(encryptedKeyStream);
    HE_MR(encryptedKeyStream);
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
            vector<long> tmp(BlockWords);
            for (int j = 0; j < BlockWords; j++)
            {
                tmp[j] = KeyStream2[i * BlockWords + j];
            }
            hera.MC(tmp);
            hera.MR(tmp);
            for (int j = 0; j < BlockWords; j++)
            {
                KeyStream2[i * BlockWords + j] = tmp[j];
            }
        }
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords,PlainBlock,nslots,Para_p))
        {
            std::cerr << "Decryption verification failed for KeyStream2 Linear Layer." << std::endl;
            return 0;
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
        return 0;
    }
    if (deflag)
    {
        hera.Sbox(KeyStream2);
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords,PlainBlock,nslots,Para_p))
        {
            std::cerr << "Decryption verification failed for KeyStream2 Sbox." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream2 Sbox." << std::endl;
    }
    start_linear = std::chrono::high_resolution_clock::now();
    // Linear Layer
    HE_MC(encryptedKeyStream);
    HE_MR(encryptedKeyStream);
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
            vector<long> tmp(BlockWords);
            for (int j = 0; j < BlockWords; j++)
            {
                tmp[j] = KeyStream2[i * BlockWords + j];
            }
            hera.MC(tmp);
            hera.MR(tmp);
            for (int j = 0; j < BlockWords; j++)
            {
                KeyStream2[i * BlockWords + j] = tmp[j];
            }
        }
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords,PlainBlock,nslots,Para_p))
        {
            std::cerr << "Decryption verification failed for KeyStream2 Linear Layer." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream2 Linear Layer." << std::endl;
    }
    // add
    start_add = std::chrono::high_resolution_clock::now();
    for (long j = 0; j < BlockWords; j++)
    {
        encryptedKeyStream[j] += encryptedRoundKeySet[Nr * BlockWords + j];
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
        for (long i = 0; i < KeyStreamWords; i++)
        {
            KeyStream2[i] = (KeyStream2[i] + RoundKeySet[Nr * KeyStreamWords + i]) % PlainMod;
        }
        if (!verifyDecryption(encryptedKeyStream, KeyStream2, secretKey, cmodulus, BlockWords,PlainBlock,nslots,Para_p))
        {
            std::cerr << "Decryption verification failed for KeyStream2 Round Key Addition." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream2 Round Key Addition." << std::endl;
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
    if (KeyStreamflag)
    {
        if (!verifyDecryption(encryptedKeyStream, KeyStream, secretKey, cmodulus, BlockPlainWords,PlainBlock,nslots,Para_p))
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
    vector<vec_long> CipherStream(BlockPlainWords);
    for (int i = 0; i < BlockPlainWords; i++)
    {
        CipherStream[i].SetLength(nslots);
                for (int j = 0; j < nslots; j++)
        {
            CipherStream[i][j] = 0;
        }
    }
    auto start_ClientOnline = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < PlainBlock; j++)
    {
        for (int i = 0; i < BlockPlainWords; i++)
        {
            CipherStream[i][j] = ((PlainStream[j * BlockPlainWords + i] + KeyStream[j * BlockPlainWords + i]) % PlainMod + PlainMod) % PlainMod;
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
    vector<Ctxt> encrypedPlainStream = encryptedKeyStream;
    // 对CipherStream进行编码
    vector<ZZX> encodedCipherStream(BlockPlainWords);
    auto start_ServerOnline = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < BlockPlainWords; i++)
    {
        cmodulus.iFFT(encodedtemp, CipherStream[i]);
        convert(encodedCipherStream[i], encodedtemp);
    }
    for (int i = 0; i < BlockPlainWords; i++)
    {
        encrypedPlainStream[i].negate();
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
    if (plainflag)
    {
        if (!verifyDecryption(encrypedPlainStream, PlainStream, secretKey, cmodulus, BlockPlainWords,PlainBlock,nslots,Para_p))
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
        filePath = "test_HERA_5.txt";
    }
    else
    {
        filePath = "../tests/test_HERA_5.txt";
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
    std::cout << "test_HERA_5.txt updated." << std::endl;
    }
    return 0;
}