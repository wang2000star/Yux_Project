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
#include "../Symmetric/Yux_p.hpp"
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
constexpr bool Rkflag = 1;        // true/1表示乘法，false/0表示加法，指的是随机向量和密钥间的操作
constexpr bool deflag = 0;        // true/1表示进行每一步解密验证，false/0表示不进行每一步解密验证
constexpr bool ompflag = 0;       // true/1表示使用OpenMP并行编码，false/0表示不使用OpenMP并行编码
constexpr bool symkeyflag = 0;    // true/1表示对称密钥同态解密验证加密，false/0表示不验证
constexpr bool KeyStreamflag = 0; // true/1表示密钥流同态解密验证，false/0表示不验证
constexpr bool plainflag = 0;     // true/1表示明文同态解密验证，false/0表示不验证
// 参数设置，paramMap[Nr-4][idx]
constexpr unsigned Nr = 9; // 轮数
constexpr long idx = 0;
constexpr unsigned Sbox_depth = 1 * Nr; // S盒深度
// 当c=2时，Qbits=1.5*bits,当c=3时，Qbits=1.5*bits - 100
// 硬编码参数值
constexpr tuple<long, long, long, long> paramMap[1][8] = {
    {// Nr = 
     // {p, log2(m), bits, c}
     {4298506241, 16, 850, 20},  // 0 *
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
constexpr long RoundConstant = 0xCD;

YuxP yuxP(PlainMod, RoundConstant); // YuxP类实例化

// 线性层HE计算
void HE_M(vector<Ctxt> &eData)
{
    // for (long j = 0; j < eData.size(); j++)
    // {
    //     eData[j].cleanUp();
    // }
    vector<Ctxt> in = eData;
    for (long j = 0; j < eData.size(); j++)
    {
        eData[j] += in[(j + 3) % 16];
        eData[j] += in[(j + 4) % 16];
        eData[j] += in[(j + 8) % 16];
        eData[j] += in[(j + 9) % 16];
        eData[j] += in[(j + 12) % 16];
        eData[j] += in[(j + 14) % 16];
    }
    // for (long j = 0; j < eData.size(); j++)
    // {
    //     eData[j].cleanUp();
    // }
}
// S盒HE计算
void HE_Sbox(vector<Ctxt> &eData)
{
    // for (int i = 0; i < eData.size(); i++)
    // {
    //     eData[i].cleanUp();
    // }
    for (int i = 0; i < 16; i += 4)
    {
        for (int j = 0; j < 2; j++)
        {
            Ctxt c0(eData[i]);
            Ctxt c1(eData[i + 1]);
            Ctxt c2(eData[i + 2]);
            Ctxt c3(eData[i + 3]);

            Ctxt temp2(c1);
            temp2.multiplyBy(c2);
            temp2 += c0;
            temp2 += c3;
            temp2.addConstant(RoundConstant);

            Ctxt temp3(c2);
            temp3.multiplyBy(c3);
            temp3 += c1;
            temp3 += temp2;
            temp3.addConstant(RoundConstant);

            eData[i] = c2;
            eData[i + 1] = c3;
            eData[i + 2] = temp2;
            eData[i + 3] = temp3;
        }
    }
    // for (int i = 0; i < eData.size(); i++)
    // {
    //     eData[i].cleanUp();
    // }
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
    vector<long> expandSymKey(NrBlockWords);
    yuxP.KeyExpansion(SymKey, expandSymKey, Nr, BlockWords);
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
    encryptSymKey(encryptedSymKey, expandSymKey, publicKey, cmodulus, nslots);
    auto end_keyEncryption = std::chrono::high_resolution_clock::now();
    double keyEncryption = std::chrono::duration<double>(end_keyEncryption - start_keyEncryption).count();
    std::cout << "expandSymKey FHE time: " << keyEncryption << "s\n";
    //  解密验证
    if (symkeyflag)
    {
        if (!verify_encryptSymKey(encryptedSymKey, expandSymKey, secretKey, cmodulus))
        {
            return 0;
        }
        std::cout << "Symmetric key encryption succeeded!" << std::endl;
    }
    // 生成随机对称明文流，只用于测试
    random_device rd;
    vector<long> PlainStream(PlainWords);
    for (int i = 0; i < PlainWords; i++)
    {
        PlainStream[i] = rd() % PlainMod;
    }
    // Generating Cipher stream
    auto start_CipherStream = std::chrono::high_resolution_clock::now();

    std::vector<vec_long> CipherStream(BlockWords);
    for (unsigned i = 0; i < BlockWords; i++)
    {
        CipherStream[i].SetLength(PlainBlock);
    }
    long block_num;
    long ir;
    std::vector<long> state(BlockWords); // 初始化 state
    std::cout << "Generating CipherStream..." << std::endl;
    uint64_t start_cycle = rdtsc();
    for (long counter = counter_begin; counter <= counter_end; counter++)
    {
        block_num = counter - counter_begin;
        for (unsigned i = 0; i < BlockWords; i++)
        {
            state[i] = (PlainStream[block_num * BlockWords + i] + expandSymKey[i]) % PlainMod;
        }
        // 逐轮进行加密
        for (unsigned r = 1; r <= Nr; r++)
        {
            yuxP.Sbox(state);   // S盒
            yuxP.Linear(state); // 线性变换
            for (unsigned i = 0; i < BlockWords; i++)
            {
                state[i] = (state[i] + expandSymKey[r * BlockWords + i]) % PlainMod;
            }
        }
        for (unsigned i = 0; i < BlockWords; i++)
        {
            CipherStream[i][block_num] = state[i];
        }
    }
    uint64_t end_cycle = rdtsc();
    auto end_cipherStream = std::chrono::high_resolution_clock::now();
    double Client_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_cipherStream - start_CipherStream).count();
    std::cout << "Encryption  total time: " << Client_time << "s\n";
    uint64_t cycle_count = end_cycle - start_cycle;
    std::cout << "Encryption total cycles: " << cycle_count << std::endl;
    // 将KeyStream同态密文写入文件
    // if (!writeEncryptedSymKey(encryptedSymKey, "Client_encryptedSymKey.bin"))
    // {
    //     return false;
    // }
    // std::cout << "Client_encryptedSymKey.bin has been written to file." << std::endl;

    //=============服务端offline阶段================
            for (int test = 0; test < 10; test++)
    {
    // 对CipherStream进行同态加密
    Ctxt tmpCtxt(*publicKey);
    vector<Ctxt> encryptedCipherStream(BlockWords, tmpCtxt);
    zz_pX encodedtemp;
    zzX encodedtemp2;
    auto start_EncryptCipher = std::chrono::high_resolution_clock::now();
    for (long i = 0; i < BlockWords; i++)
    { // encrypt the encoded key
        cmodulus.iFFT(encodedtemp, CipherStream[i]);
        convert(encodedtemp2, encodedtemp);
        publicKey->Encrypt(encryptedCipherStream[i], encodedtemp2);
    }
    auto end_EncryptCipher = std::chrono::high_resolution_clock::now();
    double EncryptCipher_time = std::chrono::duration<double>(end_EncryptCipher - start_EncryptCipher).count();
    std::cout << "EncryptCipher time: " << EncryptCipher_time << "s\n";
    int noise_budget = min_noise_budget(encryptedCipherStream);
    std::cout << "noise budget initially: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    // 生成 encryptedPlainStream
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

    vector<Ctxt> encryptedPlainStream = encryptedCipherStream;
    std::cout << "whiteround start" << std::endl;
    start_add = std::chrono::high_resolution_clock::now();
    for (long i = 0; i < BlockWords; i++)
    {
        encryptedPlainStream[i] += encryptedSymKey[Nr * BlockWords + i];
    }
    end_add = std::chrono::high_resolution_clock::now();
    Add_time += std::chrono::duration<double>(end_add - start_add).count();
    // 输出 Add_time
    std::cout << "whiteround time: " << Add_time << "s\n";
    noise_budget = min_noise_budget(encryptedPlainStream);
    std::cout << "noise budget after Whiteround: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    for (long r = Nr - 1; r >= 0; r--)
    {
        std::cout << "Round " << r << " start" << std::endl;
        start_linear = std::chrono::high_resolution_clock::now();
        // Linear Layer
        HE_M(encryptedPlainStream);
        end_linear = std::chrono::high_resolution_clock::now();
        Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
        linear_set[r] = std::chrono::duration<double>(end_linear - start_linear).count();
        noise_budget = min_noise_budget(encryptedPlainStream);
        std::cout << "noise budget after linear: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            return 0;
        }
        start_sbox = std::chrono::high_resolution_clock::now();
        // S Layer
        HE_Sbox(encryptedPlainStream);
        end_sbox = std::chrono::high_resolution_clock::now();
        Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
        sbox_set[r] = std::chrono::duration<double>(end_sbox - start_sbox).count();
        noise_budget = min_noise_budget(encryptedPlainStream);
        std::cout << "noise budget after sbox: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            return 0;
        }
        // Round Key Addition
        start_add = std::chrono::high_resolution_clock::now();
        // omp_set_num_threads(12); // 设置线程数为12
        // #pragma omp parallel for
        for (long j = 0; j < BlockWords; j++)
        {
            encryptedPlainStream[j] += encryptedSymKey[r * BlockWords + j];
        }
        end_add = std::chrono::high_resolution_clock::now();
        Add_time += std::chrono::duration<double>(end_add - start_add).count();
        noise_budget = min_noise_budget(encryptedPlainStream);
        std::cout << "noise budget after Add: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            return 0;
        }
    }
    // 输出 Add_time、Sbox_time、Linear_time
    std::cout << "RoundKey time: " << EncryptCipher_time << "s\n";
    std::cout << "Add time: " << Add_time << "s\n";
    std::cout << "Sbox time: " << Sbox_time << "s\n";
    std::cout << "Linear time: " << Linear_time << "s\n";
    // 计算总时间
    double Server_time = EncryptCipher_time + Add_time + Sbox_time + Linear_time;
    std::cout << "Server total time: " << Server_time << "s\n";
    std::cout << "sbox_timeset: " << sbox_set << endl;
    std::cout << "linear_timeset: " << linear_set << endl;
    if (plainflag)
    {
        if (!verifyDecryption(encryptedPlainStream, PlainStream, secretKey, cmodulus, BlockPlainWords, PlainBlock, nslots, Para_p))
        {
            std::cerr << "Decryption verification failed for PlainStream." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for PlainStream." << std::endl;
    }
    // 计算吞吐量,KiB/min
    double throughput = (Plainbits * 60) / (pow(2, 13) * Server_time);
    std::cout << "Server total time: " << Server_time << "s\n";
    std::cout << "Throughput: " << throughput << "KiB/min\n";
    std::string dirPath = "../tests";
    std::string filePath;
    if (!fs::exists(dirPath))
    {
        filePath = "test_Yux_14.txt";
    }
    else
    {
        filePath = "../tests/test_Yux_14.txt";
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
            << std::left << std::setw(12) << EncryptCipher_time
            << std::left << std::setw(7) << Add_time
            << std::left << std::setw(8) << Sbox_time
            << std::left << std::setw(10) << Linear_time
            << std::left << std::setw(10) << Server_time
            << std::left << std::setw(20) << throughput
            << std::left << std::setw(10) << noise_budget
            << std::endl;
    outfile.close();
    std::cout << "test_Yux_14.txt updated." << std::endl;
    }
    return 0;
}