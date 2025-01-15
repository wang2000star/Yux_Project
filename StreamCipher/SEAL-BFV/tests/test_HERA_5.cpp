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

#include "../Symmetric/HERA.hpp"
#include "../utils/tool.hpp"
#include "params.hpp"

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
constexpr unsigned Nr = 5; // 轮数
constexpr long idx = 0;
constexpr unsigned Sbox_depth = 1 * Nr; // S盒深度
// 当c=2时，Qbits=1.5*bits,当c=3时，Qbits=1.5*bits - 100
// 硬编码参数值
constexpr tuple<long, long> paramMap[1][8] = {
    {// Nr = 5
     // {p, log2(m)}
     {65537, 16},  // 0 *
     {163841, 15}, // 1
     {65537, 14},  // 2 *
     {163841, 14}, // 3
     {65537, 15},  // 4
     {163841, 15}, // 5
     {65537, 15},  // 6
     {163841, 15}} // 7
};
// p=k*m+1
//  2^10=1024,2^11=2048,2^12=4096,2^13=8192,2^14=16384,2^15=32768,2^16=65536,
// 当电脑内存有限时，log2Para_m太大，会导致内存不足，出现terminate called recursively错误，从而终止程序.
// 此外，即使正常运行，由于内存一直临界，会导致程序运行速度变慢，时间测量不准确。

constexpr long log2Para_m = get<1>(paramMap[0][idx]) - 0;
constexpr long Para_p = get<0>(paramMap[0][idx]); // plaintext prime
constexpr long Para_m = 1 << log2Para_m;          // cyclotomic polynomial
constexpr long phi_m = Para_m >> 1;               // phi(m)
constexpr long Para_r = 1;                        // Lifting [defualt = 1]
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

HERA hera(PlainMod);

// Linear transformation
void HE_MC(vector<Ciphertext> &eData, Evaluator &evaluator)
{
    vector<Ciphertext> temp = eData;
    array<int, 8> index = {0, 1, 2, 3};
    Ciphertext T2(eData[0]);
    Ciphertext T3(eData[0]);
    for (int i = 0; i < 4; i++)
    {
        int s = 4 * i;
        T2 = temp[index[0] + s];
        evaluator.add_inplace(T2, temp[index[0] + s]);
        T3 = temp[index[1] + s];
        evaluator.add_inplace(T3, temp[index[1] + s]);
        evaluator.add_inplace(T3, temp[index[1] + s]);
        eData[index[0] + s] = T2;
        evaluator.add_inplace(eData[index[0] + s], T3);
        evaluator.add_inplace(eData[index[0] + s], temp[index[2] + s]);
        evaluator.add_inplace(eData[index[0] + s], temp[index[3] + s]);

        T2 = temp[index[1] + s];
        evaluator.add_inplace(T2, temp[index[1] + s]);
        T3 = temp[index[2] + s];
        evaluator.add_inplace(T3, temp[index[2] + s]);
        evaluator.add_inplace(T3, temp[index[2] + s]);
        eData[index[1] + s] = T2;
        evaluator.add_inplace(eData[index[1] + s], T3);
        evaluator.add_inplace(eData[index[1] + s], temp[index[3] + s]);
        evaluator.add_inplace(eData[index[1] + s], temp[index[0] + s]);

        T2 = temp[index[2] + s];
        evaluator.add_inplace(T2, temp[index[2] + s]);
        T3 = temp[index[3] + s];
        evaluator.add_inplace(T3, temp[index[3] + s]);
        evaluator.add_inplace(T3, temp[index[3] + s]);
        eData[index[2] + s] = T2;
        evaluator.add_inplace(eData[index[2] + s], T3);
        evaluator.add_inplace(eData[index[2] + s], temp[index[0] + s]);
        evaluator.add_inplace(eData[index[2] + s], temp[index[1] + s]);

        T2 = temp[index[3] + s];
        evaluator.add_inplace(T2, temp[index[3] + s]);
        T3 = temp[index[0] + s];
        evaluator.add_inplace(T3, temp[index[0] + s]);
        evaluator.add_inplace(T3, temp[index[0] + s]);
        eData[index[3] + s] = T2;
        evaluator.add_inplace(eData[index[3] + s], T3);
        evaluator.add_inplace(eData[index[3] + s], temp[index[1] + s]);
        evaluator.add_inplace(eData[index[3] + s], temp[index[2] + s]);
    }
}

void HE_MR(vector<Ciphertext> &eData, Evaluator &evaluator)
{
    vector<Ciphertext> temp = eData;
    vector<int> index = {0, 4, 8, 12};
    Ciphertext T2(eData[0]);
    Ciphertext T3(eData[0]);
    for (int i = 0; i < 4; i++)
    {
        int s = i;
        T2 = temp[index[0] + s];
        evaluator.add_inplace(T2, temp[index[0] + s]);
        T3 = temp[index[1] + s];
        evaluator.add_inplace(T3, temp[index[1] + s]);
        evaluator.add_inplace(T3, temp[index[1] + s]);
        eData[index[0] + s] = T2;
        evaluator.add_inplace(eData[index[0] + s], T3);
        evaluator.add_inplace(eData[index[0] + s], temp[index[2] + s]);
        evaluator.add_inplace(eData[index[0] + s], temp[index[3] + s]);

        T2 = temp[index[1] + s];
        evaluator.add_inplace(T2, temp[index[1] + s]);
        T3 = temp[index[2] + s];
        evaluator.add_inplace(T3, temp[index[2] + s]);
        evaluator.add_inplace(T3, temp[index[2] + s]);
        eData[index[1] + s] = T2;
        evaluator.add_inplace(eData[index[1] + s], T3);
        evaluator.add_inplace(eData[index[1] + s], temp[index[3] + s]);
        evaluator.add_inplace(eData[index[1] + s], temp[index[0] + s]);

        T2 = temp[index[2] + s];
        evaluator.add_inplace(T2, temp[index[2] + s]);
        T3 = temp[index[3] + s];
        evaluator.add_inplace(T3, temp[index[3] + s]);
        evaluator.add_inplace(T3, temp[index[3] + s]);
        eData[index[2] + s] = T2;
        evaluator.add_inplace(eData[index[2] + s], T3);
        evaluator.add_inplace(eData[index[2] + s], temp[index[0] + s]);
        evaluator.add_inplace(eData[index[2] + s], temp[index[1] + s]);

        T2 = temp[index[3] + s];
        evaluator.add_inplace(T2, temp[index[3] + s]);
        T3 = temp[index[0] + s];
        evaluator.add_inplace(T3, temp[index[0] + s]);
        evaluator.add_inplace(T3, temp[index[0] + s]);
        eData[index[3] + s] = T2;
        evaluator.add_inplace(eData[index[3] + s], T3);
        evaluator.add_inplace(eData[index[3] + s], temp[index[1] + s]);
        evaluator.add_inplace(eData[index[3] + s], temp[index[2] + s]);
    }
}
void HE_M(vector<Ciphertext> &eData, Evaluator &evaluator)
{
    HE_MC(eData, evaluator);
    HE_MR(eData, evaluator);
}
// Compute the constants for Sbox,(x0,x1,x2)——>(x0,x0*x2+x1,-x0*x1+x0*x2+x2)
void HE_Sbox(vector<Ciphertext> &eData, Evaluator &evaluator, RelinKeys &relin_keys)
{
    for (int i = 0; i < eData.size(); i++)
    {
        evaluator.exponentiate_inplace(eData[i], 3, relin_keys);
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

    int SecLevel = 192;
    vector<int> bit_sizes;
    set_params(poly_modulus_degree, SecLevel, bit_sizes);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, bit_sizes));
    // parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(PlainMod);
    SEALContext context(parms, true, sec_level_type::none);
    if (SecLevel == 128)
    {
        context = SEALContext(parms, true, sec_level_type::tc128);
    }
    if (SecLevel == 192)
    {
        context = SEALContext(parms, true, sec_level_type::tc192);
    }
    if (SecLevel == 256)
    {
        context = SEALContext(parms, true, sec_level_type::tc256);
    }
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
    for (int test = 0; test < 10; test++)
    {
        std::cout << "--------------- Test = " << test << "---------------"<< std::endl;
        // Generating key stream
        std::vector<long> NonceSet(PlainBlock);
        vector<std::vector<long>> RoundKeySet(Nr + 1, vector<long>(KeyStreamWords, 0));
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
                {                     // 常规轮
                    hera.MC(state);   // MixColumns
                    hera.MR(state);   // MixRows
                    hera.Sbox(state); // S盒
                    for (unsigned i = 0; i < BlockWords; i++)
                    {
                        state[i] = (state[i] + RoundKey[i]) % PlainMod;
                    }
                }
                else
                {                     // 最后一轮
                    hera.MC(state);   // MixColumns
                    hera.MR(state);   // MixRows
                    hera.Sbox(state); // S盒
                    hera.MC(state);   // MixColumns
                    hera.MR(state);   // MixRows
                    for (unsigned i = 0; i < BlockWords; i++)
                    {
                        state[i] = (state[i] + RoundKey[i]) % PlainMod;
                    }
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
        // 编码IV
        vector<Plaintext> encodedIV(BlockWords);
        for (int i = 0; i < BlockWords; i++)
        {
            vector<long> IVi(nslots, IV[i]);
            batch_encoder.encode(IVi, encodedIV[i]);
        }
        std::cout << "Generating XOF stream..." << std::endl;
        std::vector<vector<long>> Xset(NrBlockWords, vector<long>(nslots, 0));
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

        int noise_budget = min_noise_budget(encryptedSymKey, decryptor);
        std::cout << "noise budget initially: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            return 0;
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
        vector<Ciphertext> encryptedKeyStream = encryptedSymKey;
        std::cout << "whiteround start" << std::endl;
        Plaintext encodedXseti;
        start_add = std::chrono::high_resolution_clock::now();
        for (long i = 0; i < BlockWords; i++)
        { // encrypt the encoded key
            batch_encoder.encode(Xset[i], encodedXseti);
            evaluator.multiply_plain_inplace(encryptedKeyStream[i], encodedXseti);
            evaluator.add_plain_inplace(encryptedKeyStream[i], encodedIV[i]);
        }
        end_add = std::chrono::high_resolution_clock::now();
        Add_time += std::chrono::duration<double>(end_add - start_add).count();
        // 输出 Add_time
        std::cout << "whiteround time: " << Add_time << "s\n";
        noise_budget = min_noise_budget(encryptedKeyStream, decryptor);
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
                KeyStream2[i] = (IV[i % BlockWords] + RoundKeySet[0][i]) % PlainMod;
            }
            // 使用 verifyDecryption 函数解密并验证 KeyStream2
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor, BlockWords, PlainBlock, nslots, PlainMod))
            {
                std::cerr << "Decryption verification failed for KeyStream2." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for whiteround." << std::endl;
        }
        Ciphertext encryptedRoundKeySeti;
        for (long r = 1; r < Nr; r++)
        {
            std::cout << "Round " << r << " start" << std::endl;
            start_linear = std::chrono::high_resolution_clock::now();
            // Linear Layer
            HE_M(encryptedKeyStream, evaluator);
            end_linear = std::chrono::high_resolution_clock::now();
            Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
            linear_set[r - 1] = std::chrono::duration<double>(end_linear - start_linear).count();
            noise_budget = min_noise_budget(encryptedKeyStream, decryptor);
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
                if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor, BlockWords, PlainBlock, nslots, PlainMod))
                {
                    std::cerr << "Decryption verification failed for KeyStream2 Linear Layer." << std::endl;
                    return 0;
                }
                std::cout << "Decryption verification succeeded for KeyStream2 Linear Layer." << std::endl;
            }
            start_sbox = std::chrono::high_resolution_clock::now();
            // S Layer
            HE_Sbox(encryptedKeyStream, evaluator, relin_keys);
            end_sbox = std::chrono::high_resolution_clock::now();
            Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
            sbox_set[r - 1] = std::chrono::duration<double>(end_sbox - start_sbox).count();
            noise_budget = min_noise_budget(encryptedKeyStream, decryptor);
            std::cout << "noise budget after sbox: " << noise_budget << std::endl;
            if (noise_budget <= 0)
            {
                std::cerr << "noise budget is not enough!!!" << std::endl;
                return 0;
            }
            if (deflag)
            {
                hera.Sbox(KeyStream2);
                if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor, BlockWords, PlainBlock, nslots, PlainMod))
                {
                    std::cerr << "Decryption verification failed for KeyStream2 Sbox." << std::endl;
                    return 0;
                }
                std::cout << "Decryption verification succeeded for KeyStream2 Sbox." << std::endl;
            }
            // Round Key Addition
            start_add = std::chrono::high_resolution_clock::now();
            for (long j = 0; j < BlockWords; j++)
            {
                batch_encoder.encode(Xset[r * BlockWords + j], encodedXseti);
                encryptedRoundKeySeti = encryptedSymKey[j];
                evaluator.multiply_plain_inplace(encryptedRoundKeySeti, encodedXseti);
                evaluator.add_inplace(encryptedKeyStream[j], encryptedRoundKeySeti);
            }
            end_add = std::chrono::high_resolution_clock::now();
            Add_time += std::chrono::duration<double>(end_add - start_add).count();
            noise_budget = min_noise_budget(encryptedKeyStream, decryptor);
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
                    KeyStream2[i] = (KeyStream2[i] + RoundKeySet[r][i]) % PlainMod;
                }
                if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor, BlockWords, PlainBlock, nslots, PlainMod))
                {
                    std::cerr << "Decryption verification failed for KeyStream2 Round Key Addition." << std::endl;
                    return 0;
                }
                std::cout << "Decryption verification succeeded for KeyStream2 Round Key Addition." << std::endl;
            }
        }
        // 最后一轮
        std::cout << "the last Round " << Nr << " start" << std::endl;
        start_linear = std::chrono::high_resolution_clock::now();
        // Linear Layer
        HE_M(encryptedKeyStream, evaluator);
        end_linear = std::chrono::high_resolution_clock::now();
        Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
        linear_set[Nr - 1] = std::chrono::duration<double>(end_linear - start_linear).count();
        noise_budget = min_noise_budget(encryptedKeyStream, decryptor);
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
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor, BlockWords, PlainBlock, nslots, PlainMod))
            {
                std::cerr << "Decryption verification failed for KeyStream2 Linear Layer." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream2 Linear Layer." << std::endl;
        }
        start_sbox = std::chrono::high_resolution_clock::now();
        // S Layer
        HE_Sbox(encryptedKeyStream, evaluator, relin_keys);
        // HE_Last_Sbox(encryptedKeyStream);
        end_sbox = std::chrono::high_resolution_clock::now();
        Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
        sbox_set[Nr - 1] = std::chrono::duration<double>(end_sbox - start_sbox).count();
        noise_budget = min_noise_budget(encryptedKeyStream, decryptor);
        std::cout << "noise budget after sbox: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            return 0;
        }
        if (deflag)
        {
            // yusP.Sbox_last(KeyStream2);
            hera.Sbox(KeyStream2);
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor, BlockWords, PlainBlock, nslots, PlainMod))
            {
                std::cerr << "Decryption verification failed for KeyStream2 Sbox." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream2 Sbox." << std::endl;
        }
        start_linear = std::chrono::high_resolution_clock::now();
        // Linear Layer
        HE_M(encryptedKeyStream, evaluator);
        end_linear = std::chrono::high_resolution_clock::now();
        Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
        linear_set[Nr] = std::chrono::duration<double>(end_linear - start_linear).count();
        noise_budget = min_noise_budget(encryptedKeyStream, decryptor);
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
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor, BlockWords, PlainBlock, nslots, PlainMod))
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
            batch_encoder.encode(Xset[Nr * BlockWords + j], encodedXseti);
            encryptedRoundKeySeti = encryptedSymKey[j];
            evaluator.multiply_plain_inplace(encryptedRoundKeySeti, encodedXseti);
            evaluator.add_inplace(encryptedKeyStream[j], encryptedRoundKeySeti);
        }
        end_add = std::chrono::high_resolution_clock::now();
        Add_time += std::chrono::duration<double>(end_add - start_add).count();
        noise_budget = min_noise_budget(encryptedKeyStream, decryptor);
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
                KeyStream2[i] = (KeyStream2[i] + RoundKeySet[Nr][i]) % PlainMod;
            }
            if (!verifyDecryption(encryptedKeyStream, KeyStream2, batch_encoder, decryptor, BlockWords, PlainBlock, nslots, PlainMod))
            {
                std::cerr << "Decryption verification failed for KeyStream2 Round Key Addition." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream2 Round Key Addition." << std::endl;
        }

        // 输出 XOF_time,Add_time、Sbox_time、Linear_time
        std::cout << "XOF time: " << XOF_time << "s\n";
        std::cout << "Add time: " << Add_time << "s\n";
        std::cout << "Sbox time: " << Sbox_time << "s\n";
        std::cout << "Linear time: " << Linear_time << "s\n";
        // 计算总时间
        double Server_offtime = XOF_time + Add_time + Sbox_time + Linear_time;
        std::cout << "Server offline total time: " << Server_offtime << "s\n";
        std::cout << "sbox_timeset: ";
        for (const auto &time : sbox_set)
        {
            std::cout << time << " ";
        }
        std::cout << endl;
        std::cout << "linear_timeset: ";
        for (const auto &time : linear_set)
        {
            std::cout << time << " ";
        }
        std::cout << endl;
        if (KeyStreamflag)
        {
            if (!verifyDecryption(encryptedKeyStream, KeyStream, batch_encoder, decryptor, BlockPlainWords, PlainBlock, nslots, PlainMod))
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
        vector<vector<long>> CipherStream(BlockPlainWords, vector<long>(nslots, 0));
        long jBP;
        auto start_ClientOnline = std::chrono::high_resolution_clock::now();
        uint64_t start_cycle2 = rdtsc();
        for (int j = 0; j < PlainBlock; j++)
        {
            jBP = j * BlockPlainWords;
            for (int i = 0; i < BlockPlainWords; i++)
            {
                CipherStream[i][j] = (PlainStream[jBP + i] + KeyStream[jBP + i]) % PlainMod;
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
        vector<Ciphertext> encrypedPlainStream = encryptedKeyStream;
        Plaintext encodedCipherStreami;
        auto start_ServerOnline = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < BlockPlainWords; i++)
        {
            batch_encoder.encode(CipherStream[i], encodedCipherStreami);
            evaluator.negate_inplace(encrypedPlainStream[i]);
            evaluator.add_plain_inplace(encrypedPlainStream[i], encodedCipherStreami);
        }
        auto end_ServerOnline = std::chrono::high_resolution_clock::now();
        double server_ontime = std::chrono::duration<double>(end_ServerOnline - start_ServerOnline).count();
        std::cout << "Server onine total time:" << server_ontime << "s\n";
        noise_budget = min_noise_budget(encryptedKeyStream, decryptor);
        std::cout << "noise budget after online: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            return 0;
        }
        // 同态解密验证
        if (plainflag)
        {
            if (!verifyDecryption(encrypedPlainStream, PlainStream, batch_encoder, decryptor, BlockPlainWords, PlainBlock, nslots, PlainMod))
            {
                std::cerr << "Decryption verification failed for encrypedPlainStream." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for encrypedPlainStream." << std::endl;
        }
        // 计算吞吐量,KiB/min
        double Server_totaltime = Server_offtime + server_ontime;
        double Ser_throughput = (Plainbits * 60) / (pow(2, 13) * Server_totaltime);
        std::cout << "Server total time: " << Server_totaltime << "s\n";
        std::cout << "Server Throughput: " << Ser_throughput << "KiB/min\n";
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
        outfile.close();
        std::cout << "test_HERA_5.txt updated." << std::endl;
    }
    return 0;
}