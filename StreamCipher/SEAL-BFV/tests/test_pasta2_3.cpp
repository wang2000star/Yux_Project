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
#include <math.h>

#include <SEAL-4.1/seal/seal.h>

extern "C"
{
#include "../keccak/KeccakHash.h"
}

#include "../utils/tool.hpp"
#include "../Symmetric/Pasta.hpp"
#include "params.hpp"
using namespace std;
using namespace seal;

typedef __uint128_t uint128_t;

namespace fs = std::filesystem;
/************************************************************************
  long p;          // plaintext primeplain_mod;
  long m;          // m-th cyclotomic polynomial
  long r;          // Lifting [defualt = 1]
  long bits;       // bits in the ciphertext PlainMod chain
  long c;          // columns in the key-switching matrix [default=2]
  long d;          // Degree of the field extension [default=1]
  long k;          // Security parameter [default=80]
  long s;          // Minimum number of slots [default=0]
************************************************************************/
// p^d = 1 mod m,d=1,slots=\phi(m)/d=\phi(m);m=65536=2^16,\phi(m)=2^15=32768
// 更一般的，应该有d|ord_p(m)，slots=\phi(m)/ord_p(m)
//!!!!!!!!!!!!!!!!
constexpr unsigned BlockWords = 256;      // 分组密钥字长度=KeyWords
constexpr unsigned BlockPlainWords = 128; // 明文分组字长度
constexpr double TruncRate = BlockPlainWords / (double)BlockWords;
// ===============模式设置================
constexpr bool deflag = 0;        // true/1表示进行每一步解密验证，false/0表示不进行每一步解密验证
constexpr bool symkeyflag = 0;    // true/1表示对称密钥同态解密验证加密，false/0表示不验证
constexpr bool KeyStreamflag = 0; // true/1表示密钥流同态解密验证，false/0表示不验证
constexpr bool plainflag = 0;     // true/1表示对称密文同态解密验证，false/0表示不验证
// 参数设置，paramMap[Nr-3][idx]
constexpr unsigned Nr = 3; // 轮数
constexpr long idx = 0;
constexpr unsigned Sbox_depth = 1 * Nr; // S盒深度
// 当c=2时，Qbits=1.5*bits,当c=3时，Qbits=1.5*bits - 100
// 硬编码参数值
constexpr tuple<long, long> paramMap[1][8] = {
    {
        // Nr = 4
        // {p, log2(m), bits, c}
        {65537, 15},  // 0 *
        {163841, 15}, // 1
        {65537, 14},  // 2
        {163841, 14}, // 3
        {65537, 15},  // 4
        {163841, 15}, // 5
        {65537, 15},  // 6
        {163841, 15}  // 7
    }};
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
constexpr long nslots = phi_m; // 槽数

constexpr unsigned PlainBlock = nslots - 0; // 明文分组数,应该PlainBlock<=nslots

constexpr unsigned len3 = BlockWords / 3;
// 计算 log2 的 constexpr 函数
constexpr unsigned int log2_constexpr(unsigned long long n, unsigned int p = 0)
{
    return (n <= 1) ? p : log2_constexpr(n / 2, p + 1);
}
constexpr long PlainMod = Para_p; // 明文模数

constexpr unsigned Wordbits = log2_constexpr(PlainMod - 1) + 1; // 字节比特长度=ceil(log2(PlainMod-1))

constexpr unsigned BlockSize = Wordbits * BlockWords;    // 分组比特长度=BlockWords*Wordbits
constexpr unsigned NrBlockWords = BlockWords * (Nr + 1); // Nr轮分组密钥字节长度

constexpr long PlainWords = BlockPlainWords * PlainBlock; // 明文字节长度
constexpr long Plainbits = Wordbits * PlainWords;         // 明文比特长度

constexpr unsigned NonceSize = 32;                           // Nonce比特长度
constexpr long counter_begin = 0;                            // 计数器起始值
constexpr long counter_end = PlainBlock + counter_begin - 1; // 计数器结束值

constexpr long max_prime_size = (1ULL << Wordbits) - 1;
constexpr long R_BITS = 7; // for MDS Matrix
constexpr long x_mask = (1ULL << (Wordbits - R_BITS - 2)) - 1;
constexpr long y_mask = max_prime_size >> 2;
Pasta pasta(PlainMod, BlockPlainWords, max_prime_size);

// Linear transformation
void HE_Fix_M(std::vector<Ciphertext> &encryptedKeyStream, std::vector<std::vector<long>> &fixed_mat, std::vector<std::vector<long>> &fixed_rc1,
              std::vector<std::vector<long>> &fixed_rc2, long r, Evaluator &evaluator, BatchEncoder &batch_encoder)
{
    std::vector<Ciphertext> temp = encryptedKeyStream;
    Ciphertext temp1;
    Plaintext temp2;
    vector<long> temp3(nslots);
    for (int i = 0; i < BlockPlainWords; i++)
    {
        encryptedKeyStream[i] = temp[0];
        for (int k = 0; k < nslots; k++)
        {
            temp3[k] = fixed_mat[i][0];
        }
        batch_encoder.encode(temp3, temp2);
        evaluator.multiply_plain_inplace(encryptedKeyStream[i], temp2);
        for (int j = 1; j < BlockPlainWords; j++)
        {
            temp1 = temp[j];
            for (int k = 0; k < nslots; k++)
            {
                temp3[k] = fixed_mat[i][j];
            }
            batch_encoder.encode(temp3, temp2);
            evaluator.multiply_plain_inplace(temp1, temp2);
            evaluator.add_inplace(encryptedKeyStream[i], temp1);
        }
    }
    for (int i = 0; i < BlockPlainWords; i++)
    {
        encryptedKeyStream[i + BlockPlainWords] = temp[BlockPlainWords];
        for (int k = 0; k < nslots; k++)
        {
            temp3[k] = fixed_mat[i][0];
        }
        batch_encoder.encode(temp3, temp2);
        evaluator.multiply_plain_inplace(encryptedKeyStream[i + BlockPlainWords], temp2);
        for (int j = 1; j < BlockPlainWords; j++)
        {
            temp1 = temp[BlockPlainWords + j];
            for (int k = 0; k < nslots; k++)
            {
                temp3[k] = fixed_mat[i][j];
            }
            batch_encoder.encode(temp3, temp2);
            evaluator.multiply_plain_inplace(temp1, temp2);
            evaluator.add_inplace(encryptedKeyStream[i + BlockPlainWords], temp1);
        }
    }
    for (int i = 0; i < BlockPlainWords; i++)
    {
        for (int k = 0; k < nslots; k++)
        {
            temp3[k] = fixed_rc1[r][i];
        }
        batch_encoder.encode(temp3, temp2);
        evaluator.add_plain_inplace(encryptedKeyStream[i], temp2);
        for (int k = 0; k < nslots; k++)
        {
            temp3[k] = fixed_rc2[r][i];
        }
        batch_encoder.encode(temp3, temp2);
        evaluator.add_plain_inplace(encryptedKeyStream[i + BlockPlainWords], temp2);
    }
    temp = encryptedKeyStream;
    for (int i = 0; i < BlockWords; i++)
    {
        evaluator.add_inplace(encryptedKeyStream[i], temp[i]);
        evaluator.add_inplace(encryptedKeyStream[i], temp[(i + BlockPlainWords) % BlockWords]);
    }
}
void HE_Ran_M(std::vector<Ciphertext> &encryptedKeyStream, std::vector<vector<long>> &xof_mat1, std::vector<vector<long>> &xof_mat2,
              std::vector<vector<long>> &xof_rc, Evaluator &evaluator, BatchEncoder &batch_encoder)
{
    std::vector<Ciphertext> temp = encryptedKeyStream;
    Ciphertext temp1 = temp[0];
    Plaintext encodedxof_mat1;
    Plaintext encodedxof_mat2;
    Plaintext encodedxof_rc;
    for (int i = 0; i < BlockPlainWords; i++)
    {
        encryptedKeyStream[i] = temp[0];
        batch_encoder.encode(xof_mat1[i * BlockPlainWords + 0], encodedxof_mat1);
        evaluator.multiply_plain_inplace(encryptedKeyStream[i], encodedxof_mat1);
        for (int j = 1; j < BlockPlainWords; j++)
        {
            temp1 = temp[j];
            batch_encoder.encode(xof_mat1[i * BlockPlainWords + j], encodedxof_mat1);
            evaluator.multiply_plain_inplace(temp1, encodedxof_mat1);
            evaluator.add_inplace(encryptedKeyStream[i], temp1);
        }
        batch_encoder.encode(xof_rc[i], encodedxof_rc);
        evaluator.add_plain_inplace(encryptedKeyStream[i], encodedxof_rc);
    }
    for (int i = 0; i < BlockPlainWords; i++)
    {
        encryptedKeyStream[i + BlockPlainWords] = temp[BlockPlainWords];
        batch_encoder.encode(xof_mat2[i * BlockPlainWords + 0], encodedxof_mat2);
        evaluator.multiply_plain_inplace(encryptedKeyStream[i + BlockPlainWords], encodedxof_mat2);
        for (int j = 1; j < BlockPlainWords; j++)
        {
            temp1 = temp[BlockPlainWords + j];
            batch_encoder.encode(xof_mat2[i * BlockPlainWords + j], encodedxof_mat2);
            evaluator.multiply_plain_inplace(temp1, encodedxof_mat2);
            evaluator.add_inplace(encryptedKeyStream[i + BlockPlainWords], temp1);
        }
        batch_encoder.encode(xof_rc[i + BlockPlainWords], encodedxof_rc);
        evaluator.add_plain_inplace(encryptedKeyStream[i + BlockPlainWords], encodedxof_rc);
    }
    temp = encryptedKeyStream;
    for (int i = 0; i < BlockWords; i++)
    {
        evaluator.add_inplace(encryptedKeyStream[i], temp[i]);
        evaluator.add_inplace(encryptedKeyStream[i], temp[(i + BlockPlainWords) % BlockWords]);
    }
}
// Compute the Sbox
void HE_Sbox(vector<Ciphertext> &encrypedKeyStream, Evaluator &evaluator, RelinKeys &relin_keys)
{
    vector<Ciphertext> temp = encrypedKeyStream;
    for (int i = 1; i < BlockPlainWords; i++)
    {
        evaluator.exponentiate_inplace(temp[i - 1], 2, relin_keys);
        evaluator.add_inplace(encrypedKeyStream[i], temp[i - 1]);
        evaluator.exponentiate_inplace(temp[i + BlockPlainWords - 1], 2, relin_keys);
        evaluator.add_inplace(encrypedKeyStream[i + BlockPlainWords], temp[i + BlockPlainWords - 1]);
    }
}
// Compute the last Sbox
void HE_Last_Sbox(vector<Ciphertext> &encrypedKeyStream, Evaluator &evaluator, RelinKeys &relin_keys)
{
    for (int i = 0; i < BlockWords; i++)
    {
        evaluator.exponentiate_inplace(encrypedKeyStream[i], 3, relin_keys);
    }
}
//----------------------------------------------------------------
void fixed_init_shake(Keccak_HashInstance &shake128)
{
    uint8_t seed[16] = "PASTA2_";
    seed[7] = '0' + (uint8_t)Nr;
    *((long *)(seed + 8)) = htobe64(PlainMod);

    if (SUCCESS != Keccak_HashInitialize_SHAKE128(&shake128))
        throw std::runtime_error("failed to init shake");
    if (SUCCESS != Keccak_HashUpdate(&shake128, seed, sizeof(seed) * 8))
        throw std::runtime_error("SHAKE128 update failed");
    if (SUCCESS != Keccak_HashFinal(&shake128, NULL))
        throw std::runtime_error("SHAKE128 final failed");
}
std::vector<long> get_random_yi(Keccak_HashInstance &shake128)
{
    std::vector<long> out(BlockPlainWords, 0);
    for (auto i = 0ULL; i < BlockPlainWords; i++)
    {
        bool valid = false;
        // get distinct x_i, which also imply distinct y_i
        while (!valid)
        {
            // random element of size floor(log_2(p))
            uint8_t random_bytes[sizeof(long)];
            if (SUCCESS !=
                Keccak_HashSqueeze(&shake128, random_bytes, sizeof(random_bytes) * 8))
                throw std::runtime_error("SHAKE128 squeeze failed");
            long y_i = be64toh(*((long *)random_bytes)) & y_mask;
            // check distinct x_i
            long x_i = y_i & x_mask;
            valid = true;
            for (auto j = 0ULL; j < i; j++)
            {
                if ((out[j] & x_mask) == x_i)
                {
                    valid = false;
                    break;
                }
            }
            if (valid)
                out[i] = y_i;
        }
    }
    return out;
}
long mod_inverse(long val)
{
    if (val == 0)
        throw(std::runtime_error("0 has no inverse!"));

    int64_t prev_a = 1;
    int64_t a = 0;
    long mod = PlainMod;

    while (mod != 0)
    {
        int64_t q = val / mod;
        int64_t temp = val % mod;
        val = mod;
        mod = temp;

        temp = a;
        a = prev_a - q * a;
        prev_a = temp;
    }

    if (val != 1)
        throw(std::runtime_error("value has no inverse!"));

    long res = prev_a;
    if (prev_a < 0)
        res += PlainMod;

    return res;
}
std::vector<std::vector<long>> get_random_mds_matrix(Keccak_HashInstance &shake128)
{
    std::vector<std::vector<long>> mat(BlockPlainWords, std::vector<long>(BlockPlainWords, 0));
    auto y = get_random_yi(shake128);
    auto x = y;
    for (auto &x_i : x)
        x_i &= x_mask;

    for (auto i = 0ULL; i < BlockPlainWords; i++)
    {
        for (auto j = 0ULL; j < BlockPlainWords; j++)
        {
            // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
            mat[i][j] = mod_inverse((x[i] + y[j]) % PlainMod);
        }
    }
    return mat;
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

    int SecLevel = 128;
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

    // 生成固定矩阵和轮向量
    Keccak_HashInstance shake128;

    fixed_init_shake(shake128);
    vector<vector<long>> fixed_rc1;
    vector<vector<long>> fixed_rc2;
    fixed_rc1.reserve(Nr);
    fixed_rc2.reserve(Nr);

    for (uint16_t i = 0; i < Nr; i++)
    {
        std::vector<long> rc1_(BlockPlainWords);
        std::vector<long> rc2_(BlockPlainWords);

        for (uint16_t j = 0; j < BlockPlainWords; j++)
        {
            rc1_[j] = generate_random_field_element(shake128, false, max_prime_size, PlainMod);
        }
        for (uint16_t j = 0; j < BlockPlainWords; j++)
        {
            rc2_[j] = generate_random_field_element(shake128, false, max_prime_size, PlainMod);
        }
        fixed_rc1.push_back(std::move(rc1_));
        fixed_rc2.push_back(std::move(rc2_));
    }
    std::cout << "固定向量生成成功" << std::endl;
    vector<vector<long>> fixed_mat(BlockPlainWords, vector<long>(BlockPlainWords));
    fixed_mat = get_random_mds_matrix(shake128);
    // for (int i = 0; i < BlockPlainWords; i++)
    // {
    //     for(int j=0;j<BlockPlainWords;j++)
    //     {
    //         std::cout<<fixed_mat[i][j]<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    std::cout << "固定矩阵生成成功" << std::endl;
    // 固定层生成结束
    for (int test = 0; test < 3; test++)
    {
        std::cout << "--------------- Test = " << test << "---------------"<< std::endl;
        // Generating key stream
        Keccak_HashInstance shake128_2;
        vector<vector<long>> mat_FL(BlockPlainWords, vector<long>(BlockPlainWords));
        vector<vector<long>> mat_FR(BlockPlainWords, vector<long>(BlockPlainWords));
        vector<long> rc_L(BlockPlainWords);
        vector<long> rc_R(BlockPlainWords);
        vector<long> state_L(SymKey.begin(), SymKey.begin() + BlockPlainWords);
        vector<long> state_R(SymKey.begin() + BlockPlainWords, SymKey.end());
        std::vector<long> NonceSet(PlainBlock);
        long nonce;
        long block_num;
        std::vector<long> KeyStream(BlockPlainWords * PlainBlock);
        std::cout << "Generating KeyStream..." << std::endl;
        auto start_keyStream = std::chrono::high_resolution_clock::now();
        uint64_t start_cycle1 = rdtsc();
        for (long counter = counter_begin; counter <= counter_end; counter++)
        {
            block_num = counter - counter_begin;
            nonce = generate_secure_random_int(NonceSize);
            NonceSet[block_num] = nonce;
            random_init_shake(nonce, counter, shake128_2);
            mat_FL = pasta.get_random_matrix_pasta(shake128_2);
            mat_FR = pasta.get_random_matrix_pasta(shake128_2);
            for (int i = 0; i < BlockPlainWords; i++)
            {
                rc_L[i] = generate_random_field_element(shake128_2, false, max_prime_size, PlainMod);
                rc_R[i] = generate_random_field_element(shake128_2, false, max_prime_size, PlainMod);
            }
            // 随机初始轮
            pasta.matmul(state_L, mat_FL);
            pasta.matmul(state_R, mat_FR);
            pasta.random_add_rc(state_L, state_R, rc_L, rc_R);
            pasta.mix(state_L, state_R);
            for (int r = 1; r < Nr; r++)
            {
                pasta.sbox_feistel(state_L);
                pasta.sbox_feistel(state_R);
                pasta.matmul(state_L, fixed_mat);
                pasta.matmul(state_R, fixed_mat);
                pasta.add_rc(state_L, state_R, r - 1, fixed_rc1, fixed_rc2);
                pasta.mix(state_L, state_R);
            }
            // 最后一轮
            pasta.sbox_cube(state_L);
            pasta.sbox_cube(state_R);
            pasta.matmul(state_L, fixed_mat);
            pasta.matmul(state_R, fixed_mat);
            pasta.add_rc(state_L, state_R, Nr - 1, fixed_rc1, fixed_rc2);
            pasta.mix(state_L, state_R);
            memcpy(&KeyStream[(block_num)*BlockPlainWords], state_L.data(), BlockPlainWords * sizeof(long));
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
        unsigned plain_size_square = BlockPlainWords * BlockPlainWords;
        std::vector<vector<long>> xof_mat1(plain_size_square, vector<long>(nslots, 0));
        std::vector<vector<long>> xof_mat2 = xof_mat1;
        std::vector<vector<long>> xof_rc(BlockWords, vector<long>(nslots, 0));
        Keccak_HashInstance shake128_3;
        auto start_XOF = std::chrono::high_resolution_clock::now();
        for (long counter = counter_begin; counter <= counter_end; counter++)
        {
            block_num = counter - counter_begin;
            nonce = NonceSet[block_num];
            random_init_shake(nonce, counter, shake128_3);
            mat_FL = pasta.get_random_matrix_pasta(shake128_3);
            mat_FR = pasta.get_random_matrix_pasta(shake128_3);
            for (int i = 0; i < BlockPlainWords; i++)
            {
                rc_L[i] = generate_random_field_element(shake128_2, false, max_prime_size, PlainMod);
                rc_R[i] = generate_random_field_element(shake128_2, false, max_prime_size, PlainMod);
            }
            for (int i = 0; i < BlockPlainWords; i++)
            {
                for (int j = 0; j < BlockPlainWords; j++)
                {
                    xof_mat1[i * BlockPlainWords + j][block_num] = mat_FL[i][j];
                    xof_mat2[i * BlockPlainWords + j][block_num] = mat_FR[i][j];
                }
                xof_rc[i][block_num] = rc_L[i];
                xof_rc[i + BlockPlainWords][block_num] = rc_R[i];
            }
        }
        if (PlainBlock == 1)
        {
            for (int i = 0; i < plain_size_square; i++)
            {
                for (int j = 1; j < nslots; j++)
                {
                    xof_mat1[i][j] = xof_mat1[i][0];
                    xof_mat2[i][j] = xof_mat2[i][0];
                }
            }
            for (int i = 0; i < BlockWords; i++)
            {
                for (int j = 1; j < nslots; j++)
                {
                    xof_rc[i][j] = xof_rc[i][0];
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
        // 定义Sbox_time、Linear_time
        double Sbox_time = 0, Linear_time = 0, white_time = 0;
        auto start_white = std::chrono::high_resolution_clock::now();
        auto end_white = std::chrono::high_resolution_clock::now();
        auto start_sbox = std::chrono::high_resolution_clock::now();
        auto end_sbox = std::chrono::high_resolution_clock::now();
        auto start_linear = std::chrono::high_resolution_clock::now();
        auto end_linear = std::chrono::high_resolution_clock::now();
        vector<double> sbox_set(Nr);
        vector<double> linear_set(Nr + 1);
        vector<Ciphertext> encryptedKeyStream = encryptedSymKey;

        std::cout << "random round start" << std::endl;
        start_white = std::chrono::high_resolution_clock::now();
        HE_Ran_M(encryptedKeyStream, xof_mat1, xof_mat2, xof_rc, evaluator, batch_encoder);
        end_white = std::chrono::high_resolution_clock::now();
        white_time += std::chrono::duration<double>(end_white - start_white).count();
        // 输出 random round time
        std::cout << "random round time: " << white_time << "s\n";
        noise_budget = min_noise_budget(encryptedKeyStream, decryptor);
        std::cout << "noise budget after Whiteround: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            return 0;
        }

        for (long r = 1; r < Nr; r++)
        {
            std::cout << "Round " << r << " start" << std::endl;
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
            start_linear = std::chrono::high_resolution_clock::now();
            // Linear Layer
            HE_Fix_M(encryptedKeyStream, fixed_mat, fixed_rc1, fixed_rc2, r - 1, evaluator, batch_encoder);
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
        }
        // 最后一轮
        std::cout << "the last Round " << Nr << " start" << std::endl;
        start_sbox = std::chrono::high_resolution_clock::now();
        // S Layer
        HE_Last_Sbox(encryptedKeyStream, evaluator, relin_keys);
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
        start_linear = std::chrono::high_resolution_clock::now();
        // Linear Layer
        HE_Fix_M(encryptedKeyStream, fixed_mat, fixed_rc1, fixed_rc2, Nr - 1, evaluator, batch_encoder);
        end_linear = std::chrono::high_resolution_clock::now();
        Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
        linear_set[Nr] = std::chrono::duration<double>(end_linear - start_linear).count();
        // 截断密钥流
        encryptedKeyStream.erase(encryptedKeyStream.begin() + BlockPlainWords, encryptedKeyStream.end());

        noise_budget = min_noise_budget(encryptedKeyStream, decryptor);
        std::cout << "noise budget after linear: " << noise_budget << std::endl;
        if (noise_budget <= 0)
        {
            std::cerr << "noise budget is not enough!!!" << std::endl;
            return 0;
        }
        // 输出 XOF_time,Add_time、Sbox_time、Linear_time
        std::cout << "XOF time: " << XOF_time << "s\n";
        std::cout << "whiteround time: " << white_time << "s\n";
        std::cout << "Sbox time: " << Sbox_time << "s\n";
        std::cout << "Linear time: " << Linear_time << "s\n";
        // 计算总时间
        double Server_offtime = XOF_time + white_time + Sbox_time + Linear_time;
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
        auto start_ClientOnline = std::chrono::high_resolution_clock::now();
        uint64_t start_cycle2 = rdtsc();

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
        uint64_t end_cycle2 = rdtsc();
        auto end_ClientOnline = std::chrono::high_resolution_clock::now();
        uint64_t cycle_count2 = end_cycle2 - start_cycle2;
        double Client_ontime = std::chrono::duration<double>(end_ClientOnline - start_ClientOnline).count();
        std::cout << "Client onine total time:" << std::chrono::duration<double>(end_ClientOnline - start_ClientOnline).count() << "s\n";
        double Client_totaltime = Client_offtime + Client_ontime;
        std::cout << "Client total time: " << Client_totaltime << "s\n";
        uint64_t cycle_count = cycle_count1 + cycle_count2;
        // 输出吞吐量
        double Cli_throughput = (8 * cycle_count) / Plainbits; // 单位：Cycle/Byte
        std::cout << "Client Throughput: " << Cli_throughput << "Cycle/Byte\n";
        // 服务端在线
        // 同态加密
        vector<Ciphertext> encrypedPlainStream = encryptedKeyStream;
        // 对CipherStream进行编码
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
        // for (int i = 0; i < encryptedKeyStream.size(); i++)
        // {
        //     encryptedKeyStream[i].bringToSet(encryptedKeyStream[i].naturalPrimeSet());
        // }
        // if (!verifyDecryption(encrypedPlainStream, PlainStream, secretKey, cmodulus, BlockPlainWords, PlainBlock, nslots, Para_p))
        // {
        //     std::cerr << "Decryption verification failed for encrypedPlainStream." << std::endl;
        //     return 0;
        // }
        // std::cout << "Decryption verification succeeded for encrypedPlainStream." << std::endl;
        // 计算吞吐量,KiB/min
        double Server_totaltime = Server_offtime + server_ontime;
        double Ser_throughput = (Plainbits * 60) / (pow(2, 13) * Server_totaltime);
        std::cout << "Server total time: " << Server_totaltime << "s\n";
        std::cout << "Server Throughput: " << Ser_throughput << "KiB/min\n";
        std::string dirPath = "../tests";
        std::string filePath;
        if (!fs::exists(dirPath))
        {
            filePath = "test_pasta2_3.txt";
        }
        else
        {
            filePath = "../tests/test_pasta2_3.txt";
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
        std::cout << "test_pasta2_3.txt updated." << std::endl;
    }
    return 0;
}