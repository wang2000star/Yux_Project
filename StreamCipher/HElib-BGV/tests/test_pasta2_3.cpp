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

#include <NTL/ZZX.h>
// #include <NTL/GF2X.h>
#include <helib/helib.h>
#include <helib/ArgMap.h>
#include <helib/DoubleCRT.h>

#include "helib/CModulus.h"
#include "helib/powerful.h"
#include <immintrin.h> // 包含 SIMD 指令集支持

#include "random_bit.hpp"

extern "C"
{
#include "../hybrid-HE-framework/util/keccak/KeccakHash.h"
}

#include "../utils/tool.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

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
constexpr unsigned BlockWords = 256; // 分组密钥字长度
constexpr unsigned plain_size = 128;
constexpr unsigned TruncWords = 128; // 截断字节长度
constexpr double TruncRate = TruncWords / (double)BlockWords;
// ===============模式设置================
constexpr bool Rkflag = 1;     // true/1表示乘法，false/0表示加法，指的是随机向量和密钥间的操作
constexpr bool deflag = 0;     // true/1表示进行每一步解密验证，false/0表示不进行每一步解密验证
constexpr bool ompflag = 0;    // true/1表示使用OpenMP并行编码，false/0表示不使用OpenMP并行编码
constexpr bool symkeyflag = 0; // true/1表示对称密钥同态解密验证加密，false/0表示不验证
constexpr bool plainflag = 0;  // true/1表示对称密文同态解密验证，false/0表示不验证
// 参数设置，paramMap[Nr-3][idx]
constexpr unsigned Nr = 3; // 轮数
constexpr long idx = 0;
constexpr unsigned Sbox_depth = 1 * Nr; // S盒深度
// 当c=2时，Qbits=1.5*bits,当c=3时，Qbits=1.5*bits - 100
// 硬编码参数值
constexpr tuple<long, long, long, long> paramMap[5][8] = {
    {// Nr = 3
     // {p, log2(m), bits, c}
     {65537, 16, 300, 2},   // 0 *
     {163841, 15, 240, 2},  // 1
     {65537, 14, 220, 2},   // 2 *
     {163841, 14, 230, 2},  // 3
     {65537, 15, 350, 2},   // 4
     {163841, 15, 350, 2},  // 5
     {65537, 15, 400, 2},   // 6
     {163841, 15, 400, 2}}, // 7
    {
        // Nr = 4
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
        {65537, 16, 280, 2}, // 0 *
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

constexpr long log2Para_m = get<1>(paramMap[Nr - 3][idx]) - 0;
constexpr long Para_p = get<0>(paramMap[Nr - 3][idx]);    // plaintext prime
constexpr long Para_m = 1 << log2Para_m;                  // cyclotomic polynomial
constexpr long phi_m = Para_m >> 1;                       // phi(m)
constexpr long Para_bits = get<2>(paramMap[Nr - 3][idx]); // bits in the ciphertext PlainMod chain
constexpr long Para_c = get<3>(paramMap[Nr - 3][idx]);    // columns in the key-switching matrix
constexpr long Para_r = 1;                                // Lifting [defualt = 1]
//!!!!!!!!!!!!!!!
constexpr long nslots = phi_m;             // 槽数
constexpr unsigned PlainBlock = phi_m - 0; // 明文分组数,应该PlainBlock<=phi_m
constexpr unsigned len3 = BlockWords / 3;
// 计算 log2 的 constexpr 函数
constexpr unsigned int log2_constexpr(unsigned long long n, unsigned int p = 0)
{
    return (n <= 1) ? p : log2_constexpr(n / 2, p + 1);
}
constexpr long PlainMod = Para_p;                               // 明文模数
constexpr unsigned Wordbits = log2_constexpr(PlainMod - 1) + 1; // 字节比特长度=ceil(log2(PlainMod-1))
constexpr unsigned randbits = Wordbits - 1;
constexpr unsigned BlockSize = Wordbits * BlockWords;     // 分组比特长度=BlockWords*Wordbits
constexpr unsigned NrBlockWords = BlockWords * (Nr + 1);  // Nr轮分组密钥字节长度
constexpr long PlainWords = BlockWords * PlainBlock;      // 明文字节长度
constexpr long TruncPlainWords = TruncWords * PlainBlock; // 截断后的明文字节长度
constexpr long Plainbits = Wordbits * PlainWords;         // 明文比特长度

constexpr unsigned NonceSize = 32;                           // Nonce比特长度
constexpr long counter_begin = 0;                            // 计数器起始值
constexpr long counter_end = PlainBlock + counter_begin - 1; // 计数器结束值

constexpr long max_prime_size = (1ULL << 17) - 1;
constexpr long R_BITS = 7; // for MDS Matrix
constexpr long x_mask = (1ULL << (17 - R_BITS - 2)) - 1;
constexpr long y_mask = max_prime_size >> 2;
// YusP yusP(PlainMod); // 构建明文对称加密实例

// Linear transformation
void HE_M(std::vector<Ctxt> encrypedKeyStream, std::vector<std::vector<long>> &fixed_mat, std::vector<std::vector<long>> fixed_rc1,
          std::vector<std::vector<long>> fixed_rc2, long r)
{
    std::vector<Ctxt> temp = encrypedKeyStream;
    Ctxt temp1 = temp[0];
    for (int i = 0; i < plain_size; i++)
    {
        encrypedKeyStream[i] = temp[0];
        encrypedKeyStream[i].multByConstant(static_cast<ZZX>(fixed_mat[i][0]));
        for (int j = 1; j < plain_size; j++)
        {
            temp1 = temp[j];
            temp1.multByConstant(static_cast<ZZX>(fixed_mat[i][j]));
            encrypedKeyStream[i] += temp1;
        }
    }
    for (int i = 0; i < plain_size; i++)
    {
        encrypedKeyStream[i + plain_size] = temp[plain_size];
        encrypedKeyStream[i + plain_size].multByConstant(static_cast<ZZX>(fixed_mat[i][0]));
        for (int j = 1; j < plain_size; j++)
        {
            temp1 = temp[plain_size + j];
            temp1.multByConstant(static_cast<ZZX>(fixed_mat[i][j]));
            encrypedKeyStream[i + plain_size] += temp1;
        }
    }
    for (int i = 0; i < plain_size; i++)
    {
        encrypedKeyStream[i].addConstant(to_ZZX(fixed_rc1[r][i]));
        encrypedKeyStream[i + plain_size].addConstant(to_ZZX(fixed_rc2[r][i]));
    }
    temp = encrypedKeyStream;
    for (int i = 0; i < BlockWords; i++)
    {
        encrypedKeyStream[i] += temp[i];
        encrypedKeyStream[i] += temp[(i + plain_size) % BlockWords];
    }
}
// Compute the constants for Sbox,(x0,x1,x2)——>(x0,x0*x2+x1,-x0*x1+x0*x2+x2)
void HE_Sbox(vector<Ctxt> &encrypedKeyStream)
{
    vector<Ctxt> temp = encrypedKeyStream;
    for (int i = 1; i < plain_size; i++)
    {
        temp[i - 1].square();
        encrypedKeyStream[i] += temp[i - 1];
        temp[i + plain_size - 1].square();
        encrypedKeyStream[i + plain_size] += temp[i + plain_size - 1];
    }
}
// Compute the constants for the last Sbox,(x0,x1,x2)——>(x0^3,x0*x2+x1,-x0*x1+x0*x2+x2)
void HE_Last_Sbox(vector<Ctxt> &encrypedKeyStream)
{
    for (int i = 0; i < BlockWords; i++)
    {
        encrypedKeyStream[i].cube();
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
void random_init_shake(long nonce, long block_counter, Keccak_HashInstance &shake128_)
{
    uint8_t seed[16];

    *((long *)seed) = htobe64(nonce);
    *((long *)(seed + 8)) = htobe64(block_counter);

    if (SUCCESS != Keccak_HashInitialize_SHAKE128(&shake128_))
        throw std::runtime_error("failed to init shake");
    if (SUCCESS != Keccak_HashUpdate(&shake128_, seed, sizeof(seed) * 8))
        throw std::runtime_error("SHAKE128 update failed");
    if (SUCCESS != Keccak_HashFinal(&shake128_, NULL))
        throw std::runtime_error("SHAKE128 final failed");
}
long generate_random_field_element(Keccak_HashInstance &shake128, bool allow_zero)
{
    uint8_t random_bytes[sizeof(long)];
    while (1)
    {
        if (SUCCESS !=
            Keccak_HashSqueeze(&shake128, random_bytes, sizeof(random_bytes) * 8))
            throw std::runtime_error("SHAKE128 squeeze failed");
        long ele = be64toh(*((long *)random_bytes)) & max_prime_size;
        if (!allow_zero && ele == 0)
            continue;
        if (ele < PlainMod)
            return ele;
    }
}
std::vector<long> get_random_yi(Keccak_HashInstance &shake128)
{
    std::vector<long> out(plain_size, 0);
    for (auto i = 0ULL; i < plain_size; i++)
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
    std::vector<std::vector<long>> mat(plain_size, std::vector<long>(plain_size, 0));
    auto y = get_random_yi(shake128);
    auto x = y;
    for (auto &x_i : x)
        x_i &= x_mask;

    for (auto i = 0ULL; i < plain_size; i++)
    {
        for (auto j = 0ULL; j < plain_size; j++)
        {
            // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
            mat[i][j] = mod_inverse((x[i] + y[j]) % PlainMod);
        }
    }
    return mat;
}
std::vector<std::vector<long>> get_random_matrix_pasta(Keccak_HashInstance &shake128)
{
    std::vector<std::vector<long>> mat(plain_size, std::vector<long>(plain_size));
    for (auto &m : mat[0])
    {
        m = generate_random_field_element(shake128, false);
    }
    const auto &first_row = mat[0];
    auto &prev_row = mat[0];
    for (auto i = 1ULL; i < plain_size; i++)
    {
        for (auto j = 0ULL; j < plain_size; j++)
        {
            long tmp =
                ((uint128_t)(first_row[j]) * prev_row[plain_size - 1]) % PlainMod;
            if (j)
            {
                tmp = (tmp + prev_row[j - 1]) % PlainMod;
            }
            mat[i][j] = tmp;
        }
        prev_row = mat[i];
    }
    return mat;
}
void Pasta2_sbox_cube(vector<long> &state)
{
    for (uint16_t el = 0; el < plain_size; el++)
    {
        long square = ((uint128_t)(state[el]) * state[el]) % PlainMod;
        state[el] = ((uint128_t)(square)*state[el]) % PlainMod;
    }
}
void Pasta2_sbox_feistel(vector<long> &state)
{
    vector<long> new_state(plain_size);
    new_state[0] = state[0];
    for (uint16_t el = 1; el < plain_size; el++)
    {
        long square = ((uint128_t)(state[el - 1]) * state[el - 1]) % PlainMod;
        // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
        new_state[el] = (square + state[el]) % PlainMod;
    }
    state = new_state;
}
void Pasta2_matmul(vector<long> &state, const std::vector<std::vector<long>> &matrix)
{
    vector<long> new_state(plain_size, 0);
    for (uint16_t i = 0; i < plain_size; i++)
    {
        for (uint16_t j = 0; j < plain_size; j++)
        {
            long mult = ((uint128_t)(matrix[i][j]) * state[j]) % PlainMod;
            // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
            new_state[i] = (new_state[i] + mult) % PlainMod;
        }
    }
    state = new_state;
}
void Pasta2_random_add_rc(vector<long> &state1, vector<long> &state2, vector<long> &rc1, vector<long> &rc2)
{
    for (uint16_t el = 0; el < plain_size; el++)
    {
        // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
        state1[el] = (state1[el] + rc1[el]) % PlainMod;
        state2[el] = (state2[el] + rc2[el]) % PlainMod;
    }
}
void Pasta2_add_rc(vector<long> &state1, vector<long> &state2, size_t r, vector<vector<long>> &rc1, vector<vector<long>> &rc2)
{
    for (uint16_t el = 0; el < plain_size; el++)
    {
        // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
        state1[el] = (state1[el] + rc1[r][el]) % PlainMod;
        state2[el] = (state2[el] + rc2[r][el]) % PlainMod;
    }
}
void Pasta2_mix(vector<long> &state1, vector<long> &state2)
{
    // just adding
    for (uint16_t i = 0; i < plain_size; i++)
    {
        // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
        long sum = (state1[i] + state2[i]);
        state1[i] = (state1[i] + sum) % PlainMod;
        state2[i] = (state2[i] + sum) % PlainMod;
    }
}

int main()
{
    std::cout << "Nr: " << Nr << std::endl;
    std::cout << "BlockWord: " << BlockWords << std::endl;
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
    encryptSymKey(encryptedSymKey, SymKey, publicKey, cmodulus, nslots);
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
    std::vector<long> RoundKeySet(PlainWords * (Nr + 1));
    std::vector<long> KeyStream(plain_size * PlainBlock);
    RandomBit<BlockWords * randbits> randomBit(Nr);
    long X;
    std::vector<long> RoundKey(BlockWords);
    bool bit_array[randbits];
    long nonce;
    long block_num;
    long ir;
    std::vector<long> state(BlockWords); // 初始化 state
    std::cout << "Generating KeyStream..." << std::endl;
    long start_cycle = rdtsc();

    // 生成固定矩阵和轮向量
    Keccak_HashInstance shake128;

    fixed_init_shake(shake128);
    vector<vector<long>> fixed_rc1;
    vector<vector<long>> fixed_rc2;
    fixed_rc1.reserve(Nr);
    fixed_rc2.reserve(Nr);

    for (uint16_t i = 0; i < Nr; i++)
    {
        std::vector<long> rc1_(plain_size);
        std::vector<long> rc2_(plain_size);

        for (uint16_t j = 0; j < plain_size; j++)
        {
            rc1_[j] = generate_random_field_element(shake128, false);
        }
        for (uint16_t j = 0; j < plain_size; j++)
        {
            rc2_[j] = generate_random_field_element(shake128, false);
        }
        fixed_rc1.push_back(std::move(rc1_));
        fixed_rc2.push_back(std::move(rc2_));
    } 
    std::cout<<"固定向量生成成功"<<std::endl;
    vector<vector<long>> fixed_mat(plain_size, vector<long>(plain_size));
    fixed_mat = get_random_mds_matrix(shake128);
    // for (int i = 0; i < plain_size; i++)
    // {
    //     for(int j=0;j<plain_size;j++)
    //     {
    //         std::cout<<fixed_mat[i][j]<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    std::cout<<"固定矩阵生成成功"<<std::endl;
    // 固定层生成结束

    Keccak_HashInstance shake128_2;
    vector<vector<long>> mat_FL(plain_size, vector<long>(plain_size));
    vector<vector<long>> mat_FR(plain_size, vector<long>(plain_size));
    vector<long> rc_L(plain_size);
    vector<long> rc_R(plain_size);
    vector<long> state_L(SymKey.begin(), SymKey.begin() + plain_size);
    vector<long> state_R(SymKey.begin() + plain_size, SymKey.end());
    for (long counter = counter_begin; counter <= counter_end; counter++)
    {
        nonce = generate_secure_random_int(NonceSize);
        random_init_shake(nonce, counter, shake128_2);
        mat_FL = get_random_matrix_pasta(shake128_2);
        mat_FR = get_random_matrix_pasta(shake128_2);
        for (int i = 0; i < plain_size; i++)
        {
            rc_L[i] = generate_random_field_element(shake128_2, false);
            rc_R[i] = generate_random_field_element(shake128_2, false);
        }
        // 随机初始轮
        Pasta2_matmul(state_L, mat_FL);
        Pasta2_matmul(state_R, mat_FR);
        Pasta2_random_add_rc(state_L, state_R, rc_L, rc_R);
        Pasta2_mix(state_L, state_R);
        for (int r = 1; r < Nr; r++)
        {
            Pasta2_sbox_feistel(state_L);
            Pasta2_sbox_feistel(state_R);
            Pasta2_matmul(state_L, fixed_mat);
            Pasta2_matmul(state_R, fixed_mat);
            Pasta2_add_rc(state_L, state_R, r - 1, fixed_rc1, fixed_rc2);
            Pasta2_mix(state_L, state_R);
        }
        // 最后一轮
        Pasta2_sbox_cube(state_L);
        Pasta2_sbox_cube(state_R);
        Pasta2_matmul(state_L, fixed_mat);
        Pasta2_matmul(state_R, fixed_mat);
        Pasta2_add_rc(state_L, state_R, Nr - 1, fixed_rc1, fixed_rc2);
        Pasta2_mix(state_L, state_R);
        memcpy(&KeyStream[(block_num)*plain_size], state_L.data(), plain_size * sizeof(long));
    }
    long end_cycle = rdtsc();
    auto end_keyStream = std::chrono::high_resolution_clock::now();
    double Client_offtime = std::chrono::duration_cast<std::chrono::duration<double>>(end_keyStream - start_keyStream).count();
    std::cout << "Encryption offline total time: " << Client_offtime << "s\n";
    long cycle_count = end_cycle - start_cycle;
    std::cout << "Encryption offline total cycles: " << cycle_count << std::endl;
    // 将KeyStream同态密文写入文件
    // if (!writeEncryptedSymKey(encryptedSymKey, "Client_encryptedSymKey.bin"))
    // {
    //     return false;
    // }
    // std::cout << "Client_encryptedSymKey.bin has been written to file." << std::endl;

    //=============服务端offline阶段================
    std::cout << "Generating XOF stream..." << std::endl;
    unsigned plain_size_square = plain_size * plain_size;
    std::vector<vec_long> xof_mat1(plain_size_square);
    for (int i = 0; i < plain_size_square; i++)
    {
        xof_mat1[i].SetLength(nslots);
    }
    std::vector<vec_long> xof_mat2 = xof_mat1;
    std::vector<vec_long> xof_rc(BlockWords);
    for (int i = 0; i < BlockWords; i++)
    {
        xof_rc[i].SetLength(nslots);
    }
    Keccak_HashInstance shake128_3;
    auto start_XOF = std::chrono::high_resolution_clock::now();
    for (long counter = counter_begin; counter <= counter_end; counter++)
    {
        block_num = counter - counter_begin;
        nonce = NonceSet[block_num];
        random_init_shake(nonce, counter, shake128_3);
        mat_FL = get_random_matrix_pasta(shake128_3);
        mat_FR = get_random_matrix_pasta(shake128_3);
        for (int i = 0; i < plain_size; i++)
        {
            rc_L[i] = generate_random_field_element(shake128_2, false);
            rc_R[i] = generate_random_field_element(shake128_2, false);
        }
        for (int i = 0; i < plain_size; i++)
        {
            for (int j = 0; j < plain_size; j++)
            {
                xof_mat1[i * plain_size + j][block_num] = mat_FL[i][j];
                xof_mat2[i * plain_size + j][block_num] = mat_FR[i][j];
            }
            xof_rc[i][block_num] = rc_L[i];
            xof_rc[i + plain_size][block_num] = rc_R[i];
        }
    }
    auto end_XOF = std::chrono::high_resolution_clock::now();
    double XOF_time = std::chrono::duration<double>(end_XOF - start_XOF).count();
    std::cout << "XOF stream Generation time: " << XOF_time << "s\n";

    zz_pX encodedtemp;
    std::vector<zzX> encodedxof_mat1(plain_size_square);
    std::vector<zzX> encodedxof_mat2(plain_size_square);
    std::vector<ZZX> encodedxof_rc(BlockWords);
    auto start_encode = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < plain_size_square; i++)
    {
        cmodulus.iFFT(encodedtemp, xof_mat1[i]);
        convert(encodedxof_mat1[i], encodedtemp);
        cmodulus.iFFT(encodedtemp, xof_mat2[i]);
        convert(encodedxof_mat1[i], encodedtemp);
    }
    for (int i = 0; i < BlockWords; i++)
    {
        cmodulus.iFFT(encodedtemp, xof_rc[i]);
        convert(encodedxof_rc[i], encodedtemp);
    }
    auto end_encode = std::chrono::high_resolution_clock::now();
    double Encode_time = std::chrono::duration<double>(end_encode - start_encode).count();
    std::cout << "encode time: " << Encode_time << "s\n";

    Ctxt tmpCtxt(*publicKey);
    int noise_budget = min_noise_budget(encryptedSymKey);
    std::cout << "noise budget initially: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }

    // 生成 encryptedKeyStream
    // 定义Sbox_time、Linear_time
    double Sbox_time = 0, Linear_time = 0;
    auto start_sbox = std::chrono::high_resolution_clock::now();
    auto end_sbox = std::chrono::high_resolution_clock::now();
    auto start_linear = std::chrono::high_resolution_clock::now();
    auto end_linear = std::chrono::high_resolution_clock::now();
    vector<double> sbox_set(Nr);
    vector<double> linear_set(Nr + 1);
    vector<Ctxt> encryptedKeyStream = encryptedSymKey;

    std::cout << "random round start" << std::endl;
    Ctxt temp1 = tmpCtxt;
    vector<Ctxt> temp = encryptedKeyStream;
    auto start_random = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < plain_size; i++)
    {
        encryptedKeyStream[i] = temp[0];
        encryptedKeyStream[i].multByConstant(encodedxof_mat1[i * plain_size + 0]);
        for (int j = 1; j < plain_size; j++)
        {
            temp1 = temp[j];
            temp1.multByConstant(encodedxof_mat1[i * plain_size + j]);
            encryptedKeyStream[i] += temp1;
        }
        encryptedKeyStream[i].addConstant(encodedxof_rc[i]);
    }
    for (int i = 0; i < plain_size; i++)
    {
        encryptedKeyStream[i + plain_size] = temp[plain_size];
        encryptedKeyStream[i + plain_size].multByConstant(
            encodedxof_mat1[i * plain_size + 0]);
        for (int j = 1; j < plain_size; j++)
        {
            temp1 = temp[plain_size + j];
            temp1.multByConstant(encodedxof_mat1[i * plain_size + j]);
            encryptedKeyStream[i] += temp1;
        }
        encryptedKeyStream[i + plain_size].addConstant(
            encodedxof_rc[i + plain_size]);
    }
    temp = encryptedKeyStream;
    for (int i = 0; i < BlockWords; i++)
    {
        encryptedKeyStream[i] += temp[i];
        encryptedKeyStream[i] += temp[(i + plain_size) % BlockWords];
    }
    auto end_random = std::chrono::high_resolution_clock::now();
    double random_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_random - start_random).count();
    // 输出 random round time
    std::cout << "random round time: " << random_time << "s\n";
    noise_budget = min_noise_budget(encryptedKeyStream);
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
        start_linear = std::chrono::high_resolution_clock::now();
        // Linear Layer
        HE_M(encryptedKeyStream, fixed_mat, fixed_rc1, fixed_rc2, r - 1);
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
    }
// 最后一轮
#if (1)
    std::cout << "the last Round " << Nr << " start" << std::endl;
    start_sbox = std::chrono::high_resolution_clock::now();
    // S Layer
    HE_Last_Sbox(encryptedKeyStream);
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
    start_linear = std::chrono::high_resolution_clock::now();
    // Linear Layer
    HE_M(encryptedKeyStream, fixed_mat, fixed_rc1, fixed_rc2, Nr - 1);
    end_linear = std::chrono::high_resolution_clock::now();
    Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
    linear_set[Nr] = std::chrono::duration<double>(end_linear - start_linear).count();
    //舍去后面的密文
    vector<Ctxt> TruncedencryptedKeyStream(TruncWords, tmpCtxt);
    for(int i=0;i<TruncWords;i++){
        TruncedencryptedKeyStream[i]=encryptedKeyStream[i];
    }
    noise_budget = min_noise_budget(encryptedKeyStream);
    std::cout << "noise budget after linear: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
#endif

    // 输出 XOF_time,Encode_time,Add_time、Sbox_time、Linear_time
    std::cout << "XOF time: " << XOF_time << "s\n";
    std::cout << "Encode time: " << Encode_time << "s\n";
    std::cout << "Sbox time: " << Sbox_time << "s\n";
    std::cout << "Linear time: " << Linear_time << "s\n";
    // 计算总时间
    double Server_offtime = XOF_time + Encode_time + Sbox_time + Linear_time;
    std::cout << "Server offline total time: " << Server_offtime << "s\n";
    std::cout << "sbox_timeset: " << sbox_set << endl;
    std::cout << "linear_timeset: " << linear_set << endl;
    if (plainflag)
    {
        if (!verifyDecryption(encryptedKeyStream, KeyStream, secretKey, cmodulus, TruncWords, PlainBlock, nslots, Para_p))
        {
            std::cerr << "Decryption verification failed for KeyStream." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream." << std::endl;
    }
    // 客户端在线
    // 生成随机对称明文流，只用于测试
    random_device rd;
    vector<long> PlainStream(TruncPlainWords);
    for (int i = 0; i < TruncPlainWords; i++)
    {
        PlainStream[i] = rd() % PlainMod;
    }
    // 加密
    vector<vec_long> CipherStream(TruncWords);
    for (int i = 0; i < TruncWords; i++)
    {
        CipherStream[i].SetLength(PlainBlock);
    }
    auto start_ClientOnline = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < PlainBlock; j++)
    {
        for (int i = 0; i < TruncWords; i++)
        {
            CipherStream[i][j] = ((PlainStream[j * TruncWords + i] + KeyStream[j * TruncWords + i]) % PlainMod + PlainMod) % PlainMod;
        }
    }
    auto end_ClientOnline = std::chrono::high_resolution_clock::now();
    double Client_ontime = std::chrono::duration<double>(end_ClientOnline - start_ClientOnline).count();
    std::cout << "Client onine total time:" << std::chrono::duration<double>(end_ClientOnline - start_ClientOnline).count() << "s\n";
    double Client_totaltime = Client_offtime + Client_ontime;
    std::cout << "Client total time: " << Client_totaltime << "s\n";
    // 服务端在线
    // 同态加密
    vector<Ctxt> TruncencryptedKeyStream(encryptedKeyStream.begin(), encryptedKeyStream.begin() + TruncWords);
    vector<Ctxt> encrypedPlainStream = TruncencryptedKeyStream;
    // 对CipherStream进行编码
    vector<ZZX> encodedCipherStream(TruncWords);
    auto start_ServerOnline = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < TruncWords; i++)
    {
        cmodulus.iFFT(encodedtemp, CipherStream[i]);
        convert(encodedCipherStream[i], encodedtemp);
    }
    for (int i = 0; i < TruncWords; i++)
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
    // if (!verifyDecryption(encrypedPlainStream, PlainStream, secretKey, cmodulus, TruncWords, PlainBlock, nslots, Para_p))
    // {
    //     std::cerr << "Decryption verification failed for encrypedPlainStream." << std::endl;
    //     return 0;
    // }
    // std::cout << "Decryption verification succeeded for encrypedPlainStream." << std::endl;
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
            << std::left << std::setw(10) << Para_p
            << std::left << std::setw(10) << nslots
            << std::left << std::setw(5) << Para_bits
            << std::left << std::setw(4) << Para_c
            << std::left << std::setw(6) << Qbits
            << std::fixed << std::setprecision(3)
            << std::left << std::setw(14) << SecurityLevel
            << std::left << std::setw(7) << XOF_time
            << std::left << std::setw(10) << Encode_time
            << std::left << std::setw(8) << Sbox_time
            << std::left << std::setw(10) << Linear_time
            << std::left << std::setw(10) << Server_offtime
            << std::left << std::setw(10) << server_ontime
            << std::left << std::setw(10) << Server_totaltime
            << std::left << std::setw(20) << throughput
            << std::left << std::setw(10) << noise_budget
            << std::endl;
    outfile.close();
    std::cout << "test_pasta2_3.txt updated." << std::endl;
    return 0;
}