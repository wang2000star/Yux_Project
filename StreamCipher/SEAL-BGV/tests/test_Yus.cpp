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

//#include <NTL/ZZX.h>
// #include <NTL/GF2X.h>
//#include <helib/helib.h>
//#include <helib/ArgMap.h>
//#include <helib/DoubleCRT.h>

#include <immintrin.h> // 包含 SIMD 指令集支持

#include "utils/random_bit.hpp"
#include "Symmetric/Yus_p.hpp"
#include "utils/tool.hpp"

#include <SEAL-4.1/seal/seal.h>

using namespace std;
//using namespace helib;
//using namespace NTL;
using namespace seal;
namespace fs = std::filesystem;
/************************************************************************
  uint64_t p;          // plaintext primeplain_mod;
  uint64_t m;          // m-th cyclotomic polynomial
  uint64_t r;          // Lifting [defualt = 1]
  uint64_t bits;       // bits in the ciphertext modulus chain
  uint64_t c;          // columns in the key-switching matrix [default=2]
  uint64_t d;          // Degree of the field extension [default=1]
  uint64_t k;          // Security parameter [default=80]
  uint64_t s;          // Minimum number of slots [default=0]
************************************************************************/
// p^d = 1 mod m,d=1,slots=\phi(m)/d=\phi(m);m=65536=2^16,\phi(m)=2^15=32768
// 更一般的，应该有d|ord_p(m)，slots=\phi(m)/ord_p(m)
//!!!!!!!!!!!!!!!!
constexpr unsigned BlockByte = 36; // 分组字节长度
// ===============模式设置================
static bool Rkflag = 1;     // true/1表示乘法，false/0表示加法，指的是随机向量和密钥间的操作
static bool deflag = 0;     // true/1表示进行每一步解密验证，false/0表示不进行每一步解密验证
static bool ompflag = 0;    // true/1表示使用OpenMP并行编码，false/0表示不使用OpenMP并行编码
static bool symkeyflag = 0; // true/1表示对称密钥同态解密验证加密，false/0表示不验证
static bool plainflag = 1;  // true/1表示对称密文同态解密验证，false/0表示不验证
// 参数设置，paramMap[Nr-4][idx]
static constexpr unsigned Nr = 4; // 轮数
constexpr uint64_t idx = 0;
constexpr unsigned Sbox_depth = 1 * Nr; // S盒深度
// 当c=2时，Qbits=1.5*bits,当c=3时，Qbits=1.5*bits - 100
// 硬编码参数值
constexpr tuple<uint64_t, uint64_t, uint64_t, uint64_t> paramMap[5][8] = {
    {// Nr = 4
     // {p, log2(m), bits, c}
     {65537, 16, 216, 2},   // 0 *
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
        {65537, 16, 280, 2},  // 0 *
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
        {65537, 16, 380, 2}, // 0 *
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

constexpr uint64_t log2Para_m = get<1>(paramMap[Nr - 4][idx]) - 0;
constexpr uint64_t Para_p = get<0>(paramMap[Nr - 4][idx]);    // plaintext prime
constexpr uint64_t Para_m = 1 << log2Para_m;                  // cyclotomic polynomial
constexpr uint64_t phi_m = Para_m >> 1;                       // phi(m)=nlsots
constexpr uint64_t Para_bits = get<2>(paramMap[Nr - 4][idx]); // bits in the ciphertext modulus chain
constexpr uint64_t Para_c = get<3>(paramMap[Nr - 4][idx]);    // columns in the key-switching matrix

//!!!!!!!!!!!!!!!
constexpr unsigned PlainBlock = phi_m - 0; // 明文分组数,应该PlainBlock<=phi_m
constexpr unsigned len3 = BlockByte / 3;
// 计算 log2 的 constexpr 函数
constexpr unsigned int log2_constexpr(uint64_t n, unsigned int p = 0)
{
    return (n <= 1) ? p : log2_constexpr(n / 2, p + 1);
}
constexpr uint64_t PlainMod = Para_p;                               // 明文模数
constexpr unsigned Bytebits = log2_constexpr(PlainMod - 1) + 1; // 字节比特长度=ceil(log2(PlainMod-1))
constexpr unsigned randbits = Bytebits - 1;
constexpr unsigned BlockSize = Bytebits * BlockByte;   // 分组比特长度=BlockByte*Bytebits
constexpr unsigned NrBlockByte = BlockByte * (Nr + 1); // Nr轮分组密钥字节长度
static const uint64_t PlainByte = BlockByte * PlainBlock;  // 明文字节长度
static const uint64_t Plainbits = Bytebits * PlainByte;    // 明文比特长度
static const uint64_t PlainSize = BlockSize * PlainBlock;  // 明文比特长度

static const unsigned NonceSize = 32;                           // Nonce比特长度
static const uint64_t counter_begin = 0;                            // 计数器起始值
static const uint64_t counter_end = PlainBlock + counter_begin - 1; // 计数器结束值

YusP yusP(PlainMod); // 构建明文对称加密实例

// 函数：对多项式的每个系数乘以整数 a 并取模 c
Plaintext multiplyAndMod(const Plaintext &a, uint64_t b)
{
    int len = a.coeff_count();
    Plaintext res = a;
    for (uint64_t i = 0; i < len; ++i)
    {
        res[i] = (a[i] * b) % Para_p;
    }
    return res;
}
/*
Helper function: Prints the parameters in a SEALContext.
*/
inline void print_parameters(const seal::SEALContext &context)
{
    auto &context_data = *context.key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::bfv:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::ckks:
        scheme_name = "CKKS";
        break;
    case seal::scheme_type::bgv:
        scheme_name = "BGV";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_modulus_size = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_modulus_size - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::bfv)
    {
        std::cout << "|   plain_modulus: " << context_data.parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}
// void vector_addition_avx2(const uint64_t* __restrict a, const uint64_t* __restrict b, uint64_t* __restrict result, size_t size)
// {
//     // 检查size是否为8的倍数，确保可以正确处理AVX2的256位寄存器
//     assert(size % 8 == 0);
//     __m256i va, vb, vr;
//     for (size_t i = 0; i < size; i += 8)
//     {
//         // 加载8个uint64_t整数到AVX寄存器
//         va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + i));
//         vb = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b + i));

//         // 使用AVX2指令进行向量加法
//         vr = _mm256_add_epi32(va, vb);

//         // 存储结果回内存
//         _mm256_storeu_si256(reinterpret_cast<__m256i*>(result + i), vr);
//     }
// }
int min_noise_budget(vector<Ciphertext> &eData ,Decryptor &decryptor)
{
    int min_noise = 1000;
    for (int i = 0; i < eData.size(); i++)
    {
        int noise = decryptor.invariant_noise_budget(eData[i]);
        if (noise < min_noise)
        {
            min_noise = noise;
        }
    }
    return min_noise;
}
bool writeEncryptedSymKey(const vector<Ciphertext> &encryptedSymKey, const std::string &filename)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open())
    {
        std::cerr << "Failed to open " << filename << " for writing" << std::endl;
        return false;
    }

    for (const auto &ctxt : encryptedSymKey)
    {
        ctxt.save(out);
    }

    out.close();

    return true;
}
// void encodeTo36Ciphertext(vector<ZZX> &encData, const vector<uint64_t> &data, const EncryptedArray &ea)
// {
//     uint64_t R = data.size() / PlainByte;
//     uint64_t nCiphertext = BlockByte * R;
//     uint64_t data_size = data.size();
//     uint64_t ea_size = ea.size();
//     encData.resize(nCiphertext);
// #if (opmflag)
//     omp_set_num_threads(12); // 设置线程数为12
// #pragma omp parallel for
// #endif
//     for (uint64_t i = 0; i < BlockByte; i++)
//     {
//         vector<uint64_t> slots(ea_size, 0);
//         for (uint64_t r = 0; r < R; r++)
//         {
//             for (uint64_t j = 0; j < PlainBlock; j++)
//             {
//                 uint64_t byteIdx = j * BlockByte + i + r * PlainByte;
//                 slots[j] = data[byteIdx];
//             }
//             ea.encode(encData[r * BlockByte + i], slots);
//         }
//     }
// }

// encodeTo36Ciphertext对应的解码
void decodeTo36Ciphertext(vector<uint64_t> &data, const vector<vector<uint64_t>> &encData)
{
    uint64_t R = encData.size() / BlockByte;
    uint64_t data_size = R * PlainByte;
    data.resize(data_size);
    omp_set_num_threads(16); // 设置线程数为16
#pragma omp parallel for
    for (uint64_t j = 0; j < PlainBlock; j++)
    {
        for (uint64_t r = 0; r < R; r++)
        {
            for (uint64_t i = 0; i < BlockByte; i++)
            { // i is the ciphertext number
                // j is the block number in this ctxt
                uint64_t byteIdx = j * BlockByte + i + r * PlainByte;
                data[byteIdx] = encData[r * BlockByte + i][j];
            }
        }
    }
}
// 函数：解密并验证密文是否正确，需要解码
// 函数：解密并验证密文是否正确
bool verifyDecryption36(const std::vector<Ciphertext> &encryptedVec, const vector<uint64_t> &originalVec, 
BatchEncoder &batch_encoder, Decryptor &decryptor)
{
    auto start_decrypt = std::chrono::high_resolution_clock::now();
    vector<uint64_t> decryptedVec;
    std::vector<std::vector<uint64_t>> decryptedPolys(encryptedVec.size());
    Plaintext decrypted_result;
    std::vector<uint64_t> pod_result;
    omp_set_num_threads(16); // 设置线程数为16
#pragma omp parallel for
    for (std::size_t i = 0; i < encryptedVec.size(); ++i)
    {
        decryptor.decrypt(encryptedVec[i], decrypted_result);
        batch_encoder.decode(decrypted_result, pod_result);
        decryptedPolys[i] = pod_result;
    }
    // 解码
    decodeTo36Ciphertext(decryptedVec, decryptedPolys);
    // 验证解密结果
    bool isDecryptedVecCorrect = std::equal(decryptedVec.begin(), decryptedVec.end(), originalVec.begin());
    auto end_decrypt = std::chrono::high_resolution_clock::now();
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
void encryptSymKey(vector<Ciphertext> &encryptedSymKey, vector<Ciphertext> &encryptedSymKey01, vector<Ciphertext> &encryptedSymKey02, 
const vector<vector<uint64_t>> &SymKey, BatchEncoder &batch_encoder, Encryptor &encryptor)
{
    uint64_t nslots = SymKey[0].size();
    // 加密
    Plaintext SymKey_plain;
    for (uint64_t i = 0; i < BlockByte; i++)
    { // encrypt the encoded key
        batch_encoder.encode(SymKey[i], SymKey_plain);
        encryptor.encrypt(SymKey_plain,encryptedSymKey[i]);
    }
    for (uint64_t i = 0; i < BlockByte; i += 3)
    { // encrypt the encoded key2
        vector<uint64_t> slotsData01(nslots, (SymKey[i][0] * SymKey[i + 1][0]) % PlainMod);
        vector<uint64_t> slotsData02(nslots, (SymKey[i][0] * SymKey[i + 2][0]) % PlainMod);
        int num = i / 3;
        batch_encoder.encode(slotsData01, SymKey_plain);
        encryptor.encrypt(SymKey_plain,encryptedSymKey01[num]);
        batch_encoder.encode(slotsData02, SymKey_plain);
        encryptor.encrypt(SymKey_plain,encryptedSymKey02[num]);
    }
}
bool verify_encryptSymKey(vector<Ciphertext> &encryptedSymKey, const vector<vector<uint64_t>> &SymKey,BatchEncoder &batch_encoder, Decryptor &decryptor)
{
    auto start_decrypt = std::chrono::high_resolution_clock::now();
    vector<vector<uint64_t>> decryptedSymKey = SymKey;
    bool isDecryptedSymKeyCorrect = true;
    omp_set_num_threads(12); // 设置线程数为12
    #pragma omp parallel for
    for (uint64_t i = 0; i < BlockByte; i++)
    { // encrypt the encoded key
        Plaintext decrypted_result;
        decryptor.decrypt(encryptedSymKey[i], decrypted_result);
        vector<uint64_t> pod_result;
        batch_encoder.decode(decrypted_result, pod_result);
        decryptedSymKey[i] = pod_result;
    }
    for (int i = 0; i < BlockByte; i++)
    {
        if (decryptedSymKey[i] != SymKey[i])
        {
            isDecryptedSymKeyCorrect = false;
            break;
        }
    }
    auto end_decrypt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification finished! Time: " << elapsed_seconds.count() << "s\n";
    return isDecryptedSymKeyCorrect;
}
// Linear transformation
void HE_M2(vector<Ciphertext> &eData,Evaluator &evaluator,RelinKeys &relin_keys)
{
    vector<Ciphertext> temp = eData;
    // 0,1,2,3
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
    // 4,5,6,7
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
    // 8,9,10,11
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
    // 12,13,14,15
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
    // 16,17,18,19
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
    // 20,21,22,23
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
    // 24,25,26,27
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
    // 28,29,30,31
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
    // 32,33,34,35
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

    eData[0] = temp1_2_3;
    evaluator.add_inplace(eData[0], temp4_5_6_7);
    evaluator.add_inplace(eData[0], temp9_11);
    evaluator.add_inplace(eData[0], temp[13]);
    evaluator.add_inplace(eData[0], temp16_17_19);
    evaluator.add_inplace(eData[0], temp20_22);
    evaluator.add_inplace(eData[0], temp25_26_27);
    evaluator.add_inplace(eData[0], temp29_30);
    evaluator.add_inplace(eData[0], temp33_35);
    //evaluator.relinearize_inplace(eData[0], relin_keys);

    evaluator.add_inplace(eData[1], temp[3]);
    evaluator.add_inplace(eData[1], temp4_6_7);
    evaluator.add_inplace(eData[1], temp8_9_10);
    evaluator.add_inplace(eData[1], temp[12]);
    evaluator.add_inplace(eData[1], temp17_18_19);
    evaluator.add_inplace(eData[1], temp21_22);
    evaluator.add_inplace(eData[1], temp24_27);
    evaluator.add_inplace(eData[1], temp[31]);
    evaluator.add_inplace(eData[1], temp33_34_35);
    //evaluator.relinearize_inplace(eData[1], relin_keys);

    eData[2] = temp0_3;
    evaluator.add_inplace(eData[2], temp[7]);
    evaluator.add_inplace(eData[2], temp9_10_11);
    evaluator.add_inplace(eData[2], temp13_14_15);
    evaluator.add_inplace(eData[2], temp17_18);
    evaluator.add_inplace(eData[2], temp21_22_23);
    evaluator.add_inplace(eData[2], temp[26]);
    evaluator.add_inplace(eData[2], temp[28]);
    evaluator.add_inplace(eData[2], temp32_33_34);
    //evaluator.relinearize_inplace(eData[2], relin_keys);

    eData[3] = temp0_2;
    evaluator.add_inplace(eData[3], temp4_5_6_7);
    evaluator.add_inplace(eData[3], temp8_9_10);
    evaluator.add_inplace(eData[3], temp12_14);
    evaluator.add_inplace(eData[3], temp16_19);
    evaluator.add_inplace(eData[3], temp20_22_23);
    evaluator.add_inplace(eData[3], temp[25]);
    evaluator.add_inplace(eData[3], temp28_29_30);
    evaluator.add_inplace(eData[3], temp32_33);
    //evaluator.relinearize_inplace(eData[3], relin_keys);

    evaluator.add_inplace(eData[4], temp0_1_2);
    evaluator.add_inplace(eData[4], temp6_7);
    evaluator.add_inplace(eData[4], temp9_10_11);
    evaluator.add_inplace(eData[4], temp12_13_15);
    evaluator.add_inplace(eData[4], temp20_21_22);
    evaluator.add_inplace(eData[4], temp24_25_27);
    evaluator.add_inplace(eData[4], temp[30]);
    evaluator.add_inplace(eData[4], temp[34]);
    //evaluator.relinearize_inplace(eData[4], relin_keys);

    eData[5] = temp0_1_3;
    evaluator.add_inplace(eData[5], temp[6]);
    evaluator.add_inplace(eData[5], temp[10]);
    evaluator.add_inplace(eData[5], temp12_13_14);
    evaluator.add_inplace(eData[5], temp16_17_18);
    evaluator.add_inplace(eData[5], temp20_21);
    evaluator.add_inplace(eData[5], temp24_25_26);
    evaluator.add_inplace(eData[5], temp29_31);
    evaluator.add_inplace(eData[5], temp[35]);
    //evaluator.relinearize_inplace(eData[5], relin_keys);

    eData[6] = temp0_3;
    evaluator.add_inplace(eData[6], temp5_7);
    evaluator.add_inplace(eData[6], temp8_9_10_11);
    evaluator.add_inplace(eData[6], temp12_13_15);
    evaluator.add_inplace(eData[6], temp17_19);
    evaluator.add_inplace(eData[6], temp22_23);
    evaluator.add_inplace(eData[6], temp25_26);
    evaluator.add_inplace(eData[6], temp28_31);
    evaluator.add_inplace(eData[6], temp32_33_35);
    //evaluator.relinearize_inplace(eData[6], relin_keys);

    evaluator.add_inplace(eData[7], temp1_3);
    evaluator.add_inplace(eData[7], temp4_5);
    evaluator.add_inplace(eData[7], temp9_10);
    evaluator.add_inplace(eData[7], temp12_13_14_15);
    evaluator.add_inplace(eData[7], temp16_18);
    evaluator.add_inplace(eData[7], temp[23]);
    evaluator.add_inplace(eData[7], temp24_25_27);
    evaluator.add_inplace(eData[7], temp28_30);
    evaluator.add_inplace(eData[7], temp[33]);
    //evaluator.relinearize_inplace(eData[7], relin_keys);

    eData[8] = temp2_3;
    evaluator.add_inplace(eData[8], temp4_6);
    evaluator.add_inplace(eData[8], temp[9]);
    evaluator.add_inplace(eData[8], temp13_15);
    evaluator.add_inplace(eData[8], temp16_17_19);
    evaluator.add_inplace(eData[8], temp20_21_23);
    evaluator.add_inplace(eData[8], temp24_27);
    evaluator.add_inplace(eData[8], temp28_29);
    evaluator.add_inplace(eData[8], temp32_34);
    //evaluator.relinearize_inplace(eData[8], relin_keys);

    eData[9] = temp0_2_3;
    evaluator.add_inplace(eData[9], temp[6]);
    evaluator.add_inplace(eData[9], temp8_10_11);
    evaluator.add_inplace(eData[9], temp12_13_14_15);
    evaluator.add_inplace(eData[9], temp16_18);
    evaluator.add_inplace(eData[9], temp20_22);
    evaluator.add_inplace(eData[9], temp25_26);
    evaluator.add_inplace(eData[9], temp28_29_31);
    evaluator.add_inplace(eData[9], temp34_35);
    //evaluator.relinearize_inplace(eData[9], relin_keys);

    evaluator.add_inplace(eData[10], temp[0]);
    evaluator.add_inplace(eData[10], temp4_6_7);
    evaluator.add_inplace(eData[10], temp[8]);
    evaluator.add_inplace(eData[10], temp12_13_15);
    evaluator.add_inplace(eData[10], temp16_17_18_19);
    evaluator.add_inplace(eData[10], temp[21]);
    evaluator.add_inplace(eData[10], temp26_27);
    evaluator.add_inplace(eData[10], temp28_30_31);
    evaluator.add_inplace(eData[10], temp[33]);
    //evaluator.relinearize_inplace(eData[10], relin_keys);

    eData[11] = temp[1];
    evaluator.add_inplace(eData[11], temp5_6_7);
    evaluator.add_inplace(eData[11], temp[9]);
    evaluator.add_inplace(eData[11], temp[12]);
    evaluator.add_inplace(eData[11], temp16_18_19);
    evaluator.add_inplace(eData[11], temp20_22_23);
    evaluator.add_inplace(eData[11], temp24_26_27);
    evaluator.add_inplace(eData[11], temp30_31);
    evaluator.add_inplace(eData[11], temp32_35);
    //evaluator.relinearize_inplace(eData[11], relin_keys);

    eData[12] = temp1_2_3;
    evaluator.add_inplace(eData[12], temp5_6);
    evaluator.add_inplace(eData[12], temp9_11);
    evaluator.add_inplace(eData[12], temp13_14_15);
    evaluator.add_inplace(eData[12], temp16_17_18_19);
    evaluator.add_inplace(eData[12], temp21_23);
    evaluator.add_inplace(eData[12], temp[25]);
    evaluator.add_inplace(eData[12], temp28_29_31);
    evaluator.add_inplace(eData[12], temp32_34);
    //evaluator.relinearize_inplace(eData[12], relin_keys);

    evaluator.add_inplace(eData[13], temp0_3);
    evaluator.add_inplace(eData[13], temp[7]);
    evaluator.add_inplace(eData[13], temp9_10_11);
    evaluator.add_inplace(eData[13], temp[15]);
    evaluator.add_inplace(eData[13], temp16_18_19);
    evaluator.add_inplace(eData[13], temp20_21_22);
    evaluator.add_inplace(eData[13], temp[24]);
    evaluator.add_inplace(eData[13], temp29_30_31);
    evaluator.add_inplace(eData[13], temp33_34);
    //evaluator.relinearize_inplace(eData[13], relin_keys);

    eData[14] = temp[2];
    evaluator.add_inplace(eData[14], temp[4]);
    evaluator.add_inplace(eData[14], temp8_9_10);
    evaluator.add_inplace(eData[14], temp12_15);
    evaluator.add_inplace(eData[14], temp[19]);
    evaluator.add_inplace(eData[14], temp21_22_23);
    evaluator.add_inplace(eData[14], temp25_26_27);
    evaluator.add_inplace(eData[14], temp29_30);
    evaluator.add_inplace(eData[14], temp33_34_35);
    //evaluator.relinearize_inplace(eData[14], relin_keys);

    eData[15] = temp[1];
    evaluator.add_inplace(eData[15], temp4_5_6);
    evaluator.add_inplace(eData[15], temp8_9);
    evaluator.add_inplace(eData[15], temp12_14);
    evaluator.add_inplace(eData[15], temp16_17_18_19);
    evaluator.add_inplace(eData[15], temp20_21_22);
    evaluator.add_inplace(eData[15], temp24_26);
    evaluator.add_inplace(eData[15], temp28_31);
    evaluator.add_inplace(eData[15], temp32_34_35);
    //evaluator.relinearize_inplace(eData[15], relin_keys);

    evaluator.add_inplace(eData[16], temp0_1_3);
    evaluator.add_inplace(eData[16], temp[6]);
    evaluator.add_inplace(eData[16], temp[10]);
    evaluator.add_inplace(eData[16], temp12_13_14);
    evaluator.add_inplace(eData[16], temp18_19);
    evaluator.add_inplace(eData[16], temp21_22_23);
    evaluator.add_inplace(eData[16], temp24_25_27);
    evaluator.add_inplace(eData[16], temp32_33_34);
    //evaluator.relinearize_inplace(eData[16], relin_keys);

    eData[17] = temp0_1_2;
    evaluator.add_inplace(eData[17], temp5_7);
    evaluator.add_inplace(eData[17], temp[11]);
    evaluator.add_inplace(eData[17], temp12_13_15);
    evaluator.add_inplace(eData[17], temp[18]);
    evaluator.add_inplace(eData[17], temp[22]);
    evaluator.add_inplace(eData[17], temp24_25_26);
    evaluator.add_inplace(eData[17], temp28_29_30);
    evaluator.add_inplace(eData[17], temp32_33);
    //evaluator.relinearize_inplace(eData[17], relin_keys);

    eData[18] = temp1_2;
    evaluator.add_inplace(eData[18], temp4_7);
    evaluator.add_inplace(eData[18], temp8_9_11);
    evaluator.add_inplace(eData[18], temp12_15);
    evaluator.add_inplace(eData[18], temp17_19);
    evaluator.add_inplace(eData[18], temp20_21_22_23);
    evaluator.add_inplace(eData[18], temp24_25_27);
    evaluator.add_inplace(eData[18], temp29_31);
    evaluator.add_inplace(eData[18], temp34_35);
    //evaluator.relinearize_inplace(eData[18], relin_keys);

    evaluator.add_inplace(eData[19], temp0_1_3);
    evaluator.add_inplace(eData[19], temp4_6);
    evaluator.add_inplace(eData[19], temp[9]);
    evaluator.add_inplace(eData[19], temp13_15);
    evaluator.add_inplace(eData[19], temp16_17);
    evaluator.add_inplace(eData[19], temp21_22);
    evaluator.add_inplace(eData[19], temp24_25_26_27);
    evaluator.add_inplace(eData[19], temp28_30);
    evaluator.add_inplace(eData[19], temp[35]);
    //evaluator.relinearize_inplace(eData[19], relin_keys);

    eData[20] = temp0_3;
    evaluator.add_inplace(eData[20], temp4_5);
    evaluator.add_inplace(eData[20], temp8_10);
    evaluator.add_inplace(eData[20], temp14_15);
    evaluator.add_inplace(eData[20], temp16_18);
    evaluator.add_inplace(eData[20], temp[21]);
    evaluator.add_inplace(eData[20], temp25_27);
    evaluator.add_inplace(eData[20], temp28_29_31);
    evaluator.add_inplace(eData[20], temp32_33_35);
    //evaluator.relinearize_inplace(eData[20], relin_keys);

    eData[21] = temp1_2;
    evaluator.add_inplace(eData[21], temp4_5_7);
    evaluator.add_inplace(eData[21], temp10_11);
    evaluator.add_inplace(eData[21], temp12_14_15);
    evaluator.add_inplace(eData[21], temp[18]);
    evaluator.add_inplace(eData[21], temp20_22_23);
    evaluator.add_inplace(eData[21], temp24_25_26_27);
    evaluator.add_inplace(eData[21], temp28_30);
    evaluator.add_inplace(eData[21], temp32_34);
    //evaluator.relinearize_inplace(eData[21], relin_keys);

    evaluator.add_inplace(eData[22], temp2_3);
    evaluator.add_inplace(eData[22], temp4_6_7);
    evaluator.add_inplace(eData[22], temp[9]);
    evaluator.add_inplace(eData[22], temp[12]);
    evaluator.add_inplace(eData[22], temp16_18_19);
    evaluator.add_inplace(eData[22], temp[20]);
    evaluator.add_inplace(eData[22], temp24_25_27);
    evaluator.add_inplace(eData[22], temp28_29_30_31);
    evaluator.add_inplace(eData[22], temp[33]);
    //evaluator.relinearize_inplace(eData[22], relin_keys);

    eData[23] = temp0_2_3;
    evaluator.add_inplace(eData[23], temp6_7);
    evaluator.add_inplace(eData[23], temp8_11);
    evaluator.add_inplace(eData[23], temp[13]);
    evaluator.add_inplace(eData[23], temp17_18_19);
    evaluator.add_inplace(eData[23], temp[21]);
    evaluator.add_inplace(eData[23], temp[24]);
    evaluator.add_inplace(eData[23], temp28_30_31);
    evaluator.add_inplace(eData[23], temp32_34_35);
    //evaluator.relinearize_inplace(eData[23], relin_keys);

    eData[24] = temp[1];
    evaluator.add_inplace(eData[24], temp4_5_7);
    evaluator.add_inplace(eData[24], temp8_10);
    evaluator.add_inplace(eData[24], temp13_14_15);
    evaluator.add_inplace(eData[24], temp17_18);
    evaluator.add_inplace(eData[24], temp21_23);
    evaluator.add_inplace(eData[24], temp25_26_27);
    evaluator.add_inplace(eData[24], temp28_29_30_31);
    evaluator.add_inplace(eData[24], temp33_35);
    //evaluator.relinearize_inplace(eData[24], relin_keys);

    evaluator.add_inplace(eData[25], temp[0]);
    evaluator.add_inplace(eData[25], temp5_6_7);
    evaluator.add_inplace(eData[25], temp9_10);
    evaluator.add_inplace(eData[25], temp12_15);
    evaluator.add_inplace(eData[25], temp[19]);
    evaluator.add_inplace(eData[25], temp21_22_23);
    evaluator.add_inplace(eData[25], temp[27]);
    evaluator.add_inplace(eData[25], temp28_30_31);
    evaluator.add_inplace(eData[25], temp32_33_34);
    //evaluator.relinearize_inplace(eData[25], relin_keys);

    eData[26] = temp1_2_3;
    evaluator.add_inplace(eData[26], temp5_6);
    evaluator.add_inplace(eData[26], temp9_10_11);
    evaluator.add_inplace(eData[26], temp[14]);
    evaluator.add_inplace(eData[26], temp[16]);
    evaluator.add_inplace(eData[26], temp20_21_22);
    evaluator.add_inplace(eData[26], temp24_27);
    evaluator.add_inplace(eData[26], temp[31]);
    evaluator.add_inplace(eData[26], temp33_34_35);
    //evaluator.relinearize_inplace(eData[26], relin_keys);

    eData[27] = temp0_2;
    evaluator.add_inplace(eData[27], temp4_7);
    evaluator.add_inplace(eData[27], temp8_10_11);
    evaluator.add_inplace(eData[27], temp[13]);
    evaluator.add_inplace(eData[27], temp16_17_18);
    evaluator.add_inplace(eData[27], temp20_21);
    evaluator.add_inplace(eData[27], temp24_26);
    evaluator.add_inplace(eData[27], temp28_29_30_31);
    evaluator.add_inplace(eData[27], temp32_33_34);
    //evaluator.relinearize_inplace(eData[27], relin_keys);

    evaluator.add_inplace(eData[28], temp0_1_3);
    evaluator.add_inplace(eData[28], temp8_9_10);
    evaluator.add_inplace(eData[28], temp12_13_15);
    evaluator.add_inplace(eData[28], temp[18]);
    evaluator.add_inplace(eData[28], temp[22]);
    evaluator.add_inplace(eData[28], temp24_25_26);
    evaluator.add_inplace(eData[28], temp30_31);
    evaluator.add_inplace(eData[28], temp33_34_35);
    //evaluator.relinearize_inplace(eData[28], relin_keys);

    eData[29] = temp0_1_2;
    evaluator.add_inplace(eData[29], temp4_5_6);
    evaluator.add_inplace(eData[29], temp8_9);
    evaluator.add_inplace(eData[29], temp12_13_14);
    evaluator.add_inplace(eData[29], temp17_19);
    evaluator.add_inplace(eData[29], temp[23]);
    evaluator.add_inplace(eData[29], temp24_25_27);
    evaluator.add_inplace(eData[29], temp[30]);
    evaluator.add_inplace(eData[29], temp[34]);
    //evaluator.relinearize_inplace(eData[29], relin_keys);

    eData[30] = temp0_1_3;
    evaluator.add_inplace(eData[30], temp5_7);
    evaluator.add_inplace(eData[30], temp10_11);
    evaluator.add_inplace(eData[30], temp13_14);
    evaluator.add_inplace(eData[30], temp16_19);
    evaluator.add_inplace(eData[30], temp20_21_23);
    evaluator.add_inplace(eData[30], temp24_27);
    evaluator.add_inplace(eData[30], temp29_31);
    evaluator.add_inplace(eData[30], temp32_33_34_35);
    //evaluator.relinearize_inplace(eData[30], relin_keys);

    evaluator.add_inplace(eData[31], temp0_1_2_3);
    evaluator.add_inplace(eData[31], temp4_6);
    evaluator.add_inplace(eData[31], temp[11]);
    evaluator.add_inplace(eData[31], temp12_13_15);
    evaluator.add_inplace(eData[31], temp16_18);
    evaluator.add_inplace(eData[31], temp[21]);
    evaluator.add_inplace(eData[31], temp25_27);
    evaluator.add_inplace(eData[31], temp28_29);
    evaluator.add_inplace(eData[31], temp33_34);
    //evaluator.relinearize_inplace(eData[31], relin_keys);

    eData[32] = temp1_3;
    evaluator.add_inplace(eData[32], temp4_5_7);
    evaluator.add_inplace(eData[32], temp8_9_11);
    evaluator.add_inplace(eData[32], temp12_15);
    evaluator.add_inplace(eData[32], temp16_17);
    evaluator.add_inplace(eData[32], temp20_22);
    evaluator.add_inplace(eData[32], temp26_27);
    evaluator.add_inplace(eData[32], temp28_30);
    evaluator.add_inplace(eData[32], temp[33]);
    //evaluator.relinearize_inplace(eData[32], relin_keys);

    eData[33] = temp0_1_2_3;
    evaluator.add_inplace(eData[33], temp4_6);
    evaluator.add_inplace(eData[33], temp8_10);
    evaluator.add_inplace(eData[33], temp13_14);
    evaluator.add_inplace(eData[33], temp16_17_19);
    evaluator.add_inplace(eData[33], temp22_23);
    evaluator.add_inplace(eData[33], temp24_26_27);
    evaluator.add_inplace(eData[33], temp[30]);
    evaluator.add_inplace(eData[33], temp32_34_35);
    //evaluator.relinearize_inplace(eData[33], relin_keys);

    evaluator.add_inplace(eData[34], temp0_1_3);
    evaluator.add_inplace(eData[34], temp4_5_6_7);
    evaluator.add_inplace(eData[34], temp[9]);
    evaluator.add_inplace(eData[34], temp14_15);
    evaluator.add_inplace(eData[34], temp16_18_19);
    evaluator.add_inplace(eData[34], temp[21]);
    evaluator.add_inplace(eData[34], temp[24]);
    evaluator.add_inplace(eData[34], temp28_30_31);
    evaluator.add_inplace(eData[34], temp[32]);
    //evaluator.relinearize_inplace(eData[34], relin_keys);

    eData[35] = temp[0];
    evaluator.add_inplace(eData[35], temp4_6_7);
    evaluator.add_inplace(eData[35], temp8_10_11);
    evaluator.add_inplace(eData[35], temp12_14_15);
    evaluator.add_inplace(eData[35], temp18_19);
    evaluator.add_inplace(eData[35], temp20_23);
    evaluator.add_inplace(eData[35], temp[25]);
    evaluator.add_inplace(eData[35], temp29_30_31);
    evaluator.add_inplace(eData[35], temp[33]);
    //evaluator.relinearize_inplace(eData[35], relin_keys);

    if (1)
    {
        temp = eData;
        // 0,1,2,3
    temp0_1 = temp[0];
    evaluator.add_inplace(temp0_1, temp[1]);
    temp0_2 = temp[0];
    evaluator.add_inplace(temp0_2, temp[2]);
    temp0_3 = temp[0];
    evaluator.add_inplace(temp0_3, temp[3]);
    temp1_2 = temp[1];
    evaluator.add_inplace(temp1_2, temp[2]);
    temp1_3 = temp[1];
    evaluator.add_inplace(temp1_3, temp[3]);
    temp2_3 = temp[2];
    evaluator.add_inplace(temp2_3, temp[3]);
    temp0_1_2 = temp0_1;
    evaluator.add_inplace(temp0_1_2, temp[2]);
    temp0_1_3 = temp0_1;
    evaluator.add_inplace(temp0_1_3, temp[3]);
    temp0_2_3 = temp0_2;
    evaluator.add_inplace(temp0_2_3, temp[3]);
    temp1_2_3 = temp1_2;
    evaluator.add_inplace(temp1_2_3, temp[3]);
    temp0_1_2_3 = temp0_1_2;
    evaluator.add_inplace(temp0_1_2_3, temp[3]);
    // 4,5,6,7
    temp4_5 = temp[4];
    evaluator.add_inplace(temp4_5, temp[5]);
    temp4_6 = temp[4];
    evaluator.add_inplace(temp4_6, temp[6]);
    temp4_7 = temp[4];
    evaluator.add_inplace(temp4_7, temp[7]);
    temp5_6 = temp[5];
    evaluator.add_inplace(temp5_6, temp[6]);
    temp5_7 = temp[5];
    evaluator.add_inplace(temp5_7, temp[7]);
    temp6_7 = temp[6];
    evaluator.add_inplace(temp6_7, temp[7]);
    temp4_5_6 = temp4_5;
    evaluator.add_inplace(temp4_5_6, temp[6]);
    temp4_5_7 = temp4_5;
    evaluator.add_inplace(temp4_5_7, temp[7]);
    temp4_6_7 = temp4_6;
    evaluator.add_inplace(temp4_6_7, temp[7]);
    temp5_6_7 = temp5_6;
    evaluator.add_inplace(temp5_6_7, temp[7]);
    temp4_5_6_7 = temp4_5_6;
    evaluator.add_inplace(temp4_5_6_7, temp[7]);
    // 8,9,10,11
    temp8_9 = temp[8];
    evaluator.add_inplace(temp8_9, temp[9]);
    temp8_10 = temp[8];
    evaluator.add_inplace(temp8_10, temp[10]);
    temp8_11 = temp[8];
    evaluator.add_inplace(temp8_11, temp[11]);
    temp9_10 = temp[9];
    evaluator.add_inplace(temp9_10, temp[10]);
    temp9_11 = temp[9];
    evaluator.add_inplace(temp9_11, temp[11]);
    temp10_11 = temp[10];
    evaluator.add_inplace(temp10_11, temp[11]);
    temp8_9_10 = temp8_9;
    evaluator.add_inplace(temp8_9_10, temp[10]);
    temp8_9_11 = temp8_9;
    evaluator.add_inplace(temp8_9_11, temp[11]);
    temp8_10_11 = temp8_10;
    evaluator.add_inplace(temp8_10_11, temp[11]);
    temp9_10_11 = temp9_10;
    evaluator.add_inplace(temp9_10_11, temp[11]);
    temp8_9_10_11 = temp8_9_10;
    evaluator.add_inplace(temp8_9_10_11, temp[11]);
    // 12,13,14,15
    temp12_13 = temp[12];
    evaluator.add_inplace(temp12_13, temp[13]);
    temp12_14 = temp[12];
    evaluator.add_inplace(temp12_14, temp[14]);
    temp12_15 = temp[12];
    evaluator.add_inplace(temp12_15, temp[15]);
    temp13_14 = temp[13];
    evaluator.add_inplace(temp13_14, temp[14]);
    temp13_15 = temp[13];
    evaluator.add_inplace(temp13_15, temp[15]);
    temp14_15 = temp[14];
    evaluator.add_inplace(temp14_15, temp[15]);
    temp12_13_14 = temp12_13;
    evaluator.add_inplace(temp12_13_14, temp[14]);
    temp12_13_15 = temp12_13;
    evaluator.add_inplace(temp12_13_15, temp[15]);
    temp12_14_15 = temp12_14;
    evaluator.add_inplace(temp12_14_15, temp[15]);
    temp13_14_15 = temp13_14;
    evaluator.add_inplace(temp13_14_15, temp[15]);
    temp12_13_14_15 = temp12_13_14;
    evaluator.add_inplace(temp12_13_14_15, temp[15]);
    // 16,17,18,19
    temp16_17 = temp[16];
    evaluator.add_inplace(temp16_17, temp[17]);
    temp16_18 = temp[16];
    evaluator.add_inplace(temp16_18, temp[18]);
    temp16_19 = temp[16];
    evaluator.add_inplace(temp16_19, temp[19]);
    temp17_18 = temp[17];
    evaluator.add_inplace(temp17_18, temp[18]);
    temp17_19 = temp[17];
    evaluator.add_inplace(temp17_19, temp[19]);
    temp18_19 = temp[18];
    evaluator.add_inplace(temp18_19, temp[19]);
    temp16_17_18 = temp16_17;
    evaluator.add_inplace(temp16_17_18, temp[18]);
    temp16_17_19 = temp16_17;
    evaluator.add_inplace(temp16_17_19, temp[19]);
    temp16_18_19 = temp16_18;
    evaluator.add_inplace(temp16_18_19, temp[19]);
    temp17_18_19 = temp17_18;
    evaluator.add_inplace(temp17_18_19, temp[19]);
    temp16_17_18_19 = temp16_17_18;
    evaluator.add_inplace(temp16_17_18_19, temp[19]);
    // 20,21,22,23
    temp20_21 = temp[20];
    evaluator.add_inplace(temp20_21, temp[21]);
    temp20_22 = temp[20];
    evaluator.add_inplace(temp20_22, temp[22]);
    temp20_23 = temp[20];
    evaluator.add_inplace(temp20_23, temp[23]);
    temp21_22 = temp[21];
    evaluator.add_inplace(temp21_22, temp[22]);
    temp21_23 = temp[21];
    evaluator.add_inplace(temp21_23, temp[23]);
    temp22_23 = temp[22];
    evaluator.add_inplace(temp22_23, temp[23]);
    temp20_21_22 = temp20_21;
    evaluator.add_inplace(temp20_21_22, temp[22]);
    temp20_21_23 = temp20_21;
    evaluator.add_inplace(temp20_21_23, temp[23]);
    temp20_22_23 = temp20_22;
    evaluator.add_inplace(temp20_22_23, temp[23]);
    temp21_22_23 = temp21_22;
    evaluator.add_inplace(temp21_22_23, temp[23]);
    temp20_21_22_23 = temp20_21_22;
    evaluator.add_inplace(temp20_21_22_23, temp[23]);
    // 24,25,26,27
    temp24_25 = temp[24];
    evaluator.add_inplace(temp24_25, temp[25]);
    temp24_26 = temp[24];
    evaluator.add_inplace(temp24_26, temp[26]);
    temp24_27 = temp[24];
    evaluator.add_inplace(temp24_27, temp[27]);
    temp25_26 = temp[25];
    evaluator.add_inplace(temp25_26, temp[26]);
    temp25_27 = temp[25];
    evaluator.add_inplace(temp25_27, temp[27]);
    temp26_27 = temp[26];
    evaluator.add_inplace(temp26_27, temp[27]);
    temp24_25_26 = temp24_25;
    evaluator.add_inplace(temp24_25_26, temp[26]);
    temp24_25_27 = temp24_25;
    evaluator.add_inplace(temp24_25_27, temp[27]);
    temp24_26_27 = temp24_26;
    evaluator.add_inplace(temp24_26_27, temp[27]);
    temp25_26_27 = temp25_26;
    evaluator.add_inplace(temp25_26_27, temp[27]);
    temp24_25_26_27 = temp24_25_26;
    evaluator.add_inplace(temp24_25_26_27, temp[27]);
    // 28,29,30,31
    temp28_29 = temp[28];
    evaluator.add_inplace(temp28_29, temp[29]);
    temp28_30 = temp[28];
    evaluator.add_inplace(temp28_30, temp[30]);
    temp28_31 = temp[28];
    evaluator.add_inplace(temp28_31, temp[31]);
    temp29_30 = temp[29];
    evaluator.add_inplace(temp29_30, temp[30]);
    temp29_31 = temp[29];
    evaluator.add_inplace(temp29_31, temp[31]);
    temp30_31 = temp[30];
    evaluator.add_inplace(temp30_31, temp[31]);
    temp28_29_30 = temp28_29;
    evaluator.add_inplace(temp28_29_30, temp[30]);
    temp28_29_31 = temp28_29;
    evaluator.add_inplace(temp28_29_31, temp[31]);
    temp28_30_31 = temp28_30;
    evaluator.add_inplace(temp28_30_31, temp[31]);
    temp29_30_31 = temp29_30;
    evaluator.add_inplace(temp29_30_31, temp[31]);
    temp28_29_30_31 = temp28_29_30;
    evaluator.add_inplace(temp28_29_30_31, temp[31]);
    // 32,33,34,35
    temp32_33 = temp[32];
    evaluator.add_inplace(temp32_33, temp[33]);
    temp32_34 = temp[32];
    evaluator.add_inplace(temp32_34, temp[34]);
    temp32_35 = temp[32];
    evaluator.add_inplace(temp32_35, temp[35]);
    temp33_34 = temp[33];
    evaluator.add_inplace(temp33_34, temp[34]);
    temp33_35 = temp[33];
    evaluator.add_inplace(temp33_35, temp[35]);
    temp34_35 = temp[34];
    evaluator.add_inplace(temp34_35, temp[35]);
    temp32_33_34 = temp32_33;
    evaluator.add_inplace(temp32_33_34, temp[34]);
    temp32_33_35 = temp32_33;
    evaluator.add_inplace(temp32_33_35, temp[35]);
    temp32_34_35 = temp32_34;
    evaluator.add_inplace(temp32_34_35, temp[35]);
    temp33_34_35 = temp33_34;
    evaluator.add_inplace(temp33_34_35, temp[35]);
    temp32_33_34_35 = temp32_33_34;
    evaluator.add_inplace(temp32_33_34_35, temp[35]);

        eData[0] = temp1_2_3;
    evaluator.add_inplace(eData[0], temp4_5_6_7);
    evaluator.add_inplace(eData[0], temp9_11);
    evaluator.add_inplace(eData[0], temp[13]);
    evaluator.add_inplace(eData[0], temp16_17_19);
    evaluator.add_inplace(eData[0], temp20_22);
    evaluator.add_inplace(eData[0], temp25_26_27);
    evaluator.add_inplace(eData[0], temp29_30);
    evaluator.add_inplace(eData[0], temp33_35);
    //evaluator.relinearize_inplace(eData[0], relin_keys);

    evaluator.add_inplace(eData[1], temp[3]);
    evaluator.add_inplace(eData[1], temp4_6_7);
    evaluator.add_inplace(eData[1], temp8_9_10);
    evaluator.add_inplace(eData[1], temp[12]);
    evaluator.add_inplace(eData[1], temp17_18_19);
    evaluator.add_inplace(eData[1], temp21_22);
    evaluator.add_inplace(eData[1], temp24_27);
    evaluator.add_inplace(eData[1], temp[31]);
    evaluator.add_inplace(eData[1], temp33_34_35);
    //evaluator.relinearize_inplace(eData[1], relin_keys);

    eData[2] = temp0_3;
    evaluator.add_inplace(eData[2], temp[7]);
    evaluator.add_inplace(eData[2], temp9_10_11);
    evaluator.add_inplace(eData[2], temp13_14_15);
    evaluator.add_inplace(eData[2], temp17_18);
    evaluator.add_inplace(eData[2], temp21_22_23);
    evaluator.add_inplace(eData[2], temp[26]);
    evaluator.add_inplace(eData[2], temp[28]);
    evaluator.add_inplace(eData[2], temp32_33_34);
    //evaluator.relinearize_inplace(eData[2], relin_keys);

    eData[3] = temp0_2;
    evaluator.add_inplace(eData[3], temp4_5_6_7);
    evaluator.add_inplace(eData[3], temp8_9_10);
    evaluator.add_inplace(eData[3], temp12_14);
    evaluator.add_inplace(eData[3], temp16_19);
    evaluator.add_inplace(eData[3], temp20_22_23);
    evaluator.add_inplace(eData[3], temp[25]);
    evaluator.add_inplace(eData[3], temp28_29_30);
    evaluator.add_inplace(eData[3], temp32_33);
    //evaluator.relinearize_inplace(eData[3], relin_keys);

    evaluator.add_inplace(eData[4], temp0_1_2);
    evaluator.add_inplace(eData[4], temp6_7);
    evaluator.add_inplace(eData[4], temp9_10_11);
    evaluator.add_inplace(eData[4], temp12_13_15);
    evaluator.add_inplace(eData[4], temp20_21_22);
    evaluator.add_inplace(eData[4], temp24_25_27);
    evaluator.add_inplace(eData[4], temp[30]);
    evaluator.add_inplace(eData[4], temp[34]);
    //evaluator.relinearize_inplace(eData[4], relin_keys);

    eData[5] = temp0_1_3;
    evaluator.add_inplace(eData[5], temp[6]);
    evaluator.add_inplace(eData[5], temp[10]);
    evaluator.add_inplace(eData[5], temp12_13_14);
    evaluator.add_inplace(eData[5], temp16_17_18);
    evaluator.add_inplace(eData[5], temp20_21);
    evaluator.add_inplace(eData[5], temp24_25_26);
    evaluator.add_inplace(eData[5], temp29_31);
    evaluator.add_inplace(eData[5], temp[35]);
    //evaluator.relinearize_inplace(eData[5], relin_keys);

    eData[6] = temp0_3;
    evaluator.add_inplace(eData[6], temp5_7);
    evaluator.add_inplace(eData[6], temp8_9_10_11);
    evaluator.add_inplace(eData[6], temp12_13_15);
    evaluator.add_inplace(eData[6], temp17_19);
    evaluator.add_inplace(eData[6], temp22_23);
    evaluator.add_inplace(eData[6], temp25_26);
    evaluator.add_inplace(eData[6], temp28_31);
    evaluator.add_inplace(eData[6], temp32_33_35);
    //evaluator.relinearize_inplace(eData[6], relin_keys);

    evaluator.add_inplace(eData[7], temp1_3);
    evaluator.add_inplace(eData[7], temp4_5);
    evaluator.add_inplace(eData[7], temp9_10);
    evaluator.add_inplace(eData[7], temp12_13_14_15);
    evaluator.add_inplace(eData[7], temp16_18);
    evaluator.add_inplace(eData[7], temp[23]);
    evaluator.add_inplace(eData[7], temp24_25_27);
    evaluator.add_inplace(eData[7], temp28_30);
    evaluator.add_inplace(eData[7], temp[33]);
    //evaluator.relinearize_inplace(eData[7], relin_keys);

    eData[8] = temp2_3;
    evaluator.add_inplace(eData[8], temp4_6);
    evaluator.add_inplace(eData[8], temp[9]);
    evaluator.add_inplace(eData[8], temp13_15);
    evaluator.add_inplace(eData[8], temp16_17_19);
    evaluator.add_inplace(eData[8], temp20_21_23);
    evaluator.add_inplace(eData[8], temp24_27);
    evaluator.add_inplace(eData[8], temp28_29);
    evaluator.add_inplace(eData[8], temp32_34);
    //evaluator.relinearize_inplace(eData[8], relin_keys);

    eData[9] = temp0_2_3;
    evaluator.add_inplace(eData[9], temp[6]);
    evaluator.add_inplace(eData[9], temp8_10_11);
    evaluator.add_inplace(eData[9], temp12_13_14_15);
    evaluator.add_inplace(eData[9], temp16_18);
    evaluator.add_inplace(eData[9], temp20_22);
    evaluator.add_inplace(eData[9], temp25_26);
    evaluator.add_inplace(eData[9], temp28_29_31);
    evaluator.add_inplace(eData[9], temp34_35);
    //evaluator.relinearize_inplace(eData[9], relin_keys);

    evaluator.add_inplace(eData[10], temp[0]);
    evaluator.add_inplace(eData[10], temp4_6_7);
    evaluator.add_inplace(eData[10], temp[8]);
    evaluator.add_inplace(eData[10], temp12_13_15);
    evaluator.add_inplace(eData[10], temp16_17_18_19);
    evaluator.add_inplace(eData[10], temp[21]);
    evaluator.add_inplace(eData[10], temp26_27);
    evaluator.add_inplace(eData[10], temp28_30_31);
    evaluator.add_inplace(eData[10], temp[33]);
    //evaluator.relinearize_inplace(eData[10], relin_keys);

    eData[11] = temp[1];
    evaluator.add_inplace(eData[11], temp5_6_7);
    evaluator.add_inplace(eData[11], temp[9]);
    evaluator.add_inplace(eData[11], temp[12]);
    evaluator.add_inplace(eData[11], temp16_18_19);
    evaluator.add_inplace(eData[11], temp20_22_23);
    evaluator.add_inplace(eData[11], temp24_26_27);
    evaluator.add_inplace(eData[11], temp30_31);
    evaluator.add_inplace(eData[11], temp32_35);
    //evaluator.relinearize_inplace(eData[11], relin_keys);

    eData[12] = temp1_2_3;
    evaluator.add_inplace(eData[12], temp5_6);
    evaluator.add_inplace(eData[12], temp9_11);
    evaluator.add_inplace(eData[12], temp13_14_15);
    evaluator.add_inplace(eData[12], temp16_17_18_19);
    evaluator.add_inplace(eData[12], temp21_23);
    evaluator.add_inplace(eData[12], temp[25]);
    evaluator.add_inplace(eData[12], temp28_29_31);
    evaluator.add_inplace(eData[12], temp32_34);
    //evaluator.relinearize_inplace(eData[12], relin_keys);

    evaluator.add_inplace(eData[13], temp0_3);
    evaluator.add_inplace(eData[13], temp[7]);
    evaluator.add_inplace(eData[13], temp9_10_11);
    evaluator.add_inplace(eData[13], temp[15]);
    evaluator.add_inplace(eData[13], temp16_18_19);
    evaluator.add_inplace(eData[13], temp20_21_22);
    evaluator.add_inplace(eData[13], temp[24]);
    evaluator.add_inplace(eData[13], temp29_30_31);
    evaluator.add_inplace(eData[13], temp33_34);
    //evaluator.relinearize_inplace(eData[13], relin_keys);

    eData[14] = temp[2];
    evaluator.add_inplace(eData[14], temp[4]);
    evaluator.add_inplace(eData[14], temp8_9_10);
    evaluator.add_inplace(eData[14], temp12_15);
    evaluator.add_inplace(eData[14], temp[19]);
    evaluator.add_inplace(eData[14], temp21_22_23);
    evaluator.add_inplace(eData[14], temp25_26_27);
    evaluator.add_inplace(eData[14], temp29_30);
    evaluator.add_inplace(eData[14], temp33_34_35);
    //evaluator.relinearize_inplace(eData[14], relin_keys);

    eData[15] = temp[1];
    evaluator.add_inplace(eData[15], temp4_5_6);
    evaluator.add_inplace(eData[15], temp8_9);
    evaluator.add_inplace(eData[15], temp12_14);
    evaluator.add_inplace(eData[15], temp16_17_18_19);
    evaluator.add_inplace(eData[15], temp20_21_22);
    evaluator.add_inplace(eData[15], temp24_26);
    evaluator.add_inplace(eData[15], temp28_31);
    evaluator.add_inplace(eData[15], temp32_34_35);
    //evaluator.relinearize_inplace(eData[15], relin_keys);

    evaluator.add_inplace(eData[16], temp0_1_3);
    evaluator.add_inplace(eData[16], temp[6]);
    evaluator.add_inplace(eData[16], temp[10]);
    evaluator.add_inplace(eData[16], temp12_13_14);
    evaluator.add_inplace(eData[16], temp18_19);
    evaluator.add_inplace(eData[16], temp21_22_23);
    evaluator.add_inplace(eData[16], temp24_25_27);
    evaluator.add_inplace(eData[16], temp32_33_34);
    //evaluator.relinearize_inplace(eData[16], relin_keys);

    eData[17] = temp0_1_2;
    evaluator.add_inplace(eData[17], temp5_7);
    evaluator.add_inplace(eData[17], temp[11]);
    evaluator.add_inplace(eData[17], temp12_13_15);
    evaluator.add_inplace(eData[17], temp[18]);
    evaluator.add_inplace(eData[17], temp[22]);
    evaluator.add_inplace(eData[17], temp24_25_26);
    evaluator.add_inplace(eData[17], temp28_29_30);
    evaluator.add_inplace(eData[17], temp32_33);
    //evaluator.relinearize_inplace(eData[17], relin_keys);

    eData[18] = temp1_2;
    evaluator.add_inplace(eData[18], temp4_7);
    evaluator.add_inplace(eData[18], temp8_9_11);
    evaluator.add_inplace(eData[18], temp12_15);
    evaluator.add_inplace(eData[18], temp17_19);
    evaluator.add_inplace(eData[18], temp20_21_22_23);
    evaluator.add_inplace(eData[18], temp24_25_27);
    evaluator.add_inplace(eData[18], temp29_31);
    evaluator.add_inplace(eData[18], temp34_35);
    //evaluator.relinearize_inplace(eData[18], relin_keys);

    evaluator.add_inplace(eData[19], temp0_1_3);
    evaluator.add_inplace(eData[19], temp4_6);
    evaluator.add_inplace(eData[19], temp[9]);
    evaluator.add_inplace(eData[19], temp13_15);
    evaluator.add_inplace(eData[19], temp16_17);
    evaluator.add_inplace(eData[19], temp21_22);
    evaluator.add_inplace(eData[19], temp24_25_26_27);
    evaluator.add_inplace(eData[19], temp28_30);
    evaluator.add_inplace(eData[19], temp[35]);
    //evaluator.relinearize_inplace(eData[19], relin_keys);

    eData[20] = temp0_3;
    evaluator.add_inplace(eData[20], temp4_5);
    evaluator.add_inplace(eData[20], temp8_10);
    evaluator.add_inplace(eData[20], temp14_15);
    evaluator.add_inplace(eData[20], temp16_18);
    evaluator.add_inplace(eData[20], temp[21]);
    evaluator.add_inplace(eData[20], temp25_27);
    evaluator.add_inplace(eData[20], temp28_29_31);
    evaluator.add_inplace(eData[20], temp32_33_35);
    //evaluator.relinearize_inplace(eData[20], relin_keys);

    eData[21] = temp1_2;
    evaluator.add_inplace(eData[21], temp4_5_7);
    evaluator.add_inplace(eData[21], temp10_11);
    evaluator.add_inplace(eData[21], temp12_14_15);
    evaluator.add_inplace(eData[21], temp[18]);
    evaluator.add_inplace(eData[21], temp20_22_23);
    evaluator.add_inplace(eData[21], temp24_25_26_27);
    evaluator.add_inplace(eData[21], temp28_30);
    evaluator.add_inplace(eData[21], temp32_34);
    //evaluator.relinearize_inplace(eData[21], relin_keys);

    evaluator.add_inplace(eData[22], temp2_3);
    evaluator.add_inplace(eData[22], temp4_6_7);
    evaluator.add_inplace(eData[22], temp[9]);
    evaluator.add_inplace(eData[22], temp[12]);
    evaluator.add_inplace(eData[22], temp16_18_19);
    evaluator.add_inplace(eData[22], temp[20]);
    evaluator.add_inplace(eData[22], temp24_25_27);
    evaluator.add_inplace(eData[22], temp28_29_30_31);
    evaluator.add_inplace(eData[22], temp[33]);
    //evaluator.relinearize_inplace(eData[22], relin_keys);

    eData[23] = temp0_2_3;
    evaluator.add_inplace(eData[23], temp6_7);
    evaluator.add_inplace(eData[23], temp8_11);
    evaluator.add_inplace(eData[23], temp[13]);
    evaluator.add_inplace(eData[23], temp17_18_19);
    evaluator.add_inplace(eData[23], temp[21]);
    evaluator.add_inplace(eData[23], temp[24]);
    evaluator.add_inplace(eData[23], temp28_30_31);
    evaluator.add_inplace(eData[23], temp32_34_35);
    //evaluator.relinearize_inplace(eData[23], relin_keys);

    eData[24] = temp[1];
    evaluator.add_inplace(eData[24], temp4_5_7);
    evaluator.add_inplace(eData[24], temp8_10);
    evaluator.add_inplace(eData[24], temp13_14_15);
    evaluator.add_inplace(eData[24], temp17_18);
    evaluator.add_inplace(eData[24], temp21_23);
    evaluator.add_inplace(eData[24], temp25_26_27);
    evaluator.add_inplace(eData[24], temp28_29_30_31);
    evaluator.add_inplace(eData[24], temp33_35);
    //evaluator.relinearize_inplace(eData[24], relin_keys);

    evaluator.add_inplace(eData[25], temp[0]);
    evaluator.add_inplace(eData[25], temp5_6_7);
    evaluator.add_inplace(eData[25], temp9_10);
    evaluator.add_inplace(eData[25], temp12_15);
    evaluator.add_inplace(eData[25], temp[19]);
    evaluator.add_inplace(eData[25], temp21_22_23);
    evaluator.add_inplace(eData[25], temp[27]);
    evaluator.add_inplace(eData[25], temp28_30_31);
    evaluator.add_inplace(eData[25], temp32_33_34);
    //evaluator.relinearize_inplace(eData[25], relin_keys);

    eData[26] = temp1_2_3;
    evaluator.add_inplace(eData[26], temp5_6);
    evaluator.add_inplace(eData[26], temp9_10_11);
    evaluator.add_inplace(eData[26], temp[14]);
    evaluator.add_inplace(eData[26], temp[16]);
    evaluator.add_inplace(eData[26], temp20_21_22);
    evaluator.add_inplace(eData[26], temp24_27);
    evaluator.add_inplace(eData[26], temp[31]);
    evaluator.add_inplace(eData[26], temp33_34_35);
    //evaluator.relinearize_inplace(eData[26], relin_keys);

    eData[27] = temp0_2;
    evaluator.add_inplace(eData[27], temp4_7);
    evaluator.add_inplace(eData[27], temp8_10_11);
    evaluator.add_inplace(eData[27], temp[13]);
    evaluator.add_inplace(eData[27], temp16_17_18);
    evaluator.add_inplace(eData[27], temp20_21);
    evaluator.add_inplace(eData[27], temp24_26);
    evaluator.add_inplace(eData[27], temp28_29_30_31);
    evaluator.add_inplace(eData[27], temp32_33_34);
    //evaluator.relinearize_inplace(eData[27], relin_keys);

    evaluator.add_inplace(eData[28], temp0_1_3);
    evaluator.add_inplace(eData[28], temp8_9_10);
    evaluator.add_inplace(eData[28], temp12_13_15);
    evaluator.add_inplace(eData[28], temp[18]);
    evaluator.add_inplace(eData[28], temp[22]);
    evaluator.add_inplace(eData[28], temp24_25_26);
    evaluator.add_inplace(eData[28], temp30_31);
    evaluator.add_inplace(eData[28], temp33_34_35);
    //evaluator.relinearize_inplace(eData[28], relin_keys);

    eData[29] = temp0_1_2;
    evaluator.add_inplace(eData[29], temp4_5_6);
    evaluator.add_inplace(eData[29], temp8_9);
    evaluator.add_inplace(eData[29], temp12_13_14);
    evaluator.add_inplace(eData[29], temp17_19);
    evaluator.add_inplace(eData[29], temp[23]);
    evaluator.add_inplace(eData[29], temp24_25_27);
    evaluator.add_inplace(eData[29], temp[30]);
    evaluator.add_inplace(eData[29], temp[34]);
    //evaluator.relinearize_inplace(eData[29], relin_keys);

    eData[30] = temp0_1_3;
    evaluator.add_inplace(eData[30], temp5_7);
    evaluator.add_inplace(eData[30], temp10_11);
    evaluator.add_inplace(eData[30], temp13_14);
    evaluator.add_inplace(eData[30], temp16_19);
    evaluator.add_inplace(eData[30], temp20_21_23);
    evaluator.add_inplace(eData[30], temp24_27);
    evaluator.add_inplace(eData[30], temp29_31);
    evaluator.add_inplace(eData[30], temp32_33_34_35);
    //evaluator.relinearize_inplace(eData[30], relin_keys);

    evaluator.add_inplace(eData[31], temp0_1_2_3);
    evaluator.add_inplace(eData[31], temp4_6);
    evaluator.add_inplace(eData[31], temp[11]);
    evaluator.add_inplace(eData[31], temp12_13_15);
    evaluator.add_inplace(eData[31], temp16_18);
    evaluator.add_inplace(eData[31], temp[21]);
    evaluator.add_inplace(eData[31], temp25_27);
    evaluator.add_inplace(eData[31], temp28_29);
    evaluator.add_inplace(eData[31], temp33_34);
    //evaluator.relinearize_inplace(eData[31], relin_keys);

    eData[32] = temp1_3;
    evaluator.add_inplace(eData[32], temp4_5_7);
    evaluator.add_inplace(eData[32], temp8_9_11);
    evaluator.add_inplace(eData[32], temp12_15);
    evaluator.add_inplace(eData[32], temp16_17);
    evaluator.add_inplace(eData[32], temp20_22);
    evaluator.add_inplace(eData[32], temp26_27);
    evaluator.add_inplace(eData[32], temp28_30);
    evaluator.add_inplace(eData[32], temp[33]);
    //evaluator.relinearize_inplace(eData[32], relin_keys);

    eData[33] = temp0_1_2_3;
    evaluator.add_inplace(eData[33], temp4_6);
    evaluator.add_inplace(eData[33], temp8_10);
    evaluator.add_inplace(eData[33], temp13_14);
    evaluator.add_inplace(eData[33], temp16_17_19);
    evaluator.add_inplace(eData[33], temp22_23);
    evaluator.add_inplace(eData[33], temp24_26_27);
    evaluator.add_inplace(eData[33], temp[30]);
    evaluator.add_inplace(eData[33], temp32_34_35);
    //evaluator.relinearize_inplace(eData[33], relin_keys);

    evaluator.add_inplace(eData[34], temp0_1_3);
    evaluator.add_inplace(eData[34], temp4_5_6_7);
    evaluator.add_inplace(eData[34], temp[9]);
    evaluator.add_inplace(eData[34], temp14_15);
    evaluator.add_inplace(eData[34], temp16_18_19);
    evaluator.add_inplace(eData[34], temp[21]);
    evaluator.add_inplace(eData[34], temp[24]);
    evaluator.add_inplace(eData[34], temp28_30_31);
    evaluator.add_inplace(eData[34], temp[32]);
    //evaluator.relinearize_inplace(eData[34], relin_keys);

    eData[35] = temp[0];
    evaluator.add_inplace(eData[35], temp4_6_7);
    evaluator.add_inplace(eData[35], temp8_10_11);
    evaluator.add_inplace(eData[35], temp12_14_15);
    evaluator.add_inplace(eData[35], temp18_19);
    evaluator.add_inplace(eData[35], temp20_23);
    evaluator.add_inplace(eData[35], temp[25]);
    evaluator.add_inplace(eData[35], temp29_30_31);
    evaluator.add_inplace(eData[35], temp[33]);
    //evaluator.relinearize_inplace(eData[35], relin_keys);
    }
}
// Compute the constants for Sbox,(x0,x1,x2)——>(x0,x0*x2+x1,-x0*x1+x0*x2+x2)
void HE_Sbox(vector<Ciphertext> &eData, Evaluator &evaluator,RelinKeys &relin_keys)
{
    // (x0,x1,x2)——> (x0,x0*x2+x1,-x0*x1+x0*x2+x2)
    Ciphertext c01 = eData[1];
    Ciphertext c02 = eData[2];
    //c01.multiplyBy(eData[0]);
    evaluator.multiply_inplace(c01, eData[0]);
    //evaluator.relinearize_inplace(c01, relin_keys);
    // c01*=eData[0];
    //c02.multiplyBy(eData[0]);
    evaluator.multiply_inplace(c02, eData[0]);
    //evaluator.relinearize_inplace(c02, relin_keys);
    // c02*=eData[0];

    //eData[1] += c02;
    //eData[2] -= c01;
    //eData[2] += c02;
    evaluator.add_inplace(eData[1], c02);
    //evaluator.relinearize_inplace(eData[1], relin_keys);
    evaluator.sub_inplace(eData[2], c01);
    //evaluator.relinearize_inplace(eData[2], relin_keys);
    evaluator.add_inplace(eData[2], c02);
    //evaluator.relinearize_inplace(eData[2], relin_keys);
    // omp_set_num_threads(12); // 设置线程数为12
    // #pragma omp parallel for
    for (uint64_t j = 3; j < BlockByte; j += 3)
    {
        c01 = eData[j + 1];
        //c01.multiplyBy(eData[j]);
        evaluator.multiply_inplace(c01, eData[j]);
        // std::cout << "relin_keys.size(): " << relin_keys.size() << std::endl;
        // std::cout << "c01.size(): " << c01.size() << std::endl;
        if(relin_keys.size()<c01.size()-2){
            
            std::cerr << "Error: Not enough relinearization keys." << std::endl;
            return;
        }
        //evaluator.relinearize_inplace(c01, relin_keys);

        c02 = eData[j + 2];
        //c02.multiplyBy(eData[j]);
        evaluator.multiply_inplace(c02, eData[j]);
        //evaluator.relinearize_inplace(c02, relin_keys);
        //eData[j + 1] += c02;
        evaluator.add_inplace(eData[j + 1], c02);
        //evaluator.relinearize_inplace(eData[j + 1], relin_keys);

        //eData[j + 2] -= c01;
        evaluator.sub_inplace(eData[j + 2], c01);
        //evaluator.relinearize_inplace(eData[j + 2], relin_keys);
        //eData[j + 2] += c02;
        evaluator.add_inplace(eData[j + 2], c02);
        //evaluator.relinearize_inplace(eData[j + 2], relin_keys);
    }
    // c01.cleanUp();
    // c02.cleanUp();
}

int main()
{
    std::cout << "Nr: " << Nr << std::endl;
    //=============客户端offline阶段================
    // 定义初始向量
    vector<vector<uint64_t>> IV(BlockByte,vector<uint64_t>(PlainBlock));
    for (unsigned i = 0; i < BlockByte; i++)
    {
        uint64_t t = (i+1)%PlainMod;
        for(unsigned j = 0; j < PlainBlock; j++)
        {
            IV[i][j] = t;
        }
    }
    std::cout<<"IV generated."<<std::endl;
    // Generating Public Key and encrypting the symmetric key
    auto start = std::chrono::high_resolution_clock::now();

    uint64_t p = Para_p;
    uint64_t m = Para_m;
    // uint64_t r = 1;
    // uint64_t bits = Para_bits;
    // uint64_t c = Para_c;
    // uint64_t d = 1; // slots = phi(m)/d = phi(m) = 32768 = PlainBlock
    // uint64_t k = 128;
    // uint64_t s = 1;

    // if (!m)
    //     m = FindM(k, bits, c, p, d, s, 0);

    // // Context context = ContextBuilder<BGV>()
    // //                                 .m(m)
    // //                                 .p(p)
    // //                                 .r(r)
    // //                                 .bits(bits)
    // //                                 .c(c)
    // //                                 .buildPtr();
    // shared_ptr<Context> context(ContextBuilder<BGV>()
    //                                 .m(m)
    //                                 .p(p)
    //                                 .r(r)
    //                                 .bits(bits)
    //                                 .c(c)
    //                                 .buildPtr());
    EncryptionParameters parms(scheme_type::bgv);
    size_t poly_modulus_degree = Para_m/2;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    //  vector<int> bit_sizes = {40, 40,40, 40 ,30, 30, 50,50,40,40}; // 你可以根据需要调整比特大小
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, bit_sizes));
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(PlainMod);
    SEALContext context(parms);
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

    size_t nslots = batch_encoder.slot_count();

    auto &plain_modulus = parms.plain_modulus();
    // 生成随机对称密钥
    random_device rd;
    vector<vector<uint64_t>> SymKey(BlockByte,vector<uint64_t>(PlainBlock));
    for (unsigned i = 0; i < BlockByte; i++)
    {
        uint64_t t = plain_modulus.reduce(rd());
        for (unsigned j = 0; j < PlainBlock; j++)
        {
            SymKey[i][j] = t;
        }
    }

    std::cout << "SymKey generated." << std::endl;
    //========
    // Generating symmetric key and key stream
    auto start_keyStream = std::chrono::high_resolution_clock::now();

    std::vector<uint64_t> NonceSet(PlainBlock);
    // std::vector<uint64_t> Xset(PlainByte * (Nr + 1));
    std::vector<uint64_t> RoundKeySet(PlainByte * (Nr + 1));
    std::vector<uint64_t> KeyStream(PlainByte);
    RandomBit<BlockByte * randbits> randomBit(Nr);
    std::cout << "Generating KeyStream..." << std::endl;
    for (uint64_t counter = counter_begin; counter <= counter_end; counter++)
    {
        uint64_t nonce = generate_secure_random_int(NonceSize);
        randomBit.generate_Instance_all_new(nonce, counter);
        auto &RanVecs = randomBit.roundconstants;
        // uint64_t nonce = counter;
        // std::vector<std::bitset<544>> RanVecs(Nr + 1);

        NonceSet[counter - counter_begin] = nonce;
        // 使用 std::array 代替 vector 并固定大小
        std::vector<uint64_t> state(BlockByte); // 初始化 state

        // 逐轮进行加密
        for (unsigned r = 0; r <= Nr; r++)
        {
            std::vector<uint64_t> X(BlockByte);
            std::vector<uint64_t> RoundKey(BlockByte);
            uint64_t temp;
            // 计算 X 和 RoundKey
            for (unsigned i = 0; i < BlockByte; ++i)
            {
                bool bit_array[randbits];
                for (unsigned j = 0; j < randbits; ++j)
                {
                    bit_array[j] = RanVecs[r][i * randbits + j];
                }
                BinStrToHex(bit_array, temp, randbits);
                // 强制转换为 uint64_t 类型
                X[i] = static_cast<uint64_t>(temp);
                if (Rkflag)
                {
                    RoundKey[i] = (SymKey[i][0] * X[i]) % PlainMod;
                }
                else
                {
                    RoundKey[i] = (SymKey[i][0] + X[i]) % PlainMod;
                }
            }
            // 将 X 和 RoundKey 复制到 Xset 和 RoundKeySet
            // memcpy(&Xset[PlainByte * r + BlockByte * (counter - counter_begin)], X.data(), BlockByte * sizeof(uint64_t));
            memcpy(&RoundKeySet[PlainByte * r + BlockByte * (counter - counter_begin)], RoundKey.data(), BlockByte * sizeof(uint64_t));
            if (r == 0)
            { // 初始轮
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (RoundKey[i] + IV[i][0]) % PlainMod;
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
            {                       // 最后一轮
                yusP.M36_5(state);  // 线性变换
                yusP.Sbox_5(state); // S盒
                yusP.M36_5(state);  // 线性变换
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
                memcpy(&KeyStream[(counter - counter_begin) * BlockByte], state.data(), BlockByte * sizeof(uint64_t));
            }
        }
    }
    auto end_keyStream = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds_keyStream = end_keyStream - start_keyStream;
    std::cout << "KeyStream Generation time: " << elapsed_seconds_keyStream.count() << "s\n";


    auto start_keyEncryption = std::chrono::high_resolution_clock::now();
    vector<Ciphertext> encryptedSymKey(BlockByte);
    vector<Ciphertext> encryptedSymKey01(len3);
    vector<Ciphertext> encryptedSymKey02(len3);
    encryptSymKey(encryptedSymKey, encryptedSymKey01, encryptedSymKey02, SymKey, batch_encoder, encryptor);
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
    // 离线客户端时间=KeyStream Generation time+PublicKey generation and SymKey FHE time
    double total_time_off = elapsed_seconds_keyStream.count() + elapsed_seconds_PubKey.count() + elapsed_seconds_PubKey.count() + keyEncryption;
    std::cout << "Encryption offline total time: " << total_time_off << "s\n";
    if (!writeEncryptedSymKey(encryptedSymKey, "Client_encryptedSymKey.bin"))
    {
        return false;
    }
    std::cout << "Client_encryptedSymKey.bin has been written to file." << std::endl;

    //=============服务端offline阶段================
    std::cout << "Generating XOF stream..." << std::endl;
    std::vector<std::vector<uint64_t>> Xset01(len3, std::vector<uint64_t>(PlainBlock));
    std::vector<std::vector<uint64_t>> Xset02(len3, std::vector<uint64_t>(PlainBlock));
    std::vector<std::vector<uint64_t>> Xset(NrBlockByte, std::vector<uint64_t>(PlainBlock));
    // helib::Ptxt<helib::BGV> ptxt(*context);
    // vector<helib::Ptxt<helib::BGV>> Xset(NrBlockByte,ptxt);
    // vector<helib::Ptxt<helib::BGV>> Xset01(len3,ptxt);
    // vector<helib::Ptxt<helib::BGV>> Xset02(len3,ptxt);
    auto start_XOF = std::chrono::high_resolution_clock::now();
    for (uint64_t counter = counter_begin; counter <= counter_end; counter++)
    {
        uint64_t nonce = NonceSet[counter - counter_begin];
        randomBit.generate_Instance_all_new(nonce, counter);
        auto &RanVecs = randomBit.roundconstants;
        // 逐轮进行加密
        for (unsigned r = 0; r <= Nr; r++)
        {
            std::vector<uint64_t> X(BlockByte);
            uint64_t temp;
            // 计算 X
            for (unsigned i = 0; i < BlockByte; ++i)
            {
                bool bit_array[randbits];
                for (unsigned j = 0; j < randbits; ++j)
                {
                    bit_array[j] = RanVecs[r][i * randbits + j];
                }
                BinStrToHex(bit_array, temp, randbits);
                // 强制转换为 uint64_t 类型
                //X[i] = static_cast<uint64_t>(temp);
                X[i] = temp;
                ;
            }
            // 将 X 复制到 Xset
            for (int i = 0; i < BlockByte; i++)
            {
                Xset[r * BlockByte + i][counter - counter_begin] = X[i];
            }
            if (r == 0)
            { // 初始轮
                for (int i = 0; i < BlockByte; i += 3)
                {
                    int index = i / 3;
                    Xset01[index][counter - counter_begin] = (X[i] * X[i + 1]) % PlainMod;
                    Xset02[index][counter - counter_begin] = (X[i] * X[i + 2]) % PlainMod;
                }
            }
        }
    }
    auto end_XOF = std::chrono::high_resolution_clock::now();
    double XOF_time = std::chrono::duration<double>(end_XOF - start_XOF).count();
    std::cout << "XOF stream Generation time: " << XOF_time << "s\n";

    vector<Plaintext> encodedIV(BlockByte);
    vector<Plaintext> encoded_Iv0Iv1subIv2(len3);
    vector<Plaintext> encoded_Iv0Iv2(len3);
    auto start_IV = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < len3; i ++)
    {
        batch_encoder.encode(IV[i], encodedIV[i]);
        batch_encoder.encode(IV[i + 1], encodedIV[i + 1]);
        batch_encoder.encode(IV[i + 2], encodedIV[i + 2]);
        vector<uint64_t> IV01sub2(PlainBlock,(IV[i][0]*IV[i+1][0]-IV[i+2][0])%PlainMod);
        batch_encoder.encode(IV01sub2,encoded_Iv0Iv1subIv2[i]);
        vector<uint64_t> IV02(PlainBlock,(IV[i][0]*IV[i+2][0])%PlainMod);
        batch_encoder.encode(IV02,encoded_Iv0Iv2[i]);
    }
    auto end_IV = std::chrono::high_resolution_clock::now();
    std::cout << "encodeIV time: " << std::chrono::duration<double>(end_IV - start_IV).count() << "s\n";

    vector<Plaintext> encodedXset(NrBlockByte);
    vector<Plaintext> encodedXset01(len3);
    vector<Plaintext> encodedXset02(len3);

    auto start_Xset = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < NrBlockByte; i++)
    {
        batch_encoder.encode(Xset[i], encodedXset[i]);
    }
    for (int i = 0; i < len3; i++)
    {
        batch_encoder.encode(Xset01[i], encodedXset01[i]);
        batch_encoder.encode(Xset02[i], encodedXset02[i]);
    }
    auto end_Xset = std::chrono::high_resolution_clock::now();
    double Encode_time = std::chrono::duration<double>(end_Xset - start_Xset).count();
    std::cout << "encode time: " << Encode_time << "s\n";
    std::cout << "encode time/slots:" << Encode_time/(NrBlockByte+2*len3) * pow(10,6) << "us\n";

    int noise_budget = min_noise_budget(encryptedSymKey,decryptor);
    std::cout << "noise budget initially: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    // 计算 encryptedRoundKeySet
    uint64_t eRk_len = BlockByte * (Nr + 1);
    vector<Ciphertext> encryptedRoundKeySet(eRk_len);
    for (int i = 0; i < eRk_len; i++)
    {
        encryptedRoundKeySet[i] = encryptedSymKey[i % BlockByte];
    }
    auto start_RoundKeySet_FHE = std::chrono::high_resolution_clock::now();
    if (Rkflag)
    {
        for (int i = 0; i < eRk_len; i++)
        {
            evaluator.multiply_plain_inplace(encryptedRoundKeySet[i], encodedXset[i]);
        }
    }
    else
    {
        for (int i = 0; i < eRk_len; i++)
        {
            evaluator.add_plain_inplace(encryptedRoundKeySet[i], encodedXset[i]);
        }
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
    //     if (!verifyDecryption18(encryptedRoundKeySet, RoundKeySet, secretKey, ea))
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
    auto start_roundkey = std::chrono::high_resolution_clock::now();
    auto end_roundkey = std::chrono::high_resolution_clock::now();

    vector<Ciphertext> encryptedKeyStream(encryptedRoundKeySet.begin(), encryptedRoundKeySet.begin() + BlockByte);
    std::cout << "whiteround + sbox start" << std::endl;
    Ciphertext K0Iv1R0;
    Ciphertext K1Iv0R1;
    Ciphertext K0Iv2R0;
    Ciphertext K2Iv0R2;
    vector<Ciphertext> tmp01 = encryptedSymKey01;
    vector<Ciphertext> tmp02 = encryptedSymKey02;
    Plaintext tmp;
    start_sbox = std::chrono::high_resolution_clock::now();
    for (uint64_t i = 0; i < BlockByte; i += 3) // BlockByte
    {
        int index = i / 3;
        K1Iv0R1 = encryptedSymKey[i + 1];
        //tmp = encodedIV[i] + encodedXset[i + 1];
        evaluator.multiply_plain_inplace(K1Iv0R1, multiplyAndMod(encodedXset[i + 1], IV[i][0]));
        K0Iv1R0 = encryptedSymKey[i + 0];
        evaluator.multiply_plain_inplace(K0Iv1R0, multiplyAndMod(encodedXset[i], IV[i + 1][0]));
        K2Iv0R2 = encryptedSymKey[i + 2];
        evaluator.multiply_plain_inplace(K2Iv0R2, multiplyAndMod(encodedXset[i + 2], IV[i][0]));
        K0Iv2R0 = encryptedSymKey[i + 0];
        evaluator.multiply_plain_inplace(K0Iv2R0, multiplyAndMod(encodedXset[i], IV[i + 2][0]));
        // 计算Sbox 0
        evaluator.add_plain_inplace(tmp01[index], encodedIV[i]);
        // 计算Sbox 1
        evaluator.multiply_plain_inplace(tmp02[index], encodedXset02[index]);
        evaluator.add_inplace(tmp02[index], K2Iv0R2);
        evaluator.add_inplace(tmp02[index], K0Iv2R0);
        evaluator.add_plain_inplace(tmp02[index], encoded_Iv0Iv2[index]);
        evaluator.add_plain_inplace(encryptedKeyStream[i + 1], encodedIV[i + 1]);
        evaluator.add_inplace(encryptedKeyStream[i + 1], tmp02[index]);
        // 计算Sbox 2
        evaluator.multiply_plain_inplace(tmp01[index], encodedXset01[index]);
        evaluator.add_inplace(tmp01[index], K1Iv0R1);
        evaluator.add_inplace(tmp01[index], K0Iv1R0);
        evaluator.add_plain_inplace(tmp01[index], encoded_Iv0Iv1subIv2[index]);
        evaluator.add_inplace(encryptedKeyStream[i + 2], tmp02[index]);
        evaluator.sub_inplace(encryptedKeyStream[i + 2], tmp01[index]);
    }
    end_sbox = std::chrono::high_resolution_clock::now();
    Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
    // 输出 Sbox_time
    std::cout << "SboxAndWhiteround time: " << Sbox_time << "s\n";
    noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
    std::cout << "noise budget after SboxAndWhiteround: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    // 明文密钥流
    vector<uint64_t> KeyStream2(PlainByte);
    if (deflag)
    {
        // 对IV和RoundKeySet进行异或
        for (uint64_t i = 0; i < PlainByte; i++)
        {
            KeyStream2[i] = (IV[i % BlockByte][0] + RoundKeySet[i]) % PlainMod;
        }
        // sbox
        yusP.Sbox_5(KeyStream2);
        // 使用 verifyDecryption 函数解密并验证 KeyStream
        if (!verifyDecryption36(encryptedKeyStream, KeyStream2, batch_encoder, decryptor))
        {
            std::cerr << "Decryption verification failed for KeyStream." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for whiteround." << std::endl;
    }

    for (uint64_t r = 1; r < Nr; r++)
    {
        std::cout << "Round " << r << " start" << std::endl;
        if (r > 1)
        {
            start_sbox = std::chrono::high_resolution_clock::now();
            // S Layer
            for(int j=0;j<BlockByte;j++){
                evaluator.relinearize_inplace(encryptedKeyStream[j],relin_keys);
            }
            HE_Sbox(encryptedKeyStream,evaluator,relin_keys);
            end_sbox = std::chrono::high_resolution_clock::now();
            Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
            // noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
            // std::cout << "noise budget after sbox: " << noise_budget << std::endl;
            // if (noise_budget <= 0)
            // {
            //     std::cerr << "noise budget is not enough!!!" << std::endl;
            //     return 0;
            // }
            if (deflag)
            {
                yusP.Sbox_5(KeyStream2);
                if (!verifyDecryption36(encryptedKeyStream, KeyStream2,batch_encoder, decryptor))
                {
                    std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
                    return 0;
                }
                std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
            }
        }
        start_linear = std::chrono::high_resolution_clock::now();
        // Linear Layer
        HE_M2(encryptedKeyStream,evaluator,relin_keys);
        end_linear = std::chrono::high_resolution_clock::now();
        Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
        // noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
        // std::cout << "noise budget after linear: " << noise_budget << std::endl;
        // if (noise_budget <= 0)
        // {
        //     std::cerr << "noise budget is not enough!!!" << std::endl;
        //     return 0;
        // }
        if (deflag)
        {
            for (int i = 0; i < PlainBlock; i++)
            {
                vector<uint64_t> tmp(BlockByte);
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
            if (!verifyDecryption36(encryptedKeyStream, KeyStream2, batch_encoder, decryptor))
            {
                std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
        }
        // Round Key Addition
        start_roundkey = std::chrono::high_resolution_clock::now();
        // omp_set_num_threads(12); // 设置线程数为12
        // #pragma omp parallel for
        for (uint64_t j = 0; j < BlockByte; j++)
        {
            //encryptedKeyStream[j] += encryptedRoundKeySet[r * BlockByte + j];
            evaluator.add_inplace(encryptedKeyStream[j], encryptedRoundKeySet[r * BlockByte + j]);
        }
        end_roundkey = std::chrono::high_resolution_clock::now();
        Add_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
        // noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
        // std::cout << "noise budget after Add: " << noise_budget << std::endl;
        // if (noise_budget <= 0)
        // {
        //     std::cerr << "noise budget is not enough!!!" << std::endl;
        //     return 0;
        // }
        if (deflag)
        {
            for (uint64_t i = 0; i < PlainByte; i++)
            {
                KeyStream2[i] = (KeyStream2[i] + RoundKeySet[r * PlainByte + i]) % PlainMod;
            }
            if (!verifyDecryption36(encryptedKeyStream, KeyStream2, batch_encoder, decryptor))
            {
                std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
                return 0;
            }
            std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;
        }
    }
// 最后一轮
#if (1)
    std::cout << "Round " << Nr << " start" << std::endl;
    start_linear = std::chrono::high_resolution_clock::now();
    // Linear Layer
    HE_M2(encryptedKeyStream,evaluator,relin_keys);
    end_linear = std::chrono::high_resolution_clock::now();
    Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
    // noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
    // std::cout << "noise budget after linear: " << noise_budget << std::endl;
    // if (noise_budget <= 0)
    // {
    //     std::cerr << "noise budget is not enough!!!" << std::endl;
    //     return 0;
    // }
    if (deflag)
    {
        for (int i = 0; i < PlainBlock; i++)
        {
            vector<uint64_t> tmp(BlockByte);
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
        if (!verifyDecryption36(encryptedKeyStream, KeyStream2, batch_encoder, decryptor))
        {
            std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
    }
    start_sbox = std::chrono::high_resolution_clock::now();
    // S Layer
    HE_Sbox(encryptedKeyStream,evaluator,relin_keys);
    end_sbox = std::chrono::high_resolution_clock::now();
    Sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
    // noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
    // std::cout << "noise budget after sbox: " << noise_budget << std::endl;
    // if (noise_budget <= 0)
    // {
    //     std::cerr << "noise budget is not enough!!!" << std::endl;
    //     return 0;
    // }
    if (deflag)
    {
        yusP.Sbox_5(KeyStream2);
        if (!verifyDecryption36(encryptedKeyStream, KeyStream2, batch_encoder, decryptor))
        {
            std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
    }
    start_linear = std::chrono::high_resolution_clock::now();
    // Linear Layer
    HE_M2(encryptedKeyStream,evaluator,relin_keys);
    end_linear = std::chrono::high_resolution_clock::now();
    Linear_time += std::chrono::duration<double>(end_linear - start_linear).count();
    // noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
    // std::cout << "noise budget after linear: " << noise_budget << std::endl;
    // if (noise_budget <= 0)
    // {
    //     std::cerr << "noise budget is not enough!!!" << std::endl;
    //     return 0;
    // }
    if (deflag)
    {
        for (int i = 0; i < PlainBlock; i++)
        {
            vector<uint64_t> tmp(BlockByte);
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
        if (!verifyDecryption36(encryptedKeyStream, KeyStream2, batch_encoder, decryptor))
        {
            std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
    }
    // add
    start_roundkey = std::chrono::high_resolution_clock::now();
    for (uint64_t j = 0; j < BlockByte; j++)
    {
        //encryptedKeyStream[j] += encryptedRoundKeySet[Nr * BlockByte + j];
        evaluator.add_inplace(encryptedKeyStream[j], encryptedRoundKeySet[Nr * BlockByte + j]);
    }
    end_roundkey = std::chrono::high_resolution_clock::now();
    Add_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
    noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
    // std::cout << "noise budget after Add: " << noise_budget << std::endl;
    // if (noise_budget <= 0)
    // {
    //     std::cerr << "noise budget is not enough!!!" << std::endl;
    //     return 0;
    // }
    if (deflag)
    {
        for (uint64_t i = 0; i < PlainByte; i++)
        {
            KeyStream2[i] = (KeyStream2[i] + RoundKeySet[Nr * PlainByte + i]) % PlainMod;
        }
        if (!verifyDecryption36(encryptedKeyStream, KeyStream2, batch_encoder, decryptor))
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
    std::cout << "Add time: " << Add_time << "s\n";
    std::cout << "Sbox time: " << Sbox_time << "s\n";
    std::cout << "Linear time: " << Linear_time << "s\n";
    // 计算总时间
    double total_time = XOF_time + Encode_time + RoundKey_time + Add_time + Sbox_time + Linear_time;
    std::cout << "Server offline total time: " << total_time << "s\n";
    // 计算吞吐量,KiB/min
    double throughput = (Plainbits * 60) / (pow(2, 13) * total_time);
    std::cout << "Throughput: " << throughput << "KiB/min\n";
    noise_budget = min_noise_budget(encryptedKeyStream,decryptor);
    std::cout << "noise budget finally: " << noise_budget << std::endl;
    if (noise_budget <= 0)
    {
        std::cerr << "noise budget is not enough!!!" << std::endl;
        return 0;
    }
    if (plainflag)
    {
        // for (int i = 0; i < encryptedKeyStream.size(); i++)
        // {
        //    encryptedKeyStream[i].bringToSet(encryptedKeyStream[i].naturalPrimeSet());
        // }
        if (!verifyDecryption36(encryptedKeyStream, KeyStream, batch_encoder, decryptor))
        {
            std::cerr << "Decryption verification failed for KeyStream." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream." << std::endl;
    }
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
            << std::left << std::setw(10) << p
            << std::left << std::setw(10) << nslots
            //<< std::left << std::setw(5) << bits
            //<< std::left << std::setw(4) << c
            //<< std::left << std::setw(6) << Qbits
            << std::fixed << std::setprecision(3)
            //<< std::left << std::setw(14) << SecurityLevel
            << std::left << std::setw(7) << XOF_time
            << std::left << std::setw(10) << Encode_time
            << std::left << std::setw(12) << RoundKey_time
            << std::left << std::setw(7) << Add_time
            << std::left << std::setw(8) << Sbox_time
            << std::left << std::setw(10) << Linear_time
            << std::left << std::setw(9) << total_time
            << std::left << std::setw(20) << throughput
            << std::left << std::setw(10) << noise_budget
            << std::endl;
    outfile.close();
    std::cout << "test_Yus_p_C36_ClientAndServer5_2.txt updated." << std::endl;
    return 0;
}