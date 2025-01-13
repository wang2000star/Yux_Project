#ifndef TOOL_HPP
#define TOOL_HPP

#include <iostream>
#include <vector>
#include <random>
#include <climits>
#include <fstream>
#include <chrono>

#include <SEAL-4.1/seal/seal.h>

extern "C"
{
#include "../keccak/KeccakHash.h"
}

using namespace std;
using namespace seal;
// 函数声明
long generate_secure_random_int(unsigned m);

void random_init_shake(long nonce, long block_counter, Keccak_HashInstance &shake128_);

long generate_random_field_element(Keccak_HashInstance &shake128, bool allow_zero, long max_prime_size, long PlainMod);


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

int min_noise_budget(vector<Ciphertext> &eData ,Decryptor &decryptor);

bool writeEncryptedSymKey(const vector<Ciphertext> &encryptedSymKey, const std::string &filename);

void decodeToCtxt(std::vector<long> &data, const std::vector<std::vector<long>> &encData, const long CtxtWords,const long nslots);

bool verifyDecryption(const std::vector<Ciphertext> &encryptedVec, const std::vector<long> &originalVec, BatchEncoder &batch_encoder, Decryptor &decryptor,
                      const long CtxtWords, const long PlainBlock, const long nslots, const long pmod);

void encryptSymKey(vector<Ciphertext> &encryptedSymKey, const vector<long> &SymKey, BatchEncoder &batch_encoder, Encryptor &encryptor,const long nslots);

bool verify_encryptSymKey(vector<Ciphertext> &encryptedSymKey, const vector<long> &SymKey,BatchEncoder &batch_encoder, Decryptor &decryptor);

// 声明并定义 rdtsc 函数
inline uint64_t rdtsc() {
    uint32_t lo, hi;
    __asm__ __volatile__ (
        "rdtsc"
        : "=a" (lo), "=d" (hi)
    );
    return (static_cast<uint64_t>(hi) << 32) | lo;
}

template <typename T>
void HexToBinStr(T hex, bool *bin_str, unsigned n)
{
    for (unsigned i = 0; i < n; ++i)
    {
        bin_str[i] = (hex >> i) & 1; // 位移操作替代 % 2
    }
}

// 比特串转整数
template <typename T>
void BinStrToHex(const bool *bin_str, T &dec_hex, unsigned n)
{
    dec_hex = 0; // 确保 dec_hex 从零开始
    for (unsigned i = 0; i < n; ++i)
    {
        dec_hex |= (static_cast<T>(bin_str[i]) << i); // 左移一位并添加当前比特
    }
}

#endif // TOOL_HPP