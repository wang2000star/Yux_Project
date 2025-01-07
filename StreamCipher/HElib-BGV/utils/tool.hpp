#ifndef TOOL_HPP
#define TOOL_HPP

#include <iostream>
#include <vector>
#include <random>
#include <climits>
#include <helib/helib.h>


// 函数声明
long generate_secure_random_int(unsigned m);

int min_noise_budget(std::vector<helib::Ctxt> &eData);

helib::zzX multiplyAndMod(const helib::zzX &a, long b,long pmod);

bool writeEncryptedSymKey(const std::vector<helib::Ctxt> &encryptedSymKey, const std::string &filename);

void decodeToCtxt(std::vector<long> &data, const std::vector<NTL::vec_long> &encData, const long CtxtWords,const long PlainBlock,const long nslots);

bool verifyDecryption(const std::vector<helib::Ctxt> &encryptedVec, const std::vector<long> &originalVec, const helib::SecKey &secretKey,
                      const helib::Cmodulus &cmodulus, const long CtxtWords, const long PlainBlock, const long nslots, const long pmod);

void encryptSymKey(std::vector<helib::Ctxt> &encryptedSymKey, const std::vector<long> &SymKey, std::unique_ptr<helib::PubKey> &pk, 
                   const helib::Cmodulus &cmodulus,const long nslots);

bool verify_encryptSymKey(std::vector<helib::Ctxt> &encryptedSymKey, const std::vector<long> &SymKey, const helib::SecKey &secretKey,
                          const helib::Cmodulus &cmodulus);
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