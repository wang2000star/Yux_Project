#ifndef FHETOOL_YUX2_8_HPP
#define FHETOOL_YUX2_8_HPP

#include <iostream>
#include <cstring>
#include <vector>
#include <stdint.h>
#include <string>

#include <NTL/ZZX.h>
#include <NTL/GF2X.h>
#include <helib/helib.h>
#include <helib/ArgMap.h>
#include <helib/DoubleCRT.h>

#include "params_Yux2_8.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

void encodeTo1Ctxt(NTL::Vec<NTL::ZZX>& encData, const NTL::Vec<uint8_t>& data,
                   const helib::EncryptedArrayDerived<helib::PA_GF2>& ea);

void decodeTo1Ctxt(NTL::Vec<uint8_t>& data, const NTL::Vec<NTL::ZZX>& encData,
                   const helib::EncryptedArrayDerived<helib::PA_GF2>& ea);
// 函数：解密并验证密文是否正确
bool verifyDecryption1(const vector<Ctxt> &encryptedVec, const SecKey &secretKey,
                      const EncryptedArrayDerived<PA_GF2> &ea, const Vec<uint8_t> &originalVec);

void encodeTo16Ctxt(vector<ZZX>& encData, const vector<uint8_t>& data,
		const EncryptedArrayDerived<PA_GF2>& ea);

// Decode native HE plaintext as Yux plaintext/ciphertext bytes
void decodeTo16Ctxt(vector<uint8_t>& data, const vector<ZZX>& encData,
		const EncryptedArrayDerived<PA_GF2>& ea);

bool verifyDecryption16(const vector<Ctxt> &encryptedVec, const SecKey &secretKey,
                      const EncryptedArrayDerived<PA_GF2> &ea, const vector<uint8_t> &originalVec);
                      
bool verifyDecryption_p16(const vector<Ctxt>& encryptedSymKey, const vector<long>& originalSymKey, const SecKey& secretKey, const EncryptedArray& ea);

void encodeTo16Ctxt_p(vector<vector<long>>& encData, const vector<long>& data, const EncryptedArray& ea);

void decodeTo16Ctxt_p(vector<long>& data, const vector<vector<long>>& encData, const EncryptedArray& ea);
// 通用写入函数模板
template <typename T>
bool writeToFile(const T *data, const std::string &filename, size_t length)
{
    std::ofstream out(filename);
    if (out.is_open())
    {
        for (size_t i = 0; i < length; ++i)
        {
            out << std::hex << std::setw(sizeof(T) * 2) << std::setfill('0') << static_cast<long>(data[i]) << " ";
        }
        out.close();
        return true;
    }
    else
    {
        std::cerr << "Failed to open " << filename << " for writing" << std::endl;
        return false;
    }
}

// 通用读取函数模板
template <typename T>
bool readFromFile(T *data, const std::string &filename, size_t length)
{
    std::ifstream in(filename);
    if (!in.is_open())
    {
        std::cerr << "Failed to open " << filename << " for reading" << std::endl;
        return false;
    }

    for (size_t i = 0; i < length; ++i)
    {
        long temp;
        in >> std::hex >> temp;
        data[i] = static_cast<T>(temp);
    }

    in.close();
    return true;
}

#endif //FHETOOL_HPP