#ifndef FHETOOL_YUS_P_HPP
#define FHETOOL_YUS_P_HPP

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

#include "params_Yus_p.hpp"

using namespace std;
using namespace helib;
using namespace NTL;


void encodeTo32Ctxt(vector<ZZX> &encData, const vector<long> &data, const EncryptedArray &ea);

void decodeTo32Ctxt(vector<long> &data, const vector<vector<long>> &encData,
                    const EncryptedArray &ea);

bool verifyDecryption32(const std::vector<Ctxt> &encryptedVec, const vector<long> &originalVec, const SecKey &secretKey,
                        const EncryptedArray &ea);

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