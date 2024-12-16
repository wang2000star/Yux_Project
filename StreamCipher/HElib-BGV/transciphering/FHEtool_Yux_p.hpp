#ifndef FHETOOL_YUX_P_HPP
#define FHETOOL_YUX_P_HPP

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

#include "params_Yux_p.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

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