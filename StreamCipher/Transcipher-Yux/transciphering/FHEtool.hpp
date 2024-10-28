#ifndef FHETOOL_HPP
#define FHETOOL_HPP

#include <iostream>
#include <cstring>
#include <stdint.h>
#include <string>

#include <NTL/ZZX.h>
#include <NTL/GF2X.h>
#include <helib/helib.h>
#include <helib/ArgMap.h>
#include <helib/DoubleCRT.h>

#include "params.hpp"

void encodeTo1Ctxt(NTL::Vec<NTL::ZZX>& encData, const NTL::Vec<uint8_t>& data,
                   const helib::EncryptedArrayDerived<helib::PA_GF2>& ea);

void decodeTo1Ctxt(NTL::Vec<uint8_t>& data, const NTL::Vec<NTL::ZZX>& encData,
                   const helib::EncryptedArrayDerived<helib::PA_GF2>& ea);

// 通用写入函数模板
template <typename T>
bool writeToFile(const T *data, const std::string &filename, size_t length)
{
    std::ofstream out(filename);
    if (out.is_open())
    {
        for (size_t i = 0; i < length; ++i)
        {
            out << std::hex << std::setw(sizeof(T) * 2) << std::setfill('0') << static_cast<uint64_t>(data[i]) << " ";
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
        uint64_t temp;
        in >> std::hex >> temp;
        data[i] = static_cast<T>(temp);
    }

    in.close();
    return true;
}

#endif //FHETOOL_HPP