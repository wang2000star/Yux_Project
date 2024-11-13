#ifndef HERA_HPP
#define HERA_HPP

#include <cstdint> // 如果需要使用 uint64_t 类型
#include <iostream>
#include <array>
#include <vector>
#include <cmath>

class Hera {
public:
    explicit Hera(long plainMod) : PlainMod(plainMod) {}

    void MC(std::vector<long> &A);
    void MR(std::vector<long> &A);
    void Sbox(std::vector<long> &A);

private:
    long PlainMod;
};

#endif // YUS_P_HPP