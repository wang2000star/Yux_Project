#ifndef YUS_P_HPP
#define YUS_P_HPP

#include <cstdint> // 如果需要使用 uint64_t 类型
#include <iostream>
#include <array>
#include <vector>
#include <cmath>

class YusP {
public:
    explicit YusP(long plainMod) : PlainMod(plainMod) {}

    void MC32(std::vector<long> &A);
    void MR32(std::vector<long> &A);
    void MC_MR(std::vector<long> &A);
    void MC32_2(std::vector<long> &A);
    void MR32_2(std::vector<long> &A);
    void MC64(std::vector<long> &A);
    void MR64(std::vector<long> &A);
    void Sbox(std::vector<long> &A);

private:
    long PlainMod;
};

#endif // YUS_P_HPP