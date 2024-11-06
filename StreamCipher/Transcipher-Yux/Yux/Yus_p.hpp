#ifndef YUS_P_HPP
#define YUS_P_HPP

#include <cstdint> // 如果需要使用 uint64_t 类型
#include <iostream>
#include <array>
#include <vector>
#include <cmath>

// 列混淆
void MC(std::vector<long> &A);
// 行移位
void MR(std::vector<long> &A);
// Sbox
void Sbox(std::vector<long> &A);

#endif // YUS_P_HPP