#ifndef TOOL_HPP
#define TOOL_HPP

#include <iostream>
#include <vector>
#include <random>
#include <climits>

// 函数声明
long generate_secure_random_int(unsigned m);

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