#ifndef TOOL_HPP
#define TOOL_HPP

#include <iostream>
#include <vector>
#include <random>
#include <climits>

// 函数声明
uint64_t generate_secure_random_int(unsigned m);
void HexToBinStr(uint64_t hex, bool *bin_str, unsigned n);
void BinStrToHex(const bool *bin_str, uint64_t &dec_hex, unsigned n);

#endif // TOOL_HPP