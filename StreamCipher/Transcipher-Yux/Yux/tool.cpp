#include "tool.hpp"

// 生成m比特强度的随机整数
uint64_t generate_secure_random_int(unsigned m) {
    if (m > 64) {
        throw std::invalid_argument("m exceeds the maximum bit size of uint64_t");
    }

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(0, (1ULL << m) - 1);

    return dis(gen);
}
/*
16进制转2进制，8bits，
binstr[0]是最低位，binstr[7]是最高位，
有一个性质：00~7F的最高位是0，80~FF的最高位是1
所以sboxtable前128个元素的最高位是0，后128个元素的最高位是1
*/
void HexToBinStr(uint64_t hex, bool *bin_str, unsigned n)
{
    for (unsigned i = 0; i < n; ++i)
    {
        bin_str[i] = (hex >> i) & 1; // 位移操作替代 % 2
    }
}
// 比特串转整数
void BinStrToHex(const bool *bin_str, uint64_t &dec_hex, unsigned n)
{
    dec_hex = 0; // 确保 dec_hex 从零开始
    for (unsigned i = 0; i < n; ++i)
    {
        dec_hex = (dec_hex << 1) | bin_str[i]; // 左移一位并添加当前比特
    }
}