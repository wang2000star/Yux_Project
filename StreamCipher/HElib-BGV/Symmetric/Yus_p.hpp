#ifndef YUS_P_HPP
#define YUS_P_HPP

#include <cstdint> // 如果需要使用 uint64_t 类型
#include <iostream>
#include <array>
#include <vector>
#include <cmath>

class YusP
{
public:
    explicit YusP(long plainMod) : PlainMod(plainMod) {}
    // MC32和MR32是老师第一个版本,MC_MR是合成
    // MC32_2和MR32_2是老师第二个版本
    // MC64和MR64是我对老师第一个版本的扩展
    // 第一个版本和第二个版本的s盒都是Sbox,我对第一个版本拓展也是

    void MC32(std::vector<long> &data);
    void MR32(std::vector<long> &data);
    void MC_MR(std::vector<long> &data);
    void MC32_2(std::vector<long> &data);
    void MR32_2(std::vector<long> &data);
    void MC64(std::vector<long> &data);
    void MR64(std::vector<long> &data);
    void Sbox_feistel(std::vector<long> &data);
    void MC48_3(std::vector<long> &data);
    void MR48_3(std::vector<long> &data);
    void Sbox_3(std::vector<long> &data);
    void MC64_4(std::vector<long> &data);
    void MR64_4(std::vector<long> &data);
    void Sbox_4(std::vector<long> &data);
    void M36_5(std::vector<long> &data);
    void Sbox_5(std::vector<long> &data);
    void Sbox_5_last(std::vector<long> &data);
    void MC64_6(std::vector<long> &data);
    void MR64_6(std::vector<long> &data);
    void Sbox_6(std::vector<long> &data);
    void M36_7(std::vector<long> &data);
    void Sbox_7(std::vector<long> &data);
    void M18_8(std::vector<long> &data);
    void Sbox_8(std::vector<long> &data);

    void M36(std::vector<long> &data);
    void Sbox(std::vector<long> &data);

private:
    long PlainMod;
};

#endif // YUS_P_HPP