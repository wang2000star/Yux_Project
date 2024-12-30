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
