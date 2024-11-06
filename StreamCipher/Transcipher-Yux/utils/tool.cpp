#include "tool.hpp"

// 生成m比特强度的随机整数
long generate_secure_random_int(unsigned m) {
    if (m > 64) {
        throw std::invalid_argument("m exceeds the maximum bit size of long");
    }

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<long> dis(0, (1ULL << m) - 1);

    return dis(gen);
}
