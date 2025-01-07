// class Pasta2 {
#ifndef PASTA2_3_HPP
#define PASTA2_3_HPP

#pragma once

#include <array>
#include <string>
#include <vector>
#include "KeccakHash.h"
// #
// extern "C" {
// #include "KeccakHash.h"
// }
using namespace std;

constexpr long KeyWords = 256;
constexpr long PlainWords = 128;

class Pasta2_3
{
public:
    explicit Pasta2_3(long plainMod) : PlainMod(plainMod) {}

    void mix(vector<long> &state1, vector<long> &state2);
    void matmul(std::vector<long> &state, const std::vector<std::vector<uint64_t>> &matrix);
    void add_rc(std::vector<long> &state1, std::vector<long> &state2, size_t r);
    void random_matmul(std::vector<long> &state, const std::vector<std::vector<uint64_t>> &matrix);
    void random_add_rc(std::vector<long> &state);
    void sbox_cube(std::vector<long> &state);
    void sbox_feistel(std::vector<long> &state);

private:
    long PlainMod;
};

#endif // PASTA2_3_HPP