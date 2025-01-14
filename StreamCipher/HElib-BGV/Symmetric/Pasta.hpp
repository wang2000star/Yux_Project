#ifndef PASTA_HPP
#define PASTA_HPP

#include <iostream>
#include <array>
#include <vector>
#include <cmath>

#include "../utils/tool.hpp"

extern "C"
{
#include "../keccak/KeccakHash.h"
}

using namespace std;
typedef __uint128_t uint128_t;

class Pasta
{
    public:
        explicit Pasta(long plainMod, unsigned BlockPlainWords,long max_prime_size) : PlainMod(plainMod) , BlockPlainWords(BlockPlainWords) , max_prime_size(max_prime_size) {}

        std::vector<std::vector<long>> get_random_matrix_pasta(Keccak_HashInstance &shake128);
        void sbox_cube(vector<long> &state);
        void sbox_feistel(vector<long> &state);
        void matmul(vector<long> &state, const std::vector<std::vector<long>> &matrix);
        void random_add_rc(vector<long> &state1, vector<long> &state2, vector<long> &rc1, vector<long> &rc2);
        void add_rc(vector<long> &state1, vector<long> &state2, size_t r, vector<vector<long>> &rc1, vector<vector<long>> &rc2);
        void mix(vector<long> &state1, vector<long> &state2);

    private:
        long PlainMod;
        long BlockPlainWords;
        long max_prime_size; 
};

#endif // PASTA_HPP