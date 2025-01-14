#include "Pasta.hpp"

std::vector<std::vector<long>> Pasta::get_random_matrix_pasta(Keccak_HashInstance &shake128)
{
    std::vector<std::vector<long>> mat(BlockPlainWords, std::vector<long>(BlockPlainWords));
    for (auto &m : mat[0])
    {
        m = generate_random_field_element(shake128, false, max_prime_size, PlainMod);
    }
    const auto &first_row = mat[0];
    auto &prev_row = mat[0];
    for (auto i = 1ULL; i < BlockPlainWords; i++)
    {
        for (auto j = 0ULL; j < BlockPlainWords; j++)
        {
            long tmp =
                ((uint128_t)(first_row[j]) * prev_row[BlockPlainWords - 1]) % PlainMod;
            if (j)
            {
                tmp = (tmp + prev_row[j - 1]) % PlainMod;
            }
            mat[i][j] = tmp;
        }
        prev_row = mat[i];
    }
    return mat;
}
void Pasta::sbox_cube(vector<long> &state)
{
    for (uint16_t el = 0; el < BlockPlainWords; el++)
    {
        long square = ((uint128_t)(state[el]) * state[el]) % PlainMod;
        state[el] = ((uint128_t)(square)*state[el]) % PlainMod;
    }
}
void Pasta::sbox_feistel(vector<long> &state)
{
    vector<long> new_state(BlockPlainWords);
    new_state[0] = state[0];
    for (uint16_t el = 1; el < BlockPlainWords; el++)
    {
        long square = ((uint128_t)(state[el - 1]) * state[el - 1]) % PlainMod;
        // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
        new_state[el] = (square + state[el]) % PlainMod;
    }
    state = new_state;
}
void Pasta::matmul(vector<long> &state, const std::vector<std::vector<long>> &matrix)
{
    vector<long> new_state(BlockPlainWords, 0);
    for (uint16_t i = 0; i < BlockPlainWords; i++)
    {
        for (uint16_t j = 0; j < BlockPlainWords; j++)
        {
            long mult = ((uint128_t)(matrix[i][j]) * state[j]) % PlainMod;
            // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
            new_state[i] = (new_state[i] + mult) % PlainMod;
        }
    }
    state = new_state;
}
void Pasta::random_add_rc(vector<long> &state1, vector<long> &state2, vector<long> &rc1, vector<long> &rc2)
{
    for (uint16_t el = 0; el < BlockPlainWords; el++)
    {
        // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
        state1[el] = (state1[el] + rc1[el]) % PlainMod;
        state2[el] = (state2[el] + rc2[el]) % PlainMod;
    }
}
void Pasta::add_rc(vector<long> &state1, vector<long> &state2, size_t r, vector<vector<long>> &rc1, vector<vector<long>> &rc2)
{
    for (uint16_t el = 0; el < BlockPlainWords; el++)
    {
        // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
        state1[el] = (state1[el] + rc1[r][el]) % PlainMod;
        state2[el] = (state2[el] + rc2[r][el]) % PlainMod;
    }
}
void Pasta::mix(vector<long> &state1, vector<long> &state2)
{
    // just adding
    for (uint16_t i = 0; i < BlockPlainWords; i++)
    {
        // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
        long sum = (state1[i] + state2[i]);
        state1[i] = (state1[i] + sum) % PlainMod;
        state2[i] = (state2[i] + sum) % PlainMod;
    }
}