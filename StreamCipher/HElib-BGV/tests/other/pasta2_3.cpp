#include "pasta2_3.hpp"
#include <iostream>

using namespace std;

typedef __uint128_t uint128_t;

void Pasta2_3::mix(vector<long> &state1, vector<long> &state2) {
  // just adding
  for (uint16_t i = 0; i < PlainWords; i++) {
    // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
    uint64_t sum = (state1[i] + state2[i]);
    state1[i] = (state1[i] + sum) % PlainMod;
    state2[i] = (state2[i] + sum) % PlainMod;
  }
}
//----------------------------------------------------------------

void Pasta2_3::random_add_rc(std::vector<long>& state) {
  for (uint16_t el = 0; el < PlainWords; el++) {
    // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
    state[el] = (state[el] + generate_random_field_element()) % PlainMod;
  }
}

//----------------------------------------------------------------
// requires storage of two rows in the matrix
void Pasta2_3::random_matmul(std::vector<long>& state, const std::vector<std::vector<uint64_t>>& matrix) {
  matmul(state, matrix);
  for (uint16_t el = 0; el < PlainWords; el++) {
    state[el] =
        ((uint128_t)(state[el]) * generate_random_field_element(false)) %
        PlainMod;
  }
}



void Pasta2_3::random_linear_layer(vector<long> &state1, vector<long> &state2) {
  random_matmul(state1, instance.matrix_FL);
  random_matmul(state2, instance.matrix_FR);
  random_add_rc(state1);
  random_add_rc(state2);
  mix();
}

void Pasta2_3::add_rc(std::vector<long>& state1, std::vector<long>& state2, size_t r) {
  for (uint16_t el = 0; el < PlainWords; el++) {
    // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
    state1[el] = (state1[el] + instance.rc1[r][el]) % PlainMod;
    state2[el] = (state2[el] + instance.rc2[r][el]) % PlainMod;
  }
}

//----------------------------------------------------------------


//----------------------------------------------------------------
// (2 1)(state1)
// (1 2)(state2)



//----------------------------------------------------------------
void Pasta2_3::matmul(std::vector<long>& state, const std::vector<std::vector<uint64_t>>& matrix) {
  std::vector<long> new_state;
  new_state.fill(0);

  for (uint16_t i = 0; i < PlainWords; i++) {
    for (uint16_t j = 0; j < PlainWords; j++) {
      uint64_t mult =
          ((uint128_t)(matrix[i][j]) * state[j]) % PlainMod;
      // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
      new_state[i] = (new_state[i] + mult) % PlainMod;
    }
  }
  state = new_state;
}

void Pasta2_3::fixed_linear_layer(size_t r,vector<long> &state1, vector<long> &state2,const std::vector<std::vector<uint64_t>>& matrix) {
  matmul(state1, instance.matrix);
  matmul(state2, instance.matrix);
  add_rc(state1, state2, r);
  mix();
}

//----------------------------------------------------------------

void Pasta2_3::sbox_feistel(std::vector<long>& state) {
  std::vector<long> new_state;
  new_state[0] = state[0];
  for (uint16_t el = 1; el < PlainWords; el++) {
    uint64_t square = ((uint128_t)(state[el - 1]) * state[el - 1]) % PlainMod;
    // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
    new_state[el] = (square + state[el]) % PlainMod;
  }
  state = new_state;
}

void Pasta2_3::sbox_cube(std::vector<long>& state) {
  for (uint16_t el = 0; el < PlainWords; el++) {
    uint64_t square = ((uint128_t)(state[el]) * state[el]) % PlainMod;
    state[el] = ((uint128_t)(square)*state[el]) % PlainMod;
  }
}

std::vector<uint64_t> Pasta2_3::get_random_yi( Keccak_HashInstance& shake128) const {
  std::vector<uint64_t> out(PlainWords, 0);
  for (auto i = 0ULL; i < PlainWords; i++) {
    bool valid = false;
    // get distinct x_i, which also imply distinct y_i
    while (!valid) {
      // random element of size floor(log_2(p))
      uint8_t random_bytes[sizeof(uint64_t)];
      if (SUCCESS !=
          Keccak_HashSqueeze(&shake128, random_bytes, sizeof(random_bytes) * 8))
        throw std::runtime_error("SHAKE128 squeeze failed");
      uint64_t y_i = be64toh(*((uint64_t*)random_bytes)) & y_mask;
      // check distinct x_i
      uint64_t x_i = y_i & x_mask;
      valid = true;
      for (auto j = 0ULL; j < i; j++) {
        if ((out[j] & x_mask) == x_i) {
          valid = false;
          break;
        }
      }
      if (valid) out[i] = y_i;
    }
  }
  return out;
}

//----------------------------------------------------------------
uint64_t Pasta2Instance::mod_inverse(uint64_t val) const {
  if (val == 0) throw(std::runtime_error("0 has no inverse!"));

  int64_t prev_a = 1;
  int64_t a = 0;
  uint64_t mod = modulus;

  while (mod != 0) {
    int64_t q = val / mod;
    int64_t temp = val % mod;
    val = mod;
    mod = temp;

    temp = a;
    a = prev_a - q * a;
    prev_a = temp;
  }

  if (val != 1) throw(std::runtime_error("value has no inverse!"));

  uint64_t res = prev_a;
  if (prev_a < 0) res += modulus;

  return res;
}
std::vector<std::vector<uint64_t>> Pasta2_3::get_random_mds_matrix(
    Keccak_HashInstance& shake128) const {
  std::vector<std::vector<uint64_t>> mat(PASTA2_T,
                                         std::vector<uint64_t>(PASTA2_T, 0));
  auto y = get_random_yi(shake128);
  auto x = y;
  for (auto& x_i : x) x_i &= x_mask;

  for (auto i = 0ULL; i < PASTA2_T; i++) {
    for (auto j = 0ULL; j < PASTA2_T; j++) {
      // ld(rasta_prime) ~ 60, no uint128_t for addition necessary
      mat[i][j] = mod_inverse((x[i] + y[j]) % modulus);
    }
  }
  return mat;
}

// void Pasta2_3::round(size_t r,vector<long> &state1, vector<long> &state2) {
//   if (r == PASTA2_R - 1) {
//     sbox_cube(state1);
//     sbox_cube(state2);
//   } else {
//     sbox_feistel(state1);
//     sbox_feistel(state2);
//   }
//   fixed_linear_layer(r);
// }

// void Pasta2_3::init_shake(uint64_t nonce, uint64_t block_counter) {
//   uint8_t seed[16];

//   *((uint64_t*)seed) = htobe64(nonce);
//   *((uint64_t*)(seed + 8)) = htobe64(block_counter);

//   if (SUCCESS != Keccak_HashInitialize_SHAKE128(&shake128_))
//     throw std::runtime_error("failed to init shake");
//   if (SUCCESS != Keccak_HashUpdate(&shake128_, seed, sizeof(seed) * 8))
//     throw std::runtime_error("SHAKE128 update failed");
//   if (SUCCESS != Keccak_HashFinal(&shake128_, NULL))
//     throw std::runtime_error("SHAKE128 final failed");
// }

// std::vector<long> Pasta2_3::gen_keystream(const uint64_t nonce,
//                             const uint64_t block_counter) {
//   init_shake(nonce, block_counter);

//   // init state
//   for (uint16_t i = 0; i < PlainWords; i++) {
//     state1[i] = key_[i];
//     state2[i] = key_[PlainWords + i];
//   }

//   // first random affine with mixing afterwards
//   random_linear_layer();
//   for (uint8_t r = 0; r < PASTA2_R; r++) {
//     round(r);
//   }
//   return state1;
// }

// //----------------------------------------------------------------

// std::vector<long> Pasta2_3::keystream(const uint64_t nonce, const uint64_t block_counter) {
//   return gen_keystream(nonce, block_counter);
// }
