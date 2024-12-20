#ifndef PARAMS_YUX2_8_HPP
#define PARAMS_YUX2_8_HPP

#include <stdint.h>
#include "random_bit.hpp"

#define PolyType DoubleCRT

static const unsigned BlockByte = 16;                               // 分组字节长度
static const unsigned BlockSize = 8 * BlockByte;                    // 分组比特长度
static const unsigned PlainBlock = 960;                              // 明文分组数
static const long PlainByte = BlockByte * PlainBlock;           // 明文字节长度
// PlainByte = nslots
static const long PlainSize = BlockSize * PlainBlock;           // 明文比特长度
static const unsigned NonceSize = 32;                               // Nonce比特长度
static const unsigned Nr = 5;                                       // 轮数
static const long counter_begin = 0;                            // 计数器起始值
static const long counter_end = PlainBlock + counter_begin - 1; // 计数器结束值

//static const uint8_t roundConstant = 0xCD; // x^7+x^6+x^3+x^2+1

static long mValues[][5] = { 
//{ p, phi(m),  m,   d, bits}
  { 2,  512,    771, 16, 960}, // 420, m=(3)*{257}
  { 2, 4096,   4369, 16, 360}, // m=17*(257)
  { 2, 16384, 21845, 16, 410}, // m=5*17*(257)
  { 2, 23040, 28679, 24, 700}, // m=7*17*(241) bits: 600
  { 2, 46080, 53261, 24, 1330}, // m=13*17*(241) bits: 630 
  { 2, 64512, 65281, 48, 480}  // m=97*(673)
};


#endif // PARAMS_HPP