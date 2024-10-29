#ifndef PARAMS_YUX_P_HPP
#define PARAMS_YUX_P_HPP


#include "Yux_p.hpp"

static uint64_t PlainMod = 65537;
static const long Bytebits = 16;//这里把明文字节看做16,因为2^16+1=65537                                              
constexpr uint64_t BlockByte = 16;// 分组字节长度
constexpr unsigned BlockSize = 16* BlockByte;                    // 分组比特长度

static int PlainBlock = 2;//明文分组数

static const uint64_t PlainByte = BlockByte * PlainBlock;           // 明文字节长度
// PlainByte = nslots
static const uint64_t PlainSize = BlockSize * PlainBlock;           // 明文比特长度
static const unsigned NonceSize = 32;                               // Nonce比特长度
static const unsigned Nr = 5;                                       // 轮数
static const uint64_t counter_begin = 0;                            // 计数器起始值
static const uint64_t counter_end = PlainBlock + counter_begin - 1; // 计数器结束值

static uint64_t RoundConstant = 0xCD;

static long mValues[][4] = { 
//{   p,       m,   bits}
  { 65537,  131072,  1320, 6}, // m=(3)*{257} 1250
  { 65537,  65536,  853, 17}, // m=(3)*{257}
};

#endif // PARAMS_YUX_P_HPP