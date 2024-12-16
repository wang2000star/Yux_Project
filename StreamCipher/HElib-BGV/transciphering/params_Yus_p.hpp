#ifndef PARAMS_YUS_P_HPP
#define PARAMS_YUS_P_HPP

/************
 * PlainMod = 257, Bytebits = 9, BlockSize = 32*9 = 288, PlainBlock = (Plainmod-1)/2 = 128
 * PlainByte = 65537, Bytebits = 17, BlockSize = 32*17 = 544, PlainBlock = (Plainmod-1)/2 = 32768
 */
// constexpr long PlainMod = 257;    //2^8+1
// constexpr unsigned Bytebits = 9;       // 字节比特长度=ceil(log2(PlainMod))  
constexpr long PlainMod = 65537;    //2^16+1
constexpr unsigned Bytebits = 17;       // 字节比特长度=ceil(log2(PlainMod))  

constexpr long BlockByte = 32;       // 分组字节长度
constexpr unsigned BlockSize = Bytebits * BlockByte;      // 分组比特长度=BlockByte*Bytebits
constexpr unsigned PlainBlock = (PlainMod-1)/2;//明文分组数

static const long PlainByte = BlockByte * PlainBlock;           // 明文字节长度
// PlainByte = nslots
static const long PlainSize = BlockSize * PlainBlock;           // 明文比特长度
static const unsigned NonceSize = 32;                               // Nonce比特长度
static const unsigned Nr = 6;                                       // 轮数
static const long counter_begin = 0;                            // 计数器起始值
static const long counter_end = PlainBlock + counter_begin - 1; // 计数器结束值

static long RoundConstant = 0xCD;

/************************************************************************
  long m;          // m-th cyclotomic polynomial
  long p;          // plaintext primeplain_mod;
  long r;          // Lifting [defualt = 1]
  long L;          // bits in the ciphertext modulus chain
  long c;          // columns in the key-switching matrix [default=2]
  long d;          // Degree of the field extension [default=1]
  long k;          // Security parameter [default=80]
  long s;          // Minimum number of slots [default=0]
************************************************************************/
static long mValues[][4] = { 
//{   p,       m,    L,    c}
  { 65537,  131072,  1320, 6}, 
  //{ 257,    256,     320,  23},
  { 65537,  65536,   853,  17}, 
};
// p^d = 1 mod m,d=1,slots=\phi(m)/d=\phi(m);m=65536=2^16,\phi(m)=2^15=32768

#endif // PARAMS_YUS_P_HPP