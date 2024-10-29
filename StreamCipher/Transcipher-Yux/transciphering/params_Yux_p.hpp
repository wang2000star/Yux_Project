#ifndef PARAMS_YUX_P_HPP
#define PARAMS_YUX_P_HPP

static const long PlainMod = 65537;//模数
static const long Bytebits = 16;//这里把明文字节看做16,因为2^16+1=65537                                              
static const unsigned BlockByte = 16;                               // 分组字节长度
static const unsigned BlockSize = Bytebits * BlockByte;                    // 分组比特长度
static int Nr = 5;//轮数
static int PlainBlock = 2;//明文分组数

static long mValues[][4] = { 
//{   p,       m,   bits}
  { 65537,  131072,  1320, 6}, // m=(3)*{257} 1250
  { 65537,  65536,  853, 17}, // m=(3)*{257}
};

#endif // PARAMS_YUX_P_HPP