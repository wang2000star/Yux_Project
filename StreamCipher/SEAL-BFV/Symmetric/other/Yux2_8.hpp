#ifndef YUX2_8_HPP
#define YUX2_8_HPP

#include <cstring>
#include <iostream>

// 常量定义
static const unsigned char roundConstant = 0xCD;
static const unsigned char roundConstant_4bit = 0xD; // x^3+x^2+1

// 函数声明
void addRoundKey(unsigned char state[], unsigned char RoundKey[], unsigned round);
unsigned char mul(unsigned char a, unsigned char b);
void decSboxFi(unsigned char state[], unsigned begin);
void decLinearLayer(unsigned char in[16]);
void decryption(unsigned char out[], unsigned char in[], unsigned char RoundKey[], unsigned Nr);
void encSboxFi(unsigned char state[], unsigned begin);
void encLinearLayer(unsigned char in[16]);
void encryption(unsigned char out[], unsigned char in[], unsigned char RoundKey[], unsigned Nr);
void constantForKey(unsigned char** RC, unsigned Nr);
void rotation(unsigned char *a, unsigned l, unsigned r);
unsigned KeyExpansion(unsigned char RoundKey[], unsigned Nr, unsigned blockByte, unsigned char Key[]);
void decRoundKey(unsigned char RoundKey_invert[], unsigned char RoundKey[], long Nr, long blockByte);

#endif // YUX2_8_HPP