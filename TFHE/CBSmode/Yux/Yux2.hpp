#ifndef YUX2_HPP
#define YUX2_HPP

#include <cstring>

// The round number is ROUND = 12
static const unsigned char roundConstant = 0xCD;
static const unsigned char roundConstant_4bit = 0xD; // x^3+x^2+1

// Function declarations
void addRoundKey(unsigned char state[], unsigned char RoundKey[], int round);
unsigned char mul(unsigned char a, unsigned char b);
void decSboxFi(unsigned char state[], int begin);
void decLinearLayer(unsigned char in[16]);
void decryption(unsigned char out[], unsigned char in[], unsigned char RoundKey[], int Nr);
void encSboxFi(unsigned char state[], int begin);
void encLinearLayer(unsigned char in[16]);
void encryption(unsigned char out[], unsigned char in[], unsigned char RoundKey[], int Nr);
void constantForKey(unsigned char RC[56][4], long round);
void rotation(unsigned char *a, int l, int r);
int KeyExpansion(unsigned char RoundKey[], int round, int blockByte, unsigned char Key[]);
void decRoundKey(unsigned char RoundKey_invert[], unsigned char RoundKey[], long Nr, long blockByte);

#endif // YUX2_HPP