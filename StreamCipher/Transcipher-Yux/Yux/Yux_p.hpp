#ifndef YUX_P_HPP
#define YUX_P_HPP

#include <cstdlib>
#include <cstdio>
#include <stdint.h>
#include <iostream>
#include <vector>
using namespace std;

static long roundConstant = 0xCD;

static long pmod = 65537;

static long blockByte = 16;

static long ROUND = 4;

void decSboxFi(long state[], int begin);

void decLinearLayer(long in[16]);
/*
class Yux_F_p{

  private:
    int ROUND=4;
    int BlockWords = 16;
    long Yux_p = 65537;
    long modulus = 65537;
    

    long Model_p(long state);
    void addRoundKey(long state[], long RoundKey[], int round);
    void encSboxFi(long state[], int begin);
    void encLinearLayer(long in[16]);
    
    void subtractRoundKey(long state[], long RoundKey[], int round);
    void decSboxFi(long state[], int begin);
    void decLinearLayer(long in[16]);

    void rotation(long *a, int l,int r);
    void constantForKey(long RC[56][4]);


  public:
    Yux_F_p(const int b, const int r, const long p): BlockWords(b), ROUND(r), modulus(p){}
    // YUX_F_P() = default;
    void decryption(long out[], long in[], long RoundKey[]);
    void encryption(long out[], long in[], long RoundKey[]);
    long KeyExpansion(long RoundKey[], long Key[]);
    void decRoundKey(long RoundKey_invert[], long RoundKey[]);

    vector<vector<long>> matrix = {
      { 9363, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9362, 53053, 53053, 53053, 9362, 53054, 53054, 53054 },
      { 53054, 9363, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9362, 53053, 53053, 53053, 9362, 53054, 53054 },
      { 53054, 53054, 9363, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9362, 53053, 53053, 53053, 9362, 53054 },
      { 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9362, 53053, 53053, 53053, 9362 },
      { 9362, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9362, 53053, 53053, 53053 },
      { 53053, 9362, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9362, 53053, 53053 },
      { 53053, 53053, 9362, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9362, 53053 },
      { 53053, 53053, 53053, 9362, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9362 },
      { 9362, 53053, 53053, 53053, 9362, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9363, 53054, 53054, 53054 },
      { 53054, 9362, 53053, 53053, 53053, 9362, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9363, 53054, 53054 },
      { 53054, 53054, 9362, 53053, 53053, 53053, 9362, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9363, 53054 },
      { 53054, 53054, 53054, 9362, 53053, 53053, 53053, 9362, 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9363 },
      { 9363, 53054, 53054, 53054, 9362, 53053, 53053, 53053, 9362, 53054, 53054, 53054, 9363, 53054, 53054, 53054 },
      { 53054, 9363, 53054, 53054, 53054, 9362, 53053, 53053, 53053, 9362, 53054, 53054, 53054, 9363, 53054, 53054 },
      { 53054, 53054, 9363, 53054, 53054, 53054, 9362, 53053, 53053, 53053, 9362, 53054, 53054, 53054, 9363, 53054 },
      { 53054, 53054, 53054, 9363, 53054, 53054, 53054, 9362, 53053, 53053, 53053, 9362, 53054, 53054, 53054, 9363 }
    };

};
*/

#endif // YUX_P_HPP