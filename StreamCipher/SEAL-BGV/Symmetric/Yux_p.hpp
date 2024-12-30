#ifndef YUX_P_HPP
#define YUX_P_HPP

#include <cstdlib>
#include <cstdio>
#include <stdint.h>
#include <iostream>
#include <vector>
using namespace std;

class YuxP {
public:
    explicit YuxP(long plainMod,long roundConstant) : PlainMod(plainMod),RoundConstant(roundConstant) {}

    void Sbox(vector<long> &data);
    void Sbox2(vector<long> &data);
    void Linear(vector<long> &data);

private:
    long PlainMod;
    long RoundConstant;
};


#endif // YUX_P_HPP