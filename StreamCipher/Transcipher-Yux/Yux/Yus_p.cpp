

#include "Yus_p.hpp"

using namespace std;

//static long pmod = 257;
static long pmod = 65537;

// 一维数组可以看作按照列优先存储的矩阵
const std::array<long, 64> constMC = {0, 1, 0, 1, 1, 1, 65536, 65536,
                                          1, 0, 1, 1, 1, 0, 1, 0,
                                          1, 0, 1, 0, 1, 1, 65536, 1,
                                          0, 1, 1, 1, 0, 1, 1, 1,
                                          1, 1, 1, 0, 1, 0, 1, 1,
                                          1, 0, 0, 1, 0, 1, 0, 1,
                                          1, 0, 65536, 65536, 65536, 1, 1, 65536,
                                          0, 65536, 65536, 0, 1, 0, 0, 1};

// constMR 是 constMC 的转置矩阵
const std::array<long, 64> constMR = {0, 1, 1, 0, 1, 1, 1, 0,
                                          1, 0, 0, 1, 1, 0, 0, 65536,
                                          0, 1, 1, 1, 1, 0, 65536, 65536,
                                          1, 1, 0, 1, 0, 1, 65536, 0,
                                          1, 1, 1, 0, 1, 0, 65536, 1,
                                          1, 0, 1, 1, 0, 1, 1, 0,
                                          65536, 1, 65536, 1, 1, 0, 1, 0,
                                          65536, 0, 1, 1, 1, 1, 65536, 1};

// 列混淆：A = constMC * A,A是大小为32的一维数组，按照列优先看作8x4的矩阵
void MC(vector<long> &A)
{
    std::vector<long> temp(32);
    std::array<int, 8> index = {0, 1, 2, 3, 4, 5, 6, 7};
    for (int i = 0; i < 4; i++)
    {
        int s = 8 * i;
        temp[index[0] + s] = ((A[index[1] + s] + A[index[2] + s] + A[index[4] + s] + A[index[5] + s] + A[index[6] + s]) % pmod + pmod) % pmod;
        temp[index[1] + s] = ((A[index[0] + s] + A[index[3] + s] + A[index[4] + s] - A[index[7] + s]) % pmod + pmod) % pmod;
        temp[index[2] + s] = ((A[index[1] + s] + A[index[2] + s] + A[index[3] + s] + A[index[4] + s] - A[index[6] + s] - A[index[7] + s]) % pmod + pmod) % pmod;
        temp[index[3] + s] = ((A[index[0] + s] + A[index[1] + s] + A[index[3] + s] + A[index[5] + s] - A[index[6] + s]) % pmod + pmod) % pmod;
        temp[index[4] + s] = ((A[index[0] + s] + A[index[1] + s] + A[index[2] + s] + A[index[4] + s] - A[index[6] + s] + A[index[7] + s]) % pmod + pmod) % pmod;
        temp[index[5] + s] = ((A[index[0] + s] + A[index[2] + s] + A[index[3] + s] + A[index[5] + s] + A[index[6] + s]) % pmod + pmod) % pmod;
        temp[index[6] + s] = ((A[index[1] + s] - A[index[0] + s] - A[index[2] + s] + A[index[3] + s] + A[index[4] + s] + A[index[6] + s]) % pmod + pmod) % pmod;
        temp[index[7] + s] = ((A[index[2] + s] - A[index[0] + s] + A[index[3] + s] + A[index[4] + s] + A[index[5] + s] - A[index[6] + s] + A[index[7] + s]) % pmod + pmod) % pmod; 
    }
    std::copy(temp.begin(), temp.end(), A.begin());
}

// 行移位：A = A * constMR
void MR(vector<long> &A)
{
    std::vector<long> temp(32);
    std::array<int, 8> index = {0, 1, 8, 9, 16, 17, 24, 25};
    for (int i = 0; i < 4; i++)
    {
        int s = 2 * i;
        temp[index[0] + s] = ((A[index[1] + s] + A[index[2] + s] + A[index[4] + s] + A[index[5] + s] + A[index[6] + s]) % pmod + pmod) % pmod;
        temp[index[1] + s] = ((A[index[0] + s] + A[index[3] + s] + A[index[4] + s] - A[index[7] + s]) % pmod + pmod) % pmod;
        temp[index[2] + s] = ((A[index[1] + s] + A[index[2] + s] + A[index[3] + s] + A[index[4] + s] - A[index[6] + s] - A[index[7] + s]) % pmod + pmod) % pmod;
        temp[index[3] + s] = ((A[index[0] + s] + A[index[1] + s] + A[index[3] + s] + A[index[5] + s] - A[index[6] + s]) % pmod + pmod) % pmod;
        temp[index[4] + s] = ((A[index[0] + s] + A[index[1] + s] + A[index[2] + s] + A[index[4] + s] - A[index[6] + s] + A[index[7] + s]) % pmod + pmod) % pmod;
        temp[index[5] + s] = ((A[index[0] + s] + A[index[2] + s] + A[index[3] + s] + A[index[5] + s] + A[index[6] + s]) % pmod + pmod) % pmod;
        temp[index[6] + s] = ((A[index[1] + s] - A[index[0] + s] - A[index[2] + s] + A[index[3] + s] + A[index[4] + s] + A[index[6] + s]) % pmod + pmod) % pmod;
        temp[index[7] + s] = ((A[index[2] + s] - A[index[0] + s] + A[index[3] + s] + A[index[4] + s] + A[index[5] + s] - A[index[6] + s] + A[index[7] + s]) % pmod + pmod) % pmod; 
 }
    std::copy(temp.begin(), temp.end(), A.begin());
}

void Sbox(vector<long> &A)
{
    long temp;
    for (int i = 0; i < A.size(); i += 2)
    {
        // 两次Feistel轮函数(x,y)——> (y+x^2,x)
        for (int j = 0; j < 2; j++)
        {
            temp = A[i];
            A[i] = ((A[i + 1] + A[i] * A[i]) % pmod+ pmod)%pmod;
            A[i + 1] = temp;
        }
    }
}