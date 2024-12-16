

#include "Hera.hpp"

using namespace std;

//static long PlainMod = 257;
//static long PlainMod = 65537;



// 列混淆：A = constMC * A,A是大小为32的一维数组，按照列优先看作8x4的矩阵
void Hera::MC(vector<long> &A)
{
    std::vector<long> temp(16);
    std::array<int, 8> index = {0, 1, 2, 3};
    for (int i = 0; i < 4; i++)
    {
        int s = 4 * i;
        temp[index[0] + s] = ((2*A[index[0] + s] + 3*A[index[1] + s] + A[index[2] + s] + A[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[1] + s] = ((A[index[0] + s] + 2*A[index[1] + s] + 3*A[index[2] + s] + A[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[2] + s] = ((A[index[0] + s] + A[index[1] + s] + 2*A[index[2] + s] + 3*A[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[3] + s] = ((3*A[index[0] + s] + A[index[1] + s] + A[index[2] + s] + 2*A[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
           }
    std::copy(temp.begin(), temp.end(), A.begin());
}

// 行移位：A = A * constMR
void Hera::MR(vector<long> &A)
{
    std::vector<long> temp(16);
    std::array<int, 8> index = {0, 4, 8, 12};
    for (int i = 0; i < 4; i++)
    {
        int s = i;
        temp[index[0] + s] = ((2*A[index[0] + s] + 3*A[index[1] + s] + A[index[2] + s] + A[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[1] + s] = ((A[index[0] + s] + 2*A[index[1] + s] + 3*A[index[2] + s] + A[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[2] + s] = ((A[index[0] + s] + A[index[1] + s] + 2*A[index[2] + s] + 3*A[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[3] + s] = ((3*A[index[0] + s] + A[index[1] + s] + A[index[2] + s] + 2*A[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
           }
    std::copy(temp.begin(), temp.end(), A.begin());
}

void Hera::Sbox(vector<long> &A)
{
    long temp;
    for (int i = 0; i < A.size(); i ++)
    {
        A[i] = (A[i] * A[i] * A[i]) % PlainMod;
    }
}