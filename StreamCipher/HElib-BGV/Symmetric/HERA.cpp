

#include "HERA.hpp"

using namespace std;

void HERA::MC(vector<long> &A)
{
    vector<long> temp = A;
    array<int, 8> index = {0, 1, 2, 3};
    for (int i = 0; i < 4; i++)
    {
        int s = 4 * i;
        A[index[0] + s] = ((2 * temp[index[0] + s] + 3 * temp[index[1] + s] + temp[index[2] + s] + temp[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
        A[index[1] + s] = ((temp[index[0] + s] + 2 * temp[index[1] + s] + 3 * temp[index[2] + s] + temp[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
        A[index[2] + s] = ((temp[index[0] + s] + temp[index[1] + s] + 2 * temp[index[2] + s] + 3 * temp[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
        A[index[3] + s] = ((3 * temp[index[0] + s] + temp[index[1] + s] + temp[index[2] + s] + 2 * temp[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
    }
}

void HERA::MR(vector<long> &A)
{
    vector<long> temp = A;
    array<int, 8> index = {0, 4, 8, 12};
    for (int i = 0; i < 4; i++)
    {
        int s = i;
        A[index[0] + s] = ((2 * temp[index[0] + s] + 3 * temp[index[1] + s] + temp[index[2] + s] + temp[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
        A[index[1] + s] = ((temp[index[0] + s] + 2 * temp[index[1] + s] + 3 * temp[index[2] + s] + temp[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
        A[index[2] + s] = ((temp[index[0] + s] + temp[index[1] + s] + 2 * temp[index[2] + s] + 3 * temp[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
        A[index[3] + s] = ((3 * temp[index[0] + s] + temp[index[1] + s] + temp[index[2] + s] + 2 * temp[index[3] + s]) % PlainMod + PlainMod) % PlainMod;
    }
}

void HERA::Sbox(vector<long> &A)
{
    for (int i = 0; i < A.size(); i++)
    {
        A[i] = (((A[i] * A[i]) % PlainMod) * A[i]) % PlainMod;
    }
}