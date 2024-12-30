

#include "Yus_p.hpp"

using namespace std;

// static uint64_t PlainMod = 257;
// static uint64_t PlainMod = 65537;

// 一维数组可以看作按照列优先存储的矩阵
const std::array<uint64_t, 64> constMC = {0, 1, 0, 1, 1, 1, 65536, 65536,
                                      1, 0, 1, 1, 1, 0, 1, 0,
                                      1, 0, 1, 0, 1, 1, 65536, 1,
                                      0, 1, 1, 1, 0, 1, 1, 1,
                                      1, 1, 1, 0, 1, 0, 1, 1,
                                      1, 0, 0, 1, 0, 1, 0, 1,
                                      1, 0, 65536, 65536, 65536, 1, 1, 65536,
                                      0, 65536, 65536, 0, 1, 0, 0, 1};

// constMR 是 constMC 的转置矩阵
const std::array<uint64_t, 64> constMR = {0, 1, 1, 0, 1, 1, 1, 0,
                                      1, 0, 0, 1, 1, 0, 0, 65536,
                                      0, 1, 1, 1, 1, 0, 65536, 65536,
                                      1, 1, 0, 1, 0, 1, 65536, 0,
                                      1, 1, 1, 0, 1, 0, 65536, 1,
                                      1, 0, 1, 1, 0, 1, 1, 0,
                                      65536, 1, 65536, 1, 1, 0, 1, 0,
                                      65536, 0, 1, 1, 1, 1, 65536, 1};

// 列混淆：A = constMC * A,A是大小为32的一维数组，按照列优先看作8x4的矩阵
void YusP::MC32(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp(32);
    std::array<int, 8> index = {0, 1, 2, 3, 4, 5, 6, 7};
    for (int i = 0; i < 4; i++)
    {
        int s = 8 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s;
        temp[id0] = ((A[id1] + A[id2] + A[id4] + A[id5] + A[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id1] = ((A[id0] + A[id3] + A[id4] - A[id7]) % PlainMod + PlainMod) % PlainMod;
        temp[id2] = ((A[id1] + A[id2] + A[id3] + A[id4] - A[id6] - A[id7]) % PlainMod + PlainMod) % PlainMod;
        temp[id3] = ((A[id0] + A[id1] + A[id3] + A[id5] - A[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id4] = ((A[id0] + A[id1] + A[id2] + A[id4] - A[id6] + A[id7]) % PlainMod + PlainMod) % PlainMod;
        temp[id5] = ((A[id0] + A[id2] + A[id3] + A[id5] + A[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id6] = ((A[id1] - A[id0] - A[id2] + A[id3] + A[id4] + A[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id7] = ((A[id2] - A[id0] + A[id3] + A[id4] + A[id5] - A[id6] + A[id7]) % PlainMod + PlainMod) % PlainMod;
    }
    std::copy(temp.begin(), temp.end(), A.begin());
}

// 行移位：A = A * constMR
void YusP::MR32(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp(32);
    std::array<int, 8> index = {0, 1, 8, 9, 16, 17, 24, 25};
    for (int i = 0; i < 4; i++)
    {
        int s = 2 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s;
        temp[id0] = ((A[id1] + A[id2] + A[id4] + A[id5] + A[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id1] = ((A[id0] + A[id3] + A[id4] - A[id7]) % PlainMod + PlainMod) % PlainMod;
        temp[id2] = ((A[id1] + A[id2] + A[id3] + A[id4] - A[id6] - A[id7]) % PlainMod + PlainMod) % PlainMod;
        temp[id3] = ((A[id0] + A[id1] + A[id3] + A[id5] - A[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id4] = ((A[id0] + A[id1] + A[id2] + A[id4] - A[id6] + A[id7]) % PlainMod + PlainMod) % PlainMod;
        temp[id5] = ((A[id0] + A[id2] + A[id3] + A[id5] + A[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id6] = ((A[id1] - A[id0] - A[id2] + A[id3] + A[id4] + A[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id7] = ((A[id2] - A[id0] + A[id3] + A[id4] + A[id5] - A[id6] + A[id7]) % PlainMod + PlainMod) % PlainMod;
    }
    std::copy(temp.begin(), temp.end(), A.begin());
}
void YusP::MC_MR(vector<uint64_t> &A)
{
    vector<vector<int>> M = {
        {1, 0, 0, 1, 1, 0, 0, -1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 2, 1, 1, -1, 0, 1, 1, 0, 1, 1, 1, 0},
        {0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, -1, 0, 1, 1, 0, 1, 1, 1, 0, -1, 0, 0, -1, -1, 0, 0, 1},
        {1, 1, 0, 1, 0, 1, -1, 0, 0, 1, 1, 1, 1, 0, -1, -1, 1, 2, 1, 2, 1, 1, -2, -1, 0, 1, 1, 1, 1, 0, -1, -1},
        {0, 1, 1, 1, 1, 0, -1, -1, 1, 1, 0, 1, 0, 1, -1, 0, 0, 1, 1, 1, 1, 0, -1, -1, -1, -1, 0, -1, 0, -1, 1, 0},
        {1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, -1, 1, 2, 1, 2, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, -1, 1},
        {1, 1, 1, 0, 1, 0, -1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, -1, 1, -1, 0, -1, -1, 0, -1, -1, 0},
        {-1, 0, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 0, 1, 0, -2, 1, 0, 2, 2, 1, 0, 1, -1, 1, -1, 1, 1, 0, 1, 0},
        {-1, 1, -1, 1, 1, 0, 1, 0, -1, 0, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 0, 1, 0, 1, 0, -1, -1, -1, -1, 1, -1},
        {1, 0, 0, 1, 1, 0, 0, -1, 1, 1, 1, 1, 2, 1, 1, -1, 0, 1, 1, 0, 1, 1, 1, 0, -1, -1, -1, -1, -2, -1, -1, 1},
        {1, 1, 1, 1, 2, 1, 1, -1, 1, 0, 0, 1, 1, 0, 0, -1, 1, 0, 0, 1, 1, 0, 0, -1, 0, -1, -1, 0, -1, -1, -1, 0},
        {1, 1, 0, 1, 0, 1, -1, 0, 1, 2, 1, 2, 1, 1, -2, -1, 0, 1, 1, 1, 1, 0, -1, -1, -1, -2, -1, -2, -1, -1, 2, 1},
        {1, 2, 1, 2, 1, 1, -2, -1, 1, 1, 0, 1, 0, 1, -1, 0, 1, 1, 0, 1, 0, 1, -1, 0, 0, -1, -1, -1, -1, 0, 1, 1},
        {1, 0, 1, 1, 0, 1, 1, 0, 2, 1, 2, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, -1, 1, -2, -1, -2, -1, -1, -1, 0, -1},
        {2, 1, 2, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, -1, -1, -1, 0, -1, 0, 1, -1},
        {-1, 0, 1, 1, 1, 1, -1, 1, -2, 1, 0, 2, 2, 1, 0, 1, -1, 1, -1, 1, 1, 0, 1, 0, 2, -1, 0, -2, -2, -1, 0, -1},
        {-2, 1, 0, 2, 2, 1, 0, 1, -1, 0, 1, 1, 1, 1, -1, 1, -1, 0, 1, 1, 1, 1, -1, 1, 1, -1, 1, -1, -1, 0, -1, 0},
        {1, 1, 1, 1, 2, 1, 1, -1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, -1, -1, 1, 0, -1, -1, -1},
        {0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 2, 1, 1, -1, 1, 0, 0, 1, 1, 0, 0, -1, 0, 1, 1, 0, 1, 1, 1, 0},
        {1, 2, 1, 2, 1, 1, -2, -1, 0, 1, 1, 1, 1, 0, -1, -1, 0, 1, 1, 1, 1, 0, -1, -1, 1, 0, -1, 0, -1, 1, 0, 1},
        {0, 1, 1, 1, 1, 0, -1, -1, 1, 2, 1, 2, 1, 1, -2, -1, 1, 1, 0, 1, 0, 1, -1, 0, 0, 1, 1, 1, 1, 0, -1, -1},
        {2, 1, 2, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, -1, 1, 1, 1, 1, 0, 1, 0, -1, 1, 0, -1, 0, 1, -1, 1, 2, -1},
        {1, 1, 1, 0, 1, 0, -1, 1, 2, 1, 2, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, -1, 1},
        {-2, 1, 0, 2, 2, 1, 0, 1, -1, 1, -1, 1, 1, 0, 1, 0, -1, 1, -1, 1, 1, 0, 1, 0, 0, -1, 2, 0, 0, 1, -2, 1},
        {-1, 1, -1, 1, 1, 0, 1, 0, -2, 1, 0, 2, 2, 1, 0, 1, -1, 0, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 0, 1, 0},
        {1, -1, -1, 1, 0, -1, -1, -1, 1, -1, -1, 1, 0, -1, -1, -1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0},
        {0, -1, -1, 0, -1, -1, -1, 0, 1, 1, 1, 1, 2, 1, 1, -1, 1, 1, 1, 1, 2, 1, 1, -1, 1, -1, -1, 1, 0, -1, -1, -1},
        {1, 0, -1, 0, -1, 1, 0, 1, 1, 0, -1, 0, -1, 1, 0, 1, 0, 1, 1, 1, 1, 0, -1, -1, 0, 1, 1, 1, 1, 0, -1, -1},
        {0, -1, -1, -1, -1, 0, 1, 1, 1, 2, 1, 2, 1, 1, -2, -1, 1, 2, 1, 2, 1, 1, -2, -1, 1, 0, -1, 0, -1, 1, 0, 1},
        {0, -1, 0, 1, -1, 1, 2, -1, 0, -1, 0, 1, -1, 1, 2, -1, 1, 1, 1, 0, 1, 0, -1, 1, 1, 1, 1, 0, 1, 0, -1, 1},
        {-1, -1, -1, 0, -1, 0, 1, -1, 2, 1, 2, 1, 1, 1, 0, 1, 2, 1, 2, 1, 1, 1, 0, 1, 0, -1, 0, 1, -1, 1, 2, -1},
        {0, -1, 2, 0, 0, 1, -2, 1, 0, -1, 2, 0, 0, 1, -2, 1, -1, 1, -1, 1, 1, 0, 1, 0, -1, 1, -1, 1, 1, 0, 1, 0},
        {1, -1, 1, -1, -1, 0, -1, 0, -2, 1, 0, 2, 2, 1, 0, 1, -2, 1, 0, 2, 2, 1, 0, 1, 0, -1, 2, 0, 0, 1, -2, 1}};
    std::vector<uint64_t> temp(32);
    for (int i = 0; i < 32; i++)
    {
        temp[i] = 0;
        for (int j = 0; j < 32; j++)
        {
            temp[i] += M[i][j] * A[j];
        }
        temp[i] = (temp[i] % PlainMod + PlainMod) % PlainMod;
    }
    std::copy(temp.begin(), temp.end(), A.begin());
}
/*
M=
[1,1,2,1,2,1,2,2]

[2,1,0,2,1,2,1,2]

[2,2,1,1,2,1,2,1]

[1,2,2,1,0,2,1,2]

[2,1,2,2,1,1,2,1]

[1,2,1,2,2,1,0,2]

[2,1,2,1,2,2,1,1]

[0,2,1,2,1,2,2,1]
...
[1,1,2,1,2,1,2,2]
A[id0] + A[id1] + t[2] + A[id3] + t[4] + A[id5] + t[6] + t[7];
[2,1,0,2,1,2,1,2]
t[0] + A[id1]  + t[3] + A[id4]+t[5] + A[id6] + t[7];
[2,2,1,1,2,1,2,1]
t[0] + t[1] + A[id2] + A[id3] + t[4]+A[id5] + t[6] + A[id7];
[1,2,2,1,0,2,1,2]
A[id0] + t[1] + t[2] + A[id3] + t[5] + A[id6] + t[7];
[2,1,2,2,1,1,2,1]
t[0] + A[id1] + t[2] + t[3] + A[id4]+A[id5] + t[6] + A[id7];
[1,2,1,2,2,1,0,2]
A[id0] + t[1] + A[id2] + t[3] + t[4]+A[id5] + t[7];
[2,1,2,1,2,2,1,1]
t[0] + A[id1] + t[2] + A[id3] + t[4]+t[5] + A[id6] + A[id7];
[0,2,1,2,1,2,2,1]
t[1] + A[id2] + t[3] + A[id4]+t[5] + t[6] + A[id7];

        t[0] = (A[id0]+A[id0])%PlainMod;
        t[1] = (A[id1]+A[id1])%PlainMod;
        t[2] = (A[id2]+A[id2])%PlainMod;
        t[3] = (A[id3]+A[id3])%PlainMod;
        t[4] = (A[id4]+A[id4])%PlainMod;
        t[5] = (A[id5]+A[id5])%PlainMod;
        t[6] = (A[id6]+A[id6])%PlainMod;
        t[7] = (A[id7]+A[id7])%PlainMod;

A[id0] + A[id1] + t[2] + A[id3] + t[4] + A[id5] + t[6] + t[7];
t[0] + A[id1]  + t[3] + A[id4]+t[5] + A[id6] + t[7];
t[0] + t[1] + A[id2] + A[id3] + t[4]+A[id5] + t[6] + A[id7];
A[id0] + t[1] + t[2] + A[id3] + t[5] + A[id6] + t[7];
t[0] + A[id1] + t[2] + t[3] + A[id4]+A[id5] + t[6] + A[id7];
A[id0] + t[1] + A[id2] + t[3] + t[4]+A[id5] + t[7];
t[0] + A[id1] + t[2] + A[id3] + t[4]+t[5] + A[id6] + A[id7];
t[1] + A[id2] + t[3] + A[id4]+t[5] + t[6] + A[id7];
*/
void YusP::MC32_2(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp = A;
    std::array<int, 8> index = {0, 1, 2, 3, 4, 5, 6, 7};
    std::vector<uint64_t> t = A;
    for (int i = 0; i < 4; i++)
    {
        int s = 8 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s;
        t[id0] += t[id0];
        t[id1] += t[id1];
        t[id2] += t[id2];
        t[id3] += t[id3];
        t[id4] += t[id4];
        t[id5] += t[id5];
        t[id6] += t[id6];
        t[id7] += t[id7];

        A[id0] = (A[id0] + temp[id1] + t[id2] + temp[id3] + t[id4] + temp[id5] + t[id6] + t[id7]) % PlainMod;
        A[id1] = (t[id0] + A[id1] + t[id3] + temp[id4] + t[id5] + temp[id6] + t[id7]) % PlainMod;
        A[id2] = (t[id0] + t[id1] + A[id2] + temp[id3] + t[id4] + temp[id5] + t[id6] + temp[id7]) % PlainMod;
        A[id3] = (temp[id0] + t[id1] + t[id2] + A[id3] + t[id5] + temp[id6] + t[id7]) % PlainMod;
        A[id4] = (t[id0] + temp[id1] + t[id2] + t[id3] + A[id4] + temp[id5] + t[id6] + temp[id7]) % PlainMod;
        A[id5] = (temp[id0] + t[id1] + temp[id2] + t[id3] + t[id4] + A[id5] + t[id7]) % PlainMod;
        A[id6] = (t[id0] + temp[id1] + t[id2] + temp[id3] + t[id4] + t[id5] + A[id6] + temp[id7]) % PlainMod;
        A[id7] = (t[id1] + temp[id2] + t[id3] + temp[id4] + t[id5] + t[id6] + A[id7]) % PlainMod;
    }
}

// 行移位：A = A * constMR
void YusP::MR32_2(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp = A;
    std::array<int, 8> index = {0, 1, 8, 9, 16, 17, 24, 25};
    std::vector<uint64_t> t = A;
    for (int i = 0; i < 4; i++)
    {
        int s = 2 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s;
        t[id0] += t[id0];
        t[id1] += t[id1];
        t[id2] += t[id2];
        t[id3] += t[id3];
        t[id4] += t[id4];
        t[id5] += t[id5];
        t[id6] += t[id6];
        t[id7] += t[id7];

        A[id0] = (A[id0] + temp[id1] + t[id2] + temp[id3] + t[id4] + temp[id5] + t[id6] + t[id7]) % PlainMod;
        A[id1] = (t[id0] + A[id1] + t[id3] + temp[id4] + t[id5] + temp[id6] + t[id7]) % PlainMod;
        A[id2] = (t[id0] + t[id1] + A[id2] + temp[id3] + t[id4] + temp[id5] + t[id6] + temp[id7]) % PlainMod;
        A[id3] = (temp[id0] + t[id1] + t[id2] + A[id3] + t[id5] + temp[id6] + t[id7]) % PlainMod;
        A[id4] = (t[id0] + temp[id1] + t[id2] + t[id3] + A[id4] + temp[id5] + t[id6] + temp[id7]) % PlainMod;
        A[id5] = (temp[id0] + t[id1] + temp[id2] + t[id3] + t[id4] + A[id5] + t[id7]) % PlainMod;
        A[id6] = (t[id0] + temp[id1] + t[id2] + temp[id3] + t[id4] + t[id5] + A[id6] + temp[id7]) % PlainMod;
        A[id7] = (t[id1] + temp[id2] + t[id3] + temp[id4] + t[id5] + t[id6] + A[id7]) % PlainMod;
    }
}

void YusP::MC64(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp(64);
    std::array<int, 8> index = {0, 1, 2, 3, 4, 5, 6, 7};
    for (int i = 0; i < 8; i++)
    {
        int s = 8 * i;
        temp[index[0] + s] = ((A[index[1] + s] + A[index[2] + s] + A[index[4] + s] + A[index[5] + s] + A[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[1] + s] = ((A[index[0] + s] + A[index[3] + s] + A[index[4] + s] - A[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[2] + s] = ((A[index[1] + s] + A[index[2] + s] + A[index[3] + s] + A[index[4] + s] - A[index[6] + s] - A[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[3] + s] = ((A[index[0] + s] + A[index[1] + s] + A[index[3] + s] + A[index[5] + s] - A[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[4] + s] = ((A[index[0] + s] + A[index[1] + s] + A[index[2] + s] + A[index[4] + s] - A[index[6] + s] + A[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[5] + s] = ((A[index[0] + s] + A[index[2] + s] + A[index[3] + s] + A[index[5] + s] + A[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[6] + s] = ((A[index[1] + s] - A[index[0] + s] - A[index[2] + s] + A[index[3] + s] + A[index[4] + s] + A[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[7] + s] = ((A[index[2] + s] - A[index[0] + s] + A[index[3] + s] + A[index[4] + s] + A[index[5] + s] - A[index[6] + s] + A[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
    }
    std::copy(temp.begin(), temp.end(), A.begin());
}

// 行移位：A = A * constMR
void YusP::MR64(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp(64);
    std::array<int, 8> index = {0, 8, 16, 24, 32, 40, 48, 56};
    for (int i = 0; i < 8; i++)
    {
        int s = i;
        temp[index[0] + s] = ((A[index[1] + s] + A[index[2] + s] + A[index[4] + s] + A[index[5] + s] + A[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[1] + s] = ((A[index[0] + s] + A[index[3] + s] + A[index[4] + s] - A[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[2] + s] = ((A[index[1] + s] + A[index[2] + s] + A[index[3] + s] + A[index[4] + s] - A[index[6] + s] - A[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[3] + s] = ((A[index[0] + s] + A[index[1] + s] + A[index[3] + s] + A[index[5] + s] - A[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[4] + s] = ((A[index[0] + s] + A[index[1] + s] + A[index[2] + s] + A[index[4] + s] - A[index[6] + s] + A[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[5] + s] = ((A[index[0] + s] + A[index[2] + s] + A[index[3] + s] + A[index[5] + s] + A[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[6] + s] = ((A[index[1] + s] - A[index[0] + s] - A[index[2] + s] + A[index[3] + s] + A[index[4] + s] + A[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[7] + s] = ((A[index[2] + s] - A[index[0] + s] + A[index[3] + s] + A[index[4] + s] + A[index[5] + s] - A[index[6] + s] + A[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
    }
    std::copy(temp.begin(), temp.end(), A.begin());
}
void YusP::Sbox(vector<uint64_t> &A)
{
    uint64_t temp;
    for (int i = 0; i < A.size(); i += 2)
    {
        // 两次Feistel轮函数(x,y)——> (y+x^2,x)
        for (int j = 0; j < 2; j++)
        {
            temp = A[i];
            A[i] = ((A[i + 1] + A[i] * A[i]) % PlainMod + PlainMod) % PlainMod;
            A[i + 1] = temp;
        }
    }
}
void YusP::MC48_3(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp = A;
    std::array<int, 12> index = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    for (int i = 0; i < 4; i++)
    {
        int s = 12 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s, id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s;
        A[id0] = (((temp[id1] + temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id7] + temp[id9] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        A[id1] = (((temp[id0] + temp[id1] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id9]) % PlainMod) + PlainMod) % PlainMod;
        A[id2] = (((temp[id0] + temp[id1] + temp[id2] + temp[id4] + temp[id5] + temp[id6] + temp[id10]) % PlainMod) + PlainMod) % PlainMod;
        A[id3] = (((temp[id0] + temp[id2] + temp[id4] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id9] + temp[id10]) % PlainMod) + PlainMod) % PlainMod;
        A[id4] = (((temp[id2] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        A[id5] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id9]) % PlainMod) + PlainMod) % PlainMod;
        A[id6] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        A[id7] = (((temp[id0] + temp[id1] + temp[id2] + temp[id3] + temp[id6] + temp[id7] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        A[id8] = (((temp[id0] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        A[id9] = (((temp[id0] + temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id6] + temp[id8] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        A[id10] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id9] + temp[id10]) % PlainMod) + PlainMod) % PlainMod;
        A[id11] = (((temp[id1] + temp[id3] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::MR48_3(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp = A;
    std::array<int, 12> index = {0, 1, 2, 12, 13, 14, 24, 25, 26, 36, 37, 38};

    for (int i = 0; i < 4; i++)
    {
        int s = 3 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s, id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s;
        A[id0] = (((temp[id1] + temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id7] + temp[id9] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        A[id1] = (((temp[id0] + temp[id1] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id9]) % PlainMod) + PlainMod) % PlainMod;
        A[id2] = (((temp[id0] + temp[id1] + temp[id2] + temp[id4] + temp[id5] + temp[id6] + temp[id10]) % PlainMod) + PlainMod) % PlainMod;
        A[id3] = (((temp[id0] + temp[id2] + temp[id4] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id9] + temp[id10]) % PlainMod) + PlainMod) % PlainMod;
        A[id4] = (((temp[id2] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        A[id5] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id9]) % PlainMod) + PlainMod) % PlainMod;
        A[id6] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        A[id7] = (((temp[id0] + temp[id1] + temp[id2] + temp[id3] + temp[id6] + temp[id7] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        A[id8] = (((temp[id0] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        A[id9] = (((temp[id0] + temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id6] + temp[id8] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        A[id10] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id9] + temp[id10]) % PlainMod) + PlainMod) % PlainMod;
        A[id11] = (((temp[id1] + temp[id3] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::Sbox_3(vector<uint64_t> &A)
{
    vector<uint64_t> temp = A;
    // 一次迭代，(x1,x2,x3)——> (x1,x1*x3+x2,-x1*x2+x1*x3+x3)
    for (int i = 0; i < A.size(); i += 3)
    {
        A[i + 1] = (temp[i] * temp[i + 2] + temp[i + 1]) % PlainMod;
        A[i + 2] = (((-temp[i] * temp[i + 1] + temp[i] * temp[i + 2] + temp[i + 2]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::MC64_4(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp = A;
    std::array<int, 16> index = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    for (int i = 0; i < 4; i++)
    {
        int s = 16 * i;

        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s,
            id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s,
            id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s,
            id12 = index[12] + s, id13 = index[13] + s, id14 = index[14] + s, id15 = index[15] + s;

        A[id0] = (((temp[id0] + temp[id2] + temp[id4] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id13]) % PlainMod) + PlainMod) % PlainMod;
        A[id1] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id7] + temp[id8] + temp[id10] + temp[id11] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id2] = (((temp[id0] + temp[id1] + temp[id2] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id3] = (((temp[id1] + temp[id2] + temp[id6] + temp[id7] + temp[id8] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id4] = (((temp[id0] + temp[id4] + temp[id6] + temp[id8] + temp[id9] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id5] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id11] + temp[id12] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id6] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id6] + temp[id9] + temp[id11] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id7] = (((temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        A[id8] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        A[id9] = (((temp[id0] + temp[id2] + temp[id3] + temp[id7] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id10] = (((temp[id0] + temp[id2] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id13] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id11] = (((temp[id0] + temp[id6] + temp[id9] + temp[id10] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id12] = (((temp[id0] + temp[id1] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id13] = (((temp[id3] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id11] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id14] = (((temp[id1] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id10] + temp[id11] + temp[id12] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id15] = (((temp[id2] + temp[id4] + temp[id10] + temp[id11] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::MR64_4(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp = A;
    std::array<int, 16> index = {0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51};

    for (int i = 0; i < 4; i++)
    {
        int s = 4 * i;

        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s,
            id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s,
            id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s,
            id12 = index[12] + s, id13 = index[13] + s, id14 = index[14] + s, id15 = index[15] + s;

        A[id0] = (((temp[id0] + temp[id2] + temp[id4] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id13]) % PlainMod) + PlainMod) % PlainMod;
        A[id1] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id7] + temp[id8] + temp[id10] + temp[id11] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id2] = (((temp[id0] + temp[id1] + temp[id2] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id3] = (((temp[id1] + temp[id2] + temp[id6] + temp[id7] + temp[id8] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id4] = (((temp[id0] + temp[id4] + temp[id6] + temp[id8] + temp[id9] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id5] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id11] + temp[id12] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id6] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id6] + temp[id9] + temp[id11] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id7] = (((temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        A[id8] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        A[id9] = (((temp[id0] + temp[id2] + temp[id3] + temp[id7] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id10] = (((temp[id0] + temp[id2] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id13] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id11] = (((temp[id0] + temp[id6] + temp[id9] + temp[id10] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id12] = (((temp[id0] + temp[id1] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id13] = (((temp[id3] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id11] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id14] = (((temp[id1] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id10] + temp[id11] + temp[id12] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id15] = (((temp[id2] + temp[id4] + temp[id10] + temp[id11] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::Sbox_4(vector<uint64_t> &A)
{
    vector<uint64_t> temp = A;
    // 一次迭代，(x0,x1,x2,x3)——> (x0, x0^2 + x1, -x0*x3 -x1^2 + x2, x0*x2 -x0*x3 + x1^2 + x3)
    for (int i = 0; i < A.size(); i += 4)
    {
        A[i + 1] = (temp[i] * temp[i] + temp[i + 1]) % PlainMod;
        A[i + 2] = (((-temp[i] * temp[i + 3] - temp[i + 1] * temp[i + 1] + temp[i + 2]) % PlainMod) + PlainMod) % PlainMod;
        A[i + 3] = ((temp[i] * temp[i + 2] - temp[i] * temp[i + 3] + temp[i + 1] * temp[i + 1] + temp[i + 3]) % PlainMod + PlainMod) % PlainMod;
    }
}
void YusP::M36_5(vector<uint64_t> &A)
{

    std::array<int, 36> index = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35};
    int s = 0;
    int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s,
        id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s,
        id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s,
        id12 = index[12] + s, id13 = index[13] + s, id14 = index[14] + s, id15 = index[15] + s,
        id16 = index[16] + s, id17 = index[17] + s, id18 = index[18] + s, id19 = index[19] + s,
        id20 = index[20] + s, id21 = index[21] + s, id22 = index[22] + s, id23 = index[23] + s,
        id24 = index[24] + s, id25 = index[25] + s, id26 = index[26] + s, id27 = index[27] + s,
        id28 = index[28] + s, id29 = index[29] + s, id30 = index[30] + s, id31 = index[31] + s,
        id32 = index[32] + s, id33 = index[33] + s, id34 = index[34] + s, id35 = index[35] + s;
    std::vector<uint64_t> temp = A;
    for (int i = 0; i < 2; i++)
    {
        A[id0] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id6] + temp[id7] + temp[id9] + temp[id11] + temp[id13] + temp[id16] + temp[id17] + temp[id19] + temp[id20] + temp[id22] + temp[id25] + temp[id26] + temp[id27] + temp[id29] + temp[id30] + temp[id33] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id1] = (((temp[id1] + temp[id3] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id12] + temp[id17] + temp[id18] + temp[id19] + temp[id21] + temp[id22] + temp[id24] + temp[id27] + temp[id31] + temp[id33] + temp[id34] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id2] = (((temp[id0] + temp[id3] + temp[id7] + temp[id9] + temp[id10] + temp[id11] + temp[id13] + temp[id14] + temp[id15] + temp[id17] + temp[id18] + temp[id21] + temp[id22] + temp[id23] + temp[id26] + temp[id28] + temp[id32] + temp[id33] + temp[id34]) % PlainMod) + PlainMod) % PlainMod;
        A[id3] = (((temp[id0] + temp[id2] + temp[id4] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id12] + temp[id14] + temp[id16] + temp[id19] + temp[id20] + temp[id22] + temp[id23] + temp[id25] + temp[id28] + temp[id29] + temp[id30] + temp[id32] + temp[id33]) % PlainMod) + PlainMod) % PlainMod;
        A[id4] = (((temp[id0] + temp[id1] + temp[id2] + temp[id4] + temp[id6] + temp[id7] + temp[id9] + temp[id10] + temp[id11] + temp[id12] + temp[id13] + temp[id15] + temp[id20] + temp[id21] + temp[id22] + temp[id24] + temp[id25] + temp[id27] + temp[id30] + temp[id34]) % PlainMod) + PlainMod) % PlainMod;
        A[id5] = (((temp[id0] + temp[id1] + temp[id3] + temp[id6] + temp[id10] + temp[id12] + temp[id13] + temp[id14] + temp[id16] + temp[id17] + temp[id18] + temp[id20] + temp[id21] + temp[id24] + temp[id25] + temp[id26] + temp[id29] + temp[id31] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id6] = (((temp[id0] + temp[id3] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id11] + temp[id12] + temp[id13] + temp[id15] + temp[id17] + temp[id19] + temp[id22] + temp[id23] + temp[id25] + temp[id26] + temp[id28] + temp[id31] + temp[id32] + temp[id33] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id7] = (((temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id9] + temp[id10] + temp[id12] + temp[id13] + temp[id14] + temp[id15] + temp[id16] + temp[id18] + temp[id23] + temp[id24] + temp[id25] + temp[id27] + temp[id28] + temp[id30] + temp[id33]) % PlainMod) + PlainMod) % PlainMod;
        A[id8] = (((temp[id2] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id13] + temp[id15] + temp[id16] + temp[id17] + temp[id19] + temp[id20] + temp[id21] + temp[id23] + temp[id24] + temp[id27] + temp[id28] + temp[id29] + temp[id32] + temp[id34]) % PlainMod) + PlainMod) % PlainMod;
        A[id9] = (((temp[id0] + temp[id2] + temp[id3] + temp[id6] + temp[id8] + temp[id10] + temp[id11] + temp[id12] + temp[id13] + temp[id14] + temp[id15] + temp[id16] + temp[id18] + temp[id20] + temp[id22] + temp[id25] + temp[id26] + temp[id28] + temp[id29] + temp[id31] + temp[id34] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id10] = (((temp[id0] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id12] + temp[id13] + temp[id15] + temp[id16] + temp[id17] + temp[id18] + temp[id19] + temp[id21] + temp[id26] + temp[id27] + temp[id28] + temp[id30] + temp[id31] + temp[id33]) % PlainMod) + PlainMod) % PlainMod;
        A[id11] = (((temp[id1] + temp[id5] + temp[id6] + temp[id7] + temp[id9] + temp[id12] + temp[id16] + temp[id18] + temp[id19] + temp[id20] + temp[id22] + temp[id23] + temp[id24] + temp[id26] + temp[id27] + temp[id30] + temp[id31] + temp[id32] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id12] = (((temp[id1] + temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id9] + temp[id11] + temp[id13] + temp[id14] + temp[id15] + temp[id16] + temp[id17] + temp[id18] + temp[id19] + temp[id21] + temp[id23] + temp[id25] + temp[id28] + temp[id29] + temp[id31] + temp[id32] + temp[id34]) % PlainMod) + PlainMod) % PlainMod;
        A[id13] = (((temp[id0] + temp[id3] + temp[id7] + temp[id9] + temp[id10] + temp[id11] + temp[id13] + temp[id15] + temp[id16] + temp[id18] + temp[id19] + temp[id20] + temp[id21] + temp[id22] + temp[id24] + temp[id29] + temp[id30] + temp[id31] + temp[id33] + temp[id34]) % PlainMod) + PlainMod) % PlainMod;
        A[id14] = (((temp[id2] + temp[id4] + temp[id8] + temp[id9] + temp[id10] + temp[id12] + temp[id15] + temp[id19] + temp[id21] + temp[id22] + temp[id23] + temp[id25] + temp[id26] + temp[id27] + temp[id29] + temp[id30] + temp[id33] + temp[id34] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id15] = (((temp[id1] + temp[id4] + temp[id5] + temp[id6] + temp[id8] + temp[id9] + temp[id12] + temp[id14] + temp[id16] + temp[id17] + temp[id18] + temp[id19] + temp[id20] + temp[id21] + temp[id22] + temp[id24] + temp[id26] + temp[id28] + temp[id31] + temp[id32] + temp[id34] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id16] = (((temp[id0] + temp[id1] + temp[id3] + temp[id6] + temp[id10] + temp[id12] + temp[id13] + temp[id14] + temp[id16] + temp[id18] + temp[id19] + temp[id21] + temp[id22] + temp[id23] + temp[id24] + temp[id25] + temp[id27] + temp[id32] + temp[id33] + temp[id34]) % PlainMod) + PlainMod) % PlainMod;
        A[id17] = (((temp[id0] + temp[id1] + temp[id2] + temp[id5] + temp[id7] + temp[id11] + temp[id12] + temp[id13] + temp[id15] + temp[id18] + temp[id22] + temp[id24] + temp[id25] + temp[id26] + temp[id28] + temp[id29] + temp[id30] + temp[id32] + temp[id33]) % PlainMod) + PlainMod) % PlainMod;
        A[id18] = (((temp[id1] + temp[id2] + temp[id4] + temp[id7] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id15] + temp[id17] + temp[id19] + temp[id20] + temp[id21] + temp[id22] + temp[id23] + temp[id24] + temp[id25] + temp[id27] + temp[id29] + temp[id31] + temp[id34] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id19] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id13] + temp[id15] + temp[id16] + temp[id17] + temp[id19] + temp[id21] + temp[id22] + temp[id24] + temp[id25] + temp[id26] + temp[id27] + temp[id28] + temp[id30] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id20] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id10] + temp[id14] + temp[id15] + temp[id16] + temp[id18] + temp[id21] + temp[id25] + temp[id27] + temp[id28] + temp[id29] + temp[id31] + temp[id32] + temp[id33] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id21] = (((temp[id1] + temp[id2] + temp[id4] + temp[id5] + temp[id7] + temp[id10] + temp[id11] + temp[id12] + temp[id14] + temp[id15] + temp[id18] + temp[id20] + temp[id22] + temp[id23] + temp[id24] + temp[id25] + temp[id26] + temp[id27] + temp[id28] + temp[id30] + temp[id32] + temp[id34]) % PlainMod) + PlainMod) % PlainMod;
        A[id22] = (((temp[id2] + temp[id3] + temp[id4] + temp[id6] + temp[id7] + temp[id9] + temp[id12] + temp[id16] + temp[id18] + temp[id19] + temp[id20] + temp[id22] + temp[id24] + temp[id25] + temp[id27] + temp[id28] + temp[id29] + temp[id30] + temp[id31] + temp[id33]) % PlainMod) + PlainMod) % PlainMod;
        A[id23] = (((temp[id0] + temp[id2] + temp[id3] + temp[id6] + temp[id7] + temp[id8] + temp[id11] + temp[id13] + temp[id17] + temp[id18] + temp[id19] + temp[id21] + temp[id24] + temp[id28] + temp[id30] + temp[id31] + temp[id32] + temp[id34] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id24] = (((temp[id1] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id10] + temp[id13] + temp[id14] + temp[id15] + temp[id17] + temp[id18] + temp[id21] + temp[id23] + temp[id25] + temp[id26] + temp[id27] + temp[id28] + temp[id29] + temp[id30] + temp[id31] + temp[id33] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id25] = (((temp[id0] + temp[id5] + temp[id6] + temp[id7] + temp[id9] + temp[id10] + temp[id12] + temp[id15] + temp[id19] + temp[id21] + temp[id22] + temp[id23] + temp[id25] + temp[id27] + temp[id28] + temp[id30] + temp[id31] + temp[id32] + temp[id33] + temp[id34]) % PlainMod) + PlainMod) % PlainMod;
        A[id26] = (((temp[id1] + temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id9] + temp[id10] + temp[id11] + temp[id14] + temp[id16] + temp[id20] + temp[id21] + temp[id22] + temp[id24] + temp[id27] + temp[id31] + temp[id33] + temp[id34] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id27] = (((temp[id0] + temp[id2] + temp[id4] + temp[id7] + temp[id8] + temp[id10] + temp[id11] + temp[id13] + temp[id16] + temp[id17] + temp[id18] + temp[id20] + temp[id21] + temp[id24] + temp[id26] + temp[id28] + temp[id29] + temp[id30] + temp[id31] + temp[id32] + temp[id33] + temp[id34]) % PlainMod) + PlainMod) % PlainMod;
        A[id28] = (((temp[id0] + temp[id1] + temp[id3] + temp[id8] + temp[id9] + temp[id10] + temp[id12] + temp[id13] + temp[id15] + temp[id18] + temp[id22] + temp[id24] + temp[id25] + temp[id26] + temp[id28] + temp[id30] + temp[id31] + temp[id33] + temp[id34] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id29] = (((temp[id0] + temp[id1] + temp[id2] + temp[id4] + temp[id5] + temp[id6] + temp[id8] + temp[id9] + temp[id12] + temp[id13] + temp[id14] + temp[id17] + temp[id19] + temp[id23] + temp[id24] + temp[id25] + temp[id27] + temp[id30] + temp[id34]) % PlainMod) + PlainMod) % PlainMod;
        A[id30] = (((temp[id0] + temp[id1] + temp[id3] + temp[id5] + temp[id7] + temp[id10] + temp[id11] + temp[id13] + temp[id14] + temp[id16] + temp[id19] + temp[id20] + temp[id21] + temp[id23] + temp[id24] + temp[id27] + temp[id29] + temp[id31] + temp[id32] + temp[id33] + temp[id34] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id31] = (((temp[id0] + temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id6] + temp[id11] + temp[id12] + temp[id13] + temp[id15] + temp[id16] + temp[id18] + temp[id21] + temp[id25] + temp[id27] + temp[id28] + temp[id29] + temp[id31] + temp[id33] + temp[id34]) % PlainMod) + PlainMod) % PlainMod;
        A[id32] = (((temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id15] + temp[id16] + temp[id17] + temp[id20] + temp[id22] + temp[id26] + temp[id27] + temp[id28] + temp[id30] + temp[id33]) % PlainMod) + PlainMod) % PlainMod;
        A[id33] = (((temp[id0] + temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id6] + temp[id8] + temp[id10] + temp[id13] + temp[id14] + temp[id16] + temp[id17] + temp[id19] + temp[id22] + temp[id23] + temp[id24] + temp[id26] + temp[id27] + temp[id30] + temp[id32] + temp[id34] + temp[id35]) % PlainMod) + PlainMod) % PlainMod;
        A[id34] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id6] + temp[id7] + temp[id9] + temp[id14] + temp[id15] + temp[id16] + temp[id18] + temp[id19] + temp[id21] + temp[id24] + temp[id28] + temp[id30] + temp[id31] + temp[id32] + temp[id34]) % PlainMod) + PlainMod) % PlainMod;
        A[id35] = (((temp[id0] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id11] + temp[id12] + temp[id14] + temp[id15] + temp[id18] + temp[id19] + temp[id20] + temp[id23] + temp[id25] + temp[id29] + temp[id30] + temp[id31] + temp[id33]) % PlainMod) + PlainMod) % PlainMod;
        temp = A;
    }
}
void YusP::Sbox_5(vector<uint64_t> &A)
{
    vector<uint64_t> temp = A;
    // 一次迭代，(x0,x1,x2)——>(x0,x0*x2+x1,-x0*x1+x0*x2+x2)
    for (int i = 0; i < A.size(); i += 3)
    {
        A[i + 1] = (temp[i] * temp[i + 2] + temp[i + 1]) % PlainMod;
        A[i + 2] = (((-temp[i] * temp[i + 1] + temp[i] * temp[i + 2] + temp[i + 2]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::MC64_6(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp = A;
    std::array<int, 16> index = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    for (int i = 0; i < 4; i++)
    {
        int s = 16 * i;

        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s,
            id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s,
            id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s,
            id12 = index[12] + s, id13 = index[13] + s, id14 = index[14] + s, id15 = index[15] + s;

        A[id0] = (((temp[id0] + temp[id2] + temp[id4] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id13]) % PlainMod) + PlainMod) % PlainMod;
        A[id1] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id7] + temp[id8] + temp[id10] + temp[id11] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id2] = (((temp[id0] + temp[id1] + temp[id2] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id3] = (((temp[id1] + temp[id2] + temp[id6] + temp[id7] + temp[id8] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id4] = (((temp[id0] + temp[id4] + temp[id6] + temp[id8] + temp[id9] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id5] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id11] + temp[id12] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id6] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id6] + temp[id9] + temp[id11] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id7] = (((temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        A[id8] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        A[id9] = (((temp[id0] + temp[id2] + temp[id3] + temp[id7] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id10] = (((temp[id0] + temp[id2] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id13] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id11] = (((temp[id0] + temp[id6] + temp[id9] + temp[id10] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id12] = (((temp[id0] + temp[id1] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id13] = (((temp[id3] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id11] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id14] = (((temp[id1] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id10] + temp[id11] + temp[id12] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id15] = (((temp[id2] + temp[id4] + temp[id10] + temp[id11] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::MR64_6(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp = A;
    std::array<int, 16> index = {0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51};

    for (int i = 0; i < 4; i++)
    {
        int s = 4 * i;

        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s,
            id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s,
            id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s,
            id12 = index[12] + s, id13 = index[13] + s, id14 = index[14] + s, id15 = index[15] + s;

        A[id0] = (((temp[id0] + temp[id2] + temp[id4] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id13]) % PlainMod) + PlainMod) % PlainMod;
        A[id1] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id7] + temp[id8] + temp[id10] + temp[id11] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id2] = (((temp[id0] + temp[id1] + temp[id2] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id3] = (((temp[id1] + temp[id2] + temp[id6] + temp[id7] + temp[id8] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id4] = (((temp[id0] + temp[id4] + temp[id6] + temp[id8] + temp[id9] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id5] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id11] + temp[id12] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id6] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id6] + temp[id9] + temp[id11] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id7] = (((temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        A[id8] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        A[id9] = (((temp[id0] + temp[id2] + temp[id3] + temp[id7] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id10] = (((temp[id0] + temp[id2] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id13] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id11] = (((temp[id0] + temp[id6] + temp[id9] + temp[id10] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id12] = (((temp[id0] + temp[id1] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id13] = (((temp[id3] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id11] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        A[id14] = (((temp[id1] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id10] + temp[id11] + temp[id12] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        A[id15] = (((temp[id2] + temp[id4] + temp[id10] + temp[id11] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::Sbox_6(vector<uint64_t> &A)
{
    vector<uint64_t> temp = A;
    // 一次迭代，(x0,x1,x2,x3)——> (x0, -x0*x1+x1+x0*x2+x0*x3, x1+x0*x2-x2+x0*x3, -x0*x1-x2-x0*x3+x3)

    for (int i = 0; i < A.size(); i += 4)
    {
        uint64_t t01 = temp[i] * temp[i + 1];
        uint64_t t02 = temp[i] * temp[i + 2];
        uint64_t t03 = temp[i] * temp[i + 3];
        uint64_t t02_03_1 = t02 + t03 + temp[i + 1];

        A[i + 1] = (((t02_03_1 - t01) % PlainMod) + PlainMod) % PlainMod;
        A[i + 2] = (((t02_03_1 - temp[i + 2]) % PlainMod) + PlainMod) % PlainMod;
        A[i + 3] = (((temp[i + 3] - t01 - temp[i + 2] - t03) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::M36_7(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp = A;
    A[0] = (((2 * temp[0] + 2 * temp[1] + temp[2] + 2 * temp[3] + temp[4] + temp[5] + temp[8] + temp[9] + temp[10] + temp[11] + temp[14] + 2 * temp[15] + temp[16] + 2 * temp[17] + 2 * temp[18] + 2 * temp[22] + 2 * temp[24] + temp[26] + 2 * temp[27] + temp[28] + temp[29] + temp[30] + temp[31] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
    A[1] = (((temp[0] + temp[2] + temp[3] + temp[4] + 2 * temp[6] + 2 * temp[7] + 2 * temp[9] + temp[10] + temp[11] + 2 * temp[13] + temp[14] + temp[15] + temp[17] + 2 * temp[20] + 2 * temp[21] + temp[22] + temp[23] + temp[24] + 2 * temp[25] + temp[26] + 2 * temp[27] + temp[29] + 2 * temp[30] + temp[32] + 2 * temp[33] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[2] = (((temp[1] + temp[2] + temp[3] + temp[6] + 2 * temp[8] + 2 * temp[10] + 2 * temp[11] + temp[12] + temp[16] + temp[17] + 2 * temp[20] + temp[23] + temp[26] + temp[27] + temp[29] + temp[31] + 2 * temp[32]) % PlainMod) + PlainMod) % PlainMod;
    A[3] = (((temp[1] + temp[5] + temp[6] + temp[7] + temp[9] + temp[10] + 2 * temp[12] + temp[14] + temp[15] + 2 * temp[17] + 2 * temp[19] + temp[20] + temp[21] + temp[23] + 2 * temp[25] + 2 * temp[26] + temp[28] + temp[29] + temp[30] + 2 * temp[31] + temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[4] = (((temp[1] + 2 * temp[4] + 2 * temp[5] + temp[6] + 2 * temp[7] + temp[8] + temp[9] + temp[12] + temp[13] + temp[14] + temp[15] + temp[18] + 2 * temp[19] + temp[20] + 2 * temp[21] + 2 * temp[22] + 2 * temp[26] + 2 * temp[28] + temp[30] + 2 * temp[31] + temp[32] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[5] = (((temp[0] + 2 * temp[1] + temp[3] + temp[4] + temp[6] + temp[7] + temp[8] + 2 * temp[10] + 2 * temp[11] + 2 * temp[13] + temp[14] + temp[15] + 2 * temp[17] + temp[18] + temp[19] + temp[21] + 2 * temp[24] + 2 * temp[25] + temp[26] + temp[27] + temp[28] + 2 * temp[29] + temp[30] + 2 * temp[31] + temp[33] + 2 * temp[34]) % PlainMod) + PlainMod) % PlainMod;
    A[6] = (((2 * temp[0] + temp[5] + temp[6] + temp[7] + temp[10] + 2 * temp[12] + 2 * temp[14] + 2 * temp[15] + temp[16] + temp[20] + temp[21] + 2 * temp[24] + temp[27] + temp[30] + temp[31] + temp[33] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[7] = (((temp[2] + 2 * temp[3] + temp[5] + temp[9] + temp[10] + temp[11] + temp[13] + temp[14] + 2 * temp[16] + temp[18] + temp[19] + 2 * temp[21] + 2 * temp[23] + temp[24] + temp[25] + temp[27] + 2 * temp[29] + 2 * temp[30] + temp[32] + temp[33] + temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[8] = (((temp[0] + temp[1] + temp[2] + temp[3] + temp[5] + 2 * temp[8] + 2 * temp[9] + temp[10] + 2 * temp[11] + temp[12] + temp[13] + temp[16] + temp[17] + temp[18] + temp[19] + temp[22] + 2 * temp[23] + temp[24] + 2 * temp[25] + 2 * temp[26] + 2 * temp[30] + 2 * temp[32] + temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[9] = (((temp[1] + 2 * temp[2] + temp[4] + 2 * temp[5] + temp[7] + temp[8] + temp[10] + temp[11] + temp[12] + 2 * temp[14] + 2 * temp[15] + 2 * temp[17] + temp[18] + temp[19] + 2 * temp[21] + temp[22] + temp[23] + temp[25] + 2 * temp[28] + 2 * temp[29] + temp[30] + temp[31] + temp[32] + 2 * temp[33] + temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[10] = (((temp[1] + temp[3] + 2 * temp[4] + temp[9] + temp[10] + temp[11] + temp[14] + 2 * temp[16] + 2 * temp[18] + 2 * temp[19] + temp[20] + temp[24] + temp[25] + 2 * temp[28] + temp[31] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[11] = (((temp[0] + temp[1] + temp[2] + 2 * temp[3] + temp[6] + 2 * temp[7] + temp[9] + temp[13] + temp[14] + temp[15] + temp[17] + temp[18] + 2 * temp[20] + temp[22] + temp[23] + 2 * temp[25] + 2 * temp[27] + temp[28] + temp[29] + temp[31] + 2 * temp[33] + 2 * temp[34]) % PlainMod) + PlainMod) % PlainMod;
    A[12] = (((2 * temp[0] + temp[2] + 2 * temp[3] + temp[4] + temp[5] + temp[6] + temp[7] + temp[9] + 2 * temp[12] + 2 * temp[13] + temp[14] + 2 * temp[15] + temp[16] + temp[17] + temp[20] + temp[21] + temp[22] + temp[23] + temp[26] + 2 * temp[27] + temp[28] + 2 * temp[29] + 2 * temp[30] + 2 * temp[34]) % PlainMod) + PlainMod) % PlainMod;
    A[13] = (((temp[0] + 2 * temp[1] + temp[2] + 2 * temp[3] + temp[5] + 2 * temp[6] + temp[8] + 2 * temp[9] + temp[11] + temp[12] + temp[14] + temp[15] + temp[16] + 2 * temp[18] + 2 * temp[19] + 2 * temp[21] + temp[22] + temp[23] + 2 * temp[25] + temp[26] + temp[27] + temp[29] + 2 * temp[32] + 2 * temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[14] = (((temp[2] + temp[3] + temp[5] + temp[7] + 2 * temp[8] + temp[13] + temp[14] + temp[15] + temp[18] + 2 * temp[20] + 2 * temp[22] + 2 * temp[23] + temp[24] + temp[28] + temp[29] + 2 * temp[32] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[15] = (((2 * temp[1] + 2 * temp[2] + temp[4] + temp[5] + temp[6] + 2 * temp[7] + temp[10] + 2 * temp[11] + temp[13] + temp[17] + temp[18] + temp[19] + temp[21] + temp[22] + 2 * temp[24] + temp[26] + temp[27] + 2 * temp[29] + 2 * temp[31] + temp[32] + temp[33] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[16] = (((2 * temp[2] + 2 * temp[4] + temp[6] + 2 * temp[7] + temp[8] + temp[9] + temp[10] + temp[11] + temp[13] + 2 * temp[16] + 2 * temp[17] + temp[18] + 2 * temp[19] + temp[20] + temp[21] + temp[24] + temp[25] + temp[26] + temp[27] + temp[30] + 2 * temp[31] + temp[32] + 2 * temp[33] + 2 * temp[34]) % PlainMod) + PlainMod) % PlainMod;
    A[17] = (((2 * temp[0] + 2 * temp[1] + temp[2] + temp[3] + temp[4] + 2 * temp[5] + temp[6] + 2 * temp[7] + temp[9] + 2 * temp[10] + temp[12] + 2 * temp[13] + temp[15] + temp[16] + temp[18] + temp[19] + temp[20] + 2 * temp[22] + 2 * temp[23] + 2 * temp[25] + temp[26] + temp[27] + 2 * temp[29] + temp[30] + temp[31] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
    A[18] = (((2 * temp[0] + temp[3] + temp[6] + temp[7] + temp[9] + temp[11] + 2 * temp[12] + temp[17] + temp[18] + temp[19] + temp[22] + 2 * temp[24] + 2 * temp[26] + 2 * temp[27] + temp[28] + temp[32] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
    A[19] = (((temp[0] + temp[1] + temp[3] + 2 * temp[5] + 2 * temp[6] + temp[8] + temp[9] + temp[10] + 2 * temp[11] + temp[14] + 2 * temp[15] + temp[17] + temp[21] + temp[22] + temp[23] + temp[25] + temp[26] + 2 * temp[28] + temp[30] + temp[31] + 2 * temp[33] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[20] = (((temp[0] + 2 * temp[1] + 2 * temp[2] + 2 * temp[6] + 2 * temp[8] + temp[10] + 2 * temp[11] + temp[12] + temp[13] + temp[14] + temp[15] + temp[17] + 2 * temp[20] + 2 * temp[21] + temp[22] + 2 * temp[23] + temp[24] + temp[25] + temp[28] + temp[29] + temp[30] + temp[31] + temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[21] = (((temp[1] + 2 * temp[4] + 2 * temp[5] + temp[6] + temp[7] + temp[8] + 2 * temp[9] + temp[10] + 2 * temp[11] + temp[13] + 2 * temp[14] + temp[16] + 2 * temp[17] + temp[19] + temp[20] + temp[22] + temp[23] + temp[24] + 2 * temp[26] + 2 * temp[27] + 2 * temp[29] + temp[30] + temp[31] + 2 * temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[22] = (((temp[0] + temp[1] + 2 * temp[4] + temp[7] + temp[10] + temp[11] + temp[13] + temp[15] + 2 * temp[16] + temp[21] + temp[22] + temp[23] + temp[26] + 2 * temp[28] + 2 * temp[30] + 2 * temp[31] + temp[32]) % PlainMod) + PlainMod) % PlainMod;
    A[23] = (((2 * temp[1] + 2 * temp[3] + temp[4] + temp[5] + temp[7] + 2 * temp[9] + 2 * temp[10] + temp[12] + temp[13] + temp[14] + 2 * temp[15] + temp[18] + 2 * temp[19] + temp[21] + temp[25] + temp[26] + temp[27] + temp[29] + temp[30] + 2 * temp[32] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[24] = (((temp[2] + 2 * temp[3] + temp[4] + 2 * temp[5] + 2 * temp[6] + 2 * temp[10] + 2 * temp[12] + temp[14] + 2 * temp[15] + temp[16] + temp[17] + temp[18] + temp[19] + temp[21] + 2 * temp[24] + 2 * temp[25] + temp[26] + 2 * temp[27] + temp[28] + temp[29] + temp[32] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[25] = (((2 * temp[1] + temp[2] + temp[3] + temp[5] + 2 * temp[8] + 2 * temp[9] + temp[10] + temp[11] + temp[12] + 2 * temp[13] + temp[14] + 2 * temp[15] + temp[17] + 2 * temp[18] + temp[20] + 2 * temp[21] + temp[23] + temp[24] + temp[26] + temp[27] + temp[28] + 2 * temp[30] + 2 * temp[31] + 2 * temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[26] = (((temp[0] + temp[4] + temp[5] + 2 * temp[8] + temp[11] + temp[14] + temp[15] + temp[17] + temp[19] + 2 * temp[20] + temp[25] + temp[26] + temp[27] + temp[30] + 2 * temp[32] + 2 * temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[27] = (((2 * temp[0] + temp[2] + temp[3] + 2 * temp[5] + 2 * temp[7] + temp[8] + temp[9] + temp[11] + 2 * temp[13] + 2 * temp[14] + temp[16] + temp[17] + temp[18] + 2 * temp[19] + temp[22] + 2 * temp[23] + temp[25] + temp[29] + temp[30] + temp[31] + temp[33] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
    A[28] = (((temp[0] + temp[1] + temp[2] + temp[3] + temp[6] + 2 * temp[7] + temp[8] + 2 * temp[9] + 2 * temp[10] + 2 * temp[14] + 2 * temp[16] + temp[18] + 2 * temp[19] + temp[20] + temp[21] + temp[22] + temp[23] + temp[25] + 2 * temp[28] + 2 * temp[29] + temp[30] + 2 * temp[31] + temp[32] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
    A[29] = (((2 * temp[1] + temp[2] + temp[3] + 2 * temp[5] + temp[6] + temp[7] + temp[9] + 2 * temp[12] + 2 * temp[13] + temp[14] + temp[15] + temp[16] + 2 * temp[17] + temp[18] + 2 * temp[19] + temp[21] + 2 * temp[22] + temp[24] + 2 * temp[25] + temp[27] + temp[28] + temp[30] + temp[31] + temp[32] + 2 * temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[30] = (((2 * temp[0] + 2 * temp[2] + 2 * temp[3] + temp[4] + temp[8] + temp[9] + 2 * temp[12] + temp[15] + temp[18] + temp[19] + temp[21] + temp[23] + 2 * temp[24] + temp[29] + temp[30] + temp[31] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
    A[31] = (((temp[1] + temp[2] + 2 * temp[4] + temp[6] + temp[7] + 2 * temp[9] + 2 * temp[11] + temp[12] + temp[13] + temp[15] + 2 * temp[17] + 2 * temp[18] + temp[20] + temp[21] + temp[22] + 2 * temp[23] + temp[26] + 2 * temp[27] + temp[29] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[32] = (((temp[0] + temp[1] + temp[4] + temp[5] + temp[6] + temp[7] + temp[10] + 2 * temp[11] + temp[12] + 2 * temp[13] + 2 * temp[14] + 2 * temp[18] + 2 * temp[20] + temp[22] + 2 * temp[23] + temp[24] + temp[25] + temp[26] + temp[27] + temp[29] + 2 * temp[32] + 2 * temp[33] + temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[33] = (((temp[0] + 2 * temp[2] + 2 * temp[3] + 2 * temp[5] + temp[6] + temp[7] + 2 * temp[9] + temp[10] + temp[11] + temp[13] + 2 * temp[16] + 2 * temp[17] + temp[18] + temp[19] + temp[20] + 2 * temp[21] + temp[22] + 2 * temp[23] + temp[25] + 2 * temp[26] + temp[28] + 2 * temp[29] + temp[31] + temp[32] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[34] = (((temp[2] + 2 * temp[4] + 2 * temp[6] + 2 * temp[7] + temp[8] + temp[12] + temp[13] + 2 * temp[16] + temp[19] + temp[22] + temp[23] + temp[25] + temp[27] + 2 * temp[28] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    A[35] = (((temp[1] + temp[2] + temp[3] + temp[5] + temp[6] + 2 * temp[8] + temp[10] + temp[11] + 2 * temp[13] + 2 * temp[15] + temp[16] + temp[17] + temp[19] + 2 * temp[21] + 2 * temp[22] + temp[24] + temp[25] + temp[26] + 2 * temp[27] + temp[30] + 2 * temp[31] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
}
void YusP::Sbox_7(vector<uint64_t> &A)
{
    vector<uint64_t> temp = A;
    // 一次迭代，(x0,x1,x2,x3)——> (x0, -x0*x1+x1+x0*x2+x0*x3, x1+x0*x2-x2+x0*x3, -x0*x1-x2-x0*x3+x3)

    for (int i = 0; i < A.size(); i += 4)
    {
        uint64_t t01 = temp[i] * temp[i + 1];
        uint64_t t02 = temp[i] * temp[i + 2];
        uint64_t t03 = temp[i] * temp[i + 3];
        uint64_t t02_03_1 = t02 + t03 + temp[i + 1];

        A[i + 1] = (((t02_03_1 - t01) % PlainMod) + PlainMod) % PlainMod;
        A[i + 2] = (((t02_03_1 - temp[i + 2]) % PlainMod) + PlainMod) % PlainMod;
        A[i + 3] = (((temp[i + 3] - t01 - temp[i + 2] - t03) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::M18_8(vector<uint64_t> &A)
{
    std::vector<uint64_t> temp = A;
    A[0] = (((temp[2] + 2 * temp[3] + temp[4] + temp[5] + 2 * temp[6] + 2 * temp[8] + 2 * temp[9] + temp[11] + temp[13] + 2 * temp[14] + 2 * temp[15] + temp[16]) % PlainMod) + PlainMod) % PlainMod;
    A[1] = (((2 * temp[1] + temp[2] + temp[3] + 2 * temp[5] + 2 * temp[6] + 2 * temp[7] + temp[9] + temp[10] + 2 * temp[11] + 2 * temp[12] + 2 * temp[13] + 2 * temp[14] + temp[15] + 2 * temp[16]) % PlainMod) + PlainMod) % PlainMod;
    A[2] = (((2 * temp[0] + temp[1] + temp[2] + 2 * temp[4] + temp[5] + 2 * temp[6] + 2 * temp[7] + temp[8] + temp[9] + temp[11] + temp[12] + 2 * temp[13] + 2 * temp[15] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[3] = (((2 * temp[0] + temp[1] + temp[5] + 2 * temp[6] + temp[7] + temp[8] + 2 * temp[9] + 2 * temp[11] + 2 * temp[12] + temp[14] + temp[16] + 2 * temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[4] = (((temp[0] + 2 * temp[1] + 2 * temp[4] + temp[5] + temp[6] + 2 * temp[8] + 2 * temp[9] + 2 * temp[10] + temp[12] + temp[13] + 2 * temp[14] + 2 * temp[15] + 2 * temp[16] + 2 * temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[5] = (((2 * temp[0] + temp[2] + 2 * temp[3] + temp[4] + temp[5] + 2 * temp[7] + temp[8] + 2 * temp[9] + 2 * temp[10] + temp[11] + temp[12] + temp[14] + temp[15] + 2 * temp[16]) % PlainMod) + PlainMod) % PlainMod;
    A[6] = (((temp[1] + 2 * temp[2] + 2 * temp[3] + temp[4] + temp[8] + 2 * temp[9] + temp[10] + temp[11] + 2 * temp[12] + 2 * temp[14] + 2 * temp[15] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[7] = (((2 * temp[0] + 2 * temp[1] + 2 * temp[2] + temp[3] + 2 * temp[4] + 2 * temp[7] + temp[8] + temp[9] + 2 * temp[11] + 2 * temp[12] + 2 * temp[13] + temp[15] + temp[16] + 2 * temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[8] = (((temp[0] + 2 * temp[1] + 2 * temp[3] + temp[5] + 2 * temp[6] + temp[7] + temp[8] + 2 * temp[10] + temp[11] + 2 * temp[12] + 2 * temp[13] + temp[14] + temp[15] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[9] = (((2 * temp[0] + temp[2] + temp[4] + 2 * temp[5] + 2 * temp[6] + temp[7] + temp[11] + 2 * temp[12] + temp[13] + temp[14] + 2 * temp[15] + 2 * temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[10] = (((temp[0] + temp[1] + 2 * temp[2] + 2 * temp[3] + 2 * temp[4] + 2 * temp[5] + temp[6] + 2 * temp[7] + 2 * temp[10] + temp[11] + temp[12] + 2 * temp[14] + 2 * temp[15] + 2 * temp[16]) % PlainMod) + PlainMod) % PlainMod;
    A[11] = (((temp[0] + temp[2] + temp[3] + 2 * temp[4] + 2 * temp[6] + temp[8] + 2 * temp[9] + temp[10] + temp[11] + 2 * temp[13] + temp[14] + 2 * temp[15] + 2 * temp[16] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[12] = (((2 * temp[0] + 2 * temp[2] + 2 * temp[3] + temp[5] + temp[7] + 2 * temp[8] + 2 * temp[9] + temp[10] + temp[14] + 2 * temp[15] + temp[16] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[13] = (((2 * temp[0] + 2 * temp[1] + temp[3] + temp[4] + 2 * temp[5] + 2 * temp[6] + 2 * temp[7] + 2 * temp[8] + temp[9] + 2 * temp[10] + 2 * temp[13] + temp[14] + temp[15] + 2 * temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[14] = (((2 * temp[0] + 2 * temp[1] + temp[2] + temp[3] + temp[5] + temp[6] + 2 * temp[7] + 2 * temp[9] + temp[11] + 2 * temp[12] + temp[13] + temp[14] + 2 * temp[16] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[15] = (((2 * temp[0] + temp[1] + temp[2] + 2 * temp[3] + 2 * temp[5] + 2 * temp[6] + temp[8] + temp[10] + 2 * temp[11] + 2 * temp[12] + temp[13] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[16] = (((temp[0] + 2 * temp[2] + 2 * temp[3] + 2 * temp[4] + temp[6] + temp[7] + 2 * temp[8] + 2 * temp[9] + 2 * temp[10] + 2 * temp[11] + temp[12] + 2 * temp[13] + 2 * temp[16] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    A[17] = (((2 * temp[1] + temp[2] + 2 * temp[3] + 2 * temp[4] + temp[5] + temp[6] + temp[8] + temp[9] + 2 * temp[10] + 2 * temp[12] + temp[14] + 2 * temp[15] + temp[16] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
}
void YusP::Sbox_8(vector<uint64_t> &A)
{
    vector<uint64_t> temp = A;
    // 一次迭代，(x0,x1,x2)——>(x0,x0*x2+x1,-x0*x1+x0*x2+x2)
    for (int i = 0; i < A.size(); i += 3)
    {
        A[i + 1] = (temp[i] * temp[i + 2] + temp[i + 1]) % PlainMod;
        A[i + 2] = (((-temp[i] * temp[i + 1] + temp[i] * temp[i + 2] + temp[i + 2]) % PlainMod) + PlainMod) % PlainMod;
    }
}