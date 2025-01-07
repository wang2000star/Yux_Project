

#include "Yus_p.hpp"

using namespace std;

// static long PlainMod = 257;
// static long PlainMod = 65537;

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

// 列混淆：A = constMC * data,A是大小为32的一维数组，按照列优先看作8x4的矩阵
void YusP::MC32(vector<long> &data)
{
    std::vector<long> temp(32);
    std::array<int, 8> index = {0, 1, 2, 3, 4, 5, 6, 7};
    for (int i = 0; i < 4; i++)
    {
        int s = 8 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s;
        temp[id0] = ((data[id1] + data[id2] + data[id4] + data[id5] + data[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id1] = ((data[id0] + data[id3] + data[id4] - data[id7]) % PlainMod + PlainMod) % PlainMod;
        temp[id2] = ((data[id1] + data[id2] + data[id3] + data[id4] - data[id6] - data[id7]) % PlainMod + PlainMod) % PlainMod;
        temp[id3] = ((data[id0] + data[id1] + data[id3] + data[id5] - data[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id4] = ((data[id0] + data[id1] + data[id2] + data[id4] - data[id6] + data[id7]) % PlainMod + PlainMod) % PlainMod;
        temp[id5] = ((data[id0] + data[id2] + data[id3] + data[id5] + data[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id6] = ((data[id1] - data[id0] - data[id2] + data[id3] + data[id4] + data[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id7] = ((data[id2] - data[id0] + data[id3] + data[id4] + data[id5] - data[id6] + data[id7]) % PlainMod + PlainMod) % PlainMod;
    }
    std::copy(temp.begin(), temp.end(), data.begin());
}

// 行移位：A = data * constMR
void YusP::MR32(vector<long> &data)
{
    std::vector<long> temp(32);
    std::array<int, 8> index = {0, 1, 8, 9, 16, 17, 24, 25};
    for (int i = 0; i < 4; i++)
    {
        int s = 2 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s;
        temp[id0] = ((data[id1] + data[id2] + data[id4] + data[id5] + data[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id1] = ((data[id0] + data[id3] + data[id4] - data[id7]) % PlainMod + PlainMod) % PlainMod;
        temp[id2] = ((data[id1] + data[id2] + data[id3] + data[id4] - data[id6] - data[id7]) % PlainMod + PlainMod) % PlainMod;
        temp[id3] = ((data[id0] + data[id1] + data[id3] + data[id5] - data[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id4] = ((data[id0] + data[id1] + data[id2] + data[id4] - data[id6] + data[id7]) % PlainMod + PlainMod) % PlainMod;
        temp[id5] = ((data[id0] + data[id2] + data[id3] + data[id5] + data[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id6] = ((data[id1] - data[id0] - data[id2] + data[id3] + data[id4] + data[id6]) % PlainMod + PlainMod) % PlainMod;
        temp[id7] = ((data[id2] - data[id0] + data[id3] + data[id4] + data[id5] - data[id6] + data[id7]) % PlainMod + PlainMod) % PlainMod;
    }
    std::copy(temp.begin(), temp.end(), data.begin());
}
void YusP::MC_MR(vector<long> &data)
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
    std::vector<long> temp(32);
    for (int i = 0; i < 32; i++)
    {
        temp[i] = 0;
        for (int j = 0; j < 32; j++)
        {
            temp[i] += M[i][j] * data[j];
        }
        temp[i] = (temp[i] % PlainMod + PlainMod) % PlainMod;
    }
    std::copy(temp.begin(), temp.end(), data.begin());
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
data[id0] + data[id1] + t[2] + data[id3] + t[4] + data[id5] + t[6] + t[7];
[2,1,0,2,1,2,1,2]
t[0] + data[id1]  + t[3] + data[id4]+t[5] + data[id6] + t[7];
[2,2,1,1,2,1,2,1]
t[0] + t[1] + data[id2] + data[id3] + t[4]+data[id5] + t[6] + data[id7];
[1,2,2,1,0,2,1,2]
data[id0] + t[1] + t[2] + data[id3] + t[5] + data[id6] + t[7];
[2,1,2,2,1,1,2,1]
t[0] + data[id1] + t[2] + t[3] + data[id4]+data[id5] + t[6] + data[id7];
[1,2,1,2,2,1,0,2]
data[id0] + t[1] + data[id2] + t[3] + t[4]+data[id5] + t[7];
[2,1,2,1,2,2,1,1]
t[0] + data[id1] + t[2] + data[id3] + t[4]+t[5] + data[id6] + data[id7];
[0,2,1,2,1,2,2,1]
t[1] + data[id2] + t[3] + data[id4]+t[5] + t[6] + data[id7];

        t[0] = (data[id0]+data[id0])%PlainMod;
        t[1] = (data[id1]+data[id1])%PlainMod;
        t[2] = (data[id2]+data[id2])%PlainMod;
        t[3] = (data[id3]+data[id3])%PlainMod;
        t[4] = (data[id4]+data[id4])%PlainMod;
        t[5] = (data[id5]+data[id5])%PlainMod;
        t[6] = (data[id6]+data[id6])%PlainMod;
        t[7] = (data[id7]+data[id7])%PlainMod;

data[id0] + data[id1] + t[2] + data[id3] + t[4] + data[id5] + t[6] + t[7];
t[0] + data[id1]  + t[3] + data[id4]+t[5] + data[id6] + t[7];
t[0] + t[1] + data[id2] + data[id3] + t[4]+data[id5] + t[6] + data[id7];
data[id0] + t[1] + t[2] + data[id3] + t[5] + data[id6] + t[7];
t[0] + data[id1] + t[2] + t[3] + data[id4]+data[id5] + t[6] + data[id7];
data[id0] + t[1] + data[id2] + t[3] + t[4]+data[id5] + t[7];
t[0] + data[id1] + t[2] + data[id3] + t[4]+t[5] + data[id6] + data[id7];
t[1] + data[id2] + t[3] + data[id4]+t[5] + t[6] + data[id7];
*/
void YusP::MC32_2(vector<long> &data)
{
    std::vector<long> temp = data;
    std::array<int, 8> index = {0, 1, 2, 3, 4, 5, 6, 7};
    std::vector<long> t = data;
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

        data[id0] = (data[id0] + temp[id1] + t[id2] + temp[id3] + t[id4] + temp[id5] + t[id6] + t[id7]) % PlainMod;
        data[id1] = (t[id0] + data[id1] + t[id3] + temp[id4] + t[id5] + temp[id6] + t[id7]) % PlainMod;
        data[id2] = (t[id0] + t[id1] + data[id2] + temp[id3] + t[id4] + temp[id5] + t[id6] + temp[id7]) % PlainMod;
        data[id3] = (temp[id0] + t[id1] + t[id2] + data[id3] + t[id5] + temp[id6] + t[id7]) % PlainMod;
        data[id4] = (t[id0] + temp[id1] + t[id2] + t[id3] + data[id4] + temp[id5] + t[id6] + temp[id7]) % PlainMod;
        data[id5] = (temp[id0] + t[id1] + temp[id2] + t[id3] + t[id4] + data[id5] + t[id7]) % PlainMod;
        data[id6] = (t[id0] + temp[id1] + t[id2] + temp[id3] + t[id4] + t[id5] + data[id6] + temp[id7]) % PlainMod;
        data[id7] = (t[id1] + temp[id2] + t[id3] + temp[id4] + t[id5] + t[id6] + data[id7]) % PlainMod;
    }
}

// 行移位：A = data * constMR
void YusP::MR32_2(vector<long> &data)
{
    std::vector<long> temp = data;
    std::array<int, 8> index = {0, 1, 8, 9, 16, 17, 24, 25};
    std::vector<long> t = data;
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

        data[id0] = (data[id0] + temp[id1] + t[id2] + temp[id3] + t[id4] + temp[id5] + t[id6] + t[id7]) % PlainMod;
        data[id1] = (t[id0] + data[id1] + t[id3] + temp[id4] + t[id5] + temp[id6] + t[id7]) % PlainMod;
        data[id2] = (t[id0] + t[id1] + data[id2] + temp[id3] + t[id4] + temp[id5] + t[id6] + temp[id7]) % PlainMod;
        data[id3] = (temp[id0] + t[id1] + t[id2] + data[id3] + t[id5] + temp[id6] + t[id7]) % PlainMod;
        data[id4] = (t[id0] + temp[id1] + t[id2] + t[id3] + data[id4] + temp[id5] + t[id6] + temp[id7]) % PlainMod;
        data[id5] = (temp[id0] + t[id1] + temp[id2] + t[id3] + t[id4] + data[id5] + t[id7]) % PlainMod;
        data[id6] = (t[id0] + temp[id1] + t[id2] + temp[id3] + t[id4] + t[id5] + data[id6] + temp[id7]) % PlainMod;
        data[id7] = (t[id1] + temp[id2] + t[id3] + temp[id4] + t[id5] + t[id6] + data[id7]) % PlainMod;
    }
}

void YusP::MC64(vector<long> &data)
{
    std::vector<long> temp(64);
    std::array<int, 8> index = {0, 1, 2, 3, 4, 5, 6, 7};
    for (int i = 0; i < 8; i++)
    {
        int s = 8 * i;
        temp[index[0] + s] = ((data[index[1] + s] + data[index[2] + s] + data[index[4] + s] + data[index[5] + s] + data[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[1] + s] = ((data[index[0] + s] + data[index[3] + s] + data[index[4] + s] - data[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[2] + s] = ((data[index[1] + s] + data[index[2] + s] + data[index[3] + s] + data[index[4] + s] - data[index[6] + s] - data[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[3] + s] = ((data[index[0] + s] + data[index[1] + s] + data[index[3] + s] + data[index[5] + s] - data[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[4] + s] = ((data[index[0] + s] + data[index[1] + s] + data[index[2] + s] + data[index[4] + s] - data[index[6] + s] + data[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[5] + s] = ((data[index[0] + s] + data[index[2] + s] + data[index[3] + s] + data[index[5] + s] + data[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[6] + s] = ((data[index[1] + s] - data[index[0] + s] - data[index[2] + s] + data[index[3] + s] + data[index[4] + s] + data[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[7] + s] = ((data[index[2] + s] - data[index[0] + s] + data[index[3] + s] + data[index[4] + s] + data[index[5] + s] - data[index[6] + s] + data[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
    }
    std::copy(temp.begin(), temp.end(), data.begin());
}

// 行移位：A = data * constMR
void YusP::MR64(vector<long> &data)
{
    std::vector<long> temp(64);
    std::array<int, 8> index = {0, 8, 16, 24, 32, 40, 48, 56};
    for (int i = 0; i < 8; i++)
    {
        int s = i;
        temp[index[0] + s] = ((data[index[1] + s] + data[index[2] + s] + data[index[4] + s] + data[index[5] + s] + data[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[1] + s] = ((data[index[0] + s] + data[index[3] + s] + data[index[4] + s] - data[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[2] + s] = ((data[index[1] + s] + data[index[2] + s] + data[index[3] + s] + data[index[4] + s] - data[index[6] + s] - data[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[3] + s] = ((data[index[0] + s] + data[index[1] + s] + data[index[3] + s] + data[index[5] + s] - data[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[4] + s] = ((data[index[0] + s] + data[index[1] + s] + data[index[2] + s] + data[index[4] + s] - data[index[6] + s] + data[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[5] + s] = ((data[index[0] + s] + data[index[2] + s] + data[index[3] + s] + data[index[5] + s] + data[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[6] + s] = ((data[index[1] + s] - data[index[0] + s] - data[index[2] + s] + data[index[3] + s] + data[index[4] + s] + data[index[6] + s]) % PlainMod + PlainMod) % PlainMod;
        temp[index[7] + s] = ((data[index[2] + s] - data[index[0] + s] + data[index[3] + s] + data[index[4] + s] + data[index[5] + s] - data[index[6] + s] + data[index[7] + s]) % PlainMod + PlainMod) % PlainMod;
    }
    std::copy(temp.begin(), temp.end(), data.begin());
}
void YusP::Sbox(vector<long> &data)
{
    long temp;
    for (int i = 0; i < data.size(); i += 2)
    {
        // 两次Feistel轮函数(x,y)——> (y+x^2,x)
        for (int j = 0; j < 2; j++)
        {
            temp = data[i];
            data[i] = ((data[i + 1] + data[i] * data[i]) % PlainMod + PlainMod) % PlainMod;
            data[i + 1] = temp;
        }
    }
}
void YusP::MC48_3(vector<long> &data)
{
    std::vector<long> temp = data;
    std::array<int, 12> index = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    for (int i = 0; i < 4; i++)
    {
        int s = 12 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s, id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s;
        data[id0] = (((temp[id1] + temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id7] + temp[id9] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        data[id1] = (((temp[id0] + temp[id1] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id9]) % PlainMod) + PlainMod) % PlainMod;
        data[id2] = (((temp[id0] + temp[id1] + temp[id2] + temp[id4] + temp[id5] + temp[id6] + temp[id10]) % PlainMod) + PlainMod) % PlainMod;
        data[id3] = (((temp[id0] + temp[id2] + temp[id4] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id9] + temp[id10]) % PlainMod) + PlainMod) % PlainMod;
        data[id4] = (((temp[id2] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        data[id5] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id9]) % PlainMod) + PlainMod) % PlainMod;
        data[id6] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        data[id7] = (((temp[id0] + temp[id1] + temp[id2] + temp[id3] + temp[id6] + temp[id7] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        data[id8] = (((temp[id0] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        data[id9] = (((temp[id0] + temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id6] + temp[id8] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        data[id10] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id9] + temp[id10]) % PlainMod) + PlainMod) % PlainMod;
        data[id11] = (((temp[id1] + temp[id3] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::MR48_3(vector<long> &data)
{
    std::vector<long> temp = data;
    std::array<int, 12> index = {0, 1, 2, 12, 13, 14, 24, 25, 26, 36, 37, 38};

    for (int i = 0; i < 4; i++)
    {
        int s = 3 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s, id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s;
        data[id0] = (((temp[id1] + temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id7] + temp[id9] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        data[id1] = (((temp[id0] + temp[id1] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id9]) % PlainMod) + PlainMod) % PlainMod;
        data[id2] = (((temp[id0] + temp[id1] + temp[id2] + temp[id4] + temp[id5] + temp[id6] + temp[id10]) % PlainMod) + PlainMod) % PlainMod;
        data[id3] = (((temp[id0] + temp[id2] + temp[id4] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id9] + temp[id10]) % PlainMod) + PlainMod) % PlainMod;
        data[id4] = (((temp[id2] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        data[id5] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id9]) % PlainMod) + PlainMod) % PlainMod;
        data[id6] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        data[id7] = (((temp[id0] + temp[id1] + temp[id2] + temp[id3] + temp[id6] + temp[id7] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        data[id8] = (((temp[id0] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        data[id9] = (((temp[id0] + temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id6] + temp[id8] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
        data[id10] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id9] + temp[id10]) % PlainMod) + PlainMod) % PlainMod;
        data[id11] = (((temp[id1] + temp[id3] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id11]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::Sbox_3(vector<long> &data)
{
    vector<long> temp = data;
    // 一次迭代，(x1,x2,x3)——> (x1,x1*x3+x2,-x1*x2+x1*x3+x3)
    for (int i = 0; i < data.size(); i += 3)
    {
        data[i + 1] = (temp[i] * temp[i + 2] + temp[i + 1]) % PlainMod;
        data[i + 2] = (((-temp[i] * temp[i + 1] + temp[i] * temp[i + 2] + temp[i + 2]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::MC64_4(vector<long> &data)
{
    std::vector<long> temp = data;
    std::array<int, 16> index = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    for (int i = 0; i < 4; i++)
    {
        int s = 16 * i;

        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s,
            id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s,
            id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s,
            id12 = index[12] + s, id13 = index[13] + s, id14 = index[14] + s, id15 = index[15] + s;

        data[id0] = (((temp[id0] + temp[id2] + temp[id4] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id13]) % PlainMod) + PlainMod) % PlainMod;
        data[id1] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id7] + temp[id8] + temp[id10] + temp[id11] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id2] = (((temp[id0] + temp[id1] + temp[id2] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id3] = (((temp[id1] + temp[id2] + temp[id6] + temp[id7] + temp[id8] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id4] = (((temp[id0] + temp[id4] + temp[id6] + temp[id8] + temp[id9] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id5] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id11] + temp[id12] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id6] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id6] + temp[id9] + temp[id11] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id7] = (((temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        data[id8] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        data[id9] = (((temp[id0] + temp[id2] + temp[id3] + temp[id7] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id10] = (((temp[id0] + temp[id2] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id13] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id11] = (((temp[id0] + temp[id6] + temp[id9] + temp[id10] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id12] = (((temp[id0] + temp[id1] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id13] = (((temp[id3] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id11] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id14] = (((temp[id1] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id10] + temp[id11] + temp[id12] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id15] = (((temp[id2] + temp[id4] + temp[id10] + temp[id11] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::MR64_4(vector<long> &data)
{
    std::vector<long> temp = data;
    std::array<int, 16> index = {0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51};

    for (int i = 0; i < 4; i++)
    {
        int s = 4 * i;

        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s,
            id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s,
            id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s,
            id12 = index[12] + s, id13 = index[13] + s, id14 = index[14] + s, id15 = index[15] + s;

        data[id0] = (((temp[id0] + temp[id2] + temp[id4] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id13]) % PlainMod) + PlainMod) % PlainMod;
        data[id1] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id7] + temp[id8] + temp[id10] + temp[id11] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id2] = (((temp[id0] + temp[id1] + temp[id2] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id3] = (((temp[id1] + temp[id2] + temp[id6] + temp[id7] + temp[id8] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id4] = (((temp[id0] + temp[id4] + temp[id6] + temp[id8] + temp[id9] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id5] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id11] + temp[id12] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id6] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id6] + temp[id9] + temp[id11] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id7] = (((temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        data[id8] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        data[id9] = (((temp[id0] + temp[id2] + temp[id3] + temp[id7] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id10] = (((temp[id0] + temp[id2] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id13] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id11] = (((temp[id0] + temp[id6] + temp[id9] + temp[id10] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id12] = (((temp[id0] + temp[id1] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id13] = (((temp[id3] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id11] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id14] = (((temp[id1] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id10] + temp[id11] + temp[id12] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id15] = (((temp[id2] + temp[id4] + temp[id10] + temp[id11] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::Sbox_4(vector<long> &data)
{
    vector<long> temp = data;
    // 一次迭代，(x0,x1,x2,x3)——> (x0, x0^2 + x1, -x0*x3 -x1^2 + x2, x0*x2 -x0*x3 + x1^2 + x3)
    for (int i = 0; i < data.size(); i += 4)
    {
        data[i + 1] = (temp[i] * temp[i] + temp[i + 1]) % PlainMod;
        data[i + 2] = (((-temp[i] * temp[i + 3] - temp[i + 1] * temp[i + 1] + temp[i + 2]) % PlainMod) + PlainMod) % PlainMod;
        data[i + 3] = ((temp[i] * temp[i + 2] - temp[i] * temp[i + 3] + temp[i + 1] * temp[i + 1] + temp[i + 3]) % PlainMod + PlainMod) % PlainMod;
    }
}
void YusP::M36_5(vector<long> &data)
{
    std::vector<long> temp = data;
    // for (int i = 0; i < 1; i++)// i<2
    // {
        // data[0] = (((temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7] + temp[9] + temp[11] + temp[13] + temp[16] + temp[17] + temp[19] + temp[20] + temp[22] + temp[25] + temp[26] + temp[27] + temp[29] + temp[30] + temp[33] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[1] = (((temp[1] + temp[3] + temp[4] + temp[6] + temp[7] + temp[8] + temp[9] + temp[10] + temp[12] + temp[17] + temp[18] + temp[19] + temp[21] + temp[22] + temp[24] + temp[27] + temp[31] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[2] = (((temp[0] + temp[3] + temp[7] + temp[9] + temp[10] + temp[11] + temp[13] + temp[14] + temp[15] + temp[17] + temp[18] + temp[21] + temp[22] + temp[23] + temp[26] + temp[28] + temp[32] + temp[33] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
        // data[3] = (((temp[0] + temp[2] + temp[4] + temp[5] + temp[6] + temp[7] + temp[8] + temp[9] + temp[10] + temp[12] + temp[14] + temp[16] + temp[19] + temp[20] + temp[22] + temp[23] + temp[25] + temp[28] + temp[29] + temp[30] + temp[32] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
        // data[4] = (((temp[0] + temp[1] + temp[2] + temp[4] + temp[6] + temp[7] + temp[9] + temp[10] + temp[11] + temp[12] + temp[13] + temp[15] + temp[20] + temp[21] + temp[22] + temp[24] + temp[25] + temp[27] + temp[30] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
        // data[5] = (((temp[0] + temp[1] + temp[3] + temp[6] + temp[10] + temp[12] + temp[13] + temp[14] + temp[16] + temp[17] + temp[18] + temp[20] + temp[21] + temp[24] + temp[25] + temp[26] + temp[29] + temp[31] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[6] = (((temp[0] + temp[3] + temp[5] + temp[7] + temp[8] + temp[9] + temp[10] + temp[11] + temp[12] + temp[13] + temp[15] + temp[17] + temp[19] + temp[22] + temp[23] + temp[25] + temp[26] + temp[28] + temp[31] + temp[32] + temp[33] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[7] = (((temp[1] + temp[3] + temp[4] + temp[5] + temp[7] + temp[9] + temp[10] + temp[12] + temp[13] + temp[14] + temp[15] + temp[16] + temp[18] + temp[23] + temp[24] + temp[25] + temp[27] + temp[28] + temp[30] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
        // data[8] = (((temp[2] + temp[3] + temp[4] + temp[6] + temp[9] + temp[13] + temp[15] + temp[16] + temp[17] + temp[19] + temp[20] + temp[21] + temp[23] + temp[24] + temp[27] + temp[28] + temp[29] + temp[32] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
        // data[9] = (((temp[0] + temp[2] + temp[3] + temp[6] + temp[8] + temp[10] + temp[11] + temp[12] + temp[13] + temp[14] + temp[15] + temp[16] + temp[18] + temp[20] + temp[22] + temp[25] + temp[26] + temp[28] + temp[29] + temp[31] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[10] = (((temp[0] + temp[4] + temp[6] + temp[7] + temp[8] + temp[10] + temp[12] + temp[13] + temp[15] + temp[16] + temp[17] + temp[18] + temp[19] + temp[21] + temp[26] + temp[27] + temp[28] + temp[30] + temp[31] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
        // data[11] = (((temp[1] + temp[5] + temp[6] + temp[7] + temp[9] + temp[12] + temp[16] + temp[18] + temp[19] + temp[20] + temp[22] + temp[23] + temp[24] + temp[26] + temp[27] + temp[30] + temp[31] + temp[32] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[12] = (((temp[1] + temp[2] + temp[3] + temp[5] + temp[6] + temp[9] + temp[11] + temp[13] + temp[14] + temp[15] + temp[16] + temp[17] + temp[18] + temp[19] + temp[21] + temp[23] + temp[25] + temp[28] + temp[29] + temp[31] + temp[32] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
        // data[13] = (((temp[0] + temp[3] + temp[7] + temp[9] + temp[10] + temp[11] + temp[13] + temp[15] + temp[16] + temp[18] + temp[19] + temp[20] + temp[21] + temp[22] + temp[24] + temp[29] + temp[30] + temp[31] + temp[33] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
        // data[14] = (((temp[2] + temp[4] + temp[8] + temp[9] + temp[10] + temp[12] + temp[15] + temp[19] + temp[21] + temp[22] + temp[23] + temp[25] + temp[26] + temp[27] + temp[29] + temp[30] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[15] = (((temp[1] + temp[4] + temp[5] + temp[6] + temp[8] + temp[9] + temp[12] + temp[14] + temp[16] + temp[17] + temp[18] + temp[19] + temp[20] + temp[21] + temp[22] + temp[24] + temp[26] + temp[28] + temp[31] + temp[32] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[16] = (((temp[0] + temp[1] + temp[3] + temp[6] + temp[10] + temp[12] + temp[13] + temp[14] + temp[16] + temp[18] + temp[19] + temp[21] + temp[22] + temp[23] + temp[24] + temp[25] + temp[27] + temp[32] + temp[33] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
        // data[17] = (((temp[0] + temp[1] + temp[2] + temp[5] + temp[7] + temp[11] + temp[12] + temp[13] + temp[15] + temp[18] + temp[22] + temp[24] + temp[25] + temp[26] + temp[28] + temp[29] + temp[30] + temp[32] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
        // data[18] = (((temp[1] + temp[2] + temp[4] + temp[7] + temp[8] + temp[9] + temp[11] + temp[12] + temp[15] + temp[17] + temp[19] + temp[20] + temp[21] + temp[22] + temp[23] + temp[24] + temp[25] + temp[27] + temp[29] + temp[31] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[19] = (((temp[0] + temp[1] + temp[3] + temp[4] + temp[6] + temp[9] + temp[13] + temp[15] + temp[16] + temp[17] + temp[19] + temp[21] + temp[22] + temp[24] + temp[25] + temp[26] + temp[27] + temp[28] + temp[30] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[20] = (((temp[0] + temp[3] + temp[4] + temp[5] + temp[8] + temp[10] + temp[14] + temp[15] + temp[16] + temp[18] + temp[21] + temp[25] + temp[27] + temp[28] + temp[29] + temp[31] + temp[32] + temp[33] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[21] = (((temp[1] + temp[2] + temp[4] + temp[5] + temp[7] + temp[10] + temp[11] + temp[12] + temp[14] + temp[15] + temp[18] + temp[20] + temp[22] + temp[23] + temp[24] + temp[25] + temp[26] + temp[27] + temp[28] + temp[30] + temp[32] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
        // data[22] = (((temp[2] + temp[3] + temp[4] + temp[6] + temp[7] + temp[9] + temp[12] + temp[16] + temp[18] + temp[19] + temp[20] + temp[22] + temp[24] + temp[25] + temp[27] + temp[28] + temp[29] + temp[30] + temp[31] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
        // data[23] = (((temp[0] + temp[2] + temp[3] + temp[6] + temp[7] + temp[8] + temp[11] + temp[13] + temp[17] + temp[18] + temp[19] + temp[21] + temp[24] + temp[28] + temp[30] + temp[31] + temp[32] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[24] = (((temp[1] + temp[4] + temp[5] + temp[7] + temp[8] + temp[10] + temp[13] + temp[14] + temp[15] + temp[17] + temp[18] + temp[21] + temp[23] + temp[25] + temp[26] + temp[27] + temp[28] + temp[29] + temp[30] + temp[31] + temp[33] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[25] = (((temp[0] + temp[5] + temp[6] + temp[7] + temp[9] + temp[10] + temp[12] + temp[15] + temp[19] + temp[21] + temp[22] + temp[23] + temp[25] + temp[27] + temp[28] + temp[30] + temp[31] + temp[32] + temp[33] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
        // data[26] = (((temp[1] + temp[2] + temp[3] + temp[5] + temp[6] + temp[9] + temp[10] + temp[11] + temp[14] + temp[16] + temp[20] + temp[21] + temp[22] + temp[24] + temp[27] + temp[31] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[27] = (((temp[0] + temp[2] + temp[4] + temp[7] + temp[8] + temp[10] + temp[11] + temp[13] + temp[16] + temp[17] + temp[18] + temp[20] + temp[21] + temp[24] + temp[26] + temp[28] + temp[29] + temp[30] + temp[31] + temp[32] + temp[33] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
        // data[28] = (((temp[0] + temp[1] + temp[3] + temp[8] + temp[9] + temp[10] + temp[12] + temp[13] + temp[15] + temp[18] + temp[22] + temp[24] + temp[25] + temp[26] + temp[28] + temp[30] + temp[31] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[29] = (((temp[0] + temp[1] + temp[2] + temp[4] + temp[5] + temp[6] + temp[8] + temp[9] + temp[12] + temp[13] + temp[14] + temp[17] + temp[19] + temp[23] + temp[24] + temp[25] + temp[27] + temp[30] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
        // data[30] = (((temp[0] + temp[1] + temp[3] + temp[5] + temp[7] + temp[10] + temp[11] + temp[13] + temp[14] + temp[16] + temp[19] + temp[20] + temp[21] + temp[23] + temp[24] + temp[27] + temp[29] + temp[31] + temp[32] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[31] = (((temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[6] + temp[11] + temp[12] + temp[13] + temp[15] + temp[16] + temp[18] + temp[21] + temp[25] + temp[27] + temp[28] + temp[29] + temp[31] + temp[33] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
        // data[32] = (((temp[1] + temp[3] + temp[4] + temp[5] + temp[7] + temp[8] + temp[9] + temp[11] + temp[12] + temp[15] + temp[16] + temp[17] + temp[20] + temp[22] + temp[26] + temp[27] + temp[28] + temp[30] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
        // data[33] = (((temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[6] + temp[8] + temp[10] + temp[13] + temp[14] + temp[16] + temp[17] + temp[19] + temp[22] + temp[23] + temp[24] + temp[26] + temp[27] + temp[30] + temp[32] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
        // data[34] = (((temp[0] + temp[1] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7] + temp[9] + temp[14] + temp[15] + temp[16] + temp[18] + temp[19] + temp[21] + temp[24] + temp[28] + temp[30] + temp[31] + temp[32] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
        // data[35] = (((temp[0] + temp[4] + temp[6] + temp[7] + temp[8] + temp[10] + temp[11] + temp[12] + temp[14] + temp[15] + temp[18] + temp[19] + temp[20] + temp[23] + temp[25] + temp[29] + temp[30] + temp[31] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
    // 0,1,2,3
    long temp0_1 = temp[0];
    temp0_1 += temp[1];
    long temp0_2 = temp[0];
    temp0_2 += temp[2];
    long temp0_3 = temp[0];
    temp0_3 += temp[3];
    long temp1_2 = temp[1];
    temp1_2 += temp[2];
    long temp1_3 = temp[1];
    temp1_3 += temp[3];
    long temp2_3 = temp[2];
    temp2_3 += temp[3];
    long temp0_1_2 = temp0_1;
    temp0_1_2 += temp[2];
    long temp0_1_3 = temp0_1;
    temp0_1_3 += temp[3];
    long temp0_2_3 = temp0_2;
    temp0_2_3 += temp[3];
    long temp1_2_3 = temp1_2;
    temp1_2_3 += temp[3];
    long temp0_1_2_3 = temp0_1_2;
    temp0_1_2_3 += temp[3];
    // 4,5,6,7
    long temp4_5 = temp[4];
    temp4_5 += temp[5];
    long temp4_6 = temp[4];
    temp4_6 += temp[6];
    long temp4_7 = temp[4];
    temp4_7 += temp[7];
    long temp5_6 = temp[5];
    temp5_6 += temp[6];
    long temp5_7 = temp[5];
    temp5_7 += temp[7];
    long temp6_7 = temp[6];
    temp6_7 += temp[7];
    long temp4_5_6 = temp4_5;
    temp4_5_6 += temp[6];
    long temp4_5_7 = temp4_5;
    temp4_5_7 += temp[7];
    long temp4_6_7 = temp4_6;
    temp4_6_7 += temp[7];
    long temp5_6_7 = temp5_6;
    temp5_6_7 += temp[7];
    long temp4_5_6_7 = temp4_5_6;
    temp4_5_6_7 += temp[7];
    // 8,9,10,11
    long temp8_9 = temp[8];
    temp8_9 += temp[9];
    long temp8_10 = temp[8];
    temp8_10 += temp[10];
    long temp8_11 = temp[8];
    temp8_11 += temp[11];
    long temp9_10 = temp[9];
    temp9_10 += temp[10];
    long temp9_11 = temp[9];
    temp9_11 += temp[11];
    long temp10_11 = temp[10];
    temp10_11 += temp[11];
    long temp8_9_10 = temp8_9;
    temp8_9_10 += temp[10];
    long temp8_9_11 = temp8_9;
    temp8_9_11 += temp[11];
    long temp8_10_11 = temp8_10;
    temp8_10_11 += temp[11];
    long temp9_10_11 = temp9_10;
    temp9_10_11 += temp[11];
    long temp8_9_10_11 = temp8_9_10;
    temp8_9_10_11 += temp[11];
    // 12,13,14,15
    long temp12_13 = temp[12];
    temp12_13 += temp[13];
    long temp12_14 = temp[12];
    temp12_14 += temp[14];
    long temp12_15 = temp[12];
    temp12_15 += temp[15];
    long temp13_14 = temp[13];
    temp13_14 += temp[14];
    long temp13_15 = temp[13];
    temp13_15 += temp[15];
    long temp14_15 = temp[14];
    temp14_15 += temp[15];
    long temp12_13_14 = temp12_13;
    temp12_13_14 += temp[14];
    long temp12_13_15 = temp12_13;
    temp12_13_15 += temp[15];
    long temp12_14_15 = temp12_14;
    temp12_14_15 += temp[15];
    long temp13_14_15 = temp13_14;
    temp13_14_15 += temp[15];
    long temp12_13_14_15 = temp12_13_14;
    temp12_13_14_15 += temp[15];
    // 16,17,18,19
    long temp16_17 = temp[16];
    temp16_17 += temp[17];
    long temp16_18 = temp[16];
    temp16_18 += temp[18];
    long temp16_19 = temp[16];
    temp16_19 += temp[19];
    long temp17_18 = temp[17];
    temp17_18 += temp[18];
    long temp17_19 = temp[17];
    temp17_19 += temp[19];
    long temp18_19 = temp[18];
    temp18_19 += temp[19];
    long temp16_17_18 = temp16_17;
    temp16_17_18 += temp[18];
    long temp16_17_19 = temp16_17;
    temp16_17_19 += temp[19];
    long temp16_18_19 = temp16_18;
    temp16_18_19 += temp[19];
    long temp17_18_19 = temp17_18;
    temp17_18_19 += temp[19];
    long temp16_17_18_19 = temp16_17_18;
    temp16_17_18_19 += temp[19];
    // 20,21,22,23
    long temp20_21 = temp[20];
    temp20_21 += temp[21];
    long temp20_22 = temp[20];
    temp20_22 += temp[22];
    long temp20_23 = temp[20];
    temp20_23 += temp[23];
    long temp21_22 = temp[21];
    temp21_22 += temp[22];
    long temp21_23 = temp[21];
    temp21_23 += temp[23];
    long temp22_23 = temp[22];
    temp22_23 += temp[23];
    long temp20_21_22 = temp20_21;
    temp20_21_22 += temp[22];
    long temp20_21_23 = temp20_21;
    temp20_21_23 += temp[23];
    long temp20_22_23 = temp20_22;
    temp20_22_23 += temp[23];
    long temp21_22_23 = temp21_22;
    temp21_22_23 += temp[23];
    long temp20_21_22_23 = temp20_21_22;
    temp20_21_22_23 += temp[23];
    // 24,25,26,27
    long temp24_25 = temp[24];
    temp24_25 += temp[25];
    long temp24_26 = temp[24];
    temp24_26 += temp[26];
    long temp24_27 = temp[24];
    temp24_27 += temp[27];
    long temp25_26 = temp[25];
    temp25_26 += temp[26];
    long temp25_27 = temp[25];
    temp25_27 += temp[27];
    long temp26_27 = temp[26];
    temp26_27 += temp[27];
    long temp24_25_26 = temp24_25;
    temp24_25_26 += temp[26];
    long temp24_25_27 = temp24_25;
    temp24_25_27 += temp[27];
    long temp24_26_27 = temp24_26;
    temp24_26_27 += temp[27];
    long temp25_26_27 = temp25_26;
    temp25_26_27 += temp[27];
    long temp24_25_26_27 = temp24_25_26;
    temp24_25_26_27 += temp[27];
    // 28,29,30,31
    long temp28_29 = temp[28];
    temp28_29 += temp[29];
    long temp28_30 = temp[28];
    temp28_30 += temp[30];
    long temp28_31 = temp[28];
    temp28_31 += temp[31];
    long temp29_30 = temp[29];
    temp29_30 += temp[30];
    long temp29_31 = temp[29];
    temp29_31 += temp[31];
    long temp30_31 = temp[30];
    temp30_31 += temp[31];
    long temp28_29_30 = temp28_29;
    temp28_29_30 += temp[30];
    long temp28_29_31 = temp28_29;
    temp28_29_31 += temp[31];
    long temp28_30_31 = temp28_30;
    temp28_30_31 += temp[31];
    long temp29_30_31 = temp29_30;
    temp29_30_31 += temp[31];
    long temp28_29_30_31 = temp28_29_30;
    temp28_29_30_31 += temp[31];
    // 32,33,34,35
    long temp32_33 = temp[32];
    temp32_33 += temp[33];
    long temp32_34 = temp[32];
    temp32_34 += temp[34];
    long temp32_35 = temp[32];
    temp32_35 += temp[35];
    long temp33_34 = temp[33];
    temp33_34 += temp[34];
    long temp33_35 = temp[33];
    temp33_35 += temp[35];
    long temp34_35 = temp[34];
    temp34_35 += temp[35];
    long temp32_33_34 = temp32_33;
    temp32_33_34 += temp[34];
    long temp32_33_35 = temp32_33;
    temp32_33_35 += temp[35];
    long temp32_34_35 = temp32_34;
    temp32_34_35 += temp[35];
    long temp33_34_35 = temp33_34;
    temp33_34_35 += temp[35];
    long temp32_33_34_35 = temp32_33_34;
    temp32_33_34_35 += temp[35];

    data[0] = temp1_2_3;
    data[0] += temp4_5_6_7;
    data[0] += temp9_11;
    data[0] += temp[13];
    data[0] += temp16_17_19;
    data[0] += temp20_22;
    data[0] += temp25_26_27;
    data[0] += temp29_30;
    data[0] += temp33_35;

    data[1] += temp[3];
    data[1] += temp4_6_7;
    data[1] += temp8_9_10;
    data[1] += temp[12];
    data[1] += temp17_18_19;
    data[1] += temp21_22;
    data[1] += temp24_27;
    data[1] += temp[31];
    data[1] += temp33_34_35;

    data[2] = temp0_3;
    data[2] += temp[7];
    data[2] += temp9_10_11;
    data[2] += temp13_14_15;
    data[2] += temp17_18;
    data[2] += temp21_22_23;
    data[2] += temp[26];
    data[2] += temp[28];
    data[2] += temp32_33_34;

    data[3] = temp0_2;
    data[3] += temp4_5_6_7;
    data[3] += temp8_9_10;
    data[3] += temp12_14;
    data[3] += temp16_19;
    data[3] += temp20_22_23;
    data[3] += temp[25];
    data[3] += temp28_29_30;
    data[3] += temp32_33;

    data[4] += temp0_1_2;
    data[4] += temp6_7;
    data[4] += temp9_10_11;
    data[4] += temp12_13_15;
    data[4] += temp20_21_22;
    data[4] += temp24_25_27;
    data[4] += temp[30];
    data[4] += temp[34];

    data[5] = temp0_1_3;
    data[5] += temp[6];
    data[5] += temp[10];
    data[5] += temp12_13_14;
    data[5] += temp16_17_18;
    data[5] += temp20_21;
    data[5] += temp24_25_26;
    data[5] += temp29_31;
    data[5] += temp[35];

    data[6] = temp0_3;
    data[6] += temp5_7;
    data[6] += temp8_9_10_11;
    data[6] += temp12_13_15;
    data[6] += temp17_19;
    data[6] += temp22_23;
    data[6] += temp25_26;
    data[6] += temp28_31;
    data[6] += temp32_33_35;

    data[7] += temp1_3;
    data[7] += temp4_5;
    data[7] += temp9_10;
    data[7] += temp12_13_14_15;
    data[7] += temp16_18;
    data[7] += temp[23];
    data[7] += temp24_25_27;
    data[7] += temp28_30;
    data[7] += temp[33];

    data[8] = temp2_3;
    data[8] += temp4_6;
    data[8] += temp[9];
    data[8] += temp13_15;
    data[8] += temp16_17_19;
    data[8] += temp20_21_23;
    data[8] += temp24_27;
    data[8] += temp28_29;
    data[8] += temp32_34;

    data[9] = temp0_2_3;
    data[9] += temp[6];
    data[9] += temp8_10_11;
    data[9] += temp12_13_14_15;
    data[9] += temp16_18;
    data[9] += temp20_22;
    data[9] += temp25_26;
    data[9] += temp28_29_31;
    data[9] += temp34_35;

    data[10] += temp[0];
    data[10] += temp4_6_7;
    data[10] += temp[8];
    data[10] += temp12_13_15;
    data[10] += temp16_17_18_19;
    data[10] += temp[21];
    data[10] += temp26_27;
    data[10] += temp28_30_31;
    data[10] += temp[33];

    data[11] = temp[1];
    data[11] += temp5_6_7;
    data[11] += temp[9];
    data[11] += temp[12];
    data[11] += temp16_18_19;
    data[11] += temp20_22_23;
    data[11] += temp24_26_27;
    data[11] += temp30_31;
    data[11] += temp32_35;

    data[12] = temp1_2_3;
    data[12] += temp5_6;
    data[12] += temp9_11;
    data[12] += temp13_14_15;
    data[12] += temp16_17_18_19;
    data[12] += temp21_23;
    data[12] += temp[25];
    data[12] += temp28_29_31;
    data[12] += temp32_34;

    data[13] += temp0_3;
    data[13] += temp[7];
    data[13] += temp9_10_11;
    data[13] += temp[15];
    data[13] += temp16_18_19;
    data[13] += temp20_21_22;
    data[13] += temp[24];
    data[13] += temp29_30_31;
    data[13] += temp33_34;

    data[14] = temp[2];
    data[14] += temp[4];
    data[14] += temp8_9_10;
    data[14] += temp12_15;
    data[14] += temp[19];
    data[14] += temp21_22_23;
    data[14] += temp25_26_27;
    data[14] += temp29_30;
    data[14] += temp33_34_35;

    data[15] = temp[1];
    data[15] += temp4_5_6;
    data[15] += temp8_9;
    data[15] += temp12_14;
    data[15] += temp16_17_18_19;
    data[15] += temp20_21_22;
    data[15] += temp24_26;
    data[15] += temp28_31;
    data[15] += temp32_34_35;

    data[16] += temp0_1_3;
    data[16] += temp[6];
    data[16] += temp[10];
    data[16] += temp12_13_14;
    data[16] += temp18_19;
    data[16] += temp21_22_23;
    data[16] += temp24_25_27;
    data[16] += temp32_33_34;

    data[17] = temp0_1_2;
    data[17] += temp5_7;
    data[17] += temp[11];
    data[17] += temp12_13_15;
    data[17] += temp[18];
    data[17] += temp[22];
    data[17] += temp24_25_26;
    data[17] += temp28_29_30;
    data[17] += temp32_33;

    data[18] = temp1_2;
    data[18] += temp4_7;
    data[18] += temp8_9_11;
    data[18] += temp12_15;
    data[18] += temp17_19;
    data[18] += temp20_21_22_23;
    data[18] += temp24_25_27;
    data[18] += temp29_31;
    data[18] += temp34_35;

    data[19] += temp0_1_3;
    data[19] += temp4_6;
    data[19] += temp[9];
    data[19] += temp13_15;
    data[19] += temp16_17;
    data[19] += temp21_22;
    data[19] += temp24_25_26_27;
    data[19] += temp28_30;
    data[19] += temp[35];

    data[20] = temp0_3;
    data[20] += temp4_5;
    data[20] += temp8_10;
    data[20] += temp14_15;
    data[20] += temp16_18;
    data[20] += temp[21];
    data[20] += temp25_27;
    data[20] += temp28_29_31;
    data[20] += temp32_33_35;

    data[21] = temp1_2;
    data[21] += temp4_5_7;
    data[21] += temp10_11;
    data[21] += temp12_14_15;
    data[21] += temp[18];
    data[21] += temp20_22_23;
    data[21] += temp24_25_26_27;
    data[21] += temp28_30;
    data[21] += temp32_34;

    data[22] += temp2_3;
    data[22] += temp4_6_7;
    data[22] += temp[9];
    data[22] += temp[12];
    data[22] += temp16_18_19;
    data[22] += temp[20];
    data[22] += temp24_25_27;
    data[22] += temp28_29_30_31;
    data[22] += temp[33];

    data[23] = temp0_2_3;
    data[23] += temp6_7;
    data[23] += temp8_11;
    data[23] += temp[13];
    data[23] += temp17_18_19;
    data[23] += temp[21];
    data[23] += temp[24];
    data[23] += temp28_30_31;
    data[23] += temp32_34_35;

    data[24] = temp[1];
    data[24] += temp4_5_7;
    data[24] += temp8_10;
    data[24] += temp13_14_15;
    data[24] += temp17_18;
    data[24] += temp21_23;
    data[24] += temp25_26_27;
    data[24] += temp28_29_30_31;
    data[24] += temp33_35;

    data[25] += temp[0];
    data[25] += temp5_6_7;
    data[25] += temp9_10;
    data[25] += temp12_15;
    data[25] += temp[19];
    data[25] += temp21_22_23;
    data[25] += temp[27];
    data[25] += temp28_30_31;
    data[25] += temp32_33_34;

    data[26] = temp1_2_3;
    data[26] += temp5_6;
    data[26] += temp9_10_11;
    data[26] += temp[14];
    data[26] += temp[16];
    data[26] += temp20_21_22;
    data[26] += temp24_27;
    data[26] += temp[31];
    data[26] += temp33_34_35;

    data[27] = temp0_2;
    data[27] += temp4_7;
    data[27] += temp8_10_11;
    data[27] += temp[13];
    data[27] += temp16_17_18;
    data[27] += temp20_21;
    data[27] += temp24_26;
    data[27] += temp28_29_30_31;
    data[27] += temp32_33_34;

    data[28] += temp0_1_3;
    data[28] += temp8_9_10;
    data[28] += temp12_13_15;
    data[28] += temp[18];
    data[28] += temp[22];
    data[28] += temp24_25_26;
    data[28] += temp30_31;
    data[28] += temp33_34_35;

    data[29] = temp0_1_2;
    data[29] += temp4_5_6;
    data[29] += temp8_9;
    data[29] += temp12_13_14;
    data[29] += temp17_19;
    data[29] += temp[23];
    data[29] += temp24_25_27;
    data[29] += temp[30];
    data[29] += temp[34];

    data[30] = temp0_1_3;
    data[30] += temp5_7;
    data[30] += temp10_11;
    data[30] += temp13_14;
    data[30] += temp16_19;
    data[30] += temp20_21_23;
    data[30] += temp24_27;
    data[30] += temp29_31;
    data[30] += temp32_33_34_35;

    data[31] += temp0_1_2_3;
    data[31] += temp4_6;
    data[31] += temp[11];
    data[31] += temp12_13_15;
    data[31] += temp16_18;
    data[31] += temp[21];
    data[31] += temp25_27;
    data[31] += temp28_29;
    data[31] += temp33_34;

    data[32] = temp1_3;
    data[32] += temp4_5_7;
    data[32] += temp8_9_11;
    data[32] += temp12_15;
    data[32] += temp16_17;
    data[32] += temp20_22;
    data[32] += temp26_27;
    data[32] += temp28_30;
    data[32] += temp[33];

    data[33] = temp0_1_2_3;
    data[33] += temp4_6;
    data[33] += temp8_10;
    data[33] += temp13_14;
    data[33] += temp16_17_19;
    data[33] += temp22_23;
    data[33] += temp24_26_27;
    data[33] += temp[30];
    data[33] += temp32_34_35;

    data[34] += temp0_1_3;
    data[34] += temp4_5_6_7;
    data[34] += temp[9];
    data[34] += temp14_15;
    data[34] += temp16_18_19;
    data[34] += temp[21];
    data[34] += temp[24];
    data[34] += temp28_30_31;
    data[34] += temp[32];

    data[35] = temp[0];
    data[35] += temp4_6_7;
    data[35] += temp8_10_11;
    data[35] += temp12_14_15;
    data[35] += temp18_19;
    data[35] += temp20_23;
    data[35] += temp[25];
    data[35] += temp29_30_31;
    data[35] += temp[33];
        //temp = data;
    //}
}
void YusP::Sbox_5(vector<long> &data)
{
    vector<long> temp = data;
    // 一次迭代，(x0,x1,x2)——>(x0,x0*x2+x1,-x0*x1+x0*x2+x2)
    long t02;
    for (int i = 0; i < data.size(); i += 3)
    {
        t02 = (temp[i] * temp[i + 2]) % PlainMod;
        data[i + 1] = (t02 + temp[i + 1]) % PlainMod;
        data[i + 2] = (t02 -temp[i] * temp[i + 1]  + temp[i + 2]) % PlainMod;
    }
}
void YusP::Sbox_5_last(vector<long> &data)
{
    vector<long> temp = data;
    // 一次迭代，(x0,x1,x2)——>(x0^3,x0*x2+x1,-x0*x1+x0*x2+x2)
    for (int i = 0; i < data.size(); i += 3)
    {
        data[i] = (temp[i]*temp[i]*temp[i])%PlainMod;
        data[i + 1] = (temp[i] * temp[i + 2] + temp[i + 1]) % PlainMod;
        data[i + 2] = (((-temp[i] * temp[i + 1] + temp[i] * temp[i + 2] + temp[i + 2]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::MC64_6(vector<long> &data)
{
    std::vector<long> temp = data;
    std::array<int, 16> index = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    for (int i = 0; i < 4; i++)
    {
        int s = 16 * i;

        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s,
            id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s,
            id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s,
            id12 = index[12] + s, id13 = index[13] + s, id14 = index[14] + s, id15 = index[15] + s;

        data[id0] = (((temp[id0] + temp[id2] + temp[id4] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id13]) % PlainMod) + PlainMod) % PlainMod;
        data[id1] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id7] + temp[id8] + temp[id10] + temp[id11] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id2] = (((temp[id0] + temp[id1] + temp[id2] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id3] = (((temp[id1] + temp[id2] + temp[id6] + temp[id7] + temp[id8] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id4] = (((temp[id0] + temp[id4] + temp[id6] + temp[id8] + temp[id9] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id5] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id11] + temp[id12] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id6] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id6] + temp[id9] + temp[id11] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id7] = (((temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        data[id8] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        data[id9] = (((temp[id0] + temp[id2] + temp[id3] + temp[id7] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id10] = (((temp[id0] + temp[id2] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id13] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id11] = (((temp[id0] + temp[id6] + temp[id9] + temp[id10] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id12] = (((temp[id0] + temp[id1] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id13] = (((temp[id3] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id11] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id14] = (((temp[id1] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id10] + temp[id11] + temp[id12] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id15] = (((temp[id2] + temp[id4] + temp[id10] + temp[id11] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::MR64_6(vector<long> &data)
{
    std::vector<long> temp = data;
    std::array<int, 16> index = {0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51};

    for (int i = 0; i < 4; i++)
    {
        int s = 4 * i;

        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s,
            id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s,
            id8 = index[8] + s, id9 = index[9] + s, id10 = index[10] + s, id11 = index[11] + s,
            id12 = index[12] + s, id13 = index[13] + s, id14 = index[14] + s, id15 = index[15] + s;

        data[id0] = (((temp[id0] + temp[id2] + temp[id4] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id13]) % PlainMod) + PlainMod) % PlainMod;
        data[id1] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id7] + temp[id8] + temp[id10] + temp[id11] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id2] = (((temp[id0] + temp[id1] + temp[id2] + temp[id5] + temp[id6] + temp[id7] + temp[id8] + temp[id10] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id3] = (((temp[id1] + temp[id2] + temp[id6] + temp[id7] + temp[id8] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id4] = (((temp[id0] + temp[id4] + temp[id6] + temp[id8] + temp[id9] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id5] = (((temp[id0] + temp[id3] + temp[id4] + temp[id5] + temp[id7] + temp[id11] + temp[id12] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id6] = (((temp[id1] + temp[id2] + temp[id3] + temp[id4] + temp[id5] + temp[id6] + temp[id9] + temp[id11] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id7] = (((temp[id2] + temp[id3] + temp[id5] + temp[id6] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        data[id8] = (((temp[id0] + temp[id1] + temp[id3] + temp[id4] + temp[id5] + temp[id8] + temp[id10] + temp[id12]) % PlainMod) + PlainMod) % PlainMod;
        data[id9] = (((temp[id0] + temp[id2] + temp[id3] + temp[id7] + temp[id8] + temp[id9] + temp[id11] + temp[id12] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id10] = (((temp[id0] + temp[id2] + temp[id5] + temp[id7] + temp[id8] + temp[id9] + temp[id10] + temp[id13] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id11] = (((temp[id0] + temp[id6] + temp[id9] + temp[id10] + temp[id14] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id12] = (((temp[id0] + temp[id1] + temp[id4] + temp[id5] + temp[id7] + temp[id8] + temp[id12] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id13] = (((temp[id3] + temp[id4] + temp[id6] + temp[id7] + temp[id8] + temp[id11] + temp[id12] + temp[id13] + temp[id15]) % PlainMod) + PlainMod) % PlainMod;
        data[id14] = (((temp[id1] + temp[id3] + temp[id4] + temp[id6] + temp[id9] + temp[id10] + temp[id11] + temp[id12] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
        data[id15] = (((temp[id2] + temp[id4] + temp[id10] + temp[id11] + temp[id13] + temp[id14]) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::Sbox_6(vector<long> &data)
{
    vector<long> temp = data;
    // 一次迭代，(x0,x1,x2,x3)——> (x0, -x0*x1+x1+x0*x2+x0*x3, x1+x0*x2-x2+x0*x3, -x0*x1-x2-x0*x3+x3)

    for (int i = 0; i < data.size(); i += 4)
    {
        long t01 = temp[i] * temp[i + 1];
        long t02 = temp[i] * temp[i + 2];
        long t03 = temp[i] * temp[i + 3];
        long t02_03_1 = t02 + t03 + temp[i + 1];

        data[i + 1] = (((t02_03_1 - t01) % PlainMod) + PlainMod) % PlainMod;
        data[i + 2] = (((t02_03_1 - temp[i + 2]) % PlainMod) + PlainMod) % PlainMod;
        data[i + 3] = (((temp[i + 3] - t01 - temp[i + 2] - t03) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::M36_7(vector<long> &data)
{
    std::vector<long> temp = data;
    data[0] = (((2 * temp[0] + 2 * temp[1] + temp[2] + 2 * temp[3] + temp[4] + temp[5] + temp[8] + temp[9] + temp[10] + temp[11] + temp[14] + 2 * temp[15] + temp[16] + 2 * temp[17] + 2 * temp[18] + 2 * temp[22] + 2 * temp[24] + temp[26] + 2 * temp[27] + temp[28] + temp[29] + temp[30] + temp[31] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
    data[1] = (((temp[0] + temp[2] + temp[3] + temp[4] + 2 * temp[6] + 2 * temp[7] + 2 * temp[9] + temp[10] + temp[11] + 2 * temp[13] + temp[14] + temp[15] + temp[17] + 2 * temp[20] + 2 * temp[21] + temp[22] + temp[23] + temp[24] + 2 * temp[25] + temp[26] + 2 * temp[27] + temp[29] + 2 * temp[30] + temp[32] + 2 * temp[33] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[2] = (((temp[1] + temp[2] + temp[3] + temp[6] + 2 * temp[8] + 2 * temp[10] + 2 * temp[11] + temp[12] + temp[16] + temp[17] + 2 * temp[20] + temp[23] + temp[26] + temp[27] + temp[29] + temp[31] + 2 * temp[32]) % PlainMod) + PlainMod) % PlainMod;
    data[3] = (((temp[1] + temp[5] + temp[6] + temp[7] + temp[9] + temp[10] + 2 * temp[12] + temp[14] + temp[15] + 2 * temp[17] + 2 * temp[19] + temp[20] + temp[21] + temp[23] + 2 * temp[25] + 2 * temp[26] + temp[28] + temp[29] + temp[30] + 2 * temp[31] + temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[4] = (((temp[1] + 2 * temp[4] + 2 * temp[5] + temp[6] + 2 * temp[7] + temp[8] + temp[9] + temp[12] + temp[13] + temp[14] + temp[15] + temp[18] + 2 * temp[19] + temp[20] + 2 * temp[21] + 2 * temp[22] + 2 * temp[26] + 2 * temp[28] + temp[30] + 2 * temp[31] + temp[32] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[5] = (((temp[0] + 2 * temp[1] + temp[3] + temp[4] + temp[6] + temp[7] + temp[8] + 2 * temp[10] + 2 * temp[11] + 2 * temp[13] + temp[14] + temp[15] + 2 * temp[17] + temp[18] + temp[19] + temp[21] + 2 * temp[24] + 2 * temp[25] + temp[26] + temp[27] + temp[28] + 2 * temp[29] + temp[30] + 2 * temp[31] + temp[33] + 2 * temp[34]) % PlainMod) + PlainMod) % PlainMod;
    data[6] = (((2 * temp[0] + temp[5] + temp[6] + temp[7] + temp[10] + 2 * temp[12] + 2 * temp[14] + 2 * temp[15] + temp[16] + temp[20] + temp[21] + 2 * temp[24] + temp[27] + temp[30] + temp[31] + temp[33] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[7] = (((temp[2] + 2 * temp[3] + temp[5] + temp[9] + temp[10] + temp[11] + temp[13] + temp[14] + 2 * temp[16] + temp[18] + temp[19] + 2 * temp[21] + 2 * temp[23] + temp[24] + temp[25] + temp[27] + 2 * temp[29] + 2 * temp[30] + temp[32] + temp[33] + temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[8] = (((temp[0] + temp[1] + temp[2] + temp[3] + temp[5] + 2 * temp[8] + 2 * temp[9] + temp[10] + 2 * temp[11] + temp[12] + temp[13] + temp[16] + temp[17] + temp[18] + temp[19] + temp[22] + 2 * temp[23] + temp[24] + 2 * temp[25] + 2 * temp[26] + 2 * temp[30] + 2 * temp[32] + temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[9] = (((temp[1] + 2 * temp[2] + temp[4] + 2 * temp[5] + temp[7] + temp[8] + temp[10] + temp[11] + temp[12] + 2 * temp[14] + 2 * temp[15] + 2 * temp[17] + temp[18] + temp[19] + 2 * temp[21] + temp[22] + temp[23] + temp[25] + 2 * temp[28] + 2 * temp[29] + temp[30] + temp[31] + temp[32] + 2 * temp[33] + temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[10] = (((temp[1] + temp[3] + 2 * temp[4] + temp[9] + temp[10] + temp[11] + temp[14] + 2 * temp[16] + 2 * temp[18] + 2 * temp[19] + temp[20] + temp[24] + temp[25] + 2 * temp[28] + temp[31] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[11] = (((temp[0] + temp[1] + temp[2] + 2 * temp[3] + temp[6] + 2 * temp[7] + temp[9] + temp[13] + temp[14] + temp[15] + temp[17] + temp[18] + 2 * temp[20] + temp[22] + temp[23] + 2 * temp[25] + 2 * temp[27] + temp[28] + temp[29] + temp[31] + 2 * temp[33] + 2 * temp[34]) % PlainMod) + PlainMod) % PlainMod;
    data[12] = (((2 * temp[0] + temp[2] + 2 * temp[3] + temp[4] + temp[5] + temp[6] + temp[7] + temp[9] + 2 * temp[12] + 2 * temp[13] + temp[14] + 2 * temp[15] + temp[16] + temp[17] + temp[20] + temp[21] + temp[22] + temp[23] + temp[26] + 2 * temp[27] + temp[28] + 2 * temp[29] + 2 * temp[30] + 2 * temp[34]) % PlainMod) + PlainMod) % PlainMod;
    data[13] = (((temp[0] + 2 * temp[1] + temp[2] + 2 * temp[3] + temp[5] + 2 * temp[6] + temp[8] + 2 * temp[9] + temp[11] + temp[12] + temp[14] + temp[15] + temp[16] + 2 * temp[18] + 2 * temp[19] + 2 * temp[21] + temp[22] + temp[23] + 2 * temp[25] + temp[26] + temp[27] + temp[29] + 2 * temp[32] + 2 * temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[14] = (((temp[2] + temp[3] + temp[5] + temp[7] + 2 * temp[8] + temp[13] + temp[14] + temp[15] + temp[18] + 2 * temp[20] + 2 * temp[22] + 2 * temp[23] + temp[24] + temp[28] + temp[29] + 2 * temp[32] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[15] = (((2 * temp[1] + 2 * temp[2] + temp[4] + temp[5] + temp[6] + 2 * temp[7] + temp[10] + 2 * temp[11] + temp[13] + temp[17] + temp[18] + temp[19] + temp[21] + temp[22] + 2 * temp[24] + temp[26] + temp[27] + 2 * temp[29] + 2 * temp[31] + temp[32] + temp[33] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[16] = (((2 * temp[2] + 2 * temp[4] + temp[6] + 2 * temp[7] + temp[8] + temp[9] + temp[10] + temp[11] + temp[13] + 2 * temp[16] + 2 * temp[17] + temp[18] + 2 * temp[19] + temp[20] + temp[21] + temp[24] + temp[25] + temp[26] + temp[27] + temp[30] + 2 * temp[31] + temp[32] + 2 * temp[33] + 2 * temp[34]) % PlainMod) + PlainMod) % PlainMod;
    data[17] = (((2 * temp[0] + 2 * temp[1] + temp[2] + temp[3] + temp[4] + 2 * temp[5] + temp[6] + 2 * temp[7] + temp[9] + 2 * temp[10] + temp[12] + 2 * temp[13] + temp[15] + temp[16] + temp[18] + temp[19] + temp[20] + 2 * temp[22] + 2 * temp[23] + 2 * temp[25] + temp[26] + temp[27] + 2 * temp[29] + temp[30] + temp[31] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
    data[18] = (((2 * temp[0] + temp[3] + temp[6] + temp[7] + temp[9] + temp[11] + 2 * temp[12] + temp[17] + temp[18] + temp[19] + temp[22] + 2 * temp[24] + 2 * temp[26] + 2 * temp[27] + temp[28] + temp[32] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
    data[19] = (((temp[0] + temp[1] + temp[3] + 2 * temp[5] + 2 * temp[6] + temp[8] + temp[9] + temp[10] + 2 * temp[11] + temp[14] + 2 * temp[15] + temp[17] + temp[21] + temp[22] + temp[23] + temp[25] + temp[26] + 2 * temp[28] + temp[30] + temp[31] + 2 * temp[33] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[20] = (((temp[0] + 2 * temp[1] + 2 * temp[2] + 2 * temp[6] + 2 * temp[8] + temp[10] + 2 * temp[11] + temp[12] + temp[13] + temp[14] + temp[15] + temp[17] + 2 * temp[20] + 2 * temp[21] + temp[22] + 2 * temp[23] + temp[24] + temp[25] + temp[28] + temp[29] + temp[30] + temp[31] + temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[21] = (((temp[1] + 2 * temp[4] + 2 * temp[5] + temp[6] + temp[7] + temp[8] + 2 * temp[9] + temp[10] + 2 * temp[11] + temp[13] + 2 * temp[14] + temp[16] + 2 * temp[17] + temp[19] + temp[20] + temp[22] + temp[23] + temp[24] + 2 * temp[26] + 2 * temp[27] + 2 * temp[29] + temp[30] + temp[31] + 2 * temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[22] = (((temp[0] + temp[1] + 2 * temp[4] + temp[7] + temp[10] + temp[11] + temp[13] + temp[15] + 2 * temp[16] + temp[21] + temp[22] + temp[23] + temp[26] + 2 * temp[28] + 2 * temp[30] + 2 * temp[31] + temp[32]) % PlainMod) + PlainMod) % PlainMod;
    data[23] = (((2 * temp[1] + 2 * temp[3] + temp[4] + temp[5] + temp[7] + 2 * temp[9] + 2 * temp[10] + temp[12] + temp[13] + temp[14] + 2 * temp[15] + temp[18] + 2 * temp[19] + temp[21] + temp[25] + temp[26] + temp[27] + temp[29] + temp[30] + 2 * temp[32] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[24] = (((temp[2] + 2 * temp[3] + temp[4] + 2 * temp[5] + 2 * temp[6] + 2 * temp[10] + 2 * temp[12] + temp[14] + 2 * temp[15] + temp[16] + temp[17] + temp[18] + temp[19] + temp[21] + 2 * temp[24] + 2 * temp[25] + temp[26] + 2 * temp[27] + temp[28] + temp[29] + temp[32] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[25] = (((2 * temp[1] + temp[2] + temp[3] + temp[5] + 2 * temp[8] + 2 * temp[9] + temp[10] + temp[11] + temp[12] + 2 * temp[13] + temp[14] + 2 * temp[15] + temp[17] + 2 * temp[18] + temp[20] + 2 * temp[21] + temp[23] + temp[24] + temp[26] + temp[27] + temp[28] + 2 * temp[30] + 2 * temp[31] + 2 * temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[26] = (((temp[0] + temp[4] + temp[5] + 2 * temp[8] + temp[11] + temp[14] + temp[15] + temp[17] + temp[19] + 2 * temp[20] + temp[25] + temp[26] + temp[27] + temp[30] + 2 * temp[32] + 2 * temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[27] = (((2 * temp[0] + temp[2] + temp[3] + 2 * temp[5] + 2 * temp[7] + temp[8] + temp[9] + temp[11] + 2 * temp[13] + 2 * temp[14] + temp[16] + temp[17] + temp[18] + 2 * temp[19] + temp[22] + 2 * temp[23] + temp[25] + temp[29] + temp[30] + temp[31] + temp[33] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
    data[28] = (((temp[0] + temp[1] + temp[2] + temp[3] + temp[6] + 2 * temp[7] + temp[8] + 2 * temp[9] + 2 * temp[10] + 2 * temp[14] + 2 * temp[16] + temp[18] + 2 * temp[19] + temp[20] + temp[21] + temp[22] + temp[23] + temp[25] + 2 * temp[28] + 2 * temp[29] + temp[30] + 2 * temp[31] + temp[32] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
    data[29] = (((2 * temp[1] + temp[2] + temp[3] + 2 * temp[5] + temp[6] + temp[7] + temp[9] + 2 * temp[12] + 2 * temp[13] + temp[14] + temp[15] + temp[16] + 2 * temp[17] + temp[18] + 2 * temp[19] + temp[21] + 2 * temp[22] + temp[24] + 2 * temp[25] + temp[27] + temp[28] + temp[30] + temp[31] + temp[32] + 2 * temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[30] = (((2 * temp[0] + 2 * temp[2] + 2 * temp[3] + temp[4] + temp[8] + temp[9] + 2 * temp[12] + temp[15] + temp[18] + temp[19] + temp[21] + temp[23] + 2 * temp[24] + temp[29] + temp[30] + temp[31] + temp[34]) % PlainMod) + PlainMod) % PlainMod;
    data[31] = (((temp[1] + temp[2] + 2 * temp[4] + temp[6] + temp[7] + 2 * temp[9] + 2 * temp[11] + temp[12] + temp[13] + temp[15] + 2 * temp[17] + 2 * temp[18] + temp[20] + temp[21] + temp[22] + 2 * temp[23] + temp[26] + 2 * temp[27] + temp[29] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[32] = (((temp[0] + temp[1] + temp[4] + temp[5] + temp[6] + temp[7] + temp[10] + 2 * temp[11] + temp[12] + 2 * temp[13] + 2 * temp[14] + 2 * temp[18] + 2 * temp[20] + temp[22] + 2 * temp[23] + temp[24] + temp[25] + temp[26] + temp[27] + temp[29] + 2 * temp[32] + 2 * temp[33] + temp[34] + 2 * temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[33] = (((temp[0] + 2 * temp[2] + 2 * temp[3] + 2 * temp[5] + temp[6] + temp[7] + 2 * temp[9] + temp[10] + temp[11] + temp[13] + 2 * temp[16] + 2 * temp[17] + temp[18] + temp[19] + temp[20] + 2 * temp[21] + temp[22] + 2 * temp[23] + temp[25] + 2 * temp[26] + temp[28] + 2 * temp[29] + temp[31] + temp[32] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[34] = (((temp[2] + 2 * temp[4] + 2 * temp[6] + 2 * temp[7] + temp[8] + temp[12] + temp[13] + 2 * temp[16] + temp[19] + temp[22] + temp[23] + temp[25] + temp[27] + 2 * temp[28] + temp[33] + temp[34] + temp[35]) % PlainMod) + PlainMod) % PlainMod;
    data[35] = (((temp[1] + temp[2] + temp[3] + temp[5] + temp[6] + 2 * temp[8] + temp[10] + temp[11] + 2 * temp[13] + 2 * temp[15] + temp[16] + temp[17] + temp[19] + 2 * temp[21] + 2 * temp[22] + temp[24] + temp[25] + temp[26] + 2 * temp[27] + temp[30] + 2 * temp[31] + temp[33]) % PlainMod) + PlainMod) % PlainMod;
}
void YusP::Sbox_7(vector<long> &data)
{
    vector<long> temp = data;
    // 一次迭代，(x0,x1,x2,x3)——> (x0, -x0*x1+x1+x0*x2+x0*x3, x1+x0*x2-x2+x0*x3, -x0*x1-x2-x0*x3+x3)

    for (int i = 0; i < data.size(); i += 4)
    {
        long t01 = temp[i] * temp[i + 1];
        long t02 = temp[i] * temp[i + 2];
        long t03 = temp[i] * temp[i + 3];
        long t02_03_1 = t02 + t03 + temp[i + 1];

        data[i + 1] = (((t02_03_1 - t01) % PlainMod) + PlainMod) % PlainMod;
        data[i + 2] = (((t02_03_1 - temp[i + 2]) % PlainMod) + PlainMod) % PlainMod;
        data[i + 3] = (((temp[i + 3] - t01 - temp[i + 2] - t03) % PlainMod) + PlainMod) % PlainMod;
    }
}
void YusP::M18_8(vector<long> &data)
{
    std::vector<long> temp = data;
    data[0] = (((temp[2] + 2 * temp[3] + temp[4] + temp[5] + 2 * temp[6] + 2 * temp[8] + 2 * temp[9] + temp[11] + temp[13] + 2 * temp[14] + 2 * temp[15] + temp[16]) % PlainMod) + PlainMod) % PlainMod;
    data[1] = (((2 * temp[1] + temp[2] + temp[3] + 2 * temp[5] + 2 * temp[6] + 2 * temp[7] + temp[9] + temp[10] + 2 * temp[11] + 2 * temp[12] + 2 * temp[13] + 2 * temp[14] + temp[15] + 2 * temp[16]) % PlainMod) + PlainMod) % PlainMod;
    data[2] = (((2 * temp[0] + temp[1] + temp[2] + 2 * temp[4] + temp[5] + 2 * temp[6] + 2 * temp[7] + temp[8] + temp[9] + temp[11] + temp[12] + 2 * temp[13] + 2 * temp[15] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[3] = (((2 * temp[0] + temp[1] + temp[5] + 2 * temp[6] + temp[7] + temp[8] + 2 * temp[9] + 2 * temp[11] + 2 * temp[12] + temp[14] + temp[16] + 2 * temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[4] = (((temp[0] + 2 * temp[1] + 2 * temp[4] + temp[5] + temp[6] + 2 * temp[8] + 2 * temp[9] + 2 * temp[10] + temp[12] + temp[13] + 2 * temp[14] + 2 * temp[15] + 2 * temp[16] + 2 * temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[5] = (((2 * temp[0] + temp[2] + 2 * temp[3] + temp[4] + temp[5] + 2 * temp[7] + temp[8] + 2 * temp[9] + 2 * temp[10] + temp[11] + temp[12] + temp[14] + temp[15] + 2 * temp[16]) % PlainMod) + PlainMod) % PlainMod;
    data[6] = (((temp[1] + 2 * temp[2] + 2 * temp[3] + temp[4] + temp[8] + 2 * temp[9] + temp[10] + temp[11] + 2 * temp[12] + 2 * temp[14] + 2 * temp[15] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[7] = (((2 * temp[0] + 2 * temp[1] + 2 * temp[2] + temp[3] + 2 * temp[4] + 2 * temp[7] + temp[8] + temp[9] + 2 * temp[11] + 2 * temp[12] + 2 * temp[13] + temp[15] + temp[16] + 2 * temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[8] = (((temp[0] + 2 * temp[1] + 2 * temp[3] + temp[5] + 2 * temp[6] + temp[7] + temp[8] + 2 * temp[10] + temp[11] + 2 * temp[12] + 2 * temp[13] + temp[14] + temp[15] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[9] = (((2 * temp[0] + temp[2] + temp[4] + 2 * temp[5] + 2 * temp[6] + temp[7] + temp[11] + 2 * temp[12] + temp[13] + temp[14] + 2 * temp[15] + 2 * temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[10] = (((temp[0] + temp[1] + 2 * temp[2] + 2 * temp[3] + 2 * temp[4] + 2 * temp[5] + temp[6] + 2 * temp[7] + 2 * temp[10] + temp[11] + temp[12] + 2 * temp[14] + 2 * temp[15] + 2 * temp[16]) % PlainMod) + PlainMod) % PlainMod;
    data[11] = (((temp[0] + temp[2] + temp[3] + 2 * temp[4] + 2 * temp[6] + temp[8] + 2 * temp[9] + temp[10] + temp[11] + 2 * temp[13] + temp[14] + 2 * temp[15] + 2 * temp[16] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[12] = (((2 * temp[0] + 2 * temp[2] + 2 * temp[3] + temp[5] + temp[7] + 2 * temp[8] + 2 * temp[9] + temp[10] + temp[14] + 2 * temp[15] + temp[16] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[13] = (((2 * temp[0] + 2 * temp[1] + temp[3] + temp[4] + 2 * temp[5] + 2 * temp[6] + 2 * temp[7] + 2 * temp[8] + temp[9] + 2 * temp[10] + 2 * temp[13] + temp[14] + temp[15] + 2 * temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[14] = (((2 * temp[0] + 2 * temp[1] + temp[2] + temp[3] + temp[5] + temp[6] + 2 * temp[7] + 2 * temp[9] + temp[11] + 2 * temp[12] + temp[13] + temp[14] + 2 * temp[16] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[15] = (((2 * temp[0] + temp[1] + temp[2] + 2 * temp[3] + 2 * temp[5] + 2 * temp[6] + temp[8] + temp[10] + 2 * temp[11] + 2 * temp[12] + temp[13] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[16] = (((temp[0] + 2 * temp[2] + 2 * temp[3] + 2 * temp[4] + temp[6] + temp[7] + 2 * temp[8] + 2 * temp[9] + 2 * temp[10] + 2 * temp[11] + temp[12] + 2 * temp[13] + 2 * temp[16] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
    data[17] = (((2 * temp[1] + temp[2] + 2 * temp[3] + 2 * temp[4] + temp[5] + temp[6] + temp[8] + temp[9] + 2 * temp[10] + 2 * temp[12] + temp[14] + 2 * temp[15] + temp[16] + temp[17]) % PlainMod) + PlainMod) % PlainMod;
}
void YusP::Sbox_8(vector<long> &data)
{
    vector<long> temp = data;
    // 一次迭代，(x0,x1,x2)——>(x0,x0*x2+x1,-x0*x1+x0*x2+x2)
    for (int i = 0; i < data.size(); i += 3)
    {
        data[i + 1] = (temp[i] * temp[i + 2] + temp[i + 1]) % PlainMod;
        data[i + 2] = (((-temp[i] * temp[i + 1] + temp[i] * temp[i + 2] + temp[i + 2]) % PlainMod) + PlainMod) % PlainMod;
    }
}