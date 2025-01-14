
#include <cstring>
#include "Yux_p.hpp"

void YuxP::Sbox(vector<long> &data)
{
  vector<long> temp = data;
  for (int i = 0; i < data.size(); i += 4)
  {
    for (int j = 0; j < 4; j++)
    {
      data[i] = (temp[i + 3]-temp[i+0]*temp[i+1]-temp[i+2]-RoundConstant) % PlainMod;
      data[i + 1] = temp[i+0];
      data[i + 2] = temp[i+1];
      data[i + 3] = temp[i+2];
    }
  }
}
void YuxP::Linear(vector<long> &data)
{
  vector<long> temp = data;
  // linear Layer
  for (int j = 0; j < 16; j++)
  {
    data[j] = (
               ((temp[(j) % 16] + temp[(j + 4) % 16]) * 9363) % PlainMod +
               ((temp[(j + 1) % 16] + temp[(j + 2) % 16] + temp[(j + 3) % 16] + temp[(j + 5) % 16] + temp[(j + 6) % 16] + temp[(j + 7) % 16] + temp[(j + 13) % 16] + temp[(j + 14) % 16] + temp[(j + 15) % 16]) * 53054) % PlainMod +
               ((temp[(j + 8) % 16] + temp[(j + 12) % 16]) * 9362) % PlainMod +
               ((temp[(j + 9) % 16] + temp[(j + 10) % 16] + temp[(j + 11) % 16]) * 53053) % PlainMod
              ) %PlainMod;
  }
}

void YuxP::constantForKey(vector<vector<long>> &RC)
{
  // Nr is the round number
  // Nk is the number of 64-bit Feistel works in the key 4bytes each round.

  // The first round key is the key itself. [0-4]
  int i, j, k;
  for (i = 0; i < 56; i++)
  {
    vector<long> tmp(4);
    for (j = 0; j < 4; j++)
    {
      tmp[j] = 4 * i + j + 1;
      // printf("%04x ",tmp[j]);
    }
    Sbox(tmp);
    // mcsry RC[i]=tmp;
    // memcpy(&RC[i], &tmp[0], 4);
    for (j = 0; j < 4; j++)
    {
      // printf("%04x ",tmp[j]);
      RC[i][j] = tmp[j];
    }
    // cout<< endl;
    // RC[i] = &tmp;
  }
}

// array a, length = l, <<<3
void YuxP::rotation(vector<long> &a, long l, long r)
{
  long temp[l];
  for (int i = 0; i < l; i++)
  {
    temp[i] = a[(i + r) % l];
  }
  for (int i = 0; i < l; i++)
  {
    a[i] = temp[i];
  }
}

// This function produces Nb(Nr+1) round keys.
// The round keys are used in each round to encrypt the states.
void YuxP::KeyExpansion(vector<long> &Key, vector<long> &RoundKey, long Nr, long BlockWords)
{
  // Nr is the round number
  // Nk is the number of 64-bit Feistel works in the key 4bytes each round.

  // The first round key is the key itself. [0-4]
  long Nk = BlockWords;
  for (int i = 0; i < Nk; i++)
  {
    RoundKey[i] = Key[i];
  }

  vector<vector<long>> RC(56, vector<long>(4));
  constantForKey(RC);

  int x4id = 16;
  // the rest round key is the key it self.
  for (int i = 0; i < Nr * 4; i++) //[1-9]
  {
    int x0id = i * 4;
    int x1id = i * 4 + 4;
    int x2id = i * 4 + 8;
    int x3id = i * 4 + 12;
    // x4 = x1+x2+x3
    vector<long> x4(4);
    for (long j = 0; j < 4; j++)
    {
      x4[j] = (RoundKey[x1id + j] + RoundKey[x2id + j] + RoundKey[x3id + j]) % PlainMod;
    }
    // <<<3
    rotation(x4, 4, 3);

    // x4=Sf(x4)
    Sbox(x4);

    // RK[i*4+16 ~ i*4+20] =x0+x4+RC[i]
    for (int j = 0; j < 4; j++)
    {
      x4[j] = (x4[j] + RC[i][j] + RoundKey[x0id + j]) % PlainMod;
      // printf("RC[ij]: %04x ", RC[i][j]);
      RoundKey[x4id + j] = x4[j];
    }
    x4id += 4;
  }
}
