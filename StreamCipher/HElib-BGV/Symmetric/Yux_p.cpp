
#include <cstring>
#include "Yux_p.hpp"

void YuxP::Sbox(vector<long> &data)
{
  for (int i = 0; i < data.size(); i += 4)
  {
    for (int j = 0; j < 2; j++)
    {
      long e0 = data[i];
      long e1 = data[i + 1];
      long e2 = data[i + 2];
      long e3 = data[i + 3];

      long temp2 = (e1 * e2 + e0 + e3 + RoundConstant) % PlainMod;

      long temp3 = (e2 * e3 + e1 + temp2 + RoundConstant) % PlainMod;

      data[i] = e2;
      data[i + 1] = e3;
      data[i + 2] = temp2;
      data[i + 3] = temp3;
    }
  }
}
void YuxP::Sbox2(vector<long> &data)
{
  for (int i = 0; i < data.size(); i += 4)
  {
    for (int j = 0; j < 2; j++)
    {
      long e0 = data[i];
      long e1 = data[i + 1];
      long e2 = data[i + 2];
      long e3 = data[i + 3];

      long temp2 = ((e1 + e2) * (e1 + e2) + e0 + e3 + RoundConstant) % PlainMod;

      long temp3 = ((e2 + e3) * (e2 + e3) + e1 + temp2 + RoundConstant) % PlainMod;

      data[i] = e2;
      data[i + 1] = e3;
      data[i + 2] = temp2;
      data[i + 3] = temp3;
    }
  }
}
void YuxP::Linear(vector<long> &data)
{
  int j;
  long temp[16];
  for (j = 0; j < 16; j++)
  {
    temp[j] = data[j];
  }
  // linear Layer
  for (j = 0; j < 16; j++)
  {
    data[j] = (temp[j] + temp[(j + 3) % 16] + temp[(j + 4) % 16] + temp[(j + 8) % 16] + temp[(j + 9) % 16] + temp[(j + 12) % 16] + temp[(j + 14) % 16]) % PlainMod;
  }
}