#include <cstring>
#include <iostream>
//#include "Yu2x-8.h"
// The round number is ROUND =12

static const unsigned char roundConstant = 0xCD;
static const unsigned char roundConstant_4bit = 0xD; // x^3+x^2+1

// This function adds the round key to state.
// The round key is added to the state by an XOR function.
void addRoundKey(unsigned char state[], unsigned char RoundKey[], int round)
{
  int i, blockByte = 16;
  for(i=0;i<(blockByte);i++)
    {
      int key_id = blockByte*round+i;
	    state[i] ^= RoundKey[key_id];
    }
}

unsigned char mul(unsigned char a, unsigned char b) {
    unsigned char p = 0;
    unsigned char counter;
    unsigned char hi_bit_set;
    for(counter = 0; counter < 8; counter++) {
          if((b & 1) == 1) 
                p ^= a;
          hi_bit_set = (a & 0x80);
          a <<= 1;
          if(hi_bit_set == 0x80) 
                a ^= 0x1b;          
          b >>= 1;
      }
      return p;
	}

void decSboxFi(unsigned char state[], int begin)
{
  unsigned char c0=state[begin];
  unsigned char c1=state[begin+1];
  unsigned char c2=state[begin+2];
  unsigned char c3=state[begin+3];

  unsigned char temp = mul(c1,c2) ^ c0 ^ c3 ^ roundConstant;

  state[begin] = c1;
  state[begin+1] = c2;
  state[begin+2] = c3;
  state[begin+3] = temp;
}


void decLinearLayer(unsigned char in[16])
{
    int j;  
    unsigned char temp[16];
    for(j=0;j<16;j++){
      temp[j] = in[j];
    }  
    // invert linear Layer
    for(j=0; j<16; j++){
      in[j] =    temp[j]
                ^ temp[(j+3)%16]
                ^ temp[(j+4)%16]
                ^ temp[(j+8)%16]
                ^ temp[(j+9)%16]
                ^ temp[(j+12)%16]
                ^ temp[(j+14)%16];
      }
}

// Cipher is the main function that encrypts the PlainText.
void decryption(unsigned char out[], unsigned char in[], unsigned char RoundKey[], int Nr)
{
  int i,round;
  // initial key addition
  // Add Round key before first round

  for(i=0; i<16; i++){
      out[i] = in[i];
  }

  addRoundKey(out, RoundKey, 0);

  // There will be Nr rounds.
  // The first Nr-1 rounds are identical.
  // These Nr-1 rounds are executed in the loop below.
  for(round=1;round<Nr;round++){
    // S Layer -- 4 sbox
    for(i=0; i<4; i++){
      decSboxFi(out, i*4);
      decSboxFi(out, i*4);
      decSboxFi(out, i*4);
      decSboxFi(out, i*4);
    }

    // linear Layer
    decLinearLayer(out);
    addRoundKey(out, RoundKey, round);
  }

  // The last round is given below.
  // Linear layer is not here in the last round
  {
    for(i=0; i<4; i++){
      decSboxFi(out, i*4);
      decSboxFi(out, i*4);
      decSboxFi(out, i*4);
      decSboxFi(out, i*4);
    }
    addRoundKey(out, RoundKey, Nr);
  }

}

void encSboxFi(unsigned char state[], int begin)
{
  unsigned char c0=state[begin];
  unsigned char c1=state[begin+1];
  unsigned char c2=state[begin+2];
  unsigned char c3=state[begin+3];

  unsigned char temp = mul(c0,c1) ^ c2 ^ c3 ^ roundConstant;

  state[begin] = temp;
  state[begin+1] = c0;
  state[begin+2] = c1;
  state[begin+3] = c2;
}


void encLinearLayer(unsigned char in[16])
{
    unsigned char temp[16];
    int j;
    char in_all =0;
    for(j=0;j<16;j++){
      temp[j] = in[j];
      in_all ^= in[j];
    }  
    // linear Layer
    for(j=0; j<16; j++){
      in[j] =   in_all 
                ^ temp[(j)%16]
                ^ temp[(j+4)%16]
                ^ temp[(j+9)%16]
                ^ temp[(j+10)%16]
                ^ temp[(j+11)%16];                 
      }
}


// Cipher is the main function that encrypts the PlainText.
void encryption(unsigned char out[], unsigned char in[], unsigned char RoundKey[], int Nr)
{
  int i,round;
  // initial key addition
  // Add Round key before first round

  for(i=0; i<16; i++){
      out[i] = in[i];
  }
  addRoundKey(out, RoundKey, 0);

  // There will be Nr rounds.
  // The first Nr-1 rounds are identical.
  // These Nr-1 rounds are executed in the loop below.
  for(round=1;round<Nr;round++){
    // S Layer -- 4 sbox
    for(i=0; i<4; i++){
      encSboxFi(out, i*4);
      encSboxFi(out, i*4);
      encSboxFi(out, i*4);
      encSboxFi(out, i*4);
    }
    // linear Layer
    encLinearLayer(out);
    addRoundKey(out, RoundKey, round);
  }

  // The last round is given below.
  // Linear layer is not here in the last round
  {
    for(i=0; i<4; i++){
      encSboxFi(out, i*4);
      encSboxFi(out, i*4);
      encSboxFi(out, i*4);
      encSboxFi(out, i*4);
    }
    addRoundKey(out, RoundKey, Nr);
  }
}

// This function produces Nb(Nr+1) round keys.
// The round keys are used in each round to encrypt the states.
void constantForKey(unsigned char RC[56][4], long round)
{
    // Nr is the round number
    // Nk is the number of 64-bit Feistel works in the key 4bytes each round.
    int Nr = round;
    
    // The first round key is the key itself. [0-4]
    int i,j,k;
    for(i=0;i<56;i++)
    {
      unsigned char tmp[4];
      for(j=0; j<4; j++)
      {
        tmp[j] = 4*i+j+1;
        // printf("%02x ",tmp[j]);
      }
      
      //Sf(tmp)
      encSboxFi(tmp, 0);
      encSboxFi(tmp, 0);
      encSboxFi(tmp, 0);
      encSboxFi(tmp, 0);

      // mcsry RC[i]=tmp;
      memcpy(&RC[i], &tmp[0], 4);
    }
}

// array a, length = l, <<<3
void rotation(unsigned char *a, int l,int r){
	unsigned char temp[l];
	for (int i = 0; i < l; i++){
		temp[i] = a[(i+r)%l];
	}
  for (int i =0; i< l; i++)
  {
    a[i]=temp[i];
  }
}

// This function produces Nb(Nr+1) round keys.
// The round keys are used in each round to encrypt the states.
int KeyExpansion(unsigned char RoundKey[], int round, int blockByte,  unsigned char Key[])
{
    // Nr is the round number
    // Nk is the number of 64-bit Feistel works in the key 4bytes each round.
    int Nr = round;
    int Nk = blockByte; // 128/8= 16
    
    // The first round key is the key itself. [0-4]
    int i,j,k;
    for(i=0;i<Nk;i++)
    {
        RoundKey[i]=Key[i];
    }

    unsigned char RC[56][4];
    constantForKey(RC, round);

    int x4id = 16;
    // the rest round key is the key it self.
    for(i=0; i<Nr*4; i++) //[1-9]
    { 

      int x0id = i*4;
      int x1id = x0id+4;
      int x2id = x1id+4;
      int x3id = x2id+4;
      // x4 = x1+x2+x3
      unsigned char x4[4];
      for(j=0;j<4;j++) 
      {
        x4[j] = RoundKey[x1id+j];
        x4[j] ^= RoundKey[x2id+j];
        x4[j] ^= RoundKey[x3id+j]; 
      }
      // <<<3
      rotation(x4, 4, 3);

      // x4=Sf(x4)
      encSboxFi(x4, 0);
      encSboxFi(x4, 0);
      encSboxFi(x4, 0);
      encSboxFi(x4, 0);

      // RK[i*4+16 ~ i*4+20] =x0+x4+RC[i]
      for(j=0; j<4; j++)
      {
        x4[j] ^= RoundKey[x0id+j];
        x4[j] ^= RC[i][j];
        // printf("RC[ij]: %02x", RC[i][j]);
        RoundKey[x4id+j] = x4[j];
      }
      x4id +=4;
    }
  return Nr+1;
}

//  Change to Decrypt roundkey
//  1. roundkey in reversed order
//  2. Except the first and the last roundkey, 
//     others are the roundkey with inverse transformation of the linear function.
void decRoundKey(unsigned char RoundKey_invert[], unsigned char RoundKey[], long Nr,long blockByte)
{
    int Nk = blockByte; // 128/8= 16
    if (1>Nr) return; // no data/key
    int i;
    // the first and last Decrypt round key donot need linear exchange
    long Begin = 0;
    long invertBegin = Nk*(Nr);
    memcpy(&RoundKey_invert[invertBegin], &RoundKey[Begin], Nk);

    // dec round key = L^-1(key)
    for(i = 1; i<Nr; i++)
    {
        Begin = Nk*i;
        invertBegin = Nk*(Nr-i);
        unsigned char tempKey[Nk];
        memcpy(&tempKey[0], &RoundKey[Begin], Nk);
        decLinearLayer(tempKey);
        memcpy(&RoundKey_invert[invertBegin], &tempKey[0], Nk);
    }
    // the first and last Decrypt round key donot need linear exchange
    Begin = Nk*(Nr);
    invertBegin = 0;
    memcpy(&RoundKey_invert[invertBegin], &RoundKey[Begin], Nk);

    // printf("RoundKey_invert---:\n");
    // for(int r=0;r<nRoundKeys; r++)
    // {
    //   for (int d=0; d< Nk; d++)
    //   {
    //     cout<<d;
    //     printf(". %02x ;",RoundKey_invert[r*Nk+d]);
    //   }
    //   cout<< "\n";
    // }
    // printf("RoundKey_invert---END!\n");

}

int main(){
  // 明文
    unsigned char plain[16] = {0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34};
    // 初始密钥
    unsigned char key[16] = {0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};
// 密钥扩展
#if 1
    unsigned char RoundKey[208]; // 16*13
    int ROUND = 12;              // Yux加密轮数，含最后一轮，不含白化轮
    int blockByte = 16;
    int nRoundKeys = KeyExpansion(RoundKey, ROUND, blockByte, key);
    // 将RoundKey转为二维数组RoundKey2D
    unsigned char RoundKey2D[13][16];
    for (int i = 0; i < 13; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            RoundKey2D[i][j] = RoundKey[i * 16 + j];
        }
    }
#endif

// 明文Yux加密
#if 1
    // 用轮密钥加密明文
    unsigned char state[16];
    // 令state=plain
    for (int i = 0; i < 16; i++)
    {
        state[i] = plain[i];
    }
    // 定义一个二维数组，来存储Yux每一轮加密后的state,12行16列。不是tlwe密文
    unsigned char Yux_state[3 * ROUND - 1][16];
    std::cout << "Yux加密开始" << std::endl;
    // 白化轮
    addRoundKey(state, RoundKey, 0); // 0~12
    // 将state赋值给Yux_state
    for (int i = 0; i < 16; i++)
    {
        Yux_state[0][i] = state[i];
    }
    // 11 轮常规 Yux 轮
    for (int r = 1; r < ROUND; r++)
    {
        // S Layer -- 4 sbox
        for (int i = 0; i < 4; i++)
        {
            encSboxFi(state, i * 4);
            encSboxFi(state, i * 4);
            encSboxFi(state, i * 4);
            encSboxFi(state, i * 4);
        }
        for (int i = 0; i < 16; i++)
        {
            Yux_state[3 * r - 2][i] = state[i];
        }
        // linear Layer
        encLinearLayer(state);
        for (int i = 0; i < 16; i++)
        {
            Yux_state[3 * r - 1][i] = state[i];
        }
        addRoundKey(state, RoundKey, r);
        for (int i = 0; i < 16; i++)
        {
            Yux_state[3 * r][i] = state[i];
        }
    }
    // 最后一轮
    for (int i = 0; i < 4; i++)
    {
        encSboxFi(state, i * 4);
        encSboxFi(state, i * 4);
        encSboxFi(state, i * 4);
        encSboxFi(state, i * 4);
    }
    for (int i = 0; i < 16; i++)
    {
        Yux_state[34][i] = state[i];
    }
    addRoundKey(state, RoundKey, 12);
    // 输出state
    for (int i = 0; i < 16; i++)
    {
        std::cout << std::hex << (int)state[i] << " ";
    }
    std::cout << std::endl;
    // 用轮密钥加密明文
    unsigned char state2[16];
    encryption(state2, plain, RoundKey, 12);
    // 输出state2
    for (int i = 0; i < 16; i++)
    {
        std::cout << std::hex << (int)state2[i] << " ";
    }
    std::cout << std::endl;
    unsigned char RoundKey_invert[208]; // 16*13
    decRoundKey(RoundKey_invert, RoundKey, 12, 16);
    // 对state2解密
    unsigned char state3[16];
    decryption(state3, state2, RoundKey_invert, 12);
    // 输出state3
    for (int i = 0; i < 16; i++)
    {
        std::cout << std::hex << (int)state3[i] << " ";
    }
    std::cout << std::endl;
    // 输出plain
    for (int i = 0; i < 16; i++)
    {
        std::cout << std::hex << (int)plain[i] << " ";
    }
    std::cout << std::endl;
#endif
}
// g++ -o Yux2 Yux2.cpp