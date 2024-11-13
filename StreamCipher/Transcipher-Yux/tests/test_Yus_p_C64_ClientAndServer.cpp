#include <iostream>
#include <cstring>
#include <stdint.h>
#include <chrono>
#include <vector>
#include <omp.h>
#include <array>
#include <cstdint>
#include <cstring>
#include <atomic>
#include <cmath>
#include <chrono>
#include <fstream>
#include <memory>
#include <random>
#include <climits>

#include <helib/helib.h>
#include <helib/ArgMap.h>
#include <helib/DoubleCRT.h>

#include "random_bit.hpp"
#include "Yus_p.hpp"
#include "tool.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

/************
 * PlainMod = 257, Bytebits = 9, BlockSize = 32*9 = 288, PlainBlock = (Plainmod-1)/2 = 128
 * PlainByte = 65537, Bytebits = 17, BlockSize = 32*17 = 544, PlainBlock = (Plainmod-1)/2 = 32768
 */
constexpr long PlainMod = 257;   // 2^8+1费马素数
constexpr unsigned Bytebits = 9; // 字节比特长度=ceil(log2(PlainMod))
// constexpr long PlainMod = 65537;  // 2^16+1
// constexpr unsigned Bytebits = 17; // 字节比特长度=ceil(log2(PlainMod))
constexpr long BlockByte = 64; // 分组字节长度

constexpr unsigned BlockSize = Bytebits * BlockByte; // 分组比特长度=BlockByte*Bytebits
constexpr unsigned PlainBlock = (PlainMod - 1) / 2;  // 明文分组数

static const long PlainByte = BlockByte * PlainBlock; // 明文字节长度
// PlainByte = nslots
static const long PlainSize = BlockSize * PlainBlock;           // 明文比特长度
static const unsigned NonceSize = 32;                           // Nonce比特长度
static const unsigned Nr = 4;                                   // 轮数
static const long counter_begin = 0;                            // 计数器起始值
static const long counter_end = PlainBlock + counter_begin - 1; // 计数器结束值
static bool Rkflag = 1;                                         // true/1表示乘法，false/0表示加法
static bool Deflag = 1;                                         // true/1表示进行解密，false/0表示不进行
/************************************************************************
  long p;          // plaintext primeplain_mod;
  long m;          // m-th cyclotomic polynomial
  long r;          // Lifting [defualt = 1]
  long bits;       // bits in the ciphertext modulus chain
  long c;          // columns in the key-switching matrix [default=2]
  long d;          // Degree of the field extension [default=1]
  long k;          // Security parameter [default=80]
  long s;          // Minimum number of slots [default=0]
************************************************************************/
static long mValues[][4] = {
//  {    p,      m,    bits,   c}
    {65537, 131072,    1320,   6},
    {  257,    256,     300,  2},
    //{65537,  65536,     853,  17},
};
// p^d = 1 mod m,d=1,slots=\phi(m)/d=\phi(m);m=65536=2^16,\phi(m)=2^15=32768

YusP yusP(PlainMod);

void encodeTo64Ctxt(vector<ZZX> &encData, const vector<long> &data, const EncryptedArray &ea)
{
    long R = data.size() / PlainByte;
    long nCtxt = BlockByte * R;
    long data_length = data.size();
    encData.resize(nCtxt);

    for (long r = 0; r < R; r++)
    {
        for (long i = 0; i < BlockByte; i++)
        {
            vector<long> slots(ea.size(), 0);
            for (long j = 0; j < PlainBlock; j++)
            {
                long byteIdx = j * BlockByte + i + r * PlainByte;
                if (byteIdx < data_length)
                {
                    slots[j] = data[byteIdx];
                }
            }
            ea.encode(encData[r * BlockByte + i], slots);
        }
    }
}
// encodeTo64Ctxt对应的解码
void decodeTo64Ctxt(vector<long> &data, const vector<vector<long>> &encData,
                    const EncryptedArray &ea)
{
    long R = encData.size() / BlockByte;
    long data_length = R * PlainByte;

    data.resize(data_length);

    for (long r = 0; r < R; r++)
    {
        for (long j = 0; j < PlainBlock; j++)
        {
            for (long i = 0; i < BlockByte; i++)
            { // i is the ciphertext number

                // j is the block number in this ctxt
                long byteIdx = j * BlockByte + i + r * PlainByte;
                if (byteIdx < data_length)
                {
                    data[byteIdx] = encData[r * BlockByte + i][j];
                }
            }
        }
    }
}
// 函数：解密并验证密文是否正确，需要解码
// 函数：解密并验证密文是否正确
bool verifyDecryption64(const std::vector<Ctxt> &encryptedVec, const vector<long> &originalVec, const SecKey &secretKey,
                        const EncryptedArray &ea)
{
    auto start_decrypt = std::chrono::steady_clock::now();

    vector<long> decryptedVec;
    std::vector<std::vector<long>> decryptedPolys(encryptedVec.size());

    for (std::size_t i = 0; i < encryptedVec.size(); ++i)
    {
        ea.decrypt(encryptedVec[i], secretKey, decryptedPolys[i]);
    }
    // 解码
    decodeTo64Ctxt(decryptedVec, decryptedPolys, ea);
    // 验证解密结果
    bool isDecryptedVecCorrect = true;
    long len = originalVec.size();
    for (std::size_t i = 0; i < len; ++i)
    {
        if (decryptedVec[i] != originalVec[i])
        {
            std::cout << "Decryption check failed at index " << i << ": expected " << originalVec[i]
                      << ", got " << decryptedVec[i] << std::endl;
            isDecryptedVecCorrect = false;
            break;
        }
    }

    auto end_decrypt = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification succeeded! Time: " << elapsed_seconds.count() << "s\n";

    return isDecryptedVecCorrect;
}
void encryptSymKey(vector<Ctxt> &encryptedSymKey, const vector<long> &SymKey, unique_ptr<PubKey> &pk, EncryptedArray &ea)
{
    long nslots = ea.size();
    // 加密
    encryptedSymKey.resize(BlockByte, Ctxt(*pk));
    for (long i = 0; i < BlockByte; i++)
    { // encrypt the encoded key
        vector<long> slotsData(nslots, SymKey[i]);
        ea.encrypt(encryptedSymKey[i], *pk, slotsData);
    }
}

bool verify_encryptSymKey(vector<Ctxt> &encryptedSymKey, const vector<long> &SymKey, const SecKey &secretKey, EncryptedArray &ea)
{
    vector<long> decryptedSymKey(BlockByte);
    for (long i = 0; i < BlockByte; i++)
    { // encrypt the encoded key
        vector<long> slotsData;
        ea.decrypt(encryptedSymKey[i], secretKey, slotsData);
        decryptedSymKey[i] = slotsData[0];
    }
    bool isDecryptedSymKeyCorrect = true;
    for (long i = 0; i < BlockByte; i++)
    {
        if (decryptedSymKey[i] != SymKey[i])
        {
            std::cout << "Decryption check failed at index " << i << ": expected " << SymKey[i]
                      << ", got " << decryptedSymKey[i] << std::endl;
            isDecryptedSymKeyCorrect = false;
            break;
        }
    }
    if (isDecryptedSymKeyCorrect)
    {
        std::cout << "Decryption check succeeded: Decrypted vector matches original vector." << std::endl;
    }
    else
    {
        std::cout << "Decryption check failed: Decrypted vector does not match original vector." << std::endl;
    }
    return isDecryptedSymKeyCorrect;
}
void HE_MC(vector<Ctxt> &eData)
{
    vector<Ctxt> temp = eData;
    vector<int> index = {0, 1, 2, 3, 4, 5, 6, 7};
    for (int i = 0; i < 8; i++)
    {
        int s = 8 * i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s;
        // 下面是2次赋值，36次加/减法
        eData[id0] = temp[id1];
        eData[id0] += temp[id2];
        eData[id0] += temp[id4];
        eData[id0] += temp[id5];
        eData[id0] += temp[id6];

        eData[id1] = temp[id0];
        eData[id1] += temp[id3];
        eData[id1] += temp[id4];
        eData[id1] -= temp[id7];

        eData[id2] += temp[id1];
        eData[id2] += temp[id3];
        eData[id2] += temp[id4];
        eData[id2] -= temp[id6];
        eData[id2] -= temp[id7];

        eData[id3] += temp[id0];
        eData[id3] += temp[id1];
        eData[id3] += temp[id5];
        eData[id3] -= temp[id6];

        eData[id4] += temp[id0];
        eData[id4] += temp[id1];
        eData[id4] += temp[id2];
        eData[id4] -= temp[id6];
        eData[id4] += temp[id7];

        eData[id5] += temp[id0];
        eData[id5] += temp[id2];
        eData[id5] += temp[id3];
        eData[id5] += temp[id6];

        eData[id6] += temp[id1];
        eData[id6] -= temp[id0];
        eData[id6] -= temp[id2];
        eData[id6] += temp[id3];
        eData[id6] += temp[id4];

        eData[id7] += temp[id2];
        eData[id7] -= temp[id0];
        eData[id7] += temp[id3];
        eData[id7] += temp[id4];
        eData[id7] += temp[id5];
        eData[id7] -= temp[id6];
    }
}
void HE_MR(vector<Ctxt> &eData)
{
    vector<Ctxt> temp = eData;
    vector<int> index = {0, 8, 16, 24, 32, 40, 48, 56};
    for (int i = 0; i < 8; i++)
    {
        int s = i;
        int id0 = index[0] + s, id1 = index[1] + s, id2 = index[2] + s, id3 = index[3] + s, id4 = index[4] + s, id5 = index[5] + s, id6 = index[6] + s, id7 = index[7] + s;
        // 下面是2次赋值，36次加/减法
        eData[id0] = temp[id1];
        eData[id0] += temp[id2];
        eData[id0] += temp[id4];
        eData[id0] += temp[id5];
        eData[id0] += temp[id6];

        eData[id1] = temp[id0];
        eData[id1] += temp[id3];
        eData[id1] += temp[id4];
        eData[id1] -= temp[id7];

        eData[id2] += temp[id1];
        eData[id2] += temp[id3];
        eData[id2] += temp[id4];
        eData[id2] -= temp[id6];
        eData[id2] -= temp[id7];

        eData[id3] += temp[id0];
        eData[id3] += temp[id1];
        eData[id3] += temp[id5];
        eData[id3] -= temp[id6];

        eData[id4] += temp[id0];
        eData[id4] += temp[id1];
        eData[id4] += temp[id2];
        eData[id4] -= temp[id6];
        eData[id4] += temp[id7];

        eData[id5] += temp[id0];
        eData[id5] += temp[id2];
        eData[id5] += temp[id3];
        eData[id5] += temp[id6];

        eData[id6] += temp[id1];
        eData[id6] -= temp[id0];
        eData[id6] -= temp[id2];
        eData[id6] += temp[id3];
        eData[id6] += temp[id4];

        eData[id7] += temp[id2];
        eData[id7] -= temp[id0];
        eData[id7] += temp[id3];
        eData[id7] += temp[id4];
        eData[id7] += temp[id5];
        eData[id7] -= temp[id6];
    }
}
//
void HE_MC_MR(vector<Ctxt> &eData)
{
    vector<Ctxt> temp = eData;
    vector<Ctxt> t(26, eData[0]);
    vector<int> index1 = {0, 1, 2, 3, 4, 5, 6, 7};
    for (int i = 0; i < 8; i++)
    {
        int s = 8 * i;
        int id0 = index1[0] + s, id1 = index1[1] + s, id2 = index1[2] + s, id3 = index1[3] + s, id4 = index1[4] + s, id5 = index1[5] + s, id6 = index1[6] + s, id7 = index1[7] + s;
        // 下面是26次赋值，26次加/减法
        t[0] = temp[id2];
        t[0] += temp[id4];
        t[1] = t[0];
        t[1] -= temp[id6];
        t[2] = t[1];
        t[2] += temp[id7];
        t[3] = temp[id0];
        t[3] += temp[id3];
        t[4] = t[3];
        t[4] += temp[id5];
        t[5] = t[0];
        t[5] += temp[id6];
        t[6] = temp[id1];
        t[6] += temp[id5];
        eData[id0] = t[5];
        eData[id0] += t[6];

        t[8] = t[3];
        t[8] += temp[id4];
        eData[id1] = t[8];
        eData[id1] -= temp[id7];

        t[10] = t[4];
        t[10] += t[5];
        eData[id5] = t[10];
        eData[id5] -= temp[id4];

        t[12] = temp[id1];
        t[12] += temp[id3];
        t[13] = t[1];
        t[13] += t[12];
        eData[id2] = t[13];
        eData[id2] -= temp[id7];

        t[15] = t[4];
        t[15] += temp[id1];
        eData[id3] = t[15];
        eData[id3] -= temp[id6];

        t[17] = t[2];
        t[17] += temp[id0];
        eData[id4] = t[17];
        eData[id4] += temp[id1];

        t[19] = t[12];
        t[19] -= temp[id0];
        t[20] = t[19];
        t[20] -= temp[id2];
        t[21] = t[20];
        t[21] += temp[id4];
        eData[id6] = t[21];
        eData[id6] += temp[id6];

        t[23] = temp[id3];
        t[23] += temp[id5];
        t[24] = t[23];
        t[24] -= temp[id0];
        eData[id7] = t[24];
        eData[id7] += t[2];
    }
    temp = eData;
    vector<int> index2 = {0, 8, 16, 24, 32, 40, 48, 56};
    for (int i = 0; i < 8; i++)
    {
        int s = i;
        int id0 = index2[0] + s, id1 = index2[1] + s, id2 = index2[2] + s, id3 = index2[3] + s, id4 = index2[4] + s, id5 = index2[5] + s, id6 = index2[6] + s, id7 = index2[7] + s;
        // 下面是26次赋值，26次加/减法
        t[0] = temp[id2];
        t[0] += temp[id4];
        t[1] = t[0];
        t[1] -= temp[id6];
        t[2] = t[1];
        t[2] += temp[id7];
        t[3] = temp[id0];
        t[3] += temp[id3];
        t[4] = t[3];
        t[4] += temp[id5];
        t[5] = t[0];
        t[5] += temp[id6];
        t[6] = temp[id1];
        t[6] += temp[id5];
        eData[id0] = t[5];
        eData[id0] += t[6];

        t[8] = t[3];
        t[8] += temp[id4];
        eData[id1] = t[8];
        eData[id1] -= temp[id7];

        t[10] = t[4];
        t[10] += t[5];
        eData[id5] = t[10];
        eData[id5] -= temp[id4];

        t[12] = temp[id1];
        t[12] += temp[id3];
        t[13] = t[1];
        t[13] += t[12];
        eData[id2] = t[13];
        eData[id2] -= temp[id7];

        t[15] = t[4];
        t[15] += temp[id1];
        eData[id3] = t[15];
        eData[id3] -= temp[id6];

        t[17] = t[2];
        t[17] += temp[id0];
        eData[id4] = t[17];
        eData[id4] += temp[id1];

        t[19] = t[12];
        t[19] -= temp[id0];
        t[20] = t[19];
        t[20] -= temp[id2];
        t[21] = t[20];
        t[21] += temp[id4];
        eData[id6] = t[21];
        eData[id6] += temp[id6];

        t[23] = temp[id3];
        t[23] += temp[id5];
        t[24] = t[23];
        t[24] -= temp[id0];
        eData[id7] = t[24];
        eData[id7] += t[2];
    }
}
// Compute the constants for Sbox
void HE_Sbox(vector<Ctxt> &eData)
{
    // #pragma omp parallel for
    for (long j = 0; j < BlockByte; j += 2)
    {
        Ctxt c0(eData[j]);
        Ctxt c1(eData[j + 1]);

        for (int i = 0; i < 2; i++)
        {
            Ctxt temp = c0;
            // c0.frobeniusAutomorph(1);
            c0.square();
            c0 += c1;
            c1 = temp;
        }
        // #pragma omp critical
        {
            eData[j] = c0;
            eData[j + 1] = c1;
        }
    }
}

int main()
{
    std::cout << "Nr=" << Nr << std::endl;
    //=============客户端offline阶段================
    // 定义初始向量
    vector<long> IV(BlockByte);
    for (unsigned i = 0; i < BlockByte; i++)
    {
        IV[i] = i + 1;
    }
    // 生成随机对称密钥
    GF2X rnd;
    int Bytebitsdiv8 = ceil(Bytebits / 8);
    vector<uint8_t> SymKey0(Bytebitsdiv8 * BlockByte);
    random(rnd, 8 * SymKey0.size());
    BytesFromGF2X(SymKey0.data(), rnd, SymKey0.size());
    vector<long> SymKey(BlockByte);
    for (unsigned i = 0; i < BlockByte; i++)
    {
        SymKey[i] = 0;
        for (unsigned j = 0; j < Bytebitsdiv8; j++)
        {
            SymKey[i] += (SymKey0[Bytebitsdiv8 * i + j] << (8 * j));
        }
        SymKey[i] %= PlainMod;
    }

    std::cout << "SymKey generated." << std::endl;
    //========
    // Generating symmetric key and key stream
    auto start_keyStream = std::chrono::steady_clock::now();

    std::vector<long> NonceSet(PlainBlock);
    std::vector<long> Xset(PlainByte * (Nr + 1));
    std::vector<long> RoundKeySet(PlainByte * (Nr + 1));
    std::vector<long> KeyStream(PlainByte);

    long total_tasks = counter_end - counter_begin + 1;
    std::atomic<long> completed_tasks(0);
    // 定义进度打印的粒度，例如每完成 1% 的任务打印一次
    long progress_step = total_tasks / 100;
    if (progress_step == 0)
        progress_step = 1; // 防止除零
    RandomBit<BlockSize> randomBit(Nr);
#pragma omp parallel for firstprivate(randomBit)
    for (long counter = counter_begin; counter <= counter_end; counter++)
    {
        long nonce = generate_secure_random_int(NonceSize);
        randomBit.generate_Instance_all_new(nonce, counter);
        auto &RanVecs = randomBit.roundconstants;
        // long nonce = counter;
        // std::vector<std::bitset<544>> RanVecs(Nr + 1);

        NonceSet[counter - counter_begin] = nonce;
        // 使用 std::array 代替 vector 并固定大小
        std::vector<long> state(BlockByte); // 初始化 state

        // 逐轮进行加密
        for (unsigned r = 0; r <= Nr; r++)
        {
            std::vector<long> X(BlockByte);
            std::vector<long> RoundKey(BlockByte);
            uint64_t temp;
            // 计算 X 和 RoundKey
            for (unsigned i = 0; i < BlockByte; ++i)
            {
                bool bit_array[Bytebits];
                for (unsigned j = 0; j < Bytebits; ++j)
                {
                    bit_array[j] = RanVecs[r][i * Bytebits + j];
                }
                BinStrToHex(bit_array, temp, Bytebits);
                // 强制转换为 long 类型
                X[i] = static_cast<long>(temp % PlainMod);
                if (Rkflag)
                {
                    RoundKey[i] = (SymKey[i] * X[i]) % PlainMod;
                }
                else
                {
                    RoundKey[i] = (SymKey[i] + X[i]) % PlainMod;
                }
            }

            // 将 X 和 RoundKey 复制到 Xset 和 RoundKeySet
            memcpy(&Xset[PlainByte * r + BlockByte * (counter - counter_begin)], X.data(), BlockByte * sizeof(long));
            memcpy(&RoundKeySet[PlainByte * r + BlockByte * (counter - counter_begin)], RoundKey.data(), BlockByte * sizeof(long));
            if (r == 0)
            { // 初始轮
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (RoundKey[i] + IV[i]) % PlainMod;
                }
            }
            else if (r < Nr)
            {                     // 常规轮
                yusP.MC64(state); // 列混淆
                yusP.MR64(state); // 行移位
                yusP.Sbox(state); // S盒
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
            }
            else
            {                     // 最后一轮
                yusP.MC64(state); // 列混淆
                yusP.MR64(state); // 行移位
                yusP.Sbox(state); // S盒
                yusP.MC64(state); // 列混淆
                yusP.MR64(state); // 行移位
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
                memcpy(&KeyStream[(counter - counter_begin) * BlockByte], state.data(), BlockByte * sizeof(long));
            }
        }

        // 更新已完成的任务数
        long local_completed = ++completed_tasks;

        // 定期打印进度
        if (local_completed % progress_step == 0 || local_completed == total_tasks)
        {
#pragma omp critical
            {
                std::cout << "Progress: " << (local_completed * 100) / total_tasks << "% completed.\r" << std::flush;
            }
        }
    }
    std::cout << std::endl; // 完成后换行

    auto end_keyStream = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_keyStream = end_keyStream - start_keyStream;
    std::cout << "KeyStream Generation time: " << elapsed_seconds_keyStream.count() << "s\n";

    // Generating Public Key and encrypting the symmetric key

    int idx = 1;
    long p = mValues[idx][0];
    long m = mValues[idx][1]; // m=65536,phi(m)=32768=2^15
    long r = 1;
    long bits = mValues[idx][2];
    long c = mValues[idx][3];
    long d = 1; // slots = phi(m)/d = phi(m) = 32768 = PlainBlock
    long k = 128;
    long s = 1;

    if (!m)
        m = FindM(k, bits, c, p, d, s, 0);
    auto start = std::chrono::steady_clock::now();
    shared_ptr<Context> context(ContextBuilder<BGV>()
                                    .m(m)
                                    .p(p)
                                    .r(r)
                                    .bits(bits)
                                    .c(c)
                                    .buildPtr());
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_context = end - start;
    std::cout << "Context build time: " << elapsed_seconds_context.count() << "s\n";

    auto start_PubKey = std::chrono::steady_clock::now();
    SecKey secretKey(*context);
    secretKey.GenSecKey();
    unique_ptr<PubKey> publicKey = std::make_unique<helib::PubKey>(secretKey);
    helib::EncryptedArray ea(context->getEA());
    long nslots = ea.size();
    if (nslots != PlainBlock)
    {
        std::cerr << "nslots != PlainBlock" << std::endl;
        return false;
    }
    std::cout << "p=" << p << std::endl;
    std::cout << "m=" << m << std::endl;
    std::cout << "nslots=" << nslots << std::endl;
    std::cout << "bits=" << bits << std::endl;
    std::cout << "c=" << c << std::endl;
    auto end_PubKey = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_PubKey = end_PubKey - start_PubKey;
    std::cout << "PublicKey generation time: " << elapsed_seconds_PubKey.count() << "s\n";
    return 0;

    auto start_keyEncryption = std::chrono::steady_clock::now();
    vector<Ctxt> encryptedSymKey;
    encryptSymKey(encryptedSymKey, SymKey, publicKey, ea);
    auto end_keyEncryption = std::chrono::steady_clock::now();
    double keyEncryption = std::chrono::duration<double>(end_keyEncryption - start_keyEncryption).count();
    std::cout << "SymKey FHE time: " << keyEncryption << "s\n";

    // 解密验证
    if (!verify_encryptSymKey(encryptedSymKey, SymKey, secretKey, ea))
    {
        return 0;
    }
    std::cout << "Symmetric key encryption succeeded!" << std::endl;

    // 离线客户端时间=KeyStream Generation time+PublicKey generation and SymKey FHE time
    double total_time_off = elapsed_seconds_keyStream.count() + elapsed_seconds_PubKey.count() + keyEncryption;
    std::cout << "Encryption offline total time: " << total_time_off << "s\n";
    //=============服务端offline阶段================
    // 计算 encryptedRoundKeySet
    auto start_RoundKeySet_FHE = std::chrono::steady_clock::now();
    vector<Ctxt> encryptedRoundKeySet;
    Ctxt tmpCtxt(*publicKey);
    long eRk_len = BlockByte * (Nr + 1);
    encryptedRoundKeySet.resize(eRk_len, tmpCtxt);

    for (int i = 0; i < eRk_len; i++)
    {
        encryptedRoundKeySet[i] = encryptedSymKey[i % BlockByte];
    }

    vector<ZZX> encodedXset;
    encodeTo64Ctxt(encodedXset, Xset, ea); // encode as HE plaintext
    if (Rkflag)
    {
        for (int i = 0; i < eRk_len; i++)
        {
            encryptedRoundKeySet[i].multByConstant(encodedXset[i]);
        }
    }
    else
    {
        for (int i = 0; i < eRk_len; i++)
        {
            encryptedRoundKeySet[i].addConstant(encodedXset[i]);
        }
    }
    auto end_RoundKeySet_FHE = std::chrono::steady_clock::now();
    // 使用 verifyDecryption 函数解密并验证 RoundKeySet
    if (!verifyDecryption64(encryptedRoundKeySet, RoundKeySet, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for RoundKeySet." << std::endl;
        return 0;
    }
    std::cout << "Decryption verification succeeded for RoundKeySet." << std::endl;

    std::chrono::duration<double> elapsed_seconds_RoundKeySet_FHE = end_RoundKeySet_FHE - start_RoundKeySet_FHE;
    std::cout << "RoundKeySet FHE succeeded! Time: " << elapsed_seconds_RoundKeySet_FHE.count() << "s\n";
    // // 测试赋值和加法时间
    // Ctxt tt = encryptedRoundKeySet[0];
    // auto start_assign = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < 10000; i++)
    // {
    //     tt = encryptedRoundKeySet[i%eRk_len];
    // }
    // auto end_assign = std::chrono::high_resolution_clock::now();
    // std::cout << "assign time: " << std::chrono::duration<double>(end_assign - start_assign).count() << "s\n";
    // // 测试addConstant时间
    // tt = encryptedRoundKeySet[0];
    // auto start_addConstant = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < 10000; i++)
    // {
    //     tt+=encryptedRoundKeySet[i%eRk_len];
    // }
    // auto end_addConstant = std::chrono::high_resolution_clock::now();
    // std::cout << "addConstant time: " << std::chrono::duration<double>(end_addConstant - start_addConstant).count() << "s\n";
    // // 测试addConstant时间
    // vector<Ctxt> tt(eRk_len, tmpCtxt);
    // for (int i = 0; i < eRk_len; i++)
    // {
    //     tt[i] = encryptedRoundKeySet[i];
    // }
    // auto start_addConstant = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < eRk_len; i++)
    // {
    //     tt[i].addConstant(encodedXset[i]);
    // }
    // auto end_addConstant = std::chrono::high_resolution_clock::now();
    // std::cout << "addConstant time: " << std::chrono::duration<double>(end_addConstant - start_addConstant).count() << "s\n";
    // //测试multByConstant时间
    // vector<Ctxt> tt2(eRk_len, tmpCtxt);
    // for (int i = 0; i < eRk_len; i++)
    // {
    //     tt2[i] = encryptedRoundKeySet[i];
    // }
    // auto start_multByConstant = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < eRk_len; i++)
    // {
    //     tt2[i].multByConstant(encodedXset[i]);
    // }
    // auto end_multByConstant = std::chrono::high_resolution_clock::now();
    // std::cout << "multByConstant time: " << std::chrono::duration<double>(end_multByConstant - start_multByConstant).count() << "s\n";
    // // 测试multiplyBy时间
    // vector<Ctxt> tt3(eRk_len, tmpCtxt);
    // for (int i = 0; i < eRk_len; i++)
    // {
    //     tt3[i] = encryptedRoundKeySet[i];
    // }
    // auto start_multiplyBy = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < eRk_len; i++)
    // {
    //     tt3[i]*=tt3[i];
    // }
    // auto end_multiplyBy = std::chrono::high_resolution_clock::now();
    // std::cout << "multiplyBy time: " << std::chrono::duration<double>(end_multiplyBy - start_multiplyBy).count() << "s\n";
    // // 测试square时间
    // vector<Ctxt> tt4(eRk_len, tmpCtxt);
    // for (int i = 0; i < eRk_len; i++)
    // {
    //     tt4[i] = encryptedRoundKeySet[i];
    // }
    // auto start_square = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < eRk_len; i++)
    // {
    //     tt4[i].square();
    // }
    // auto end_square = std::chrono::high_resolution_clock::now();
    // std::cout << "square time: " << std::chrono::duration<double>(end_square - start_square).count() << "s\n";
    // 测试+=时间
    // vector<Ctxt> tt5(encryptedRoundKeySet);
    // auto start_plus = std::chrono::high_resolution_clock::now();
    // for (int i = 0; i < eRk_len; i++)
    // {
    //     tt5[i]+=tt5[i];
    // }
    // auto end_plus = std::chrono::high_resolution_clock::now();
    // std::cout << "+= time: " << std::chrono::duration<double>(end_plus - start_plus).count() << "s\n";

    // 生成 encryptedKeyStream
    // 定义roundkey_time、sbox_time、linear_layer_time
    double sbox_time = 0, linear_layer_time = 0, roundkey_time = 0;

    vector<Ctxt> encryptedKeyStream;
    encryptedKeyStream.resize(BlockByte, tmpCtxt);
    std::copy(encryptedRoundKeySet.begin(), encryptedRoundKeySet.begin() + BlockByte, encryptedKeyStream.begin());

    vector<long> expandedIV(BlockByte * PlainBlock);
    for (long j = 0; j < PlainBlock; j++)
    {
        memcpy(&expandedIV[BlockByte * j], IV.data(), BlockByte * sizeof(long));
    }
    // 对expanded进行simd编码，这样会返回nRoundKeys个多项式数组即encoded，nRoundKeys=encoded.length()
    vector<ZZX> encoded_expandedIV;
    encodeTo64Ctxt(encoded_expandedIV, expandedIV, ea); // encode as HE plaintext

    std::cout << "whiteround start" << std::endl;
    auto start_roundkey = std::chrono::high_resolution_clock::now();
    for (long i = 0; i < BlockByte; i++)
    { // encrypt the encoded key
        encryptedKeyStream[i].addConstant(encoded_expandedIV[i]);
    }
    auto end_roundkey = std::chrono::high_resolution_clock::now();
    roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
    // 输出 roundkey_time
    std::cout << "whiteround time: " << roundkey_time << "s\n";
    // 明文密钥流
    vector<long> KeyStream2(PlainByte);
    // 对IV和RoundKeySet进行异或
    for (long i = 0; i < PlainByte; i++)
    {
        KeyStream2[i] = (expandedIV[i] + RoundKeySet[i]) % PlainMod;
    }
    // 使用 verifyDecryption 函数解密并验证 KeyStream
    if (!verifyDecryption64(encryptedKeyStream, KeyStream2, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for KeyStream." << std::endl;
        return 0;
    }
    std::cout << "Decryption verification succeeded for whiteround." << std::endl;

    auto start_sbox = std::chrono::high_resolution_clock::now();
    auto start_linear = std::chrono::high_resolution_clock::now();
    auto end_sbox = std::chrono::high_resolution_clock::now();
    auto end_linear = std::chrono::high_resolution_clock::now();

    for (long r = 1; r < Nr; r++)
    {
        std::cout << "Round " << r << " start" << std::endl;
        start_linear = std::chrono::high_resolution_clock::now();
        // #pragma omp parallel for
        // MC Layer
        HE_MC(encryptedKeyStream);
        //  MR Layer
        HE_MR(encryptedKeyStream);
        // // MC Layer+MR Layer
        // HE_MC_MR(encryptedKeyStream);
        end_linear = std::chrono::high_resolution_clock::now();
        linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
        for (int i = 0; i < PlainBlock; i++)
        {
            vector<long> tmp(BlockByte);
            for (int j = 0; j < BlockByte; j++)
            {
                tmp[j] = KeyStream2[i * BlockByte + j];
            }
            yusP.MC64(tmp);
            yusP.MR64(tmp);
            for (int j = 0; j < BlockByte; j++)
            {
                KeyStream2[i * BlockByte + j] = tmp[j];
            }
        }
        if (!verifyDecryption64(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
        start_sbox = std::chrono::high_resolution_clock::now();
        // S Layer
        HE_Sbox(encryptedKeyStream);
        end_sbox = std::chrono::high_resolution_clock::now();
        sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();

        yusP.Sbox(KeyStream2);
        if (!verifyDecryption64(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
        start_roundkey = std::chrono::high_resolution_clock::now();
        // #pragma omp parallel for
        for (long j = 0; j < BlockByte; j++)
        {
            encryptedKeyStream[j] += encryptedRoundKeySet[r * BlockByte + j];
        }
        end_roundkey = std::chrono::high_resolution_clock::now();
        roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
        //
        for (long i = 0; i < PlainByte; i++)
        {
            KeyStream2[i] = (KeyStream2[i] + RoundKeySet[r * PlainByte + i]) % PlainMod;
        }
        if (!verifyDecryption64(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;
    }

    // 最后一轮
    std::cout << "Round " << Nr << " start" << std::endl;
    start_linear = std::chrono::high_resolution_clock::now();
    // MC Layer
    HE_MC(encryptedKeyStream);
    // MR Layer
    HE_MR(encryptedKeyStream);
    // // MC Layer+MR Layer
    // HE_MC_MR(encryptedKeyStream);
    end_linear = std::chrono::high_resolution_clock::now();
    linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
    for (int i = 0; i < PlainBlock; i++)
    {
        vector<long> tmp(BlockByte);
        for (int j = 0; j < BlockByte; j++)
        {
            tmp[j] = KeyStream2[i * BlockByte + j];
        }
        yusP.MC64(tmp);
        yusP.MR64(tmp);
        for (int j = 0; j < BlockByte; j++)
        {
            KeyStream2[i * BlockByte + j] = tmp[j];
        }
    }
    if (!verifyDecryption64(encryptedKeyStream, KeyStream2, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
        return 0;
    }
    std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
    start_sbox = std::chrono::high_resolution_clock::now();
    // S Layer
    HE_Sbox(encryptedKeyStream);
    end_sbox = std::chrono::high_resolution_clock::now();
    sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
    yusP.Sbox(KeyStream2);
    if (!verifyDecryption64(encryptedKeyStream, KeyStream2, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
        return 0;
    }
    std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
    start_linear = std::chrono::high_resolution_clock::now();
    // MC Layer
    HE_MC(encryptedKeyStream);
    // MR Layer
    HE_MR(encryptedKeyStream);
    // // MC Layer+MR Layer
    // HE_MC_MR(encryptedKeyStream);
    end_linear = std::chrono::high_resolution_clock::now();
    linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
    for (int i = 0; i < PlainBlock; i++)
    {
        vector<long> tmp(BlockByte);
        for (int j = 0; j < BlockByte; j++)
        {
            tmp[j] = KeyStream2[i * BlockByte + j];
        }
        yusP.MC64(tmp);
        yusP.MR64(tmp);
        for (int j = 0; j < BlockByte; j++)
        {
            KeyStream2[i * BlockByte + j] = tmp[j];
        }
    }
    if (!verifyDecryption64(encryptedKeyStream, KeyStream2, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
        return 0;
    }
    std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
    // add
    start_roundkey = std::chrono::high_resolution_clock::now();
    for (long j = 0; j < BlockByte; j++)
    {
        encryptedKeyStream[j] += encryptedRoundKeySet[Nr * BlockByte + j];
    }
    end_roundkey = std::chrono::high_resolution_clock::now();
    roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
    for (long i = 0; i < PlainByte; i++)
    {
        KeyStream2[i] = (KeyStream2[i] + RoundKeySet[Nr * PlainByte + i]) % PlainMod;
    }
    if (!verifyDecryption64(encryptedKeyStream, KeyStream2, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
        return 0;
    }
    std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;

    // 输出 roundkey_time、sbox_time、linear_layer_time
    std::cout << "RoundKey time: " << roundkey_time << "s\n";
    std::cout << "Sbox time: " << sbox_time << "s\n";
    std::cout << "Linear Layer time: " << linear_layer_time << "s\n";
    // 计算总时间
    double total_time = roundkey_time + sbox_time + linear_layer_time + elapsed_seconds_RoundKeySet_FHE.count();
    std::cout << "Server offline total time: " << total_time << "s\n";
    return 0;
}