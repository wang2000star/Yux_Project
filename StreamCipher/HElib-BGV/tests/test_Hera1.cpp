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
#include "Hera.hpp"
#include "tool.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

/************
 * PlainMod = 257, Bytebits = 9, BlockSize = 32*9 = 288, PlainBlock = (Plainmod-1)/2 = 128
 * PlainByte = 65537, Bytebits = 17, BlockSize = 32*17 = 544, PlainBlock = (Plainmod-1)/2 = 32768
 */
// constexpr long PlainMod = 257;    //2^8+1
// constexpr unsigned Bytebits = 9;       // 字节比特长度=ceil(log2(PlainMod))
constexpr long PlainMod = 65537;  // 2^16+1
constexpr unsigned Bytebits = 17; // 字节比特长度=ceil(log2(PlainMod))
constexpr long BlockByte = 16;    // 分组字节长度

constexpr unsigned BlockSize = Bytebits * BlockByte; // 分组比特长度=BlockByte*Bytebits
constexpr unsigned PlainBlock = (PlainMod - 1) / 4;  // 明文分组数

static const long PlainByte = BlockByte * PlainBlock; // 明文字节长度
// PlainByte = nslots
static const long PlainSize = BlockSize * PlainBlock;           // 明文比特长度
static const unsigned NonceSize = 32;                           // Nonce比特长度
static const unsigned Nr = 4;                                   // 轮数
static const long counter_begin = 0;                            // 计数器起始值
static const long counter_end = PlainBlock + counter_begin - 1; // 计数器结束值
static bool Rkflag = 1;                                         // true/1表示乘法，false/0表示加法

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
    {65537, 131072, 1320, 6},
    //{  257,    256,     500,  17},
    {65537, 32768, 300, 2},
};
// p^d = 1 mod m,d=1,slots=\phi(m)/d=\phi(m);m=65536=2^16,\phi(m)=2^15=32768

Hera hera(PlainMod);

void encodeTo16Ctxt(vector<ZZX> &encData, const vector<long> &data, const EncryptedArray &ea)
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
// encodeTo16Ctxt对应的解码
void decodeTo16Ctxt(vector<long> &data, const vector<vector<long>> &encData,
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
bool verifyDecryption16(const std::vector<Ctxt> &encryptedVec, const vector<long> &originalVec, const SecKey &secretKey,
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
    decodeTo16Ctxt(decryptedVec, decryptedPolys, ea);
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
    encryptedSymKey.resize(32, Ctxt(*pk));
    for (long i = 0; i < 32; i++)
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
    vector<Ctxt> temp(BlockByte, eData[0]);
    vector<int> index = {0, 1, 2, 3};
    Ctxt T2(eData[0]);
    Ctxt T3(eData[0]);
    for (int i = 0; i < 4; i++)
    {
        int s = 4 * i;
        T2 = eData[index[0] + s];
        T2 += eData[index[0] + s];
        T3 = eData[index[1] + s];
        T3 += eData[index[1] + s];
        T3 += eData[index[1] + s];
        temp[index[0] + s] = T2;
        temp[index[0] + s] += T3;
        temp[index[0] + s] += eData[index[2] + s];
        temp[index[0] + s] += eData[index[3] + s];

        T2 = eData[index[1] + s];
        T2 += eData[index[1] + s];
        T3 = eData[index[2] + s];
        T3 += eData[index[2] + s];
        T3 += eData[index[2] + s];
        temp[index[1] + s] = T2;
        temp[index[1] + s] += T3;
        temp[index[1] + s] += eData[index[3] + s];
        temp[index[1] + s] += eData[index[0] + s];

        T2 = eData[index[2] + s];
        T2 += eData[index[2] + s];
        T3 = eData[index[3] + s];
        T3 += eData[index[3] + s];
        T3 += eData[index[3] + s];
        temp[index[2] + s] = T2;
        temp[index[2] + s] += T3;
        temp[index[2] + s] += eData[index[0] + s];
        temp[index[2] + s] += eData[index[1] + s];

        T2 = eData[index[3] + s];
        T2 += eData[index[3] + s];
        T3 = eData[index[0] + s];
        T3 += eData[index[0] + s];
        T3 += eData[index[0] + s];
        temp[index[3] + s] = T2;
        temp[index[3] + s] += T3;
        temp[index[3] + s] += eData[index[1] + s];
        temp[index[3] + s] += eData[index[2] + s];
    }
    std::copy(temp.begin(), temp.end(), eData.begin());
}
void HE_MR(vector<Ctxt> &eData)
{
    vector<Ctxt> temp(BlockByte, eData[0]);
    vector<int> index = {0, 4, 8, 12};
    Ctxt T2(eData[0]);
    Ctxt T3(eData[0]);
    for (int i = 0; i < 4; i++)
    {
        int s = i;
        T2 = eData[index[0] + s];
        T2 += eData[index[0] + s];
        T3 = eData[index[1] + s];
        T3 += eData[index[1] + s];
        T3 += eData[index[1] + s];
        temp[index[0] + s] = T2;
        temp[index[0] + s] += T3;
        temp[index[0] + s] += eData[index[2] + s];
        temp[index[0] + s] += eData[index[3] + s];

        T2 = eData[index[1] + s];
        T2 += eData[index[1] + s];
        T3 = eData[index[2] + s];
        T3 += eData[index[2] + s];
        T3 += eData[index[2] + s];
        temp[index[1] + s] = T2;
        temp[index[1] + s] += T3;
        temp[index[1] + s] += eData[index[3] + s];
        temp[index[1] + s] += eData[index[0] + s];

        T2 = eData[index[2] + s];
        T2 += eData[index[2] + s];
        T3 = eData[index[3] + s];
        T3 += eData[index[3] + s];
        T3 += eData[index[3] + s];
        temp[index[2] + s] = T2;
        temp[index[2] + s] += T3;
        temp[index[2] + s] += eData[index[0] + s];
        temp[index[2] + s] += eData[index[1] + s];

        T2 = eData[index[3] + s];
        T2 += eData[index[3] + s];
        T3 = eData[index[0] + s];
        T3 += eData[index[0] + s];
        T3 += eData[index[0] + s];
        temp[index[3] + s] = T2;
        temp[index[3] + s] += T3;
        temp[index[3] + s] += eData[index[1] + s];
        temp[index[3] + s] += eData[index[2] + s];
    }
    std::copy(temp.begin(), temp.end(), eData.begin());
}
// Compute the constants for Sbox
void HE_Sbox(vector<Ctxt> &eData)
{
    // #pragma omp parallel for
    for (long j = 0; j < BlockByte; j++)
    {
        Ctxt temp(eData[j]);
        temp.multiplyBy(eData[j]);
        temp.multiplyBy(eData[j]);
        eData[j] = temp;
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
            { // 常规轮

                hera.MR(state);   // 行移位
                hera.MC(state);   // 列混淆
                hera.Sbox(state); // S盒
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
            }
            else
            {                     // 最后一轮
                hera.MR(state);   // 行移位
                hera.MC(state);   // 列混淆
                hera.Sbox(state); // S盒
                hera.MR(state);   // 再次行移位
                hera.MC(state);   // 再次列混淆
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
    auto start_PubKey = std::chrono::steady_clock::now();

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

    shared_ptr<Context> context(ContextBuilder<BGV>()
                                    .m(m)
                                    .p(p)
                                    .r(r)
                                    .bits(bits)
                                    .c(c)
                                    .buildPtr());

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
    encodeTo16Ctxt(encodedXset, Xset, ea); // encode as HE plaintext
    for (int i = 0; i < eRk_len; i++)
    {
        if (Rkflag)
        {
            encryptedRoundKeySet[i].multByConstant(encodedXset[i]);
        }
        else
        {
            encryptedRoundKeySet[i].addConstant(encodedXset[i]);
        }
    }
    auto end_RoundKeySet_FHE = std::chrono::steady_clock::now();
    // 使用 verifyDecryption 函数解密并验证 RoundKeySet
    if (!verifyDecryption16(encryptedRoundKeySet, RoundKeySet, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for RoundKeySet." << std::endl;
        return 0;
    }
    std::cout << "Decryption verification succeeded for RoundKeySet." << std::endl;

    std::chrono::duration<double> elapsed_seconds_RoundKeySet_FHE = end_RoundKeySet_FHE - start_RoundKeySet_FHE;
    std::cout << "RoundKeySet FHE succeeded! Time: " << elapsed_seconds_RoundKeySet_FHE.count() << "s\n";

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
    encodeTo16Ctxt(encoded_expandedIV, expandedIV, ea); // encode as HE plaintext

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
    if (!verifyDecryption16(encryptedKeyStream, KeyStream2, secretKey, ea))
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
        //  MR Layer
        HE_MR(encryptedKeyStream);
        // MC Layer
        HE_MC(encryptedKeyStream);
        end_linear = std::chrono::high_resolution_clock::now();
        linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
        for (int i = 0; i < PlainBlock; i++)
        {
            vector<long> tmp(BlockByte);
            for (int j = 0; j < BlockByte; j++)
            {
                tmp[j] = KeyStream2[i * BlockByte + j];
            }
            hera.MR(tmp);
            hera.MC(tmp);
            for (int j = 0; j < BlockByte; j++)
            {
                KeyStream2[i * BlockByte + j] = tmp[j];
            }
        }
        if (!verifyDecryption16(encryptedKeyStream, KeyStream2, secretKey, ea))
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

        hera.Sbox(KeyStream2);
        if (!verifyDecryption16(encryptedKeyStream, KeyStream2, secretKey, ea))
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
        if (!verifyDecryption16(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
            return 0;
        }
        std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;
    }

    // 最后一轮
    std::cout << "Round " << Nr << " start" << std::endl;
    start_linear = std::chrono::high_resolution_clock::now();
    // MR Layer
    HE_MR(encryptedKeyStream);
    // MC Layer
    HE_MC(encryptedKeyStream);
    end_linear = std::chrono::high_resolution_clock::now();
    linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
    for (int i = 0; i < PlainBlock; i++)
    {
        vector<long> tmp(BlockByte);
        for (int j = 0; j < BlockByte; j++)
        {
            tmp[j] = KeyStream2[i * BlockByte + j];
        }
        hera.MR(tmp);
        hera.MC(tmp);
        for (int j = 0; j < BlockByte; j++)
        {
            KeyStream2[i * BlockByte + j] = tmp[j];
        }
    }
    if (!verifyDecryption16(encryptedKeyStream, KeyStream2, secretKey, ea))
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
    hera.Sbox(KeyStream2);
    if (!verifyDecryption16(encryptedKeyStream, KeyStream2, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
        return 0;
    }
    std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
    start_linear = std::chrono::high_resolution_clock::now();
    // MR Layer
    HE_MR(encryptedKeyStream);
    // MC Layer
    HE_MC(encryptedKeyStream);
    end_linear = std::chrono::high_resolution_clock::now();
    linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
    for (int i = 0; i < PlainBlock; i++)
    {
        vector<long> tmp(BlockByte);
        for (int j = 0; j < BlockByte; j++)
        {
            tmp[j] = KeyStream2[i * BlockByte + j];
        }
        hera.MR(tmp);
        hera.MC(tmp);
        for (int j = 0; j < BlockByte; j++)
        {
            KeyStream2[i * BlockByte + j] = tmp[j];
        }
    }
    if (!verifyDecryption16(encryptedKeyStream, KeyStream2, secretKey, ea))
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
    if (!verifyDecryption16(encryptedKeyStream, KeyStream2, secretKey, ea))
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