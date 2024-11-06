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

#include "params_Yus_p.hpp"
#include "Yus_p.hpp"
#include "tool.hpp"
using namespace std;
using namespace helib;
using namespace NTL;

void encodeTo32Ctxt(vector<ZZX> &encData, const vector<long> &data, const EncryptedArray &ea)
{
    long R = data.size() / PlainByte;
    std::cout << "R: " << R << std::endl;
    long nCtxt = BlockByte * R;
    std::cout << "nCtxt: " << nCtxt << std::endl;
    long data_length = data.size();
    encData.resize(nCtxt);
    std::cout << "data_length: " << data_length << std::endl;
    std::cout << "ea.size: " << ea.size() << std::endl;

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
// encodeTo32Ctxt对应的解码
void decodeTo32Ctxt(vector<long> &data, const vector<vector<long>> &encData,
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
bool verifyDecryption32(const std::vector<Ctxt> &encryptedVec, const vector<long> &originalVec, const SecKey &secretKey,
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
    decodeTo32Ctxt(decryptedVec, decryptedPolys, ea);
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
    vector<long> decryptedSymKey(32);
    for (long i = 0; i < 32; i++)
    { // encrypt the encoded key
        vector<long> slotsData;
        ea.decrypt(encryptedSymKey[i], secretKey, slotsData);
        decryptedSymKey[i] = slotsData[0];
    }
    bool isDecryptedSymKeyCorrect = true;
    for (long i = 0; i < 32; i++)
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
    vector<Ctxt> temp(32, eData[0]);
    vector<int> index = {0, 1, 2, 3, 4, 5, 6, 7};
    for (int i = 0; i < 4; i++)
    {
        int s = 8 * i;
        temp[index[0] + s] = eData[index[1] + s];
        temp[index[0] + s] += eData[index[2] + s];
        temp[index[0] + s] += eData[index[4] + s];
        temp[index[0] + s] += eData[index[5] + s];
        temp[index[0] + s] += eData[index[6] + s];

        temp[index[1] + s] = eData[index[0] + s];
        temp[index[1] + s] += eData[index[3] + s];
        temp[index[1] + s] += eData[index[4] + s];
        temp[index[1] + s] -= eData[index[7] + s];

        temp[index[2] + s] = eData[index[1] + s];
        temp[index[2] + s] += eData[index[2] + s];
        temp[index[2] + s] += eData[index[3] + s];
        temp[index[2] + s] += eData[index[4] + s];
        temp[index[2] + s] -= eData[index[6] + s];
        temp[index[2] + s] -= eData[index[7] + s];

        temp[index[3] + s] = eData[index[0] + s];
        temp[index[3] + s] += eData[index[1] + s];
        temp[index[3] + s] += eData[index[3] + s];
        temp[index[3] + s] += eData[index[5] + s];
        temp[index[3] + s] -= eData[index[6] + s];

        temp[index[4] + s] = eData[index[0] + s];
        temp[index[4] + s] += eData[index[1] + s];
        temp[index[4] + s] += eData[index[2] + s];
        temp[index[4] + s] += eData[index[4] + s];
        temp[index[4] + s] -= eData[index[6] + s];
        temp[index[4] + s] += eData[index[7] + s];

        temp[index[5] + s] = eData[index[0] + s];
        temp[index[5] + s] += eData[index[2] + s];
        temp[index[5] + s] += eData[index[3] + s];
        temp[index[5] + s] += eData[index[5] + s];
        temp[index[5] + s] += eData[index[6] + s];

        temp[index[6] + s] = eData[index[1] + s];
        temp[index[6] + s] -= eData[index[0] + s];
        temp[index[6] + s] -= eData[index[2] + s];
        temp[index[6] + s] += eData[index[3] + s];
        temp[index[6] + s] += eData[index[4] + s];
        temp[index[6] + s] += eData[index[6] + s];

        temp[index[7] + s] = eData[index[2] + s];
        temp[index[7] + s] -= eData[index[0] + s];
        temp[index[7] + s] += eData[index[3] + s];
        temp[index[7] + s] += eData[index[4] + s];
        temp[index[7] + s] += eData[index[5] + s];
        temp[index[7] + s] -= eData[index[6] + s];
        temp[index[7] + s] += eData[index[7] + s];
    }
    std::copy(temp.begin(), temp.end(), eData.begin());
}
void HE_MR(vector<Ctxt> &eData)
{
    vector<Ctxt> temp(32, eData[0]);
    vector<int> index = {0, 1, 8, 9, 16, 17, 24, 25};
    for (int i = 0; i < 4; i++)
    {
        int s = 2 * i;
        temp[index[0] + s] = eData[index[1] + s];
        temp[index[0] + s] += eData[index[2] + s];
        temp[index[0] + s] += eData[index[4] + s];
        temp[index[0] + s] += eData[index[5] + s];
        temp[index[0] + s] += eData[index[6] + s];

        temp[index[1] + s] = eData[index[0] + s];
        temp[index[1] + s] += eData[index[3] + s];
        temp[index[1] + s] += eData[index[4] + s];
        temp[index[1] + s] -= eData[index[7] + s];

        temp[index[2] + s] = eData[index[1] + s];
        temp[index[2] + s] += eData[index[2] + s];
        temp[index[2] + s] += eData[index[3] + s];
        temp[index[2] + s] += eData[index[4] + s];
        temp[index[2] + s] -= eData[index[6] + s];
        temp[index[2] + s] -= eData[index[7] + s];

        temp[index[3] + s] = eData[index[0] + s];
        temp[index[3] + s] += eData[index[1] + s];
        temp[index[3] + s] += eData[index[3] + s];
        temp[index[3] + s] += eData[index[5] + s];
        temp[index[3] + s] -= eData[index[6] + s];

        temp[index[4] + s] = eData[index[0] + s];
        temp[index[4] + s] += eData[index[1] + s];
        temp[index[4] + s] += eData[index[2] + s];
        temp[index[4] + s] += eData[index[4] + s];
        temp[index[4] + s] -= eData[index[6] + s];
        temp[index[4] + s] += eData[index[7] + s];

        temp[index[5] + s] = eData[index[0] + s];
        temp[index[5] + s] += eData[index[2] + s];
        temp[index[5] + s] += eData[index[3] + s];
        temp[index[5] + s] += eData[index[5] + s];
        temp[index[5] + s] += eData[index[6] + s];

        temp[index[6] + s] = eData[index[1] + s];
        temp[index[6] + s] -= eData[index[0] + s];
        temp[index[6] + s] -= eData[index[2] + s];
        temp[index[6] + s] += eData[index[3] + s];
        temp[index[6] + s] += eData[index[4] + s];
        temp[index[6] + s] += eData[index[6] + s];

        temp[index[7] + s] = eData[index[2] + s];
        temp[index[7] + s] -= eData[index[0] + s];
        temp[index[7] + s] += eData[index[3] + s];
        temp[index[7] + s] += eData[index[4] + s];
        temp[index[7] + s] += eData[index[5] + s];
        temp[index[7] + s] -= eData[index[6] + s];
        temp[index[7] + s] += eData[index[7] + s];
    }
    std::copy(temp.begin(), temp.end(), eData.begin());
}
// Compute the constants for Sbox
void HE_Sbox(vector<Ctxt> &eData)
{
    int begin;
    for (long j = 0; j < 16; j++)
    {
        begin = j * 2;
        Ctxt c0(eData[begin]);
        Ctxt c1(eData[begin + 1]);

        for (int i = 0; i < 2; i++)
        {
            Ctxt temp = c0;
            // c0.frobeniusAutomorph(1);
            c0.multiplyBy(c0);
            c0 += c1;
            c1 = temp;
        }

        eData[begin] = c0;
        eData[begin + 1] = c1;
    }
}

int main()
{
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
                RoundKey[i] = (SymKey[i] + X[i]) % PlainMod;
                // RoundKey[i] = (SymKey[i] * X[i]) % PlainMod;
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

                MR(state);   // 行移位
                MC(state);   // 列混淆
                Sbox(state); // S盒
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = (state[i] + RoundKey[i]) % PlainMod;
                }
            }
            else
            {                // 最后一轮
                MR(state);   // 行移位
                MC(state);   // 列混淆
                Sbox(state); // S盒
                MR(state);   // 再次行移位
                MC(state);   // 再次列混淆
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
    long m = mValues[idx][1]; // m=65536,phi(m)=32768=2^15
    long p = mValues[idx][0];
    long r = 1;
    long L = mValues[idx][2];
    long c = mValues[idx][3];
    long d = 1; // slots = phi(m)/d = phi(m) = 32768 = PlainBlock
    long k = 128;
    long s = 1;

    if (!m)
        m = FindM(k, L, c, p, d, s, 0);
    printf("Selected m = %ld\n", m);

    shared_ptr<Context> context(ContextBuilder<BGV>()
                                    .m(m)
                                    .p(p)
                                    .r(r)
                                    .bits(L)
                                    .c(c)
                                    .buildPtr());

    printf("Security");
    SecKey secretKey(*context);
    secretKey.GenSecKey();
    unique_ptr<PubKey> publicKey = std::make_unique<helib::PubKey>(secretKey);
    helib::EncryptedArray ea(context->getEA());
    long nslots = ea.size();
    printf("nslots = %ld\n", nslots);
    if (nslots != PlainBlock)
    {
        std::cerr << "nslots != PlainBlock" << std::endl;
        return false;
    }
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
    encodeTo32Ctxt(encodedXset, Xset, ea); // encode as HE plaintext
    for (int i = 0; i < eRk_len; i++)
    {
        encryptedRoundKeySet[i].addConstant(encodedXset[i]);
    }
    auto end_RoundKeySet_FHE = std::chrono::steady_clock::now();
    // 使用 verifyDecryption 函数解密并验证 RoundKeySet
    if (!verifyDecryption32(encryptedRoundKeySet, RoundKeySet, secretKey, ea))
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
    encodeTo32Ctxt(encoded_expandedIV, expandedIV, ea); // encode as HE plaintext

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
    if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
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
            MR(tmp);
            MC(tmp);
            for (int j = 0; j < BlockByte; j++)
            {
                KeyStream2[i * BlockByte + j] = tmp[j];
            }
        }
        if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
            return 0;
        }
        start_sbox = std::chrono::high_resolution_clock::now();
        // S Layer
        HE_Sbox(encryptedKeyStream);
        end_sbox = std::chrono::high_resolution_clock::now();
        sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();

        Sbox(KeyStream2);
        // 输出前32个KeyStream2
        for (int i = 0; i < 32; i++)
        {
            std::cout << KeyStream2[i] << " ";
        }
        std::cout << std::endl;
        if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
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
        if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
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
        MR(tmp);
        MC(tmp);
        for (int j = 0; j < BlockByte; j++)
        {
            KeyStream2[i * BlockByte + j] = tmp[j];
        }
    }
    if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
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
    Sbox(KeyStream2);
    if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
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
        MR(tmp);
        MC(tmp);
        for (int j = 0; j < BlockByte; j++)
        {
            KeyStream2[i * BlockByte + j] = tmp[j];
        }
    }
    if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
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
    Sbox(KeyStream2);
    if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
        return 0;
    }
    std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
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
    if (!verifyDecryption32(encryptedKeyStream, KeyStream2, secretKey, ea))
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
    double total_time = roundkey_time + sbox_time + linear_layer_time;
    std::cout << "Encryption total time: " << total_time << "s\n";
    return 0;
}