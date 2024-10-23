#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <climits>
#include "random_bit.hpp"
#include "tool.hpp"
#include "Yux2_8.hpp"

int main()
{

    // 初始化常量
    const unsigned blocksize = 16; // 分组长度
    const unsigned Nr = 6;         // 轮数
    // 初始密钥
    unsigned char key[16] = {0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};
    // 初始IV
    unsigned char iv[16];
    for (unsigned i = 0; i < 16; i++)
    {
        iv[i] = i + 1;
    }
    unsigned plain_block = 24;//明文分组数
    unsigned char PlainStream[16 * plain_block];
    for (unsigned i = 0; i < 16 * plain_block; i++)
    {
        PlainStream[i] = i + 1;
    }
    // 创建RandomBit对象
    RandomBit<8 * blocksize> randomBit(Nr);
    auto &RanVecs = randomBit.roundconstants;
    // 生成Nonce
    unsigned m = 64; // 生成m比特强度的随机整数, m <= 64
    uint64_t counter_begin = 0;//
    uint64_t counter_end = plain_block + counter_begin-1;
    uint64_t counter;
    
    uint64_t nonce_set[plain_block];
    std::vector<std::array<std::array<unsigned char, 16>, Nr + 1>> X_set;
    std::vector<std::array<std::array<unsigned char, 16>, Nr + 1>> RoundKey_set;
    
    unsigned char KeyStream[16 * plain_block];
    unsigned char CipherStream[16 * plain_block];

    for (counter = counter_begin; counter <= counter_end; counter++)
    {
        uint64_t nonce = generate_secure_random_int(m);
        nonce_set[counter - counter_begin] = nonce;
        randomBit.generate_Instance_all_new(nonce, counter); // nonce, counter

        // 为LinMatrices创建别名
        // auto &RanMats = randomBit.LinMatrices;
        // 为roundconstants创建别名

        uint64_t temp;
        std::array<std::array<unsigned char, 16>, Nr + 1> X;
        // 将RanVecs转换为X，比特转字节
        for (unsigned r = 0; r <= Nr; r++)
        {
            for (unsigned i = 0; i < 16; ++i)
            {
                bool bit_array[8];
                for (unsigned j = 0; j < 8; ++j)
                {
                    bit_array[j] = RanVecs[r][i * 8 + j];
                }
                BinStrToHex(bit_array, temp, 8);
                X[r][i] = static_cast<unsigned char>(temp);
            }
        }
        // 把X存入X_set
        X_set.push_back(X);

        // 生成轮密钥
        std::array<std::array<unsigned char, 16>, Nr + 1> RoundKey;
        for (unsigned r = 0; r <= Nr; r++)
        {
            for (unsigned i = 0; i < 16; i++)
            {
                RoundKey[r][i] = mul(key[i], X[r][i]);
            }
        }
        // 把RoundKey存入RoundKey_set
        RoundKey_set.push_back(RoundKey);

        // the first round
        std::array<unsigned char, 16> state;

        for (unsigned i = 0; i < 16; i++)
        {
            state[i] = RoundKey[0][i] ^ iv[i];
        }
        // Nr-1 轮常规 Yux 轮
        for (unsigned r = 1; r < Nr; r++)
        {
            // S Layer -- 4 sbox
            for (unsigned i = 0; i < 4; i++)
            {
                encSboxFi(state.data(), i * 4);
                encSboxFi(state.data(), i * 4);
                encSboxFi(state.data(), i * 4);
                encSboxFi(state.data(), i * 4);
            }
            // linear Layer
            encLinearLayer(state.data());
            // addRoundKey
            for (unsigned i = 0; i < 16; i++)
            {
                state[i] ^= RoundKey[r][i];
            }
        }
        // 最后一轮
        for (unsigned i = 0; i < 4; i++)
        {
            encSboxFi(state.data(), i * 4);
            encSboxFi(state.data(), i * 4);
            encSboxFi(state.data(), i * 4);
            encSboxFi(state.data(), i * 4);
        }
        for (unsigned i = 0; i < 16; i++)
        {
            state[i] ^= RoundKey[Nr][i];
        }
        for (unsigned i = 0; i < 16; i++)
        {
            KeyStream[(counter - counter_begin) * 16 + i] = state[i];
        }
    }

    // 输出keystream
    printf("KeyStream :\n");
    for (unsigned i = 0; i < 16 * plain_block; i++)
    {
        printf("%02x ", KeyStream[i]);
    }
    std::cout << std::endl;
    // 用最终的keystream对明文进行加密
    for (unsigned i = 0; i < plain_block; i++)
    {
        for (unsigned j = 0; j < 16; j++)
        {
            CipherStream[i * 16 + j] = PlainStream[i * 16 + j] ^ KeyStream[i * 16 + j];
        }
    }
    // 输出密文
    printf("cipher:\n");
    for (unsigned i = 0; i < 16 * plain_block; i++)
    {
        printf("%02x ", CipherStream[i]);
    }
    std::cout << std::endl;
    // 用Keystream对密文进行解密
    unsigned char PlainStream2[16 * plain_block];
    for (unsigned i = 0; i < plain_block; i++)
    {
        for (unsigned j = 0; j < 16; j++)
        {
            PlainStream2[i * 16 + j] = CipherStream[i * 16 + j] ^ KeyStream[i*16 + j];
        }
    }
    // 检查解密是否正确
    for (unsigned i = 0; i < 16*plain_block; i++)
    {
        if (PlainStream[i] != PlainStream2[i])
        {
            std::cout << "Decryption failed!" << std::endl;
            return 1;
        }
    }
    std::cout << "Decryption verification succeeded!" << std::endl;
    //下面对key进行BGV加密，对X_set进行BGV加密，对RoundKey_set进行BGV加密
    // 然后生成BGV下Keystream,即FHE(KeyStream)
    // 最后用FHE(KeyStream)和密文cipher进行BGV异或，得到FHE(plain)
    return 0;
}