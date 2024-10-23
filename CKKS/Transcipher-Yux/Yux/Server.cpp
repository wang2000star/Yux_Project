#include "Server.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

bool Server_offline()
{    
// 生成初始向量
#if 1
    std::array<unsigned char, BlockByte> IV;
    for (unsigned i = 0; i < BlockByte; i++)
    {
        IV[i] = i + 1;
    }
#endif
// 读取Nonce.txt，生成Xset
#if 1
    std::vector<std::array<std::array<unsigned char, BlockByte>, Nr + 1>> Xset;
    // 创建RandomBit对象
    RandomBit<BlockSize> randomBit(Nr);
    // 为roundconstants创建别名
    auto &RanVecs = randomBit.roundconstants;

    // 读取Client_NonceSet.txt
    std::array<uint64_t, PlainBlock> NonceSet;

    // 从文件Client_NonceSet.txt读取NonceSet
    std::ifstream in("Client_NonceSet.txt");
    if (in.is_open())
    {
        for (unsigned i = 0; i < PlainBlock; i++)
        {
            in >> std::hex >> NonceSet[i];
        }
        in.close();
    }
    else
    {
        std::cerr << "Failed to open Client_NonceSet.txt for reading" << std::endl;
    }

    for (uint64_t counter = counter_begin; counter <= counter_end; counter++)
    {
        uint64_t nonce = NonceSet[counter - counter_begin];
        randomBit.generate_Instance_all_new(nonce, counter); // nonce, counter

        uint64_t temp;
        std::array<std::array<unsigned char, BlockByte>, Nr + 1> X;
        // 将RanVecs转换为X，比特转字节
        for (unsigned r = 0; r <= Nr; r++)
        {
            for (unsigned i = 0; i < BlockByte; ++i)
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
        Xset.push_back(X);
    }
#endif
// 读取Client_publickey，生成publicKey
#if 1
    // 设置
    long idx = 3; // 0
    long c = 9;
    bool packed = true;
    if (idx > 5)
        idx = 5;
    long p = mValues[idx][0];
    long m = mValues[idx][2];
    long bits = mValues[idx][4];

    Context context(ContextBuilder<BGV>()
                        .m(m)
                        .p(p)
                        .r(1)
                        .c(c)
                        .bits(bits)
                        .build());
    // initialize context
    PubKey publicKey(context);
    std::ifstream in2("Client_publickey", std::ios::binary);
    if (!in2.is_open())
    {
        std::cerr << "Failed to open " << "Client_publickey" << " for reading" << std::endl;
        throw std::runtime_error("Failed to open public key file");
    }
    publicKey.readFrom(in2, context);
    in2.close();
#endif
// 读取Client_encryptedSymKey.txt，生成encryptedSymKey
#if 1
    std::vector<Ctxt> encrypted_SymKey;
    std::ifstream in3("Client_encryptedSymKey.txt", std::ios::binary);
    if (in3.is_open())
    {
        while (in3.peek() != EOF)
        {
            Ctxt ctxt(publicKey);
            ctxt.read(in3);
            encrypted_SymKey.push_back(ctxt);
        }
        in3.close();
    }
    else
    {
        std::cerr << "Failed to open Client_encryptedSymKey.txt for reading" << std::endl;
    }
#endif
// 重点1：对Xset进行BGV公钥加密,得到encryptedXset
#if 1
    
    std::vector<std::vector<Ctxt>> encryptedXset;

#endif
// 重点2：让encryptedXset的每一个元素分别和encryptedSymKey进行BGV同态异或，得到encrypted_RoundKeySet
#if 1
    
    std::vector<Ctxt> encrypted_RoundKeySet;
#endif:
// 重点3：对IV进行BGV公钥加密,得到encryptedIV
#if 1
    std::vector<std::vector<Ctxt>> encryptedIV;
#endif
// 重点4：对roundConstants进行BGV公钥加密，得到encryptedRoundConstants
#if 1
    std::vector<std::vector<Ctxt>> encryptedRoundConstants;
#endif
/* 重点5：利用encrypted_RoundKeySet和encryptedIV，执行Yux加密/解密算法，得到encryptedKeyStream
   涉及的运算有比特异或、模多项式乘法、移位
   操作分为S盒、线性层、轮密钥加，共Nr轮，最后一轮没有线性层，有白化轮密钥加。
   核心点，使用simd打包批处理技术，提高单位时间的吞吐量
*/
#if 1
    std::vector<Ctxt> encryptedKeyStream;
#endif
// 将encryptedKeyStream写入到文件Server_encryptedKeyStream.txt
#if 1
    std::ofstream out5("Server_encryptedKeyStream.txt", std::ios::binary);
    if (out5.is_open())
    {
        for (const auto &ctxt : encryptedKeyStream)
        {
            ctxt.writeTo(out5);
        }
        out5.close();
    }
    else
    {
        std::cerr << "Failed to open Server_encryptedKeyStream.txt for writing" << std::endl;
    }
#endif
return true;
}

bool Server_online()
{
// 读取Client_publickey，生成publicKey
#if 1
    // 设置
    long idx = 3; // 0
    long c = 9;
    bool packed = true;
    if (idx > 5)
        idx = 5;
    long p = mValues[idx][0];
    long m = mValues[idx][2];
    long bits = mValues[idx][4];

    Context context(ContextBuilder<BGV>()
                        .m(m)
                        .p(p)
                        .r(1)
                        .c(c)
                        .bits(bits)
                        .build());
    // initialize context
    PubKey publicKey(context);
    std::ifstream in2("Client_publickey", std::ios::binary);
    if (!in2.is_open())
    {
        std::cerr << "Failed to open " << "Client_publickey" << " for reading" << std::endl;
        throw std::runtime_error("Failed to open public key file");
    }
    publicKey.readFrom(in2, context);
    in2.close();
#endif
// 读取Client_CipherStream.txt，生成cipherStream
#if 1
    std::vector<Ctxt> cipherStream;
    std::ifstream in4("Client_CipherStream.txt", std::ios::binary);
    if (in4.is_open())
    {
        while (in4.peek() != EOF)
        {
            Ctxt ctxt(publicKey);
            ctxt.read(in4);
            cipherStream.push_back(ctxt);
        }
        in4.close();
    }
    else
    {
        std::cerr << "Failed to open Client_CipherStream.txt for reading" << std::endl;
    }
#endif
// 重点6：对cipherStream进行BGV公钥加密,得到encryptedCipherStream
#if 1 
    std::vector<Ctxt> encryptedCipherStream;
#endif
// 重点7：对encryptedCipherStream和encryptedKeyStream进行BGV同态异或，得到encryptedPlainStream
#if 1
    std::vector<Ctxt> encryptedPlainStream;
#endif
// 将encryptedPlainStream写入到文件Server_encryptedPlainStream.txt
#if 1
    std::ofstream out6("Server_encryptedPlainStream.txt", std::ios::binary);
    if (out6.is_open())
    {
        for (const auto &ctxt : encryptedPlainStream)
        {
            ctxt.writeTo(out6);
        }
        out6.close();
    }
    else
    {
        std::cerr << "Failed to open Server_encryptedPlainStream.txt for writing" << std::endl;
    }
#endif
return true;
}