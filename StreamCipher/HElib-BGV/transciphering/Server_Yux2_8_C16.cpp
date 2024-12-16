#include "Server_Yux2_8_C16.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

bool readEncryptedSymKey(vector<Ctxt> &encryptedSymKey, const std::string &filename, const PubKey &publicKey)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open())
    {
        std::cerr << "Failed to open " << filename << " for reading" << std::endl;
        return false;
    }
    encryptedSymKey.clear();
    while (in.peek() != EOF)
    {
        Ctxt ctxt(publicKey);
        ctxt.read(in);
        encryptedSymKey.push_back(ctxt);
    }

    in.close();
    return true;
}
// 读取加密密钥流
bool readEncryptedKeyStream(vector<Ctxt> &encryptedKeyStream, const std::string &filename, const PubKey &publicKey)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open())
    {
        std::cerr << "Failed to open " << filename << " for reading" << std::endl;
        return false;
    }
    encryptedKeyStream.clear();
    while (in.peek() != EOF)
    {
        Ctxt ctxt(publicKey);
        ctxt.read(in);
        encryptedKeyStream.push_back(ctxt);
    }

    in.close();
    return true;
}

void encryptPlainStream(vector<Ctxt> &encryptedPlainStream, vector<uint8_t> &CipherStream, vector<Ctxt> &encryptedKeyStream, const PubKey &hePK,
                        const EncryptedArrayDerived<PA_GF2> &ea)
{
    vector<ZZX> encoded;
    encodeTo16Ctxt(encoded, CipherStream, ea); // encode as HE plaintext
                                               // 对encoded进行同态加密，得到eKey,
    Ctxt tmpCtxt(hePK);
    encryptedPlainStream.resize(encoded.size(), tmpCtxt); // allocate space

#pragma omp parallel for
    for (long i = 0; i < (long)encryptedPlainStream.size(); i++) // encrypt the encoded key
    {
        encryptedPlainStream[i] = encryptedKeyStream[i];
        encryptedPlainStream[i].addConstant(encoded[i]);
    }
}

void buildSboxConstant(Ctxt &encA,
                       const EncryptedArrayDerived<PA_GF2> &ea)
{
    // char --> GF2X --> ZZX -->Ctxt
    GF2X polyConstant;
    GF2XFromBytes(polyConstant, &roundConstant, 1);
    // cout << "----Round Constant: " << polyConstant << "  \n";
    vector<GF2X> slots(ea.size(), polyConstant);
    ZZX ZZXConstant;
    ea.encode(ZZXConstant, slots);
    encA.DummyEncrypt(ZZXConstant);
}

void buildallOne(Ctxt &allOne,
                 const EncryptedArrayDerived<PA_GF2> &ea)
{
    // char --> GF2X --> ZZX -->Ctxt
    unsigned char One = 0x1;
    GF2X polyConstant;
    GF2XFromBytes(polyConstant, &One, 1);
    // cout << "----Round Constant: " << polyConstant << "  \n";
    vector<GF2X> slots(ea.size(), polyConstant);
    ZZX ZZXConstant;
    ea.encode(ZZXConstant, slots);
    allOne.DummyEncrypt(ZZXConstant);
}
// Compute the constants for Sbox
void decSboxFunc2(vector<Ctxt> &eData, long begin, Ctxt &encA, Ctxt &allOne)
{
    Ctxt c0 = eData[begin];
    Ctxt c1 = eData[begin + 1];
    Ctxt c2 = eData[begin + 2];
    Ctxt c3 = eData[begin + 3];
    /*

    Ctxt temp(c1);
    temp.multiplyBy(c2);
    temp += c0;
    temp += c3;
    temp += encA;
    eData[begin] = c1;
    eData[begin + 1] = c2;
    eData[begin + 2] = c3;
    eData[begin + 3] = temp;*/

    Ctxt temp2 = c2;
    Ctxt temp3 = c3;
    Ctxt temp0 = allOne;

#pragma omp parallel for
    for (long i = 0; i < 2; i++)
    {
        switch (i)
        {
        case 0:
            temp2.multiplyBy(c1);
            temp2 += c0;
            temp2 += c3;
            temp2 += encA;
            break;
        case 1:
            temp3 += c1;
            temp0 += c2;
            temp3.multiplyBy(temp0);
            temp3 += c0;
            break;
        }
    }

    /*
        Ctxt temp1(c1);
        temp1.multiplyBy(c2);
        temp1 += c0;
        temp1 += c3;
        Ctxt temp2(temp1);
        temp2 += encA;

        Ctxt temp3(c2);
        temp3.multiplyBy(c3);
        temp3 += c1;
        temp3 += temp1;
    */
    eData[begin] = c2;
    eData[begin + 1] = c3;
    eData[begin + 2] = temp2;
    eData[begin + 3] = temp3;
}

void mydecSboxFunc(vector<Ctxt> &eData, long begin, Ctxt &encA, const EncryptedArrayDerived<PA_GF2> &ea)
{
    Ctxt c0 = eData[begin];
    Ctxt c1 = eData[begin + 1];
    Ctxt c2 = eData[begin + 2];
    Ctxt c3 = eData[begin + 3];
    Ctxt temp(c1);
#pragma omp parallel for
    for (int i = 0; i < 2; i++)
    {
        switch (i)
        {
        case 0:
            temp += c2;
            // 输出temp的长度
            //  移位操作实现 x^2
            temp.frobeniusAutomorph(1);
            break;
        case 1:
            c0 += c3;
            c0 += encA;
            break;
        }
    }
    temp += c0;
    eData[begin] = c1;
    eData[begin + 1] = c2;
    eData[begin + 2] = c3;
    eData[begin + 3] = temp;
}

void decSboxFunc(vector<Ctxt> &eData, long begin, Ctxt &encA, const EncryptedArrayDerived<PA_GF2> &ea)
{
    Ctxt c0 = eData[begin];
    Ctxt c1 = eData[begin + 1];
    Ctxt c2 = eData[begin + 2];
    Ctxt c3 = eData[begin + 3];
    Ctxt temp(c1);
    temp.multiplyBy(c2);
    temp += c0;
    temp += c3;
    temp += encA;
    eData[begin] = c1;
    eData[begin + 1] = c2;
    eData[begin + 2] = c3;
    eData[begin + 3] = temp;
}
// Server offline
namespace Server_Yux2_8_C16
{
    bool Server_offline()
    {
        // 生成初始向量 IV
        vector<uint8_t> IV(BlockByte);
        for (unsigned i = 0; i < BlockByte; i++)
        {
            IV[i] = i + 1;
        }

        // 读取 Client_NonceSet.txt，生成 Xset
        Vec<long> NonceSet(INIT_SIZE, PlainBlock);
        if (!readFromFile<long>(NonceSet.data(), "Client_NonceSet.txt", PlainBlock))
        {
            std::cerr << "Failed to open Client_NonceSet.txt for reading" << std::endl;
            return false;
        }
        std::cout << "NonceSet read succeeded!\n";
        vector<uint8_t> Xset(PlainByte * (Nr + 1));
        auto start_Xset = std::chrono::steady_clock::now();
        RandomBit<BlockSize> randomBit(Nr);
#pragma omp parallel for firstprivate(randomBit)
        for (long counter = counter_begin; counter <= counter_end; counter++)
        {
            long nonce = NonceSet[counter - counter_begin];
            randomBit.generate_Instance_all_new(nonce, counter);
            auto &RanVecs = randomBit.roundconstants;

            for (unsigned r = 0; r <= Nr; r++)
            {
                uint64_t temp = 0;
                std::vector<uint8_t> X(BlockByte);

                for (unsigned i = 0; i < BlockByte; ++i)
                {
                    bool bit_array[8];
                    for (unsigned j = 0; j < 8; ++j)
                    {
                        bit_array[j] = RanVecs[r][i * 8 + j];
                    }
                    BinStrToHex(bit_array, temp, 8);
                    X[i] = static_cast<uint8_t>(temp);
                }
                memcpy(&Xset[PlainByte * r + BlockByte * (counter - counter_begin)], X.data(), BlockByte);
            }
        }
        auto end_Xset = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds_Xset = end_Xset - start_Xset;
        std::cout << "Xset generation succeeded! Time: " << elapsed_seconds_Xset.count() << "s\n";

        // 读取公钥
        std::ifstream inContext("Client_context", std::ios::binary);
        if (!inContext.is_open())
        {
            std::cerr << "Failed to open Client_context for reading" << std::endl;
            throw std::runtime_error("Failed to open context file");
        }
        Context context = Context::readFrom(inContext);
        inContext.close();

        static const uint8_t YuxPolyBytes[] = {0x1B, 0x1};
        const GF2X YuxPoly = GF2XFromBytes(YuxPolyBytes, 2);
        EncryptedArrayDerived<PA_GF2> ea(context, YuxPoly, context.getAlMod());

        std::ifstream inPublicKey("Client_publickey", std::ios::binary);
        if (!inPublicKey.is_open())
        {
            std::cerr << "Failed to open Client_publickey for reading" << std::endl;
            throw std::runtime_error("Failed to open public key file");
        }
        PubKey publicKey = PubKey::readFrom(inPublicKey, context);
        inPublicKey.close();

        // 读取私钥
        std::ifstream inSecretKey("Client_secretkey", std::ios::binary);
        if (!inSecretKey.is_open())
        {
            std::cerr << "Failed to open Client_secretkey for reading" << std::endl;
            throw std::runtime_error("Failed to open secret key file");
        }
        SecKey secretKey = SecKey::readFrom(inSecretKey, context);
        inSecretKey.close();

        // 读取加密后的对称密钥
        vector<Ctxt> encryptedSymKey;
        if (!readEncryptedSymKey(encryptedSymKey, "Client_encryptedSymKey.bin", publicKey))
        {
            std::cerr << "Failed to open Client_encryptedSymKey.bin for reading" << std::endl;
            return false;
        }
        std::cout << "encryptedSymKey read succeeded!\n";

        // 读取Client_RoundKeySet.txt，生成 RoundKeySet
        Vec<uint8_t> RoundKeySet(INIT_SIZE, PlainByte * (Nr + 1));
        if (!readFromFile<uint8_t>(RoundKeySet.data(), "Client_RoundKeySet.txt", PlainByte * (Nr + 1)))
        {
            std::cerr << "Failed to open Client_RoundKeySet.txt for reading" << std::endl;
            return false;
        }
        std::cout << "RoundKeySet read succeeded!\n";

        vector<Ctxt> encryptedRoundKeySet;
        Ctxt tmpCtxt(publicKey);
        encryptedRoundKeySet.resize(BlockByte * (Nr + 1), tmpCtxt);
        for (int i = 0; i < BlockByte * (Nr + 1); i++)
        {
            encryptedRoundKeySet[i] = encryptedSymKey[i % BlockByte];
        }
        vector<ZZX> encodedXset;
        encodeTo16Ctxt(encodedXset, Xset, ea); // encode as HE plaintext
        long encodedXset_length = BlockByte * (Nr + 1);
        // 计算 encryptedRoundKeySet
        auto start_RoundKeySet_FHE = std::chrono::steady_clock::now();           
        // omp_set_num_threads(Nr+1);
        //#pragma omp parallel for
        for (long i = 0; i < encodedXset_length; i++) // encrypt the encoded key
        {
            encryptedRoundKeySet[i].multByConstant(encodedXset[i]);
        }
        auto end_RoundKeySet_FHE = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds_RoundKeySet_FHE = end_RoundKeySet_FHE - start_RoundKeySet_FHE;
        std::cout << "RoundKeySet FHE succeeded! Time: " << elapsed_seconds_RoundKeySet_FHE.count() << "s\n";
        // 使用 verifyDecryption 函数解密并验证 RoundKeySet
        // if (!verifyDecryption16(encryptedRoundKeySet, secretKey, ea, RoundKeySet))
        // {
        //     std::cerr << "Decryption verification failed for RoundKeySet." << std::endl;
        //     return false;
        // }
        // std::cout << "Decryption verification succeeded for RoundKeySet." << std::endl;

        // 加密 RoundConstants
        auto start_RoundConstants_FHE = std::chrono::steady_clock::now();
        Ctxt encA(ZeroCtxtLike, encryptedSymKey[0]);
        buildSboxConstant(encA, ea);
        // Ctxt allOne(ZeroCtxtLike, encryptedSymKey[0]);
        // buildallOne(allOne, ea);
        auto end_RoundConstants_FHE = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds_RoundConstants_FHE = end_RoundConstants_FHE - start_RoundConstants_FHE;
        std::cout << "RoundConstants FHE succeeded! Time: " << elapsed_seconds_RoundConstants_FHE.count() << "s\n";

        // 生成 encryptedKeyStream
        // 定义roundkey_time、sbox_time、linear_layer_time
        double sbox_time = 0, linear_layer_time = 0, roundkey_time = 0;

        vector<Ctxt> encryptedKeyStream;
        encryptedKeyStream.resize(BlockByte, tmpCtxt);
        std::copy(encryptedRoundKeySet.begin(), encryptedRoundKeySet.begin() + BlockByte, encryptedKeyStream.begin());
        vector<uint8_t> expandedIV(BlockByte * PlainBlock);
        for (long j = 0; j < PlainBlock; j++)
        {
            memcpy(&expandedIV[BlockByte * j], &IV[0], BlockByte);
        }
        // 对expanded进行simd编码，这样会返回nRoundKeys个多项式数组即encoded，nRoundKeys=encoded.length()
        vector<ZZX> encoded_expandedIV;
        encodeTo16Ctxt(encoded_expandedIV, expandedIV, ea); // encode as HE plaintext

        std::cout << "whiteround start" << std::endl;
        auto start_roundkey = std::chrono::high_resolution_clock::now();     
        //#pragma omp parallel for
        for (long i = 0; i < BlockByte; i++) // encrypt the encoded key
        {
            encryptedKeyStream[i].addConstant(encoded_expandedIV[i]);
        }
        auto end_roundkey = std::chrono::high_resolution_clock::now();
        roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
        // 输出 roundkey_time
        std::cout << "whiteround time: " << roundkey_time << "s\n";
        // 明文密钥流
        // Vec<uint8_t> KeyStream(INIT_SIZE, PlainByte);
        // // 对IV和RoundKeySet进行异或
        // for (long i = 0; i < PlainByte; i++)
        // {
        //     KeyStream[i] = IV[i % BlockByte] ^ RoundKeySet[i];
        // }

        // // 使用 verifyDecryption 函数解密并验证 KeyStream
        // if (!verifyDecryption16(encryptedKeyStream, secretKey, ea, KeyStream))
        // {
        //     // std::cerr << "Decryption verification failed for KeyStream." << std::endl;
        //     return false;
        // }
        // std::cout << "Decryption verification succeeded for whiteround." << std::endl;

        // 创建副本
        // vector<Ctxt> tempCtxt = encryptedKeyStream;

        // 测试
        Ctxt temp = encryptedKeyStream[0];
        auto start_mul = std::chrono::high_resolution_clock::now();
        for (int i=0;i<10;i++){
            temp.multiplyBy(temp);
        }
        auto end_mul = std::chrono::high_resolution_clock::now();
        temp = encryptedKeyStream[0];
        std::cout << "mul time: " << std::chrono::duration<double>(end_mul - start_mul).count() << "s\n";
        // auto start_multByConstant = std::chrono::steady_clock::now();
        // for (long i = 0; i < 100; i++) // encrypt the encoded key
        // {
        //     temp.multByConstant(encodedXset[0]);
        // }
        // auto end_multByConstant = std::chrono::steady_clock::now();
        // std::cout << "multByConstant time: " << std::chrono::duration<double>(end_multByConstant - start_multByConstant).count() << "s\n";
        auto start_fro = std::chrono::high_resolution_clock::now();
        for (int i=0;i<10;i++){
            temp.frobeniusAutomorph(1);
        }
        auto end_fro = std::chrono::high_resolution_clock::now();
        std::cout << "fro time: " << std::chrono::duration<double>(end_fro - start_fro).count() << "s\n";
        // Ctxt temp = encryptedKeyStream[0];
        // auto start_add = std::chrono::high_resolution_clock::now();
        // for (int i=0;i<96;i++){
        //     temp += temp;
        // }
        // auto end_add = std::chrono::high_resolution_clock::now();
        // std::cout << "add time: " << std::chrono::duration<double>(end_add - start_add).count() << "s\n";

        // auto start_addconstant = std::chrono::high_resolution_clock::now();
        // for (long i = 0; i < 96; i++) // encrypt the encoded key
        // {
        //     temp.addConstant(encoded_expandedIV[0]);
        // }
        // auto end_addconstant = std::chrono::high_resolution_clock::now();
        // std::cout << "addConstant time: " << std::chrono::duration<double>(end_addconstant - start_addconstant).count() << "s\n";
        
        auto start_sbox = std::chrono::high_resolution_clock::now();
        auto start_linear = std::chrono::high_resolution_clock::now();
        auto end_sbox = std::chrono::high_resolution_clock::now();
        auto end_linear = std::chrono::high_resolution_clock::now();

        vector<Ctxt> in;
        in.resize(encryptedKeyStream.size(), Ctxt(ZeroCtxtLike, encryptedKeyStream[0]));
        for (long r = 1; r < Nr; r++)
        {
            std::cout << "Round " << r << " start" << std::endl;
            start_sbox = std::chrono::high_resolution_clock::now();
            // S Layer - 4 sbox
            #pragma omp parallel for
            for (long j = 0; j < 4; j++)
            {
                // for (int step = 0; step < 2; step++)
                // {
                //     decSboxFunc2(encryptedKeyStream, 4 * j, encA, allOne);
                // }
                for (int step = 0; step < 4; step++)
                {
                    mydecSboxFunc(encryptedKeyStream, 4 * j, encA, ea);
                }
            }
            end_sbox = std::chrono::high_resolution_clock::now();
            sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
            // for (unsigned t = 0; t < PlainBlock; t++)
            // {
            //     Vec<uint8_t> key(INIT_SIZE, BlockByte);
            //     memcpy(key.data(), &KeyStream[BlockByte * t], BlockByte);
            //     for (unsigned k = 0; k < 4; k++)
            //     {
            //         decSboxFi(key.data(), k * 4);
            //         decSboxFi(key.data(), k * 4);
            //         decSboxFi(key.data(), k * 4);
            //         decSboxFi(key.data(), k * 4);
            //     }
            //     memcpy(&KeyStream[BlockByte * t], key.data(), BlockByte);
            // }
            // return to natural PrimeSet to save memery
            // 创建 encryptedKeyStream[j] 的副本
            // for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
            // {
            //     tempCtxt[j] = encryptedKeyStream[j];
            //     tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
            // }
            // // 使用 verifyDecryption 函数解密并验证 KeyStream
            // if (!verifyDecryption16(tempCtxt, secretKey, ea, KeyStream))
            // {
            //     std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
            //     return false;
            // }
            // std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
            start_linear = std::chrono::high_resolution_clock::now();
            std::copy(encryptedKeyStream.begin(), encryptedKeyStream.end(), in.begin());
//#pragma omp parallel for
            for (long j = 0; j < BlockByte; j++)
            {
                encryptedKeyStream[j] += in[(j + 3) % 16];
                encryptedKeyStream[j] += in[(j + 4) % 16];
                encryptedKeyStream[j] += in[(j + 8) % 16];
                encryptedKeyStream[j] += in[(j + 9) % 16];
                encryptedKeyStream[j] += in[(j + 12) % 16];
                encryptedKeyStream[j] += in[(j + 14) % 16];
            }
            end_linear = std::chrono::high_resolution_clock::now();
            linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
            //
            // for (unsigned t = 0; t < PlainBlock; t++)
            // {
            //     Vec<uint8_t> key(INIT_SIZE, BlockByte);
            //     memcpy(key.data(), &KeyStream[BlockByte * t], BlockByte);
            //     decLinearLayer(key.data());
            //     memcpy(&KeyStream[BlockByte * t], key.data(), BlockByte);
            // }
            // return to natural PrimeSet to save memery
            // 创建 encryptedKeyStream[j] 的副本
            // for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
            // {
            //     tempCtxt[j] = encryptedKeyStream[j];
            //     tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
            // }
            // // 使用 verifyDecryption 函数解密并验证 KeyStream
            // if (!verifyDecryption16(tempCtxt, secretKey, ea, KeyStream))
            // {
            //     std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
            //     return false;
            // }
            // std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
            start_roundkey = std::chrono::high_resolution_clock::now();
//#pragma omp parallel for
for (long j = 0; j < BlockByte; j++)
            {
                encryptedKeyStream[j] += encryptedRoundKeySet[r * BlockByte+j];
            }
            end_roundkey = std::chrono::high_resolution_clock::now();
            roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
            //
            // for (long t = 0; t < PlainByte; t++)
            // {
            //     KeyStream[t] ^= RoundKeySet[r * PlainByte + t];
            // }
            // return to natural PrimeSet to save memery
            // 创建 encryptedKeyStream[j] 的副本
            // for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
            // {
            //     tempCtxt[j] = encryptedKeyStream[j];
            //     tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
            // }
            // 使用 verifyDecryption 函数解密并验证 KeyStream
            // if (!verifyDecryption16(tempCtxt, secretKey, ea, KeyStream))
            // {
            //     std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
            //     return false;
            // }
            // std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;
        }
        // 最后一轮
        std::cout << "th last round start" << std::endl;
        start_sbox = std::chrono::high_resolution_clock::now();
// S Layer - 4 sbox
#pragma omp parallel for
        for (long j = 0; j < 4; j++)
        {
            // for (int step = 0; step < 2; step++)
            // {
            //     decSboxFunc2(encryptedKeyStream, 4 * j, encA, allOne);
            // }
            for (int step = 0; step < 4; step++)
            {
                mydecSboxFunc(encryptedKeyStream, 4 * j, encA, ea);
            }
        }
        end_sbox = std::chrono::high_resolution_clock::now();
        sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
        // for (unsigned t = 0; t < PlainBlock; t++)
        // {
        //     Vec<uint8_t> key(INIT_SIZE, BlockByte);
        //     memcpy(key.data(), &KeyStream[BlockByte * t], BlockByte);
        //     for (unsigned k = 0; k < 4; k++)
        //     {
        //         decSboxFi(key.data(), k * 4);
        //         decSboxFi(key.data(), k * 4);
        //         decSboxFi(key.data(), k * 4);
        //         decSboxFi(key.data(), k * 4);
        //     }
        //     memcpy(&KeyStream[BlockByte * t], key.data(), BlockByte);
        // }
        // return to natural PrimeSet to save memery
        // 创建 encryptedKeyStream[j] 的副本
        // for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
        // {
        //     tempCtxt[j] = encryptedKeyStream[j];
        //     tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
        // }
        // 使用 verifyDecryption 函数解密并验证 KeyStream
        // if (!verifyDecryption16(tempCtxt, secretKey, ea, KeyStream))
        // {
        //     std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
        //     return false;
        // }
        // std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
        start_roundkey = std::chrono::high_resolution_clock::now();
//#pragma omp parallel for
        for (long j = 0; j < BlockByte; j++)
        {
            encryptedKeyStream[j] += encryptedRoundKeySet[Nr * BlockByte + j];
        }
        end_roundkey = std::chrono::high_resolution_clock::now();
        roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
        //
        // for (long t = 0; t < PlainByte; t++)
        // {
        //     KeyStream[t] ^= RoundKeySet[Nr * PlainByte + t];
        // }
        // return to natural PrimeSet to save memery
        // 创建 encryptedKeyStream[j] 的副本
        // for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
        // {
        //     tempCtxt[j] = encryptedKeyStream[j];
        //     tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
        // }
        // 使用 verifyDecryption 函数解密并验证 KeyStream
        // if (!verifyDecryption16(tempCtxt, secretKey, ea, KeyStream))
        // {
        //     std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
        //     return false;
        // }
        // std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;

        // 读取Client_KeyStream.txt，生成KeyStream
        vector<uint8_t> KeyStream(PlainByte);
        if (!readFromFile<uint8_t>(KeyStream.data(), "Client_KeyStream.txt", PlainByte))
        {
            return false;
        }
        std::cout << "KeyStream read succeeded!\n";
        // 使用 verifyDecryption 函数解密并验证 KeyStream
        if (!verifyDecryption16(encryptedKeyStream, secretKey, ea, KeyStream))
        {
            std::cerr << "Decryption verification failed for KeyStream." << std::endl;
            return false;
        }
        std::cout << "Decryption verification succeeded for KeyStream." << std::endl;

        double elapsed_seconds_KeyStream_FHE = roundkey_time + sbox_time + linear_layer_time;
        std::cout << "Round key addition total time: " << roundkey_time << "s\n";
        std::cout << "Sbox total time: " << sbox_time << "s\n";
        std::cout << "Linear layer total time: " << linear_layer_time << "s\n";
        std::cout << "KeyStream FHE Total Time: " << elapsed_seconds_KeyStream_FHE << "s\n";

        for (long i = 0; i < (long)encryptedKeyStream.size(); i++)
        {
            encryptedKeyStream[i].bringToSet(encryptedKeyStream[i].naturalPrimeSet());
        }
        // 写入 encryptedKeyStream 到文件
        std::ofstream outKeyStream("Server_encryptedKeyStream.bin", std::ios::binary);
        if (outKeyStream.is_open())
        {
            for (const auto &ctxt : encryptedKeyStream)
            {
                ctxt.writeTo(outKeyStream);
            }
            outKeyStream.close();
        }
        else
        {
            std::cerr << "Failed to open Server_encryptedKeyStream.bin for writing" << std::endl;
        }
        // 输出RoundKeySetFHE时间 + RoundConstantsFHE时间 + KeyStreamFHE时间
        double total_time_off = elapsed_seconds_RoundKeySet_FHE.count() + elapsed_seconds_RoundConstants_FHE.count() + elapsed_seconds_KeyStream_FHE;
        std::cout << "Decryption offline time (RoundKey+Constant+KeyStream): " << total_time_off << "s\n";
        return true;
    }

    bool Server_online()
    {
        // 读取Client_publickey，生成publicKey
#if 1
        // 从文件中读取上下文
        std::ifstream inContext("Client_context", std::ios::binary);
        if (!inContext.is_open())
        {
            std::cerr << "Failed to open Client_context for reading" << std::endl;
            throw std::runtime_error("Failed to open context file");
        }
        Context context = Context::readFrom(inContext);
        inContext.close();

        // 从文件中读取公钥
        std::ifstream inPublicKey("Client_publickey", std::ios::binary);
        if (!inPublicKey.is_open())
        {
            std::cerr << "Failed to open Client_publickey for reading" << std::endl;
            throw std::runtime_error("Failed to open public key file");
        }
        PubKey publicKey = PubKey::readFrom(inPublicKey, context);
        inPublicKey.close();
        printf("Public key read succeeded!\n");
        static const uint8_t YuxPolyBytes[] = {0x1B, 0x1}; // X^8+X^4+X^3+X+1
        const GF2X YuxPoly = GF2XFromBytes(YuxPolyBytes, 2);
        EncryptedArrayDerived<PA_GF2> ea(context, YuxPoly, context.getAlMod());
#endif

        // 读取Client_CipherStream.txt，生成cipherStream
#if 1
        vector<uint8_t> CipherStream(PlainByte);
        if (!readFromFile<uint8_t>(CipherStream.data(), "Client_CipherStream.txt", PlainByte))
        {
            return false;
        }
        std::cout << "CipherStream read succeeded!\n";
#endif

        // 读取Server_encryptedKeyStream.bin，生成encryptedKeyStream
#if 1
        vector<Ctxt> encryptedKeyStream;
        if (!readEncryptedKeyStream(encryptedKeyStream, "Server_encryptedKeyStream.bin", publicKey))
        {
            return false;
        }
        printf("Encrypted Key Stream read succeeded!\n");
#endif

        // 重点6：得到encryptedPlainStream
        vector<Ctxt> encryptedPlainStream;
        Ctxt tmpCtxt(publicKey);
        encryptedPlainStream.resize(BlockByte, tmpCtxt);

        auto start_encryptedPlainStream = std::chrono::steady_clock::now();
        encryptPlainStream(encryptedPlainStream, CipherStream, encryptedKeyStream, publicKey, ea);
        auto end_encryptedPlainStream = std::chrono::steady_clock::now();
        double elapsed_seconds_CipherStream_FHE = std::chrono::duration<double>(end_encryptedPlainStream - start_encryptedPlainStream).count();
        std::cout << "encryptedPlainStream FHE succeeded! Time: " << elapsed_seconds_CipherStream_FHE << "s\n";
// 读取CLient_PlainStream.txt，生成PlainStream
#if 1
        vector<uint8_t> PlainStream(PlainByte);
        if (!readFromFile<uint8_t>(PlainStream.data(), "Client_PlainStream.txt", PlainByte))
        {
            return false;
        }
        std::cout << "PlainStream read succeeded!\n";
        // 从文件中读取私钥
        std::ifstream inSecretKey("Client_secretkey", std::ios::binary);
        if (!inSecretKey.is_open())
        {
            std::cerr << "Failed to open Client_secretkey for reading" << std::endl;
            throw std::runtime_error("Failed to open secret key file");
        }
        SecKey secretKey = SecKey::readFrom(inSecretKey, context);
        inSecretKey.close();
        printf("Secret key read succeeded!\n");
#endif

        for (long i = 0; i < (long)encryptedPlainStream.size(); i++)
        {
            encryptedPlainStream[i].bringToSet(encryptedPlainStream[i].naturalPrimeSet());
        }
        if (!verifyDecryption16(encryptedPlainStream, secretKey, ea, PlainStream))
        {
            std::cerr << "Decryption verification failed for PlainStream." << std::endl;
            return false;
        }
        std::cout << "-----------------Check Server_decryption succeeded!-----------------" << std::endl;
        // 将encryptedPlainStream写入到文件Server_encryptedPlainStream.bin
#if 1
        std::ofstream out6("Server_encryptedPlainStream.bin", std::ios::binary);
        if (out6.is_open())
        {
            for (const auto &ctxt : encryptedPlainStream)
            {
                ctxt.writeTo(out6);
            }
            out6.close();
            printf("Encrypted PlainStream successfully written to Server_encryptedPlainStream.bin\n");
        }
        else
        {
            std::cerr << "Failed to open Server_encryptedPlainStream.bin for writing" << std::endl;
            return false;
        }
#endif

        // 输出
        double total_time_on = elapsed_seconds_CipherStream_FHE;
        std::cout << "Decryption online time: " << total_time_on << "s\n";
        return true;
    }
} // namespace Server_Yux2_8_C16