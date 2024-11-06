#include "Server_Yus_p_C32.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

static double total_time_off;
static double total_time_on;

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
        temp[index[0] + s] += eData[index[5]];
        temp[index[0] + s] += eData[index[6]];

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
        temp[index[3] + s] -= eData[index[5] + s];
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
        temp[index[0] + s] += eData[index[5]+s];
        temp[index[0] + s] += eData[index[6]+s];
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
        temp[index[3] + s] -= eData[index[5] + s];
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
            //c0.frobeniusAutomorph(1);
            c0.multiplyBy(c0);
            c0 += c1;
            c1 = temp;
        }

        eData[begin] = c0;
        eData[begin + 1] = c1;
    }
}
// Server offline
namespace Server_Yus_p_C32
{
    bool Server_offline()
    {
        // 定义开始时间
        auto start = std::chrono::steady_clock::now();

        // 生成初始向量 IV
        vector<long> IV(BlockByte);
        for (unsigned i = 0; i < BlockByte; i++)
        {
            IV[i] = i + 1;
        }

        // 读取 Client_NonceSet.txt
        vector<long> NonceSet(PlainBlock);
        if (!readFromFile<long>(NonceSet.data(), "Client_NonceSet.txt", PlainBlock))
        {
            std::cerr << "Failed to open Client_NonceSet.txt for reading" << std::endl;
            return false;
        }
        std::cout << "NonceSet read succeeded!\n";

        // 读取 Client_Xset.txt，生成 Xset
        vector<long> Xset(PlainByte * (Nr + 1));
        if (!readFromFile<long>(Xset.data(), "Client_Xset.txt", PlainByte * (Nr + 1)))
        {
            std::cerr << "Failed to open Client_Xset.txt for reading" << std::endl;
            return false;
        }
        std::cout << "Xset read succeeded!\n";

        // 读取context
        std::ifstream inContext("Client_context", std::ios::binary);
        if (!inContext.is_open())
        {
            std::cerr << "Failed to open Client_context for reading" << std::endl;
            throw std::runtime_error("Failed to open context file");
        }
        // 从文件中读取 Context 对象
        auto context = helib::Context::readFrom(inContext);
        inContext.close();
        std::cout << "Context read succeeded!\n";

        helib::EncryptedArray ea(context.getEA());

        std::ifstream inPublicKey("Client_publickey", std::ios::binary);
        if (!inPublicKey.is_open())
        {
            std::cerr << "Failed to open Client_publickey for reading" << std::endl;
            throw std::runtime_error("Failed to open public key file");
        }
        // 从文件中读取公钥
        auto publicKey = std::make_unique<helib::PubKey>(context);
        inPublicKey.close();
        std::cout << "Public key read succeeded!\n";

#if 1
        // 从文件中读取私钥
        std::ifstream inSecretKey("Client_secretkey", std::ios::binary);
        if (!inSecretKey.is_open())
        {
            std::cerr << "Failed to open Client_secretkey for reading" << std::endl;
            throw std::runtime_error("Failed to open secret key file");
        }
        SecKey secretKey = SecKey::readFrom(inSecretKey, context);
        inSecretKey.close();
        std::cout << "Secret key read succeeded!\n";
#endif

        // 读取对称密钥
        vector<Ctxt> encryptedSymKey;
        if (!readEncryptedSymKey(encryptedSymKey, "Client_encryptedSymKey.bin", *publicKey))
        {
            return false;
        }
        std::cout << "Symmetric key read succeeded!\n";
        vector<long> slotsData;
        ea.decrypt(encryptedSymKey[0], secretKey, slotsData);
        std::cout << "decryptedSymKey[0] = " << slotsData[0] << std::endl;
        slotsData.clear();
        encryptedSymKey[0].square();
        ea.decrypt(encryptedSymKey[0], secretKey, slotsData);
        std::cout << "decryptedSymKey[0]^2 = " << slotsData[0] << std::endl;

        // 读取Client_RoundKeySet.txt，生成 RoundKeySet
        vector<long> RoundKeySet(PlainByte * (Nr + 1));
        if (!readFromFile<long>(RoundKeySet.data(), "Client_RoundKeySet.txt", PlainByte * (Nr + 1)))
        {
            std::cerr << "Failed to open Client_RoundKeySet.txt for reading" << std::endl;
            return false;
        }
        std::cout << "RoundKeySet read succeeded!\n";

        vector<Ctxt> encryptedRoundKeySet;
        Ctxt tmpCtxt(*publicKey);
        encryptedRoundKeySet.resize(BlockByte * (Nr + 1), tmpCtxt);
        for (int i = 0; i < BlockByte * (Nr + 1); i++)
        {
            encryptedRoundKeySet[i] = encryptedSymKey[i % BlockByte];
        }

        // 创建副本
        // vector<Ctxt> temp;
        // temp.resize(encryptedSymKey.size(), tmpCtxt);
        // for (long i = 0; i < temp.size(); i++)
        // {
        //     temp[i] = encryptedSymKey[i];
        //     //temp[i].bringToSet(temp[i].naturalPrimeSet());
        // }
        // //测试
        // temp[0] += temp[0];
        // // 对temp[0]进行解密
        // vector<long> temp0;
        // ea.decrypt(temp[0], secretKey, temp0);
        // std::cout << "temp[0] = " << temp0[0] << std::endl;
        // std::cout << encryptedSymKey.size() << std::endl;


        vector<ZZX> encodedXset;
        encodeTo32Ctxt(encodedXset, Xset, ea); // encode as HE plaintext
        long encodedXset_length = BlockByte * (Nr + 1);
        // 计算 encryptedRoundKeySet
        auto start_RoundKeySet_FHE = std::chrono::steady_clock::now();
        // omp_set_num_threads(Nr+1);
        // #pragma omp parallel for
        for (long i = 0; i < encodedXset_length; i++) // encrypt the encoded key
        {
            encryptedRoundKeySet[i].addConstant(encodedXset[i]);
            //encryptedRoundKeySet[i].multByConstant(encodedXset[i]);
        }
        auto end_RoundKeySet_FHE = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds_RoundKeySet_FHE = end_RoundKeySet_FHE - start_RoundKeySet_FHE;
        std::cout << "RoundKeySet FHE succeeded! Time: " << elapsed_seconds_RoundKeySet_FHE.count() << "s\n";
        // 使用 verifyDecryption 函数解密并验证 RoundKeySet

        
        // if (!verifyDecryption32(temp, RoundKeySet,secretKey, ea))
        // {
        //     std::cerr << "Decryption verification failed for RoundKeySet." << std::endl;
        //     return false;
        // }
        // std::cout << "Decryption verification succeeded for RoundKeySet." << std::endl;

        // 生成 encryptedKeyStream
        // 定义roundkey_time、sbox_time、linear_layer_time
        double sbox_time = 0, linear_layer_time = 0, roundkey_time = 0;

        vector<Ctxt> encryptedKeyStream;
        encryptedKeyStream.resize(BlockByte, tmpCtxt);
        std::copy(encryptedRoundKeySet.begin(), encryptedRoundKeySet.begin() + BlockByte, encryptedKeyStream.begin());

        vector<long> expandedIV(BlockByte * PlainBlock);
        for (long j = 0; j < PlainBlock; j++)
        {
            memcpy(&expandedIV[BlockByte * j], &IV[0], BlockByte);
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
        // vector<long> KeyStream(PlainByte);
        // // 对IV和RoundKeySet进行异或
        // for (long i = 0; i < PlainByte; i++)
        // {
        //     KeyStream[i] = IV[i % BlockByte] ^ RoundKeySet[i];
        // }

        // // 使用 verifyDecryption 函数解密并验证 KeyStream
        // if (!verifyDecryption32(encryptedKeyStream, secretKey, ea, KeyStream))
        // {
        //     // std::cerr << "Decryption verification failed for KeyStream." << std::endl;
        //     return false;
        // }
        // std::cout << "Decryption verification succeeded for whiteround." << std::endl;

        // 创建副本
        // vector<Ctxt> tempCtxt = encryptedKeyStream;

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
            start_sbox = std::chrono::high_resolution_clock::now();
// S Layer - 4 sbox

                HE_Sbox(encryptedKeyStream);
            
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

            start_roundkey = std::chrono::high_resolution_clock::now();
            // #pragma omp parallel for
            for (long j = 0; j < BlockByte; j++)
            {
                encryptedKeyStream[j] += encryptedRoundKeySet[r * BlockByte + j];
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
        start_linear = std::chrono::high_resolution_clock::now();
        // #pragma omp parallel for
        //  MR Layer
        HE_MR(encryptedKeyStream);
        // MC Layer
        HE_MC(encryptedKeyStream);
        end_linear = std::chrono::high_resolution_clock::now();
        linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
        start_sbox = std::chrono::high_resolution_clock::now();
// S Layer - 16 sbox
        HE_Sbox(encryptedKeyStream);
        
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
        start_linear = std::chrono::high_resolution_clock::now();
        // #pragma omp parallel for
        //  MR Layer
        HE_MR(encryptedKeyStream);
        // MC Layer
        HE_MC(encryptedKeyStream);
        end_linear = std::chrono::high_resolution_clock::now();
        linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
        start_roundkey = std::chrono::high_resolution_clock::now();
        // #pragma omp parallel for
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
        vector<long> KeyStream(PlainByte);
        if (!readFromFile<long>(KeyStream.data(), "Client_KeyStream.txt", PlainByte))
        {
            return false;
        }
        std::cout << "KeyStream read succeeded!\n";
        // 使用 verifyDecryption 函数解密并验证 KeyStream
        if (!verifyDecryption32(encryptedKeyStream, KeyStream, secretKey, ea))
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
        double total_time_off = elapsed_seconds_RoundKeySet_FHE.count() + elapsed_seconds_KeyStream_FHE;
        std::cout << "Decryption offline time (RoundKey+Constant+KeyStream): " << total_time_off << "s\n";
        return true;
    }

    bool Server_online()
    {
        // 定义开始时间
        auto start = std::chrono::steady_clock::now();
        // 读取Client_publickey，生成publicKey
#if 1
        // 从文件中读取上下文
        std::ifstream inContext("Client_context", std::ios::binary);
        if (!inContext.is_open())
        {
            std::cerr << "Failed to open Client_context for reading" << std::endl;
            throw std::runtime_error("Failed to open context file");
        }
        auto context = helib::Context::readFrom(inContext);
        inContext.close();
        helib::EncryptedArray ea(context.getEA());
        // 从文件中读取公钥
        std::ifstream inPublicKey("Client_publickey", std::ios::binary);
        if (!inPublicKey.is_open())
        {
            std::cerr << "Failed to open Client_publickey for reading" << std::endl;
            throw std::runtime_error("Failed to open public key file");
        }
        // 创建一个空的 PubKey 对象
        auto publicKey = std::make_unique<helib::PubKey>(context);

        // 从文件中读取 PubKey 对象
        publicKey->readFrom(inPublicKey, context);
        inPublicKey.close();
        printf("Public key read succeeded!\n");
#endif

        // 读取Client_CipherStream.txt，生成cipherStream
#if 1
        vector<long> CipherStream(PlainByte);
        if (!readFromFile<long>(CipherStream.data(), "Client_CipherStream.txt", PlainByte))
        {
            return false;
        }
        std::cout << "CipherStream read succeeded!\n";
#endif

        // 读取Server_encryptedKeyStream.bin，生成encryptedKeyStream
#if 1
        vector<Ctxt> encryptedKeyStream;
        if (!readEncryptedKeyStream(encryptedKeyStream, "Server_encryptedKeyStream.bin", *publicKey))
        {
            return false;
        }
        printf("Encrypted Key Stream read succeeded!\n");
#endif

        // 重点6：得到encryptedPlainStream
        vector<Ctxt> encryptedPlainStream;
        Ctxt tmpCtxt(*publicKey);
        encryptedPlainStream.resize(BlockByte, tmpCtxt);
        std::copy(encryptedKeyStream.begin(), encryptedKeyStream.begin() + BlockByte, encryptedPlainStream.begin());
        vector<ZZX> encoded;
        encodeTo32Ctxt(encoded, CipherStream, ea);
        auto start_encryptedPlainStream = std::chrono::steady_clock::now();

#pragma omp parallel for
        for (long i = 0; i < (long)encryptedPlainStream.size(); i++) // encrypt the encoded key
        {
            encryptedPlainStream[i].addConstant(encoded[i]);
        }
        auto end_encryptedPlainStream = std::chrono::steady_clock::now();
        double elapsed_seconds_CipherStream_FHE = std::chrono::duration<double>(end_encryptedPlainStream - start_encryptedPlainStream).count();
        std::cout << "encryptedPlainStream FHE succeeded! Time: " << elapsed_seconds_CipherStream_FHE << "s\n";
// 读取CLient_PlainStream.txt，生成PlainStream
#if 1
        vector<long> PlainStream(PlainByte);
        if (!readFromFile<long>(PlainStream.data(), "Client_PlainStream.txt", PlainByte))
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
        if (!verifyDecryption32(encryptedPlainStream, PlainStream, secretKey, ea))
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
} // namespace Server_Yus_p_C32