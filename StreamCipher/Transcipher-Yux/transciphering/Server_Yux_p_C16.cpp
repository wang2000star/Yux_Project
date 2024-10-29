#include "Server_Yux_p_C16.hpp"

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
// encryptXset函数既可以用于加密Xset, 也可以用于加密CipherStream
void encryptXset(vector<Ctxt>& encryptedXset, vector<uint64_t>& Xset, unique_ptr<PubKey>& he_pk, EncryptedArray& ea)
{   
    //编码
    vector<vector<long>> encoded;
    encodeTo16Ctxt_p(encoded, Xset, ea);
    //加密
    encryptedXset.resize(encoded.size(), Ctxt(*he_pk));
    for (long i=0; i<encoded.size(); i++){ // encrypt the encoded key
      ea.encrypt(encryptedXset[i], *he_pk, encoded[i]);
    }

}

void encryptIV(vector<Ctxt>& eKey, vector<uint64_t>& Key, unique_ptr<PubKey>& he_pk, EncryptedArray& ea)
{   
    //将key拷贝exKey
    vector<uint64_t> exKey(Key.size()*PlainBlock);
    for (size_t i = 0; i < Key.size(); i++)
    {
        for (size_t j = 0; j < PlainBlock; j++)
        {
            exKey[i*PlainBlock+j] = Key[i];
        }
    }
    //编码
    vector<vector<long>> encoded;
    encodeTo16Ctxt_p(encoded, exKey, ea);
    //加密
    eKey.resize(encoded.size(), Ctxt(*he_pk));
    for (long i=0; i<encoded.size(); i++){ // encrypt the encoded key
      ea.encrypt(eKey[i], *he_pk, encoded[i]);
    }
}

void buildRoundConstant(Ctxt& encA, const helib::EncryptedArray& ea, long roundConstant)
{
  // long --> ZZX -->Ctxt 
  vector<long> slots(ea.size(), roundConstant);
  ZZX ZZXConstant;
  ea.encode(ZZXConstant, slots);
  encA.DummyEncrypt(ZZXConstant);
}
// Compute the constants for Sbox
void decSboxFunc2(vector<Ctxt> &eData, long begin, Ctxt &encA)
{
    Ctxt c0(eData[begin]);
    Ctxt c1(eData[begin + 1]);
    Ctxt c2(eData[begin + 2]);
    Ctxt c3(eData[begin + 3]);
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

    eData[begin] = c2;
    eData[begin + 1] = c3;
    eData[begin + 2] = temp2;
    eData[begin + 3] = temp3;
}
// Server offline
namespace Server_Yux_p_C16
{
bool Server_offline()
{
    // 定义开始时间
    auto start = std::chrono::steady_clock::now();

    // 生成初始向量 IV
    vector<uint64_t> IV(BlockByte);
    for (unsigned i = 0; i < BlockByte; i++)
    {
        IV[i] = i + 1;
    }

    // 读取 Client_NonceSet.txt，生成 Xset
    vector<uint64_t> Xset(PlainByte * (Nr + 1));
    RandomBit<BlockSize> randomBit(Nr);
    auto &RanVecs = randomBit.roundconstants;
    vector<uint64_t> NonceSet(PlainBlock);
    if (!readFromFile<uint64_t>(NonceSet.data(), "Client_NonceSet.txt", PlainBlock))
    {
        std::cerr << "Failed to open Client_NonceSet.txt for reading" << std::endl;
        return false;
    }

    auto start_Xset = std::chrono::steady_clock::now();
    for (uint64_t counter = counter_begin; counter <= counter_end; counter++)
    {
        uint64_t nonce = NonceSet[counter - counter_begin];
        randomBit.generate_Instance_all_new(nonce, counter);
        for (unsigned r = 0; r <= Nr; r++)
        {
            Vec<uint64_t> X(INIT_SIZE, BlockByte);
            uint64_t temp;
            for (unsigned i = 0; i < BlockByte; ++i)
            {
                bool bit_array[Bytebits];
                for (unsigned j = 0; j < Bytebits; ++j)
                {
                    bit_array[j] = RanVecs[r][i * Bytebits + j];
                }
                BinStrToHex(bit_array, temp, Bytebits);
                X[i] = static_cast<uint64_t>(temp);
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
    // 创建一个空的 Context 对象
    auto context = std::make_shared<helib::Context>();
    context->readFrom(inContext);
    inContext.close();

    helib::EncryptedArray ea(context->getEA());

    std::ifstream inPublicKey("Client_publickey", std::ios::binary);
    if (!inPublicKey.is_open())
    {
        std::cerr << "Failed to open Client_publickey for reading" << std::endl;
        throw std::runtime_error("Failed to open public key file");
    }
     // 创建一个空的 PubKey 对象
    auto publicKey = std::make_unique<helib::PubKey>(*context);

    // 从文件中读取 PubKey 对象
    publicKey->readFrom(inPublicKey, *context);
    inPublicKey.close();

#if 1
    // 从文件中读取私钥
    std::ifstream inSecretKey("Client_secretkey", std::ios::binary);
    if (!inSecretKey.is_open())
    {
        std::cerr << "Failed to open Client_secretkey for reading" << std::endl;
        throw std::runtime_error("Failed to open secret key file");
    }
    SecKey secretKey = SecKey::readFrom(inSecretKey, *context);
    inSecretKey.close();
#endif

    // 读取对称密钥
    vector<Ctxt> encryptedSymKey;
    if (!readEncryptedSymKey(encryptedSymKey, "Client_encryptedSymKey.bin", *publicKey))
    {
        return false;
    }

    // 对 Xset 进行加密
    auto start_Xset_FHE = std::chrono::steady_clock::now();
    vector<Ctxt> encryptedXset;
    encryptXset(encryptedXset, Xset, publicKey, ea);
    auto end_Xset_FHE = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_Xset_FHE = end_Xset_FHE - start_Xset_FHE;
    std::cout << "Xset FHE succeeded! Time: " << elapsed_seconds_Xset_FHE.count() << "s\n";

    // 使用 verifyDecryption 函数解密并验证 Xset
    if (!verifyDecryption_p16(encryptedXset, Xset, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for Xset." << std::endl;
        return false;
    }
    std::cout << "Decryption verification succeeded for Xset." << std::endl;
    // 读取Client_RoundKeySet.txt，生成 RoundKeySet
    vector<uint64_t> RoundKeySet(PlainByte * (Nr + 1));
    if (!readFromFile<uint64_t>(RoundKeySet.data(), "Client_RoundKeySet.txt", PlainByte * (Nr + 1)))
    {
        std::cerr << "Failed to open Client_RoundKeySet.txt for reading" << std::endl;
        return false;
    }
    // 计算 encrypted_RoundKeySet
    auto start_RoundKeySet_FHE = std::chrono::steady_clock::now();
    vector<Ctxt> encrypted_RoundKeySet;
    Ctxt tmpCtxt(*publicKey);
    encrypted_RoundKeySet.resize(encryptedXset.size(), tmpCtxt);

    // vector<Ctxt> encrypted_RoundKeySet(encryptedXset.size(), Ctxt(ZeroCtxtLike, encryptedXset[0]));
    //omp_set_num_threads(16);
    //#pragma omp parallel for
    for (int r = 0; r < BlockByte * (Nr + 1); r++)
    {
        encrypted_RoundKeySet[r] = encryptedXset[r];
        // encrypted_RoundKeySet[r] += encryptedSymKey[0];
        encrypted_RoundKeySet[r].multiplyBy(encryptedSymKey[r % BlockByte]);
    }
    auto end_RoundKeySet_FHE = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_RoundKeySet_FHE = end_RoundKeySet_FHE - start_RoundKeySet_FHE;
    std::cout << "RoundKeySet FHE succeeded! Time: " << elapsed_seconds_RoundKeySet_FHE.count() << "s\n";
    // 使用 verifyDecryption 函数解密并验证 RoundKeySet
    if (!verifyDecryption_p16(encrypted_RoundKeySet, RoundKeySet, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for RoundKeySet." << std::endl;
        return false;
    }
    {
        std::cerr << "Decryption verification failed for RoundKeySet." << std::endl;
        return false;
    }
    std::cout << "Decryption verification succeeded for RoundKeySet." << std::endl;
    // 加密 IV
    auto start_IV_FHE = std::chrono::steady_clock::now();
    vector<Ctxt> encryptedIV;
    encryptIV(encryptedIV, IV, publicKey, ea);
    auto end_IV_FHE = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_IV_FHE = end_IV_FHE - start_IV_FHE;
    std::cout << "IV FHE succeeded! Time: " << elapsed_seconds_IV_FHE.count() << "s\n";
    // 使用 verifyDecryption 函数解密并验证 IV
    vector<uint64_t> expandedIV(PlainByte*BlockByte);
    for (unsigned i = 0; i < BlockByte; i++)
    {
        for (unsigned j = 0; j < PlainByte; j++)
        {
            expandedIV[i * PlainByte + j] = IV[i];
        }
    }
    if (!verifyDecryption_p16(encryptedIV, expandedIV, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for IV." << std::endl;
        return false;
    }
    std::cout << "Decryption verification succeeded for IV." << std::endl;

    // 加密 RoundConstants
    auto start_RoundConstants_FHE = std::chrono::steady_clock::now();
    Ctxt encA(ZeroCtxtLike, encryptedSymKey[0]);
    buildRoundConstant(encA, ea, RoundConstant);
    auto end_RoundConstants_FHE = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_RoundConstants_FHE = end_RoundConstants_FHE - start_RoundConstants_FHE;
    std::cout << "RoundConstants FHE succeeded! Time: " << elapsed_seconds_RoundConstants_FHE.count() << "s\n";

    // 生成 encryptedKeyStream
    // 定义roundkey_time、sbox_time、linear_layer_time
    double sbox_time = 0, linear_layer_time = 0, roundkey_time = 0;

    vector<Ctxt> encryptedKeyStream;
    encryptedKeyStream.resize(encryptedIV.size(), tmpCtxt);
    std::cout << "whiteround start" << std::endl;
    auto start_roundkey = std::chrono::steady_clock::now();
    for (long j = 0; j < (long)encryptedIV.size(); j++)
    {
        encryptedKeyStream[j] = encryptedIV[j];
        encryptedKeyStream[j] += encrypted_RoundKeySet[j]; // 初始轮密钥加
    }

    auto end_roundkey = std::chrono::steady_clock::now();
    roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
    // 明文密钥流
    vector<uint64_t> KeyStream(PlainByte);
    // 对IV和RoundKeySet进行模加
    for (long i = 0; i < PlainByte; i++)
    {
        KeyStream[i] = (IV[i % BlockByte] + RoundKeySet[i])%PlainMod;
    }

    // 使用 verifyDecryption 函数解密并验证 KeyStream
    if (!verifyDecryption_p16(encryptedKeyStream, KeyStream, secretKey, ea))
    {
        // std::cerr << "Decryption verification failed for KeyStream." << std::endl;
        return false;
    }
    std::cout << "Decryption verification succeeded for whiteround." << std::endl;

    // 创建副本
    vector<Ctxt> tempCtxt = encryptedKeyStream;
    for (long r = 1; r < Nr; r++)
    {
        std::cout << "Round " << r << " start" << std::endl;
        auto start_sbox = std::chrono::steady_clock::now();
        // S Layer - 4 sbox
        for (long j = 0; j < 4; j++)
        {
            for (int step = 0; step < 2; step++)
            {
                decSboxFunc2(encryptedKeyStream, 4 * j, encA);
            }
        }
        auto end_sbox = std::chrono::steady_clock::now();
        sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
        for (unsigned t = 0; t < PlainBlock; t++)
        {
            vector<uint64_t> key(BlockByte);
            memcpy(key.data(), &KeyStream[BlockByte * t], BlockByte);
            for (unsigned k = 0; k < 4; k++)
            {
                decSboxFi(key.data(), k * 4);
                decSboxFi(key.data(), k * 4);
                decSboxFi(key.data(), k * 4);
                decSboxFi(key.data(), k * 4);
            }
            memcpy(&KeyStream[BlockByte * t], key.data(), BlockByte);
        }
        // return to natural PrimeSet to save memery
        // 创建 encryptedKeyStream[j] 的副本
        for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
        {
            tempCtxt[j] = encryptedKeyStream[j];
            tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
        }
        // 使用 verifyDecryption 函数解密并验证 KeyStream
        if (!verifyDecryption_p16(tempCtxt, KeyStream, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
            return false;
        }
        std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
        auto start_linear = std::chrono::steady_clock::now();
        vector<Ctxt> in;
        in.resize(encryptedKeyStream.size(), Ctxt(ZeroCtxtLike, encryptedKeyStream[0]));
        for (long j = 0; j < encryptedKeyStream.size(); j++)
        {
            in[j] = encryptedKeyStream[j];
        }
        for (long j = 0; j < encryptedKeyStream.size(); j++)
        {
            encryptedKeyStream[j] += in[(j + 3) % 16];
            encryptedKeyStream[j] += in[(j + 4) % 16];
            encryptedKeyStream[j] += in[(j + 8) % 16];
            encryptedKeyStream[j] += in[(j + 9) % 16];
            encryptedKeyStream[j] += in[(j + 12) % 16];
            encryptedKeyStream[j] += in[(j + 14) % 16];
        }
        auto end_linear = std::chrono::steady_clock::now();
        linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
        //
        for (unsigned t = 0; t < PlainBlock; t++)
        {
            vector<uint64_t> key(BlockByte);
            memcpy(key.data(), &KeyStream[BlockByte * t], BlockByte);
            decLinearLayer(key.data());
            memcpy(&KeyStream[BlockByte * t], key.data(), BlockByte);
        }
        // return to natural PrimeSet to save memery
        // 创建 encryptedKeyStream[j] 的副本
        for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
        {
            tempCtxt[j] = encryptedKeyStream[j];
            tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
        }
        // 使用 verifyDecryption 函数解密并验证 KeyStream
        if (!verifyDecryption_p16(tempCtxt, KeyStream, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
            return false;
        }
        std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
        auto start_roundkey = std::chrono::steady_clock::now();
        for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
        {
            encryptedKeyStream[j] += encrypted_RoundKeySet[r * BlockByte + j];
        }

        auto end_roundkey = std::chrono::steady_clock::now();
        roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
        //
        for (long t = 0; t < PlainByte; t++)
        {
            KeyStream[t] ^= RoundKeySet[r * PlainByte + t];
        }
        // return to natural PrimeSet to save memery
        // 创建 encryptedKeyStream[j] 的副本
        for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
        {
            tempCtxt[j] = encryptedKeyStream[j];
            tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
        }
        // 使用 verifyDecryption 函数解密并验证 KeyStream
        if (!verifyDecryption_p16(tempCtxt, KeyStream, secretKey, ea))
        {
            std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
            return false;
        }

        std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;
    }
    // 最后一轮
    std::cout << "th last round start" << std::endl;
    auto start_sbox = std::chrono::steady_clock::now();
    // S Layer - 4 sbox
    for (long j = 0; j < 4; j++)
    {
        for (int step = 0; step < 2; step++)
        {
            decSboxFunc2(encryptedKeyStream, 4 * j, encA);
        }
    }
    auto end_sbox = std::chrono::steady_clock::now();
    sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
    for (unsigned t = 0; t < PlainBlock; t++)
    {
        vector<uint64_t> key(BlockByte);
        memcpy(key.data(), &KeyStream[BlockByte * t], BlockByte);
        for (unsigned k = 0; k < 4; k++)
        {
            decSboxFi(key.data(), k * 4);
            decSboxFi(key.data(), k * 4);
            decSboxFi(key.data(), k * 4);
            decSboxFi(key.data(), k * 4);
        }
        memcpy(&KeyStream[BlockByte * t], key.data(), BlockByte);
    }
    // return to natural PrimeSet to save memery
    // 创建 encryptedKeyStream[j] 的副本
    for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
    {
        tempCtxt[j] = encryptedKeyStream[j];
        tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
    }
    // 使用 verifyDecryption 函数解密并验证 KeyStream
    if (!verifyDecryption_p16(tempCtxt, KeyStream, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
        return false;
    }

    std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
    start_roundkey = std::chrono::steady_clock::now();
    for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
    {
        encryptedKeyStream[j] += encrypted_RoundKeySet[Nr * BlockByte + j];
    }

    end_roundkey = std::chrono::steady_clock::now();
    roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
    //
    for (long t = 0; t < PlainByte; t++)
    {
        KeyStream[t] ^= RoundKeySet[Nr * PlainByte + t];
    }
    // return to natural PrimeSet to save memery
    // 创建 encryptedKeyStream[j] 的副本
    for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
    {
        tempCtxt[j] = encryptedKeyStream[j];
        tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
    }
    // 使用 verifyDecryption 函数解密并验证 KeyStream
    if (!verifyDecryption_p16(tempCtxt, KeyStream, secretKey, ea))
    {
        std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
        return false;
    }
    std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;

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

    // 计算结束时间
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Server offline time: " << elapsed_seconds.count() << "s\n";

    // 输出Xse生成时间+XsetFHE时间+RoundKeySetFHE时间+IVFHE时间+RoundConstantsFHE时间+KeyStreamFHE时间
    total_time_off = elapsed_seconds_Xset.count() + elapsed_seconds_Xset_FHE.count() + elapsed_seconds_RoundKeySet_FHE.count() + elapsed_seconds_IV_FHE.count() + elapsed_seconds_RoundConstants_FHE.count() + elapsed_seconds_KeyStream_FHE;

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
    // 创建一个空的 Context 对象
    auto context = std::make_shared<helib::Context>();
    context->readFrom(inContext);
    inContext.close();

    // 从文件中读取公钥
    std::ifstream inPublicKey("Client_publickey", std::ios::binary);
    if (!inPublicKey.is_open())
    {
        std::cerr << "Failed to open Client_publickey for reading" << std::endl;
        throw std::runtime_error("Failed to open public key file");
    }
    // 创建一个空的 PubKey 对象
    auto publicKey = std::make_unique<helib::PubKey>(*context);

    // 从文件中读取 PubKey 对象
    publicKey->readFrom(inPublicKey, *context);
    inPublicKey.close();
    printf("Public key read succeeded!\n");
#endif

    // 读取Client_CipherStream.txt，生成cipherStream
#if 1
    vector<uint64_t> CipherStream(PlainByte);
    if (!readFromFile<uint64_t>(CipherStream.data(), "Client_CipherStream.txt", PlainByte))
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

    // 重点6：对cipherStream进行BGV公钥加密，得到encryptedCipherStream
#if 1
    // 定义CipherStream FHE开始时间
    auto start_CipherStream_FHE = std::chrono::steady_clock::now();
    // 定义ea
    helib::EncryptedArray ea(*context);

    vector<Ctxt> encryptedCipherStream;
    encryptXset(encryptedCipherStream, CipherStream, publicKey, ea);
    // 定义CipherStream FHE结束时间
    auto end_CipherStream_FHE = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_CipherStream_FHE = end_CipherStream_FHE - start_CipherStream_FHE;
    std::cout << "CipherStream FHE succeeded! Time: " << elapsed_seconds_CipherStream_FHE.count() << "s\n";
#endif

    // 重点7：对encryptedCipherStream和encryptedKeyStream进行BGV同态异或，得到encryptedPlainStream
#if 1
    // 定义PlainStream FHE开始时间
    auto start_PlainStream_FHE = std::chrono::steady_clock::now();
    vector<Ctxt> encryptedPlainStream;
    Ctxt tmpCtxt(*publicKey);
    encryptedPlainStream.resize(encryptedCipherStream.size(), tmpCtxt);
    for (long j = 0; j < (long)encryptedCipherStream.size(); j++)
    {
        encryptedPlainStream[j] = encryptedCipherStream[j];
        encryptedPlainStream[j] += encryptedKeyStream[j];
    }
    // 定义PlainStream FHE结束时间
    auto end_PlainStream_FHE = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_PlainStream_FHE = end_PlainStream_FHE - start_PlainStream_FHE;
    std::cout << "PlainStream FHE succeeded! Time: " << elapsed_seconds_PlainStream_FHE.count() << "s\n";
#endif

    // 重点8：对encryptedPlainStream进行解密，得到PlainStream
// 读取CLient_PlainStream.txt，生成PlainStream
#if 1
    vector<uint64_t> PlainStream(PlainByte);
    if (!readFromFile<uint64_t>(PlainStream.data(), "Client_PlainStream.txt", PlainByte))
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
    SecKey secretKey = SecKey::readFrom(inSecretKey, *context);
    inSecretKey.close();
    printf("Secret key read succeeded!\n");
#endif

    for (long i = 0; i < (long)encryptedPlainStream.size(); i++)
    {
        encryptedPlainStream[i].bringToSet(encryptedPlainStream[i].naturalPrimeSet());
    }
    if (!verifyDecryption_p16(encryptedPlainStream, PlainStream, secretKey, ea))
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
    // 定义结束时间
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Server online time: " << elapsed_seconds.count() << "s\n";

    // 输出CipherStream FHE时间+PlainStream FHE时间
    total_time_on = elapsed_seconds_CipherStream_FHE.count() + elapsed_seconds_PlainStream_FHE.count();
    std::cout << "Decryption offline time (Xset+RoundKey+IV+Constant+KeyStream): " << total_time_off << "s\n";
    std::cout << "Decryption online time(CipherStream+PlainStream): " << total_time_on << "s\n";
    double total_time = total_time_off + total_time_on;
    std::cout << "Decryption total time(offline+online): " << total_time << "s\n";
    return true;
}
}// namespace Server_Yux2_8_C16