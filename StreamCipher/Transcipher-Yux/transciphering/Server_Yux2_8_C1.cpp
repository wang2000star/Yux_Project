#include "Server_Yux2_8_C1.hpp"

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
// encryptXset函数既可以用于加密Xset，也可以用于加密CipherStream
void encryptXset(vector<Ctxt> &encryptedXset, Vec<uint8_t> &Xset, const PubKey &hePK,
                 const EncryptedArrayDerived<PA_GF2> &ea)
{
    // Xset长度
    long length_s = PlainByte * (Nr + 1);
    // 明文槽数/单个分组的字节数，也就是每个Ctxt多项式打包表示了多少个分组
    long blocksPerCtxt = ea.size() / BlockByte;

    // 对expanded进行simd编码，这样会返回nRoundKeys个多项式数组即encoded，nRoundKeys=encoded.length()
    Vec<ZZX> encoded;
    encodeTo1Ctxt(encoded, Xset, ea); // encode as HE plaintext
    {
        Ctxt tmpCtxt(hePK);
        encryptedXset.resize(encoded.length(), tmpCtxt);
    } // allocate space

    for (long i = 0; i < (long)encryptedXset.size(); i++) // encrypt the encoded key
    {
        hePK.Encrypt(encryptedXset[i], encoded[i]);
    }
}
// encryptSymIV函数既可以用于加密IV，也可以用于加密SymKey
void encryptSymIV(vector<Ctxt> &encryptedIV, Vec<uint8_t> &roundKeySchedule, const PubKey &hePK,
                  const EncryptedArrayDerived<PA_GF2> &ea)
{
    // roundKeySchedules是总轮密钥，已经拓展了，按照加密顺序。
    // 密钥总轮数
    long nRoundKeys = 1; // 如果每轮密钥都相同，这里就是1；如果每轮密钥都不同，这里就是总轮数Nr+1
    // 单轮密钥长度，这里就是字节数
    long round_key_length = BlockByte;
    // 总轮密钥长度，这里还是字节数，等于单轮密钥长度x密钥总轮数
    long length_s = round_key_length * nRoundKeys;
    // 明文槽数/单个分组的字节数，也就是每个Ctxt多项式打包表示了多少个分组
    long blocksPerCtxt = ea.size() / BlockByte;
    // Expand the key-schedule, copying each round key blocksPerCtxt times
    Vec<uint8_t> expanded(INIT_SIZE, nRoundKeys * blocksPerCtxt * BlockByte);
    for (long i = 0; i < nRoundKeys; i++)
    {
        uint8_t *roundKey = &roundKeySchedule[16 * i];
        for (long j = 0; j < blocksPerCtxt; j++)
            memcpy(&expanded[16 * (i * blocksPerCtxt + j)], roundKey, 16);
    }
    // 对expanded进行simd编码，这样会返回nRoundKeys个多项式数组即encoded，nRoundKeys=encoded.length()
    Vec<ZZX> encoded;
    encodeTo1Ctxt(encoded, expanded, ea); // encode as HE plaintext
                                          // 对encoded进行同态加密，得到eKey,
    {
        Ctxt tmpCtxt(hePK);
        encryptedIV.resize(encoded.length(), tmpCtxt);
    } // allocate space

    for (long i = 0; i < (long)encryptedIV.size(); i++) // encrypt the encoded key
    {
        hePK.Encrypt(encryptedIV[i], encoded[i]);
    }
}
void buildRoundConstant(Ctxt &encA,
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
// Compute the constants for Sbox
void buildLinEnc(vector<PolyType> &encLinTran,
                 const EncryptedArrayDerived<PA_GF2> &ea)
{
    // encLinTran[0]: The constants only have nonzero entires in their slots corresponding
    // to bytes 3,7,11,15 of each blocks
    // encLinTran[1]: The constants only have zero entires in their slots corresponding
    // to bytes 0,4,8,12 of each blocks, others is 1;

    Vec<uint8_t> bytes(INIT_SIZE, ea.size());
    long blocksPerCtxt = ea.size() / 16;
    Vec<ZZX> tmp;

    memset(bytes.data(), 0, bytes.length());
    /*
      void Transcipher1::*memset(void Transcipher1::*s, int ch, size_t n);
      函数解释：将s中前n个字节 （typedef unsigned int size_t）用 ch 替换并返回 s
      讲bytes设置为0
    */
    for (long j = 0; j < blocksPerCtxt; j++)
    {
        uint8_t *bptr = &bytes[16 * j];
        bptr[3] = bptr[7] = bptr[11] = bptr[15] = 1;
    }
    encodeTo1Ctxt(tmp, bytes, ea);
    encLinTran[0] = tmp[0];

    memset(bytes.data(), 1, bytes.length());
    for (long j = 0; j < blocksPerCtxt; j++)
    {
        uint8_t *bptr = &bytes[16 * j];
        bptr[3] = bptr[7] = bptr[11] = bptr[15] = 0;
    }
    encodeTo1Ctxt(tmp, bytes, ea);
    encLinTran[1] = tmp[0];
}
// Compute the constants for Sbox
void buildLinEnc2(vector<PolyType> &encLinTran,
                  const EncryptedArrayDerived<PA_GF2> &ea)
{
    // encLinTran[0]: The constants only have nonzero entires in their slots corresponding
    // to bytes 3,7,11,15 of each blocks // 0001000100010001
    // encLinTran[1]: The constants only have nonzero entires in their slots corresponding
    // to bytes 2,6,10,14 of each blocks // 0010001000100010
    // encLinTran[1]: The constants only have zero entires in their slots corresponding
    // to bytes 01,45,89,1213 of each blocks, others is 1; //1100110011001100

    Vec<uint8_t> bytes(INIT_SIZE, ea.size());
    long blocksPerCtxt = ea.size() / 16;
    Vec<ZZX> tmp;

    memset(bytes.data(), 0, bytes.length());
    /*
      void Transcipher1::*memset(void Transcipher1::*s, int ch, size_t n);
      函数解释：将s中前n个字节 （typedef unsigned int size_t）用 ch 替换并返回 s
      讲bytes设置为0
    */
    for (long j = 0; j < blocksPerCtxt; j++)
    {
        uint8_t *bptr = &bytes[16 * j];
        bptr[3] = bptr[7] = bptr[11] = bptr[15] = 1;
    }
    encodeTo1Ctxt(tmp, bytes, ea);
    encLinTran[0] = tmp[0];

    memset(bytes.data(), 0, bytes.length());
    for (long j = 0; j < blocksPerCtxt; j++)
    {
        uint8_t *bptr = &bytes[16 * j];
        bptr[2] = bptr[6] = bptr[10] = bptr[14] = 1;
    }
    encodeTo1Ctxt(tmp, bytes, ea);
    encLinTran[1] = tmp[0];

    memset(bytes.data(), 1, bytes.length());
    for (long j = 0; j < blocksPerCtxt; j++)
    {
        uint8_t *bptr = &bytes[16 * j];
        bptr[3] = bptr[7] = bptr[11] = bptr[15] = 0;
        bptr[2] = bptr[6] = bptr[10] = bptr[14] = 0;
    }
    encodeTo1Ctxt(tmp, bytes, ea);
    encLinTran[2] = tmp[0];
}
// Sbox function
void decSboxFunc(Ctxt &c, vector<PolyType> encLinTran, Ctxt &encA, const EncryptedArrayDerived<PA_GF2> &ea)
{
    // The basic rotation amount along the 1st dimension
    long rotAmount = ea.getContext().getZMStar().OrderOf(0) / 16;

    c.cleanUp();
    Ctxt c1(c), c2(c), c3(c), c4(c);
    ea.rotate1D(c1, 0, 1 * rotAmount);
    ea.rotate1D(c2, 0, 2 * rotAmount);
    ea.rotate1D(c3, 0, 3 * rotAmount);
    ea.rotate1D(c4, 0, 15 * rotAmount);
    c1.cleanUp();
    c2.cleanUp();
    c3.cleanUp();
    c4.cleanUp();

    c1.multiplyBy(c2);
    c += c1;
    c += c3;
    c += encA;

    const PolyType &p4 = encLinTran[0];   // 0001000100010001
    const PolyType &p123 = encLinTran[1]; // 1110111011101110

    c.multByConstant(p4);
    c4.multByConstant(p123);
    c += c4;

    c.cleanUp();
}
// Sbox function
void decSboxFunc2(Ctxt &c, vector<PolyType> encLinTran, Ctxt &encA, const EncryptedArrayDerived<PA_GF2> &ea)
{
    // The basic rotation amount along the 1st dimension
    long rotAmount = ea.getContext().getZMStar().OrderOf(0) / 16;

    const PolyType &p3 = encLinTran[0];  // 0001000100010001
    const PolyType &p2 = encLinTran[1];  // 0010001000100010
    const PolyType &p01 = encLinTran[2]; // 1100110011001100

    c.cleanUp();
    Ctxt c1(c), c2(c), c15(c);
    ea.shift1D(c1, 0, 1 * rotAmount);
    ea.shift1D(c2, 0, 2 * rotAmount);
    ea.rotate1D(c15, 0, 15 * rotAmount);
    c1.cleanUp();
    c2.cleanUp();
    c15.cleanUp();

    c1.multiplyBy(c);
    c1 += c2;
    c1 += encA;

    Ctxt y3(c1);

    y3.multByConstant(p3);

    c15 += c1;
    c15.multByConstant(p2);
    Ctxt y2(c15);
    ea.shift1D(c15, 0, 1 * rotAmount);
    c15.cleanUp();
    y3 += c15;

    ea.rotate1D(c, 0, 14 * rotAmount);
    c.cleanUp();
    c.multByConstant(p01);
    c += y2;
    c += y3;
    c.cleanUp();
}
// Linear layer
void Linear_function(Ctxt &c, const EncryptedArrayDerived<PA_GF2> &ea)
{
    // The basic rotation amount along the 1st dimension
    long rotAmount = ea.getContext().getZMStar().OrderOf(0) / 16;

    c.cleanUp();
    // 循环左移   3  4  8 9 12 14
    // 即循环右移 13 12 8 7  4  2
    Ctxt c3(c), c4(c), c8(c), c9(c), c12(c), c14(c);
    ea.rotate1D(c3, 0, 13 * rotAmount);
    ea.rotate1D(c4, 0, 12 * rotAmount);
    ea.rotate1D(c8, 0, 8 * rotAmount);
    ea.rotate1D(c9, 0, 7 * rotAmount);
    ea.rotate1D(c12, 0, 4 * rotAmount);
    ea.rotate1D(c14, 0, 2 * rotAmount);

    c3.cleanUp();
    c4.cleanUp();
    c8.cleanUp();
    c9.cleanUp();
    c12.cleanUp();
    c14.cleanUp();

    c += c3;
    c += c4;
    c += c8;
    c += c9;
    c += c12;
    c += c14;
    c.cleanUp();
}

// Server offline
namespace Server_Yux2_8_C1
{
bool Server_offline()
{
    // 定义开始时间
    auto start = std::chrono::steady_clock::now();

    // 生成初始向量 IV
    Vec<uint8_t> IV(INIT_SIZE, BlockByte);
    for (unsigned i = 0; i < BlockByte; i++)
    {
        IV[i] = i + 1;
    }

    // 读取 Client_NonceSet.txt，生成 Xset
    Vec<uint8_t> Xset(INIT_SIZE, PlainByte * (Nr + 1));
    RandomBit<BlockSize> randomBit(Nr);
    auto &RanVecs = randomBit.roundconstants;
    Vec<uint64_t> NonceSet(INIT_SIZE, PlainBlock);
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
            Vec<uint8_t> X(INIT_SIZE, BlockByte);
            uint64_t temp;
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
#endif

    // 读取对称密钥
    vector<Ctxt> encryptedSymKey;
    if (!readEncryptedSymKey(encryptedSymKey, "Client_encryptedSymKey.bin", publicKey))
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
    if (!verifyDecryption1(encryptedXset, secretKey, ea, Xset))
    {
        std::cerr << "Decryption verification failed for Xset." << std::endl;
        return false;
    }
    std::cout << "Decryption verification succeeded for Xset." << std::endl;
    // 读取Client_RoundKeySet.txt，生成 RoundKeySet
    Vec<uint8_t> RoundKeySet(INIT_SIZE, PlainByte * (Nr + 1));
    if (!readFromFile<uint8_t>(RoundKeySet.data(), "Client_RoundKeySet.txt", PlainByte * (Nr + 1)))
    {
        std::cerr << "Failed to open Client_RoundKeySet.txt for reading" << std::endl;
        return false;
    }
    // 计算 encrypted_RoundKeySet
    auto start_RoundKeySet_FHE = std::chrono::steady_clock::now();
    vector<Ctxt> encrypted_RoundKeySet;
    Ctxt tmpCtxt(publicKey);
    encrypted_RoundKeySet.resize(encryptedXset.size(), tmpCtxt);

    // vector<Ctxt> encrypted_RoundKeySet(encryptedXset.size(), Ctxt(ZeroCtxtLike, encryptedXset[0]));
    for (int r = 0; r < Nr + 1; r++)
    {
        encrypted_RoundKeySet[r] = encryptedXset[r];
        // encrypted_RoundKeySet[r] += encryptedSymKey[0];
        encrypted_RoundKeySet[r].multiplyBy(encryptedSymKey[0]);
    }
    auto end_RoundKeySet_FHE = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_RoundKeySet_FHE = end_RoundKeySet_FHE - start_RoundKeySet_FHE;
    std::cout << "RoundKeySet FHE succeeded! Time: " << elapsed_seconds_RoundKeySet_FHE.count() << "s\n";
    // 使用 verifyDecryption 函数解密并验证 RoundKeySet
    if (!verifyDecryption1(encrypted_RoundKeySet, secretKey, ea, RoundKeySet))
    {
        std::cerr << "Decryption verification failed for RoundKeySet." << std::endl;
        return false;
    }
    std::cout << "Decryption verification succeeded for RoundKeySet." << std::endl;
    // 加密 IV
    auto start_IV_FHE = std::chrono::steady_clock::now();
    vector<Ctxt> encryptedIV;
    encryptSymIV(encryptedIV, IV, publicKey, ea);
    auto end_IV_FHE = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_IV_FHE = end_IV_FHE - start_IV_FHE;
    std::cout << "IV FHE succeeded! Time: " << elapsed_seconds_IV_FHE.count() << "s\n";
    // 使用 verifyDecryption 函数解密并验证 IV
    Vec<uint8_t> expandedIV(INIT_SIZE, PlainByte);
    for (long i = 0; i < 1; i++)
    {
        for (long j = 0; j < PlainBlock; j++)
            memcpy(&expandedIV[16 * (i * PlainBlock + j)], &IV[BlockByte * i], BlockByte);
    }
    if (!verifyDecryption1(encryptedIV, secretKey, ea, expandedIV))
    {
        std::cerr << "Decryption verification failed for IV." << std::endl;
        return false;
    }
    std::cout << "Decryption verification succeeded for IV." << std::endl;
    // 加密 RoundConstants
    auto start_RoundConstants_FHE = std::chrono::steady_clock::now();
    Ctxt encA(ZeroCtxtLike, encryptedSymKey[0]);
    buildRoundConstant(encA, ea);
    vector<PolyType> encLinTran(3, DoubleCRT(ea.getContext(), ea.getContext().fullPrimes()));
    buildLinEnc2(encLinTran, ea);
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
        encryptedKeyStream[j] += encrypted_RoundKeySet[0]; // 初始轮密钥加
    }
    auto end_roundkey = std::chrono::steady_clock::now();
    roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
    // 明文密钥流
    Vec<uint8_t> KeyStream(INIT_SIZE, PlainByte);
    // 对IV和RoundKeySet进行异或
    for (long i = 0; i < PlainByte; i++)
    {
        KeyStream[i] = IV[i % BlockByte] ^ RoundKeySet[i];
    }

    // 使用 verifyDecryption 函数解密并验证 KeyStream
    if (!verifyDecryption1(encryptedKeyStream, secretKey, ea, KeyStream))
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
        for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
        {
            auto start_sbox = std::chrono::steady_clock::now();
            for (long step = 0; step < 2; step++)
            {
                decSboxFunc2(encryptedKeyStream[j], encLinTran, encA, ea);
            }
            auto end_sbox = std::chrono::steady_clock::now();
            sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
            for (unsigned t = 0; t < PlainBlock; t++)
            {
                Vec<uint8_t> key(INIT_SIZE, BlockByte);
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
            tempCtxt[j] = encryptedKeyStream[j];
            tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
            // 使用 verifyDecryption 函数解密并验证 KeyStream
            if (!verifyDecryption1(tempCtxt, secretKey, ea, KeyStream))
            {
                std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
                return false;
            }
            std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
            auto start_linear = std::chrono::steady_clock::now();
            Linear_function(encryptedKeyStream[j], ea);
            auto end_linear = std::chrono::steady_clock::now();
            linear_layer_time += std::chrono::duration<double>(end_linear - start_linear).count();
            //
            for (unsigned t = 0; t < PlainBlock; t++)
            {
                Vec<uint8_t> key(INIT_SIZE, BlockByte);
                memcpy(key.data(), &KeyStream[BlockByte * t], BlockByte);
                decLinearLayer(key.data());
                memcpy(&KeyStream[BlockByte * t], key.data(), BlockByte);
            }
            // return to natural PrimeSet to save memery
            // 创建 encryptedKeyStream[j] 的副本
            tempCtxt[j] = encryptedKeyStream[j];
            tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
            // 使用 verifyDecryption 函数解密并验证 KeyStream
            if (!verifyDecryption1(tempCtxt, secretKey, ea, KeyStream))
            {
                std::cerr << "Decryption verification failed for KeyStream Linear Layer." << std::endl;
                return false;
            }
            std::cout << "Decryption verification succeeded for KeyStream Linear Layer." << std::endl;
            auto start_roundkey = std::chrono::steady_clock::now();
            encryptedKeyStream[j] += encrypted_RoundKeySet[r];
            auto end_roundkey = std::chrono::steady_clock::now();
            roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
            //
            for (long t = 0; t < PlainByte; t++)
            {
                KeyStream[t] ^= RoundKeySet[r * PlainByte + t];
            }
            // return to natural PrimeSet to save memery
            // 创建 encryptedKeyStream[j] 的副本
            tempCtxt[j] = encryptedKeyStream[j];
            tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
            // 使用 verifyDecryption 函数解密并验证 KeyStream
            if (!verifyDecryption1(tempCtxt, secretKey, ea, KeyStream))
            {
                std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
                return false;
            }
            std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;
        }
    }
    // 最后一轮
    std::cout << "th last round start" << std::endl;
    for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
    {
        auto start_sbox = std::chrono::steady_clock::now();
        for (long step = 0; step < 2; step++)
        {
            decSboxFunc2(encryptedKeyStream[j], encLinTran, encA, ea);
        }
        auto end_sbox = std::chrono::steady_clock::now();
        sbox_time += std::chrono::duration<double>(end_sbox - start_sbox).count();
        for (unsigned t = 0; t < PlainBlock; t++)
        {
            Vec<uint8_t> key(INIT_SIZE, BlockByte);
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
        tempCtxt[j] = encryptedKeyStream[j];
        tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
        // 使用 verifyDecryption 函数解密并验证 KeyStream
        if (!verifyDecryption1(tempCtxt, secretKey, ea, KeyStream))
        {
            std::cerr << "Decryption verification failed for KeyStream Sbox." << std::endl;
            return false;
        }
        std::cout << "Decryption verification succeeded for KeyStream Sbox." << std::endl;
        auto start_roundkey = std::chrono::steady_clock::now();
        encryptedKeyStream[j] += encrypted_RoundKeySet[Nr];
        auto end_roundkey = std::chrono::steady_clock::now();
        roundkey_time += std::chrono::duration<double>(end_roundkey - start_roundkey).count();
        //
        for (long t = 0; t < PlainByte; t++)
        {
            KeyStream[t] ^= RoundKeySet[Nr * PlainByte + t];
        }
        // 创建 encryptedKeyStream[j] 的副本
        tempCtxt[j] = encryptedKeyStream[j];
        tempCtxt[j].bringToSet(tempCtxt[j].naturalPrimeSet());
        // 使用 verifyDecryption 函数解密并验证 KeyStream
        if (!verifyDecryption1(tempCtxt, secretKey, ea, KeyStream))
        {
            std::cerr << "Decryption verification failed for KeyStream Round Key Addition." << std::endl;
            return false;
        }
        std::cout << "Decryption verification succeeded for KeyStream Round Key Addition." << std::endl;
    }
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
#endif

    // 读取Client_CipherStream.txt，生成cipherStream
#if 1
    Vec<uint8_t> CipherStream(INIT_SIZE, PlainByte);
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

    // 重点6：对cipherStream进行BGV公钥加密，得到encryptedCipherStream
#if 1
    // 定义CipherStream FHE开始时间
    auto start_CipherStream_FHE = std::chrono::steady_clock::now();
    static const uint8_t YuxPolyBytes[] = {0x1B, 0x1}; // X^8+X^4+X^3+X+1
    const GF2X YuxPoly = GF2XFromBytes(YuxPolyBytes, 2);
    EncryptedArrayDerived<PA_GF2> ea(context, YuxPoly, context.getAlMod());
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
    Ctxt tmpCtxt(publicKey);
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
    Vec<uint8_t> PlainStream(INIT_SIZE, PlainByte);
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
    if (!verifyDecryption1(encryptedPlainStream, secretKey, ea, PlainStream))
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
} // namespace Server_Yux2_8_C1