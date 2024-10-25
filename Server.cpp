#include "Server.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

// encryptXset函数既可以用于加密Xset，也可以用于加密CipherStream
void encryptXset(vector<Ctxt> &encryptedXset, Vec<uint8_t> &Xset, const PubKey &hePK,
                 const EncryptedArrayDerived<PA_GF2> &ea)
{
    // Xset长度
    long length_s = BlockByte * (Nr + 1) * PlainBlock;
    // 明文槽数/单个分组的字节数，也就是每个Ctxt多项式打包表示了多少个分组
    long blocksPerCtxt = ea.size() / BlockByte;

    // 对expanded进行simd编码，这样会返回nRoundKeys个多项式数组即encoded，nRoundKeys=encoded.length()
    Vec<ZZX> encoded;
    encodeTo1Ctxt(encoded, Xset, ea); // encode as HE plaintext
    // 对encoded进行同态加密，得到eKey,
    {
        Ctxt tmpCtxt(hePK);
        encryptedXset.resize(encoded.length(), tmpCtxt);
    } // allocate space
    for (long i = 0; i < (long)encryptedXset.size(); i++) // encrypt the encoded key
        hePK.Encrypt(encryptedXset[i], encoded[i]);
}
// encryptSymIV函数既可以用于加密IV，也可以用于加密SymKey
void encryptSymIV(vector<Ctxt> &encryptedSymKey, Vec<uint8_t> &roundKeySchedule, const PubKey &hePK,
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
        encryptedSymKey.resize(encoded.length(), tmpCtxt);
    } // allocate space
    for (long i = 0; i < (long)encryptedSymKey.size(); i++) // encrypt the encoded key
        hePK.Encrypt(encryptedSymKey[i], encoded[i]);
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
bool Server_offline()
{
// 生成初始向量
#if 1
    Vec<uint8_t> IV(INIT_SIZE, BlockByte);
    for (unsigned i = 0; i < BlockByte; i++)
    {
        IV[i] = i + 1;
    }
#endif
// 读取Nonce.txt，生成Xset
#if 1
    Vec<uint8_t> Xset(INIT_SIZE, BlockByte * (Nr + 1) * PlainBlock);
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

        // 将RanVecs转换为X，比特转字节
        for (unsigned r = 0; r <= Nr; r++)
        {
            Vec<uint8_t> X(INIT_SIZE, BlockByte);
            for (unsigned i = 0; i < BlockByte; ++i)
            {
                bool bit_array[8];
                for (unsigned j = 0; j < 8; ++j)
                {
                    bit_array[j] = RanVecs[r][i * 8 + j];
                }
                BinStrToHex(bit_array, temp, 8);
                X[i] = static_cast<unsigned char>(temp);
            }
            // 把X存入X_set
            memcpy(&Xset[PlainByte * r + BlockByte * (counter - counter_begin)], X.data(), BlockByte);
        }
    }
    printf("Xset generation succeeded!\n");
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
    printf("Public key read succeeded!\n");
#endif
// 读取Client_encryptedSymKey.txt，生成encryptedSymKey
#if 1
    vector<Ctxt> encryptedSymKey;
    std::ifstream in3("Client_encryptedSymKey.txt", std::ios::binary);
    if (in3.is_open())
    {
        while (in3.peek() != EOF)
        {
            Ctxt ctxt(publicKey);
            ctxt.read(in3);
            encryptedSymKey.push_back(ctxt);
        }
        in3.close();
    }
    else
    {
        std::cerr << "Failed to open Client_encryptedSymKey.txt for reading" << std::endl;
    }
    printf("SymKey read succeeded!\n");
#endif
// 重点1：对Xset进行BGV公钥加密,得到encryptedXset
#if 1
    static const uint8_t YuxPolyBytes[] = {0x1B, 0x1}; // X^8+X^4+X^3+X+1
    const GF2X YuxPoly = GF2XFromBytes(YuxPolyBytes, 2);
    EncryptedArrayDerived<PA_GF2> ea(context, YuxPoly, context.getAlMod());
    vector<Ctxt> encryptedXset;
    encryptXset(encryptedXset, Xset, publicKey, ea);
    printf("Xset FHE succeeded!\n");
#endif
// 重点2：让encryptedXset的每一个元素分别和encryptedSymKey进行BGV同态异或，得到encrypted_RoundKeySet
#if 1
    vector<Ctxt> encrypted_RoundKeySet;
    encrypted_RoundKeySet.resize(encryptedXset.size());
    // 把encryptedXset的每一个元素copy到encrypted_RoundKeySet
    for (int r = 0; r < Nr + 1; r++)
    {
        encrypted_RoundKeySet[r] = encryptedXset[r];
        encrypted_RoundKeySet[r] *= encryptedSymKey[0];
    }
#endif
// 重点3：对IV进行BGV公钥加密,得到encryptedIV
#if 1
    // 生成同态加密的IVm
    vector<Ctxt> encryptedIV;
    encryptSymIV(encryptedIV, IV, publicKey, ea);

#endif
// 重点4：对roundConstants进行BGV公钥加密，得到encryptedRoundConstants
#if 1
    Ctxt encA(ZeroCtxtLike, encryptedSymKey[0]);
    buildRoundConstant(encA, ea); // Sbox里面使用的常数
    vector<PolyType> encLinTran;
    encLinTran.resize(3, DoubleCRT(ea.getContext(), ea.getContext().fullPrimes()));
    buildLinEnc2(encLinTran, ea);
#endif
/* 重点5：利用encrypted_RoundKeySet和encryptedIV，执行Yux加密/解密算法，得到encryptedKeyStream
   涉及的运算有比特异或、模多项式乘法、移位
   操作分为S盒、线性层、轮密钥加，共Nr轮，最后一轮没有线性层，有白化轮密钥加。
   核心点，使用simd打包批处理技术，提高单位时间的吞吐量
*/
#if 1
    vector<Ctxt> encryptedKeyStream;
    encryptedKeyStream.resize(encryptedIV.size());
    // 白化轮密钥加
    for (long j = 0; j < (long)encryptedIV.size(); j++)
    {

        encryptedKeyStream[j] = encryptedIV[j];
        encryptedKeyStream[j] += encrypted_RoundKeySet[0]; // initial key addition
    }

    // 第1到Nr-1轮
    for (long i = 1; i < Nr; i++)
    {
        for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
        {
            // S Layer
            for (long step = 0; step < 2; step++)
            {
                decSboxFunc2(encryptedKeyStream[j], encLinTran, encA, ea);
            }
            // Linear Layer
            Linear_function(encryptedKeyStream[j], ea);
            // Add round key
            encryptedKeyStream[j] += encrypted_RoundKeySet[i];
        }
    }
    // The last round is given below.
    // Linear layer is not here in the last round
    for (long j = 0; j < (long)encryptedKeyStream.size(); j++)
    {
        // S Layer
        for (long step = 0; step < 2; step++)
        {
            decSboxFunc2(encryptedKeyStream[j], encLinTran, encA, ea);
        }
        // Add round key
        encryptedKeyStream[j] += encrypted_RoundKeySet[Nr];
    }

    cout << "enc Finish! \n";
    // return to natural PrimeSet to save memery
    for (int i = 0; i < encryptedKeyStream.size(); i++)
        encryptedKeyStream[i].bringToSet(encryptedKeyStream[i].naturalPrimeSet());

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
    Vec<uint8_t> CipherStream(INIT_SIZE, BlockByte * PlainBlock);
    std::ifstream in4("Client_CipherStream.txt", std::ios::binary);
    if (in4.is_open())
    {
        while (in4.peek() != EOF)
        {
            CipherStream.append(in4.get());
        }
        in4.close();
    }
    else
    {
        std::cerr << "Failed to open Client_CipherStream.txt for reading" << std::endl;
    }
#endif
// 读取Server_encryptedKeyStream.txt，生成encryptedKeyStream
#if 1
    vector<Ctxt> encryptedKeyStream;
    std::ifstream in5("Server_encryptedKeyStream.txt", std::ios::binary);
    if (in5.is_open())
    {
        while (in5.peek() != EOF)
        {
            Ctxt ctxt(publicKey);
            ctxt.read(in5);
            encryptedKeyStream.push_back(ctxt);
        }
        in5.close();
    }
    else
    {
        std::cerr << "Failed to open Server_encryptedKeyStream.txt for reading" << std::endl;
    }
#endif
// 重点6：对cipherStream进行BGV公钥加密,得到encryptedCipherStream
#if 1
    static const uint8_t YuxPolyBytes[] = {0x1B, 0x1}; // X^8+X^4+X^3+X+1
    const GF2X YuxPoly = GF2XFromBytes(YuxPolyBytes, 2);
    EncryptedArrayDerived<PA_GF2> ea(context, YuxPoly, context.getAlMod());
    vector<Ctxt> encryptedCipherStream;
    encryptXset(encryptedCipherStream, CipherStream, publicKey, ea);

#endif
// 重点7：对encryptedCipherStream和encryptedKeyStream进行BGV同态异或，得到encryptedPlainStream
#if 1
    vector<Ctxt> encryptedPlainStream;
    encryptedPlainStream.resize(encryptedCipherStream.size());
    for (long j = 0; j < (long)encryptedCipherStream.size(); j++)
    {
        encryptedPlainStream[j] = encryptedCipherStream[j];
        encryptedPlainStream[j] += encryptedKeyStream[j];
    }
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