#include "Client.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

// 通用写入函数模板
template <typename T>
bool writeToFile(const Vec<T> &data, const string &filename, size_t length)
{
    std::ofstream out(filename);
    if (out.is_open())
    {
        for (size_t i = 0; i < length; ++i)
        {
            out << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(data[i]) << " ";
        }
        out.close();
        return true;
    }
    else
    {
        std::cerr << "Failed to open " << filename << " for writing" << std::endl;
        return false;
    }
}

// 通用读取函数模板
template <typename T>
bool readFromFile(Vec<T> &data, const string &filename, size_t length)
{
    std::ifstream in(filename);
    if (!in.is_open())
    {
        std::cerr << "Failed to open " << filename << " for reading" << std::endl;
        return false;
    }

    // 获取文件大小
    in.seekg(0, std::ios::end);
    std::streamsize size = in.tellg();
    in.seekg(0, std::ios::beg);

    // 计算预期的文件大小
    std::streamsize expected_size = length * (sizeof(T) * 2 + 1) - 1; // 每个字节占用2个字符和1个空格，最后一个字节没有空格

    if (size != expected_size)
    {
        std::cerr << "File size does not match expected length. Expected: " << expected_size << ", Got: " << size << std::endl;
        return false;
    }

    for (size_t i = 0; i < length; ++i)
    {
        int temp;
        in >> std::hex >> temp;
        data[i] = static_cast<T>(temp);
    }

    in.close();
    return true;
}

// run the Yux key-expansion and then encrypt the expanded key.
void encryptSymKey(vector<Ctxt> &encryptedSymKey, Vec<uint8_t> &roundKeySchedule, const PubKey &hePK,
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

bool Client_offline()
{
    // 定义对称密钥SymKey,生成随机数集合nonce_set,Xset,RoundKeySet和密钥流KeyStream
#if 1
    // 生成初始向量
    vector<uint8_t> iv(BlockByte);
    for (unsigned i = 0; i < BlockByte; i++)
    {
        iv[i] = i + 1;
    }
    // 生成对称密钥
    GF2X rnd;
    Vec<uint8_t> SymKey(INIT_SIZE, BlockByte); // 8*BlockByte
    random(rnd, 8 * SymKey.length());
    BytesFromGF2X(SymKey.data(), rnd, SymKey.length());
    // 定义随机数集合nonce_set，随机向量Xset,轮密钥RoundKeySet和密钥流KeyStream
    Vec<uint64_t> NonceSet(INIT_SIZE, PlainBlock);
    Vec<uint8_t> Xset(INIT_SIZE, PlainByte * (Nr + 1));
    Vec<uint8_t> RoundKeySet(INIT_SIZE, PlainByte * (Nr + 1));
    Vec<uint8_t> KeyStream(INIT_SIZE, PlainByte);

    // 创建RandomBit对象
    RandomBit<BlockSize> randomBit(Nr);
    // 为roundconstants创建别名
    auto &RanVecs = randomBit.roundconstants;
    // 定义counter_begin和counter_end
    uint64_t counter_begin = 0;                            // 计数器起始值
    uint64_t counter_end = PlainBlock + counter_begin - 1; // 计数器结束值

    for (uint64_t counter = counter_begin; counter <= counter_end; counter++)
    {
        uint64_t nonce = generate_secure_random_int(NonceSize);
        NonceSet[counter - counter_begin] = nonce;
        randomBit.generate_Instance_all_new(nonce, counter); // nonce, counter
        // 将RanVecs转换为X，比特转字节
        Vec<uint8_t> state(INIT_SIZE, BlockByte); // 单组密钥流状态
        for (unsigned r = 0; r <= Nr; r++)
        {
            Vec<uint8_t> X(INIT_SIZE, BlockByte);
            Vec<uint8_t> RoundKey(INIT_SIZE, BlockByte);
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
                RoundKey[i] = mul(SymKey[i], X[i]);
            } // 把X存入X_set
            memcpy(&Xset[PlainByte * r + BlockByte * (counter - counter_begin)], X.data(), BlockByte);
            memcpy(&RoundKeySet[BlockByte * r + BlockByte * (counter - counter_begin)], RoundKey.data(), BlockByte);
            if (r == 0)
            {
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] = RoundKey[i] ^ iv[i];
                }
            }
            else if (r < Nr)
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
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] ^= RoundKey[i];
                }
            }
            else
            {
                for (unsigned i = 0; i < 4; i++)
                {
                    encSboxFi(state.data(), i * 4);
                    encSboxFi(state.data(), i * 4);
                    encSboxFi(state.data(), i * 4);
                    encSboxFi(state.data(), i * 4);
                }
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    state[i] ^= RoundKey[i];
                }
                memcpy(&KeyStream[(counter - counter_begin) * BlockByte], state.data(), BlockByte);
            }
        }
    }

    // 将Nonce以字节形式写入到文件Client_NonceSet.txt
    if (!writeToFile(NonceSet, "Client_NonceSet.txt", PlainBlock))
    {
        return false;
    }
    // 将X_set写入到文件Client_Xset.txt
    if (!writeToFile(Xset, "Client_Xset.txt", PlainByte * (Nr + 1)))
    {
        return false;
    }
    // 将RoundKey_set写入到文件Client_RoundKeySet.txt
    if (!writeToFile(RoundKeySet, "Client_RoundKeySet.txt", PlainByte * (Nr + 1)))
    {
        return false;
    }
    // 将KeyStream以字节形式写入到文件Client_KeyStream.txt
    if (!writeToFile(KeyStream, "Client_KeyStream.txt", PlainByte))
    {
        return false;
    }

    printf("KeyStream Generation Succeeded!\n");
#endif
    // 同态加密初始对称密钥SymKey,得到encryptedSymKey
    // 设置
    long idx = 3; // 0
    // amap.arg("sz", idx, "parameter-sets: toy=0 through huge=5");
    long c = 9;
    // amap.arg("c", c, "number of columns in the key-switching matrices");
    bool packed = true;
    // amap.arg("packed", packed, "use packed bootstrapping");
    // amap.parse(argc, argv);
    if (idx > 5)
        idx = 5;
    // lm: I have set the parameters of idx = 0, 4
    long p = mValues[idx][0];
    //  long phim = mValues[idx][1];
    long m = mValues[idx][2];
    // bits是？
    long bits = mValues[idx][4];
    long phi_m = mValues[idx][1];
    // p^d=1(mod m)
    long d = mValues[idx][3];
    long nslots = phi_m / d;
    if (d % 8 != 0)
    {
        throw std::logic_error("d is not a multiple of 8");
    }
    if (nslots < PlainByte)
    {
        throw std::logic_error("nslots is too small");
    }
    Context context(ContextBuilder<BGV>()
                        .m(m)
                        .p(p)
                        .r(1)
                        .c(c)
                        .bits(bits)
                        .build());

    // initialize context
    //----生成同态加密的公私钥，
    SecKey secretKey(context);
    // construct a secret key structure associated with the context
    PubKey &publicKey = secretKey;
    // an "upcast": SecKey is a subclass of PubKey
    secretKey.GenSecKey();
    addSome1DMatrices(secretKey);

    // 将公私钥写入到文件Client_publickey
    std::ofstream out1("Client_publickey", std::ios::binary);
    if (!out1.is_open())
    {
        std::cerr << "Failed to open " << "Client_publickey" << " for writing" << std::endl;
        return false;
    }
    publicKey.writeTo(out1);
    out1.close();
    // 将私钥写入到文件Client_secretkey
    std::ofstream out2("Client_secretkey", std::ios::binary);
    if (!out2.is_open())
    {
        std::cerr << "Failed to open " << "Client_secretkey" << " for writing" << std::endl;
        return false;
    }
    secretKey.writeTo(out2);
    out2.close();

    // 构造ea,ea.size=PlainByte=nslots,可以简单理解为打包环境设置
    static const uint8_t YuxPolyBytes[] = {0x1B, 0x1}; // X^8+X^4+X^3+X+1
    const GF2X YuxPoly = GF2XFromBytes(YuxPolyBytes, 2);
    EncryptedArrayDerived<PA_GF2> ea(context, YuxPoly, context.getAlMod());

    // 生成同态加密的初始对称密钥SymKey,打包状态
    std::cout << "SymKey encrypted start!" << std::endl;
    vector<Ctxt> encryptedSymKey;
    encryptSymKey(encryptedSymKey, SymKey, publicKey, ea);
    std::cout << "SymKey  encrypted cceeded!" << std::endl;

    // 将encryptedSymKey以字节形式写入到文件Client_encryptedSymKey.txt
    std::ofstream out3("Client_encryptedSymKey.txt");
    if (out3.is_open())
    {
        for (long i = 0; i < (long)encryptedSymKey.size(); i++)
        {
            std::ostringstream oss;
            encryptedSymKey[i].writeTo(oss);
            std::string str = oss.str();
            for (char c : str)
            {
                out3 << std::hex << std::setw(2) << std::setfill('0') << (int)(unsigned char)c << " ";
            }
        }
        out3.close();
    }
    else
    {
        std::cerr << "Failed to open Client_encryptedSymKey.txt for writing" << std::endl;
        return false;
    }
    return true;
}

bool Client_online()
{
    // 生成随机明文
    GF2X rnd;
    Vec<uint8_t> PlainStream(INIT_SIZE, PlainByte); // 8*10
    random(rnd, 8 * PlainStream.length());
    BytesFromGF2X(PlainStream.data(), rnd, PlainByte);

    // 将PlainStream以字节形式写入到文件Client_PlainStream.txt
    if (!writeToFile(PlainStream, "Client_PlainStream.txt", PlainByte))
    {
        return false;
    }

    // 从文件Client_KeyStream.txt读取KeyStream,注意判断长度是否等于PlainByte
    Vec<uint8_t> KeyStream(INIT_SIZE, PlainByte);
    if (!readFromFile(KeyStream, "Client_KeyStream.txt", PlainByte))
    {
        return false;
    }

    // 生成密文
    Vec<uint8_t> CipherStream(INIT_SIZE, PlainByte);
    for (unsigned i = 0; i < PlainByte; i++)
    {
        CipherStream[i] = PlainStream[i] ^ KeyStream[i];
    }
    // 检查加密解密是否正确
    for (unsigned i = 0; i < PlainByte; i++)
    {
        if (PlainStream[i] != (CipherStream[i] ^ KeyStream[i]))
        {
            std::cout << "Plain encryption failed!" << std::endl;
            return false;
        }
    }
    std::cout << "Plain encryption succeeded!" << std::endl;

    // 将CipherStream以字节形式写入到文件Client_CipherStream.txt
    if (!writeToFile(CipherStream, "Client_CipherStream.txt", PlainByte))
    {
        return false;
    }

    return true;
}