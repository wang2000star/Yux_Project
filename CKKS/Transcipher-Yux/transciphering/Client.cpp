#include "Client.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

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

    vector<vector<uint8_t>> Xset;
    Xset.resize(PlainBlock);
    vector<std::array<std::array<unsigned char, BlockByte>, Nr + 1>> Xset;
    std::vector<std::array<std::array<unsigned char, BlockByte>, Nr + 1>> RoundKeySet;
    std::array<uint64_t, PlainBlock> NonceSet; // Nonce集合
    std::array<unsigned char, BlockByte * PlainBlock> KeyStream;

    // 创建RandomBit对象
    RandomBit<BlockSize> randomBit(Nr);
    // 为roundconstants创建别名
    auto &RanVecs = randomBit.roundconstants;

    // 定义counter_begin和counter_end
    uint64_t counter_begin = 0;
    uint64_t counter_end = PlainBlock - 1;

    for (uint64_t counter = counter_begin; counter <= counter_end; counter++)
    {
        uint64_t nonce = generate_secure_random_int(NonceSize);
        NonceSet[counter - counter_begin] = nonce;
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

        // 生成轮密钥
        std::array<std::array<unsigned char, BlockByte>, Nr + 1> RoundKey;
        for (unsigned r = 0; r <= Nr; r++)
        {
            for (unsigned i = 0; i < BlockByte; i++)
            {
                RoundKey[r][i] = mul(SymKey[i], X[r][i]);
            }
        }
        // 把RoundKey存入RoundKey_set
        RoundKeySet.push_back(RoundKey);

        // the first round
        std::array<unsigned char, BlockByte> state;

        for (unsigned i = 0; i < BlockByte; i++)
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
            for (unsigned i = 0; i < BlockByte; i++)
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
        for (unsigned i = 0; i < BlockByte; i++)
        {
            state[i] ^= RoundKey[Nr][i];
        }
        for (unsigned i = 0; i < BlockByte; i++)
        {
            KeyStream[(counter - counter_begin) * BlockByte + i] = state[i];
        }
    }

    // 将KeyStream以字节形式写入到文件Client_KeyStream.txt
    std::ofstream out("Client_KeyStream.txt");
    if (out.is_open())
    {
        for (unsigned i = 0; i < PlainByte; i++)
        {
            out << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(KeyStream[i]) << " ";
        }
        out.close();
    }
    else
    {
        std::cerr << "Failed to open Client_KeyStream.txt for writing" << std::endl;
        return false;
    }

    // 将Nonce以字节形式写入到文件Client_NonceSet.txt
    std::ofstream out2("Client_NonceSet.txt");
    if (out2.is_open())
    {
        for (unsigned i = 0; i < PlainBlock; i++)
        {
            out2 << std::hex << std::setw(2) << std::setfill('0') << NonceSet[i] << " ";
        }
        out2.close();
    }
    else
    {
        std::cerr << "Failed to open Client_NonceSet.txt for writing" << std::endl;
        return false;
    }
    // 将X_set写入到文件Client_Xset.txt
    std::ofstream out3("Client_Xset.txt");
    if (out3.is_open())
    {
        for (const auto &x : Xset)
        {
            for (const auto &round : x)
            {
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    out3 << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(round[i]) << " ";
                }
                out3 << std::endl;
            }
        }
        out3.close();
    }
    else
    {
        std::cerr << "Failed to open Client_Xset.txt for writing" << std::endl;
        return false;
    }

    // 将RoundKey_set写入到文件Client_RoundKeySet.txt
    std::ofstream out4("Client_RoundKeySet.txt");
    if (out4.is_open())
    {
        for (const auto &roundKey : RoundKeySet)
        {
            for (const auto &round : roundKey)
            {
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    out4 << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(round[i]) << " ";
                }
                out4 << std::endl;
            }
        }
        out4.close();
    }
    else
    {
        std::cerr << "Failed to open Client_RoundKeySet.txt for writing" << std::endl;
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

    long bits = mValues[idx][4];
    long phi_m = mValues[idx][1];
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
    std::ofstream out5("Client_publickey", std::ios::binary);
    if (!out5.is_open())
    {
        std::cerr << "Failed to open " << "Client_publickey" << " for writing" << std::endl;
        return false;
    }
    publicKey.writeTo(out5);
    out5.close();
    // 将私钥写入到文件Client_secretkey
    std::ofstream out6("Client_secretkey", std::ios::binary);
    if (!out6.is_open())
    {
        std::cerr << "Failed to open " << "Client_secretkey" << " for writing" << std::endl;
        return false;
    }
    secretKey.writeTo(out6);
    out6.close();

    static const uint8_t YuxPolyBytes[] = {0x1B, 0x1}; // X^8+X^4+X^3+X+1
    const GF2X YuxPoly = GF2XFromBytes(YuxPolyBytes, 2);
    EncryptedArrayDerived<PA_GF2> ea(context, YuxPoly, context.getAlMod());
    // 生成同态加密的初始对称密钥SymKey,打包状态
    vector<Ctxt> encryptedSymKey;
    encryptSymKey(encryptedSymKey, SymKey, publicKey, ea);
    std::cout << "SymKey FHE succeeded!" << std::endl;

    // 将encryptedSymKey以字节形式写入到文件Client_encryptedSymKey.txt
    std::ofstream out7("Client_encryptedSymKey.txt");
    if (out7.is_open())
    {
        for (long i = 0; i < (long)encryptedSymKey.size(); i++)
        {
            std::ostringstream oss;
            encryptedSymKey[i].writeTo(oss);
            std::string str = oss.str();
            for (char c : str)
            {
                out7 << std::hex << std::setw(2) << std::setfill('0') << (int)(unsigned char)c << " ";
            }
        }
        out7.close();
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
    std::ofstream out1("Client_PlainStream.txt");
    if (out1.is_open())
    {
        for (unsigned i = 0; i < PlainByte; i++)
        {
            out1 << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(PlainStream[i]) << " ";
        }
        out1.close();
    }
    else
    {
        std::cerr << "Failed to open Client_PlainStream.txt for writing" << std::endl;
        return false;
    }
    // 从文件KeyStream.txt读取KeyStream
    // 从文件Client_KeyStream.txt读取KeyStream
    std::ifstream in("Client_KeyStream.txt");
    unsigned char KeyStream[PlainByte];
    if (in.is_open())
    {
        unsigned int byte;
        for (unsigned i = 0; i < PlainByte; i++)
        {
            in >> std::hex >> byte;
            KeyStream[i] = static_cast<unsigned char>(byte);
        }
        in.close();
    }
    else
    {
        std::cerr << "Failed to open Client_KeyStream.txt for reading" << std::endl;
        return false;
    }
    Vec<uint8_t> CipherStream(INIT_SIZE, PlainByte);
    // 密文
    for (unsigned i = 0; i < PlainByte; i++)
    {
        CipherStream[i] = PlainStream[i] ^ KeyStream[i];
    }
    // 检查解密是否正确
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
    std::ofstream out7("Client_CipherStream.txt");
    if (out7.is_open())
    {
        for (unsigned i = 0; i < PlainByte; i++)
        {
            out7 << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(CipherStream[i]) << " ";
        }
        out7.close();
    }
    else
    {
        std::cerr << "Failed to open Client_CipherStream.txt for writing" << std::endl;
        return false;
    }

    return true;
}