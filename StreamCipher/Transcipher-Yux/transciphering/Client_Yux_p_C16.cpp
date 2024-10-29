#include "Client_Yux_p_C16.hpp"

using namespace std;
using namespace helib;
using namespace NTL;




void encryptSymKey(vector<Ctxt>& eKey, vector<uint64_t>& Key, unique_ptr<PubKey>& he_pk, EncryptedArray& ea)
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

bool writeEncryptedSymKey(const vector<Ctxt> &encryptedSymKey, const std::string &filename)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open())
    {
        std::cerr << "Failed to open " << filename << " for writing" << std::endl;
        return false;
    }

    for (const auto &ctxt : encryptedSymKey)
    {
        ctxt.writeTo(out);
    }

    out.close();

    return true;
}

namespace C1ient_Yux_p_C16{
bool Client_offline()
{
    printf("Client offline 16 start!\n");
    auto start = std::chrono::steady_clock::now();

    // Generating symmetric key and key stream
    auto start_keyStream = std::chrono::steady_clock::now();

    vector<uint8_t> iv(BlockByte);
    for (unsigned i = 0; i < BlockByte; i++)
    {
        iv[i] = i + 1;
    }

    GF2X rnd;
    Vec<uint8_t> SymKey0(INIT_SIZE, 2*BlockByte);
    random(rnd, 8 * SymKey0.length());
    BytesFromGF2X(SymKey0.data(), rnd, SymKey0.length());

    vector<uint64_t> SymKey(BlockByte);
    // SymKey0的两个字节合并成一个16位的字，即一个uint16_t，存入SymKey
    for (unsigned i = 0; i < BlockByte; i++)
    {
        SymKey[i] = (SymKey0[2 * i] << 8) | SymKey0[2 * i + 1];
    }


    Vec<uint64_t> NonceSet(INIT_SIZE, PlainBlock);
    Vec<uint64_t> Xset(INIT_SIZE, PlainByte * (Nr + 1));
    Vec<uint64_t> RoundKeySet(INIT_SIZE, PlainByte * (Nr + 1));
    Vec<uint64_t> KeyStream(INIT_SIZE, PlainByte);

    RandomBit<BlockSize> randomBit(Nr);
    auto &RanVecs = randomBit.roundconstants;

    uint64_t counter_begin = 0;
    uint64_t counter_end = PlainBlock + counter_begin - 1;

    for (uint64_t counter = counter_begin; counter <= counter_end; counter++)
    {
        uint64_t nonce = generate_secure_random_int(NonceSize);
        NonceSet[counter - counter_begin] = nonce;
        randomBit.generate_Instance_all_new(nonce, counter);

        Vec<uint64_t> state(INIT_SIZE, BlockByte);
        for (unsigned r = 0; r <= Nr; r++)
        {
            Vec<uint64_t> X(INIT_SIZE, BlockByte);
            Vec<uint64_t> RoundKey(INIT_SIZE, BlockByte);
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
                // RoundKey[i] = (SymKey[i]+X[i])%PlainMod;
                RoundKey[i] = (SymKey[i]*X[i])%PlainMod;
            }
            memcpy(&Xset[PlainByte * r + BlockByte * (counter - counter_begin)], X.data(), BlockByte);
            memcpy(&RoundKeySet[PlainByte * r + BlockByte * (counter - counter_begin)], RoundKey.data(), BlockByte);
            if (r == 0)
            {
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    // state[i] = (RoundKey[i] + iv[i]) % PlainMod; 
                    state[i] = (RoundKey[i] * iv[i])%PlainMod;
                }
            }
            else if (r < Nr)
            {
                for (unsigned i = 0; i < 4; i++)
                {
                    decSboxFi(state.data(), i * 4);
                    decSboxFi(state.data(), i * 4);
                    decSboxFi(state.data(), i * 4);
                    decSboxFi(state.data(), i * 4);
                }
                decLinearLayer(state.data());
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    // state[i] = (state[i] + RoundKey[i]) % PlainMod;
                    state[i] = (state[i]+RoundKey[i])%PlainMod;
                }
            }
            else
            {
                for (unsigned i = 0; i < 4; i++)
                {
                    decSboxFi(state.data(), i * 4);
                    decSboxFi(state.data(), i * 4);
                    decSboxFi(state.data(), i * 4);
                    decSboxFi(state.data(), i * 4);
                }
                for (unsigned i = 0; i < BlockByte; i++)
                {
                    //state[i] = (state[i] + RoundKey[i]) % PlainMod;
                    state[i] = (state[i] * RoundKey[i]) % PlainMod;
                }
                memcpy(&KeyStream[(counter - counter_begin) * BlockByte], state.data(), BlockByte);
            }
        }
    }
    auto end_keyStream = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_keyStream = end_keyStream - start_keyStream;
    std::cout << "KeyStream Generation time: " << elapsed_seconds_keyStream.count() << "s\n";

    if (!writeToFile(NonceSet.data(), "Client_NonceSet.txt", PlainBlock) ||
        !writeToFile(Xset.data(), "Client_Xset.txt", PlainByte * (Nr + 1)) ||
        !writeToFile(RoundKeySet.data(), "Client_RoundKeySet.txt", PlainByte * (Nr + 1)) ||
        !writeToFile(KeyStream.data(), "Client_KeyStream.txt", PlainByte))
    {
        return false;
    }

    // generate FHE sk and pk；use it Encrypt the symmetric key
    auto start_keyEncryption = std::chrono::steady_clock::now();

    int idx = 1;
    if (idx >1) {
      idx = 1;
    }

    long m = mValues[idx][1];
    long p = mValues[idx][0];
    long r = 1;
    long L = mValues[idx][2];
    long c = mValues[idx][3];
    long d = 1;
    long k = 128;
    long s = 1;

    if (!m) m = FindM(k, L, c, p, d, s, 0);

    shared_ptr<Context> context(ContextBuilder<BGV>()
                                             .m(m)
                                             .p(p)
                                             .r(r)
                                             .bits(L)
                                             .c(c)
                                             .buildPtr());
    SecKey secretKey(*context);
    secretKey.GenSecKey();
    unique_ptr<PubKey> publicKey = std::make_unique<helib::PubKey>(secretKey);
    helib::EncryptedArray ea(context->getEA());
    long nslots = ea.size();
    printf("nslots = %ld\n", nslots);
    

    std::ofstream outContext("Client_context", std::ios::binary);
    if (!outContext.is_open())
    {
        std::cerr << "Failed to open Client_context for writing" << std::endl;
        return false;
    }
    context->writeTo(outContext);
    outContext.close();

    std::ofstream outPublicKey("Client_publickey", std::ios::binary);
    if (!outPublicKey.is_open())
    {
        std::cerr << "Failed to open Client_publickey for writing" << std::endl;
        return false;
    }
    publicKey->writeTo(outPublicKey);
    outPublicKey.close();

    std::ofstream outSecretKey("Client_secretkey", std::ios::binary);
    if (!outSecretKey.is_open())
    {
        std::cerr << "Failed to open Client_secretkey for writing" << std::endl;
        return false;
    }
    secretKey.writeTo(outSecretKey);
    outSecretKey.close();

    vector<Ctxt> encryptedSymKey;
    encryptSymKey(encryptedSymKey, SymKey, publicKey, ea);

    auto end_keyEncryption = std::chrono::steady_clock::now();
    double keyEncryption = std::chrono::duration<double>(end_keyEncryption - start_keyEncryption).count();
    std::cout << "PublicKey generation and SymKey FHE time: " << keyEncryption << "s\n";


    // 解密验证
    // Verify the decryption of the encrypted symmetric key
    std::cout << "Verifying decryption of encrypted symmetric key..." << std::endl;

    vector<uint64_t> expanded(PlainByte*BlockByte);
    // 分别把每个SymKey的元素赋值PlainByte次
    for (unsigned i = 0; i < BlockByte; i++)
    {
        for (unsigned j = 0; j < PlainByte; j++)
        {
            expanded[i * PlainByte + j] = SymKey[i];
        }
    }
    // 解密验证
    bool decryptionCorrect = verifyDecryption_p16(encryptedSymKey, expanded, secretKey, ea);
    if (!decryptionCorrect)
    {
        std::cerr << "Decryption verification failed for SymKey." << std::endl;
        return false;
    }
    if (!writeEncryptedSymKey(encryptedSymKey, "Client_encryptedSymKey.bin"))
    {
        return false;
    }

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Client offline total time: " << elapsed_seconds.count() << "s\n";

    // 离线客户端时间=KeyStream Generation time+PublicKey generation and SymKey FHE time
    double total_time_off = elapsed_seconds_keyStream.count() + keyEncryption;
    std::cout << "Encryption offline total time: " << total_time_off << "s\n";

    return true;
}

bool Client_online()
{
    auto start = std::chrono::steady_clock::now();

    GF2X rnd;
    Vec<uint8_t> PlainStream0(INIT_SIZE, 2*BlockByte*PlainBlock);
    random(rnd, 8 * PlainStream0.length());
    BytesFromGF2X(PlainStream0.data(), rnd, PlainStream0.length());

    vector<uint64_t> PlainStream(BlockByte*PlainBlock);
    for (unsigned i = 0; i < PlainBlock; i++)
    {
        for (unsigned j = 0; j < BlockByte; j++)
        {
            PlainStream[i * BlockByte + j] = (PlainStream0[2 * (i * BlockByte + j)] << 8) | PlainStream0[2 * (i * BlockByte + j) + 1];
        }
    }

    if (!writeToFile<uint64_t>(PlainStream.data(), "Client_PlainStream.txt", PlainByte))
    {
        return false;
    }

    vector<uint64_t> KeyStream(PlainByte);
    if (!readFromFile<uint64_t>(KeyStream.data(), "Client_KeyStream.txt", PlainByte))
    {
        return false;
    }

    vector<uint64_t> CipherStream(PlainByte);
    auto start_encryption = std::chrono::steady_clock::now();
    for (unsigned i = 0; i < PlainByte; i++)
    {
        CipherStream[i] = PlainStream[i] ^ KeyStream[i];
    }
    auto end_encryption = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_encryption = end_encryption - start_encryption;

    for (unsigned i = 0; i < PlainByte; i++)
    {
        if (PlainStream[i] != (CipherStream[i] ^ KeyStream[i]))
        {
            std::cout << "Plain encryption failed!" << std::endl;
            return false;
        }
    }
    std::cout << "Plain encryption succeeded!" << std::endl;

    if (!writeToFile(CipherStream.data(), "Client_CipherStream.txt", PlainByte))
    {
        return false;
    }

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Client online total time: " << elapsed_seconds.count() << "s\n";

    std::cout << "Encryption online time: " << elapsed_seconds_encryption.count() << "s\n";
    return true;
}
} // namespace C1ient_Yux_p_C16