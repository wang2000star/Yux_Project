#include "Client.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

// run the Yux key-expansion and then encrypt the expanded key.
void encryptSymKey(vector<Ctxt> &encryptedSymKey, Vec<uint8_t> &roundKeySchedule, const PubKey &hePK,
                   const EncryptedArrayDerived<PA_GF2> &ea)
{
    auto start_encryptSymKey = std::chrono::steady_clock::now();
    
    long nRoundKeys = 1;
    long round_key_length = BlockByte;
    long length_s = round_key_length * nRoundKeys;
    long blocksPerCtxt = ea.size() / BlockByte;

    Vec<uint8_t> expanded(INIT_SIZE, nRoundKeys * blocksPerCtxt * BlockByte);
    for (long i = 0; i < nRoundKeys; i++)
    {
        uint8_t *roundKey = &roundKeySchedule[16 * i];
        for (long j = 0; j < blocksPerCtxt; j++)
            memcpy(&expanded[16 * (i * blocksPerCtxt + j)], roundKey, 16);
    }

    Vec<ZZX> encoded;
    encodeTo1Ctxt(encoded, expanded, ea);

    {
        Ctxt tmpCtxt(hePK);
        encryptedSymKey.resize(encoded.length(), tmpCtxt);
    }

    for (long i = 0; i < (long)encryptedSymKey.size(); i++)
        hePK.Encrypt(encryptedSymKey[i], encoded[i]);

    auto end_encryptSymKey = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_encryptSymKey = end_encryptSymKey - start_encryptSymKey;
    std::cout << "SymKey encryption time: " << elapsed_seconds_encryptSymKey.count() << "s\n";
}

bool writeEncryptedSymKey(const vector<Ctxt>& encryptedSymKey, const std::string& filename)
{
    auto start_writeSymKey = std::chrono::steady_clock::now();
    
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open())
    {
        std::cerr << "Failed to open " << filename << " for writing" << std::endl;
        return false;
    }

    for (const auto& ctxt : encryptedSymKey)
    {
        ctxt.writeTo(out);
    }

    out.close();

    auto end_writeSymKey = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_writeSymKey = end_writeSymKey - start_writeSymKey;
    std::cout << "Write encrypted SymKey to file time: " << elapsed_seconds_writeSymKey.count() << "s\n";

    return true;
}
bool verifySymKeyDecryption(const vector<Ctxt>& encryptedSymKey, const SecKey& secretKey,
                      const EncryptedArrayDerived<PA_GF2>& ea, const Vec<uint8_t>& originalSymKey)
{
    auto start_decrypt = std::chrono::steady_clock::now();

    // 解密密钥
    Vec<ZZX> decryptedPolys(INIT_SIZE, encryptedSymKey.size());
    for (long i = 0; i < encryptedSymKey.size(); ++i)
    {
        secretKey.Decrypt(decryptedPolys[i], encryptedSymKey[i]);
    }

    // 解码解密后的多项式到明文字节
    Vec<uint8_t> decryptedSymKey(INIT_SIZE, originalSymKey.length());
    decodeTo1Ctxt(decryptedSymKey, decryptedPolys, ea);

    // 验证解密结果是否与原始对称密钥匹配
    bool isDecryptionCorrect = true;
    for (long i = 0; i < originalSymKey.length(); ++i)
    {
        if (decryptedSymKey[i] != originalSymKey[i])
        {
            std::cout << "Decryption check failed at index " << i << ": expected " << (int)originalSymKey[i]
                      << ", got " << (int)decryptedSymKey[i] << std::endl;
            isDecryptionCorrect = false;
        }
    }

    if (isDecryptionCorrect)
    {
        std::cout << "Decryption check succeeded: Decrypted symmetric key matches original key." << std::endl;
    }
    else
    {
        std::cout << "Decryption check failed: Decrypted symmetric key does not match original key." << std::endl;
    }

    auto end_decrypt = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_decrypt = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification time: " << elapsed_seconds_decrypt.count() << "s\n";

    return isDecryptionCorrect;
}

bool Client_offline()
{
    auto start = std::chrono::steady_clock::now();

    // Generating symmetric key and key stream
    auto start_keyStream = std::chrono::steady_clock::now();

    vector<uint8_t> iv(BlockByte);
    for (unsigned i = 0; i < BlockByte; i++) {
        iv[i] = i + 1;
    }

    GF2X rnd;
    Vec<uint8_t> SymKey(INIT_SIZE, BlockByte);
    random(rnd, 8 * SymKey.length());
    BytesFromGF2X(SymKey.data(), rnd, SymKey.length());

    Vec<uint64_t> NonceSet(INIT_SIZE, PlainBlock);
    Vec<uint8_t> Xset(INIT_SIZE, PlainByte * (Nr + 1));
    Vec<uint8_t> RoundKeySet(INIT_SIZE, PlainByte * (Nr + 1));
    Vec<uint8_t> KeyStream(INIT_SIZE, PlainByte);

    RandomBit<BlockSize> randomBit(Nr);
    auto &RanVecs = randomBit.roundconstants;

    uint64_t counter_begin = 0;
    uint64_t counter_end = PlainBlock + counter_begin - 1;

    for (uint64_t counter = counter_begin; counter <= counter_end; counter++) {
        uint64_t nonce = generate_secure_random_int(NonceSize);
        NonceSet[counter - counter_begin] = nonce;
        randomBit.generate_Instance_all_new(nonce, counter);

        Vec<uint8_t> state(INIT_SIZE, BlockByte);
        for (unsigned r = 0; r <= Nr; r++) {
            Vec<uint8_t> X(INIT_SIZE, BlockByte);
            Vec<uint8_t> RoundKey(INIT_SIZE, BlockByte);
            uint64_t temp;
            for (unsigned i = 0; i < BlockByte; ++i) {
                bool bit_array[8];
                for (unsigned j = 0; j < 8; ++j) {
                    bit_array[j] = RanVecs[r][i * 8 + j];
                }
                BinStrToHex(bit_array, temp, 8);
                X[i] = static_cast<uint8_t>(temp);
                //RoundKey[i] = SymKey[i] ^ X[i];
                RoundKey[i] = mul(SymKey[i], X[i]);
            }
            memcpy(&Xset[PlainByte * r + BlockByte * (counter - counter_begin)], X.data(), BlockByte);
            memcpy(&RoundKeySet[PlainByte * r + BlockByte * (counter - counter_begin)], RoundKey.data(), BlockByte);
            if (r == 0) {
                for (unsigned i = 0; i < BlockByte; i++) {
                    state[i] = RoundKey[i] ^ iv[i];
                }
            } else if (r < Nr) {
                for (unsigned i = 0; i < 4; i++) {
                    decSboxFi(state.data(), i * 4);
                    decSboxFi(state.data(), i * 4);
                    decSboxFi(state.data(), i * 4);
                    decSboxFi(state.data(), i * 4);
                }
                decLinearLayer(state.data());
                for (unsigned i = 0; i < BlockByte; i++) {
                    state[i] ^= RoundKey[i];
                }
            } else {
                for (unsigned i = 0; i < 4; i++) {
                    decSboxFi(state.data(), i * 4);
                    decSboxFi(state.data(), i * 4);
                    decSboxFi(state.data(), i * 4);
                    decSboxFi(state.data(), i * 4);
                }
                for (unsigned i = 0; i < BlockByte; i++) {
                    state[i] ^= RoundKey[i];
                }
                memcpy(&KeyStream[(counter - counter_begin) * BlockByte], state.data(), BlockByte);
            }
        }
    }
    auto end_keyStream = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_keyStream = end_keyStream - start_keyStream;
    std::cout << "Symmetric Key and KeyStream Generation time: " << elapsed_seconds_keyStream.count() << "s\n";

    // Write KeyStream to files
    auto start_fileWrite = std::chrono::steady_clock::now();

    if (!writeToFile(NonceSet.data(), "Client_NonceSet.txt", PlainBlock) ||
        !writeToFile(Xset.data(), "Client_Xset.txt", PlainByte * (Nr + 1)) ||
        !writeToFile(RoundKeySet.data(), "Client_RoundKeySet.txt", PlainByte * (Nr + 1)) ||
        !writeToFile(KeyStream.data(), "Client_KeyStream.txt", PlainByte))
    {
        return false;
    }

    auto end_fileWrite = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_fileWrite = end_fileWrite - start_fileWrite;
    std::cout << "Write KeyStream to file time: " << elapsed_seconds_fileWrite.count() << "s\n";

    // Encrypt the symmetric key
    auto start_keyEncryption = std::chrono::steady_clock::now();

    long idx = 3;
    long c = 9;
    bool packed = true;
    if (idx > 5) idx = 5;
    long p = mValues[idx][0];
    long m = mValues[idx][2];
    long bits = mValues[idx][4];
    long phi_m = mValues[idx][1];
    long d = mValues[idx][3];
    long nslots = phi_m / d;
    if (d % 8 != 0 || nslots < PlainByte) {
        throw std::logic_error("Invalid parameters for encryption.");
    }

    Context context(ContextBuilder<BGV>().m(m).p(p).r(1).c(c).bits(bits).build());
    SecKey secretKey(context);
    PubKey &publicKey = secretKey;
    secretKey.GenSecKey();
    addSome1DMatrices(secretKey);

    std::ofstream outContext("Client_context", std::ios::binary);
    if (!outContext.is_open()) {
        std::cerr << "Failed to open Client_context for writing" << std::endl;
        return false;
    }
    context.writeTo(outContext);
    outContext.close();

    std::ofstream outPublicKey("Client_publickey", std::ios::binary);
    if (!outPublicKey.is_open()) {
        std::cerr << "Failed to open Client_publickey for writing" << std::endl;
        return false;
    }
    publicKey.writeTo(outPublicKey);
    outPublicKey.close();

    std::ofstream outSecretKey("Client_secretkey", std::ios::binary);
    if (!outSecretKey.is_open()) {
        std::cerr << "Failed to open Client_secretkey for writing" << std::endl;
        return false;
    }
    secretKey.writeTo(outSecretKey);
    outSecretKey.close();

    static const uint8_t YuxPolyBytes[] = {0x1B, 0x1};
    const GF2X YuxPoly = GF2XFromBytes(YuxPolyBytes, 2);
    EncryptedArrayDerived<PA_GF2> ea(context, YuxPoly, context.getAlMod());

    vector<Ctxt> encryptedSymKey;
    encryptSymKey(encryptedSymKey, SymKey, publicKey, ea);
    // 解密验证
    // Verify the decryption of the encrypted symmetric key
    std::cout << "Verifying decryption of encrypted symmetric key..." << std::endl;
   
    Vec<uint8_t> expanded(INIT_SIZE, PlainByte);
    for (long i = 0; i < 1; i++)
    {
        for (long j = 0; j <  PlainBlock; j++)
            memcpy(&expanded[16 * (i * PlainBlock + j)], &SymKey[BlockByte * i], BlockByte);
    }
    bool decryptionCorrect = verifySymKeyDecryption(encryptedSymKey, secretKey, ea, expanded);
    if (!decryptionCorrect)
    {
        std::cerr << "Decryption verification failed for SymKey." << std::endl;
        return false;
    }

    if (!writeEncryptedSymKey(encryptedSymKey, "Client_encryptedSymKey.bin"))
    {
        return false;
    }

    auto end_keyEncryption = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_keyEncryption = end_keyEncryption - start_keyEncryption;
    std::cout << "Symmetric key encryption and write to file time: " << elapsed_seconds_keyEncryption.count() << "s\n";

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Client offline total time: " << elapsed_seconds.count() << "s\n";

    return true;
}

bool Client_online()
{
    auto start = std::chrono::steady_clock::now();

    GF2X rnd;
    Vec<uint8_t> PlainStream(INIT_SIZE, PlainByte);
    random(rnd, 8 * PlainStream.length());
    BytesFromGF2X(PlainStream.data(), rnd, PlainByte);

    if (!writeToFile<uint8_t>(PlainStream.data(), "Client_PlainStream.txt", PlainByte))
    {
        return false;
    }

    Vec<uint8_t> KeyStream(INIT_SIZE, PlainByte);
    if (!readFromFile<uint8_t>(KeyStream.data(), "Client_KeyStream.txt", PlainByte))
    {
        return false;
    }

    Vec<uint8_t> CipherStream(INIT_SIZE, PlainByte);
    for (unsigned i = 0; i < PlainByte; i++) {
        CipherStream[i] = PlainStream[i] ^ KeyStream[i];
    }

    for (unsigned i = 0; i < PlainByte; i++) {
        if (PlainStream[i] != (CipherStream[i] ^ KeyStream[i])) {
            std::cout << "Plain encryption failed!" << std::endl;
            return false;
        }
    }
    std::cout << "Plain encryption succeeded!" << std::endl;

    if (!writeToFile(CipherStream.data(), "Client_CipherStream.txt", PlainByte)) {
        return false;
    }

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Client online total time: " << elapsed_seconds.count() << "s\n";

    return true;
}