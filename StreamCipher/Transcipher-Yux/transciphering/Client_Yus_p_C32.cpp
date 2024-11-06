#include "Client_Yus_p_C32.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

void encryptSymKey(vector<Ctxt> &encryptedSymKey, vector<long> &SymKey, unique_ptr<PubKey> &pk, EncryptedArray &ea)
{
    long nslots = ea.size();
    // 加密
    encryptedSymKey.resize(32, Ctxt(*pk));
    for (long i = 0; i < 32; i++)
    { // encrypt the encoded key
        vector<long> slotsData(nslots, SymKey[i]);
        ea.encrypt(encryptedSymKey[i], *pk, slotsData);
    }
}

bool verify_encryptSymKey(vector<Ctxt> &encryptedSymKey, vector<long> &SymKey, const SecKey &secretKey, EncryptedArray &ea)
{
    vector<long> decryptedSymKey(32);
    for (long i = 0; i < 32; i++)
    { // encrypt the encoded key
        vector<long> slotsData;
        ea.decrypt(encryptedSymKey[i], secretKey, slotsData);
        decryptedSymKey[i] = slotsData[0];
    }
    bool isDecryptedSymKeyCorrect = true;
    for (long i = 0; i < 32; i++)
    {
        if (decryptedSymKey[i] != SymKey[i])
        {
            std::cout << "Decryption check failed at index " << i << ": expected " << SymKey[i]
                      << ", got " << decryptedSymKey[i] << std::endl;
            isDecryptedSymKeyCorrect = false;
            break;
        }
    }
    if (isDecryptedSymKeyCorrect)
    {
        std::cout << "Decryption check succeeded: Decrypted vector matches original vector." << std::endl;
    }
    else
    {
        std::cout << "Decryption check failed: Decrypted vector does not match original vector." << std::endl;
    }
    return isDecryptedSymKeyCorrect;
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

namespace C1ient_Yus_p_C32
{
    bool Client_offline()
    {
        // 定义初始向量
        vector<long> iv(BlockByte);
        for (unsigned i = 0; i < BlockByte; i++)
        {
            iv[i] = i + 1;
        }
        // 生成随机对称密钥
        GF2X rnd;
        int Bytebitsdiv8 = ceil(Bytebits / 8);
        vector<uint8_t> SymKey0(Bytebitsdiv8 * BlockByte);
        random(rnd, 8 * SymKey0.size());
        BytesFromGF2X(SymKey0.data(), rnd, SymKey0.size());
        vector<long> SymKey(BlockByte);
        for (unsigned i = 0; i < BlockByte; i++)
        {
            SymKey[i] = 0;
            for (unsigned j = 0; j < Bytebitsdiv8; j++)
            {
                SymKey[i] += (SymKey0[Bytebitsdiv8 * i + j] << (8 * j));
            }
            SymKey[i] %= PlainMod;
        }

        std::cout << "SymKey generated." << std::endl;
        // 将对称密钥写入文件
        if (!writeToFile(SymKey.data(), "Client_SymKey.txt", BlockByte))
        {
            return false;
        }
        std::cout << "Client_SymKey has been written to file." << std::endl;
        //========
        // Generating symmetric key and key stream
        auto start_keyStream = std::chrono::steady_clock::now();

        std::vector<long> NonceSet(PlainBlock);
        std::vector<long> Xset(PlainByte * (Nr + 1));
        std::vector<long> RoundKeySet(PlainByte * (Nr + 1));
        std::vector<long> KeyStream(PlainByte);

        RandomBit<BlockSize> randomBit(Nr+1);

        long total_tasks = counter_end - counter_begin + 1;
        std::atomic<long> completed_tasks(0);

        // 定义进度打印的粒度，例如每完成 1% 的任务打印一次
        long progress_step = total_tasks / 100;
        if (progress_step == 0)
            progress_step = 1; // 防止除零

#pragma omp parallel for firstprivate(randomBit)
        for (long counter = counter_begin; counter <= counter_end; counter++)
        {
            // long nonce = generate_secure_random_int(NonceSize);          
            // randomBit.generate_Instance_all_new(nonce, counter);
            // auto &RanVecs = randomBit.roundconstants;
            long nonce = counter; 
            std::vector<std::bitset<544>> RanVecs(Nr + 1);

            NonceSet[counter - counter_begin] = nonce;
            // 使用 std::array 代替 vector 并固定大小
            std::vector<long> state; // 初始化 state
            std::vector<long> X;
            std::vector<long> RoundKey;
            bool bit_array[Bytebits];
            uint64_t temp;

            // 逐轮进行加密
            for (unsigned r = 0; r <= Nr; r++)
            {
                // 计算 X 和 RoundKey
                for (unsigned i = 0; i < BlockByte; ++i)
                {
                    for (unsigned j = 0; j < Bytebits; ++j)
                    {
                        bit_array[j] = RanVecs[r][i * Bytebits + j];
                    }
                    BinStrToHex(bit_array, temp, Bytebits);
                    X[i] = temp % PlainMod;
                    RoundKey[i] = (SymKey[i] * X[i]) % PlainMod;
                }

                // 将 X 和 RoundKey 复制到 Xset 和 RoundKeySet
                memcpy(&Xset[PlainByte * r + BlockByte * (counter - counter_begin)], X.data(), BlockByte * sizeof(long));
                memcpy(&RoundKeySet[PlainByte * r + BlockByte * (counter - counter_begin)], RoundKey.data(), BlockByte * sizeof(long));

                if (r == 0)
                { // 初始轮
                    for (unsigned i = 0; i < BlockByte; i++)
                    {
                        state[i] = (RoundKey[i] + iv[i]) % PlainMod;
                    }
                }
                else if (r < Nr)
                {                // 常规轮
                    MC(state);   // 行移位
                    MR(state);   // 列混淆
                    Sbox(state); // S盒
                    for (unsigned i = 0; i < BlockByte; i++)
                    {
                        state[i] = (state[i] + RoundKey[i]) % PlainMod;
                    }
                }
                else
                {                // 最后一轮
                    MC(state);   // 行移位
                    MR(state);   // 列混淆
                    Sbox(state); // S盒
                    MC(state);   // 再次行移位
                    MR(state);   // 再次列混淆
                    for (unsigned i = 0; i < BlockByte; i++)
                    {
                        state[i] = (state[i] + RoundKey[i]) % PlainMod;
                    }
                    memcpy(&KeyStream[(counter - counter_begin) * BlockByte], state.data(), BlockByte * sizeof(long));
                }
            }

            // 更新已完成的任务数
            long local_completed = ++completed_tasks;

            // 定期打印进度
            if (local_completed % progress_step == 0 || local_completed == total_tasks)
            {
#pragma omp critical
                {
                    std::cout << "Progress: " << (local_completed * 100) / total_tasks << "% completed.\r" << std::flush;
                }
            }
        }
        std::cout << std::endl; // 完成后换行

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
        std::cout << "Client_NonceSet, Client_Xset, Client_RoundKeySet, Client_KeyStream have been written to files." << std::endl;
        // Generating Public Key and encrypting the symmetric key
        auto start_PubKey = std::chrono::steady_clock::now();

        int idx = 1;
        long m = mValues[idx][1]; // m=65536,phi(m)=32768=2^15
        long p = mValues[idx][0];
        long r = 1;
        long L = mValues[idx][2];
        long c = mValues[idx][3];
        long d = 1; // slots = phi(m)/d = phi(m) = 32768 = PlainBlock
        long k = 128;
        long s = 1;

        if (!m)
            m = FindM(k, L, c, p, d, s, 0);
        printf("Selected m = %ld\n", m);

        shared_ptr<Context> context(ContextBuilder<BGV>()
                                        .m(m)
                                        .p(p)
                                        .r(r)
                                        .bits(L)
                                        .c(c)
                                        .buildPtr());

        printf("Security");
        SecKey secretKey(*context);
        secretKey.GenSecKey();
        unique_ptr<PubKey> publicKey = std::make_unique<helib::PubKey>(secretKey);
        helib::EncryptedArray ea(context->getEA());
        long nslots = ea.size();
        printf("nslots = %ld\n", nslots);
        if (nslots != PlainBlock)
        {
            std::cerr << "nslots != PlainBlock" << std::endl;
            return false;
        }
        auto end_PubKey = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds_PubKey = end_PubKey - start_PubKey;
        std::cout << "PublicKey generation time: " << elapsed_seconds_PubKey.count() << "s\n";
        std::ofstream outContext("Client_context", std::ios::binary);
        if (!outContext.is_open())
        {
            std::cerr << "Failed to open Client_context for writing" << std::endl;
            return false;
        }
        context->writeTo(outContext);
        outContext.close();
        std::cout << "Client_context has been written to file." << std::endl;

        std::ofstream outPublicKey("Client_publickey", std::ios::binary);
        if (!outPublicKey.is_open())
        {
            std::cerr << "Failed to open Client_publickey for writing" << std::endl;
            return false;
        }
        publicKey->writeTo(outPublicKey);
        outPublicKey.close();
        std::cout << "Client_publickey has been written to file." << std::endl;

        std::ofstream outSecretKey("Client_secretkey", std::ios::binary);
        if (!outSecretKey.is_open())
        {
            std::cerr << "Failed to open Client_secretkey for writing" << std::endl;
            return false;
        }
        secretKey.writeTo(outSecretKey);
        outSecretKey.close();
        std::cout << "Client secretkey has been written to file." << std::endl;

        auto start_keyEncryption = std::chrono::steady_clock::now();
        vector<Ctxt> encryptedSymKey;
        encryptSymKey(encryptedSymKey, SymKey, publicKey, ea);
        auto end_keyEncryption = std::chrono::steady_clock::now();
        double keyEncryption = std::chrono::duration<double>(end_keyEncryption - start_keyEncryption).count();
        std::cout << "SymKey FHE time: " << keyEncryption << "s\n";

        // 解密验证
        if (!verify_encryptSymKey(encryptedSymKey, SymKey, secretKey, ea))
        {
            return false;
        }
        std::cout << "Symmetric key encryption succeeded!" << std::endl;
        if (!writeEncryptedSymKey(encryptedSymKey, "Client_encryptedSymKey.bin"))
        {
            return false;
        }
        std::cout << "Client_encryptedSymKey.bin has been written to file." << std::endl;

        vector<long> slotsData;
        ea.decrypt(encryptedSymKey[0], secretKey, slotsData);
        std::cout << "decryptedSymKey[0] = " << slotsData[0] << std::endl;
        std::cout << SymKey[0]  << std::endl;
        for(int i=0;i<10;i++){
        encryptedSymKey[0].square();
        ea.decrypt(encryptedSymKey[0], secretKey, slotsData);
        std::cout << "decryptedSymKey[0]^2 = " << slotsData[0] << std::endl;
        std::cout<< (SymKey[0]*SymKey[0])%PlainMod << std::endl;
        }
        // 离线客户端时间=KeyStream Generation time+PublicKey generation and SymKey FHE time
        double total_time_off = elapsed_seconds_keyStream.count() + elapsed_seconds_PubKey.count() + keyEncryption;
        std::cout << "Encryption offline total time: " << total_time_off << "s\n";

        return true;
    }

    bool Client_online()
    {
        // 生成随机明文流
        GF2X rnd;
        int Bytebitsdiv8 = ceil(Bytebits / 8);
        vector<uint8_t> PlainStream0(Bytebitsdiv8 * PlainByte);
        random(rnd, 8 * PlainStream0.size());
        BytesFromGF2X(PlainStream0.data(), rnd, PlainStream0.size());
        vector<long> PlainStream(BlockByte);
        for (unsigned i = 0; i < BlockByte; i++)
        {
            PlainStream[i] = 0;
            for (unsigned j = 0; j < Bytebitsdiv8; j++)
            {
                PlainStream[i] += (PlainStream0[Bytebitsdiv8 * i + j] << (8 * j));
            }
            PlainStream[i] %= PlainMod;
        }

        if (!writeToFile<long>(PlainStream.data(), "Client_PlainStream.txt", PlainByte))
        {
            return false;
        }

        vector<long> KeyStream(PlainByte);
        if (!readFromFile<long>(KeyStream.data(), "Client_KeyStream.txt", PlainByte))
        {
            return false;
        }

        vector<long> CipherStream(PlainByte);
        auto start_encryption = std::chrono::steady_clock::now();
        for (unsigned i = 0; i < PlainByte; i++)
        {
            CipherStream[i] = (PlainStream[i] + KeyStream[i]) % PlainMod;
        }
        auto end_encryption = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds_encryption = end_encryption - start_encryption;
        std::cout << "Plain encryption succeeded!" << std::endl;

        if (!writeToFile(CipherStream.data(), "Client_CipherStream.txt", PlainByte))
        {
            return false;
        }
        std::cout << "Client_CipherStream has been written to file." << std::endl;
        std::cout << "Encryption online time: " << elapsed_seconds_encryption.count() << "s\n";
        return true;
    }
} // namespace C1ient_Yus_p_C32