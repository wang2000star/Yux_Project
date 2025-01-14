#include "tool.hpp"

// 生成m比特强度的随机整数
long generate_secure_random_int(unsigned m) {
    if (m > 64) {
        throw std::invalid_argument("m exceeds the maximum bit size of long");
    }

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<long> dis(0, (1ULL << m) - 1);

    return dis(gen);
}
void random_init_shake(long nonce, long block_counter, Keccak_HashInstance &shake128_)
{
    uint8_t seed[16];

    *((long *)seed) = htobe64(nonce);
    *((long *)(seed + 8)) = htobe64(block_counter);

    if (SUCCESS != Keccak_HashInitialize_SHAKE128(&shake128_))
        throw std::runtime_error("failed to init shake");
    if (SUCCESS != Keccak_HashUpdate(&shake128_, seed, sizeof(seed) * 8))
        throw std::runtime_error("SHAKE128 update failed");
    if (SUCCESS != Keccak_HashFinal(&shake128_, NULL))
        throw std::runtime_error("SHAKE128 final failed");
}
long generate_random_field_element(Keccak_HashInstance &shake128, bool allow_zero, long max_prime_size, long PlainMod)
{
    uint8_t random_bytes[sizeof(long)];
    while (1)
    {
        if (SUCCESS !=
            Keccak_HashSqueeze(&shake128, random_bytes, sizeof(random_bytes) * 8))
            throw std::runtime_error("SHAKE128 squeeze failed");
        long ele = be64toh(*((long *)random_bytes)) & max_prime_size;
        if (!allow_zero && ele == 0)
            continue;
        if (ele < PlainMod)
            return ele;
    }
}

int min_noise_budget(vector<Ciphertext> &eData ,Decryptor &decryptor)
{
    int min_noise = 1000;
    long size = eData.size();
    for (int i = 0; i < size; i++)
    {
        int noise = decryptor.invariant_noise_budget(eData[i]);
        if (noise < min_noise)
        {
            min_noise = noise;
        }
    }
    return min_noise;
}


bool writeEncryptedSymKey(const vector<Ciphertext> &encryptedSymKey, const std::string &filename)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open())
    {
        std::cerr << "Failed to open " << filename << " for writing" << std::endl;
        return false;
    }

    for (const auto &ctxt : encryptedSymKey)
    {
        ctxt.save(out);
    }

    out.close();

    return true;
}

void decodeToCtxt(std::vector<long> &data, const std::vector<std::vector<long>> &encData, const long CtxtWords,const long nslots)
{
    long R = encData.size() / CtxtWords;
    long AllByte = CtxtWords * nslots;
    long data_size = R * AllByte;
    data.resize(data_size);
    long byteIdx;
    long rAB;
    long jCB;
    long rCB;
    //    omp_set_num_threads(16); // 设置线程数为16
    // #pragma omp parallel for
    for (long j = 0; j < nslots; j++)
    {
        jCB = j * CtxtWords;
        for (long r = 0; r < R; r++)
        {
            rAB = r * AllByte;
            rCB = r * CtxtWords;
            for (long i = 0; i < CtxtWords; i++)
            { // i is the ciphertext number
                // j is the block number in this ctxt
                byteIdx = jCB + i + rAB;
                data[byteIdx] = encData[rCB + i][j];
            }
        }
    }
}

// 函数：解密并验证密文是否正确，需要解码
bool verifyDecryption(const std::vector<Ciphertext> &encryptedVec, const std::vector<long> &originalVec, BatchEncoder &batch_encoder, Decryptor &decryptor,
                      const long CtxtWords, const long PlainBlock, const long nslots, const long pmod)
{
    std::vector tempVec = originalVec;
    for (int i = 0; i < originalVec.size(); i++)
    {
        tempVec[i] = (tempVec[i] + pmod) % pmod;
    }
    int size = encryptedVec.size();
    auto start_decrypt = std::chrono::high_resolution_clock::now();
    std::vector<long> decryptedVec = originalVec;
    Plaintext decrypted_result;
    std::vector<long> pod_result;
    std::vector<std::vector<long>> decryptedPolys(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        decryptor.decrypt(encryptedVec[i], decrypted_result);
        batch_encoder.decode(decrypted_result, pod_result);
        decryptedPolys[i] = pod_result;
    }
    decodeToCtxt(decryptedVec, decryptedPolys, CtxtWords, nslots);
    for (int i = 0; i < decryptedVec.size(); i++)
    {
        decryptedVec[i] = (decryptedVec[i] + pmod) % pmod;
    }
    // 验证解密结果
    bool isDecryptedVecCorrect = std::equal(decryptedVec.begin(), decryptedVec.begin()+PlainBlock*CtxtWords, tempVec.begin());
    auto end_decrypt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification finished! Time: " << elapsed_seconds.count() << "s\n";
    // 如果解密结果不正确，输出错误的位置
    if (!isDecryptedVecCorrect)
    {
        for (size_t i = 0; i < CtxtWords*2; i++)
        {
            if (decryptedVec[i] != tempVec[i])
            {
                std::cout << "Error at position " << i << ": " << decryptedVec[i] << " != " << originalVec[i] << std::endl;
                // break;
            }
        }
    }
    return isDecryptedVecCorrect;
}

void encryptSymKey(vector<Ciphertext> &encryptedSymKey, const vector<long> &SymKey, BatchEncoder &batch_encoder, Encryptor &encryptor,const long nslots)
{
    long BlockWords = SymKey.size();
    // 加密
    for (long i = 0; i < BlockWords; i++)
    { // encrypt the encoded key
        vector<long> SymKey_temp(nslots, SymKey[i]);
        Plaintext SymKey_plain;
        batch_encoder.encode(SymKey_temp, SymKey_plain);
        encryptor.encrypt(SymKey_plain,encryptedSymKey[i]);
    }
}

bool verify_encryptSymKey(vector<Ciphertext> &encryptedSymKey, const vector<long> &SymKey,BatchEncoder &batch_encoder, Decryptor &decryptor)
{
    auto start_decrypt = std::chrono::high_resolution_clock::now();
    long BlockWords = SymKey.size();
    std::vector<long> decryptedSymKey = SymKey;
    Plaintext decrypted_result;
    vector<long> pod_result;
    //     omp_set_num_threads(16); // 设置线程数为16
    // #pragma omp parallel for
    for (long i = 0; i < BlockWords; i++)
    { // encrypt the encoded key
        decryptor.decrypt(encryptedSymKey[i], decrypted_result);
        batch_encoder.decode(decrypted_result, pod_result);
        decryptedSymKey[i] = pod_result[0];
    }
    bool isDecryptedSymKeyCorrect = std::equal(SymKey.begin(), SymKey.end(), decryptedSymKey.begin());
    auto end_decrypt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification finished! Time: " << elapsed_seconds.count() << "s\n";
    return isDecryptedSymKeyCorrect;
}