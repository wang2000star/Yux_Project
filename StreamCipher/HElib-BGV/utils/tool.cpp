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


int min_noise_budget(std::vector<helib::Ctxt> &eData)
{
    int min_noise = 10000;
    int noise;
    for (int i = 0; i < eData.size(); i++)
    {
        noise = eData[i].bitCapacity();
        if (noise < min_noise)
        {
            min_noise = noise;
        }
    }
    return min_noise;
}
// 函数：对多项式的每个系数乘以整数 a 并取模 c
helib::zzX multiplyAndMod(const helib::zzX &a, long b,long pmod)
{
    helib::zzX res;
    res.SetLength(helib::lsize(a));
    int len = helib::lsize(a);
    for (long i = 0; i < len; ++i)
    {
        res[i] = (a[i] * b) % pmod;
    }
    return res;
}

bool writeEncryptedSymKey(const std::vector<helib::Ctxt> &encryptedSymKey, const std::string &filename)
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

void decodeToCtxt(std::vector<long> &data, const std::vector<NTL::vec_long> &encData, const long CtxtWords,const long PlainBlock,const long nslots)
{
    long R = encData.size() / CtxtWords;
    long AllByte = CtxtWords * PlainBlock;
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
bool verifyDecryption(const std::vector<helib::Ctxt> &encryptedVec, const std::vector<long> &originalVec, const helib::SecKey &secretKey,
                      const helib::Cmodulus &cmodulus, const long CtxtWords, const long PlainBlock, const long nslots, const long pmod)
{
    std::vector tempVec = originalVec;
    for (int i = 0; i < originalVec.size(); i++)
    {
        tempVec[i] = (tempVec[i] + pmod) % pmod;
    }
    int size = encryptedVec.size();
    auto start_decrypt = std::chrono::high_resolution_clock::now();
    std::vector<long> decryptedVec = originalVec;
    std::vector<NTL::vec_long> decryptedPolys(size);
    std::vector<NTL::ZZX> Polys(size);
    //     omp_set_num_threads(16); // 设置线程数为16
    // #pragma omp parallel for
    for (std::size_t i = 0; i < size; ++i)
    {
        secretKey.Decrypt(Polys[i], encryptedVec[i]);
        cmodulus.FFT(decryptedPolys[i], Polys[i]);
    }
    decodeToCtxt(decryptedVec, decryptedPolys, CtxtWords, PlainBlock, nslots);
    // 验证解密结果
    bool isDecryptedVecCorrect = std::equal(decryptedVec.begin(), decryptedVec.end(), tempVec.begin());
    auto end_decrypt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification finished! Time: " << elapsed_seconds.count() << "s\n";
    // 如果解密结果不正确，输出错误的位置
    if (!isDecryptedVecCorrect)
    {
        for (size_t i = 0; i < CtxtWords; i++)
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

void encryptSymKey(std::vector<helib::Ctxt> &encryptedSymKey, const std::vector<long> &SymKey, std::unique_ptr<helib::PubKey> &pk, const helib::Cmodulus &cmodulus,const long nslots)
{
    NTL::vec_long slotsData;
    slotsData.SetLength(nslots);
    long BlockWords = SymKey.size();
    encryptedSymKey.resize(BlockWords, helib::Ctxt(*pk));
    NTL::zz_pX temp;
    NTL::ZZX encodedData;
    for (long i = 0; i < BlockWords; i++)
    { // encrypt the encoded key
        for (long j = 0; j < nslots; j++)
        {
            slotsData[j] = SymKey[i];
        }
        cmodulus.iFFT(temp, slotsData);
        conv(encodedData, temp);
        pk->Encrypt(encryptedSymKey[i], encodedData);
    }
}

bool verify_encryptSymKey(std::vector<helib::Ctxt> &encryptedSymKey, const std::vector<long> &SymKey, const helib::SecKey &secretKey,
                          const helib::Cmodulus &cmodulus)
{
    auto start_decrypt = std::chrono::high_resolution_clock::now();
    long BlockWords = SymKey.size();
    std::vector<long> decryptedSymKey = SymKey;
    NTL::ZZX encodedSymKey;
    NTL::vec_long slotsData;
    //     omp_set_num_threads(16); // 设置线程数为16
    // #pragma omp parallel for
    for (long i = 0; i < BlockWords; i++)
    { // encrypt the encoded key
        secretKey.Decrypt(encodedSymKey, encryptedSymKey[i]);
        cmodulus.FFT(slotsData, encodedSymKey);
        decryptedSymKey[i] = slotsData[0];
    }
    bool isDecryptedSymKeyCorrect = std::equal(SymKey.begin(), SymKey.end(), decryptedSymKey.begin());
    auto end_decrypt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification finished! Time: " << elapsed_seconds.count() << "s\n";
    return isDecryptedSymKeyCorrect;
}