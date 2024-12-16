#include "FHEtool_Yus_p.hpp"

using namespace std;
using namespace helib;
using namespace NTL;


void encodeTo32Ctxt(vector<ZZX> &encData, const vector<long> &data, const EncryptedArray &ea)
{
    long R = data.size()/PlainByte;
    std::cout << "R: " << R << std::endl;
    long nCtxt = BlockByte * R;
    std::cout << "nCtxt: " << nCtxt << std::endl;
    long data_length = data.size();
    encData.resize(nCtxt);
    std::cout << "data_length: " << data_length << std::endl;
    std::cout << "ea.size: " << ea.size() << std::endl;
    
    for (long r = 0; r < R; r++)
    {
        for (long i = 0; i < BlockByte; i++)
        {
            vector<long> slots(ea.size(), 0);
            for (long j = 0; j < PlainBlock; j++)
            {
                long byteIdx = j * BlockByte + i + r * PlainByte;
                if (byteIdx < data_length)
                {
                    slots[j] = data[byteIdx];
                }
            }
            ea.encode(encData[r * BlockByte + i], slots);
        }
    }
}
// encodeTo32Ctxt对应的解码
void decodeTo32Ctxt(vector<long> &data, const vector<vector<long>> &encData,
                    const EncryptedArray &ea)
{
  long R = encData.size() / BlockByte;
  std::cout << "R: " << R << std::endl;
  long data_length = R * PlainByte;

  data.resize(data_length);

  for (long r = 0; r < R; r++)
  {
    for (long j = 0; j < PlainBlock; j++){
    for (long i = 0; i < BlockByte; i++)
    { // i is the ciphertext number
      
       // j is the block number in this ctxt
        long byteIdx = j * BlockByte + i + r * PlainByte;
        if (byteIdx < data_length)
        {
          data[byteIdx] = encData[r * BlockByte + i][j];
        }
      
    }
  }
}
}
// 函数：解密并验证密文是否正确，需要解码
// 函数：解密并验证密文是否正确
bool verifyDecryption32(const std::vector<Ctxt> &encryptedVec, const vector<long> &originalVec, const SecKey &secretKey,
                        const EncryptedArray &ea)
{
    auto start_decrypt = std::chrono::steady_clock::now();

    vector<long> decryptedVec;
    std::vector<std::vector<long>> decryptedPolys(encryptedVec.size());

    for (std::size_t i = 0; i < encryptedVec.size(); ++i)
    {
        ea.decrypt(encryptedVec[i], secretKey, decryptedPolys[i]);
    }
    // 解码
    decodeTo32Ctxt(decryptedVec, decryptedPolys, ea);
    // 验证解密结果
    bool isDecryptedVecCorrect = true;
    long len = originalVec.size();
    for (std::size_t i = 0; i < len; ++i)
    {
        if (decryptedVec[i] != originalVec[i])
        {
            std::cout << "Decryption check failed at index " << i << ": expected " << originalVec[i]
                      << ", got " << decryptedVec[i] << std::endl;
            isDecryptedVecCorrect = false;
            break;
        }
    }

    auto end_decrypt = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification succeeded! Time: " << elapsed_seconds.count() << "s\n";

    return isDecryptedVecCorrect;
}