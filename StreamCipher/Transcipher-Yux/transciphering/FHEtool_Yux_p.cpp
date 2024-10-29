#include "FHEtool_Yux_p.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

void encodeTo16Ctxt_p(vector<vector<long>>& encData, const vector<uint64_t>& data, const EncryptedArray& ea)
{
    long R = divc(data.size(), PlainByte);
    long nCtxt = BlockByte * R;

    // 调整 encData 的大小并初始化为 0
    encData.resize(nCtxt, vector<long>(ea.size(), 0));
  for (long r=0;r<R;r++){
    for (long i=0;i<BlockByte;i++){
      vector<long> slots(ea.size(), 0);
      for (long j=0;j<PlainBlock;j++){
        long byteIdx = j*BlockByte + i + r*PlainByte;
        if (byteIdx < data.size()){
          slots[j] = data[byteIdx];
        }
      }
      memcpy(encData[r*BlockByte+i].data(), slots.data(), slots.size()*sizeof(long));
    }
  }
}
// encodeTo16Ctxt_p对应的解码
void decodeTo16Ctxt_p(vector<uint64_t>& data, const vector<vector<long>>& encData, const EncryptedArray& ea)
{
    long nCtxt = encData.size();
    long R = nCtxt / BlockByte;
    data.clear();
    data.resize(R * PlainByte, 0);

    for (long r = 0; r < R; r++) {
        for (long i = 0; i < BlockByte; i++) {
            vector<long> slots(ea.size(), 0);
            memcpy(slots.data(), encData[r * BlockByte + i].data(), slots.size() * sizeof(long));
            for (long j = 0; j < PlainBlock; j++) {
                long byteIdx = j * BlockByte + i + r * PlainByte;
                if (byteIdx < data.size()) {
                    data[byteIdx] = slots[j];
                }
            }
        }
    }
}

// 函数：解密并验证密文是否正确，需要解码
bool verifyDecryption_p16(const vector<Ctxt>& encryptedSymKey, const vector<uint64_t>& originalSymKey, const SecKey& secretKey, const EncryptedArray& ea)
{
    // 解密加密的对称密钥
    vector<vector<long>> decryptedSymKey(encryptedSymKey.size());
    for (size_t i = 0; i < encryptedSymKey.size(); ++i)
    {
        ea.decrypt(encryptedSymKey[i], secretKey, decryptedSymKey[i]);
    }

    // 验证解密后的对称密钥是否与原始对称密钥匹配
    bool isMatch = true;
    for (size_t i = 0; i < originalSymKey.size(); ++i)
    {
        bool slotMatch = false;
        for (const auto& slotValue : decryptedSymKey[i])
        {
            if (slotValue == originalSymKey[i])
            {
                slotMatch = true;
                break;
            }
        }
        if (!slotMatch)
        {
            isMatch = false;
            std::cerr << "Mismatch at index " << i << ": expected " << originalSymKey[i] << ", got ";
            for (const auto& slotValue : decryptedSymKey[i])
            {
                std::cerr << slotValue << " ";
            }
            std::cerr << std::endl;
        }
    }

    return isMatch;
}