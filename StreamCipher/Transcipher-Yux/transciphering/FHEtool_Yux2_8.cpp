#include "FHEtool_Yux2_8.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

// Encode plaintext/ciphertext bytes as native HE plaintext
// 编码明文/密文字节为本机HE明文
void encodeTo1Ctxt(Vec<ZZX>& encData, const Vec<uint8_t>& data,
		                            const EncryptedArrayDerived<PA_GF2>& ea)
{
  long nAllBlocks = divc(data.length(),BlockByte); // ceil( data.length()/BlockByte )
  long blocksPerCtxt = ea.size() / BlockByte;  // = nSlots/BlockByte
  long nCtxt = divc(nAllBlocks, blocksPerCtxt);
  // We encode blocksPerCtxt = n/BlockByte blocks in the slots of one ctxt.
  encData.SetLength(nCtxt);

  for (long i=0; i<nCtxt; i++) {         // i is the cipehrtext number
    // Copy the bytes into Hypercube<GF2X>'es to be used for encoding
    vector<GF2X> slots(ea.size(), GF2X::zero());
    for (long j=0; j<blocksPerCtxt; j++) { // j is the block number in this ctxt
      long blockShift = (i*blocksPerCtxt +j)*BlockByte;  // point to block
      for (long k=0; k<BlockByte; k++) {         // k is the byte number in this block
        long byteIdx= blockShift+ k;      // column orded within block
        if (byteIdx < data.length()) {
          long slotIdx = j + k*blocksPerCtxt;
          //ong slotIdx = k + j*BlockByte;
          GF2XFromBytes(slots[slotIdx], &data[byteIdx], 1);// copy byte as poly
        }
      }
    }
    ea.encode(encData[i], slots);
  }
}

// Decode native HE plaintext/ciphertext as Yux plaintext/ciphertext bytes
// 将本机HE明文/密文解码为Yux明文/密文字节，就是把多项式转换为字节数组
// 
void decodeTo1Ctxt(Vec<uint8_t>& data, const Vec<ZZX>& encData,
		const EncryptedArrayDerived<PA_GF2>& ea)
{
  // Check the size of the data array
  long nBytes = encData.length() * ea.size(); // total number of input bytes
  if (data.length()<=0 || data.length()>nBytes)
    data.SetLength(nBytes);
  long nAllBlocks = divc(data.length(),16);       // ceil( data.length()/16 )
  long blocksPerCtxt = ea.size() / 16;        // = nSlots/16
  long nCtxt = divc(nAllBlocks, blocksPerCtxt);   // <= encData.length()

  // We encode blocksPerCtxt = n/16 blocks in the slots of one ctxt.

  vector<GF2X> slots;
  for (long i=0; i<nCtxt; i++) {         // i is the cipehrtext number
    ea.decode(slots, encData[i]);
    for (long j=0; j<blocksPerCtxt; j++) { // j is the block number in this ctxt
      long blockShift = (i*blocksPerCtxt +j)*16;  // point to block
      for (long k=0; k<16; k++) {         // k is the byte number in this block
        long byteIdx= blockShift +k;      // column orded within block
        if (byteIdx < data.length()) {
          long slotIdx = j + k*blocksPerCtxt;
          //long slotIdx = k + j*BlockByte;
          BytesFromGF2X(&data[byteIdx], slots[slotIdx], 1);// copy poly as byte
        }
      }
    }
  }
}

// 函数：解密并验证密文是否正确
bool verifyDecryption1(const vector<Ctxt> &encryptedVec, const SecKey &secretKey,
                      const EncryptedArrayDerived<PA_GF2> &ea, const Vec<uint8_t> &originalVec)
{
    // 定义解密开始时间
    auto start_decrypt = std::chrono::steady_clock::now();

    Vec<uint8_t> decryptedVec(INIT_SIZE, originalVec.length());
    Vec<ZZX> decryptedPolys(INIT_SIZE, encryptedVec.size());

    // 使用 secretKey 逐个解密 encryptedVec
    for (long i = 0; i < encryptedVec.size(); ++i)
    {
        secretKey.Decrypt(decryptedPolys[i], encryptedVec[i]);
    }

    // 解码多项式到明文字节
    decodeTo1Ctxt(decryptedVec, decryptedPolys, ea);

    // 验证解密结果是否与原始明文一致
    bool isDecryptionCorrect = true;
    for (long i = 0; i < originalVec.length(); ++i)
    {
        if (decryptedVec[i] != originalVec[i])
        {
            std::cout << "Decryption check failed at index " << i << ": expected " << (int)originalVec[i]
                      << ", got " << (int)decryptedVec[i] << std::endl;
            isDecryptionCorrect = false;
            break;
        }
    }

    if (isDecryptionCorrect)
    {
        std::cout << "Decryption check succeeded: Decrypted vector matches original vector." << std::endl;
    }
    else
    {
        std::cout << "Decryption check failed: Decrypted vector does not match original vector." << std::endl;
    }

    // 定义解密结束时间
    auto end_decrypt = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_decrypt = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification time: " << elapsed_seconds_decrypt.count() << "s\n";

    return isDecryptionCorrect;
}

// Encode  plaintext/ciphertext bytes as native HE plaintext
// packing
void encodeTo16Ctxt(Vec<ZZX>& encData, const Vec<uint8_t>& data,
		const EncryptedArrayDerived<PA_GF2>& ea)
{
  long R = divc(data.length(), PlainByte);
  long nCtxt = BlockByte * R;
  encData.SetLength(nCtxt);
  for (long r = 0; r < R; r++)
  {
    for (long i = 0; i < BlockByte; i++)
    {
      vector<GF2X> slots(ea.size(), GF2X::zero());
      for (long j = 0; j < PlainBlock; j++)
      { // j is the block number in this ctxt
        long byteIdx = j * BlockByte + i + r * PlainByte;
        if (byteIdx < data.length())
        {
          GF2XFromBytes(slots[j], &data[byteIdx], 1);
        }
      }
      ea.encode(encData[r*BlockByte+i], slots);
    }  
  }
}

// Decode native HE plaintext as Yux plaintext/ciphertext bytes
void decodeTo16Ctxt2(Vec<uint8_t>& data, const Vec<ZZX>& encData,
		const EncryptedArrayDerived<PA_GF2>& ea)
{
  // Check the size of the data array
  long nAllBytes = encData.length() * ea.size(); // total number of input bytes
  if (data.length()<=0 || data.length()>nAllBytes)
    data.SetLength(nAllBytes);

  // 一个分组有16个字节
  long nBytes = BlockByte;
  long nAllBlocks = divc(data.length(), nBytes); // ceil( data.length()/16 ) =(a + b - 1)/b
  // long blocksPerCtxt = ea.size() / 16;  // = nSlots/16
  long nCtxt = nBytes;

  vector<GF2X> slots;
  for (long i=0; i<nCtxt; i++) {
    ea.decode(slots, encData[i]);
    for(long j=0; j<nAllBlocks; j++) {
      long slotIdx = j;
      long byteIdx = i + nBytes*j;
      if (byteIdx < data.length()) {
        BytesFromGF2X(&data[byteIdx], slots[slotIdx], 1);
      }
    }
  }
}
void decodeTo16Ctxt(Vec<uint8_t>& data, const Vec<ZZX>& encData,
                      const EncryptedArrayDerived<PA_GF2>& ea)
{
    long R = encData.length() / BlockByte;
    long nBytes = R * PlainByte;

    data.SetLength(nBytes);

    for (long r = 0; r < R; r++)
    {
        for (long i = 0; i < BlockByte; i++)
        { // i is the ciphertext number
            vector<GF2X> slots;
            ea.decode(slots, encData[r * BlockByte + i]);

            for (long j = 0; j < PlainBlock; j++)
            { // j is the block number in this ctxt
                long byteIdx = j * BlockByte + i + r * PlainByte;
                if (byteIdx < data.length())
                {
                    BytesFromGF2X(&data[byteIdx], slots[j], 1); // copy poly as byte
                }
            }
        }
    }
}
// 函数：解密并验证密文是否正确
bool verifyDecryption16(const vector<Ctxt> &encryptedVec, const SecKey &secretKey,
                      const EncryptedArrayDerived<PA_GF2> &ea, const Vec<uint8_t> &originalVec)
{
    // 定义解密开始时间
    auto start_decrypt = std::chrono::steady_clock::now();

    Vec<uint8_t> decryptedVec(INIT_SIZE, originalVec.length());
    Vec<ZZX> decryptedPolys(INIT_SIZE, encryptedVec.size());

    // 使用 secretKey 逐个解密 encryptedVec
    for (long i = 0; i < encryptedVec.size(); ++i)
    {
        secretKey.Decrypt(decryptedPolys[i], encryptedVec[i]);
    }

    // 解码多项式到明文字节
    decodeTo16Ctxt(decryptedVec, decryptedPolys, ea);

    // 验证解密结果是否与原始明文一致
    bool isDecryptionCorrect = true;
    for (long i = 0; i < originalVec.length(); ++i)
    {
        if (decryptedVec[i] != originalVec[i])
        {
            std::cout << "Decryption check failed at index " << i << ": expected " << (int)originalVec[i]
                      << ", got " << (int)decryptedVec[i] << std::endl;
            isDecryptionCorrect = false;
            break;
        }
    }

    if (isDecryptionCorrect)
    {
        std::cout << "Decryption check succeeded: Decrypted vector matches original vector." << std::endl;
    }
    else
    {
        std::cout << "Decryption check failed: Decrypted vector does not match original vector." << std::endl;
    }

    // 定义解密结束时间
    auto end_decrypt = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_decrypt = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification time: " << elapsed_seconds_decrypt.count() << "s\n";

    return isDecryptionCorrect;
}
