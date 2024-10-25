#include "FHEtool.hpp"

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
          BytesFromGF2X(&data[byteIdx], slots[slotIdx], 1);// copy poly as byte
        }
      }
    }
  }
}