#include <iostream>
#include "Client.hpp"
#include "Server.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

// 读取加密密钥流
extern bool readEncryptedKeyStream(vector<Ctxt>& encryptedKeyStream, const std::string& filename, const PubKey& publicKey);

int main()
{
    // 调用 Client_offline 方法
#if 1
    std::cout << "==========================================" << std::endl;
    std::cout << "           Client Offline Start           " << std::endl;
    std::cout << "==========================================" << std::endl;

    if (Client_offline())
    {
        std::cout << "==========================================" << std::endl;
        std::cout << "         Client Offline Succeeded         " << std::endl;
        std::cout << "==========================================" << std::endl;
    }
    else
    {
        std::cout << "==========================================" << std::endl;
        std::cout << "          Client Offline Failed           " << std::endl;
        std::cout << "==========================================" << std::endl;
        return 0;
    }
#endif
    // 调用 Server_offline 方法
#if 1
    std::cout << "---------------------Server_offline start!------------------------" << std::endl;
    if (Server_offline())
    {
        std::cout << "-----------------------Server_offline succeeded!------------------" << std::endl;
    }
    else
    {
        std::cout << "------------------------Server_offline failed!----------------------" << std::endl;
        return 0;
    }
#endif
    // 调用 Client_online 方法
#if 1

    std::cout << "==========================================" << std::endl;
    std::cout << "           Client Online Start            " << std::endl;
    std::cout << "==========================================" << std::endl;

    if (Client_online())
    {
        std::cout << "==========================================" << std::endl;
        std::cout << "         Client Online Succeeded          " << std::endl;
        std::cout << "==========================================" << std::endl;
    }
    else
    {
        std::cout << "==========================================" << std::endl;
        std::cout << "          Client Online Failed            " << std::endl;
        std::cout << "==========================================" << std::endl;
        return 0;
    }
#endif

    // 调用 Server_online 方法
#if 1
    std::cout << "----------------------Server_online start!---------------------" << std::endl;
    if (Server_online())
    {
        std::cout << "-----------------------Server_online succeeded!-------------------" << std::endl;
    }
    else
    {
        std::cout << "---------------------Server_online failed!-------------------------" << std::endl;
        return 0;
    }
#endif

// ===============检查服务端是否正确解密================
// BGV预设
#if 1
     // 从文件中读取上下文
    std::ifstream inContext("Client_context", std::ios::binary);
    if (!inContext.is_open()) {
        std::cerr << "Failed to open Client_context for reading" << std::endl;
        throw std::runtime_error("Failed to open context file");
    }
    Context context = Context::readFrom(inContext);
    inContext.close();
    
    // 从文件中读取公钥
    std::ifstream inPublicKey("Client_publickey", std::ios::binary);
    if (!inPublicKey.is_open()) {
        std::cerr << "Failed to open Client_publickey for reading" << std::endl;
        throw std::runtime_error("Failed to open public key file");
    }
    PubKey publicKey = PubKey::readFrom(inPublicKey, context);
    inPublicKey.close();
    // 从文件中读取私钥
    std::ifstream inSecretKey("Client_secretkey", std::ios::binary);
    if (!inSecretKey.is_open()) {
        std::cerr << "Failed to open Client_secretkey for reading" << std::endl;
        throw std::runtime_error("Failed to open secret key file");
    }
    SecKey secretKey = SecKey::readFrom(inSecretKey, context);
    inSecretKey.close();

    uint8_t YuxPolyBytes[] = {0x1B, 0x1}; // X^8+X^4+X^3+X+1
    const GF2X YuxPoly = GF2XFromBytes(YuxPolyBytes, 2);
    EncryptedArrayDerived<PA_GF2> ea(context, YuxPoly, context.getAlMod());
#endif
// 从Server_encryptedPlainStream.txt读取encryptedPlainStream
#if 1
    vector<Ctxt> encryptedPlainStream;
    if (!readEncryptedKeyStream(encryptedPlainStream, "Server_encryptedPlainStream.bin", publicKey))
    {
        return false;
    }
#endif
// 从Client_plainStream.txt读取PlainStream
#if 1
    Vec<uint8_t> PlainStream(INIT_SIZE, PlainByte);
    std::ifstream in4("Client_PlainStream.txt", std::ios::binary);
    if (in4.is_open())
    {
        in4.read((char *)PlainStream.data(), PlainByte);
        in4.close();
    }
    else
    {
        std::cerr << "Failed to open Client_PlainStream.txt for reading" << std::endl;
        return false;
    }
#endif
// 解密encryptedPlainStream得到PlainStream2，然后与PlainStream比较
#if 1
    std::cout << "-----------------Check Server_decryption start!-----------------" << std::endl;
    for (int i = 0; i < encryptedPlainStream.size(); i++)
        encryptedPlainStream[i].bringToSet(encryptedPlainStream[i].naturalPrimeSet());
        
    Vec<uint8_t> decryptedPlainStream(INIT_SIZE, PlainByte);
    Vec<ZZX> poly(INIT_SIZE, encryptedPlainStream.size());
    for (long i = 0; i < poly.length(); i++)
    {
        secretKey.Decrypt(poly[i], encryptedPlainStream[i]);
    }
    decodeTo1Ctxt(decryptedPlainStream, poly, ea);

    // 验证解密后的密文是否与原始密文相同
    for (long i = 0; i < PlainByte; i++)
    {
        if (decryptedPlainStream[i] != PlainStream[i])
        {
            std::cout << "Check Server_decryption failed!" << std::endl;
            return 0;
        }
    }
    std::cout << "-----------------Check Server_decryption succeeded!-----------------" << std::endl;
#endif
    return 0;
}