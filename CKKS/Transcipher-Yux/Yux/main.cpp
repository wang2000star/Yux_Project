#include <iostream>
#include "Client.hpp"
#include "Server.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

int main()
{
    // 调用 Client_offline 方法
#if 1
    if (Client_offline())
    {
        std::cout << "Client_offline succeeded!" << std::endl;
    }
    else
    {
        std::cout << "Client_offline failed!" << std::endl;
    }
#endif

    // 调用 Client_online 方法
#if 1
    if (Client_online())
    {
        std::cout << "Client_online succeeded!" << std::endl;
    }
    else
    {
        std::cout << "Client_online failed!" << std::endl;
    }
#endif
    // 调用 Server_offline 方法
#if 0
    if (Server_offline())
    {
        std::cout << "Server_offline succeeded!" << std::endl;
    }
    else
    {
        std::cout << "Server_offline failed!" << std::endl;
    }
#endif
    // 调用 Server_online 方法
#if 0
    if (Server_online())
    {
        std::cout << "Server_online succeeded!" << std::endl;
    }
    else
    {
        std::cout << "Server_online failed!" << std::endl;
    }
#endif

// ===============检查服务端是否正确解密================
#if 0
// BGV预设
#if 1
    long idx = 3; // 0
    long c = 9;
    bool packed = true;
    if (idx > 5)
        idx = 5;
    long p = mValues[idx][0];
    long m = mValues[idx][2];
    long bits = mValues[idx][4];

    Context context(ContextBuilder<BGV>()
                        .m(m)
                        .p(p)
                        .r(1)
                        .c(c)
                        .bits(bits)
                        .build());
    // initialize context
    PubKey publicKey(context);
    SecKey secretKey(context);
    uint8_t YuxPolyBytes[] = {0x1B, 0x1}; // X^8+X^4+X^3+X+1
    const GF2X YuxPoly = GF2XFromBytes(YuxPolyBytes, 2);
    EncryptedArrayDerived<PA_GF2> ea(context, YuxPoly, context.getAlMod());
#endif
// 从Client_publickey读取publicKey, 从Client_secretKey读取secretKey
#if 1
    std::ifstream in2("Client_publickey", std::ios::binary);
    if (!in2.is_open())
    {
        std::cerr << "Failed to open Client_publickey for reading" << std::endl;
        return false;
    }
    publicKey.readFrom(in2, context);
    in2.close();
    std::ifstream in4("Client_secretKey", std::ios::binary);
    if (!in4.is_open())
    {
        std::cerr << "Failed to open Client_secretKey for reading" << std::endl;
        return false;
    }
    secretKey.readFrom(in4, context);
    in4.close();
#endif
// 从Server_encryptedPlainStream.txt读取encryptedPlainStream
#if 1
    std::vector<Ctxt> encryptedPlainStream;
    std::ifstream in("Server_encryptedPlainStream.txt", std::ios::binary);
    if (in.is_open())
    {
        while (in.peek() != EOF)
        {
            Ctxt ctxt(publicKey);
            ctxt.read(in);
            encryptedPlainStream.push_back(ctxt);
        }
        in.close();
    }
    else
    {
        std::cerr << "Failed to open Server_encryptedPlainStream.txt for reading" << std::endl;
        return false;
    }
#endif
// 从Client_plainStream.txt读取PlainStream
#if 1
    Vec<uint8_t> PlainStream(INIT_SIZE, PlainByte);
    std::ifstream in3("Client_plainStream.txt", std::ios::binary);
    if (in3.is_open())
    {
        in3.read((char *)PlainStream.data(), PlainByte);
        in3.close();
    }
    else
    {
        std::cerr << "Failed to open Client_plainStream.txt for reading" << std::endl;
        return false;
    }
#endif
// 解密encryptedPlainStream得到PlainStream2，然后与PlainStream比较
#if 1
    Vec<uint8_t> decryptedPlainStream(INIT_SIZE, PlainByte);
    for (long i = 0; i < PlainByte; i++)
    {
        decryptedPlainStream[i] = 0;
    }

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
            std::cout << "PlainStream FHE failed!" << std::endl;
            break;
        }
    }
#endif
#endif
    return 0;
}