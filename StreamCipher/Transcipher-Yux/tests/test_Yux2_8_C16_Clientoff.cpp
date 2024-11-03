#include <iostream>
#include "Client_Yux2_8_C16.hpp"
#include "Server_Yux2_8_C16.hpp"

using namespace std;
using namespace helib;
using namespace NTL;
using namespace C1ient_Yux2_8_C16;
using namespace Server_Yux2_8_C16;

int main()
{
    
    // 输出组长度
    std::cout << "BlockByte: " << BlockByte << std::endl;
    // 输出明文分组长度
    std::cout << "PlainBlock: " << PlainBlock << std::endl;
    // 输出明文长度
    std::cout << "PlainByte: " << PlainByte << std::endl;
    // 输出轮数
    std::cout << "Nr: " << Nr << std::endl;
    // 调用 Client_offline 方法
#if 1
    std::cout << "==========================================" << std::endl;
    std::cout << "           Client Offline Start           " << std::endl;
    std::cout << "==========================================" << std::endl;

    if (C1ient_Yux2_8_C16::Client_offline())
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

    return 0;
}