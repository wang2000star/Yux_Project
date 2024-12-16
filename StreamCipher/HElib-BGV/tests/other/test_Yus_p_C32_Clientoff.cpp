#include <iostream>
#include "Client_Yus_p_C32.hpp"
//#include "Server_Yus_p_C32.hpp"

using namespace std;
using namespace helib;
using namespace NTL;
using namespace C1ient_Yus_p_C32;
//using namespace Server_Yus_p_C32;

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

    if (C1ient_Yus_p_C32::Client_offline())
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