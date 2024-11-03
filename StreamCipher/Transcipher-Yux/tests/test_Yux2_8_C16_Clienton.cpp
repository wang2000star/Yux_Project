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
    // 调用 Client_online 方法
#if 1

    std::cout << "==========================================" << std::endl;
    std::cout << "           Client Online Start            " << std::endl;
    std::cout << "==========================================" << std::endl;
    if (C1ient_Yux2_8_C16::Client_online())
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

    return 0;
}