#include <iostream>
#include "Client_Yus_p_C32.hpp"
#include "Server_Yus_p_C32.hpp"

using namespace std;
using namespace helib;
using namespace NTL;
using namespace C1ient_Yus_p_C32;
using namespace Server_Yus_p_C32;

int main()
{
    
    // 调用 Client_online 方法
#if 1

    std::cout << "==========================================" << std::endl;
    std::cout << "           Client Online Start            " << std::endl;
    std::cout << "==========================================" << std::endl;
    if (C1ient_Yus_p_C32::Client_online())
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