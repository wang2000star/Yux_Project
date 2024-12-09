#include <iostream>
#include "Client_Yux_p_C16.hpp"
#include "Server_Yux_p_C16.hpp"

using namespace std;
using namespace helib;
using namespace NTL;
using namespace C1ient_Yux_p_C16;
using namespace Server_Yux_p_C16;

int main()
{

    // 调用 Server_online 方法
#if 1
    std::cout << "----------------------Server_online start!---------------------" << std::endl;
    if (Server_Yux_p_C16::Server_online())
    {
        std::cout << "-----------------------Server_online succeeded!-------------------" << std::endl;
    }
    else
    {
        std::cout << "---------------------Server_online failed!-------------------------" << std::endl;
        return 0;
    }
#endif

    return 0;
}