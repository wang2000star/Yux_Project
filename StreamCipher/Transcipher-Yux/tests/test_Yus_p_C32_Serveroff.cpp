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
   
    // 调用 Server_offline 方法
#if 1
    std::cout << "---------------------Server_offline start!------------------------" << std::endl;
    if (Server_Yus_p_C32::Server_offline())
    {
        std::cout << "-----------------------Server_offline succeeded!------------------" << std::endl;
    }
    else
    {
        std::cout << "------------------------Server_offline failed!----------------------" << std::endl;
        return 0;
    }
#endif

    return 0;
}