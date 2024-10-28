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
    auto start_test_time = std::chrono::steady_clock::now();
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
    // 调用 Server_offline 方法
#if 1
    std::cout << "---------------------Server_offline start!------------------------" << std::endl;
    if (Server_Yux2_8_C16::Server_offline())
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

    // 调用 Server_online 方法
#if 1
    std::cout << "----------------------Server_online start!---------------------" << std::endl;
    if (Server_Yux2_8_C16::Server_online())
    {
        std::cout << "-----------------------Server_online succeeded!-------------------" << std::endl;
    }
    else
    {
        std::cout << "---------------------Server_online failed!-------------------------" << std::endl;
        return 0;
    }
#endif
    auto end_test_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_test_time = end_test_time - start_test_time;
    std::cout << "Test time: " << elapsed_seconds_test_time.count() << "s\n";

    return 0;
}