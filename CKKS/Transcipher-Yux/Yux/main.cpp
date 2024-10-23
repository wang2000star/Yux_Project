#include <iostream>
#include "Client.hpp"
#include "Server.hpp"

using namespace std;
using namespace helib;
using namespace NTL;

int main() {
    // 调用 Client_offline 方法
#if 1
    if (Client_offline()) {
        std::cout << "Client_offline succeeded!" << std::endl;  
    } else {
        std::cout << "Client_offline failed!" << std::endl;
    }
#endif

    // 调用 Client_online 方法
#if 0
    if (Client_online()) {
        std::cout << "Client_online succeeded!" << std::endl;       
    } else {
        std::cout << "Client_online failed!" << std::endl;
    }
#endif
    // 调用 Server_offline 方法
#if 0
    if(Server_offline())
    {
        std::cout << "Server_offline succeeded!" << std::endl;
    }
    else
    {
        std::cout << "Server_offline failed!" << std::endl;
    }
#endif
    // 调用 Server_online 方法
#if 1
    if(Server_online())
    {
        std::cout << "Server_online succeeded!" << std::endl;
    }
    else
    {
        std::cout << "Server_online failed!" << std::endl;
    }
    if(Client_Check()){
        std::cout << "Client_Check succeeded!" << std::endl;
    }
    else
    {
        std::cout << "Client_Check failed!" << std::endl;
    }
#endif
    return 0;
}