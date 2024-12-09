#ifndef SERVER_YUS_P_C32_HPP
#define SERVER_YUS_P_C32_HPP


#include <iostream>
#include <cstring>
#include <stdint.h>
#include <chrono>
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <vector>

#include <helib/helib.h>
#include <helib/ArgMap.h>
#include <helib/DoubleCRT.h>

#include "random_bit.hpp"
#include "tool.hpp"

#include "FHEtool_Yus_p.hpp"
#include "Yus_p.hpp"
#include "params_Yus_p.hpp"

namespace Server_Yus_p_C32
{
    bool Server_offline();
    bool Server_online();
}

#endif // CLIENT_HPP