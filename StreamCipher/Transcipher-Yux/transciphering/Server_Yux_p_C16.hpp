#ifndef SERVER_YUX_P_C16_HPP
#define SERVER_YUX_P_C16_HPP


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

#include "FHEtool.hpp"
#include "random_bit.hpp"
#include "tool.hpp"
#include "Yux_p.hpp"

#include "params_Yux_p.hpp"

namespace Server_Yux_p_C16
{
    bool Server_offline();
    bool Server_online();
}

#endif // CLIENT_HPP