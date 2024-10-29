#ifndef CLIENT_YUX_P_HPP
#define CLIENT_YUX_P_HPP



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

#include "params_Yux_p.hpp"

#include "FHEtool_Yux_p.hpp"
#include "random_bit.hpp"
#include "tool.hpp"
#include "Yux_p.hpp"

namespace C1ient_Yux_p_C16{
    bool Client_offline();
    bool Client_online();
}

#endif // CLIENT_YUX_P_HPP