#ifndef CLIENT_YUS_P_C32_HPP
#define CLIENT_YUX_P_C32_HPP



#include <iostream>
#include <cstring>
#include <stdint.h>
#include <chrono>
#include <vector>
#include <omp.h>
#include <array>
#include <cstdint>
#include <cstring>
#include <atomic>
#include <cmath>
#include <chrono>
#include <fstream>
#include <memory>

#include <helib/helib.h>
#include <helib/ArgMap.h>
#include <helib/DoubleCRT.h>

#include "random_bit.hpp"
#include "tool.hpp"

#include "params_Yus_p.hpp"
#include "FHEtool_Yus_p.hpp"
#include "Yus_p.hpp"

namespace C1ient_Yus_p_C32{
    bool Client_offline();
    bool Client_online();
}

#endif // CLIENT_YUX_P_HPP