#ifndef CLIENT_YUX2_8_C1_HPP
#define CLIENT_YUX2_8_C1_HPP

#include <iostream>
#include <cstring>
#include <stdint.h>
#include <vector>
#include <array>
#include <random>
#include <climits>
#include <fstream>
#include <chrono>

#include "random_bit.hpp"
#include "tool.hpp"
#include "Yux2_8.hpp"

#include "params_Yux2_8.hpp"
#include "FHEtool.hpp"

namespace C1ient_Yux2_8_C1{
    bool Client_offline();
    bool Client_online();
}

#endif // CLIENT_YUX2_8_C1_HPP