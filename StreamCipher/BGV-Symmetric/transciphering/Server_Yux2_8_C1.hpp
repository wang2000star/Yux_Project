#ifndef SERVER_YUX2_8_C1_HPP
#define SERVER_YUX2_8_C1_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include <stdint.h>
#include <cstring>
#include <cstdint>
#include <array>
#include <fstream>
#include <chrono>

#include "random_bit.hpp"
#include "tool.hpp"
#include "Yux2_8.hpp"

#include "params_Yux2_8.hpp"
#include "FHEtool_Yux2_8.hpp"
#include "Client_Yux2_8_C1.hpp"

namespace Server_Yux2_8_C1{
    bool Server_offline();
    bool Server_online();
}
#endif // CLIENT_HPP