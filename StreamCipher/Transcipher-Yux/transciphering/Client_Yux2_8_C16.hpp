#ifndef CLIENT_YUX2_8_C16_HPP
#define CLIENT_YUX2_8_C16_HPP

#include <iostream>
#include <cstring>
#include <stdint.h>
#include <vector>
#include <array>
#include <random>
#include <climits>
#include <fstream>
#include <chrono>
#include <omp.h> // 包含 OpenMP 头文件，OpenMP 是一种并行编程模型

#include "random_bit.hpp"
#include "tool.hpp"
#include "Yux2_8.hpp"

#include "params_Yux2_8.hpp"
#include "FHEtool.hpp"

namespace C1ient_Yux2_8_C16{
    bool Client_offline();
    bool Client_online();
}

#endif // CLIENT_YUX2_8_C16_HPP