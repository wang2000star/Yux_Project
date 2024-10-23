#ifndef CLIENT_HPP
#define CLIENT_HPP

#include <iostream>
#include <cstring>
#include <stdint.h>
#include <vector>
#include <array>
#include <random>
#include <climits>
#include <fstream>

#include "random_bit.hpp"
#include "tool.hpp"
#include "Yux2_8.hpp"
#include "params.hpp"
#include "FHEtool.hpp"

bool Client_offline();
bool Client_online();
bool Client_Check();
#endif // CLIENT_HPP