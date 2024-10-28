#ifndef SERVER_HPP
#define SERVER_HPP

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

#include "params.hpp"
#include "FHEtool.hpp"
#include "Client.hpp"


bool Server_offline();
bool Server_online();

#endif // CLIENT_HPP