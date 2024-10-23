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

#include "params.hpp"
#include "tool.hpp"
#include "random_bit.hpp"
#include "FHEtool.hpp"
#include "Client.hpp"
#include "Yux2_8.hpp"

bool Server_offline();
bool Server_online();

#endif // CLIENT_HPP