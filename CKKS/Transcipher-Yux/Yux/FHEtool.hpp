#ifndef FHETOOL_HPP
#define FHETOOL_HPP

#include <iostream>
#include <cstring>
#include <stdint.h>

#include <NTL/ZZX.h>
#include <NTL/GF2X.h>
#include <helib/helib.h>
#include <helib/ArgMap.h>
#include <helib/DoubleCRT.h>

#include "params.hpp"

void encodeTo1Ctxt(NTL::Vec<NTL::ZZX>& encData, const NTL::Vec<uint8_t>& data,
                   const helib::EncryptedArrayDerived<helib::PA_GF2>& ea);

void decodeTo1Ctxt(NTL::Vec<uint8_t>& data, const NTL::Vec<NTL::ZZX>& encData,
                   const helib::EncryptedArrayDerived<helib::PA_GF2>& ea);

#endif // FHETOOL_HPP