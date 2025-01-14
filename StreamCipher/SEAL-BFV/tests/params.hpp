#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <vector>
#include <iostream>

void set_params(long poly_modulus_degree, int SecLevel, std::vector<int> &bit_sizes)
{
    if (poly_modulus_degree == 16384)
    {
        if (SecLevel == 128)
        {
        bit_sizes = {53, 53, 53, 53, 53, 53, 53, 53};
        // logq = 424, maximum depth = 10
        }
        if (SecLevel == 192)
        {
        bit_sizes = {35,35,35,35,35,35,35,35};
        }
    }
    if (poly_modulus_degree == 32768)
    {
        if (SecLevel = 192){
        bit_sizes = {42, 41, 41, 41, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42};
        // logq = 585, maximum depth = 14
        }
    }
    if (poly_modulus_degree == 65536)
    {
        SecLevel = 256;
        bit_sizes = {44, 43, 43, 43, 43, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44};
        // logq = 920, maximum depth = 23
    }
}
#endif