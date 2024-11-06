#ifndef RANDOM_BIT_HPP
#define RANDOM_BIT_HPP

#include <vector>
#include <bitset>
#include <string>
#include <cstdint>

typedef long tKeccakLane;

template <unsigned BlockSize>
class RandomBit
{
public:
    // 构造函数
    RandomBit(unsigned rounds);

    // 类型定义
    typedef std::bitset<BlockSize> block;
    typedef std::bitset<BlockSize> keyblock;

    // 成员函数声明
    unsigned rank_of_Matrix(const std::vector<block> &matrix);
    unsigned LFSR86540(unsigned char *LFSR);
    void KeccakF1600_StatePermute(void *state);
    bool getrandbit();
    block getrandblock();
    void initrand(const long nonce, const long counter);
    void generate_Instance_all_new(const long nonce, const long counter);

    // 成员变量
    unsigned rounds; // Number of rounds
    std::vector<std::vector<block>> LinMatrices;
    std::vector<block> roundconstants;
    unsigned squeezedbits;
    unsigned char shakestate[200];
    const unsigned shakerate = 1088;
};

// 构造函数实现
template <unsigned BlockSize>
RandomBit<BlockSize>::RandomBit(unsigned r) : rounds(r)
{
    // 初始化其他成员变量
    squeezedbits = 0;
    std::fill(std::begin(shakestate), std::end(shakestate), 0);
}

// 成员函数实现
template <unsigned BlockSize>
unsigned RandomBit<BlockSize>::rank_of_Matrix(const std::vector<block> &matrix)
{
    std::vector<block> mat; // Copy of the matrix
    for (auto u : matrix)
    {
        mat.push_back(u);
    }
    unsigned size = mat[0].size();
    // Transform to upper triangular matrix
    unsigned row = 0;
    for (unsigned col = 1; col <= size; ++col)
    {
        if (!mat[row][size - col])
        {
            unsigned r = row;
            while (r < mat.size() && !mat[r][size - col])
            {
                ++r;
            }
            if (r >= mat.size())
            {
                continue;
            }
            else
            {
                auto temp = mat[row];
                mat[row] = mat[r];
                mat[r] = temp;
            }
        }
        for (unsigned i = row + 1; i < mat.size(); ++i)
        {
            if (mat[i][size - col])
                mat[i] ^= mat[row];
        }
        ++row;
        if (row == size)
            break;
    }
    return row;
}

// Shake taken from Keccak Code package
#define ROL64(a, offset) ((((long)a) << offset) ^ (((long)a) >> (64 - offset)))
#define fun(x, y) ((x) + 5 * (y))
#define readLane(x, y) (((tKeccakLane *)state)[fun(x, y)])
#define writeLane(x, y, lane) (((tKeccakLane *)state)[fun(x, y)]) = (lane)
#define XORLane(x, y, lane) (((tKeccakLane *)state)[fun(x, y)]) ^= (lane)

template <unsigned BlockSize>
unsigned RandomBit<BlockSize>::LFSR86540(unsigned char *LFSR)
{
    unsigned result = ((*LFSR) & 0x01) != 0;
    if (((*LFSR) & 0x80) != 0)
        /* Primitive polynomial over GF(2): x^8+x^6+x^5+x^4+1 */
        (*LFSR) = ((*LFSR) << 1) ^ 0x71;
    else
        (*LFSR) <<= 1;
    return result;
}

template <unsigned BlockSize>
void RandomBit<BlockSize>::KeccakF1600_StatePermute(void *state)
{
    unsigned round, x, y, j, t;
    unsigned char LFSRstate = 0x01;

    for (round = 0; round < 24; round++)
    {
        { /* === θ step (see [Keccak Reference, Section 2.3.2]) === */
            tKeccakLane C[5], D;
            /* Compute the parity of the columns */
            for (x = 0; x < 5; x++)
            {
                C[x] = readLane(x, 0) ^ readLane(x, 1) ^ readLane(x, 2) ^ readLane(x, 3) ^ readLane(x, 4);
            }
            for (x = 0; x < 5; x++)
            {
                /* Compute the θ effect for a given column */
                D = C[(x + 4) % 5] ^ ROL64(C[(x + 1) % 5], 1);
                /* Add the θ effect to the whole column */
                for (y = 0; y < 5; y++)
                {
                    XORLane(x, y, D);
                }
            }
        }

        { /* === ρ and π steps (see [Keccak Reference, Sections 2.3.3 and 2.3.4]) === */
            tKeccakLane current, temp;
            /* Start at coordinates (1 0) */
            x = 1;
            y = 0;
            current = readLane(x, y);
            /* Iterate over ((0 1)(2 3))^t * (1 0) for 0 ≤ t ≤ 23 */
            for (t = 0; t < 24; t++)
            {
                /* Compute the rotation constant r = (t+1)(t+2)/2 */
                unsigned r = ((t + 1) * (t + 2) / 2) % 64;
                /* Compute ((0 1)(2 3)) * (x y) */
                unsigned Y = (2 * x + 3 * y) % 5;
                x = y;
                y = Y;
                /* Swap current and state(x,y), and rotate */
                temp = readLane(x, y);
                writeLane(x, y, ROL64(current, r));
                current = temp;
            }
        }

        { /* === χ step (see [Keccak Reference, Section 2.3.1]) === */
            tKeccakLane temp[5];
            for (y = 0; y < 5; y++)
            {
                /* Take a copy of the plane */
                for (x = 0; x < 5; x++)
                    temp[x] = readLane(x, y);
                /* Compute χ on the plane */
                for (x = 0; x < 5; x++)
                    writeLane(x, y, temp[x] ^ ((~temp[(x + 1) % 5]) & temp[(x + 2) % 5]));
            }
        }

        { /* === ι step (see [Keccak Reference, Section 2.3.5]) === */
            for (j = 0; j < 7; j++)
            {
                unsigned bitPosition = (1 << j) - 1; /* 2^j-1 */
                if (LFSR86540(&LFSRstate))
                {
                    XORLane(0, 0, (tKeccakLane)1 << bitPosition);
                }
            }
        }
    }
}
// End of Keccak Code package

template <unsigned BlockSize>
bool RandomBit<BlockSize>::getrandbit()
{
    bool tmp = 0;
    tmp = (bool)((shakestate[squeezedbits / 8] >> (squeezedbits % 8)) & 1);
    squeezedbits++;
    if (squeezedbits == shakerate)
    {
        KeccakF1600_StatePermute(shakestate);
        squeezedbits = 0;
    }
    return tmp;
}

template <unsigned BlockSize>
typename RandomBit<BlockSize>::block RandomBit<BlockSize>::getrandblock()
{
    block tmp = 0;
    for (unsigned i = 0; i < BlockSize; ++i)
        tmp[i] = getrandbit();
    return tmp;
}

// Shake like rand functions
template <unsigned BlockSize>
void RandomBit<BlockSize>::initrand(const long nonce, const long counter)
{

    for (unsigned i = 0; i < 200; ++i)
        shakestate[i] = 0;

    for (unsigned i = 0; i < 8; i++)
        shakestate[i] = (rounds >> (8 * i)) & 0xff;

    for (unsigned i = 8; i < 16; i++)
        shakestate[i] = (BlockSize >> (8 * i)) & 0xff;

    for (unsigned i = 24; i < 32; i++)
        shakestate[i] = (nonce >> (8 * i)) & 0xff;

    for (unsigned i = 40; i < 48; i++)
        shakestate[i] = (counter >> (8 * i)) & 0xff;

    shakestate[32] = 0x1F;

    shakestate[shakerate / 8 - 1] ^= 0x80;

    KeccakF1600_StatePermute(shakestate);
    squeezedbits = 0;
}

template <unsigned BlockSize>
void RandomBit<BlockSize>::generate_Instance_all_new(const long nonce, const long counter)
{
    // Initialize RNG state
    initrand(nonce, counter);
    // Create LinMatrices and invLinMatrices
    LinMatrices.clear();
    for (unsigned r = 0; r <= rounds; ++r)
    {
        // Create matrix
        std::vector<block> mat;
        // Fill matrix with random bits
        do
        {
            mat.clear();
            for (unsigned i = 0; i < BlockSize; ++i)
            {
                mat.push_back(getrandblock());
            }
            // Repeat if matrix is not invertible
        } while (rank_of_Matrix(mat) != BlockSize);
        LinMatrices.push_back(mat);
    }

    // Create roundconstants
    roundconstants.clear();
    for (unsigned r = 0; r <= rounds; ++r)
    {
        roundconstants.push_back(getrandblock());
    }

    return;
}

#endif // RANDOM_BIT_HPP