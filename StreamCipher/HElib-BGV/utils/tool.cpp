#include "tool.hpp"

// 生成m比特强度的随机整数
long generate_secure_random_int(unsigned m) {
    if (m > 64) {
        throw std::invalid_argument("m exceeds the maximum bit size of long");
    }

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<long> dis(0, (1ULL << m) - 1);

    return dis(gen);
}
void random_init_shake(long nonce, long block_counter, Keccak_HashInstance &shake128_)
{
    uint8_t seed[16];

    *((long *)seed) = htobe64(nonce);
    *((long *)(seed + 8)) = htobe64(block_counter);

    if (SUCCESS != Keccak_HashInitialize_SHAKE128(&shake128_))
        throw std::runtime_error("failed to init shake");
    if (SUCCESS != Keccak_HashUpdate(&shake128_, seed, sizeof(seed) * 8))
        throw std::runtime_error("SHAKE128 update failed");
    if (SUCCESS != Keccak_HashFinal(&shake128_, NULL))
        throw std::runtime_error("SHAKE128 final failed");
}
long generate_random_field_element(Keccak_HashInstance &shake128, bool allow_zero, long max_prime_size, long PlainMod)
{
    uint8_t random_bytes[sizeof(long)];
    while (1)
    {
        if (SUCCESS !=
            Keccak_HashSqueeze(&shake128, random_bytes, sizeof(random_bytes) * 8))
            throw std::runtime_error("SHAKE128 squeeze failed");
        long ele = be64toh(*((long *)random_bytes)) & max_prime_size;
        if (!allow_zero && ele == 0)
            continue;
        if (ele < PlainMod)
            return ele;
    }
}

int min_noise_budget(std::vector<helib::Ctxt> &eData)
{
    int min_noise = 10000;
    int noise;
    for (int i = 0; i < eData.size(); i++)
    {
        noise = eData[i].bitCapacity();
        if (noise < min_noise)
        {
            min_noise = noise;
        }
    }
    return min_noise;
}
// 函数：对多项式的每个系数乘以整数 a 并取模 c
helib::zzX multiplyAndMod(const helib::zzX &a, long b,long pmod)
{
    helib::zzX res;
    res.SetLength(helib::lsize(a));
    int len = helib::lsize(a);
    for (long i = 0; i < len; ++i)
    {
        res[i] = (a[i] * b) % pmod;
    }
    return res;
}

bool writeEncryptedSymKey(const std::vector<helib::Ctxt> &encryptedSymKey, const std::string &filename)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open())
    {
        std::cerr << "Failed to open " << filename << " for writing" << std::endl;
        return false;
    }

    for (const auto &ctxt : encryptedSymKey)
    {
        ctxt.writeTo(out);
    }

    out.close();

    return true;
}

void decodeToCtxt(std::vector<long> &data, const std::vector<NTL::vec_long> &encData, const long CtxtWords,const long nslots)
{
    long R = encData.size() / CtxtWords;
    long AllByte = CtxtWords * nslots;
    long data_size = R * AllByte;
    data.resize(data_size);
    long byteIdx;
    long rAB;
    long jCB;
    long rCB;
    //    omp_set_num_threads(16); // 设置线程数为16
    // #pragma omp parallel for
    for (long j = 0; j < nslots; j++)
    {
        jCB = j * CtxtWords;
        for (long r = 0; r < R; r++)
        {
            rAB = r * AllByte;
            rCB = r * CtxtWords;
            for (long i = 0; i < CtxtWords; i++)
            { // i is the ciphertext number
                // j is the block number in this ctxt
                byteIdx = jCB + i + rAB;
                data[byteIdx] = encData[rCB + i][j];
            }
        }
    }
}

// 函数：解密并验证密文是否正确，需要解码
bool verifyDecryption(const std::vector<helib::Ctxt> &encryptedVec, const std::vector<long> &originalVec, const helib::SecKey &secretKey,
                      const helib::Cmodulus &cmodulus, const long CtxtWords, const long PlainBlock, const long nslots, const long pmod)
{
    std::vector tempVec = originalVec;
    for (int i = 0; i < originalVec.size(); i++)
    {
        tempVec[i] = (tempVec[i] + pmod) % pmod;
    }
    int size = encryptedVec.size();
    std::cout << "size: " << size << std::endl;
    auto start_decrypt = std::chrono::high_resolution_clock::now();
    std::vector<long> decryptedVec = originalVec;
    std::vector<NTL::vec_long> decryptedPolys(size);
    std::vector<NTL::ZZX> Polys(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        secretKey.Decrypt(Polys[i], encryptedVec[i]);
        cmodulus.FFT(decryptedPolys[i], Polys[i]);
    }
    decodeToCtxt(decryptedVec, decryptedPolys, CtxtWords, nslots);
    // 验证解密结果
    bool isDecryptedVecCorrect = std::equal(decryptedVec.begin(), decryptedVec.begin()+PlainBlock*CtxtWords, tempVec.begin());
    auto end_decrypt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification finished! Time: " << elapsed_seconds.count() << "s\n";
    // 如果解密结果不正确，输出错误的位置
    if (!isDecryptedVecCorrect)
    {
        for (size_t i = 0; i < CtxtWords; i++)
        {
            if (decryptedVec[i] != tempVec[i])
            {
                std::cout << "Error at position " << i << ": " << decryptedVec[i] << " != " << originalVec[i] << std::endl;
                // break;
            }
        }
    }
    return isDecryptedVecCorrect;
}

void encryptSymKey(std::vector<helib::Ctxt> &encryptedSymKey, const std::vector<long> &SymKey, std::unique_ptr<helib::PubKey> &pk, const helib::Cmodulus &cmodulus,const long nslots)
{
    NTL::vec_long slotsData;
    slotsData.SetLength(nslots);
    long BlockWords = SymKey.size();
    encryptedSymKey.resize(BlockWords, helib::Ctxt(*pk));
    NTL::zz_pX temp;
    NTL::ZZX encodedData;
    for (long i = 0; i < BlockWords; i++)
    { // encrypt the encoded key
        for (long j = 0; j < nslots; j++)
        {
            slotsData[j] = SymKey[i];
        }
        cmodulus.iFFT(temp, slotsData);
        conv(encodedData, temp);
        pk->Encrypt(encryptedSymKey[i], encodedData);
    }
}

bool verify_encryptSymKey(std::vector<helib::Ctxt> &encryptedSymKey, const std::vector<long> &SymKey, const helib::SecKey &secretKey,
                          const helib::Cmodulus &cmodulus)
{
    auto start_decrypt = std::chrono::high_resolution_clock::now();
    long BlockWords = SymKey.size();
    std::vector<long> decryptedSymKey = SymKey;
    NTL::ZZX encodedSymKey;
    NTL::vec_long slotsData;
    //     omp_set_num_threads(16); // 设置线程数为16
    // #pragma omp parallel for
    for (long i = 0; i < BlockWords; i++)
    { // encrypt the encoded key
        secretKey.Decrypt(encodedSymKey, encryptedSymKey[i]);
        cmodulus.FFT(slotsData, encodedSymKey);
        decryptedSymKey[i] = slotsData[0];
    }
    bool isDecryptedSymKeyCorrect = std::equal(SymKey.begin(), SymKey.end(), decryptedSymKey.begin());
    auto end_decrypt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_decrypt - start_decrypt;
    std::cout << "Decryption and verification finished! Time: " << elapsed_seconds.count() << "s\n";
    return isDecryptedSymKeyCorrect;
}
std::string get_cpu_model()
{
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line))
    {
        if (line.find("model name") != std::string::npos)
        {
            return line.substr(line.find(":") + 2);
        }
    }
    return "Unknown CPU model";
}
// std::string get_memory_info() {
//     std::ifstream meminfo("/proc/meminfo");
//     std::string line;
//     std::ostringstream meminfo_stream;
//     while (std::getline(meminfo, line)) {
//         if (line.find("MemTotal") != std::string::npos || line.find("MemFree") != std::string::npos) {
//             meminfo_stream << line << std::endl;
//         }
//     }
//     return meminfo_stream.str();
// }

// std::string get_memory_info() {
//     std::ifstream meminfo("/proc/meminfo");
//     std::string line;
//     std::ostringstream meminfo_stream;
//     long mem_total = 0;
//     long swap_total = 0;

//     while (std::getline(meminfo, line)) {
//         if (line.find("MemTotal") != std::string::npos) {
//             mem_total = std::stol(line.substr(line.find(":") + 1));
//         }
//         if (line.find("SwapTotal") != std::string::npos) {
//             swap_total = std::stol(line.substr(line.find(":") + 1));
//         }
//     }

//     meminfo_stream << "MemTotal: " << mem_total << " kB" << std::endl;
//     meminfo_stream << "SwapTotal: " << swap_total << " kB" << std::endl;

//     return meminfo_stream.str();
// }
std::string get_memory_info()
{
    std::ifstream meminfo("/proc/meminfo");
    std::string line;
    std::ostringstream meminfo_stream;
    long mem_total_kb = 0;
    long swap_total_kb = 0;

    while (std::getline(meminfo, line))
    {
        if (line.find("MemTotal") != std::string::npos)
        {
            mem_total_kb = std::stol(line.substr(line.find(":") + 1));
        }
        if (line.find("SwapTotal") != std::string::npos)
        {
            swap_total_kb = std::stol(line.substr(line.find(":") + 1));
        }
    }

    double mem_total_gib = mem_total_kb / 1024.0 / 1024.0;
    double swap_total_gib = swap_total_kb / 1024.0 / 1024.0;

    meminfo_stream << "MemTotal: " << mem_total_gib << " GiB" << std::endl;
    meminfo_stream << "SwapTotal: " << swap_total_gib << " GiB" << std::endl;

    return meminfo_stream.str();
}

std::string get_os_version()
{
    std::ifstream osrelease("/etc/os-release");
    std::string line;
    std::ostringstream osrelease_stream;
    while (std::getline(osrelease, line))
    {
        if (line.find("PRETTY_NAME") != std::string::npos)
        {
            osrelease_stream << line.substr(line.find("=") + 1) << std::endl;
        }
    }
    return osrelease_stream.str();
}

std::string to_lower(const std::string &str) {
    std::string lower_str = str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return lower_str;
}

std::string get_environment_info() {
    std::string environment = "Unknown";

    // Check for WSL
    std::ifstream version("/proc/version");
    std::string version_line;
    if (std::getline(version, version_line)) {
        std::cout << "Version line: " << version_line << std::endl; // 调试信息
        std::string lower_version_line = to_lower(version_line);
        if (lower_version_line.find("microsoft") != std::string::npos) {
            //std::cout << "Detected Microsoft in version line" << std::endl; // 调试信息
            if (lower_version_line.find("wsl2") != std::string::npos) {
                //std::cout << "Detected WSL2 in version line" << std::endl; // 调试信息
                return "WSL2";
            } else {
                return "WSL";
            }
        }
    }

    // Check for virtual machine
    std::ifstream dmi("/sys/class/dmi/id/product_name");
    std::string dmi_line;
    if (std::getline(dmi, dmi_line)) {
        //std::cout << "DMI line: " << dmi_line << std::endl; // 调试信息
        std::string lower_dmi_line = to_lower(dmi_line);
        if (lower_dmi_line.find("virtualbox") != std::string::npos || lower_dmi_line.find("vmware") != std::string::npos) {
            return "Virtual Machine";
        }
    }

    // If not WSL or VM, assume it's a regular Linux system
    return "Linux";
}