
#include <random>
#include"randutils.hpp"

static uint64_t seed = 328575958951598690;
randutils::seed_seq_fe128 seeder{uint32_t(seed),uint32_t(seed >> 32)};
std::mt19937 mt19937Engine(seeder);
std::uniform_real_distribution<> rDist(0.0, 1.0);

inline int mod(int a, int b)
{
    int ret = a % b;
    if (ret<0)
        ret+=b;
    return ret;
}

template<typename T> double meanOf(const T& array){
    double mean = 0;
    for (auto& val: array){
        mean+=val;
    }
    return mean/((double)array.size());
}

template<typename T> double variance(const T& array){
    double mean = meanOf(array);
    double variance = 0.0;
    for(auto& val : array){
        auto diff = val - mean;
        variance += diff * diff;
    }
    return variance/((double)array.size());
}

double sample_exponential(double lambda) {
    double u = rDist(mt19937Engine);  // Generate a uniform random number
    return -std::log(1 - u) / lambda;  // Inverse transform sampling
}