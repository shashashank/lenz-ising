// This file creates an executable for getting resetting statistics without simulation
#include <cmath>
#include <iostream>
#include <array>
#include <vector>
#include <fstream>
#include "randutils.hpp"


static uint64_t seed = 328575958951598690;
randutils::seed_seq_fe128 seeder{uint32_t(seed),uint32_t(seed >> 32)};
std::mt19937 mt19937Engine(seeder);
const int n = 1000000, stepSize = 1000, maxVal =10000;

double energy(double t, double eta, double eta0, double gamma, int terms=100){
    double exp_term = std::exp(-2 * t);

    double total_sum = 0.0;
    for (int l = 1; l < terms + 1; l++){
        double term1 = (std::pow(eta0, l) - std::pow(eta, l)) *
         (std::cyl_bessel_i(1-l, 2.0*gamma*t) - std::cyl_bessel_i(1+l, 2.0*gamma*t));
        total_sum += term1;
    }
    return -(eta + exp_term * total_sum);
}

double magnetisation(double t, double gamma, double m0=0.992){
    return m0*std::exp(-(1-gamma)*t);
}

void binArray(std::vector<int>& array, int bins=50){
    std::vector<int> hist(bins); 
    hist = {0};
    for (int i = 0; i < array.size(); i++){
        hist[array[i]] += 1;
    }
    array.resize(bins);
    array = hist;
}

int main(void){
    std::cout << "here" << std::endl;
    double eta, eta0, gamma;
    std::array<double, maxVal*stepSize> x, energyArray, magArray;
    printf("here\n");
    std::vector<int> time_between_calls1(n);
    printf("here\n");
    std::array<double, 100> hist;
    printf("here\n");
    std::ofstream histFile;
    histFile.open("histFile.txt");
    double valTmp = 0;

    printf("here\n");
    for (auto i: x)
    {
        i = valTmp;
        valTmp+=1/stepSize;
    }
    printf("here\n");

    eta0=std::tanh(0); // infinite temperature
    for(double T = 0.5; T < 5.5; T+=0.5)
    {
        eta=std::tanh(1/T); gamma=std::tanh(2/T);
        for (int i = 0; i < stepSize*maxVal; i++){
            energyArray[i] = energy(x[i], eta, eta0, gamma, 100);
            magArray[i] = magnetisation(x[i], gamma, 0.992);
        }
        for(double r = 0.5; r < 1; r+=0.01){
            std::exponential_distribution<> expDis(r);
            for (int i = 0; i < n; i++)
            {
                int tmp = static_cast<int>(expDis(mt19937Engine)*stepSize);
                if (tmp >= maxVal){
                    time_between_calls1[i] = energy(tmp, eta, eta0, gamma, 100);
                }else{
                    time_between_calls1[i] = energyArray[tmp];
                }
            }
            // binArray(time_between_calls1);
            for (int i = 0; i < 100; i++)
            {
                histFile << time_between_calls1[i] << " ";
            }
            histFile << "\n";
        }
    }
}