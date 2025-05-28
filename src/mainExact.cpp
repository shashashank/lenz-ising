// This file creates an executable for getting resetting statistics without simulation
#include <cmath>
#include <iostream>
#include <array>
#include <vector>
#include <fstream>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <omp.h>
#include <iomanip>
#include "randutils.hpp"
#include "processing.h"
#include <gsl/gsl_histogram.h>


static uint64_t seed = 328575958951598690;
randutils::seed_seq_fe128 seeder{uint32_t(seed),uint32_t(seed >> 32)};
std::mt19937 mt19937Engine(seeder);
const int n = 10000000, stepSize = 10000, maxVal = 100;

double energy(double t, double eta, double eta0, double gamma, int terms=100){
    double total_sum = 0.0, eta0Mult = 1.0, etaMult = 1.0;
    for (int l = 1; l < terms + 1; l++){
        eta0Mult *= eta0; etaMult *= eta;
        double term1 = (eta0Mult-etaMult)*
        (std::cyl_bessel_i(l-1,2.0*gamma*t)-std::cyl_bessel_i(1+l,2.0*gamma*t));
        total_sum += term1;
    }
    return -(eta + std::exp(-2 * t) * total_sum);
}

double magnetisation(double t, double gamma, double m0=0.992){
    return m0*std::exp(-(1-gamma)*t);
}

template <typename T> 
void binArray(std::vector<T>& array, FILE* stream ,int bins=50){
    gsl_histogram * h = gsl_histogram_alloc (bins);
    T maxValArr = *std::max_element(array.begin(), array.end());
    T minValArr = *std::min_element(array.begin(), array.end());
    gsl_histogram_set_ranges_uniform (h, minValArr, maxValArr);
    for (size_t i = 0; i < n; i++){
        gsl_histogram_increment(h, array[i]);
    }
    gsl_histogram_fprintf (stream, h, "%g", "%g");
    gsl_histogram_free (h);
}

template<typename T>
void writeToFile(const std::vector<T>& array, FILE* stream){

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < array.size(); i++){
        fprintf(stream, "%f\n", array[i]);
    }
}

int main(void){
    double eta, eta0, gamma;
    std::vector<double> energyArray(maxVal*stepSize), magArray(maxVal);
    std::vector<double> statisticsEnergy(n);
    std::ofstream rawFile; rawFile.open("rawFile");
    // FILE * histFile = fopen("histFile.txt", "w");
    eta0=0.984064; // infinite temperature
    // eta0=1.0; // 0 temperature
    typedef std::chrono::high_resolution_clock Clock;
    auto t0=Clock::now();
    omp_set_num_threads(6);
    omp_set_dynamic(0);
    for(double T = 0.5; T < 5.5; T+=0.5){
        eta=std::tanh(1.0/T); gamma=std::tanh(2.0/T);
#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < maxVal*stepSize; i++){
            double x = static_cast<double>(i)/stepSize;
            energyArray[i] = energy(x, eta, eta0, gamma, 100);
            // energyArray[i] = magnetisation(x, gamma, 0.992);
            // magArray[i] = magnetisation(x, gamma, 0.992);
        }
        printf("done energy calc T=%f\n", T);
        for(double r = 0.1; r < 2.0; r+=0.01){
            std::exponential_distribution<> expDis(r);
            for (size_t i = 0; i < n; i++){
                statisticsEnergy[i] = expDis(mt19937Engine);
            }
#pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < n; i++){
                statisticsEnergy[i] = (statisticsEnergy[i] >= maxVal)? energy(statisticsEnergy[i], eta, eta0, gamma, 100)
                : energyArray[static_cast<int>(statisticsEnergy[i]*stepSize)];
                // statisticsEnergy[i] = energy(statisticsEnergy[i], eta, eta0, gamma, 100);                
                // statisticsEnergy[i] = (statisticsEnergy[i] >= maxVal)? magnetisation(statisticsEnergy[i], gamma, 0.992)
                // : energyArray[static_cast<int>(statisticsEnergy[i]*stepSize)];
            }
            // binArray(statisticsEnergy, histFile, 100);
            // calculate standard deviation
            rawFile << std::fixed << std::setprecision(10) << std::sqrt(variance(statisticsEnergy)) << "\t";
        }
        rawFile << std::endl;
    }
    rawFile.close();
    // fclose(histFile);
    auto t1=Clock::now();
    std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count() << " milliseconds\n";
    std::cout << std::endl;
}