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
#include <iomanip>
#include "randutils.hpp"
#include "processing.h"
#include <gsl/gsl_histogram.h>

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


template <typename T> 
double histVariance(std::vector<T>& array, const int bins=50){
    gsl_histogram * h = gsl_histogram_alloc (bins);
    T maxValArr = *std::max_element(array.begin(), array.end());
    T minValArr = *std::min_element(array.begin(), array.end());
    gsl_histogram_set_ranges_uniform (h, minValArr, maxValArr);
    for (int i = 0; i < n; i++){
        gsl_histogram_increment(h, array[i]);
    }
    // calculate standard deviation
    std::vector<double> binValues(bins);
    for (int i = 0; i < bins; i++){
        binValues[i] = gsl_histogram_get(h, i);
    }
    double mean = 0.0, variance = 0.0;
    for (int i = 0; i < bins; i++){
        mean += binValues[i];
    }
    mean /= bins;
    for (int i = 0; i < bins; i++){
        variance += (binValues[i] - mean) * (binValues[i] - mean);
    }
    variance /= bins;
    gsl_histogram_free (h);
    return variance;
}


double magnetisation(double t, double gamma, double m0=0.992){
    return m0*std::exp(-(1-gamma)*t);
}

int main(int argc, char **argv){

    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h")){
        // To do: write help
    }


    double T;
    const std::string &Tstring = input.getCmdOption("-T");
    if (!Tstring.empty()){
        T = std::stod(Tstring);
        // std::cout << "temperature T = " << T << "\n";
    }else{
        std::cerr << "temperature not provided" << std::endl;
        exit(0);
    }

    double m0;
    const std::string &m0string = input.getCmdOption("-m0");
    if (!m0string.empty()){
        m0 = std::stod(m0string);
        // std::cout << "reset mag. m0 = " << m0 << "\n";
    }else{
        std::cerr << "reset magnetisation not provided" << std::endl;
        exit(0);
    }

    double r;
    const std::string &rstring = input.getCmdOption("-r");
    if (!rstring.empty()){
        r = std::stod(rstring);
        // std::cout << "reset rate r = " << r << "\n";
    }else{
        std::cerr << "resetting rate not provided" << std::endl;
        exit(0);
    }

    static uint64_t seed;
    const std::string &seedstring = input.getCmdOption("-s");
    if (!seedstring.empty()){
        seed = std::stoull(seedstring);
    }else{
        seed = 328575958951598690;
    }

    randutils::seed_seq_fe128 seeder{uint32_t(seed),uint32_t(seed >> 32)};
    std::mt19937 mt19937Engine(seeder);

    double eta, eta0, gamma;

    std::exponential_distribution<> expDis(r);
    std::vector<double> energyArray(maxVal*stepSize);//, magArray(maxVal);
    std::vector<double> statisticsEnergy(n);//, statisticsMag(n);
    std::ofstream rawFile; rawFile.open("rawFile"+std::to_string(r), std::ios_base::app);
    std::ofstream histFile; histFile.open("histFile"+std::to_string(r), std::ios_base::app);
    eta0=m0*m0; // infinite temperature
    // eta0=1.0; // 0 temperature
    eta=std::tanh(1.0/T); gamma=std::tanh(2.0/T);
    for (int i = 0; i < maxVal*stepSize; i++){
        double x = static_cast<double>(i)/stepSize;
        energyArray[i] = energy(x, eta, eta0, gamma, 100);
        // energyArray[i] = magnetisation(x, gamma, 0.992);
    }
    printf("done energy calc T=%f\n", T);
    for (size_t i = 0; i < n; i++){
        statisticsEnergy[i] = expDis(mt19937Engine);
    }
    for (int i = 0; i < n; i++){
        statisticsEnergy[i] = (statisticsEnergy[i] >= maxVal)? energy(statisticsEnergy[i], eta, eta0, gamma, 100)
        : energyArray[static_cast<int>(statisticsEnergy[i]*stepSize)];
    }
    
    // histogram and calculate standard deviation
    double std = std::sqrt(histVariance(statisticsEnergy, 100));
    histFile << std::fixed << std::setprecision(10) << std << "\n";

    // standard deviation without histogram
    rawFile << std::fixed << std::setprecision(10) << std::sqrt(variance(statisticsEnergy)) << "\n";
}