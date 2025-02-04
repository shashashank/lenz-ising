// 1D Ising Resetting for energy

#include"ising.h"
#include<iostream>
#include<fstream>
#include<string>
#include<chrono>
#include <iomanip>
// #include<gsl_randist.h>

int main(int argc, char *argv[]){
    std::ofstream data;
    data.open("energy.txt");
    const int L = 100; // L**2 = 10000
    isingLattice lattice(L);
    double T = 3.5, beta = 1.0/T;
    double m0 = 0.992;
    double r = 3.0;
    std::exponential_distribution<> expDis(r);
    int totalSteps = 0;
    while (totalSteps<300000){

        // resetting/initialising the lattice
        lattice.initialise(0.0);
        for (size_t i = 0; i < 1000; i++){
            lattice.glauber1DimSweep(beta);
        }

        // iters from exp dist scaled to lattice size and spin flips attempted
        int stepsItr =  static_cast<int>(expDis(mt19937Engine)*L*L);
        lattice.glauber1dInterval(beta, stepsItr);
        data << std::fixed<<std::setprecision(10)<<lattice.energy1D()<<"\n";
        totalSteps+=stepsItr/(L*L);
    }
    data.close();
    
}