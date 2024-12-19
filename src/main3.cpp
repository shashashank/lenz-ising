#include"ising.h"
#include<iostream>
#include<fstream>
#include<string>
#include<chrono>
#include <iomanip>
// #include<gsl_randist.h>

int main(int argc, char *argv[]){
    std::ofstream data, config, times;
    data.open("data.txt");
    times.open("times.txt");
    // config.open("config.txt");
    const int L = 100; // L**2 = 10.000
    isingLattice lattice(L);
    double T = 3.5, beta = 1.0/T;
    double m0 = 0.992;
    double r = 0.1175;
    std::exponential_distribution<> expDis(r);
    // data<<"#MAGNETISATION\tENERGY\n";
    for (size_t i = 0; i < 10000; i++)
    {
        // steps = gsl_ran_exponential(mt19937Engine,1/r);
        int stepsItr =  static_cast<int>(expDis(mt19937Engine));
        lattice.initialise(m0);
        for (size_t j = 0; j < stepsItr; j++)
        {
            lattice.glauber1DimSweep(beta);
        }
        times << stepsItr << "\n";
        data << lattice.magnetisation() << "\n";
    }
    data.close();
    times.close();
}