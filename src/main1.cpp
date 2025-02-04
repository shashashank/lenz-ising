// This file creates an executable for simulating the ising model for statistics without resetting

#include"ising.h"
#include<iostream>
#include<fstream>
#include<string>
#include<chrono>
#include <iomanip>

int main(int argc, char *argv[]){

    const int thmTime=50000, simTime = 100, steps=10;
    const int L = 100; // N = L**2 = 10000
    const int nTrials = 100;
    double betaInitial=1.0, betaFinal = 1.0/1.1;

    std::ofstream data;    
    std::vector<double> energy(steps*simTime), mag(steps*simTime);
    // std::vector<double> times(301); int count;
    double specHeat, sus;
    typedef std::chrono::high_resolution_clock Clock;
    using FpMilliseconds = std::chrono::duration<float, std::chrono::milliseconds::period>;
    for (int j =0; j<nTrials; j++){
        isingLattice lattice(L);

        // initialization
        lattice.initialise(0);

        //thermalization
        for (int k = 0; k < thmTime; k++){
            lattice.glauber1DimSweep(betaInitial);
        }

        // dynamics and data collection
        for (int k = 0; k < simTime; k++){
            for (size_t i = 0; i < steps; i++){
                energy[k*steps+i] = lattice.energy1D();
                mag[k*steps+i] = lattice.magnetisation();
                lattice.glauber1dInterval(betaFinal, L*L/steps);
            }
        }
        // writing to disk
        data.open("data"+std::to_string(j));
        data<< "#MAGNETISATION\tENERGY\n";
        for (size_t i = 0; i < steps*simTime; i++)
        {
            data << std::fixed<<std::setprecision(10)<<mag[i]<<"\t";
            data << std::fixed<<std::setprecision(10)<<energy[i]<<"\n";
        }
        data.close();
        std::cout << "Fininshed with " << j+1 << " trials" <<std::endl;
    }
}