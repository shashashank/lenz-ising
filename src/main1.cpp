// This file creates an executable for simulating the ising model for statistics without resetting

#include"ising.h"
#include<iostream>
#include<fstream>
#include<string>
#include<chrono>
#include <iomanip>

int main(int argc, char *argv[]){
    std::ofstream data;
    const int thmTime=10000, simTime = 100, steps=10;
    std::vector<double> energy(steps*simTime), mag(steps*simTime);
    // std::vector<double> times(301); int count;
    double specHeat, sus;
    typedef std::chrono::high_resolution_clock Clock;
    using FpMilliseconds = std::chrono::duration<float, std::chrono::milliseconds::period>;
    const int L = 100; // L**2 = 10.000
    double betaInitial=1.0/10.0, betaFinal = 1.0/1.5;
    for (int j =0; j<100; j++){
        isingLattice lattice(L);

        // initialization
        lattice.initialise(0);

        //thermalization
        for (int k = 0; k < thmTime; k++)
        {
            // lattice.kawasaki1DimSweep(beta);
            lattice.glauber1DimSweep(betaInitial);
        }

        // dynamics and data collection
        for (int k = 0; k < simTime; k++)
        {
            for (size_t i = 0; i < steps; i++)
            {
                energy[k*steps+i] = lattice.energy1D();
                mag[k*steps+i] = lattice.magnetisation();
                lattice.glauber1dInterval(betaFinal, L*L/steps);
            }
            // energy[k] = lattice.energy1D();
            // mag[k] = lattice.magnetisation();
            // lattice.glauber1DimSweep(betaFinal);
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