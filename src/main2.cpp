#include"ising.h"
#include<iostream>
#include<fstream>
#include<string>
#include<chrono>
#include <iomanip>

int main(int argc, char *argv[]){
    std::ofstream data;
    int thmTime = 10000, simTime = 100000;
    std::vector<double> energy(simTime), mag(simTime);
    std::vector<double> times(301); int count;
    double specHeat, sus;
    typedef std::chrono::high_resolution_clock Clock;
    using FpMilliseconds = std::chrono::duration<float, std::chrono::milliseconds::period>;
    for(int L = 4; L<=128; L*=2){
        isingLattice lattice(L*L);
        count = 0;
        data.open("data"+std::to_string(L));
        for (double beta = 0.3; beta <=0.6 ; beta+=0.001)
        {
            // initialization
            lattice.initialise(0.992);
            // thermalization
            auto t0=Clock::now();
            for (size_t k = 0; k < thmTime; k++)
            {
                lattice.glauber2DimSweep(beta);
            }
            auto t1=Clock::now();
            //data collection
            for (size_t k = 0; k < simTime; k++)
            {
                lattice.glauber2DimSweep(beta);
                energy[k] = lattice.energy2D();
                mag[k] = lattice.magnetisation();
            }
            // data writing
            data << "#SUSCEPTIBILITY\tSPECIFICHEAT\n";
            sus = beta*L*L*variance(mag);
            specHeat = beta*beta*L*L*variance(energy);
            data << beta << "\t" << sus << "\t" << specHeat << "\n";
            times[count++] = FpMilliseconds(t1-t0).count(); // time in ms
        }
        data.close();
        // std::cout << "Time taken: " << FpMilliseconds(t1-t0).count() << " ms\n";
        std::cout<< "Avg time for "<< thmTime <<" thm sweeps (over T) for L="<< L <<" is " << meanOf(times) <<"ms" <<std::endl;
    }
}