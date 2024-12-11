#include"ising.h"
#include<iostream>
#include<fstream>
#include<string>
#include<chrono>

int main(int argc, char *argv[]){
    std::ofstream data;
    std::vector<double> energy(1000), mag(1000), times(31);
    double specHeat, sus; int count;
    typedef std::chrono::high_resolution_clock Clock;
    using FpMilliseconds = std::chrono::duration<float, std::chrono::milliseconds::period>;
    for(int L = 4; L<=128; L*=2){
        isingLattice lattice(L);    
        data.open("data"+std::to_string(L));
        data << "SUSCEPTIBILITY\tSPECIFIC HEAT\n";
        count = 0;
        for (double beta = 0.3; beta <= 0.6; beta+=0.01)
        {
            lattice.initialize(0.5);
            lattice.boltzmann2D(beta);
            // thermalization
            auto t0=Clock::now();
            for (size_t k = 0; k < 1000; k++)
            {            
                lattice.metropolis2DimSweep(beta);
            }
            auto t1=Clock::now();
            //data collection
            for (int k = 0; k < 1000; k++)
            {
                lattice.metropolis2DimSweep(beta);
                energy[k] = lattice.energy2D();
                mag[k] = lattice.magnetisation();
            }
            sus = beta*L*L*variance(mag);
            specHeat = beta*beta*L*L*variance(energy);
            data << sus << "\t" << specHeat << "\n";
            times[count++] = FpMilliseconds(t1-t0).count(); // time in ms
        }
        data.close();
        std::cout << count << std::endl;
        std::cout << "Finished with " << L<< std::endl;
        std::cout<< "Average time for 1000 sweeps over all temps for L="<< L <<" is " << meanOf(times) <<std::endl;
    }
    data.close();
}