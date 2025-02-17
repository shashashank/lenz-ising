// This file creates an executable for simulating the 1d ising model for statistics without resetting

#include"ising.h"
#include<chrono>
#include <iomanip>

int main(int argc, char *argv[]){

    const int thmTime=50000, simTime = 100, steps=1;
    const int N = 1000;
    // const int nTrials = 100;
    // double betaInitial=1.0, betaFinal = 1.0/1.1;

    std::ofstream data;    
    std::vector<double> energy(steps*simTime), mag(steps*simTime);
    // std::vector<double> times(301); int count;
    double specHeat, sus;
    typedef std::chrono::high_resolution_clock Clock;
    using FpMilliseconds = std::chrono::duration<float, std::chrono::milliseconds::period>;

    isingLattice lattice(N);

    data.open("data");
    data<< "#BETA\tMAGNETISATION\tENERGY\n";

    // for (int j =0; j<nTrials; j++){
    //     isingLattice lattice(N);
    for (double j =0; j<10; j+=0.1){
        double beta = 1.0/j;

        // initialization
        lattice.initialise(0.5);

        //thermalization
        for (int k = 0; k < thmTime; k++){
            lattice.glauber1DimSweep(1.0/j);
        }

        // dynamics and data collection
        for (size_t k = 0; k < simTime; k++){
            energy[k] = lattice.energy1D();
            mag[k] = lattice.magnetisation();
            lattice.glauber1dInterval(beta, N);    
        }
        
        // for (int k = 0; k < simTime; k++){
        //     for (size_t i = 0; i < steps; i++){
        //         energy[k*steps+i] = lattice.energy1D();
        //         mag[k*steps+i] = lattice.magnetisation();
        //         lattice.glauber1dInterval(betaFinal, N/steps);
        //     }
        // }

        // writing to disk
        // data.open("data"+std::to_string(j));
        // // data<< "#BETA\tMAGNETISATION\tENERGY\n";
        for (size_t i = 0; i < steps*simTime; i++)
        {
            data << std::fixed<<std::setprecision(10)<<1.0/j<<"\t";
            data << std::fixed<<std::setprecision(10)<<mag[i]<<"\t";
            data << std::fixed<<std::setprecision(10)<<energy[i]<<"\n";
        }
        // data.close();
        std::cout << "Fininshed with " << j*10 << " trials" <<std::endl;
    }
    data.close();
}