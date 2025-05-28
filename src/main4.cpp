// 1D Ising Resetting for energy

#include"ising.h"
#include<chrono>
#include<iomanip>
#include<omp.h>
#include<format>
// #include<gsl_randist.h>

int main(int argc, char **argv){
    int N = 1000;
    isingLattice lattice(N);
    double energyTmp, magTmp;
    std::ofstream data;
    data.open("data.txt");
    for(int t = 0; t < 10000; t++){

        lattice.initialise(0.992);
        // lattice.glauber1dInterval(0, 10000);

        for (int i = 0; i < 1000; i++){
            lattice.glauber1dInterval(1, N/10);
            // energyTmp = lattice.energy1D();
            magTmp = lattice.magnetisation();
            // data << std::fixed<<std::setprecision(10) << energyTmp << "\t";
            data << std::fixed<<std::setprecision(10) << magTmp    << "\n";
        }
    }
}





    // InputParser input(argc, argv);
    // if(input.cmdOptionExists("-h")){
    //     // To do: write help
    // }

    // // const std::string &filename = input.getCmdOption("-f");
    // // if (!filename.empty()){
    // //     // Do interesting things ...
    // // }

    // double T;
    // const std::string &Tstring = input.getCmdOption("-T");
    // if (!Tstring.empty()){
    //     T = std::stod(Tstring);
    //     std::cout << "T = " << T << std::endl;
    // }else{
    //     std::cerr << "T not provided" << std::endl;
    //     exit(0);
    // }

    // int t=10000;
    // const std::string &tstring = input.getCmdOption("-t");
    // if (!tstring.empty()){
    //     t = std::stoi(tstring);
    //     std::cout << "t = " << t << std::endl;
    // }else{
    //     std::cerr << "t not provided" << std::endl;
    //     exit(0);
    // }

    // double r;
    // const std::string &rstring = input.getCmdOption("-r");
    // if (!rstring.empty()){
    //     r = std::stod(rstring);
    //     std::cout << "r = " << r << std::endl;
    // }else{
    //     std::cerr << "r not provided" << std::endl;
    //     exit(0);
    // }

    // int N;
    // const std::string &Nstring = input.getCmdOption("-N");
    // if (!Nstring.empty()){
    //     N = std::stoi(Nstring);
    //     std::cout << "N = " << N << std::endl;
    // }else{
    //     N = 10000;
    // }

    // typedef std::chrono::high_resolution_clock Clock;
    // std::ofstream data;
    // data.open(std::format("energy{0}_{1}.txt", T, r));
    // // const int N = 10000;
    // const double beta = 1.0/T;
    // std::exponential_distribution<> expDis(r);
    // auto t0=Clock::now();
    // isingLattice lattice(N);
    // for(int totalSteps = 0; totalSteps < t; totalSteps++){
    
    //     // iters from exp dist scaled to lattice size and spin flips attempted
    //     int stepsItr =  static_cast<int>(expDis(mt19937Engine)*N);

    //     lattice.initialise(0.992); // initialising the lattice
    //     for (size_t i = 0; i < 1000; i++){
    //         lattice.glauber1DimSweep(100);
    //     }

    //     lattice.glauber1dInterval(beta, stepsItr);

    //     double energy = lattice.energy1D();
    //     data << std::fixed<<std::setprecision(10)<<energy<<"\n";
        
    // }
    // auto t1=Clock::now();
    // std::chrono::duration<double> dts = (t1-t0); // time in s
    // std::cout << "Time in s: "<< dts.count() << std::endl;
    // data.close();
