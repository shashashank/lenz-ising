// 1D Ising Resetting for magnetisation

#include"ising.h"
#include<chrono>
#include<iomanip>
#include<omp.h>
// #include<gsl_randist.h>

int main(int argc, char **argv){

    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h")){
        // To do: write help
    }

    // const std::string &filename = input.getCmdOption("-f");
    // if (!filename.empty()){
    //     // Do interesting things ...
    // }

    double T;
    const std::string &Tstring = input.getCmdOption("-T");
    if (!Tstring.empty()){
        T = std::stod(Tstring);
        std::cout << "T = " << T << std::endl;
    }else{
        std::cerr << "T not provided" << std::endl;
        exit(0);
    }

    int t=10000;
    const std::string &tstring = input.getCmdOption("-t");
    if (!tstring.empty()){
        t = std::stoi(tstring);
        std::cout << "t = " << t << std::endl;
    }else{
        std::cerr << "t not provided" << std::endl;
        exit(0);
    }

    double m0;
    const std::string &m0string = input.getCmdOption("-m0");
    if (!m0string.empty()){
        m0 = std::stod(m0string);
        std::cout << "m0 = " << m0 << std::endl;
    }else{
        std::cerr << "m0 not provided" << std::endl;
        exit(0);
    }

    double r;
    const std::string &rstring = input.getCmdOption("-r");
    if (!rstring.empty()){
        r = std::stod(rstring);
        std::cout << "r = " << r << std::endl;
    }else{
        std::cerr << "r not provided" << std::endl;
        exit(0);
    }

    typedef std::chrono::high_resolution_clock Clock;
    std::ofstream data;
    data.open("magnetisation"+std::to_string(r)+".txt");
    const int N = 10000;
    // double T = 3.5;
    const double beta = 1.0/T;
    // double m0 = 0.992;
    // double r = 3.0;
    std::exponential_distribution<> expDis(r);
    // int totalSteps = 0;
    // while (totalSteps<10000){

    omp_set_num_threads(4);
    omp_set_dynamic(0);

    auto t0=Clock::now();
    #pragma omp parallel for schedule(dynamic) shared(data)
    for(int totalSteps = 0; totalSteps < t; totalSteps++){
        
        // iters from exp dist scaled to lattice size and spin flips attempted
        int stepsItr =  static_cast<int>(expDis(mt19937Engine)*N);

        isingLattice lattice(N);
        lattice.initialise(m0); // initialising the lattice           
        lattice.glauber1dInterval(beta, stepsItr);

        double mag = lattice.magnetisation();

        #pragma omp critical
        {
            data << std::fixed<<std::setprecision(10) << mag <<"\n";
        }
    }
    auto t1=Clock::now();
    std::chrono::duration<double> dts = (t1-t0); // time in s
    std::cout << "Time in s: "<< dts.count() << std::endl;
    data.close();
    
}