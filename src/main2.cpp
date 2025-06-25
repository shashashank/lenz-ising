// 2D Ising Resetting for magnetisation and energy

#include"ising2d.h"
#include<chrono>
#include<iomanip>

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

    int totalSteps=10000;
    const std::string &tstring = input.getCmdOption("-t");
    if (!tstring.empty()){
        totalSteps = std::stoi(tstring);
    }
    // std::cout << "no. of trials t = " << t << "\n";

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

    int L;
    const std::string &Nstring = input.getCmdOption("-L");
    if (!Nstring.empty()){
        L = std::stoi(Nstring);
    }else{
        L = 100;
    }
    // std::cout << "no. of spins N = " << N << "\n";

    static uint64_t seed;
    const std::string &seedstring = input.getCmdOption("-s");
    if (!seedstring.empty()){
        seed = std::stoull(seedstring);
    }else{
        seed = 328575958951598690;
    }
    // std::cout << "seed = " << seed << "\n";

    std::ofstream reset;
    reset.open("resetData"+std::to_string(r)+".txt");
    reset.setf(std::ios::fixed); reset.precision(10);
    const double beta = 1.0/T;
    
    std::exponential_distribution<> expDis(r);
    randutils::seed_seq_fe128 seeder{uint32_t(seed),uint32_t(seed >> 32)};
    std::mt19937 mt19937Engine(seeder);
    isingLattice2D lattice(L, &mt19937Engine);
    
    double energyTmp, magTmp;
    for(int t = 0; t < totalSteps; t++){
        // iters from exp dist scaled to lattice size and spin flips attempted
        int stepsItr =  static_cast<int>(expDis(mt19937Engine)*L*L);
        lattice.initialise(m0); // initialising the lattice
        lattice.glauberInterval(beta, stepsItr);
        energyTmp = lattice.energy();
        magTmp = lattice.magnetisation();
        reset << energyTmp << "\t" << magTmp << "\n";
    }
    reset.close();
}