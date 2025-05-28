// 1D Ising Resetting for magnetisation

#include"ising.h"
#include<chrono>
#include<iomanip>
// #include<omp.h>
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


    int N;
    const std::string &Nstring = input.getCmdOption("-N");
    if (!Nstring.empty()){
        N = std::stoi(Nstring);
    }else{
        N = 10000;
    }
    std::cout << "N = " << N << std::endl;

    // double r;
    // const std::string &rstring = input.getCmdOption("-r");
    // if (!rstring.empty()){
    //     r = std::stod(rstring);
    //     std::cout << "r = " << r << std::endl;
    // }else{
    //     std::cerr << "r not provided" << std::endl;
    //     exit(0);
    // }

    typedef std::chrono::high_resolution_clock Clock;
    // double T = 3.5;
    const double beta = 1.0/T;
    // double m0 = 0.992;
    // double r = 3.0;
    // int totalSteps = 0;
    // while (totalSteps<10000){

    // omp_set_num_threads(4);
    // omp_set_dynamic(0);

    auto t0=Clock::now();
    isingLattice lattice(N);
    double maxE = 0, minE=0, eTmp;
    double maxM = 0, minM=0, mTmp;
    std::vector<double> energies(t), mags(t);
    // #pragma omp parallel for schedule(dynamic) shared(data)
    for(int totalSteps = 0; totalSteps < t; totalSteps++){
        
        // iters from exp dist scaled to lattice size and spin flips attempted
        // int stepsItr =  static_cast<int>(expDis(mt19937Engine)*N);

        // isingLattice lattice(N);
        lattice.initialise(m0); // initialising the lattice           
        // lattice.glauber1dInterval(beta, stepsItr);

        // double mag = lattice.magnetisation();

        // #pragma omp critical
        // {
        //     data << std::fixed<<std::setprecision(10) << mag <<"\n";
        // }
        eTmp = lattice.energy1D();
        mTmp = lattice.magnetisation(); 
        energies[totalSteps] = eTmp;
        mags[totalSteps] = mTmp;
        if (totalSteps==0){
            maxE = eTmp;
            minE = eTmp;
            minM = mTmp;
            maxM = mTmp;
        }
        if (eTmp>maxE){
            maxE = eTmp;
        }
        if (eTmp<minE){
            minE = eTmp;
        }

        if (mTmp>maxM){
            maxM = mTmp;
        }

        if (mTmp<minM){
            minM = mTmp;
        }

    }
    std::cout<<"Energy"<<std::endl;
    std::cout << "max = " << maxE << ";\t";
    std::cout << "min = " << minE << ";\t";
    std::cout << "avg = " << std::accumulate(energies.begin(), energies.end(), 0.0)/energies.size() << ";\t";
    std::cout << "var = " << std::accumulate(energies.begin(), energies.end(), 0.0, 
        [&](const double previous, const double current){
            return previous + (current - std::accumulate(energies.begin(), energies.end(), 0.0)/energies.size())*(current - std::accumulate(energies.begin(), energies.end(), 0.0)/energies.size());
        })/energies.size() << std::endl;
    std::cout<<"Magnetisation"<<std::endl;
    std::cout << "max = " << maxM << ";\t";
    std::cout << "min = " << minM << ";\t";
    std::cout << "avg = " << std::accumulate(mags.begin(), mags.end(), 0.0)/mags.size() << ";\t";
    std::cout << "var = " << std::accumulate(mags.begin(), mags.end(), 0.0, 
        [&](const double previous, const double current){
            return previous + (current - std::accumulate(mags.begin(), mags.end(), 0.0)/mags.size())*(current - std::accumulate(mags.begin(), mags.end(), 0.0)/mags.size());
        })/mags.size() << std::endl;
    // auto t1=Clock::now();
    // std::chrono::duration<double> dts = (t1-t0); // time in s
    // std::cout << "Time in s: "<< dts.count() << std::endl;
    // data.close();
    
}