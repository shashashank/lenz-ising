#ifndef ISING_H
#define ISING_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <random>
#include"randutils.hpp"

static uint64_t seed = 328575958941598690;
randutils::seed_seq_fe128 seeder{uint32_t(seed),uint32_t(seed >> 32)};
std::mt19937 mt19937Engine(seeder);
std::uniform_real_distribution<> rDist(0.0, 1.0);

template<typename T> class spin{
    private:
        // bool state:1;
        T state;
    public:
        spin(bool state){
            this->state = state;
        }
        spin(){
            this->state = 1;
        }
        bool getState(void){
            return this->state;
        }
        void flip(void){
            this->state = !this->state;
        }
        spin upSpin(void) {return spin(true);}
        spin downSpin(void) {return spin(false);}
        spin& operator=(const spin &rhs){
            this->state = rhs.state;
            return *this;
        }
        int operator+(const spin &rhs){
            return this->state + rhs.state;
        }
        int operator*(const spin &rhs){
            return this->state * rhs.state;
        }
        int operator-(const spin &rhs){
            return this->state - rhs.state;
        }
        // void operator-(const spin &rhs){
        //     this->state.flip();
        // }
};

class isingLattice{
    private:
        const long unsigned int N, nSqrd;
        std::vector<int> lattice;
        double probArray[5]={0};
        // std::uniform_real_distribution<> rDist(0.0, 1.0), latticeCoordDist(0, nSqrd);
    public:
        isingLattice(long unsigned int N): N(N), nSqrd(N*N){
            lattice.resize(nSqrd);
            for (size_t i = 0; i < nSqrd; i++) lattice[i] = 1;
        };
        void initialize(double p){
            for(size_t i=0; i<nSqrd; ++i)
                lattice[i] = (rDist(mt19937Engine) < p)? -1 : 1;
        }
        void boltzmann2D(double b){
            for (int i = 2; i < 5; ++i) probArray[i] = exp(-2*b*i);
        };
        void metropolis2DimSweep(double);
        void metropolis1DimSweep(double);
        void metropolis1DimSweep(double, double);
        void printLattice(std::ostream &data);
        void glauber1DimSweep(double, double);
        double magnetisation(void) const;
        double energy1D(void) const;
        double energy2D(void) const;
};

double isingLattice::energy1D (void) const{
    double e = 0;
    for (int i = 0; i < nSqrd; ++i)
        e += lattice[i]*(lattice[((i==0)?nSqrd-1:i-1)]+lattice[((i==nSqrd-1)?0:i+1)]);
    return -e/nSqrd;
}

double isingLattice::energy2D(void) const{
    double e = 0.0, a, b, c, d;
    for(int x=0; x<N; ++x){
        for (int y = 0; y < N; y++){
            a = lattice[N*y+((x==0)?N-1:x-1)];
            b = lattice[N*y+((x==N-1)?0:x+1)];
            c = lattice[N*((y==0)?N-1:y-1)+x];
            d = lattice[N*((y==N-1)?0:y+1)+x];
            e += lattice[N*y+x]*(a+b+c+d);
        }
    }
    return -e/((double) nSqrd);
}

double isingLattice::magnetisation(void) const{
    double m=0;
    for (int i = 0; i < nSqrd; ++i) m += lattice[i];
    return m/nSqrd;
}

void isingLattice::printLattice(std::ostream &data){
    for(size_t i=0; i<nSqrd; ++i){
        switch (lattice[i]){
            case 1:
                data.write("1", 1);
                break;
            case -1:
                data.write("0", 1);
                break;
            default:
                std::cerr << "Error: lattice value not 1 or -1" << std::endl;
                exit(0);
        }
    }
    data.write("\n", 1);
}

void isingLattice::metropolis2DimSweep(double b){
    int ide, a, c, d;
    for(int x=0; x<N; ++x){
        for (int y = 0; y < N; y++){
            a = lattice[N*y+((x==0)?N-1:x-1)];
            b = lattice[N*y+((x==N-1)?0:x+1)];
            c = lattice[N*((y==0)?N-1:y-1)+x];
            d = lattice[N*((y==N-1)?0:y+1)+x];
            ide = lattice[N*y+x]*(a+b+c+d);
            if (ide <=0 || rDist(mt19937Engine) < probArray[ide]){
                lattice[N*y+x] = -lattice[N*y+x];
            }
        }
    }
}

void isingLattice::metropolis1DimSweep(double b){
    int ide;
    double boltz = exp(-4*b);
    for(int x=0; x<nSqrd; ++x){
        ide = lattice[x]*(lattice[((x==0)?nSqrd-1:x-1)]+lattice[((x==nSqrd-1)?0:x+1)]);
        if (ide <=0 || rDist(mt19937Engine) < boltz){
            lattice[x] = -lattice[x];
        }
    }
}

void isingLattice::metropolis1DimSweep(double b, double h){
    int ide;
    double boltz = exp(-4*b);
    for(int x=0; x<nSqrd; ++x){
        ide = lattice[x]*(lattice[((x==0)?nSqrd-1:x-1)]+lattice[((x==nSqrd-1)?0:x+1)] + h);
        if (ide <=0 || rDist(mt19937Engine) < exp(-2*b*ide)){
            lattice[x] = -lattice[x];
        }
    }
}

void isingLattice::glauber1DimSweep(double b, double h){
    int ide;
    for(int x=0; x<nSqrd; ++x){
        ide = lattice[x]*(lattice[((x==0)?nSqrd-1:x-1)]+lattice[((x==nSqrd-1)?0:x+1)] + h);
        if (rDist(mt19937Engine) <= 1/(1+exp(-2*ide))){
            lattice[x] = -lattice[x];
        }
    }
}

template<typename T> double meanOf(const T& array){
    double mean = 0;
    for (auto& val: array){
        mean+=val;
    }
    return mean/((double)array.size());
}

template<typename T> double variance(const T& array){
    double mean = meanOf(array);
    double variance = 0.0;
    for(auto& val : array){
        auto diff = val - mean;
        variance += diff * diff;
    }
    return variance/((double)array.size());
}

// void updateGlauber(const size_t n, const double b){
//     double r, W, k;
//     int deltaE, S;

//     S = lattice[mod(n+1,nSqrd)] + lattice[mod(n+N,nSqrd)] + lattice[mod(n-1,nSqrd)] + lattice[mod(n-N,nSqrd)] + h;
//     deltaE = 2*(lattice[n])*S;   //-ve taken care of by S

//     r = rDist(gen);
//     k = exp(-b*deltaE);
//     W = k/(1.0+k);

//     if (r <= W){
//         lattice[n] = -lattice[n];
//     }
// }

#endif