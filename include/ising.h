#ifndef ISING_H
#define ISING_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include"processing.h"
#include"randutils.hpp"

class isingLattice{
    private:
        const int N;
        int *lattice; // lattice
        int *nn; // neighbours
        double probArray[3]={0};
        std::uniform_int_distribution<> lDist;
        std::uniform_real_distribution<> rDist;
        std::mt19937 mt19937Engine;
    public:
        isingLattice(int n, std::mt19937 rng): N(n), lattice(new int[n]), mt19937Engine(rng){
            rDist = std::uniform_real_distribution<>(0.0, 1.0);
            lDist = std::uniform_int_distribution<>(0, N-1);
        };
        void initialise(double m0){
            for(int i=0; i<N; ++i)
                lattice[i] = (rDist(mt19937Engine) < (1+m0)/2.0)? 1 : -1;
        }
        void boltz1dGlauber(double b, double h=0.0){
            for (int i = 0; i < 3; ++i) probArray[i] = 1.0/(1.0+exp(2.0*b*((i*2.0)-2.0+h)));
        };
        void metropolis1DimSweep(double);
        void metropolis1DimSweep(double, double);
        void printLattice(std::ostream &data);
        void glauber1DimSweep(double, double);
        void glauber1dInterval(double,int);
        void kawasaki1DimSweep(double, double);
        double magnetisation(void) const;
        double energy1D(double) const;
        void writeConfig(std::ostream &data) const;
        ~isingLattice() {
            // lattice.clear();
            delete[] lattice;
        };
};

void isingLattice::writeConfig(std::ostream &data) const{
    for (int i = 0; i < N; ++i)
        data << lattice[i];
    data <<"\n";
    // data << std::endl;
}

double isingLattice::energy1D(double h=0.0) const{
    double e = 0;
    for (int i = 0; i < N; ++i){
        e += lattice[i]*(lattice[mod(i-1,N)] + lattice[mod(i+1,N)]);
    }
    return -e/((double) 2.0*N);
}

double isingLattice::magnetisation(void) const{
    double m=0.0;
    for (int i = 0; i < N; ++i) m += lattice[i];
    return m/N;
}

void isingLattice::printLattice(std::ostream &data){
    for(int i=0; i<N; ++i){
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

void isingLattice::metropolis1DimSweep(double beta){
    int ide;
    double boltz = exp(-4*beta);
    for(int i=0; i<N; ++i){
        int x = lDist(mt19937Engine);
        ide = lattice[x]*(lattice[mod(x-1,N)]+lattice[mod(x-1,N)]);
        if (ide <=0 || rDist(mt19937Engine) < boltz){
            lattice[x] = -lattice[x];
        }
    }
}

void isingLattice::metropolis1DimSweep(double beta, double h){
    double ide;
    for(int i=0; i<N; ++i){
        int x = lDist(mt19937Engine);
        ide = lattice[x]*(lattice[mod(x-1,N)]+lattice[mod(x+1,N)] + h);
        if (ide <=0 || rDist(mt19937Engine) < exp(-2*beta*ide)){
            lattice[x] = -lattice[x];
        }
    }
}

void isingLattice::glauber1DimSweep(double beta, double h=0.0){
    int ide;
    boltz1dGlauber(beta, h);
    for(int i=0; i<N; ++i){
        int x = lDist(mt19937Engine);
        ide = static_cast<int>(lattice[x]*(lattice[mod(x-1,N)] + lattice[mod(x+1,N)])/2.0+1.0);
        if (rDist(mt19937Engine) <= probArray[ide]){
            lattice[x] = -lattice[x];
        }
    }
}

// glauber random update over a given number of sites
void isingLattice::glauber1dInterval(double beta, int nSites){
    int ide;
    boltz1dGlauber(beta, 0.0);
    for(int i=0; i<nSites; ++i){
        int x = lDist(mt19937Engine);
        ide = static_cast<int>(lattice[x]*(lattice[mod(x-1,N)] + lattice[mod(x+1,N)])/2.0+1.0);
        // if (rDist(mt19937Engine) <= 1.0/(1.0+exp(2.0*ide*beta))){
        if (rDist(mt19937Engine) <= probArray[ide]){
            lattice[x] = -lattice[x];
        }
    }
}

// 1d kawasaki random update
void isingLattice::kawasaki1DimSweep(double beta, double h=0.0){
    double ide, boltz = exp(-beta*4); int tmp, x, y;
    for(int i=0; i<N; ++i){
        x = lDist(mt19937Engine); y = mod(x+1,N);
        if (lattice[x]==lattice[y]) continue;
        ide = lattice[x]*lattice[mod(x-1,N)] + lattice[y]*lattice[mod(y+1,N)];
        if (rDist(mt19937Engine) <= 1.0/(1.0+exp(2.0*ide*beta))){
        // if (ide<=0 || rDist(mt19937Engine) < boltz){
            tmp = lattice[x];
            lattice[x] = lattice[y];
            lattice[y] = tmp;
        }
    }
}

#endif