#ifndef ISING_H
#define ISING_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include"processing.h"

class isingLattice{
    private:
        const int N, nSqrd;
        // std::vector<int> lattice;
        int* lattice;
        double probArray[5]={0};
        // std::uniform_real_distribution<> rDist(0.0, 1.0), latticeCoordDist(0, nSqrd);
    public:
        isingLattice(int N): N(N), nSqrd(N*N), lattice(new int[N*N]){
            // lattice.resize(nSqrd);
        };
        void initialise(double m0){
            for(size_t i=0; i<nSqrd; ++i)
                lattice[i] = (rDist(mt19937Engine) < (1+m0)/2.0)? 1 : -1;
        }
        void boltz2dMetro(double b){
            for (int i = 2; i < 5; i+=2) probArray[i] = exp(-2*b*i);
        };
        void boltz2dGlauber(double b, double h=0){
            for (int i = 0; i < 5; ++i) probArray[i] = 1/(1+exp(-2.0*b*((i*2.0)-4.0+h)));
        };
        void metropolis2DimSweep(double);
        void metropolis1DimSweep(double);
        void metropolis1DimSweep(double, double);
        void printLattice(std::ostream &data);
        void glauber1DimSweep(double, double);
        void glauber1dInterval(double,int);
        void glauber1dIntervalTyp(double,int, int);
        void glauber2DimSweep(double, double);
        void kawasaki1DimSweep(double, double);
        double magnetisation(void) const;
        double energy1D(double) const;
        double energy2D(double) const;
        double siteEnergy2D(int, int, double) const;
        void writeConfig(std::ostream &data) const;
        ~isingLattice() {
            // lattice.clear();
            delete[] lattice;
        };
};

void isingLattice::writeConfig(std::ostream &data) const{
    for (int i = 0; i < nSqrd; ++i)
        data << lattice[i];
    data <<"\n";
    // data << std::endl;
}

double isingLattice::energy1D (double h=0.0) const{
    double e = 0;
    for (int i = 0; i < nSqrd; ++i)
        // e += lattice[i]*(lattice[((i==0)?nSqrd-1:i-1)]+lattice[((i==nSqrd-1)?0:i+1)]);
        e += lattice[i]*(lattice[mod(i-1,nSqrd)] + lattice[mod(i+1,nSqrd)]);
    return -e/((double) 2.0*nSqrd);
}

double isingLattice::energy2D(double h=0.0) const{
    double e = 0.0;
    for(int x=0; x<N; ++x){
        for (int y = 0; y < N; y++)
            e = siteEnergy2D(x,y,h);
    }
    return e/((double) 2.0*nSqrd);
}

double isingLattice::magnetisation(void) const{
    double m=0.0;
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

// multiplies with a -1
double isingLattice::siteEnergy2D(int x, int y, double h=0.0) const{
    double energy = 0.0;
    // a = lattice[N*y+((x==0)?N-1:x-1)];
    // b = lattice[N*y+((x==N-1)?0:x+1)];
    // c = lattice[N*((y==0)?N-1:y-1)+x];
    // d = lattice[N*((y==N-1)?0:y+1)+x];
    // ide = lattice[N*y+x]*(a+b+c+d);
    energy += lattice[mod(y,N)*N +   mod(x+1,N)];
    energy += lattice[mod(y+1,N)*N + mod(x,N)];
    energy += lattice[mod(y,N)*N +   mod(x-1,N)];
    energy += lattice[mod(y-1,N)*N + mod(x,N)];
    energy = (energy + h)*(-lattice[N*y+x]);
    return energy;
}

void isingLattice::metropolis2DimSweep(double beta){
    int ide;
    std::uniform_int_distribution<> lDist(0, nSqrd);
    // to do use l dist
    boltz2dMetro(beta);
    for(int x=0; x<N; ++x){
        for (int y = 0; y < N; y++){
            ide = (int) -siteEnergy2D(x, y);
            if (ide <=0 || rDist(mt19937Engine) < probArray[ide]){
            // if(ide<=0 || rDist(mt19937Engine) <  exp(-2*beta*ide)){
                lattice[N*y+x] = -lattice[N*y+x];
            }
        }
    }
}

void isingLattice::metropolis1DimSweep(double beta){
    int ide;

    std::uniform_int_distribution<> lDist(0, nSqrd);
    // to do use ldist
    double boltz = exp(-4*beta);
    for(int x=0; x<nSqrd; ++x){
        ide = lattice[x]*(lattice[((x==0)?nSqrd-1:x-1)]+lattice[((x==nSqrd-1)?0:x+1)]);
        if (ide <=0 || rDist(mt19937Engine) < boltz){
            lattice[x] = -lattice[x];
        }
    }
}

void isingLattice::metropolis1DimSweep(double beta, double h){
    double ide;

    std::uniform_int_distribution<> lDist(0, nSqrd);
    for(int x=0; x<nSqrd; ++x){
        ide = lattice[x]*(lattice[((x==0)?nSqrd-1:x-1)]+lattice[((x==nSqrd-1)?0:x+1)] + h);
        if (ide <=0 || rDist(mt19937Engine) < exp(-2*beta*ide)){
            lattice[x] = -lattice[x];
        }
    }
}

void isingLattice::glauber1DimSweep(double beta, double h=0.0){
    double ide;
    std::uniform_int_distribution<> lDist(0, nSqrd);
    for(int i=0; i<nSqrd; ++i){
        int x = lDist(mt19937Engine);
        // ide = lattice[x]*(lattice[((x==0)?nSqrd-1:x-1)]+lattice[((x==nSqrd-1)?0:x+1)] + h);
        ide = lattice[x]*(lattice[mod(x-1,nSqrd)] + lattice[mod(x+1,nSqrd)]);
        if (rDist(mt19937Engine) <= 1.0/(1.0+exp(2.0*ide*beta))){
            lattice[x] = -lattice[x];
        }
    }
}

void isingLattice::glauber1dInterval(double beta, int end){
    double ide; int count =0;
    std::uniform_int_distribution<> lDist(0, nSqrd);
    for(int i=0; i<end; ++i){
        int x = lDist(mt19937Engine);
        ide = lattice[x]*(lattice[mod(x-1,nSqrd)] + lattice[mod(x+1,nSqrd)]);
        if (rDist(mt19937Engine) <= 1.0/(1.0+exp(2.0*ide*beta))){
            lattice[x] = -lattice[x];
        }
    }
}

void isingLattice::glauber1dIntervalTyp(double beta, int start, int finish){
    double ide; int count =0;
    for(int x=start; x<finish; ++x){
        ide = lattice[x]*(lattice[mod(x-1,nSqrd)] + lattice[mod(x+1,nSqrd)]);
        if (rDist(mt19937Engine) <= 1.0/(1.0+exp(2.0*ide*beta))){
            lattice[x] = -lattice[x];
        }
    }
}


void isingLattice::kawasaki1DimSweep(double beta, double h=0.0){
    double ide, boltz = exp(-beta*4); int tmp, x, y;
    std::uniform_int_distribution<> discr(0, nSqrd-1);
    for(int i=0; i<nSqrd; ++i){
        x = discr(mt19937Engine); y = mod(x+1,nSqrd);
        if (lattice[x]==lattice[y]){continue;}
        ide = lattice[x]*lattice[mod(x-1,nSqrd)] + lattice[y]*lattice[mod(y+1,nSqrd)];
        if (rDist(mt19937Engine) <= 1.0/(1.0+exp(2.0*ide*beta))){
        // if (ide<=0 || rDist(mt19937Engine) < boltz){
            tmp = lattice[x];
            lattice[x] = lattice[y];
            lattice[y] = tmp;
        }
    }
}

void isingLattice::glauber2DimSweep(double beta, double h=0.0){
    std::uniform_int_distribution<> lDist(0, nSqrd);
    // to do use lDist and not typewriter
    for(int x=0; x<N; ++x){
        for(int y=0; y<N; ++y){
            if (rDist(mt19937Engine) <= 1.0/(1.0+exp(-2.0*siteEnergy2D(x, y, h)*beta))){
                lattice[N*y+x] = -lattice[N*y+x];
            }
        }
    }
}

#endif