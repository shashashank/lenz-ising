#ifndef ISING_H
#define ISING_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include"processing.h"
#include"randutils.hpp"

static uint64_t seed = 328575958951598690;
randutils::seed_seq_fe128 seeder{uint32_t(seed),uint32_t(seed >> 32)};
std::mt19937 mt19937Engine(seeder);
std::uniform_real_distribution<> rDist(0.0, 1.0);

class isingLattice{
    private:
        const int N, dim;
        int *lattice; // lattice
        int *nn, *nnn, *plaquette; // neighbours
        double probArray[5]={0};
        // std::uniform_real_distribution<> rDist(0.0, 1.0), latticeCoordDist(0, nSqrd);
        void neighbourGen(int dim){
            if (dim==2){
                int L = static_cast<int>(sqrt(N));
                for (int i = 0; i < N; ++i)
                {
                    nn[4*i+0] = (mod(i+1,L)==0) ? i-L+1 : i+1;
                    nn[4*i+1] = (mod(i-1,L)==0) ? i+L-1 : i-1;
                    nn[4*i+2] = mod(i-L,N);
                    nn[4*i+3] = mod(i+L,N);
                }
            }
            if (dim==3){
                int L = static_cast<int>(pow(N, 1.0/3.0));
                int L2 = L*L;
                for (int i = 0; i < N; ++i)
                {
                    nn[6*i+0] = (mod(i+1,L)==0) ? i-L+1 : i+1;
                    nn[6*i+1] = (mod(i-1,L)==0) ? i+L-1 : i-1;
                    nn[6*i+2] = mod(i-L,N);
                    nn[6*i+3] = mod(i+L,N);
                    nn[6*i+4] = mod(i-L2,N);
                    nn[6*i+5] = mod(i+L2,N);
                }
                for (int i = 0; i < N; i++)
                {
                    nnn[8*i+0] = nn[6*nn[6*i+0]+2];
                    nnn[8*i+1] = nn[6*nn[6*i+0]+2];
                    nnn[8*i+2] = nn[6*nn[6*i+1]+3];
                    nnn[8*i+3] = nn[6*nn[6*i+1]+3];
                    nnn[8*i+4] = nn[6*nn[6*i+0]+4];
                    nnn[8*i+5] = nn[6*nn[6*i+1]+4];
                    nnn[8*i+6] = nn[6*nn[6*i+0]+5];
                    nnn[8*i+7] = nn[6*nn[6*i+1]+5];
                }

                for (int i = 0; i < N; i++)
                {
                    plaquette[14*i+0] = nn[6*i+0];
                    plaquette[14*i+1] = nn[6*i+1];
                    plaquette[14*i+2] = nn[6*i+2];
                    plaquette[14*i+3] = nn[6*i+3];
                    plaquette[14*i+4] = nn[6*i+4];
                    plaquette[14*i+5] = nn[6*i+5];
                    plaquette[14*i+6] = nnn[8*i+0];
                    plaquette[14*i+7] = nnn[8*i+1];
                    plaquette[14*i+8] = nnn[8*i+2];
                    plaquette[14*i+9] = nnn[8*i+3];
                    plaquette[14*i+10] = nnn[8*i+4];
                    plaquette[14*i+11] = nnn[8*i+5];
                    plaquette[14*i+12] = nnn[8*i+6];
                    plaquette[14*i+13] = nnn[8*i+7];
                }
                
                

            }
            
        }
    public:
        isingLattice(int n): N(n), dim(1),lattice(new int[n]){
            nn = nullptr;
            nnn = nullptr;
            plaquette = nullptr;
        };

        isingLattice(int n, int d): N(n), dim(d), lattice(new int[n]){
            nn = new int[2*d*n];
            // nnn = new int[4*(d-1)*n];
            if (d==3){
                nnn = new int[8*n];
                plaquette = new int[14*n];
            }else{
                nnn = nullptr;
                plaquette = nullptr;
            }
            neighbourGen(d);
        };

        void initialise(double m0){
            for(size_t i=0; i<N; ++i)
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
        double energy3D(double) const;
        double siteEnergy2D(int, double) const;
        double siteEnergy3D(int, double) const;
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

double isingLattice::energy2D(double h=0.0) const{
    double e = 0.0;
    for(int x=0; x< N; ++x)
        e = siteEnergy2D(x,h);
    return -e/((double) 2.0*N);
}

double isingLattice::energy3D(double h=0.0) const{
    double e = 0.0;
    for(int x=0; x< N; ++x)
        e = siteEnergy3D(x,h);
    return -e/((double) 2.0*N);
}

double isingLattice::magnetisation(void) const{
    double m=0.0;
    for (int i = 0; i < N; ++i) m += lattice[i];
    return m/N;
}

void isingLattice::printLattice(std::ostream &data){
    for(size_t i=0; i<N; ++i){
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
double isingLattice::siteEnergy2D(int x, double h=0.0) const{
    double energy = 0.0;
    energy += lattice[nn[x+0]];
    energy += lattice[nn[x+1]];
    energy += lattice[nn[x+2]];
    energy += lattice[nn[x+3]];
    energy = (energy + h)*(-lattice[x]);
    return energy;
}

double isingLattice::siteEnergy3D(int x, double h=0.0) const{
    double energy = 0.0;
    energy += lattice[nn[x+0]];
    energy += lattice[nn[x+1]];
    energy += lattice[nn[x+2]];
    energy += lattice[nn[x+3]];
    energy += lattice[nn[x+4]];
    energy += lattice[nn[x+5]];
    energy = (energy + h)*(-lattice[x]);
    return energy;
}

void isingLattice::metropolis2DimSweep(double beta){
    int ide;
    std::uniform_int_distribution<> lDist(0, N);
    // to do use l dist
    boltz2dMetro(beta);
    for(int x=0; x<N; ++x){
        ide = (int) -siteEnergy2D(x);
        if (ide <=0 || rDist(mt19937Engine) < probArray[ide]){
            lattice[x] = -lattice[x];
        }
    }
}

void isingLattice::metropolis1DimSweep(double beta){
    int ide;

    std::uniform_int_distribution<> lDist(0, N);
    // to do use ldist
    double boltz = exp(-4*beta);
    for(int x=0; x<N; ++x){
        ide = lattice[x]*(lattice[((x==0)?N-1:x-1)]+lattice[((x==N-1)?0:x+1)]);
        if (ide <=0 || rDist(mt19937Engine) < boltz){
            lattice[x] = -lattice[x];
        }
    }
}

void isingLattice::metropolis1DimSweep(double beta, double h){
    double ide;

    std::uniform_int_distribution<> lDist(0, N);
    for(int x=0; x<N; ++x){
        ide = lattice[x]*(lattice[((x==0)?N-1:x-1)]+lattice[((x==N-1)?0:x+1)] + h);
        if (ide <=0 || rDist(mt19937Engine) < exp(-2*beta*ide)){
            lattice[x] = -lattice[x];
        }
    }
}

void isingLattice::glauber1DimSweep(double beta, double h=0.0){
    double ide;
    std::uniform_int_distribution<> lDist(0, N);
    for(int i=0; i<N; ++i){
        int x = lDist(mt19937Engine);
        // ide = lattice[x]*(lattice[((x==0)?nSqrd-1:x-1)]+lattice[((x==nSqrd-1)?0:x+1)] + h);
        ide = lattice[x]*(lattice[mod(x-1,N)] + lattice[mod(x+1,N)]);
        if (rDist(mt19937Engine) <= 1.0/(1.0+exp(2.0*ide*beta))){
            lattice[x] = -lattice[x];
        }
    }
}

// glauber random update over a given number of sites
void isingLattice::glauber1dInterval(double beta, int nSites){
    double ide; int count =0;
    std::uniform_int_distribution<> lDist(0, N);
    for(int i=0; i<nSites; ++i){
        int x = lDist(mt19937Engine);
        ide = lattice[x]*(lattice[mod(x-1,N)] + lattice[mod(x+1,N)]);
        if (rDist(mt19937Engine) <= 1.0/(1.0+exp(2.0*ide*beta))){
            lattice[x] = -lattice[x];
        }
    }
}

// 1d glauber typewriter update over an interval from start to finish
void isingLattice::glauber1dIntervalTyp(double beta, int start, int finish){
    double ide; int count =0;
    for(int x=start; x<finish; ++x){
        ide = lattice[x]*(lattice[mod(x-1,N)] + lattice[mod(x+1,N)]);
        if (rDist(mt19937Engine) <= 1.0/(1.0+exp(2.0*ide*beta))){
            lattice[x] = -lattice[x];
        }
    }
}

// 1d kawasaki random update
void isingLattice::kawasaki1DimSweep(double beta, double h=0.0){
    double ide, boltz = exp(-beta*4); int tmp, x, y;
    std::uniform_int_distribution<> discr(0, N-1);
    for(int i=0; i<N; ++i){
        x = discr(mt19937Engine); y = mod(x+1,N);
        if (lattice[x]==lattice[y]){continue;}
        ide = lattice[x]*lattice[mod(x-1,N)] + lattice[y]*lattice[mod(y+1,N)];
        if (rDist(mt19937Engine) <= 1.0/(1.0+exp(2.0*ide*beta))){
        // if (ide<=0 || rDist(mt19937Engine) < boltz){
            tmp = lattice[x];
            lattice[x] = lattice[y];
            lattice[y] = tmp;
        }
    }
}

void isingLattice::glauber2DimSweep(double beta, double h=0.0){
    std::uniform_int_distribution<> lDist(0, N);
    // to do use lDist and not typewriter
    for(int x=0; x<N; ++x){
        if (rDist(mt19937Engine) <= 1.0/(1.0+exp(-2.0*siteEnergy2D(x, h)*beta))){
            lattice[x] = -lattice[x];
        }
    }
}

#endif