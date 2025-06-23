#ifndef ISING2D_H
#define ISING2D_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include"processing.h"
#include"randutils.hpp"

class isingLattice2D{
    private:
        const int L, N;
        int *lattice; // lattice
        int *nn; // neighbours
        double probArray[5]={0};
        std::uniform_int_distribution<> lDist;
        std::uniform_real_distribution<> rDist;
        std::mt19937 *mt19937Engine;
        void neighboutList(void){
            for (int i = 0; i < N; ++i)
            {
                nn[(4*i)+0] = (mod(i,L)==0) ?  i+L-1 : i-1;
                nn[(4*i)+1] = (mod(i+1,L)==0) ? i+1-L : i+1;
                nn[(4*i)+2] = mod(i-L,N);
                nn[(4*i)+3] = mod(i+L,N);
            }
        };
    public:
        isingLattice2D(int l, std::mt19937 *rng): L(l), N(l*l), mt19937Engine(rng){
            rDist = std::uniform_real_distribution<>(0.0, 1.0);
            lDist = std::uniform_int_distribution<>(0, N-1);
            lattice = new int[N];
            nn = new int[4*N];
            neighboutList();
        };

        void initialise(double m0){
            for(int i=0; i<N; ++i)
                lattice[i] = (rDist(*mt19937Engine) < (1.0+m0)/2.0)? 1 : -1;
        }
        void boltzGlauber(double b, double h=0.0){
            for (int i = 0; i < 5; ++i) probArray[i] = 1.0/(1.0+exp(2.0*b*((i*2.0)-4.0+h)));
        };
        void metropolisSweep(double, double);
        void metropolisPUDSweep(double, double);
        void glauberSweep(double, double);
        void glauberInterval(double, int);
        double magnetisation(void) const;
        double energy(double) const;
        double siteEnergy(int, double) const;
        double siteEnergyPUD(int, double) const;
        void writeConfig(std::ostream &data) const;
        // void flipSeriesOfSites(int, int);
        ~isingLattice2D() {
            delete[] lattice;
            delete[] nn;
        };
};

void isingLattice2D::writeConfig(std::ostream &data) const{
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

double isingLattice2D::magnetisation(void) const{
    double m=0.0;
    for (int i = 0; i < N; ++i) m += lattice[i];
    return m/N;
}

double isingLattice2D::siteEnergy(int x, double J=1.0) const{
    double energy = lattice[nn[(4*x)+0]] + lattice[nn[(4*x)+1]]
                + lattice[nn[(4*x)+2]] + lattice[nn[(4*x)+3]];
    return J*energy;
}

double isingLattice2D::energy(double k=0.0) const{
    double e = 0.0;
    for(int x=0; x< N; ++x){
        e += lattice[x]*siteEnergy(x);
    }
    return -e/(2*N);
}

double isingLattice2D::siteEnergyPUD(int x, double J) const{
    if (mod(x,2)==0) J=1.0;
    double energy = lattice[nn[(4*x)+0]] + lattice[nn[(4*x)+1]]
                +J*(lattice[nn[(4*x)+2]] + lattice[nn[(4*x)+3]]);
    return energy;
}


void isingLattice2D::metropolisSweep(double beta, double J){
    double deltaE;
    for(int i=0; i<N; ++i){
        int x = lDist(*mt19937Engine);
        deltaE = 2.0*lattice[x]*siteEnergy(x, J);
        if (deltaE <=0 || rDist(*mt19937Engine) < exp(-beta*deltaE)){
            lattice[x] = -lattice[x];
        }
    }
}

void isingLattice2D::metropolisPUDSweep(double beta, double J){
    double deltaE;
    for(int i=0; i<N; ++i){
        deltaE = 2.0*lattice[i]*siteEnergyPUD(i, J);
        if (deltaE <=0 || rDist(*mt19937Engine) < exp(-beta*deltaE)){
            lattice[i] = -lattice[i];
        }
    }
}

void isingLattice2D::glauberSweep(double beta, double h=0.0){
    int ide, x;
    boltzGlauber(beta, h);
    for(int i=0; i<N; ++i){
        x = lDist(*mt19937Engine);
        ide = static_cast<int>(lattice[x]*siteEnergy(x)/2.0+2.0);
        // ide = lattice[x]*siteEnergy(x);
        // if (rDist(*mt19937Engine) <= 1.0/(1.0+exp(2.0*ide*beta))){
        if (rDist(*mt19937Engine) <= probArray[ide]){
            lattice[x] = -lattice[x];
        }
    }
}

// glauber random update over a given number of sites
void isingLattice2D::glauberInterval(double beta, int nSites){
    int ide, x;
    boltzGlauber(beta, 0.0);
    for(int i=0; i<nSites; ++i){
        x = lDist(*mt19937Engine);
        ide = static_cast<int>(lattice[x]*siteEnergy(x)/2.0+2.0);
        // ide = lattice[x]*siteEnergy(x);
        // if (rDist(*mt19937Engine) <= 1.0/(1.0+exp(2.0*ide*beta))){
        if (rDist(*mt19937Engine) <= probArray[ide]){
            lattice[x] = -lattice[x];
        }
    }
}

#endif