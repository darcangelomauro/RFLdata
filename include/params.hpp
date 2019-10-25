#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <iostream>
#include <string>

// Struct that packs the simulation
// parameters all together (hmc and mmc).
// 
//
// L    =   # of integration steps in leapfrog
// 
// AR   =   target acceptance rate
// 
// dAR  =   randomize dt in interval [dt_min, dt_max]
//          such that dt_min gives AR+dAR and
//          dt_max gives AR-dAR
//
//
// iter_duav    =   # of dual averaging iterations
// 
// iter         =   # of iterations per sample
//
// samples      =   # of samples to be collected
//
//
// Explanation of MC mode:
//
// fix          =   HMC with fixed stepsize
//
// rand         =   HMC with random stepsize
//
// mmc          =   MMC
//
struct Simul_params
{
    // Geometric parameters
    int p;
    int q;
    int dim;

    // Integration length (HMC only)
    int L;

    // Acceptance rate
    double AR;
    double dAR=0;

    // Simulation parameters
    int iter_duav;
    int iter;
    int samples;

    // MC mode
    std::string mode;

    // Control string
    std::string control;
};

// Function to read simulation parameters from stream
bool read_init_stream(std::istream&, struct Simul_params&);

// Checks that the necessary parameters are there
bool params_validity(const struct Simul_params&);

#endif

