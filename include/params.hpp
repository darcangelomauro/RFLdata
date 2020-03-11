#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <iostream>
#include <string>

// Struct that packs the simulation
// parameters all together (hmc and mmc).
// 
// L    =   # of integration steps in leapfrog
// 
// AR   =   target acceptance rate
// 
// iter_duav    =   # of dual averaging iterations
//
// iter_therm   =   # of thermalization iterations
// 
// iter_sim     =   # of iterations per sample
//
// samples      =   # of samples to be collected
//
//
// Available MC modes:
// hmc
// mmc
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

    // Simulation parameters
    int iter_duav;
    int iter_therm;
    int iter_sim;
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

