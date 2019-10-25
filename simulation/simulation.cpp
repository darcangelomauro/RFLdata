#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sys/stat.h>
#include <ctime>
#include <armadillo>
#include <gsl/gsl_rng.h>
#include "geometry.hpp"
#include "utils.hpp"
#include "params.hpp"

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    // ARGUMENT LIST:
    // The first argument is time(NULL) (it will be used as global seed).
    // The second argument is the job rank (it will be used as local seed).
    // The third argument is the g2 value.
    // The fourth argument is the path to the working directory
    if(argc != 5)
    {
        cerr << "Error: need to pass global seed, local rank, g2 value, and path to working directory as arguments." << endl;
        return 0;
    }
    // First two arguments are converted to unsigned long, the third to double.
    unsigned long global_seed = stoul(argv[1]);
    unsigned long local_rank = stoul(argv[2]);
    double g2 = stod(argv[3]);


    // STRING OPERATIONS FOR FILE NAMES
    string path = argv[4];
    string path_g2 = path + "/" + cc_to_name(g2);
    string path_g2_rank = path_g2 + "/" + argv[2];
    string init_filename = path_g2_rank + "/init.tmp"; 


    // RNG:
    // Initialize random number generator with global+local
    // (to avoid identical initialization among different jobs).
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd2);
    gsl_rng_set(engine, global_seed+local_rank);
    clog << "RNG seed: " << global_seed+local_rank << endl << endl;
    
    
    //********* BEGIN PARAMETER INITIALIZATION **********//
    struct Simul_params sm;
    ifstream in_init;
    
    // Check that parameter reading is successful
    in_init.open(init_filename);
    if(!read_init_stream(in_init, sm))
    {
        cerr << "Error: couldn't read file " + init_filename << endl;
        return 1;
    }
    in_init.close();

    // Log what has been read
    clog << "File " + init_filename + " contains the following parameters:" << endl;
    clog << sm.control << endl;
    
    // Check that parameters are consistent with simulation type
    if(!params_validity(sm))
    {
        cerr << "Error: file " + init_filename + " does not contain the necessary parameters." << endl;
        return 1;
    }
    //********* END PARAMETER INITIALIZATION **********//



    //********* BEGIN MONTE CARLO **********//

    // Create geometry
    Geom24 G(sm.p, sm.q, sm.dim, g2);
    G.shuffle(engine);
    clog << G << endl;

    // Open output files
    ofstream out_s, out_hl;
    string s_filename = path_g2_rank + "/" + data_to_name(sm.p, sm.q, sm.dim, g2, "GEOM") + "_S.txt";
    string hl_filename = path_g2_rank + "/" + data_to_name(sm.p, sm.q, sm.dim, g2, "GEOM") + "_HL.txt";
    out_s.open(s_filename);
    out_hl.open(hl_filename);

    if(sm.mode == "fix")
    {
        // DUAL AVERAGING
        double dt = 0.005;
        clog << "Duav averaging start timestamp: " << time(NULL) << endl;
        G.HMC_duav(sm.L, dt, sm.iter_duav, engine, sm.AR);
        clog << "Dual averaging end timestamp: " << time(NULL) << endl;

        // SIMULATION
        clog << "Simulation start timestamp: " << time(NULL) << endl;
        double ar = 0;
        for(int i=0; i<sm.samples; ++i)
        {
            ar += G.HMC(sm.L, dt, sm.iter, engine);
            G.print_S(out_s);
            G.print_HL(out_hl);
        }
        clog << "Simulation end timestamp: " << time(NULL) << endl;
        clog << "Integration step: " << dt << endl;
        clog << "Acceptance rate: " << ar/sm.samples << endl;
    }

    else if(sm.mode == "rand")
    {
        // DUAL AVERAGING
        double dt_min = 0.005;
        double dt_max = 0.005;
        clog << "Duav averaging start timestamp: " << time(NULL) << endl;
        G.HMC_duav(sm.L, dt_min, sm.iter_duav, engine, sm.AR + sm.dAR);
        G.HMC_duav(sm.L, dt_max, sm.iter_duav, engine, sm.AR - sm.dAR);
        clog << "Dual averaging end timestamp: " << time(NULL) << endl;

        // SIMULATION
        clog << "Simulation start timestamp: " << time(NULL) << endl;
        double ar = 0;
        for(int i=0; i<sm.samples; ++i)
        {
            ar += G.HMC(sm.L, dt_min, dt_max, sm.iter, engine);
            G.print_S(out_s);
            G.print_HL(out_hl);
        }
        clog << "Simulation end timestamp: " << time(NULL) << endl;
        clog << "Integration step: " << dt_min << " " << dt_max << endl;
        clog << "Acceptance rate: " << ar/sm.samples << endl;
    }
    
    else if(sm.mode == "mmc")
    {
        // DUAL AVERAGING
        double scale = 0.005;
        clog << "Duav averaging start timestamp: " << time(NULL) << endl;
        G.MMC_duav(scale, sm.iter_duav, engine, sm.AR);
        clog << "Dual averaging end timestamp: " << time(NULL) << endl;

        // SIMULATION
        clog << "Simulation start timestamp: " << time(NULL) << endl;
        double ar = 0;
        for(int i=0; i<sm.samples; ++i)
        {
            ar += G.MMC(scale, sm.iter, engine);
            G.print_S(out_s);
            G.print_HL(out_hl);
        }
        clog << "Simulation end timestamp: " << time(NULL) << endl;
        clog << "Metropolis scale: " << scale << endl;
        clog << "Acceptance rate: " << ar/sm.samples << endl;
    }

    else
        cerr << "Error: mode not recognized" << endl;


    // Close output files
    out_s.close();
    out_hl.close();
    
    //********* END MONTE CARLO **********//


    // Free memory
    gsl_rng_free(engine);

    return 0;
}
