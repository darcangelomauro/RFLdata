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
    // The first argument is the seed.
    // The second argument is the g2 value.
    // The third argument is the path to the work directory
    if(argc != 4)
    {
        cerr << "Error: need to pass seed, g2 value, and path to work directory as arguments." << endl;
        return 0;
    }
    // First argument (seed) is converted to unsigned long, the second (g2) to double.
    unsigned long seed = stoul(argv[1]);
    double g2 = stod(argv[2]);


    // STRING OPERATIONS FOR FILE NAMES
    string path_lvl0 = argv[3];
    string path_lvl1 = path_lvl0 + "/" + cc_to_name(g2);
    string init_filename = path_lvl0 + "/init.txt"; 


    // RNG:
    // Initialize random number generator
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd2);
    gsl_rng_set(engine, seed);
    clog << "RNG seed: " << seed << endl << endl;
    
    
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

    // CREATE GEOMETRY
    Geom24 G(sm.p, sm.q, sm.dim, g2);
    clog << G << endl;


    // OPEN OUTPUT FILES:
    // start.txt stores the final H,L configuration
    // scalar.txt stores whatever scalar parameter was computed with dual averaging
    ofstream out_start, out_scalar;
    string start_filename = path_lvl1 + "/" + "start.txt";
    string scalar_filename = path_lvl1 + "/" + "scalar.txt";
    out_start.open(start_filename);
    out_scalar.open(scalar_filename);

    if(sm.mode == "fix")
    {
        double dt = 0.005;
        clog << "Duav averaging start timestamp: " << time(NULL) << endl;
        G.shuffle(engine);
        G.HMC_duav(sm.L, dt, sm.iter_duav, engine, sm.AR);
        clog << "Dual averaging end timestamp: " << time(NULL) << endl;
        clog << "Integration step: " << dt << endl;

        out_scalar << dt;
        G.print_HL(out_start);
    }

    else if(sm.mode == "rand")
    {
        double dt_min = 0.005;
        double dt_max = 0.005;
        clog << "Duav averaging start timestamp: " << time(NULL) << endl;
        G.shuffle(engine);
        G.HMC_duav(sm.L, dt_min, sm.iter_duav, engine, sm.AR + sm.dAR);
        G.shuffle(engine);
        G.HMC_duav(sm.L, dt_max, sm.iter_duav, engine, sm.AR - sm.dAR);
        clog << "Dual averaging end timestamp: " << time(NULL) << endl;
        clog << "Integration step: " << dt_min << " " << dt_max << endl;
        
        out_scalar << dt_min << " " << dt_max;
        G.print_HL(out_start);
    }
    
    else if(sm.mode == "mmc")
    {
        double scale = 0.005;
        clog << "Duav averaging start timestamp: " << time(NULL) << endl;
        G.shuffle(engine);
        G.MMC_duav(scale, sm.iter_duav, engine, sm.AR);
        clog << "Dual averaging end timestamp: " << time(NULL) << endl;
        clog << "Metropolis scale: " << scale << endl;
        
        out_scalar << scale;
        G.print_HL(out_start);
    }

    else
        cerr << "Error: mode not recognized" << endl;


    // Close output files
    out_start.close();
    out_scalar.close();
    
    //********* END MONTE CARLO **********//


    // Free memory
    gsl_rng_free(engine);

    return 0;
}
