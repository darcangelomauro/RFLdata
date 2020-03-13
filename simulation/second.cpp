#include <iostream>
#include <fstream>
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
    // The fourth argument is the job index
    if(argc != 5)
    {
        cerr << "Error: need to pass seed, g2 value, path to work directory, and job index as arguments." << endl;
        return 0;
    }
    // First argument (seed) is converted to unsigned long, the second (g2) to double.
    unsigned long seed = stoul(argv[1]);
    double g2 = stod(argv[2]);


    // STRING OPERATIONS FOR FILE NAMES
    string path_lvl0 = argv[3];
    string path_lvl1 = path_lvl0 + "/" + cc_to_name(g2);
    string path_lvl2 = path_lvl1 + "/" + argv[4];
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

    // OPEN OUTPUT FILES
    ofstream out_s, out_hl;
    string s_filename = path_lvl2 + "/" + data_to_name(sm.p, sm.q, sm.dim, g2, "GEOM") + "_S.txt";
    string hl_filename = path_lvl2 + "/" + data_to_name(sm.p, sm.q, sm.dim, g2, "GEOM") + "_HL.txt";
    out_s.open(s_filename);
    out_hl.open(hl_filename);
    


    // OPEN INPUT FILES ======================
    ifstream in_scalar, in_start;
    
    string scalar_filename = path_lvl1 + "/" + "scalar.txt";
    string start_filename = path_lvl1 + "/" + "start.txt";
    
    in_scalar.open(scalar_filename);
    if(!in_scalar)
    {
        cerr << "Error: no scalar file was found" << endl;
        return 1;
    }

    in_start.open(start_filename);
    if(!in_start)
    {
        cerr << "Error: no start file was found" << endl;
        return 1;
    }
    // OPEN INPUT FILES ======================



    if(sm.mode == "hmc")
    {
        // Read input data
        G.read_mat(in_start);
        double dt;
        in_scalar >> dt;
        in_start.close();
        in_scalar.close();

        // Start simulation
        clog << "Simulation start timestamp: " << time(NULL) << endl;
        double ar = 0;
        for(int i=0; i<sm.samples; ++i)
        {
            ar += G.HMC(sm.L, dt, sm.iter_sim, engine);
            G.print_S(out_s);
            G.print_HL(out_hl);
        }
        clog << "Simulation end timestamp: " << time(NULL) << endl;
        clog << "Integration step: " << dt << endl;
        clog << "Acceptance rate: " << ar/sm.samples << endl;
    }
    
    else if(sm.mode == "mmc")
    {
        // Read input data
        G.read_mat(in_start);
        double scale;
        in_scalar >> scale;
        in_start.close();
        in_scalar.close();

        // Start simulation
        clog << "Simulation start timestamp: " << time(NULL) << endl;
        double ar = 0;
        for(int i=0; i<sm.samples; ++i)
        {
            ar += G.MMC(scale, sm.iter_sim, engine);
            G.adjust();
            G.print_S(out_s);
            G.print_HL(out_hl);
        }
        clog << "Simulation end timestamp: " << time(NULL) << endl;
        clog << "Metropolis scale: " << scale << endl;
        clog << "Acceptance rate: " << ar/sm.samples << endl;
    }

    else
        cerr << "Error: mode not recognized" << endl;


    // Close files
    out_s.close();
    out_hl.close();
    
    //********* END MONTE CARLO **********//


    // Free memory
    gsl_rng_free(engine);

    return 0;
}
