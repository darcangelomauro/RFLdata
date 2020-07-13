#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <armadillo>
#include <gsl/gsl_rng.h>
#include "geometry.hpp"
#include "utils.hpp"

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));
    
    
    // create geometry
    cout << "insert p, q, dim, g2" << endl;
    Geom24 G(cin); 
    G.shuffle(engine);
    
    int p = G.get_p();
    int q = G.get_q();
    int dim = G.get_dim();
    double g2 = G.get_g2();

    string name = data_to_name(p, q, dim, g2, "GEOM");
    ofstream out_S(name + "_S.txt");
    ofstream out_HL(name+ "_HL.txt");
    out_S.precision(16);
    out_HL.precision(16);
    
        
    // thermalize first
    double tgt = 0.8;
    double dt = 0.000001;
    G.HMC_duav(10, dt, 10000, engine, tgt, "leapfrog");
    cout << "dual averaging complete" << endl;
    cout << "dual averaged dt: " << dt << endl;
    double ar = G.HMC(10, dt, 10000, engine, "leapfrog");
    cout << "thermalization complete" << endl;
    cout << "acceptance rate: " << ar << endl;
    
    
    for(int i=1; i<1000; ++i)
    {
        G.HMC(10, dt, 1000, engine, "leapfrog");
        G.print_S(out_S);
        G.print_HL(out_HL);
    }

    out_S.close();
    out_HL.close();
}
