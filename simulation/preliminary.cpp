#include <iostream>
#include <ctime>
#include <chrono>
#include <armadillo>
#include <gsl/gsl_rng.h>
#include "geometry.hpp"
#include "utils.hpp"
#include "params.hpp"

using namespace std;
using namespace std::chrono;
using namespace arma;

int main(int argc, char** argv)
{
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));
    
    
    // create geometry
    cout << "insert p, q, dim, g2" << endl;
    Geom24 G(cin); 
    G.shuffle(engine);

    string mode;
    cout << "choose mode (hmc or mmc)" << endl;
    cin >> mode;


    if(mode == "hmc")
    {
        int Nt, iter;
        double dt;
        cout << "insert Nt, dt, iter" << endl;
        cin >> Nt >> dt >> iter;

        auto start = high_resolution_clock::now();
        G.HMC_duav(Nt, dt, iter, engine, 0.65);
        auto end = high_resolution_clock::now();
        
        auto duration = duration_cast<milliseconds>(end-start);

        cout << "Time: " << duration.count()/1000. << " sec" << endl;
    }
    else if(mode == "mmc")
    {
        int iter;
        double scale;
        cout << "insert scale, iter" << endl;
        cin >> scale >> iter;

        auto start = high_resolution_clock::now();
        G.MMC_duav(scale, iter, engine, 0.232);
        auto end = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(end-start);
        
        cout << "Time: " << duration.count()/1000. << " sec" << endl;
    }
    else
        cout << "mode not recognized" << endl;

    // Free memory
    gsl_rng_free(engine);

    return 0;
}
