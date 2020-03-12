#include <iostream>
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

        clock_t start = clock();
        G.HMC_duav(Nt, dt, iter, engine, 0.65);
        clock_t end = clock();

        cout << "Time: " << double(end-start)/CLOCKS_PER_SEC << " sec" << endl;
    }
    else if(mode == "mmc")
    {
        int iter;
        double scale;
        cout << "insert scale, iter" << endl;
        cin >> scale >> iter;

        clock_t start = clock();
        G.MMC_duav(scale, iter, engine, 0.232);
        clock_t end = clock();

        cout << "Time: " << double(end-start)/CLOCKS_PER_SEC << " sec" << endl;
    }
    else
        cout << "mode not recognized" << endl;

    // Free memory
    gsl_rng_free(engine);

    return 0;
}
