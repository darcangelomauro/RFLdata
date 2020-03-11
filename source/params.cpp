#include <iostream>
#include <string>
#include "params.hpp"

using namespace std;

bool read_init_stream(istream& in, struct Simul_params& sm)
{
    bool success = false;
    if(in)
    {
        string temp;
        
        while(in >> temp)
        {
            if(temp == "p:")
            {
                sm.control += temp;
                in >> sm.p;
            }
            else if(temp == "q:")
            {
                sm.control += temp;
                in >> sm.q;
            }
            else if(temp == "dim:")
            {
                sm.control += temp;
                in >> sm.dim;
            }
            else if(temp == "L:")
            {
                sm.control += temp;
                in >> sm.L;
            }
            else if(temp == "AR:")
            {
                sm.control += temp;
                in >> sm.AR;
            }
            else if(temp == "iter_duav:")
            {
                sm.control += temp;
                in >> sm.iter_duav;
            }
            else if(temp == "iter_therm:")
            {
                sm.control += temp;
                in >> sm.iter_therm;
            }
            else if(temp == "iter_sim:")
            {
                sm.control += temp;
                in >> sm.iter_sim;
            }
            else if(temp == "samples:")
            {
                sm.control += temp;
                in >> sm.samples;
            }
            else if(temp == "mode:")
            {
                sm.control += temp;
                in >> sm.mode;
            }
        }

        success = true;

    }

    return success;
}

bool params_validity(const struct Simul_params& sm)
{
    if(sm.control.find("p:") == std::string::npos)
    {
        cerr << "p not found" << endl;
        return 0;
    }
    if(sm.control.find("q:") == std::string::npos)
    {
        cerr << "q not found" << endl;
        return 0;
    }
    if(sm.control.find("dim:") == std::string::npos)
    {
        cerr << "dim not found" << endl;
        return 0;
    }
    
    if(sm.control.find("mode:") == std::string::npos)
    {
        cerr << "mode not found" << endl;
        return 0;
    }
    if(sm.control.find("iter_duav:") == std::string::npos)
    {
        cerr << "iter_duav not found" << endl;
        return 0;
    }
    if(sm.control.find("iter_therm:") == std::string::npos)
    {
        cerr << "iter_therm not found" << endl;
        return 0;
    }
    if(sm.control.find("iter_sim:") == std::string::npos)
    {
        cerr << "iter_sim not found" << endl;
        return 0;
    }

    if(sm.mode == "hmc")
    {
        if(sm.control.find("L:") == std::string::npos)
        {
            cerr << "L not found" << endl;
            return 0;
        }
        if(sm.control.find("AR:") == std::string::npos)
        {
            cerr << "AR not found" << endl;
            return 0;
        }
    }
    
    else if(sm.mode == "mmc")
    {
        if(sm.control.find("AR:") == std::string::npos)
        {
            cerr << "AR not found" << endl;
            return 0;
        }
    }

    return 1;
}




