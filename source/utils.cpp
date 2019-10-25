#include <sstream>
#include <iostream>
#include <algorithm>
#include <ctime>
#include "utils.hpp"
#include "geometry.hpp"
#include "clifford.hpp"

using namespace std;


string data_to_name(const int& p, const int& q, const int& dim, const double& g2, const string& prefix)
{
    string sp = to_string(p);
    string sq = to_string(q);
    string sdim = to_string(dim);
    
    ostringstream osg2;
    osg2 << g2;
    string sg2 = osg2.str();
    replace(sg2.begin(), sg2.end(), '.', 'd');
    replace(sg2.begin(), sg2.end(), '-', 'n');
    
    return prefix + "p" + sp + "q" + sq + "dim" + sdim + "g" + sg2;
}

string cc_to_name(const double& g2)
{
    ostringstream osg2;
    osg2 << g2;
    string sg2 = osg2.str();
    replace(sg2.begin(), sg2.end(), '.', 'd');
    replace(sg2.begin(), sg2.end(), '-', 'n');
    
    return sg2;
}

void name_to_data(const string& s, int& p, int& q, int& dim, double& g2, const string& prefix)
{
    // generate a clean string s1 from s that
    // ignores everything before the first appearence
    // of prefix (if found)
    size_t start = s.find(prefix);
    string s1;
    if(start != string::npos)
        s1 = s.substr(start);
    else
        s1 = s;

    // extract data from string.
    // data is formatted as follows:
    //
    // GEOMp[int]q[int]dim[int]g[char][int]d[int][garbage].txt
    // 
    // where [int] denotes an integer numerical value and [char] a single character.
    // the last part g[char][int]d[int] is a double written as
    // integer part and decimal part separated by a letter d, with a 'n'
    // in front if g2 is negative.
    
    string s_p = s1.substr(s1.find("p")+1, s1.find("q")-1);
    string s_q = s1.substr(s1.find("q")+1, s1.find("d")-s1.find("q")-1);
    string s_dim = s1.substr(s1.find("m")+1, s1.find("g")-s1.find("m")-1);
    string s_g = s1.substr(s1.find("g")+1, s1.find(".")-s1.find("g")-1);
    replace(s_g.begin(), s_g.end(), 'd', '.');
    replace(s_g.begin(), s_g.end(), 'n', '-');

    p = stoi(s_p);
    q = stoi(s_q);
    dim = stoi(s_dim);
    g2 = stod(s_g);
}

int n_meas(const int& n_tot, const int& gap)
{
    int res = 0;

    for(int i=0; i<n_tot; ++i)
    {
        if( !(i%gap) )
            ++res;
    }

    return res;
}
