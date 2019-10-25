#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>

std::string data_to_name(const int&, const int&, const int&, const double&, const std::string&);
std::string cc_to_name(const double&);
void name_to_data(const std::string&, int&, int&, int&, double&, const std::string&);
int n_meas(const int& n_tot, const int& n_gap);

#endif

