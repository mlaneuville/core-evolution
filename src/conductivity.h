#include <string>
#include <cmath>

using namespace std;

// EOS parameters - need references
const double T0 = 1812;
const double rho0 = 7010;
const double K0 = 130e9;
const double Kp0 = 4;
const double alpha = 1e-5 / exp(-0.400485 * 135 / 130);
const double alphap = 130 * log(2.) / (360 - 135);
const double cp = 750; // Gubbins et al 2003


double murnaghan(double P, double T);
double volume_ratio(double P, double T);
double conductivity(double P, double T);
double diffusivity(double P, double T);
double integrand(double z);
double resistivity_temperature(double T, double f);
double resistivity_volume(double f, string element);
void print_kmap(void);
