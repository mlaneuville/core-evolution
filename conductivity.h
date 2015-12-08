#include <string>
#include <cmath>

using namespace std;

double murnaghan(double P, double T);
double volume_ratio(double P, double T);
double conductivity(double P, double T);
double integrand(double z);
double resistivity_temperature(double T, double f);
double resistivity_volume(double f, string element);
void print_kmap(void);