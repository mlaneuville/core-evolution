#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <stdlib.h>

#include <cmath>
#include <cassert>

using namespace std;

double murnaghan(double P, double T);
double volume_ratio(double P, double T);
double conductivity(double P, double T);
double integrand(double z);
double resistivity_temperature(double T, double f);
double resistivity_volume(double f, string element);
