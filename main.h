#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <stdlib.h>

#include <cmath>
#include <cassert>

using namespace std;

const double rho = 11e3; 

double PI = 3.14159265359;
double G = 6.67e-11;
double Ma = 1e6*365*24*3600.;

double tmax = 1000*Ma; // duration of the run - will move to command line argument at some point

int num_points = 200; // grid size

double R = 3480e3; // core size [m]
double k0 = 1e-5/pow(R,2); // thermal diffusivity scaled with core radius

double dx = 1./num_points;
double dt = 0.25*pow(dx,2)/k0;

// what is the duration between snapshots?
double snapshot = 10*Ma;
