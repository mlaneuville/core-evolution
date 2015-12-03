#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <stdlib.h>

#include <cmath>
#include <cassert>

using namespace std;

int num_points = 100;

double R = 3480e3;
double k0 = 1e-5/pow(R,2);
double dx = 1./num_points;
double dt = 0.25*pow(dx,2)/k0;

double Ma = 1e6*365*24*3600.;
double tmax = 1000*Ma;

double snapshot = tmax/3.;
