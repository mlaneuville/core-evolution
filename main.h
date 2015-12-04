#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <stdlib.h>

#include <cmath>
#include <cassert>

using namespace std;

double PI = 3.14159265359;
double G = 6.67e-11;

int num_points = 200;

double R = 3480e3;
double k0 = 1e-5/pow(R,2);
double dx = 1./num_points;
double dt = 0.25*pow(dx,2)/k0;

double Ma = 1e6*365*24*3600.;
double tmax = 1000*Ma;

// the denominator determines how many snapshot will be written
double snapshot = tmax/3.;
