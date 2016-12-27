#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <stdlib.h>

#include <cmath>
#include <cassert>

using namespace std;

string run_name = "testrun";
string body = "earth";
string subfolder = "";
bool full_convecting = true;

bool constant_diff = false;
double constant_diff_value = 1e-5;

double tbl_conductivity = 5; // W/m/K
double mu = 1e5; // kinematic viscosity
double T_mantle = 4500;

double PI = 3.14159265359;
double G = 6.67e-11;
double Ma = 1e6 * 365 * 24 * 3600.;

int num_points = 200; // grid size
double dx = 1. / num_points;

double tmax = 1000 * Ma; // duration of the run - will move to command line argument at some point
double snapshot = 10 * Ma;
double dt = 0.01 * pow(dx, 2) / 1e-5; // will be updated as soon as actual diffusivity is computed

double R = 3480e3; // core size [m]
double TBL = 100e3; // thermal boundary layer thickness
double shell_volume = 4 * PI * (pow(R, 3) - pow(R - R *dx, 3)) / 3; //  volume of the top shell, for dt calculation
