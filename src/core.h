#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <stdlib.h>

#include <cmath>
#include <cassert>

using namespace std;

class Core
{
    string run_name;
    string body;

    double tmax; // duration of the run - will move to command line argument at some point

    int num_points; // grid size
    bool constant_diff;
    double constant_diff_value;
    double mu; // kinematic viscosity
    double cp;

    double YEAR = 365*24*3600;
    double MA = 1e6*YEAR;
    double G = 6.67e-11;

    double dx = 1/num_points;
    double dt = 0.01 * pow(dx, 2) / 1e-5; // will be updated as soon as actual diffusivity is computed
    double R; // core size [m]
    //  volume of the top shell, for dt calculation
    double shell_volume;

    bool full_convecting = true;
    int polynom_order;
    double *cT, *cP, *cg; // polynomial coefficients for temperature, pressure, gravity
    double *gravity, *pressure; // profiles; won't change with time

    double kmax; // tracks largest diffusivity to control dt
    double *T, *T_new;
    double *K;
    double heatflow;

    void read_profiles(string);

    double get_diffusivity(int);
    double calculate_heat_flow(void);

    double gradient_forward(int);
    double gradient_backward(int);
    double gradient_adiabat(int);

    int is_convective(int);

    double thermal_diffusion(int);

public:
    double get_tcmb(void);
    double gettimestep(void);
    void runstep(double);
    void load_config(string);
    void output(double);
    void run();
    void write_params_to_file(FILE *);
    void initialize(void);
};
