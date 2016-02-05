#include "yaml-cpp/yaml.h"
#include "main.h"
#include "conductivity.h"
// This file is created & removed by the makefile. It is used to keep track 
// of the revision number used to generate a given datafile.
#include "revision.h"

// ------------------------------------------------------------------- \\

class Simulation
{
    int polynom_order;
    double *cT, *cP, *cg; // polynomial coefficients for temperature, pressure, gravity
    double *gravity, *pressure; // profiles; won't change with time

    double kmax; // tracks largest diffusivity to control dt
    double *T, *T_new;
    double *K;

    void write_params_to_file(FILE*);
    void read_profiles(string);
    void initialize(void);
    void iterate(double);

    double get_diffusivity(int);
    double calculate_heat_flow(void);

    double gradient_forward(int);
    double gradient_backward(int);
    double gradient_adiabat(int); 

    int is_convective(int);

    double thermal_diffusion(int);

public:
    void load_config(string);
    void run();
};
   
// ------------------------------------------------------------------- \\

void Simulation::read_profiles(string s)
// read polynome coefficients from dat/polynomes_prefix.dat
// this assumes that the coefficients are listed with higher
// order first.
{
    cout << "Reading profiles from " << s << "..." << endl;
    string line;
    polynom_order = -1; // header line should not be counted
    ifstream data_file (s.c_str());

    if (data_file.is_open())
    {
        // count lines to get highest order
        while (getline(data_file, line)) { polynom_order++; }
        cT = (double*)malloc(polynom_order*sizeof(double));
        cP = (double*)malloc(polynom_order*sizeof(double));
        cg = (double*)malloc(polynom_order*sizeof(double));
        cout << " ... polynom order = " << polynom_order << endl;

        // the file pointer has to be set back to origin before reading the data
        data_file.clear();
        data_file.seekg(0);

        int i=0;
        while (getline(data_file, line))
        {
            // not very elegant but i don't know how to split a line better
            istringstream iss(line);
            string sub1, sub2, sub3;
            iss >> sub1 >> sub2 >> sub3;

            if (sub1 != "#") // we don't want python commented line
            {
                cT[i] = stof(sub1);
                cP[i] = stof(sub2);
                cg[i] = stof(sub3);
                i++;
            }
        }
    } else {
        cout << "Error reading data file! Please use {earth, mars, venus or vesta} as argument." << endl;
        exit(0);
    }
}

int Simulation::is_convective(int i)
// check if the local gradient is super-adiabatic and returns state
{
    if (0.5*(gradient_forward(i)+gradient_backward(i)) > gradient_adiabat(i)) { return 1; }
    else { return 0; }
}

double Simulation::get_diffusivity(int i)
// Returns diffusivity scaled by the core radius.
{
    if (constant_diff) return constant_diff_value/pow(R,2);
    else return diffusivity(pressure[i], T[i])/pow(R,2);
}

void Simulation::initialize(void)
{
    gravity = (double*)malloc(num_points*sizeof(double));
    pressure = (double*)malloc(num_points*sizeof(double));
    T = (double*)malloc(num_points*sizeof(double));
    T_new = (double*)malloc(num_points*sizeof(double));
    K = (double*)malloc(num_points*sizeof(double));

    for (int i=0; i<num_points; i++)
    {
        double rad = i*dx*R;
        T_new[i] = 0.;
        gravity[i] = 0.;
        pressure[i] = 0.;
        // use the polynomial coefficients to initialize the different profiles
        // gravity and pressure won't change with time.
        for(int j=0; j<polynom_order; j++) 
        {
            T_new[i] += cT[j]*pow(rad, polynom_order-j-1);
            gravity[i] += cg[j]*pow(rad, polynom_order-j-1);
            pressure[i] += cP[j]*pow(rad, polynom_order-j-1);
        }
        T[i] = T_new[i];
        K[i] = get_diffusivity(i);
    }
    kmax = get_diffusivity(num_points-1);
}

// Short-hand functions to compute temperature gradients.
double Simulation::gradient_forward(int x) {return (T[x]-T[x+1])/dx;}
double Simulation::gradient_backward(int x) {return (T[x-1]-T[x])/dx;}

double Simulation::gradient_adiabat(int x) 
// Computes local adiabatic gradient
{
    double grad = alpha*gravity[x]*T[x]/cp;
    return grad*R; // has to be normalized by R_core
}

double Simulation::calculate_heat_flow(void)
// calculate CMB heat flow from core temperature, mantle temperature
// and thermal boundary layer thickness.
{
    return tbl_conductivity*(T[num_points-1]-T_mantle)/TBL*4*PI*pow(R,2);
}

double Simulation::thermal_diffusion(int x)
// Computes the right hand side.
{
    double r;
    double rw, re; // west, east radii
    double kw, ke; // west, east thermal diffusivity
    double qw, qe; // west, east prefactor
    double fw, fe; // west, east effective diffusivity
    double gw, ge; // west, east temperature gradient

    r = x*dx;
    rw = r - 0.5*dx;
    re = r + 0.5*dx;

    // this is not used at x=0 or x=num_points-1 (see below)
    kw = 2*K[x-1]*K[x]/(K[x-1]+K[x]);
    ke = 2*K[x+1]*K[x]/(K[x+1]+K[x]);

    gw = gradient_backward(x);
    ge = gradient_forward(x);

    qw = kw*gw;
    qe = ke*ge;

    fw = 0;
    fe = 0;
 
    if (is_convective(x))
    { 
        // effective diffusivity from Kimura et al. 2009
        double Taw, Tae;
        Taw = 0.5*(gradient_adiabat(x)+gradient_adiabat(x-1));
        Tae = 0.5*(gradient_adiabat(x+1)+gradient_adiabat(x));

        double dT_backward = max(0., gw - Taw)/R;
        double dT_forward = max(0., ge - Tae)/R;
        double dist = (num_points-x)*dx*R;

        double kc_backward = min(1e-3, alpha*gravity[x]*pow(dist,4)/(18*mu)*dT_backward)/pow(R,2);
        double kc_forward = min(1e-3, alpha*gravity[x]*pow(dist,4)/(18*mu)*dT_forward)/pow(R,2);
        double temp = max(kc_backward, kc_forward);
        if (temp > kmax) kmax = temp;

        fw = kc_backward*dT_backward*R; // all non-dimensional
        fe = kc_forward*dT_forward*R;
    } else {
        // if at least a point is not convective, we are not "full_convecting"
        if (x != 0) full_convecting = false;
    }

    // no flux at r=0
    if (x==0) return 6*K[0]*(T[1]-T[0])/pow(dx,2);

    // heat flow at r=rcmb; relatively arbitrary for now, but depends on temperature
    if (x==num_points-1) 
    {
        double heat = calculate_heat_flow()/(4*PI*pow(R,3)*11e6*K[x]);
        return (2*x*T[x-1]-2*x*T[x]-(1+x)*2*dx*heat)*K[x]/x/pow(dx,2);
    }

    return ((qw+fw)*pow(rw,2) - (qe+fe)*pow(re,2))/dx/pow(r,2);
}

void Simulation::iterate(double time)
{
    double dT;

    for (int i=0; i<num_points; i++)
    {
        dT = dt*thermal_diffusion(i);
        T_new[i] = T[i] + dT;
    }

    double kmax = 0;
    for (int i=0; i<num_points; i++)
    {
        T[i] = T_new[i];
        K[i] = get_diffusivity(i); 
        if (K[i] > kmax) kmax = K[i];
    }
    // making sure we use a proper timestep for the new diffusivity distribution
    dt = 0.01*pow(dx,2)/kmax;

}

void Simulation::write_params_to_file(FILE *f)
{
    fprintf(f, "# code revision: %s\n", revision.c_str());
    fprintf(f, "# target body: %s\n", body.c_str());
    fprintf(f, "# num_points: %d\n", num_points);

    fprintf(f, "# core radius: %.3e\n", R);
    fprintf(f, "# TBL thickness: %.3e\n", TBL);

    fprintf(f, "# TBL conductivity: %.f\n", tbl_conductivity);
    fprintf(f, "# Kinematic visc: %.3e\n", mu);
    fprintf(f, "# Mantle temperature: %.f\n", T_mantle);

    fprintf(f, "# Constant diffusivity: %d\n", int(constant_diff));
    if (constant_diff)
        fprintf(f, "# Diffusivity value: %.3e\n", constant_diff_value);

    fprintf(f, "#\n");
    return;
}

void Simulation::run()
// argument "prefix" will be prepended to the datafiles names.
// the number of snapshots to output is set in main.h
{
    double time = 0;
    int last_out = 0;

    cout << "Initializes simulation using revision " << revision << "..." << endl;

    ostringstream fname1;
    fname1 << "dat/polynomes_" << body << ".dat";
    read_profiles(fname1.str());

    initialize();

    ostringstream fname2;
    fname2 << "out/" << run_name << "-output.txt";
    FILE *f = fopen(fname2.str().c_str(), "w"); // make sure we don't append to an old file
    write_params_to_file(f);
    fclose(f);

    cout << "Iterates..." << endl;
    while (time < tmax)
    {
        iterate(time);
        time += dt;

        // is it time to output a snapshot?
        if (time >= snapshot*(1+last_out))
        {
            cout << " ... snapshot at t = " << time/Ma << " Ma" << endl;
            FILE *f = fopen(fname2.str().c_str(), "a");

            for (int x=0;x<num_points;x++) fprintf(f, "%.9g %.9g %.9g %.9g %.9g %.9g %d\n", x*dx, time/Ma, T[x], pow(R,2)*K[x], gradient_adiabat(x)/R, calculate_heat_flow()/1e12, is_convective(x));
            fprintf(f,"\n");
            fclose(f);

            last_out++;
        }

        if (full_convecting)
        {
            cout << "Full convecting!" << endl;
            exit(0);
        }
        full_convecting = true; // starts true, then checks if at least one point is not
    }
    cout << "Done!..." << endl;
}
    
void Simulation::load_config(string fname)
{
    YAML::Node config = YAML::LoadFile(fname);
    run_name = config["run_name"].as<string>();
    body = config["body"].as<string>();

    num_points = config["num_points"].as<int>();
    dx = 1./num_points;

    R = config["core_radius"].as<int>()*1e3;
    TBL = config["tbl_thickness"].as<double>()*1e3;

    snapshot = config["snapshot"].as<int>()*Ma;
    tmax = config["tmax"].as<double>()*Ma;

    tbl_conductivity = config["tbl_conductivity"].as<double>();
    mu = config["kinematic_visc"].as<double>();
    T_mantle = config["mantle_temperature"].as<double>();

    constant_diff = config["constant_diff"].as<bool>();
    constant_diff_value = config["constant_diff_value"].as<double>();
}

int main(int argc, char **argv)
{
    Simulation s;
    s.load_config("config.yaml");
    s.run();
}
