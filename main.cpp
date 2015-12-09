#include "main.h"
#include "conductivity.h"
// This file is created & removed by the makefile. It is used to keep track 
// of the revision number used to generate a given datafile.
#include "revision.h" 

// ------------------------------------------------------------------- \\

class Simulation
{
    double *T, *T_new;
    double *K;

    void initialize(void);
    void iterate(double);

    double get_diffusivity(int);

    double gradient_forward(int);
    double gradient_backward(int);
    double gradient_adiabat(int); 

    int is_convective(int);

    double thermal_diffusion(int);

public:
    void run(string);
};
   
// ------------------------------------------------------------------- \\

int Simulation::is_convective(int i)
// check if the local gradient is super-adiabatic and returns state
{
    if (0.5*(gradient_forward(i)+gradient_backward(i)) > gradient_adiabat(i)) { return 1; }
    else { return 0; }
}

double Simulation::get_diffusivity(int i)
// Returns diffusivity scaled by the core radius.
// TODO: diffusivity takes into account variable density, but the estimation
// of pressure from the radius does not so far.
{
    double pressure = 363.85e9 - 2*PI*pow(i*dx,2)*G*pow(rho,2)/3;
    return diffusivity(pressure, T[i])/pow(R,2);
}

void Simulation::initialize(void)
{
    T = (double*)malloc(num_points*sizeof(double));
    T_new = (double*)malloc(num_points*sizeof(double));
    K = (double*)malloc(num_points*sizeof(double));

    for (int i=0; i<num_points; i++)
    {
        // initialize with a linear gradient from 2000 to 4500 K
        T_new[i] = 2000 + i*2500./num_points;
        T[i] = T_new[i];
        K[i] = get_diffusivity(i);
    }
}

// Short-hand functions to compute temperature gradients.
double Simulation::gradient_forward(int x) {return (T[x]-T[x+1])/dx;}
double Simulation::gradient_backward(int x) {return (T[x-1]-T[x])/dx;}

double Simulation::gradient_adiabat(int x) 
// Computes local adiabatic gradient
{
    // memo: alpha is defined in conductivity.h 
    double gravity = 4*PI*G*rho*R*x/num_points/3; // see Eq. 6 of Labrosse 2015 for higher order
    double grad = alpha*gravity*T[x]/cp;
    return grad/R; // has to be normalized by R_core
}

double Simulation::thermal_diffusion(int x)
// Computes the right hand side.
{
    double radius = dx*x;
    double next_radius = dx*(x+1);
    double prev_radius = dx*(x-1);  

    if (is_convective(x))
    { 
        // effective diffusivity here
    }
    

    // no flux at r=0
    if (x==0) return 6*K[0]*(T[1]-T[0])/pow(dx,2);

    // fixed heat flow at r=rcmb
    if (x==num_points-1) 
    {
        double heat = 1e13/(4*PI*pow(R,3)*11e6*K[x]);
        return (2*x*T[x-1]-2*x*T[x]-(1+x)*2*dx*heat)*K[x]/x/pow(dx,2);
    }

    double diffusion;
    diffusion = K[x]*((x+1)*T[x+1]-2*x*T[x]+(x-1)*T[x-1])/pow(dx,2)/x;
    diffusion += (K[x+1]-K[x-1])*(T[x+1]-T[x-1])/4/pow(dx,2);

    return diffusion;
}

void Simulation::iterate(double time)
{
    double dT;

    for (int i=0; i<num_points; i++)
    {
        dT = dt*thermal_diffusion(i);
        T_new[i] = T[i] + dT;
    }

    double dTa = gradient_adiabat(num_points-10);
    
    double kmax = 0;
    for (int i=0; i<num_points; i++)
    {
        T[i] = T_new[i];
        K[i] = get_diffusivity(i); 
        if (K[i] > kmax) kmax = K[i];
    }
    // making sure we use a proper timestep for the new diffusivity distribution
    dt = 0.25*pow(dx,2)/kmax;

}

void Simulation::run(string prefix)
// argument "prefix" will be prepended to the datafiles names.
// the number of snapshots to output is set in main.h
{
    double time = 0;
    int last_out = 0;
    string prepend = "output-";

    cout << "Initializes simulation using revision " << revision << "..." << endl;
    initialize();

    ostringstream fname;
    fname << prefix << "output-" << last_out << ".txt";
    FILE *f = fopen(fname.str().c_str(), "w"); // make sure we don't append to an old file
    fclose(f);

    cout << "Iterates..." << endl;
    while (time < tmax)
    {
        iterate(time);
        time += dt;

        // is it time to output a snapshot?
        if (time >= snapshot*(1+last_out))
        {
            FILE *f = fopen(fname.str().c_str(), "a");

            for (int x=0;x<num_points;x++) fprintf(f, "%.9g %.9g %.9g %.9g %d\n", x*dx, time/Ma, T[x], K[x]/k0, is_convective(x));
            fprintf(f,"\n");
            fclose(f);

            last_out++;
        }
    }
    cout << "Done!..." << endl;
}

int main(int argc, char **argv)
{
    // uncomment to produce kmap.txt which shows k = k(T,P)
    //print_kmap();

    string prefix = "";
    if (argc != 1) prefix = argv[1] + string("-");

    Simulation s;
    s.run(prefix);
}
