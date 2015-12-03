#include "main.h"
#include "conductivity.h"
#include "revision.h"

class Simulation
{
    double *T, *T_new;
    double *K;

    void initialize(void);
    void iterate(double);

    double get_diffusivity(int);

    double gradient_forward(int);
    double gradient_backward(int);

    double thermal_diffusion(int);

public:
    void run(string);
    
};

double Simulation::get_diffusivity(int i)
{
    double pressure = 100e9; // TODO: of course this will need to change
    double rhocp = 7e6; // TODO: change with real value at some point 
    return conductivity(pressure, T[i])/rhocp/pow(R,2);
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
    // the CMB is kept at 2500 K
    T_new[num_points-1] = 2500;
    T[num_points-1] = 2500;
}

double Simulation::gradient_forward(int x) {return (T[x]-T[x+1])/dx;}
double Simulation::gradient_backward(int x) {return (T[x-1]-T[x])/dx;}

double Simulation::thermal_diffusion(int x)
{
    double radius = dx*x;
    double next_radius = dx*(x+1);
    double prev_radius = dx*(x-1);

    // no flux at r=0
    if (x==0) return 6*K[0]*(T[1]-T[0])/pow(dx,2);
    // fixed temperature at r=rcmb
    if (x==num_points-1) return 0.;

    // k = k(r,T)
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

    for (int i=0; i<num_points; i++)
    {
        T[i] = T_new[i];
        K[i] = get_diffusivity(i); 
    }


}

void Simulation::run(string prefix)
{
    double time = 0;
    int last_out = 0;
    string prepend = "output-";

    cout << "Initializes simulation using revision " << revision << "..." << endl;
    initialize();

    cout << "Iterates..." << endl;
    while (time < tmax)
    {
        iterate(time);
        time += dt;

        if (time >= snapshot*(1+last_out))
        {
            ostringstream fname;
            fname << prefix << "output-" << last_out << ".txt";
            FILE *f = fopen(fname.str().c_str(), "w");

            for (int x=0;x<num_points;x++) fprintf(f, "%.9g %.9g %.9g %.9g\n", x*dx, time/Ma, T[x], K[x]/k0);
            fprintf(f,"\n");
            fclose(f);

            last_out++;
        }
    }
    cout << "Done!..." << endl;
}

int main(int argc, char **argv)
{
    /*
    FILE *f = fopen("kmap.txt", "w");

    double pressure, temperature;
    for (int x=0; x<=100; x++)
    {
        pressure = 130e9 + x*320e9/100;
        for (int y=0; y<=100; y++) { temperature = 2000 + y*3000./100; fprintf(f, "%.9g %.9g %.9g\n", pressure, temperature, conductivity(pressure,temperature)); }
        fprintf(f,"\n");
    }
    fclose(f);
    */

    string prefix = "";
    if (argc != 1) prefix = argv[1] + string("-");

    Simulation s;
    s.run(prefix);
}
