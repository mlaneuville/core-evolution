#include "main.h"
#include "revision.h"

class Simulation
{
    double *T, *T_new;
    double *K; // not used for now, constant K

    void initialize(void);
    void iterate(double);

    double get_diffusivity(int);

    double gradient_forward(int);
    double gradient_backward(int);

    double thermal_diffusion(int);

public:
    void run();
    
};

double Simulation::get_diffusivity(int i)
{
    // k(T) = k0/2 @ 2500K, K0 at 4500K, linear in between
    // (purely arbitrary, I know)
//    return k0;
    return 0.5*k0 + (T[i]-2500)/2000.*0.5*k0;
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

    // constant diffusivity
    // return k0*((x+1)*T[x+1]-2*x*T[x]+(x-1)*T[x-1])/pow(dx,2)/x;
    
    // k = k(r,T)
    double diffusion = K[x]*((x+1)*T[x+1]-2*x*T[x]+(x-1)*T[x-1])/pow(dx,2)/x;
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

void Simulation::run(void)
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
            fname << "output-" << last_out << ".txt";
            FILE *f = fopen(fname.str().c_str(), "w");

            for (int x=0;x<num_points;x++) fprintf(f, "%.9g %.9g %.9g\n", x*dx, time/Ma, T[x]);
            fprintf(f,"\n");
            fclose(f);

            last_out++;
        }
    }
    cout << "Done!..." << endl;
}

int main(int argc, char **argv)
{
    Simulation s;
    s.run();
}
