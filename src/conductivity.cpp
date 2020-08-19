#include "conductivity.h"

double murnaghan(double P, double T)
// Uses Murnaghan EOS to provide rho(P,T)
{
    double rho = rho0 * pow(1 + Kp0 * P / K0, 1. / Kp0);
    double dtemp = exp(-alpha * exp(-alphap * P / K0) * (T - T0));
    return rho * dtemp;
}

double volume_ratio(double P, double T)
// V/V0 = rho0/rho
// note: STP conditions should be computed only once
{
    return murnaghan(0, 300) / murnaghan(P, T);
}

double diffusivity(double P, double T)
// Returns diffusivity using rho(P,T) and constant cp
{
    return conductivity(P, T) / murnaghan(P, T) / cp;
}

double conductivity(double P, double T)
// From Gomi et al. (2013); Eq 11-13
{
    double concentration_Si = 22.5;
    double f = volume_ratio(P, T);

    // ideal
    double ideal;
    ideal = resistivity_temperature(T, f) * resistivity_volume(f, "Fe");
    ideal += resistivity_volume(f, "Si") * concentration_Si;

    // saturation
    double saturation = 1.68e-6 * pow(f, 1. / 3);

    // total
    double L = 2.44e-8;
    return L * T * (1. / ideal + 1. / saturation);
}

double integrand(double z)
{
    return pow(z, 5) / (exp(z) - 1) / (1 - exp(-z));
}


double simpsons(double(*f)(double x), double a, double b, int n)
// integral using simpsons algorithm
{
    double h = (b - a) / n;
    double x;
    double r;
    char m = 0;
    double s = 0.0;

    for(x = a; x <= b; x += h) {
        r = f(x);
        if(x == a || x == b) {
            s += r;
        } else {
            m = !m;
            s += r * (m + 1) * 2.0;
        }
    }
    return s * (h / 3.0);
}


double resistivity_temperature(double T, double f)
// Eq 1 from Gomi et al 2013
{
    double theta0 = 417;
    double gamma = 1.52;

    double theta = theta0 * exp(-gamma * log(f)); // debye temperature
    double integral = simpsons(integrand, 1e-5, theta / T, 100);

    return pow(T / theta, 5) * integral;
}

double resistivity_volume(double f, string element)
// Eq 9 from Gomi et al 2013
// f = V/V0
{
    double F1, F2, F3;
    if(element.compare("Fe") == 0) {
        // for hcp-iron
        F1 = 5.26e-9; // ohm m
        F2 = 1.24;
        F3 = -3.21;
    }
    if(element.compare("Si") == 0) {
        // silicium
        F1 = 3.77e-8;
        F2 = 1.48;
        F3 = -3.1;
    }

    return F1 * pow(F2 - f, F3);
}

void print_kmap(void)
// Produces a data file with k = k(T,P) - mostly for debugging purposes
{
    FILE *f = fopen("kmap.txt", "w");

    double pressure, temperature;
    for(int x = 0; x <= 100; x++) {
        pressure = 130e9 + x * 320e9 / 100;
        for(int y = 0; y <= 100; y++) {
            temperature = 2000 + y * 3000. / 100;
            fprintf(f, "%.9g %.9g %.9g\n", pressure, temperature, conductivity(pressure, temperature));
        }
        fprintf(f, "\n");
    }
    fclose(f);
}
