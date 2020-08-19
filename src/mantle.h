#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

class Mantle
{
    string run_name;
    string body;

    double dt, idtcore;
    double tlithos, tbottom, tref, tsurface, tmantle, tcore; // K
    double rmars, rlithos, rcore;
    double dcrust, dlithos; // km
    double rhocpcore, rhocpman; // J/K/m3
    double vcore, vmantle; // m3
    double acore; // m2
    double rhom;
    double alpha, g, gc, cm, kc, km, diffusivity, fc, eta0;
    double ra, rai, racrit, raicrit;
    double R, A, DELTA;
    double qc, qr, ql;

    void updaterayleighs(void);
    double viscosity(double);
    double botboundarylayer(void);
    double topboundarylayer(void);

    void gettlithos(void);
    void gettbottom(void);

    double coreheatflow(void);
    double lithosheatflow(void);
    double radioactive(double);

    double lithosdiffusion(void);

    double energycore(double);
    double energymantle(double, double, double);
    double energylithos(double);

public:
    double get_tcmb(void);
    double getcoreheatflow(void);
    void runstep(double, double, double);
    void runtests(void);
    void initialize(double);
    void runsim(void);
    void output(double);
    void load_config(string);
};
