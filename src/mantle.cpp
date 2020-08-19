#include "yaml-cpp/yaml.h"
#include "mantle.h"
#include "fd1d_heat_steady.hpp"

bool DEBUG = false;
double YEAR = 365*24*3600;

// ----------------------------------------

void Mantle::load_config(string fname)
{
    YAML::Node config = YAML::LoadFile(fname);
    run_name = config["run_name"].as<string>();
    body = config["body"].as<string>();

    YAML::Node bodyconf = config[body];
    rmars = bodyconf["radius"].as<double>();
    rcore = bodyconf["core_radius"].as<double>();

    YAML::Node mantleconf = config["mantle"];
    rhom = mantleconf["density"].as<double>();
    cm = mantleconf["specific_heat"].as<double>();
    km = mantleconf["cond_mantle"].as<double>();
    kc = mantleconf["cond_crust"].as<double>();
    alpha = mantleconf["thermal_expansivity"].as<double>();
    g = mantleconf["gravity"].as<double>();
    eta0 = mantleconf["ref_viscosity"].as<double>();
    idtcore = mantleconf["initial_dtcore"].as<double>();
}

void Mantle::initialize(double tc)
{
    // initial conditions
    tcore = tc;
    tbottom = tc-idtcore;
    tmantle = tbottom/(1+alpha*g*(rmars-rcore)/cm)-600;
    tlithos = tmantle-100;
    rlithos = 3300e3;
    tsurface = 250;
    dt = 3.0381e10; //1e5*YEAR;

    // model setting
    acore = 4*pow(rcore, 2)*M_PI;
    vcore = 4*pow(rcore, 3)*M_PI/3;
    rhocpcore = 7200*840;
    rhocpman = 3400*1000; //rhom*cm;
    gc = 1.5; // ?
    DELTA = 2.21;
    diffusivity = 1e-6; //km/(rhom*cm);
    fc = 1;
    A = 3e5;
    R = 8.3144;

    tref = 1600;
    racrit = 450;
    ra = 1e6;
    rai = 1e6;
    raicrit = 450;

    return;
}

double Mantle::viscosity(double temp)
// Mantle viscosity for a given temperature.
{
    return eta0*exp(A*(tref-temp)/R/tref/temp);
}

void Mantle::updaterayleighs(void)
// Update the different Rayleigh numbers.
{
    double dT;
    double prefactor = alpha*rhom*g/diffusivity/viscosity(tmantle);

    dT = tmantle - tlithos + tcore - tbottom;
    ra = prefactor*dT*pow(rlithos-rcore, 3);

    dT = tmantle - tsurface + tcore - tbottom;
    rai = prefactor*dT*pow(rmars-rcore, 3);
    raicrit = 0.28*pow(rai, 0.21);

    if (DEBUG) { cout << "ra = " << ra << ", rai = " << rai << ", raic = " << raicrit << endl; }
    return;
}

void Mantle::gettlithos(void)
// Temperature at the base of the lithosphere.
{
    tlithos = tmantle - DELTA*R*pow(tmantle, 2)/A;
    return;
}

void Mantle::gettbottom(void)
// Temperature at the base of the mantle.
{
    double tbltop = topboundarylayer();
    double tblbot = botboundarylayer();
    assert(tbltop > 0);
    //assert(tblbot > 0);
    tblbot = max(0., tblbot);
    double dR = rlithos - rcore - tbltop - tblbot;
    tbottom = tmantle*(1 + alpha*g*dR/cm);
    return;
}

double Mantle::topboundarylayer(void)
// Upper boundary layer thickness according to Turcotte and Schubert.
{
    return (rlithos - rcore)*pow(racrit/ra, 1./3);
}

double Mantle::botboundarylayer(void)
// Lower boundary layer thickness according to Turcotte and Schubert.
{
    double nc = viscosity(0.5*(tcore+tbottom));
    return pow(diffusivity*fc*nc*raicrit/(alpha*rhom*gc*(tcore-tbottom)), 1./3);
}

double Mantle::lithosheatflow(void)
// Computes heat flow to the base of the lithosphere for a given mantle state.
{
    double tbltop = topboundarylayer();
    double flux = km*(tmantle - tlithos)/tbltop;
    assert(flux > 0);
    return flux;
}

double Mantle::radioactive(double time)
// in W/m3
{
    double t = 4.5e9 - time/YEAR;

    double cu = 16e-9;
    double cth = 56e-9;
    double ck = 305e-6;

    double cst;
    cst = cu*0.9928*9.46e-5*exp(t*log(2.0)/4.47e9);
    cst += cu*0.0071*5.69e-4*exp(t*log(2.0)/0.704e9);
    cst += cth*2.64e-5*exp(t*log(2.0)/14.0e9);
    cst += ck*1.19e-4*2.92e-5*exp(t*log(2.0)/1.25e9);

    //double tau = 2.5e9*YEAR;
    //double flux = exp(-time/tau)*8.5e-8;
    double flux = 3400.*cst;

    assert(flux > 0);
    return flux;
}

double Mantle::getcoreheatflow(void)
{
    return coreheatflow()*4*M_PI*pow(rcore, 2);
}

double Mantle::coreheatflow(void)
// Computes core heat flow for a given mantle state.
{
    double tblbot = botboundarylayer();
    //double flux = km*(tcore - tbottom)/tblbot;
    double flux = km*(tcore - 2000)/50e3;
    //assert(flux > 0);
    return max(0., flux);
}

double k2(double x) {return 3.;}
double f2(double x) {return 0.;}

double Mantle::lithosdiffusion(void)
// Computes the thermal gradient at the base of the lithosphere. Boundary
// conditions are tlithos and tsurface, respectively. The conductivity changes
// between mantle and crust, so does radioactive heat sources content.
{
    double a, b;
    double ua, ub;
    double *x, *u;
    int n = 51;

    double dx = (rmars - rlithos)/n;
    a = tlithos;
    b = tsurface;

    ua = tlithos;
    ub = tsurface;

    x = r8vec_even(n, a, b);
    u = fd1d_heat_steady(n, a, b, ua, ub, k2, f2, x);
    double q = 3*(u[0] - u[1])/dx;
    delete x;
    delete u;

    return q;
}

double Mantle::energylithos(double ql)
// Given heat flux from the convecting mantle, computes lithospheric growth.
//  - needs ql in W
//  - returns dlithos in km
{
    double condflux = lithosdiffusion();
    double dd = (ql-condflux)/rhocpman/(tmantle-tlithos);
    return dd;
}

double Mantle::energycore(double qc)
// Given a core mantle boundary heat flow, compute core cooling.
//  - needs qc in W
//  - returns dT in K
{
    double dT = qc*acore/(rhocpcore*vcore);
    return dT;
}

double Mantle::energymantle(double qc, double ql, double qr)
// Conservation of energy in the mantle.
{
    vmantle = 4*M_PI*(pow(rlithos, 3) - pow(rcore, 3))/3;
    double prefactor = rhocpman*vmantle;
    return (-ql*4*M_PI*pow(rlithos, 2) + qc*acore + qr*vmantle)/prefactor;
}

// ----------------------------------------

void Mantle::runtests(void)
// Some tests to make sure things behave as expected.
{
    cout << "=== RUNNING TESTS ===" << endl;
    assert(viscosity(1600.) == eta0);
    assert(viscosity(1700.) < eta0);
    assert(viscosity(1500.) > eta0);
    cout << "=== TESTS COMPLETED ===" << endl;
}


void print_header(void) {
    cout << "time,";
    cout << "tlithos,";
    cout << "tmantle,";
    cout << "tbottom,";
    cout << "tcore,";
    cout << "qc,";
    cout << "ql,";
    cout << "qr,";
    cout << "rlithos,";
    cout << "bbl,";
    cout << "tbl" << endl;
}

double Mantle::get_tcmb(void) { return tcore; }

// ----------------------------------------

void Mantle::runstep(double time, double ts, double tc) {

    dt = ts;
    tcore = tc;
    double dtcore;
    double dtman, drlithos;

    // compute fluxes for this timestep prep
    qc = coreheatflow();
    ql = lithosheatflow();
    qr = radioactive(time);

    // core update
    dtcore = energycore(qc)*dt;
    tcore -= dtcore;
    assert(tcore > 0);

    // mantle update
    dtman = energymantle(qc, ql, qr)*dt;
    tmantle += dtman;
    assert(tmantle > 0);

    // lithosphere update
    drlithos = energylithos(ql)*dt;
    rlithos += drlithos;
    assert(rlithos > rcore);
    assert(rlithos < rmars);

    // update variables
    updaterayleighs();
    gettlithos();
    gettbottom();

    return;
}

void Mantle::output(double time)
{
    ostringstream fname;
    fname << "out/" << run_name << "-mantle-output.txt";

    FILE *f = fopen(fname.str().c_str(), "a");
    fprintf(f, "%.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",
            time, tlithos, tmantle, tbottom, tcore, qc*dt, ql*dt, qr*vmantle,
            rlithos, botboundarylayer(), topboundarylayer());
    fclose(f);
}

void Mantle::runsim(void)
{
    double TMAX = 4.5e9*YEAR;

    print_header();

    int i = 0;
    while(i*dt < TMAX) {
        runstep(i*dt, dt, tcore);
        i++;
    }

    return;
}

// ----------------------------------------
