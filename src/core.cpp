#include "yaml-cpp/yaml.h"
#include "core.h"
#include "conductivity.h"

// ------------------------------------------------------------------- \\

void Core::load_config(string fname)
{
    if(!ifstream(fname)) {
        cout << "Config file not found. Please copy default.yaml to config.yaml to start." << endl;
        exit(0);
    }

    YAML::Node config = YAML::LoadFile(fname);
    run_name = config["run_name"].as<string>();
    body = config["body"].as<string>();
    tmax = config["tmax"].as<double>()*MA;

    YAML::Node coreconfig = config["core"];
    num_points = coreconfig["num_points"].as<int>();
    dx = 1./num_points;

    YAML::Node bodyconfig = config[body];
    R = bodyconfig["core_radius"].as<double>();

    mu = coreconfig["kinematic_visc"].as<double>();
    cp = coreconfig["specific_heat"].as<double>();
    constant_diff = coreconfig["constant_diff"].as<bool>();
    constant_diff_value = coreconfig["constant_diff_value"].as<double>();
}

void Core::read_profiles(string s)
// read polynome coefficients from dat/polynomes_prefix.dat
// this assumes that the coefficients are listed with higher
// order first.
{
    cout << "Reading profiles from " << s << "..." << endl;
    string line;
    polynom_order = -1; // header line should not be counted
    ifstream data_file(s.c_str());

    if(data_file.is_open()) {
        // count lines to get highest order
        while(getline(data_file, line)) {
            polynom_order++;
        }
        cT = (double *)malloc(polynom_order * sizeof(double));
        cP = (double *)malloc(polynom_order * sizeof(double));
        cg = (double *)malloc(polynom_order * sizeof(double));
        cout << " ... polynom order = " << polynom_order << endl;

        // the file pointer has to be set back to origin before reading the data
        data_file.clear();
        data_file.seekg(0);

        int i = 0;
        while(getline(data_file, line)) {
            // not very elegant but i don't know how to split a line better
            istringstream iss(line);
            string sub1, sub2, sub3;
            iss >> sub1 >> sub2 >> sub3;

            if(sub1 != "#") { // we don't want python commented line
                cT[i] = stof(sub1);
                cP[i] = stof(sub2);
                cg[i] = stof(sub3);
                i++;
            }
        }
        cout << "Done reading from file" << endl;
    } else {
        cout << "Error reading data file! Please use {earth, mars, venus or vesta} as argument." << endl;
        exit(0);
    }
}

void Core::initialize(void)
{
    shell_volume = 4*M_PI*(pow(R, 3)-pow(R - R*dx, 3))/3;

    ostringstream fname1;
    fname1 << "dat/polynomes_" << body << ".dat";
    read_profiles(fname1.str());

    cout << "Initializing..." << endl;
    gravity = (double *)malloc(num_points * sizeof(double));
    pressure = (double *)malloc(num_points * sizeof(double));
    T = (double *)malloc(num_points * sizeof(double));
    T_new = (double *)malloc(num_points * sizeof(double));
    K = (double *)malloc(num_points * sizeof(double));

    for(int i = 0; i < num_points; i++) {
        double rad = i * dx * R;
        T_new[i] = 0.;
        gravity[i] = 0.;
        pressure[i] = 0.;
        // use the polynomial coefficients to initialize the different profiles
        // gravity and pressure won't change with time.
        for(int j = 0; j < polynom_order; j++) {
            T_new[i] += cT[j] * pow(rad, polynom_order - j - 1);
            gravity[i] += cg[j] * pow(rad, polynom_order - j - 1);
            pressure[i] += cP[j] * pow(rad, polynom_order - j - 1);
        }
        T[i] = T_new[i];
        T[i] = 2880;
        T_new[i] = 2880;
        K[i] = get_diffusivity(i);
    }
    kmax = get_diffusivity(num_points - 1);
    cout << "Done initializing!" << endl;
}

int Core::is_convective(int i)
// check if the local gradient is super-adiabatic and returns state
{
    if(0.5 * (gradient_forward(i) + gradient_backward(i)) > gradient_adiabat(i)) {
        return 1;
    } else {
        return 0;
    }
}

double Core::get_diffusivity(int i)
// Returns diffusivity scaled by the core radius.
{
    if(constant_diff) {
        return constant_diff_value / pow(R, 2);
    } else {
        return diffusivity(pressure[i], T[i]) / pow(R, 2);
    }
}

double Core::get_tcmb(void) { return T[num_points-1]; }
double Core::gettimestep(void) { return dt; }

// Short-hand functions to compute temperature gradients.
double Core::gradient_forward(int x) { return (T[x] - T[x + 1]) / dx; }
double Core::gradient_backward(int x) { return (T[x - 1] - T[x]) / dx; }

double Core::gradient_adiabat(int x)
// Computes local adiabatic gradient
{
    double grad = alpha * gravity[x] * T[x] / cp;
    return grad * R; // has to be normalized by R_core
}

double Core::calculate_heat_flow(void)
// calculate CMB heat flow from core temperature, mantle temperature
// and thermal boundary layer thickness.
{
    //return mantle.getcoreheatflow();
    return heatflow;
}

double Core::thermal_diffusion(int x)
// Computes the right hand side.
{
    double r;
    double rw, re; // west, east radii
    double kw, ke; // west, east thermal diffusivity
    double qw, qe; // west, east prefactor
    double fw, fe; // west, east effective diffusivity
    double gw, ge; // west, east temperature gradient

    r = x * dx;
    rw = r - 0.5 * dx;
    re = r + 0.5 * dx;

    // this is not used at x=0 or x=num_points-1 (see below)
    kw = 2 * K[x - 1] * K[x] / (K[x - 1] + K[x]);
    ke = 2 * K[x + 1] * K[x] / (K[x + 1] + K[x]);

    gw = gradient_backward(x);
    ge = gradient_forward(x);

    qw = kw * gw;
    qe = ke * ge;

    fw = 0;
    fe = 0;

    if(is_convective(x)) {
        // effective diffusivity from Kimura et al. 2009
        double Taw, Tae;
        Taw = 0.5 * (gradient_adiabat(x) + gradient_adiabat(x - 1));
        Tae = 0.5 * (gradient_adiabat(x + 1) + gradient_adiabat(x));

        double dT_backward = max(0., gw - Taw) / R;
        double dT_forward = max(0., ge - Tae) / R;
        double dist = (num_points - x) * dx * R;

        double kc_backward = min(1e-3, alpha*gravity[x]*pow(dist, 4)/(18*mu)*dT_backward)/pow(R, 2);
        double kc_forward = min(1e-3, alpha*gravity[x]*pow(dist, 4)/(18*mu)*dT_forward)/pow(R, 2);
        double temp = max(kc_backward, kc_forward);
        if(temp > kmax) {
            kmax = temp;
        }

        fw = kc_backward * dT_backward * R; // all non-dimensional
        fe = kc_forward * dT_forward * R;
    } else {
        // if at least a point is not convective, we are not "full_convecting"
        if(x != 0) {
            full_convecting = false;
        }
    }

    // no flux at r=0
    if(x == 0) {
        return 6 * K[0] * (T[1] - T[0]) / pow(dx, 2);
    }

    // heat flow at r=rcmb; relatively arbitrary for now, but depends on temperature
    if(x == num_points - 1) {
        // TODO: generalize for different bodies
        double heat = calculate_heat_flow() / (4 * M_PI * pow(R, 3) * 11e6 * K[x]);
        return (2 * x * T[x - 1] - 2 * x * T[x] - (1 + x) * 2 * dx * heat) * K[x] / x / pow(dx, 2);
    }

    return ((qw + fw) * pow(rw, 2) - (qe + fe) * pow(re, 2)) / dx / pow(r, 2);
}

void Core::runstep(double qc)
{
    heatflow = qc;
    double dT;

    // making sure we use a proper timestep for the new diffusivity distribution
    dt = 0.01 * pow(dx, 2) / kmax;
    // the heat flow from the top shell to the core should not be more than 10x
    // typical conducting heat flow
    double dt_boundary = dt;

    if(calculate_heat_flow() > 1e8)
        // this condition is valid only when the heat flow is large enough
        // as calculate_heat_flow() tends towards 0, dt would become infinite.
    {
        dt_boundary = 10 * 5e3 * cp * 1e-2 * shell_volume / calculate_heat_flow();
    }

    if(dt_boundary < dt) {
        dt = dt_boundary;
    }
    assert(dt > 0);

    for(int i = 0; i < num_points; i++) {
        dT = dt * thermal_diffusion(i);
        T_new[i] = T[i] + dT;
    }

    double kmax = 0;
    for(int i = 0; i < num_points; i++) {
        T[i] = T_new[i];
        K[i] = get_diffusivity(i);
        if(K[i] > kmax) {
            kmax = K[i];
        }
    }
}

void Core::write_params_to_file(FILE *f)
{
    //fprintf(f, "# 'code revision':'%s'\n", revision.c_str());
    fprintf(f, "# 'target body': '%s'\n", body.c_str());
    fprintf(f, "# 'num_points': %d\n", num_points);

    fprintf(f, "# 'core radius': %.3e\n", R);

    fprintf(f, "# 'Kinematic visc': %.3e\n", mu);

    fprintf(f, "# 'Constant diffusivity': %d\n", int(constant_diff));
    if(constant_diff) {
        fprintf(f, "# 'Diffusivity value': %.3e\n", constant_diff_value);
    }

    fprintf(f, "#\n");
    return;
}

void Core::output(double time)
{
    ostringstream fname;
    fname << "out/" << run_name << "-core-output.txt";

    FILE *f = fopen(fname.str().c_str(), "a");

    for(int x = 0; x < num_points; x++) {
        fprintf(f, "%.9g %.9g %.9g %.9g %.9g %.9g %d\n", R*x*dx/1e3,
                time, T[x], pow(R, 2)*K[x], gradient_adiabat(x)/R,
                heatflow/1e12, is_convective(x));
    }
    fprintf(f, "\n");
    fclose(f);
}
void Core::run()
// argument "prefix" will be prepended to the datafiles names.
// the number of snapshots to output is set in main.h
{
    double time = 0;
    int last_out = 0;

    //cout << "Initializes simulation using revision " << revision << "..." << endl;


    initialize();

    //ostringstream fname2;
    //fname2 << "out/" << subfolder << "/" << run_name << "-output.txt";
    //// TODO: this will segfault if subfolder does not exist!
    //FILE *f = fopen(fname2.str().c_str(), "w"); // make sure we don't append to an old file
    //write_params_to_file(f);
    //fclose(f);

    //cout << "Initialize mantle..." << endl;
    //mantle.initialize(T[num_points-1]);

    cout << "Iterates..." << endl;
    while(time < tmax) {
        //mantle.runstep(time, dt, T[num_points-1]);
        runstep(time);
        time += dt;

        // is it time to output a snapshot?
        //if(time >= snapshot * (1 + last_out)) {
        //    cout << " ... snapshot at t = " << time / Ma << " Ma" << endl;
        //    FILE *f = fopen(fname2.str().c_str(), "a");

        //    for(int x = 0; x < num_points; x++) {
        //        fprintf(f, "%.9g %.9g %.9g %.9g %.9g %.9g %d\n", R*x*dx/1e3,
        //                time/Ma, T[x], pow(R, 2)*K[x], gradient_adiabat(x)/R,
        //                //mantle.getcoreheatflow()/1e12, is_convective(x));
        //                1e12, is_convective(x));
        //    }
        //    fprintf(f, "\n");
        //    fclose(f);
        //
        //    last_out++;
        //}

        if(full_convecting) {
            cout << "Full convecting!" << endl;
            exit(0);
        }
        full_convecting = true; // starts true, then checks if at least one point is not
    }
    cout << "Done!..." << endl;
}
