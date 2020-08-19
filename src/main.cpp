#include "yaml-cpp/yaml.h"
#include "main.h"
#include "../revision.h"

void Simulation::load_config(string fname)
{
    if(!ifstream(fname)) {
        cout << "Config file not found. Please copy default.yaml to ";
        cout << "config.yaml to start." << endl;
        exit(1);
    }

    core.load_config(fname);
    mantle.load_config(fname);

    YAML::Node config = YAML::LoadFile(fname);
    string run_name = config["run_name"].as<string>();

    ostringstream fn;
    fn<< "out/" << run_name << "-core-output.txt";
    FILE *f = fopen(fn.str().c_str(), "w"); // make sure we don't append to an old file
    core.write_params_to_file(f);
    fclose(f);

    return;
}

void Simulation::initialize(void)
{
    core.initialize();
    //mantle.initialize(core.get_tcmb());
    mantle.initialize(2250);

    return;
}

double Simulation::get_tcmb(void) {
    if(coupled) { return core.get_tcmb(); }
    else { return mantle.get_tcmb(); }
}

void Simulation::output(void)
{
    mantle.output(time);
    if(coupled) { core.output(time); }
}

void Simulation::run(void)
{
    double YEAR = 365*24*3600;
    double tmax = 500e6*YEAR;
    double dt = 1e4*YEAR;
    int last_out = 0;
    double snapshot = 1e7*YEAR;

    time = 0;
    coupled = true;
    while(time < tmax) {
        // base mantle evolution model
        mantle.runstep(time, dt, get_tcmb());
        time += dt;

        // if coupled, run core evolution based on qcmb
        if(coupled) {
            core.runstep(mantle.getcoreheatflow());
            dt = core.gettimestep();
        }

        // produce output every dt = snapshot
        if(time >= snapshot*(1+last_out)) {
            cout << " ... snapshot at " << time/YEAR << endl;
            output();
            last_out++;
        }
    }
}

int main(int argc, char **argv){
    Simulation s;

    s.load_config("config.yaml");
    s.initialize();
    s.run();
    return 0;
}
