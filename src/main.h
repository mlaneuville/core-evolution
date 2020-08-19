#include <string>
#include <sstream>

#include "core.h"
#include "mantle.h"

class Simulation
{
    Mantle mantle;
    Core core;

    double time;
    bool coupled;

    void output(void);
    double get_tcmb(void);

public:
    void load_config(string);
    void initialize(void);
    void run(void);
};
