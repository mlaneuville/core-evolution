Solves the diffusion equation in a sphere with variable diffusivity.

### How-to

* Run `make` to prepare the executable file.
* Then `./a.out <optional_run_name>` will run the simulation.
* The output will be `<optional_run_name->output-XXX.txt`, where XXX is output number.
* Each file contains a temperature and diffusivity profile at a given time.

The number of outputs can be set in `main.h`

A way to visualize a temperature profile in gnuplot for instance is `plot "output-0.txt" u 1:3 w l`.


### plot.py

to run it: `python plot.py N name`
* N is the number of plot you want (by default, it would be only the
first file output-0.txt which will be used)
* name is the `<optional_run_name>` (without the dash) if used to
obtain the .txt files.

so for
`./a.out test`
with 3 outputs (the default parameters in main.h)
you should use
`python plot.py 3 test`
