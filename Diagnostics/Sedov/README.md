# Sedov diagnostics

Process a Sedov problem to produce rho, u and p as a function of r,
for comparison to the analytic solution.

## Building & running

The n-dimensional diagnostic can be built by executing `make
DIM=n`. This will produce the executable `sedov_nd.exe`. To run, the
executable must be provided with a list of arguments. For all
problems, you must provide the name of the plotfile to be analyzed and
the output file name:

```
./sedov_nd.ex -p plotfile_name -s output_name
```
