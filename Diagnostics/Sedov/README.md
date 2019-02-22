# Sedov diagnostics

Process a Sedov problem to produce rho, u and p as a function of r, for
comparison to the analytic solution.

## Building & running

The n-dimensional diagnostic can be built by executing `make DIM=n`. This will
produce the executable `sedov_nd.exe`. To run, the executable must be provided
with a list of arguments. For all problems, you must provide the name of the
plotfile to be analyzed and the name of the slicefile where the results shall
be output:
```
./sedov_nd.exe -p plotfile_name -s slicefile_name
```

Additional arguments depend on the problem:

- **1d**: no additional arguments are required
- **2d cylindrical**: extra arguents `--xctr x` and `--yctr y` giving the coordinates
of the domain center (x,y) can be provided (both default to 0.0)
- **2d spherical**: the argument `--sphr` *must* be provided to indicate a spherical
problem, and the `--yctr y` argument can be provided to give the y coordinate of
the domain center (defaults to 0.0)
- **3d cylindrical**: extra arguents `--xctr x` and `--yctr y` giving the coordinates
of the domain center (x,y) can be provided (both default to 0.0)
- **3d spherical**: the argument `--sphr` *must* be provided to indicate a spherical
problem, and extra arguents `--xctr x`, `--yctr y` and `--zctr z` giving the coordinates
of the domain center (x,y,z) can be provided (all default to 0.0)
