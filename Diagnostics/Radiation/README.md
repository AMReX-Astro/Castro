# Radiation diagnostics

This directory contains diagnostics for Castro's radiation tests.
These are:

- `gaussian_pulse`: Process a 2-d gaussian radiation pulse.
- `lgt_frnt1d`: Process a 1-d (radiative) sedov problem to produce rho, u, and
    p as a function of r, for comparison to the analytic solution.
- `rad_shock`: extract a 1-d slice of the data of a radiative shock test (all
        variables or a single variable) along the specified coordinate direction
        from a plotfile.  The plotfile can be 1-, 2-, or 3-d.
- `rad_source`: Analysis routine for the radiation source test.  Take a list of
    files and print out (rho e) and the total radiation energy density in the
    first zone as a function of time.
- `rad_sphere`: Print out the radiation quantities at a specified distance from
    the origin.  This is written for the 1-d radiating sphere problem.
- `rhd_shocktube`: Analysis routine for RHD_shocktube.

## Building

To build one of the radiation diagnostics, specify the executable name
as the target, e.g.

```
   make DIM=1 rad_sphere.ex
```

Take care that you compile with the correct dimension (i.e. set `DIM`
equal to the same value it had for the code that generated the plotfile
to be investigated).

## Running

Command line arguments must be passed to the executables to provide the plotfile(s)
to be anaylzed and problem parameters:

- `gaussian_pulse`: `-p plotfile_name` (or `--plotfile`) to provide the plotfile,
    `-s slicefile_name` (or `--slicefile`) to provide the name of the file to
    output the results, and `--xctr x`, `--yctr y` to provide the coordinates
    of the center of the domain (x,y) (if not provided, both default to 0.0).
- `lgt_frnt1d`: `-p plotfile_name` (or `--plotfile`) to provide the plotfile,
    `-s slicefile_name` (or `--slicefile`) to provide the name of the file to
    output the results.
- `rad_shock`: `-p plotfile_name` (or `--plotfile`) to provide the plotfile,
    `-s slicefile_name` (or `--slicefile`) to provide the name of the file to
    output the results, and `--idir d` for the direction d along which to take
    the slice (where d is an integer 1-3 corresponding to the x-z directions,
    and defaults to 1).
- `rad_source`: the arguments are assumed to be a list of the plotfiles to be analyzed
- `rad_sphere`:  `-p plotfile_name` (or `--plotfile`) to provide the plotfile,
    `-g groupfile_name` (or `--groupfile`) to provide the name of the group file
    (this is output at runtime by the `RadSphere` problem), and `-r radius` (or
    `--radius`) to provide the radius at which to print out the
    radiation quantities.
- `rhd_shocktube`: `-p plotfile_name` (or `--plotfile`) to provide the plotfile,
    `-g groupfile_name` (or `--groupfile`) to provide the name of the group file
    (this is output at runtime by the `RadSphere` problem), and
    `-s slicefile_name` (or `--slicefile`) to provide the name of the file to
    output the results.
