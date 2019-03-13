# Sedov diagnostic routines

The diagnostic routines in `Diagnostics/DustCollapse` are used to construct
the data which appears in Figure 12 in the first CASTRO paper.

The analytic maximum density has been computed from the initial density
and the analytic r(t) using conservation of mass; this is hard-wired
into the `main.cpp` file.

The radius of the star in each plotfile is computed by first computing the
radial average of density, then finding the radius at which
the density equals half of the analytic maximum density.

## Building and running

Typing 'make DIM=n' will build the diagnostic routine for the n-dimensional problem.

To run the 1d problem, run
```
./dustcollapse_1d.exe plotfile1 plotfile2
```
i.e. run the executable followed by a list of the plotfiles to be analyzed.
For the 2d and 3d problems, additional argument(s) can be passed in *before* the
list of plotfiles to give the coordinates of the domain center. If these are not
provided, then they default to 0.0. For the 2d problem, this would be
```
./dustcollapse_2d.exe --xctr x --yctr y plotfile1 plotfile2
```
and similarly for 3d
```
./dustcollapse_3d.exe --xctr x --yctr y --zctr z plotfile1 plotfile2
```
For the 2d and 3d problems, it is also possible to print the profile to file
by providing the argument `--profile`. This will create the profile file
`prof.profile` in the plotfile's directory.

## Analytic solution

To make the analytic profile, type
```
gfortran analytic.f90
```
then
```
a.out > analytic.txt
```
## Plotting

To use gnuplot to make Figure 12, use `DustCollapse.gp` once you have created the files
`analytic.txt`, `results_1d.txt`, `results_2d.txt` and `results_3d.txt`
