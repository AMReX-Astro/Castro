# Problem Setup

This is a pure diffusion problem (no hydro).  It uses the explicit
diffusion solver in Castro to diffuse a Gaussian thermal pulse.
Because diffusion in Castro is incorporated in the energy equation, we
are solving:

   d(ρe)/dt = ∇·(k ∇T)

where L is the Laplacian and k is the thermal conductivity.  Here we
assume that k is constant, but Castro does not require that
assumption.  For this problem, we take ρ = constant, and using a
gamma-law EOS, we have e = c_v T, so this can take the form of a
diffusion equation as:

dT/dt = k/(ρ c_v) ∇^2 T = D ∇^2 T

where D is the constant diffusion coefficient.

Because we are doing the diffusion explicitly, there is a constraint
on the timestep,

dt < 0.5 dx^2/D

For constant diffusion coefficient, there is an analytic solution: the
diffusion of a Gaussian remains Gaussian, with the amplitude
decreasing and the width increasing with time.


# Testing

The problem creates a derived variable, `analytic`, which is the
analytic solution at the time of the plotfile output.  This allows you
to compare the current solution to the analytic solution to assess the
error.

It also will use the Castro problem-specific post-simulation hooks (in
`Prob.cpp`) to output the L-inf norm of the error (numerical vs
analyic solution) at the end of the simulation.  This can be used
for convergence testing.


## 2-d cylindrical with AMR

This uses the 2nd order accurate predictor-corrector formulation of
diffusion that is used with the CTU hydrodynamics solver.  A test of
the diffusion in 2-d cylindrical coordinates, with 2 levels of
refinement can be run as:

```
./Castro2d.gnu.ex inputs.2d.cyl amr.n_cell=64  128
./Castro2d.gnu.ex inputs.2d.cyl amr.n_cell=128 256
./Castro2d.gnu.ex inputs.2d.cyl amr.n_cell=256 512
```

At the end, each run will report the norm of the error against the
analytic solution, giving:

```
 base resolution      L-inf error
 (64 , 128)          0.0003707056645
 (128, 256)          9.414571162e-05
 (256, 512)          2.437072009e-05
```


## SDC-4 in 1-d

A convergence test of the 4th-order SDC algorithm can be run as:

```
./Castro1d.gnu.TRUESDC.ex inputs.1d castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=64
./Castro1d.gnu.TRUESDC.ex inputs.1d castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=128
./Castro1d.gnu.TRUESDC.ex inputs.1d castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=256
./Castro1d.gnu.TRUESDC.ex inputs.1d castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=512
```

Note: this is Cartesian, not spherical (since we don't have spherical
implemented to 4th order).  The norm of the error output at the end as
a function of resolution is:

```
 resolution     L-inf error
 64           8.639542243e-05
128           5.84073812e-06
256           3.725743019e-07
512           2.340631822e-08
```


# Non-constant Conductivity

There is no analytic solution for non-constant conductivity, so we can
only do convergence testing by varying the resolution.  In this
manner, the error reported at the end of a run is meaningless, since
it is comparing against the analytic solution for constant
conductivity.


## Second-order predictor-corrector algorithm

The powerlaw conductivity is simply k = k0 T^ν.  To build the test with this
conductivity we do:

```
 make DIM=1 CONDUCTIVITY_DIR=powerlaw -j 20
```

Tests can be run as:

```
 ./Castro1d.gnu.ex inputs.1d.powerlaw amr.n_cell=64
 mv diffuse_plt00048 diffuse_64
 ./Castro1d.gnu.ex inputs.1d.powerlaw
 mv diffuse_plt00191 diffuse_128
 ./Castro1d.gnu.ex inputs.1d.powerlaw amr.n_cell=256
 mv diffuse_plt00761 diffuse_256
```

Then the error can be measured using the RichardsonConvergenceTest
tool in `amrex/Tools/C_util/Convergence` as:

```
RichardsonConvergenceTest1d.gnu.ex coarFile=diffuse_64 mediFile=diffuse_128 fineFile=diffuse_256
```

This gives:

```
Level  L1 norm of Error in Each Component
-----------------------------------------------
Warning: BoxArray lengths are not the same at level 0
  0    Level  L1 norm of Error in Each Component
-----------------------------------------------
Warning: BoxArray lengths are not the same at level 0
  0    \begin{table}[p]
\begin{center}
\begin{tabular}{|cccc|} \hline
Variable & $e_{4h \rightarrow 2h}$ & Order & $e_{2h \rightarrow h}$\\
\hline
density&         0.000000e+00 & ------------ &0.000000e+00 \\
xmom&         0.000000e+00 & ------------ &0.000000e+00 \\
ymom&         0.000000e+00 & ------------ &0.000000e+00 \\
zmom&         0.000000e+00 & ------------ &0.000000e+00 \\
rho_E&         3.479414e-04 & 2.012958966 & 8.620750e-05 \\
rho_e&         3.479414e-04 & 2.012958966 & 8.620750e-05 \\
Temp&         3.479414e-04 & 2.012958966 & 8.620750e-05 \\
rho_X&         0.000000e+00 & ------------ &0.000000e+00 \\
```

(some bits were edited out)

e.g. we see second-order convergence in the temperature


## 4th order SDC

This is built as the previous test:

```
make DIM=1 CONDUCTIVITY_DIR=powerlaw -j 20
```

and then run as:

```
./Castro1d.gnu.ex inputs.1d.powerlaw castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=64
mv diffuse_plt00048 diffuse_64
./Castro1d.gnu.ex inputs.1d.powerlaw castro.time_integration_method=2 castro.sdc_order=4
mv diffuse_plt00190 diffuse_128
./Castro1d.gnu.ex inputs.1d.powerlaw castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=256
mv diffuse_plt00760 diffuse_256
./Castro1d.gnu.ex inputs.1d.powerlaw castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=512
mv diffuse_plt03040 diffuse_512

RichardsonConvergenceTest1d.gnu.ex coarFile=diffuse_64 mediFile=diffuse_128 fineFile=diffuse_256 > convergence_diffusion.1d.lo.sdc4.out
RichardsonConvergenceTest1d.gnu.ex coarFile=diffuse_128 mediFile=diffuse_256 fineFile=diffuse_512 > convergence_diffusion.1d.hi.sdc4.out
```

This gives (for the lower resolution runs):

```
Level  L1 norm of Error in Each Component
-----------------------------------------------
Warning: BoxArray lengths are not the same at level 0
  0    Level  L1 norm of Error in Each Component
-----------------------------------------------
Warning: BoxArray lengths are not the same at level 0
  0    \begin{table}[p]
\begin{center}
\begin{tabular}{|cccc|} \hline
Variable & $e_{4h \rightarrow 2h}$ & Order & $e_{2h \rightarrow h}$\\
\hline
density&         0.000000e+00 & ------------ &0.000000e+00 \\
xmom&         0.000000e+00 & ------------ &0.000000e+00 \\
ymom&         0.000000e+00 & ------------ &0.000000e+00 \\
zmom&         0.000000e+00 & ------------ &0.000000e+00 \\
rho_E&         1.111626e-05 & 3.948910124 & 7.198104e-07 \\
rho_e&         1.111626e-05 & 3.948910124 & 7.198104e-07 \\
Temp&         1.063477e-05 & 3.952987539 & 6.866892e-07 \\
rho_X&         0.000000e+00 & ------------ &0.000000e+00 \\
pressure&         7.410837e-06 & 3.948910124 & 4.798736e-07 \\
```

e.g. we see fourth-order convergence in the temperature


## 4th order SDC (2-d)

This is built as the previous test:

```
make DIM=2 CONDUCTIVITY_DIR=powerlaw -j 20 USE_MPI=TRUE
```

and then run as:

```
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.2d.powerlaw castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=64 64
mv diffuse_plt00039 diffuse_2d_64
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.2d.powerlaw castro.time_integration_method=2 castro.sdc_order=4
mv diffuse_plt00157 diffuse_2d_128
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.2d.powerlaw castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=256 256
mv diffuse_plt00626 diffuse_2d_256
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.2d.powerlaw castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=512 512
mv diffuse_plt02504 diffuse_2d_512

RichardsonConvergenceTest2d.gnu.ex coarFile=diffuse_2d_64 mediFile=diffuse_2d_128 fineFile=diffuse_2d_256 > convergence_diffusion.2d.lo.sdc4.out
RichardsonConvergenceTest2d.gnu.ex coarFile=diffuse_2d_128 mediFile=diffuse_2d_256 fineFile=diffuse_2d_512 > convergence_diffusion.2d.hi.sdc4.out
```

This gives (for the lower resolution runs):

```
Level  L1 norm of Error in Each Component
-----------------------------------------------
Warning: BoxArray lengths are not the same at level 0
  0    Level  L1 norm of Error in Each Component
-----------------------------------------------
  0    \begin{table}[p]
\begin{center}
\begin{tabular}{|cccc|} \hline
Variable & $e_{4h \rightarrow 2h}$ & Order & $e_{2h \rightarrow h}$\\
\hline
density&         0.000000e+00 & ------------ &0.000000e+00 \\
xmom&         0.000000e+00 & ------------ &0.000000e+00 \\
ymom&         0.000000e+00 & ------------ &0.000000e+00 \\
zmom&         0.000000e+00 & ------------ &0.000000e+00 \\
rho_E&         1.902161e-06 & 3.957610923 & 1.224299e-07 \\
rho_e&         1.902161e-06 & 3.957610923 & 1.224299e-07 \\
Temp&         1.770452e-06 & 3.966033724 & 1.132894e-07 \\
```

e.g. we see fourth-order convergence in the temperature
