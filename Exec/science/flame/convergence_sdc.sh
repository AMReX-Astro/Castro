#!/bin/bash

# echo the commands
set -x

DIM=1
EXEC=./Castro${DIM}d.gnu.MPI.ex

RUNPARAMS="castro.sdc_order=2 castro.time_integration_method=2 castro.ppm_type=0 castro.sdc_solve_for_rhoe=1 castro.sdc_solver_tol=1.e-8 castro.sdc_solver=2 castro.sdc_use_analytic_jac=0"

mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=1024 amr.plot_file=flame_1024_sdc_plt &> sdc_1024.out
mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=2048 amr.plot_file=flame_2048_sdc_plt &> sdc_2048.out
mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=4096 amr.plot_file=flame_4096_sdc_plt &> sdc_4096.out



