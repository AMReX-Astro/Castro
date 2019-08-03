#!/bin/bash

# echo the commands
set -x

DIM=1
EXEC=./Castro${DIM}d.gnu.MPI.ex

RUNPARAMS="
castro.time_integration_method=2
castro.sdc_order=4
castro.limit_fourth_order=1
castro.use_reconstructed_gamma1=1
castro.sdc_solve_for_rhoe=1
castro.sdc_solver_tol_dens=1.e-8
castro.sdc_solver_tol_spec=1.e-8
castro.sdc_solver_tol_ener=1.e-6
castro.sdc_solver=2
castro.sdc_use_analytic_jac=0"

#mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=1024 amr.plot_file=flame_1024_sdc4_plt &> sdc4_1024.out
mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=2048 amr.plot_file=flame_2048_sdc4_plt &> sdc4_2048.out
mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=4096 amr.plot_file=flame_4096_sdc4_plt &> sdc4_4096.out
mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=8192 amr.plot_file=flame_8192_sdc4_plt &> sdc4_8192.out


