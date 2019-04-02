#!/bin/bash

# echo the commands
set -x

DIM=1
EXEC=./Castro${DIM}d.gnu.MPI.ex

RUNPARAMS="castro.sdc_order=4 castro.time_integration_method=2 castro.limit_fourth_order=1 castro.use_reconstructed_gamma1=1 castro.sdc_solve_for_rhoe=1 castro.sdc_solver_tol=1.e-8 castro.sdc_solver=2 castro.sdc_use_analytic_jac=0"

mpiexec -n 8 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=512


