#!/bin/bash

# echo the commands
set -x

DIM=1
EXEC=./Castro${DIM}d.gnu.MPI.ex

RUNPARAMS="
castro.sdc_order=2
castro.time_integration_method=2
castro.sdc_solve_for_rhoe=1
castro.sdc_solver_tol_dens=1.e-8
castro.sdc_solver_tol_spec=1.e-6
castro.sdc_solver_tol_ener=1.e-5
castro.sdc_solver_atol=1.e-5
castro.sdc_solver=2
castro.sdc_use_analytic_jac=0
castro.sdc_solver_relax_factor=5"

mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=2048 castro.dtnuc_e=0.1 amrex.fpe_trap_invalid=1 stop_time=1.5e-4
# castro.init_shrink=0.01 castro.cfl=0.1
# max_step=25 





