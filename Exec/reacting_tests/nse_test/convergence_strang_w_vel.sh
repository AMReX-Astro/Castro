#!/bin/bash

# echo the commands
set -x

DIM=2
EXEC=./Castro${DIM}d.gnu.MPI.ex

RUNPARAMS="
problem.u0=1.e8
problem.v0=1.e8
"
mpiexec -n 4 ${EXEC} inputs.64 ${RUNPARAMS} >& /dev/null
mpiexec -n 4 ${EXEC} inputs.128 ${RUNPARAMS} >& /dev/null
mpiexec -n 4 ${EXEC} inputs.256 ${RUNPARAMS} >& /dev/null

RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=nse_test_64_plt00125 mediFile=nse_test_128_plt00250 fineFile=nse_test_256_plt00500  >& nse_convergence_strang_lo.out

mpiexec -n 4 ${EXEC} inputs.512 ${RUNPARAMS} >& /dev/null

RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=nse_test_128_plt00250 mediFile=nse_test_256_plt00500 fineFile=nse_test_512_plt01000 >& nse_convergence_strang_hi.out


