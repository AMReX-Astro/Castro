#!/bin/bash

# echo the commands
set -x

DIM=1
EXEC=./Castro${DIM}d.gnu.MPI.TRUESDC.ex

RUNPARAMS="
castro.sdc_order=2
castro.time_integration_method=2
castro.ppm_type=0
castro.sdc_solver=1
castro.use_retry=0"

mpiexec -n 8 ${EXEC} inputs.64 ${RUNPARAMS} &> 64.out
mpiexec -n 16 ${EXEC} inputs.128 ${RUNPARAMS} &> 128.out
mpiexec -n 16 ${EXEC} inputs.256 ${RUNPARAMS} &> 256.out

RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=react_converge_64_plt00301 mediFile=react_converge_128_plt00601 fineFile=react_converge_256_plt01201 > convergence.${DIM}d.lo.sdc2.out

mpiexec -n 16 ${EXEC} inputs.512 ${RUNPARAMS} &> 512.out

RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=react_converge_128_plt00601 mediFile=react_converge_256_plt01201 fineFile=react_converge_512_plt02401 > convergence.${DIM}d.hi.sdc2.out
