#!/bin/bash

# echo the commands
set -x

DIM=2
EXEC=./Castro${DIM}d.gnu.MPI.ex

RUNPARAMS="castro.mol_order=4 castro.time_integration_method=1"

mpiexec -n 4 ${EXEC} inputs.64 ${RUNPARAMS} amr.plot_file=acoustic_pulse_64_rk4_plt &> 64.out
mpiexec -n 4 ${EXEC} inputs.128 ${RUNPARAMS} amr.plot_file=acoustic_pulse_128_rk4_plt &> 128.out
mpiexec -n 4 ${EXEC} inputs.256 ${RUNPARAMS} amr.plot_file=acoustic_pulse_256_rk4_plt &> 256.out
mpiexec -n 4 ${EXEC} inputs.512 ${RUNPARAMS} amr.plot_file=acoustic_pulse_512_rk4_plt &> 512.out

RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_64_rk4_plt00101 mediFile=acoustic_pulse_128_rk4_plt00201 fineFile=acoustic_pulse_256_rk4_plt00401 > convergence.${DIM}d.lo.rk4.out
RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_128_rk4_plt00201 mediFile=acoustic_pulse_256_rk4_plt00401 fineFile=acoustic_pulse_512_rk4_plt00801 > convergence.${DIM}d.hi.rk4.out


