#!/bin/bash

# echo the commands
set -x

DIM=2
EXEC=./Castro${DIM}d.gnu.MPI.ex

RUNPARAMS="castro.fourth_order=1 castro.mol_order=4 castro.time_integration_method=1"

mpiexec -n 4 ${EXEC} inputs.2d.64 ${RUNPARAMS} amr.plot_file=acoustic_pulse_64_rk4_plt &> 64.out
mpiexec -n 4 ${EXEC} inputs.2d.128 ${RUNPARAMS} amr.plot_file=acoustic_pulse_128_rk4_plt &> 128.out
mpiexec -n 4 ${EXEC} inputs.2d.256 ${RUNPARAMS} amr.plot_file=acoustic_pulse_256_rk4_plt &> 256.out
mpiexec -n 4 ${EXEC} inputs.2d.512 ${RUNPARAMS} amr.plot_file=acoustic_pulse_512_rk4_plt &> 512.out

RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_64_rk4_plt00081 mediFile=acoustic_pulse_128_rk4_plt00161 fineFile=acoustic_pulse_256_rk4_plt00321 > convergence.${DIM}d.lo.rk4.out
RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_128_rk4_plt00161 mediFile=acoustic_pulse_256_rk4_plt00321 fineFile=acoustic_pulse_512_rk4_plt00641 > convergence.${DIM}d.hi.rk4.out
