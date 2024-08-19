#!/bin/bash

# echo the commands
set -x

DIM=3
EXEC=./Castro${DIM}d.gnu.MPI.ex

RUNPARAMS="
castro.ppm_type=1
castro.time_integration_method=0
problem.init_as_1d=1
geometry.prob_hi=1 0.125 0.125"

mpiexec -n 8 ${EXEC} inputs.64 amr.n_cell=64 8 8 ${RUNPARAMS} amr.plot_file=acoustic_pulse_64_ppm_plt &> 64.out
mpiexec -n 8 ${EXEC} inputs.128 amr.n_cell=128 16 16  ${RUNPARAMS} amr.plot_file=acoustic_pulse_128_ppm_plt &> 128.out
mpiexec -n 8 ${EXEC} inputs.256 amr.n_cell=256 32 32 ${RUNPARAMS} amr.plot_file=acoustic_pulse_256_ppm_plt &> 256.out

RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_64_ppm_plt00081 mediFile=acoustic_pulse_128_ppm_plt00161 fineFile=acoustic_pulse_256_ppm_plt00321 > convergence.${DIM}d.lo.ppm.x.out

#mpiexec -n 4 ${EXEC} inputs.512 ${RUNPARAMS} amr.plot_file=acoustic_pulse_512_ppm_plt &> 512.out

#RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_128_ppm_plt00161 mediFile=acoustic_pulse_256_ppm_plt00321 fineFile=acoustic_pulse_512_ppm_plt00641 > convergence.${DIM}d.hi.ppm.out

RUNPARAMS="
castro.ppm_type=1
castro.time_integration_method=0
problem.init_as_1d=2
geometry.prob_hi=0.125 1 0.125"

mpiexec -n 8 ${EXEC} inputs.64 amr.n_cell=8 64 8 ${RUNPARAMS} amr.plot_file=acoustic_pulse_64_ppm_plt &> 64.out
mpiexec -n 8 ${EXEC} inputs.128 amr.n_cell=16 128 16  ${RUNPARAMS} amr.plot_file=acoustic_pulse_128_ppm_plt &> 128.out
mpiexec -n 8 ${EXEC} inputs.256 amr.n_cell=32 256 32 ${RUNPARAMS} amr.plot_file=acoustic_pulse_256_ppm_plt &> 256.out

RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_64_ppm_plt00081 mediFile=acoustic_pulse_128_ppm_plt00161 fineFile=acoustic_pulse_256_ppm_plt00321 > convergence.${DIM}d.lo.ppm.y.out

#mpiexec -n 4 ${EXEC} inputs.512 ${RUNPARAMS} amr.plot_file=acoustic_pulse_512_ppm_plt &> 512.out

#RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_128_ppm_plt00161 mediFile=acoustic_pulse_256_ppm_plt00321 fineFile=acoustic_pulse_512_ppm_plt00641 > convergence.${DIM}d.hi.ppm.out


RUNPARAMS="
castro.ppm_type=1
castro.time_integration_method=0
problem.init_as_1d=3
geometry.prob_hi=0.125 0.125 1"

mpiexec -n 8 ${EXEC} inputs.64 amr.n_cell=8 8 64 ${RUNPARAMS} amr.plot_file=acoustic_pulse_64_ppm_plt &> 64.out
mpiexec -n 8 ${EXEC} inputs.128 amr.n_cell=16 16 128 ${RUNPARAMS} amr.plot_file=acoustic_pulse_128_ppm_plt &> 128.out
mpiexec -n 8 ${EXEC} inputs.256 amr.n_cell=32 32 256 ${RUNPARAMS} amr.plot_file=acoustic_pulse_256_ppm_plt &> 256.out

RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_64_ppm_plt00081 mediFile=acoustic_pulse_128_ppm_plt00161 fineFile=acoustic_pulse_256_ppm_plt00321 > convergence.${DIM}d.lo.ppm.z.out

#mpiexec -n 4 ${EXEC} inputs.512 ${RUNPARAMS} amr.plot_file=acoustic_pulse_512_ppm_plt &> 512.out

#RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_128_ppm_plt00161 mediFile=acoustic_pulse_256_ppm_plt00321 fineFile=acoustic_pulse_512_ppm_plt00641 > convergence.${DIM}d.hi.ppm.out


