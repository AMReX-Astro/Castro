#!/bin/bash

# echo the commands
set -x

EXEC=./Castro2d.gnu.MPI.TRUESDC.ex

CONV_TOOL=RichardsonConvergenceTest2d.gnu.ex

mpiexec -n 8 ${EXEC} inputs_2d.32 amr.plot_file=bubble_32_sdc4_plt >& 32.out
mpiexec -n 8 ${EXEC} inputs_2d.64 amr.plot_file=bubble_64_sdc4_plt >& 64.out
mpiexec -n 8 ${EXEC} inputs_2d.128 amr.plot_file=bubble_128_sdc4_plt >& 128.out

${CONV_TOOL} coarFile=bubble_32_sdc4_plt00334 mediFile=bubble_64_sdc4_plt00667 fineFile=bubble_128_sdc4_plt01334 > sdc_converge.lo.out

mpiexec -n 8 ${EXEC} inputs_2d.256 amr.plot_file=bubble_256_sdc4_plt >& 256.out

${CONV_TOOL} coarFile=bubble_64_sdc4_plt00667 mediFile=bubble_128_sdc4_plt01334 fineFile=bubble_256_sdc4_plt02667 > sdc_converge.mid.out

mpiexec -n 8 ${EXEC} inputs_2d.512 amr.plot_file=bubble_512_sdc4_plt >& 512.out

${CONV_TOOL} coarFile=bubble_128_sdc4_plt01334 mediFile=bubble_256_sdc4_plt02667 fineFile=bubble_512_sdc4_plt05334 > sdc_converge.hi.out
