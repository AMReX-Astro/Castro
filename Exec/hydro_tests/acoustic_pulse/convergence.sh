#!/bin/bash

mpiexec -n 8 ./Castro2d.gnu.MPI.ex inputs.2d.64 castro.fourth_order=1 castro.mol_order=4 castro.time_integration_method=1
mpiexec -n 8 ./Castro2d.gnu.MPI.ex inputs.2d.128 castro.fourth_order=1 castro.mol_order=4 castro.time_integration_method=1
mpiexec -n 8 ./Castro2d.gnu.MPI.ex inputs.2d.256 castro.fourth_order=1 castro.mol_order=4 castro.time_integration_method=1
mpiexec -n 8 ./Castro2d.gnu.MPI.ex inputs.2d.512 castro.fourth_order=1 castro.mol_order=4 castro.time_integration_method=1

RichardsonConvergenceTest2d.gnu.ex coarFile=acoustic_pulse_128_plt00161 mediFile=acoustic_pulse_256_plt00321 fineFile=acoustic_pulse_512_plt00641 mediError=med.out coarError=coar.out

