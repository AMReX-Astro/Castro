#!/bin/bash

mpiexec -n 4 ./Castro2d.gnu.MPI.ex inputs.128 >& /dev/null
mpiexec -n 4 ./Castro2d.gnu.MPI.ex inputs.256 >& /dev/null
mpiexec -n 4 ./Castro2d.gnu.MPI.ex inputs.512 >& /dev/null

RichardsonConvergenceTest2d.gnu.ex coarFile=nse_test_128_plt00160 mediFile=nse_test_256_plt00320 fineFile=nse_test_512_plt00640 >& nse_convergence.out

