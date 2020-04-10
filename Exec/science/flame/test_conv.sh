#!/bin/bash

rm -rf plt*

# ./Castro1d.gnu.MPI.TRUESDC.ex inputs.1d.sdc.map.avg amr.n_cell=4096 max_step=1
# mv plt00000 plt00000.4096
# ./Castro1d.gnu.MPI.TRUESDC.ex inputs.1d.sdc.map.avg amr.n_cell=2048 max_step=1
# mv plt00000 plt00000.2048
# ./Castro1d.gnu.MPI.TRUESDC.ex inputs.1d.sdc.map.avg amr.n_cell=1024 max_step=1
# mv plt00000 plt00000.1024

./Castro1d.gnu.MPI.TRUESDC.ex inputs.1d.sdc.conv.512
./Castro1d.gnu.MPI.TRUESDC.ex inputs.1d.sdc.conv.1024
./Castro1d.gnu.MPI.TRUESDC.ex inputs.1d.sdc.conv.2048
#./Castro1d.gnu.MPI.TRUESDC.ex inputs.1d.sdc.conv.4096

RichardsonConvergenceTest1d.gnu.ex coarFile=plt_512_00050 mediFile=plt_1024_00100 fineFile=plt_2048_00200
#RichardsonConvergenceTest1d.gnu.ex coarFile=plt_1024_00100 mediFile=plt_2048_00200 fineFile=plt_4096_00400
#RichardsonConvergenceTest1d.gnu.ex coarFile=plt_1024_00000 mediFile=plt_2048_00000 fineFile=plt_4096_00000
