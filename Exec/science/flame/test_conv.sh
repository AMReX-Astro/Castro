#!/bin/bash

rm -rf plt*

./Castro1d.gnu.MPI.ex inputs.1d.sdc.map.avg amr.n_cell=4096 max_step=1
mv plt00000 plt00000.4096
./Castro1d.gnu.MPI.ex inputs.1d.sdc.map.avg amr.n_cell=2048 max_step=1
mv plt00000 plt00000.2048
./Castro1d.gnu.MPI.ex inputs.1d.sdc.map.avg amr.n_cell=1024 max_step=1
mv plt00000 plt00000.1024
RichardsonConvergenceTest1d.gnu.ex coarFile=plt00000.1024 mediFile=plt00000.2048 fineFile=plt00000.4096
