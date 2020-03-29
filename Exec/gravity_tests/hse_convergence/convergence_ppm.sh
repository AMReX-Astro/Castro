#!/bin/bash

./Castro1d.gnu.MPI.ex inputs.ppm.64 >& 64.out
./Castro1d.gnu.MPI.ex inputs.ppm.128 >& 128.out
./Castro1d.gnu.MPI.ex inputs.ppm.256 >& 256.out
./Castro1d.gnu.MPI.ex inputs.ppm.512 >& 512.out

fextrema.gnu.ex -v magvel hse_64_plt00667 > ppm.converge.out
fextrema.gnu.ex -v magvel hse_128_plt01334 >> ppm.converge.out
fextrema.gnu.ex -v magvel hse_256_plt02667 >> ppm.converge.out
fextrema.gnu.ex -v magvel hse_512_plt05334 >> ppm.converge.out

