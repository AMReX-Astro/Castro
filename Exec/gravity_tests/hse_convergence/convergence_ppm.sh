#!/bin/sh

nohup ./Castro1d.gnu.MPI.ex inputs.ppm.64 >& 64.out &
nohup ./Castro1d.gnu.MPI.ex inputs.ppm.128 >& 128.out &
nohup ./Castro1d.gnu.MPI.ex inputs.ppm.256 >& 256.out &
nohup ./Castro1d.gnu.MPI.ex inputs.ppm.512 >& 512.out &
