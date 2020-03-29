#!/bin/bash

ofile=ppm.converge.out

./Castro1d.gnu.MPI.ex inputs.ppm.64 >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} > ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.128 >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.256 >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.512 >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} >> ${ofile}



