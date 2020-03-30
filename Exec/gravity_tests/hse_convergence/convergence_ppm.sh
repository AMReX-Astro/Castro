#!/bin/bash


## ppm

ofile=ppm.converge.out

./Castro1d.gnu.MPI.ex inputs.ppm.64 >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel > ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.128 >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.256 >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.512 >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}


## ppm + grav_source_type = 4

ofile=ppm-grav4.converge.out

RUNPARAMS="
castro.grav_source_type=4
"""

./Castro1d.gnu.MPI.ex inputs.ppm.64 ${RUNPARAMS} >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel > ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.128 ${RUNPARAMS} >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.256 ${RUNPARAMS} >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.512 ${RUNPARAMS} >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}


## ppm + reflect

ofile=ppm-reflect.converge.out

RUNPARAMS="
castro.lo_bc=3
castro.hi_bc=3
"""

./Castro1d.gnu.MPI.ex inputs.ppm.64 ${RUNPARAMS} >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel > ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.128 ${RUNPARAMS} >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.256 ${RUNPARAMS} >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.512 ${RUNPARAMS} >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}



## plm

ofile=plm.converge.out

RUNPARAMS="
castro.ppm_type=0
"""

./Castro1d.gnu.MPI.ex inputs.ppm.64 ${RUNPARAMS} >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel > ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.128 ${RUNPARAMS} >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.256 ${RUNPARAMS} >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.512 ${RUNPARAMS} >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}


## plm + well-balanced

ofile=plm-wellbalanced.converge.out

RUNPARAMS="
castro.ppm_type=0
castro.plm_well_balanced=1
castro.use_pslope=0
"""

./Castro1d.gnu.MPI.ex inputs.ppm.64 ${RUNPARAMS} >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel > ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.128 ${RUNPARAMS} >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.256 ${RUNPARAMS} >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.512 ${RUNPARAMS} >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}



## plm + reflect

ofile=plm-reflect.converge.out

RUNPARAMS="
castro.ppm_type=0
castro.lo_bc=3
castro.hi_bc=3
"""

./Castro1d.gnu.MPI.ex inputs.ppm.64 ${RUNPARAMS} >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel > ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.128 ${RUNPARAMS} >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.256 ${RUNPARAMS} >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.512 ${RUNPARAMS} >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}



## plm + reflect + nopslope

ofile=plm-reflect-nopslope.converge.out

RUNPARAMS="
castro.ppm_type=0
castro.lo_bc=3
castro.hi_bc=3
castro.use_pslope=0
"""

./Castro1d.gnu.MPI.ex inputs.ppm.64 ${RUNPARAMS} >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel > ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.128 ${RUNPARAMS} >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.256 ${RUNPARAMS} >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.512 ${RUNPARAMS} >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}



## plm + reflect + well-balanced

ofile=plm-reflect-wellbalanced.converge.out

RUNPARAMS="
castro.ppm_type=0
castro.lo_bc=3
castro.hi_bc=3
castro.plm_well_balanced=1
castro.use_pslope=0
"""

./Castro1d.gnu.MPI.ex inputs.ppm.64 ${RUNPARAMS} >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel > ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.128 ${RUNPARAMS} >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.256 ${RUNPARAMS} >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

./Castro1d.gnu.MPI.ex inputs.ppm.512 ${RUNPARAMS} >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1` 
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}




