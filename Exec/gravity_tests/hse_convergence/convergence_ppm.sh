#!/bin/bash


EXEC=./Castro1d.gnu.MPI.ex

## ppm

ofile=ppm.converge.out

${EXEC} inputs.ppm.64 >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel > ${ofile}

${EXEC} inputs.ppm.128 >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

${EXEC} inputs.ppm.256 >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

${EXEC} inputs.ppm.512 >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}


## ppm + hse reflect

ofile=ppm-hsereflect.converge.out

RUNPARAMS="
castro.hse_interp_temp=1
castro.hse_reflect_vels=1
"""

${EXEC} inputs.ppm.64 ${RUNPARAMS} >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel > ${ofile}

${EXEC} inputs.ppm.128 ${RUNPARAMS} >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

${EXEC} inputs.ppm.256 ${RUNPARAMS} >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

${EXEC} inputs.ppm.512 ${RUNPARAMS} >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}


## ppm + grav_source_type = 4

ofile=ppm-grav4.converge.out

RUNPARAMS="
castro.grav_source_type=4
"""

${EXEC} inputs.ppm.64 ${RUNPARAMS} >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel > ${ofile}

${EXEC} inputs.ppm.128 ${RUNPARAMS} >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

${EXEC} inputs.ppm.256 ${RUNPARAMS} >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

${EXEC} inputs.ppm.512 ${RUNPARAMS} >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}


## ppm + reflect

ofile=ppm-reflect.converge.out

RUNPARAMS="
castro.lo_bc=3
castro.hi_bc=3
"""

${EXEC} inputs.ppm.64 ${RUNPARAMS} >& 64.out
pfile=`ls -t | grep -i hse_64_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel > ${ofile}

${EXEC} inputs.ppm.128 ${RUNPARAMS} >& 128.out
pfile=`ls -t | grep -i hse_128_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

${EXEC} inputs.ppm.256 ${RUNPARAMS} >& 256.out
pfile=`ls -t | grep -i hse_256_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

${EXEC} inputs.ppm.512 ${RUNPARAMS} >& 512.out
pfile=`ls -t | grep -i hse_512_plt | head -1`
fextrema.gnu.ex -v magvel ${pfile} | grep -i magvel >> ${ofile}

