#!/bin/bash


EXEC=./Castro1d.gnu.MPI.TRUESDC.ex


## sdc-2 + reflect

ofile=sdc2-reflect.converge.out

RUNPARAMS="
castro.time_integration_method=2
castro.sdc_order=2
castro.ppm_type=0
castro.use_pslope=1
castro.limit_fourth_order=1
castro.use_reconstructed_gamma1=1
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


## sdc-2 + ppm + reflect

ofile=sdc2-ppm-reflect.converge.out

RUNPARAMS="
castro.time_integration_method=2
castro.sdc_order=2
castro.ppm_type=1
castro.limit_fourth_order=1
castro.use_reconstructed_gamma1=1
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


## sdc-2 + ppm

ofile=sdc2-ppm.converge.out

RUNPARAMS="
castro.time_integration_method=2
castro.sdc_order=2
castro.ppm_type=1
castro.limit_fourth_order=1
castro.use_reconstructed_gamma1=1
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


## sdc-2 + HSE

ofile=sdc2.converge.out

RUNPARAMS="
castro.time_integration_method=2
castro.sdc_order=2
castro.ppm_type=0
castro.use_pslope=1
castro.limit_fourth_order=1
castro.use_reconstructed_gamma1=1
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



## sdc-4 + reflect

ofile=sdc4-reflect.converge.out

RUNPARAMS="
castro.time_integration_method=2
castro.sdc_order=4
castro.limit_fourth_order=1
castro.use_reconstructed_gamma1=1
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

