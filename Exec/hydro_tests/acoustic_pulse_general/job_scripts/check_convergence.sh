#!/bin/sh

# SDC 4th order
CONV_TOOL=RichardsonConvergenceTest2d.gnu.sandybridge.ex

${CONV_TOOL} coarFile=acoustic_pulse_64_sdc4_plt00101 mediFile=acoustic_pulse_128_sdc4_plt00201 fineFile=acoustic_pulse_256_sdc4_plt00401 > convergence.2d.lo.sdc4.out
${CONV_TOOL} coarFile=acoustic_pulse_128_sdc4_plt00201 mediFile=acoustic_pulse_256_sdc4_plt00401 fineFile=acoustic_pulse_512_sdc4_plt00801 > convergence.2d.hi.sdc4.out


