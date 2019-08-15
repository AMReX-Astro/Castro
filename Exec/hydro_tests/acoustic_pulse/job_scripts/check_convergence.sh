#!/bin/sh

# SDC 4th order
CONV_TOOL=RichardsonConvergenceTest2d.gnu.sandybridge.ex

${CONV_TOOL} coarFile=acoustic_pulse_64_sdc4_plt00081 mediFile=acoustic_pulse_128_sdc4_plt00161 fineFile=acoustic_pulse_256_sdc4_plt00321 > convergence.2d.lo.sdc4.out
${CONV_TOOL} coarFile=acoustic_pulse_128_sdc4_plt00161 mediFile=acoustic_pulse_256_sdc4_plt00321 fineFile=acoustic_pulse_512_sdc4_plt00641 > convergence.2d.hi.sdc4.out

