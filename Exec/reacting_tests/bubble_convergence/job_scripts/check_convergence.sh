#!/bin/sh

# SDC 4th order
CONV_TOOL=RichardsonConvergenceTest2d.gnu.sandybridge.ex

${CONV_TOOL} coarFile=bubble_64_sdc4_plt00667 mediFile=bubble_128_sdc4_plt01334 fineFile=bubble_256_sdc4_plt02667 > sdc_converge.lo.out
${CONV_TOOL} coarFile=bubble_128_sdc4_plt01334 mediFile=bubble_256_sdc4_plt02667 fineFile=bubble_512_sdc4_plt05334 > sdc_converge.hi.out


