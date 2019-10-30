#!/bin/sh

# Strang

RichardsonConvergenceTest2d.gnu.sandybridge.ex coarFile=react_converge_64_strang_plt00301 mediFile=react_converge_128_strang_plt00601 fineFile=react_converge_256_strang_plt01201 > convergence.2d.lo.strang.out

RichardsonConvergenceTest2d.gnu.sandybridge.ex coarFile=react_converge_128_strang_plt00601 mediFile=react_converge_256_strang_plt01201 fineFile=react_converge_512_strang_plt02401 > convergence.2d.hi.strang.out


# SDC 2nd order

RichardsonConvergenceTest2d.gnu.sandybridge.ex coarFile=react_converge_64_sdc_plt00301 mediFile=react_converge_128_sdc_plt00601 fineFile=react_converge_256_sdc_plt01201 > convergence.2d.lo.sdc.out

RichardsonConvergenceTest2d.gnu.sandybridge.ex coarFile=react_converge_128_sdc_plt00601 mediFile=react_converge_256_sdc_plt01201 fineFile=react_converge_512_sdc_plt02401 > convergence.2d.hi.sdc.out


# SDC 4th order

RichardsonConvergenceTest2d.gnu.sandybridge.ex coarFile=react_converge_64_sdc4_plt00301 mediFile=react_converge_128_sdc4_plt00601 fineFile=react_converge_256_sdc4_plt01201 > convergence.2d.lo.sdc4.out

RichardsonConvergenceTest2d.gnu.sandybridge.ex coarFile=react_converge_128_sdc4_plt00601 mediFile=react_converge_256_sdc4_plt01201 fineFile=react_converge_512_sdc4_plt02401 > convergence.2d.hi.sdc4.out

