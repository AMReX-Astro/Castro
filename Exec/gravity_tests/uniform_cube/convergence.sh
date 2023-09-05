#!/bin/bash

EXEC=./Castro3d.gnu.ex

${EXEC} inputs amr.n_cell=16 16 16 amr.plot_file=cube_plt_16_ &> cube_16.out
error_16=$(grep "Error" cube_16.out | awk '{print $3}')
echo "ncell = 16 error =" $error_16

${EXEC} inputs amr.n_cell=32 32 32 amr.plot_file=cube_plt_32_ &> cube_32.out
error_32=$(grep "Error" cube_32.out | awk '{print $3}')
echo "ncell = 32 error =" $error_32

${EXEC} inputs amr.n_cell=64 64 64 amr.plot_file=cube_plt_64_ &> cube_64.out
error_64=$(grep "Error" cube_64.out | awk '{print $3}')
echo "ncell = 64 error =" $error_64

echo "Average convergence rate =" $(echo "0.5 * (sqrt($error_16 / $error_32) + sqrt($error_32 / $error_64))" | bc -l)
