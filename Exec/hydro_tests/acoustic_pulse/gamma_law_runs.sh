#!/bin/sh


make realclean 2>1 > /dev/null
make EOS_DIR=gamma_law -j 4 2>1 > /dev/null
status=$?

if [ $status -ne 0 ]; then
    echo "make failed"
    exit
fi

echo "running the CTU version..."
./convergence.sh

echo "running 2nd order SDC..."
./convergence_sdc.sh

echo "running 4th order SDC..."
./convergence_sdc4.sh
