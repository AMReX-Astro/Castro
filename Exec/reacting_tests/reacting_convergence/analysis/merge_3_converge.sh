#!/bin/bash

# merge 2 sets of outputs from create_pretty_table for when we want to
# have 5 resolutions all compared together

python create_pretty_tables.py --simple convergence.2d.lo.strang.ppm.out convergence.2d.hi.strang.ppm.out > _lo.out
python create_pretty_tables.py --simple convergence.2d.hi.strang.ppm.out convergence.2d.vhi.strang.ppm.out | awk '{print $5 "      " $6}' > _hi.out

paste _lo.out _hi.out

