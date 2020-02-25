#!/bin/bash

./Castro2d.gnu.ex inputs-sod-x amr.plot_file=sod_sdc_plt castro.sdc_order=4 castro.time_integration_method=2 amr.n_cell=128 8 8 geometry.prob_hi=1 0.0625 0.0625 amr.max_level=0
./Castro2d.gnu.ex inputs-sod-x amr.plot_file=sod_ppm_plt amr.n_cell=128 8 8 geometry.prob_hi=1 0.0625 0.0625 amr.max_level=0

./compare_sdc.py --ppm sod_ppm_plt????? --sdc sod_sdc_plt????? -o sod_sdc.png


./Castro2d.gnu.ex inputs-test2-x amr.plot_file=test2_sdc_plt castro.sdc_order=4 castro.time_integration_method=2 amr.n_cell=128 8 8 geometry.prob_hi=1 0.0625 0.0625 amr.max_level=0
./Castro2d.gnu.ex inputs-test2-x amr.plot_file=test2_ppm_plt  amr.n_cell=128 8 8 geometry.prob_hi=1 0.0625 0.0625 amr.max_level=0

./compare_sdc.py --ppm test2_ppm_plt????? --sdc test2_sdc_plt????? -o test2_sdc.png


./Castro2d.gnu.ex inputs-test3-x amr.plot_file=test3_sdc_plt castro.sdc_order=4 castro.time_integration_method=2 amr.n_cell=128 8 8 geometry.prob_hi=1 0.0625 0.0625 amr.max_level=0
./Castro2d.gnu.ex inputs-test3-x amr.plot_file=test3_ppm_plt amr.n_cell=128 8 8 geometry.prob_hi=1 0.0625 0.0625 amr.max_level=0

./compare_sdc.py --ppm test3_ppm_plt????? --sdc test3_sdc_plt????? -o test3_sdc.png




