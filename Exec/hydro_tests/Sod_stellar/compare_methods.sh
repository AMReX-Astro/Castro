#!/bin/bash

./Castro2d.gnu.ex inputs-test1-helm amr.plot_file=test1_sdc_plt castro.sdc_order=4 castro.time_integration_method=2 castro.fourth_order=1
./Castro2d.gnu.ex inputs-test1-helm amr.plot_file=test1_ppm_plt

./compare_sdc.py --ppm test1_ppm_plt????? --sdc test1_sdc_plt????? -o test1_sdc.png

./Castro2d.gnu.ex inputs-test2-helm amr.plot_file=test2_sdc_plt castro.sdc_order=4 castro.time_integration_method=2 castro.fourth_order=1
./Castro2d.gnu.ex inputs-test2-helm amr.plot_file=test2_ppm_plt

./compare_sdc.py --ppm test2_ppm_plt????? --sdc test2_sdc_plt????? -o test2_sdc.png

./Castro2d.gnu.ex inputs-test3-helm amr.plot_file=test3_sdc_plt castro.sdc_order=4 castro.time_integration_method=2 castro.fourth_order=1
./Castro2d.gnu.ex inputs-test3-helm amr.plot_file=test3_ppm_plt

./compare_sdc.py --ppm test3_ppm_plt????? --sdc test3_sdc_plt????? -o test3_sdc.png

./Castro2d.gnu.ex inputs-test4-helm amr.plot_file=test4_sdc_plt castro.sdc_order=4 castro.time_integration_method=2 castro.fourth_order=1
./Castro2d.gnu.ex inputs-test4-helm amr.plot_file=test4_ppm_plt

./compare_sdc.py --ppm test4_ppm_plt????? --sdc test4_sdc_plt????? -o test4_sdc.png

