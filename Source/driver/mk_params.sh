#!/bin/sh

./parse_castro_params.py -m meth_params.template ./_cpp_parameters
PWD=`pwd`
cd ../../Docs/runtime_parameters/
./rp.py > runtime_parameters.tex
cd ${PWD}
