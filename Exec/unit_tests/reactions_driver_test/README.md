# Model file

The model file should have the same number of zone entries as there are in the grid.

So, e.g. for a 16 x 16 x 16 cell grid, the model file should have 16^3 entries.

You can use `make_react_model.py` to generate a new model file.

# Building on groot

make -j 4 COMP=PGI CUDA_ARCH=60 COMPILE_CUDA_PATH=/usr/local/cuda-9.2 USE_CUDA=TRUE USE_MPI=FALSE USE_OMP=FALSE USE_ACC=FALSE NETWORK_DIR=ignition_simple EOS_DIR=helmholtz INTEGRATOR_DIR=VODE DIMENSION_AGNOSTIC=TRUE
