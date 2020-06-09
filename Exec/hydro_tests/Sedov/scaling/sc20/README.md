Sedov scaling test
==================

This directory contains a Sedov benchmark problem for testing weak scaling.
The problem is the Sedov-Taylor blast wave, a very commonly solved hydrodynamics
problem which is often used for hydrodynamics benchmarking. This setup was
designed for use on the OLCF Summit supercomputer.

To build and run, have the following modules loaded:

```
spectrum-mpi/10.3.1.2-20200121
cuda/10.1.243
pgi/19.10
```

`make -j16` will then build the executable `Castro3d.pgi.TPROF.MPI.CUDA.ex`.
Then run `bash run_scaling.sh` to submit a series of training scripts in the
`scaling_results` subdirectory. Each script will correspond to inputs files with
a different problem size. `run_script.sh` is the template, which has text replacement
done on various parameters before it is submitted for each of the cases. Note that
you will need to ensure the project flag (`-P`) is valid for you if running on Summit.

To generate a scaling plot based on the results in the `scaling_results` directory,
simply run the `plot.py` script, which depends only on matplotlib and numpy. It is
recommended to have the `python/3.6.6-anaconda3-5.3.0` module loaded.
