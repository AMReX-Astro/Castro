White dwarf merger collision science test
=========================================

This directory includes three setups of white dwarf mergers. Each setup is
the same except for the spatial resolution. Two equal mass white dwarfs
(0.64 solar masses each) collide with each other nearly head-on. The collision
is expected to result in a thermonuclear ignition and detonation.

To build the problem on OLCF's Summit supercomputer, have the following modules loaded:

```
spectrum-mpi/10.3.1.2-20200121
cuda/10.1.243
pgi/19.10
```

Then you can do `make -j16` to compile, producing the executable
`Castro3d.pgi.TPROF.MPI.CUDA.ex`. The build process will also symlink the file
`helm_table.dat` into your directory. To run the science problem, copy the executable
and the `helm_table.dat` file into each of the three subdirectories, `collision_ml0`,
`collision_ml1`, and `collision_ml2` (which are the low, medium, and high resolution
cases, respectively). Each directory has a `run_script.sh` which can be submitted on
Summit using `bsub run_script.sh`. The runs will not all complete in a single job
submission, so when you need to resume a job, add `amr.restart = chk?????` to the end
of the execution line in `run_script.sh` where `chk?????` is the last checkpoint file
that was generated (numerically highest).

Results can be visualized, for example from the `smallplt?????` files, using the `yt`
package, which can be installed with `conda install -c conda-forge yt`. For example,
you can use the `vol_render_density` function in the wdmerger analysis script,
https://github.com/AMReX-Astro/wdmerger/blob/master/analysis/wdmerger.py. An example
invocation would be

```
python -c 'import wdmerger; wdmerger.vol_render_density("collision_ml0.png", "smallplt00595", zoom_factor=0.55, annotate_max_T=True)'
```

which would generate a 3D volume render of the results in plotfile `smallplt00595`, and
would include a dot at the location of the peak temperature.
