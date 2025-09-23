The Verification directory contains the analytic solutions for the
Sedov problem, as computed using Frank Timmes' exact Sedov routine,
http://www.cococubed.com/codes/sedov/sedov3.f

To compare Castro output, we need to radially bin the plotfile data to
get density, total velocity, and pressure as a function of radius.
This can be accomplished using main.cpp routine in
Castro/Diagnostics/Sedov

For the 2-d Cylindrical explosion in Cartesian coordinates:

  * run:
   ```
   ./Castro2d.gnu.MPI.ex inputs.2d.cyl_in_cartcoords
   ```

  * create the angle-averaged profile:

    compile Castro/Diagnostics/Sedov/fsedov2d_cyl_in_cartcoords.f90 as:
    ```
    make programs=fsedov2d_cyl_in_cartcoords
    ```
    and then run it as:
    ```
    ./fsedov2d_cyl_in_cartcoords.Linux.gfortran.exe -p ../plt00144 -s sedov2d.out
    ```

  * plot:

    A Gnuplot script to plot the data overtop the analytic solution is
    provided, simply execute
    ```
    gnuplot sedov2d.gp
    ```
    in the Sedov/Verification/ directory (where the above slices should have been
    output).


When running inputs.2d.sph_in_cylcoords or inputs.3d.sph, use sedov.gp instead
of sedov2d.gp

NOTES:

  * To get good agreement with the analytic solution for the Sedov
    problem, it is necessary to start off with a small timestep.  This
    is accomplished by setting castro.init_shrink = 0.1 in the inputs
    file.

  * Subsampling is critical to ensure that the perturbation is as
    spherical as possible.

  * It is important to use a lot of 'steps' in the sedov3.f routine to
    get a good analytic solution.
