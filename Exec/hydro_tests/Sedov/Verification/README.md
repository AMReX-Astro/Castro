The Verification directory contains the analytic solutions for the
Sedov problem, as computed using Frank Timmes' exact Sedov routine,
http://www.cococubed.com/codes/sedov/sedov3.f

To compare Castro output, we need to radially bin the plotfile data to
get density, total velocity, and pressure as a function of radius.
This can be accomplished using main.cpp routine in
Castro/Diagnostics/Sedov

Cylindrical explosion is available in 1-d cylindrical and
2-d Cartesian coordinates using inputs.1d.cyl and inputs.2d.cyl_in_cartcoords,
respectively.

Spherical explosion is available in 1-d spherical, 2-d spherical,
2-d cylindrical and 3-d Cartesian coordinates by using inputs.1d.sph,
inputs.2d.sph_in_sphcoords, inputs.2d.sph_in_cylcoords, and inputs.3d.sph,
respectively.

For the 2-d Cylindrical explosion in Cartesian coordinates:
  * compile program:
    ```
    make DIM=2
    ```
  * run:
    ```
    ./Castro2d.gnu.MPI.ex inputs.2d.cyl_in_cartcoords
    ```

  * create the angle-averaged profile:

    compile main.cpp in Castro/Diagnostics/Sedov as:
    ```
    make DIM=2
    ```
    and then run it as:
    ```
    ./sedov_2d.ex -p ../sedov_2d_cyl_in_cart_plt00159 -s sedov_2d_cyl_in_cart.out
    ```

  * plot:

    A Gnuplot script to plot the data overtop the analytic solution is
    provided, simply execute
    ```
    gnuplot sedov_cyl.gp
    ```
    in the Sedov/Verification/ directory (where the above slices should have been
    output). The output is then saved in sedov_cyl.eps.
    When working with a spherical problem, using sedov_sph.gp instead.


NOTES:

  * To get good agreement with the analytic solution for the Sedov
    problem, it is necessary to start off with a small timestep.  This
    is accomplished by setting castro.init_shrink = 0.1 in the inputs
    file.

  * Subsampling is critical to ensure that the perturbation is as
    spherical as possible.

  * It is important to use a lot of 'steps' in the sedov3.f routine to
    get a good analytic solution.
