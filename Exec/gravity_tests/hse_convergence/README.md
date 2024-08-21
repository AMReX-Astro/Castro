# hse_convergence

This is meant to be a simple 1-d test for assessing the convergence of
hydro + gravity in maintaining HSE.  Convergence can be measure either
with the RichardsonConvergenceTest tool or by looking at the max |U|
in the plotfiles.

To run this problem, use one of the convergence scripts:

  * `convergence_plm.sh` :

    this runs CTU + PLM using:
    1. the default HSE BCs and `use_pslope`
    2. the HSE BCs with reflection and `use_pslope`
    3. reflect BCs instead of HSE BCs without `use_pslope`
    4. reflect BCs with `use_pslope`

    These tests show that the best results (by far) come from
    `use_pslope=1` and reflecting BCs

  * convergence_ppm.sh :

    this runs CTU + PPM in a similar set of configurations as PLM above
    1. the default HSE BCs
    2. HSE BCs with reflection
    3. reflecting BCs
    4. reflecting BCs with `use_pslope`

    These tests show that the best results (by far) come from
    reflecting BCs with `use_pslope=1`, just like the PLM case.

  * convergence_sdc.sh :

    this uses the `TRUE_SDC` integration, with the following variations:
    1. SDC-2 + PLM and reflecting BCs
    2. SDC-2 + PPM and reflecting BCs
    3. SDC-2 + PLM with HSE BCs
    4. SDC-2 + PPM with HSE BCs
    5. SDC-4 + reflect

    These tests show that the PLM + reflect (which uses the
    well-balanced `use_pslope`) and the SDC-4 + reflect give the lowest
    errors and expected (or better) convergence.
