# hse_convergence

This is meant to be a simple 1-d test for assessing the convergence of
hydro + gravity in maintaining HSE.  Convergence can be measure either
with the RichardsonConvergenceTest tool or by looking at the max |U|
in the plotfiles.

To run this problem, use one of the convergence scripts:

  * ``convergence_plm.sh`` :

    this runs CTU + PLM using the default HSE BCs and default
    use_pslope, then with reflect BCs, then without use_pslope, and
    finally runs with reflect instead of HSE BCs.

    These tests show that the best results come from HSE BCs + reflect vel

  * convergence_ppm.sh :

    this runs CTU + PPM in a similar set of configurations as PLM above

    These tests show that the best results come from HSE BCs + reflect vel

  * convergence_sdc.sh :

    this uses the TRUE_SDC integation, first with SDC-2 + PLM  and reflecting BCs,
    the SDC-2 + PPM and reflecting BCs, then the same but HSE BCs, and finally
    SDC-4 + reflect

    These tests show that the PLM + reflect (which uses the
    well-balanced use_pslope) and the SDC-4 + reflect give the lowest
    errors and expected (or better) convergence:


