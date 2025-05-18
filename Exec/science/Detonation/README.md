# `Detonation`

A simple (carbon or helium) detonation.  The reaction network should
be set through the `GNUmakefile`.

> [!NOTE]
> By default this setup uses the simplified-SDC coupling of hydro
> and reactions, since this gives the best results for this
> problem.

This sets up a domain with a uniform density (`dens`).  A large
temperature (`T_l`) is placed in the left side of the domain, and an
ambient temperature (`T_r`) is in the right.  The parameter `center_T`
specifies where, as a fraction of the domain length, to put the
interface, while `width_T` determines how wide the transition region
is, as a fraction of the domain length. Finally, the parameters
`cfrac`, `nfrac`, and `ofrac`, are the initial C12, N14, O16 fraction,
and the remaining material is He4.

The left side and right side can optionally approach each other with a
non-zero infall velocity, `vel`.

When run, the large temperature in the left side instantly flashes the
fuel, generating a lot of energy and a large overpressure.  The
signature of a propagating detonation is that a narrow energy
generation zone will keep pace with the rightward propagating
shockwave (as seen in the pressure field).

Note: if the domain is too small, then the burning will decouple from
the shock wave, and you will not get a detonation.

## Inputs files

Some important inputs files

* `inputs-det-x.nse` : this produces a nice detonation that gets hot
  enough for the ash to be in NSE.

* `inputs-det-x.subch_base` : this produces a He4 detonation using the
  same initial condition in the shell burning stage of the subchandra
  problem.

* `inputs-collision` : this is the inputs file used for the 1D collision
  simulations from Katz & Zingale, 2019, ApJ, 874, 169

Inputs files used in the regression test suite:

* `inputs-det-x.test`

* `inputs-collision.testsuite`

* `inputs-det-x.simplified_sdc`

* `inputs-det-x.regrid`


## Publications

This problem setup has been used in the following papers:

* *Numerical Stability of Detonations in White Dwarf Simulations*

  Max P. Katz & M. Zingake, 2019, ApJ, 874, 169 (DOI:
  10.3847/1538-4357/ab0c00)

* *A Framework for Exploring Nuclear Physics Sensitivity in Numerical Simulations*

  Zhi Chen, Eric T. Johnson, Max Katz, Alexander Smith Clark, Brendan
  Boyd, & Michael Zingale, 2024, Journal of Physics: Conference Series,
  2742, 1, 012021 (DOI: 10.1088/1742-6596/2742/1/012021)

