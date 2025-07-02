# `massive_star`

This setup models the convective burning shells in a massive star
leading up to core collapse.  Some details:

* The initial model is from MESA and is put onto a uniform grid using
  the routines in the AMReX Astrophysics initial model repo
  (https://github.com/amrex-astro/initial_models) in the
  ``massive_star`` directory.

* We use simplified-SDC together with aprox19 and an NSE table
  (generated from pynucastro)

* The Castro ``drive_initial_convection`` functionality is used to
  establish the initial convective velocity field.

## Publications

This problem setup was used in:

* *Strong Coupling of Hydrodynamics and Reactions in Nuclear
  Statistical Equilibrium for Modeling Convection in Massive Stars*

  Zingale, Michael, Chen, Zhi, Johnson, Eric T., Katz, Max P., &
  Clark, Alexander Smith 2024, The Astrophysical Journal, 977, 1,
  p. 30
