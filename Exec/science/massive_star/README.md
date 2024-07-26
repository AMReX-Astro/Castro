# massive star

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
