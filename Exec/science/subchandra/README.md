# subchandra

This is the double detonation Type Ia supernova setup.  It takes a 1-d
initial model representing a sub-Chandra WD and a surface He layer and
initializes a small perturbation at the pole in the He layer.

Initial models for this setup can be created via
https://github.com/AMReX-Astro/initial_models in the sub_chandra/ directory.

This setup has been used in the following papers:

* Zingale et al. 2022
  (https://ui.adsabs.harvard.edu/abs/2022ApJ...936....6Z/abstract)

  This paper introduced the simplified-SDC time-integration method and
  used the double detonation setup as a test problem.

* Zingale et al. 2024
  (https://ui.adsabs.harvard.edu/abs/2023arXiv230901802Z/abstract)

  This paper showed that the simplified-SDC method can model double
  detonations without cutting the timestep or limiting the reaction
  rates.  It also explored disabling burning in shocks for this
  problem.
