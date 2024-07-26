# flame_wave

This is the XRB flame setup.  It has been used for several papers,
including:

* Eiden et al. 2020
  (https://ui.adsabs.harvard.edu/abs/2020ApJ...894....6E)

  This modeled pure He flames and used "boosted" flames, which were
  artificially sped up.

* Harpole et al. 2021
  (https://ui.adsabs.harvard.edu/abs/2021ApJ...912...36H/abstract)

  This also modeled pure He flames, but without any boosting.  A
  variety of rotation rates were explored.  The inputs files for this
  are in `inputs_He/`

* Zingale et al. 2023
  (https://ui.adsabs.harvard.edu/abs/2023ApJ...952..160Z/abstract)

  This explored both 2D and 3D pure He flames.  Another change was in
  the lower boundary condition -- a reflecting boundary was used an
  `castro.use_pslope` was enabled.  The inputs files for this are in
  `inputs_He/` and also available on Zenodo:
  https://zenodo.org/record/7692201

* Chen et al. 2023
  (https://ui.adsabs.harvard.edu/abs/2023ApJ...955..128C/abstract)

  This looked at 2D pure He flames and explored different reaction
  networks.  The inputs files are available on Zenodo:
  https://zenodo.org/record/8117761

