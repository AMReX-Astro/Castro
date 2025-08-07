# Problem Description

This is a simple 1-d flame problem.  It sets up a fuel and ash region
with a smooth transition in temperature and composition.  The density
is computed such that the domain is at constant pressure initially.
Thermal diffusion heats the fuel to the point of ignition and a flame
ignites and propagates to the right.



# Usage

The different inputs files setup different fuels or solvers.

* `inputs.He` : a pure He flame at a density of 2.e6 g/cc, with 80 cm
  resolution.  It takes about 0.01 s for the flame to get established
  and has a speed of about 1.e5 cm/s.

* `inputs.H_He` : a mixed H/He flame appropriate for X-ray bursts.
  This has 20% H and 69% He (with the remainder O14/O15).  This uses
  1280 cm resolution, takes about 0.2 s to get established, and has a
  speed of about 1.6e5 cm/s.

  Note: The more H relative to He, the slower the flame.

  Note: the `NETWORK_DIR` will need to be set for this to include
  a network that has O14, O15, like `CNO_extras`.

* `inputs.C` : a pure C flame at a density of 5.e8 g/cc, appropriate
  for SN Ia.  This is a very small domain and the resolution
  requirements are very dependent on the density---this uses a resolution
  of about 2.e-5 cm.  It takes about 5.e-9 s to get established and has a speed
  of about 6e6 cm/s.

## Analysis

There are several scripts in the `analysis/` subdirectory.  These
include:

* `profiles.py` : plot the flame structure (T, enuc, and velocity)
  from a sequence of plotfiles, evenly spaced in time.  There are
  several options---do `./profiles.py -h` to see them.

  To change the interval between plotfiles plotted, use ``--dt``,
  like:

  ```
  ./profiles.py --dt 0.01
  ```

  This will loop over all plotfiles and plot the ones corresponding
  to 0s, 0.01s, 0.02s, ...

  Other options allow you to zoom in on the structure.

* `snapshot.py` : for a single plotfile, show the temperature,
  energy generation, and composition as a function of position.
  This will focus on only the most abundant nuclei.
  
# Publications

The SDC inputs files were used for a demonstration of a 4th order
accurate flame in

* *Improved Coupling of Hydrodynamics and Nuclear Reactions via
  Spectral Deferred Corrections*

  Zingale, M., Katz, M. P., Bell, J. B., Minion, M. L., Nonaka, A. J.,
  & Zhang, W. 2019, ApJ, 886, 2, p. 105 (DOI:
  10.3847/1538-4357/ab4e1d)
