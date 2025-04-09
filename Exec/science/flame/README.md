# Problem Description

This is a simple 1-d flame problem.  It sets up a fuel and ash region
with a smooth transition in temperature and composition.  The density
is computed such that the domain is at constant pressure initially.
Thermal diffusion heats the fuel to the point of ignition and a flame
ignites and propagates to the right.



# Usage

The different inputs files setup different fuels or solvers.

* `inputs.He` : a pure He flame at a density of 2.e6 g/cc

* `inputs.H_He` : a mixed H/He flame appropriate for X-ray bursts.



A time series of flame profiles can be made using the
`analysis/profiles.py` script.  This has options to set the dt to skip
between plotfiles and to zoom in on the domain.

# Publications

The SDC inputs files were used for a demonstration of a 4th order
accurate flame in

* Improved Coupling of Hydrodynamics and Nuclear Reactions via
  Spectral Deferred Corrections Zingale, M., Katz, M. P., Bell, J. B.,
  Minion, M. L., Nonaka, A. J., & Zhang, W. 2019, ApJ, 886, 2, p. 105
