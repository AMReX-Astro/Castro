# `Gradient Detonation`

This is a similar setup to `Detonation`, except that we want to study
the impact of composition gradient on shock detection by applying
a linear interpolation function to the composition.

A simple (carbon or helium) detonation.  The reaction network should
be set through the GNUmakefile.

This sets up a domain with a uniform density (dens).  A large
temperature (T_ign) is placed in the accretion region with temperature T_l on the left side of the domain, and an ambient temperature (T_r) representing the WD core is in the right.  The transition from envelope composition to core composition occurs at the center of the domain, with the width of said region specified by delta. Finally, the parameter cfrac_core, nfrac_core, and ofrac_core, are the initial C12, N14, O16 fraction in the core, and the remaining material is He4. This is the same for cfrac_env...etc, which represents the material in the envelope.

The left side and right side can optionally approach each other with a
non-zero infall velocity, vel.

When run, the large temperature in the left side instantly flashes the
fuel, generating a lot of energy and a large overpressure.  The
signature of a propagating detonation is that a narrow energy
generation zone will keep pace with the rightward propagating
shockwave (as seen in the pressure field).

Note: if the domain is too small, then the burning will decouple from
the shock wave, and you will not get a detonation.
