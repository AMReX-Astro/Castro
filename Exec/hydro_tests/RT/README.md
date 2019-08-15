# Problem Description

This is a single mode Rayleigh-Taylor instability problem.  The
instability is triggered by perturbing the interface between the dense
and light fluid.

This problem was used in the Castro I paper to compare different PPM
types.


# Particles

This can be run with particles by building as:

```
make USE_PARTICLES=TRUE
```

and running with the `inputs_2d.particles` inputs file.  This will
assign particles according to the ``particle_file`` locations.

These particles are moved with the fluid velocity each time step, and
if they leave the domain then they just "disappear."

This test is primarily designed to test the ability of Castro to
handle particles.
