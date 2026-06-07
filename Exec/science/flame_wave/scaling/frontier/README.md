New scaling study done on June 6, 2026.

For the Microphysics repo, I am using this PR:
https://github.com/AMReX-Astro/Microphysics/pull/1998

Tests were run with several integrators (VODE, Rosenbrock, and RKC) and two
different inputs files (one for He and one for H/He).

* He (`inputs.He.25cm.static.1000Hz`)

  * we try two networks: `iso7` and `ase`

* H/He (`inputs.H_He.64cm_x_16cm.static.1000Hz`)

  * we do one network: `cno-he-burn-34am`

  * Rosenbrock does not fit on < 256 nodes (out of memory)

  * RKC does not work at all (lots of retries)

  * VODE with a 32-bit Jacobian seems to work really well
