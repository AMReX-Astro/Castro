# Castro
*an adaptive mesh, astrophysical radiation hydrodynamics simulation code*

`Castro` is an adaptive-mesh compressible radiation hydrodynamics
code for astrophysical flows.  `Castro` supports a general equation of
state, full Poisson gravity, and reactive flows, and is parallelized
with MPI + OpenMP.

More information on Castro can be found here:

http://boxlib-codes.github.io/Castro/


## Getting Started

To build `Castro`, you need a copy of the `BoxLib` library:

https://github.com/BoxLib-Codes/BoxLib.git

(optionally) for Radiation, you'll want to get the
[CastroRadiation](https://github.com/BoxLib-Codes/CastroRadiation) source too.

There is a User's Guide in `Castro/Docs/` (type `make` to build
from LaTeX source) that will guide you through running your first
problem.  A PDF of the User's Guide can be found at

http://bender.astro.sunysb.edu/Castro/staging/Castro/Docs/CastroUserGuide.pdf


## Call Graph

A doxygen-generated call graph for `Castro` is available here:

http://bender.astro.sunysb.edu/Castro/staging/Castro/html/


## Development Model:

New features are committed to the `development` branch.  Nightly
regression testing is used to ensure that no answers change (or if
they do, that the changes were expected).  No changes should ever
be pushed directly into `master`.

On the first workday of each month, we perform a merge of
`development` into `master`, in coordination with `BoxLib`, `Maestro`,
and `Microphysics`.  For this merge to take place, we need to be
passing the regression tests.  To accommodate this need, we close the
merge window into `development` a few days before the merge day.
While the merge window is closed, only bug fixes should be pushed into
`development`.  Once the merge from `development` -> `master` is done,
the merge window reopens.


## Mailing list

You can subscribe to the castro-help mailing list at google groups:

https://groups.google.com/forum/#!forum/castro-help
