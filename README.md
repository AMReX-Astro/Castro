# Castro
*an adaptive mesh, astrophysical radiation hydrodynamics simulation code*

`Castro` is an adaptive-mesh compressible radiation hydrodynamics
code for astrophysical flows.  `Castro` supports a general equation of
state, full Poisson gravity, and reactive flows, and is parallelized
with MPI + OpenMP.

More information on Castro can be found here:

http://amrex-astro.github.io/Castro/


## Getting Started

You can download Castro and all the necessary dependencies (AMReX and
the Microphysics repo) using the script:

https://raw.githubusercontent.com/AMReX-Astro/Castro/development/Util/scripts/get_castro.sh

Running this will clone the repositories and create a file
`castro_exports.sh` which you can source or copy to your `.bashrc`
to define the necessary environment variables.

The User's Guide in `Castro/Docs/` (type `make` to build from
LaTeX source) will guide you through running your first problem.

A PDF of the User's Guide can be found at

http://bender.astro.sunysb.edu/Castro/staging/Castro/Docs/CastroUserGuide.pdf


## Call Graph

A doxygen-generated call graph for `Castro` is available here:

http://bender.astro.sunysb.edu/Castro/staging/Castro/html/


## Development Model:

New features are committed to the `development` branch.  Nightly
regression testing is used to ensure that no answers change (or if
they do, that the changes were expected).  No changes should ever be
pushed directly into `master`.

On the first workday of each month, we perform a merge of
`development` into `master`, in coordination with `AMReX`, `Maestro`,
and `Microphysics`.  For this merge to take place, we need to be
passing the regression tests.  To accommodate this need, we close the
merge window into `development` a few days before the merge day.
While the merge window is closed, only bug fixes should be pushed into
`development`.  Once the merge from `development` -> `master` is done,
the merge window reopens.


## Mailing list

You can subscribe to the castro-help mailing list at google groups:

https://groups.google.com/forum/#!forum/castro-help
