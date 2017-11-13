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

Development generally follows the following ideas:

  * New features are committed to the `development` branch.

    Nightly regression testing is used to ensure that no answers
    change (or if they do, that the changes were expected).

    If a change is critical, we can cherry-pick the commit from
    `development` to `master`.

  * Contributions are welcomed from anyone.  *Any contributions that
    have the potential to change answers should be done via pull
    requests.*   A pull request should be generated from your fork of
    Castro and target the `development` branch.  (If you mistakenly
    target `master`, we can change it for you.)

    Please add a line to `CHANGES` summarizing your change if it
    is a bug fix or new feature.  Reference the PR or issue as
    appropriate. Additionally, if your change fixes a bug (or if
    you find a bug but do not fix it), and there is no current
    issue describing the bug, please file a separate issue describing
    the bug, regardless of how significant the bug is. If possible,
    in both the `CHANGES` file and the issue, please cite the pull
    request numbers or git commit hashes where the problem was
    introduced and fixed, respectively.

    If there are a number of small commits making up the PR, we may
    wish to squash commits upon merge to have a clean history.
    *Please ensure that your PR title and first post are descriptive,
    since these will be used for a squashed commit message.*

  * On the first workday of each month, we perform a merge of
    `development` into `master`, in coordination with `AMReX`,
    `Maestro`, and `Microphysics`.  For this merge to take place, we
    need to be passing the regression tests.

    To accommodate this need, we close the merge window into
    `development` a few days before the merge day.  While the merge
    window is closed, only bug fixes should be pushed into
    `development`.  Once the merge from `development` -> `master` is
    done, the merge window reopens.


## Core Developers

People who make a number of substantive contributions will be named
"core developers" of Castro.  The criteria for becoming a core
developer are flexible, but generally involve one of the following:

  * 10 non-merge commits to `Castro/Source/` or `Castro/Docs/` or one
    of the problems that is not your own science problem *or*

  * addition of a new algorithm / module  *or*

  * substantial input into the code design process or testing

Core developers will be recognized in the following ways:

  * invited to the group's slack team

  * listed in the User's Guide and website as a core developer

  * invited to co-author general code papers / proceedings describing
    Castro, its performance, etc.  (Note: science papers will always
    be left to the science leads to determine authorship).

If a core developer is inactive for 3 years, we may reassess their
status as a core developer.



## Mailing list

You can subscribe to the castro-help mailing list at google groups:

https://groups.google.com/forum/#!forum/castro-help
