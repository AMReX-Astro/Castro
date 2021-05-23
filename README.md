[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2301848.svg)](https://doi.org/10.5281/zenodo.2301848)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02513/status.svg)](https://doi.org/10.21105/joss.02513)
[![AMReX](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io)
[![yt-project](https://img.shields.io/static/v1?label="works%20with"&message="yt"&color="blueviolet")](https://yt-project.org)
[![github pages](https://github.com/AMReX-Astro/Castro/workflows/github%20pages/badge.svg)](https://github.com/AMReX-Astro/Castro/actions?query=workflow%3A%22github+pages%22)

![Castro](https://github.com/AMReX-Astro/Castro/blob/development/Util/logo/castro_logo_hot_200.png)

*an adaptive mesh, astrophysical radiation hydrodynamics simulation code*

`Castro` is an adaptive-mesh compressible radiation / MHD / hydrodynamics
code for astrophysical flows.  `Castro` supports a general equation of
state, full Poisson gravity, and reactive flows, and is parallelized
with MPI + OpenMP for CPUs and MPI + CUDA for GPUs.

More information on Castro can be found here:

http://amrex-astro.github.io/Castro/


## Getting Started

The "Getting Started" section of the User's Guide walks you
through running your first problem:

https://amrex-astro.github.io/Castro/docs/getting_started.html

This will have you clone Castro and its dependencies (AMReX and
StarKiller Microphysics),

The User's Guide in written in re-structured text using Sphinx, with
the source in `Castro/Docs/`, and is built automatically
from the `development` branch.

## Running at Supercomputer Centers

Documentation for running the AMReX Astrophysics codes at popular
supercomputing centers can be found at:
https://amrex-astro.github.io/workflow/

## Call Graph

A doxygen-generated call graph for `Castro` is available here:

http://bender.astro.sunysb.edu/Castro/staging/Castro/html/


## Development Model:

Development generally follows the following ideas:

  * New features are committed to the `development` branch.

    Nightly regression testing is used to ensure that no answers
    change (or if they do, that the changes were expected).

    If a change is critical, we can cherry-pick the commit from
    `development` to `main`.

  * Contributions are welcomed from anyone in the form of a pull
    request from your fork of Castro, targeting the `development`
    branch. (If you mistakenly target `main`, we can change it
    for you.)

    Please add a line to `CHANGES.md` summarizing your change if it
    is a bug fix or new feature.  Reference the PR or issue as
    appropriate. Additionally, if your change fixes a bug (or if
    you find a bug but do not fix it), and there is no current
    issue describing the bug, please file a separate issue describing
    the bug, regardless of how significant the bug is. If possible,
    in both the `CHANGES.md` file and the issue, please cite the pull
    request numbers or git commit hashes where the problem was
    introduced and fixed, respectively.

    We will squash commits upon merge to have a clean history.
    *Please ensure that your PR title and and the PR summary field are
    descriptive, since these will be used for a squashed commit message.*

  * On the first workday of each month, we perform a merge of
    `development` into `main`, in coordination with `AMReX`,
    `Maestro`, and `Microphysics`.  For this merge to take place, we
    need to be passing the regression tests.

    To accommodate this need, we close the merge window into
    `development` a few days before the merge day.  While the merge
    window is closed, only bug fixes should be pushed into
    `development`.  Once the merge from `development` -> `main` is
    done, the merge window reopens.


## Core Developers

People who make a number of substantive contributions will be named
"core developers" of Castro.  The criteria for becoming a core
developer are flexible, but generally involve one of the following:

  * 10 non-merge commits to `Castro/Source/` or `Castro/Docs/`
    or one of the problems that is not your own science problem *or*

  * addition of a new algorithm / module  *or*

  * substantial input into the code design process or testing

Core developers will be recognized in the following ways:

  * invited to the group's slack team

  * listed in the User's Guide and website as a core developer

  * listed in the author list on the Zenodo DOI for the project
    (as given in the .zenodo.json file)

  * invited to co-author general code papers / proceedings describing
    Castro, its performance, etc.  (Note: science papers will always
    be left to the science leads to determine authorship).

If a core developer is inactive for 3 years, we may reassess their
status as a core developer.



## Mailing list

You can subscribe to the castro-help mailing list at google groups:

https://groups.google.com/forum/#!forum/castro-help
