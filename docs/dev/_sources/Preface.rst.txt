*******
Preface
*******

Welcome to the Castro User’s Guide!

In this User’s Guide we describe how to download and run Castro, a
massively parallel code that solves the multicomponent compressible
hydrodynamic equations for astrophysical flows including self-gravity,
nuclear reactions and radiation. Castro uses an Eulerian grid and
incorporates adaptive mesh refinement (AMR). Our approach to AMR uses
a nested hierarchy of logically-rectangular grids with simultaneous
refinement in both space and time, utilizing the
AMReX library.

The core algorithms in Castro are described in a series of papers:

  * *CASTRO: A New Compressible Astrophysical Solver. I. Hydrodynamics
    and Self-gravity*, A. S. Almgren, V. E. Beckner, J. B. Bell,
    M. S. Day, L. H. Howell, C. C. Joggerst, M. J. Lijewski,
    A. Nonaka, M. Singer, & M. Zingale, 2010, ApJ, 715, 1221
    http://dx.doi.org/10.1088/0004-637X/715/2/1221

  * *CASTRO: A New Compressible Astrophysical Solver. II. Gray
    Radiation
    Hydrodynamics*, W. Zhang, L. Howell, A. Almgren, A. Burrows,
    & J. Bell, 2011, ApJS, 196, 20
    http://dx.doi.org/10.1088/0067-0049/196/2/20

  * *CASTRO: A New Compressible Astrophysical Solver. III. Multigroup
    Radiation
    Hydrodynamics*, W. Zhang, L. Howell, A. Almgren, A. Burrows, J. Dolence,
    & J. Bell, 2013, ApJS, 204, 7
    http://dx.doi.org/10.1088/0067-0049/204/1/7

Improvements to the gravity solver and rotation were described in:

  * *Double White Dwarf Mergers on Adaptive Meshes I. Methodology and
    Code
    Verification*, M. P. Katz, M. Zingale, A. C. Calder, F. D. Swesty,
    A. S. Almgren, W. Zhang 2016, ApJ, 819, 94.
    http://dx.doi.org/10.3847/0004-637X/819/2/94

The development of AMReX library is led by the
Center for Computational Sciences and Engineering / Lawrence Berkeley
National Laboratory. Castro development is done collaboratively,
including the CCSE and Stony Brook University.

Castro *core developers* are those who have made substantial
contributions to the code. The process for becoming a core developer
is described in the `README.md <https://github.com/AMReX-Astro/Castro/blob/master/README.md>`_ in the Castro root directory.
Current Castro core developers are:

  * Ann Almgren
  * Maria G. Barrios Sazo
  * John Bell
  * Vince Beckner
  * Marc Day
  * Alice Harpole
  * Max Katz
  * Mike Lijewski
  * Chris Malone
  * Andy Nonaka
  * Don Willcox
  * Weiqun Zhang
  * Michael Zingale

All Castro development takes place on the project’s github
page, https://github.com/AMReX-Astro/Castro

External contributions are welcomed. Fork the Castro repo, modify your
local copy, and issue a pull-request to the AMReX-Astro/Castro
project. Further guidelines are given in the `README.md
<https://github.com/AMReX-Astro/Castro/blob/master/README.md>`_ file.

To get help, subscribe to the *castro-help* google group mailing list:
https://groups.google.com/forum/#!forum/castro-help

Acknowledging and Citing Castro
===============================

If you use Castro in your research, we would appreciate it if you
cited the relevant code papers describing its design, features, and
testing. A list of these can be found in the `CITATION
<https://github.com/AMReX-Astro/Castro/blob/master/CITATION>`_ file in
the root ``Castro/`` directory.

The development Castro is supported by the science application
interests of the contributors. There is a lot of effort behind the
scenes: testing, optimization, development of new features, bug
fixing, ..., that is often done under the radar. Nevertheless,
we are happy to volunteer our time to help new users come up to speed
with Castro. When significant new development / debugging for you
application is provided by a member of the Castro development
community, we would appreciate consideration of inviting the
developer(s) for co-authorship on any science paper that results.

