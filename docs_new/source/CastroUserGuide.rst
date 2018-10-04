.. raw:: latex

   \frontmatter

|  
| CASTRO

--------------

| 

--------------

|  
| User’s Guide

.. raw:: latex

   \vfill

.. raw:: latex

   \shorttoc{Chapter Listing}{0}

.. raw:: latex

   \setcounter{tocdepth}{2}

.. raw:: latex

   \tableofcontents

.. raw:: latex

   \clearpage

.. raw:: latex

   \listoffigures

.. raw:: latex

   \addcontentsline{toc}{chapter}{list of figures}

.. raw:: latex

   \clearpage

.. raw:: latex

   \listoftables

.. raw:: latex

   \addcontentsline{toc}{chapter}{list of tables}

.. raw:: latex

   \clearpage

Preface
=======

.. raw:: latex

   \markboth{\chaptername
   \ \thechapter.\ Preface}{}

.. raw:: latex

   \addcontentsline{toc}{chapter}{preface}

Welcome to the Castro User’s Guide!

In this User’s Guide we describe how to download and run Castro, a
massively parallel code that solves the multicomponent compressible
hydrodynamic equations for astrophysical flows including self-gravity,
nuclear reactions and radiation. Castro uses an Eulerian grid and
incorporates adaptive mesh refinement (AMR). Our approach to AMR uses
a nested hierarchy of logically-rectangular grids with simultaneous
refinement in both space and time, utilizing the
AMReX library [1]_.

The core algorithms in Castro are described in a series of papers:

-  *CASTRO: A New Compressible Astrophysical Solver. I. Hydrodynamics and Self-gravity*,
   A. S. Almgren, V. E. Beckner, J. B. Bell, M. S. Day, L. H. Howell, C. C. Joggerst, M. J. Lijewski,
   A. Nonaka, M. Singer, & M. Zingale, 2010, ApJ, 715, 1221
   http://dx.doi.org/10.1088/0004-637X/715/2/1221

-  *CASTRO: A New Compressible Astrophysical Solver. II. Gray Radiation Hydrodynamics*,
   W. Zhang, L. Howell, A. Almgren, A. Burrows, & J. Bell, 2011, ApJS, 196, 20
   http://dx.doi.org/10.1088/0067-0049/196/2/20

-  *CASTRO: A New Compressible Astrophysical Solver. III. Multigroup Radiation Hydrodynamics*,
   W. Zhang, L. Howell, A. Almgren, A. Burrows, J. Dolence, & J. Bell, 2013, ApJS, 204, 7
   http://dx.doi.org/10.1088/0067-0049/204/1/7

Improvements to the gravity solver and rotation were described in:

-  *Double White Dwarf Mergers on Adaptive Meshes I. Methodology
   and Code Verification,*
   M. P. Katz, M. Zingale, A. C. Calder, F. D. Swesty, A. S. Almgren, W. Zhang
   2016, ApJ, 819, 94.
   http://dx.doi.org/10.3847/0004-637X/819/2/94

The development of AMReX library is led by the
Center for Computational Sciences and Engineering / Lawrence Berkeley
National Laboratory. Castro development is done collaboratively,
including the CCSE and Stony Brook University.

Castro *core developers* are those who have made substantial
contributions to the code. The process for becoming a core developer
is described in the README.md in the Castro root directory.
Current Castro core developers are:

    Ann Almgren
    Maria G. Barrios Sazo
    John Bell
    Vince Beckner
    Marc Day
    Max Katz
    Mike Lijewski
    Chris Malone
    Andy Nonaka
    Don Willcox
    Weiqun Zhang
    Michael Zingale

| All Castro development takes place on the project’s github
  page
| https://github.com/AMReX-Astro/Castro
| External contributions are welcomed. Fork the Castro repo, modify
  your local copy, and issue a pull-request to the
  AMReX-Astro/Castro project. Further guidelines are given in the
  README.md file.

To get help, subscribe to the *castro-help* google group mailing list:
https://groups.google.com/forum/#!forum/castro-help

Acknowledging and Citing Castro
-------------------------------

If you use Castro in your research, we would appreciate it if you
cited the relevant code papers describing its design, features, and
testing. A list of these can be found in the
`CITATION <https://github.com/AMReX-Astro/Castro/blob/master/CITATION>`__ file in the root Castro/ directory.

The development Castro is supported by the science application
interests of the contributors. There is a lot of effort behind the
scenes: testing, optimization, development of new features, bug
fixing, :math:`\ldots`, that is often done under the radar. Nevertheless,
we are happy to volunteer our time to help new users come up to speed
with Castro. When significant new development / debugging for you
application is provided by a member of the Castro development
community, we would appreciate consideration of inviting the
developer(s) for co-authorship on any science paper that results.

.. raw:: latex

   \clearpage

.. raw:: latex

   \mainmatter

Introduction
============

Introduction to Castro
----------------------

Castro is a adaptive mesh, radiation hydrodynamics code that is
designed to model astrophysical reacting flows on massively parallel
computers.

The major capabilities:

-  1-, 2-, and 3-dimensional unsplit, 2nd-order hydrodynamics

-  multigroup flux-limited diffusion radiation hydrodynamics

-  adaptive mesh refinement with subcycling; jumps of 2x and 4x between levels

-  arbitrary equation of state (gamma-law and stellar EOSes are bundled)

-  general nuclear reaction networks

-  explicit thermal diffusion

-  full Poisson gravity (with isolated boundary conditions)

-  rotation (in the co-rotating frame) in 2-d axisymmetric and 3-d.

-  parallelization via MPI + OpenMP

Units and Conventions
---------------------

Castro works in CGS units unless otherwise specified.
Table \ `[table:units] <#table:units>`__ shows some of the common symbols / names used
throughout the code documentation and papers.

.. raw:: latex

   \centering

.. table:: [table:units] Common quantities and units.

   +-----------------------+-----------------------+-----------------------+
   | name                  | units                 | description           |
   +=======================+=======================+=======================+
   | :math:`t`             | s                     | time                  |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\rho`          | :math:`\mathrm{g~cm^{ | mass density          |
   |                       | -3}}`                 |                       |
   +-----------------------+-----------------------+-----------------------+
   | :math:`{\bf u}`       | :math:`\mathrm{cm~s^{ | velocity vector       |
   |                       | -1}}`                 |                       |
   +-----------------------+-----------------------+-----------------------+
   | :math:`p`             | :math:`\mathrm{dyn~cm | pressure              |
   |                       | ^{-2}}`               |                       |
   +-----------------------+-----------------------+-----------------------+
   | :math:`{\bf g}`       | :math:`\mathrm{cm~s^{ | gravitational         |
   |                       | -2}}`                 | acceleration          |
   +-----------------------+-----------------------+-----------------------+
   | :math:`{\bf S}`       | varies                | source term           |
   +-----------------------+-----------------------+-----------------------+
   | :math:`E`             | :math:`\mathrm{erg~g^ | specific total energy |
   |                       | {-1}}`                |                       |
   +-----------------------+-----------------------+-----------------------+
   | :math:`e`             | :math:`\mathrm{erg~g^ | specific internal     |
   |                       | {-1}}`                | energy                |
   +-----------------------+-----------------------+-----------------------+
   | :math:`T`             | :math:`K`             | temperature           |
   +-----------------------+-----------------------+-----------------------+
   | :math:`{k_\mathrm{th} | :math:`\mathrm{erg~cm | thermal conductivity  |
   | }`                    | ^{-1}~s^{-1}~K~{-1}}` |                       |
   +-----------------------+-----------------------+-----------------------+
   | :math:`X_k`           | –                     | mass fraction of      |
   |                       |                       | species :math:`k`     |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\dot\omega_k`  | :math:`\mathrm{s^{-1} | species creation rate |
   |                       | }`                    | (from reactions)      |
   +-----------------------+-----------------------+-----------------------+

Physical constants, again using the CGS system are available
in Castro/constants/constants_cgs.f90

Getting Started
===============

Downloading the Code
--------------------

Castro is built on top of the AMReX framework. In order to run
Castro  you must download two separate git modules.

.. raw:: latex

   \vspace{.1in}

First, make sure that git is installed on your machine—we recommend version 1.7.x or higher.

.. raw:: latex

   \vspace{.1in}

#. Clone/fork the AMReX repository from the AMReX-Codes
   github page (https://github.com/AMReX-Codes/amrex/). To
   clone via the command line, simply type:

   ::

       git clone https://github.com/AMReX-Codes/amrex.git

   Alternately, if you have a github account with your
   machine’s SSH-keys registered, you can do:

   ::

       git clone ssh://git@github.com/AMReX-Codes/amrex.git

   This will create a directory called amrex/ on your machine.

   You will want to periodically update AMReX by typing

   ::

       git pull

   in the amrex/ directory.

   Note: actively development is done on the development branch
   in each repo, and merged into the master branch periodically.
   If you wish to use the Castro development branch, then you
   should also switch to the development branch for AMReX.

#. Set the environment variable, AMREX_HOME, on your
   machine to point to the path name where you have put AMReX.
   You can add this to your .bashrc as:

   ::

       export AMREX_HOME={\em /path/to/amrex/}

   where you replace ``/path/to/amrex/`` will the full path to the
   amrex/ directory.

#. Clone/fork the Castro repository from the same
   github organization as above, using either HTTP access:

   ::

       git clone https://github.com/AMReX-Astro/Castro.git

   or SSH access if you have it enabled:

   ::

       git clone ssh://git@github.com:/AMReX-Astro/Castro.git

   Or, as above, you can download a ZIP file of the code from
   `our main github page <https://github.com/AMReX-Astro>`__,
   by clicking on the Castro link.

   As with AMReX, development on Castro is done in the
   development branch, so you should work there if you want
   the latest source.

#. We recommend setting the CASTRO_HOME environment
   variable to point to the path name where you have put Castro.
   Add the following to your .bashrc:

   ::

       export CASTRO_HOME="/path/to/Castro/"

#. (optional) An additional repository, Microphysics.git is
   available at the starkiller-astro github page. This add
   additional reaction networks and EOSes and can be cloned following
   the same procedure as above [2]_:

   ::

       git clone https://github.com/starkiller-astro/Microphysics.git

   or via SSH as

   ::

       git clone ssh://git@github.com:/starkiller-astro/Microphysics.git

   To access the Microphysics routines, set the MICROPHYSICS_HOME
   environment variable to point to the Microphysics/ directory.

Building the Code
-----------------

In Castro each different problem setup is stored in its own
sub-directory under Castro/Exec/. You build the
Castro executable in the problem sub-directory. Here we’ll
build the Sedov problem:

#. From the directory in which you checked out the Castro git repo,
   type

   ::

       cd Castro/Exec/hydro_tests/Sedov

   This will put you into a directory in which you can run the Sedov
   problem in 1-d, 2-d or 3-d.

#. In Sedov/, edit the GNUmakefile, and set

   -  This is the dimensionality—here we pick 2-d.

   -  This is the set of compilers. gnu are a good default
      choice (this will use g++ and gfortran. You can
      also choose pgi and intel for example.

      If you want to try other compilers than the GNU suite and they
      don’t work, please let us know.

   -  This disabled debugging checks and results in a more
      optimized executable.

   -  This turns off parallelization via MPI. Set it to TRUE to
      build with MPI—this requires that you have the MPI library
      installed on your machine. In this case, the build system will
      need to know about your MPI installation. This can be done by
      editing the makefiles in the AMReX tree, but the default
      fallback is to look for the standard MPI wrappers (e.g. 
      mpic++ and mpif90) to do the build.

#. Now type make.

   The resulting executable will look something like
   Castro2d.Linux.gnu.ex, which means this is a 2-d version
   of the code, made on a Linux machine, with COMP = gnu.

Running the Code
----------------

#. Castro takes an input file that overrides the runtime parameter defaults.
   The code is run as:

   ::

       Castro2d.Linux.gcc.gfortran.ex inputs.2d.cyl_in_cartcoords

   This will run the 2-d cylindrical Sedov problem in Cartesian (:math:`x`-:math:`y`
   coordinates). You can see other possible options, which should be
   clear by the names of the inputs files.

#. You will notice that running the code generates directories that
   look like plt00000/, plt00020/, etc, and chk00000/,
   chk00020/, etc. These are “plotfiles” and “checkpoint”
   files. The plotfiles are used for visualization, the checkpoint
   files are used for restarting the code.

Visualization of the Results
----------------------------

There are several options for visualizing the data. The popular
VisIt package supports the AMReX file format natively, as does the
yt python package [3]_. The standard tool used within the
AMReX-community is Amrvis, which we demonstrate here. Amrvis is available on github.

#. Get Amrvis:

   ::

       git clone https://github.com/AMReX-Codes/Amrvis

   Then cd into Amrvis/, edit the GNUmakefile there
   to set DIM = 2, and again set COMP to compilers that
   you have. Leave DEBUG = FALSE.

   Type make to build, resulting in an executable that
   looks like amrvis2d...ex.

   If you want to build amrvis with DIM = 3, you must first
   download and build volpack:

   ::

       git clone https://ccse.lbl.gov/pub/Downloads/volpack.git

   Then cd into volpack/ and type make.

   Note: Amrvis requires the OSF/Motif libraries and headers. If you don’t have these
   you will need to install the development version of motif through your package manager.
   On most Linux distributions, the motif library is provided by the
   openmotif package, and its header files (like Xm.h) are provided
   by openmotif-devel. If those packages are not installed, then use the
   package management tool to install them, which varies from
   distribution to distribution, but is straightforward.
   lesstif gives some functionality and will allow you to build the amrvis executable,
   but Amrvis may not run properly.

   You may then want to create an alias to amrvis2d, for example

   ::

       alias amrvis2d /tmp/Amrvis/amrvis2d...ex

   where /tmp/Amrvis/amrvis2d...ex is the full path and name of the Amrvis executable.

#. Configure Amrvis:

   Copy the amrvis.defaults file to your home directory (you can
   rename it to .amrvis.defaults if you wish). Then edit the
   file, and change the palette line to point to the full
   path/filename of the Palette file that comes with Amrvis.

#. Visualize:

   Return to the Castro/Exec/hydro_tests/Sedov directory. You should
   have a number of output files, including some in the form pltXXXXX,
   where XXXXX is a number corresponding to the timestep the file
   was output.
   amrvis2d *filename* to see a single plotfile, or amrvis2d -a
   \*plt\*, which will animate the sequence of plotfiles.

   Try playing
   around with this—you can change which variable you are
   looking at, select a region and click “Dataset” (under View)
   in order to look at the actual numbers, etc. You can also export the
   pictures in several different formats under "File/Export".

   Some users have found that Amrvis does not work properly under X
   with the proprietary Nvidia graphics driver. A fix for this is
   provided in the FAQ (§ `5 <#ch:faq:vis>`__)—this is due to the default
   behavior of the DAC in mappuing colors.

   Note: yt is a great alternative to using Amrvis for visualization,
   and understands Castro plotfiles well.

   Please know that we do have a number of conversion routines to other
   formats (such as matlab), but it is hard to describe them all. If you
   would like to display the data in another format, please let us know
   (again, asalmgren@lbl.gov) and we will point you to whatever we have
   that can help.

You have now completed a brief introduction to Castro.

Other Distributed Problem Setups
--------------------------------

There are a number of standard problem setups that come with Castro.
These can be used as a starting point toward writing your own setup.
We organize these into subdirectories by broad type (radiation, hydro,
gravity, etc.): The standard categories and *some* of the included
problems are:

-  gravity_tests:

   -  DustCollapse:

      A pressureless cloud collapse that is a standard test problem for
      gravity. An analytic solution that describes the radius of the
      sphere as a function of time is found in Colgate and
      White :raw-latex:`\cite{colgwhite}`. This problem is also found in the FLASH
      User’s Guide.

   -  hydrostatic_adjust:

      Model a 1-d stellar atmosphere (plane-parallel or
      spherical/self-gravitating) and dump energy in via an analytic
      heat source and watch the atmosphere’s hydrostatic state adjust in
      response. This is the counterpart to the Maestro 
      test_basestate unit test.

-  hydro_tests:

   -  double_bubble:

      Initialize 1 or 2 bubbles in a stratified atmosphere (isothermal
      or isentropic) and allow for the bubbles to have the same or a
      different :math:`\gamma` from one another / the background atmosphere.
      This uses the multigamma EOS.

      An analogous problem is implemented in Maestro.

   -  HCBubble:

   -  KH:

      A Kelvin-Helmholtz shear instability problem.

   -  oddeven:

      A grid-aligned shock hitting a very small density perturbation.
      This demonstrates the odd-even decoupling problem discussed in
      :raw-latex:`\cite{quirk1997}`. This setup serves to test the
      castro.hybrid_riemann option to hydrodynamics.

   -  reacting_bubble:

      A reacting bubble in a stratified white dwarf atmosphere. This
      problem was featured in the Maestro reaction
      paper :raw-latex:`\cite{maestro:III}`.

   -  RT:

      A single-model Rayleigh-Taylor instability problem.

   -  RT_particles:

   -  Sedov:

      The standard Sedov-Taylor blast wave problem. This setup was used
      in the first Castro paper :raw-latex:`\cite{castro_I}`.

   -  Sod:

      A one-dimensional shock tube setup, including the classic Sod
      problem. This setup was used in the original Castro paper.

   -  Sod_stellar:

      A version of the Sod shock tube for the general stellar equation
      of state. This setup and the included inputs files was used
      in :raw-latex:`\cite{zingalekatz}`.

   -  toy_convect:

      A simple nova-like convection problem with an external heating
      source. This problem shows how to use the model parser to
      initialize a 1-d atmosphere on the Castro grid, incorporate a
      custom tagging routine, sponge the fluid above the atmosphere, and
      write a custom diagnostics routine.

      A Maestro version of this problem setup also exists.

-  radiation_tests:

-  science:

-  unit_tests:

Inputs Files
============

The Castro executable uses two inputs files at runtime to set and alter the
behavior of the algorithm and initial conditions.

The main inputs file, typically named inputs
is used to set BoxLib parameters and the control flow in the
C portions of the Castro code. Each parameter here has a
namespace (like amr.\ *optionname* or castro.
*optionname*). Parameters set here are read using the
BoxLib ParmParse class infrastructure.

The second inputs file, typically named probin is used by the Fortran code that initializes the problem
setup. It is read at problem initialization (via a Fortran
namelist) and the problem-specific quantities are stored in a
Fortran module defined in the problem’s
probdata.f90 file.

Only the inputs file is specified on the commandline. The
associated probin file is specified in the inputs file
using the amr.probin_file parameter, e.g.,

::

    amr.probin_file = my_special_probin

for example, has the Fortran code read a file called my_special_probin.

Working with probin Files
-------------------------

There are three different Fortran namelists that can be defined in the
probin file:

-  &fortin is the main namelist read by the problem’s probinit
   subroutine in the Prob_?d.f90 file.

-  &extern is used to set different microphysics options

-  &tagging is used to get the parameters (defined in )
   that affect how we tag for refinement.

Common inputs Options
---------------------

**Important**: because the inputs file is handled by the C portion of
the code, any quantities you specify in scientific notation, must take the
form 1.e5 and not 1.d5—the ‘d’ specifier is not recognized.

Additionally, note that in Castro, all quantities are in CGS units.

Problem Geometry
~~~~~~~~~~~~~~~~

The geometry namespace is used by BoxLib to define the
computational domain. The main parameters here are:

-  : physical location of low corner of the
   domain (type: Real; must be set)

   Note: a number is needed for each dimension in the problem

-  : physical location of high corner of the
   domain (type: Real; must be set)

   Note: a number is needed for each dimension in the problem

-  : coordinate system, 0 = Cartesian,
   1 = :math:`r`-:math:`z` (2-d only), 2 = spherical (1-d only) (must be set)

-  : is the domain periodic in this direction?
   0 if false, 1 if true (default: 0 0 0)

   Note: an integer is needed for each dimension in the problem

-  castro.center: physical location of problem center on the
   domain (type: Real; default: 0.0 0.0 0.0). The problem
   center is used for gravity, rotation, and some other quantities.
   This is not necessarily the geometric center of the domain—often
   you should choose it to coincide with the center of mass of your
   system. See § \ `[soft:prob_params] <#soft:prob_params>`__ for more details.

   Note: a number is needed for each dimension in the problem

As an example, the following:

::

    geometry.prob_lo = 0 0 0
    geometry.prob_hi = 1.e8 2.e8 2.e8 
    geometry.coord_sys = 0 
    geometry.is_periodic = 0 1 0 
    castro.center = 5.e7 1.e8 1.e8

This defines the domain to run from :math:`(0,0,0)` at the lower left to
:math:`(10^8,\, 2\times 10^8,\, 2\times 10^8)` at the upper right in physical
space, specifies a Cartesian geometry, and makes the domain periodic
in the :math:`y`-direction only. The problem center is set to be halfway in
between the lower left and upper right corners.

Domain Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~

Boundary conditions are specified using integer keys that are interpreted
by BoxLib. The runtime parameters that we use are:

-  castro.lo_bc: boundary type of each low face (must be set)

-  castro.hi_bc: boundary type of each high face (must be set)

The valid boundary types are:

+-------------------------+------------------+--+--+
| 0 – Interior / Periodic | 3 – Symmetry     |  |  |
+-------------------------+------------------+--+--+
| 1 – Inflow              | 4 – Slip Wall    |  |  |
+-------------------------+------------------+--+--+
| 2 – Outflow             | 5 – No Slip Wall |  |  |
+-------------------------+------------------+--+--+

Note: castro.lo_bc and castro.hi_bc must be
consistent with geometry.is_periodic—if the domain is
periodic in a particular direction then the low and high bc’s must be
set to 0 for that direction.

As an example, the following:

::

    castro.lo_bc = 1 4 0 
    castro.hi_bc = 2 4 0 

    geometry.is_periodic = 0 0 1

This defines a problem with inflow (1) in the low-\ :math:`x` direction,
outflow (2) in the high-\ :math:`x` direction, slip wall (4) on
the low and high :math:`y`-faces, and periodic in the :math:`z`-direction.
See § \ `4.2 <#soft:phys_bcs>`__ for more information.

Resolution
~~~~~~~~~~

The grid resolution is specified by defining the resolution at the
coarsest level (level 0) and the number of refinement levels and
factor of refinement between levels. The relevant parameters are:

-  : number of cells in each direction at the
   coarsest level (Integer :math:`> 0`; must be set)

-  : number of levels of refinement above the
   coarsest level (Integer :math:`\geq 0`; must be set)

-  : ratio of coarse to fine grid spacing
   between subsequent levels (2 or 4; must be set)

-  : how often (in terms of number of steps)
   to regrid (Integer; must be set)

-  : should we regrid immediately
   after restarting? (0 or 1; default: 0)

Note: if amr.max_level = 0 then you do not need to set
amr.ref_ratio or amr.regrid_int.

Some examples:

::

    amr.n_cell = 32 64 64

would define the domain to have 32 cells in the :math:`x`-direction, 64 cells
in the :math:`y`-direction, and 64 cells in the :math:`z`-direction *at the
coarsest level*. (If this line appears in a 2D inputs file then the
final number will be ignored.)

::

    amr.max_level = 2 

would allow a maximum of 2 refined levels in addition to the coarse
level. Note that these additional levels will only be created only if
the tagging criteria are such that cells are flagged as needing
refinement. The number of refined levels in a calculation must be
:math:`\leq` amr.max_level, but can change in time and need not
always be equal to amr.max_level.

::

    amr.ref_ratio = 2 4 

would set factor of 2 refinement between levels 0 and 1, and factor of 4
refinement between levels 1 and 2. Note that you must have at least
amr.max_level values of amr.ref_ratio (Additional values
may appear in that line and they will be ignored).

::

    amr.regrid_int = 2 2

tells the code to regrid every 2 steps. Thus in this example, new
level 1 grids will be created every 2 level-0 time steps, and new
level 2 grids will be created every 2 level-1 time steps. If
amr.regrid_int :math:`<` 0 for any level, then regridding starting at that
level will be disabled. If amr.regrid_int = -1 only, then we
never regrid for any level. Note that this is not compatible with
amr.regrid_on_restart = 1.

Regridding
~~~~~~~~~~

The details of the regridding strategy are described in
§ \ `1 <#sec:tagging>`__; here we cover how the input parameters can
control the gridding.

As described later, the user defines Fortran subroutines which tag
individual cells at a given level if they need refinement. This list
of tagged cells is sent to a grid generation routine, which uses the
Berger-Rigoutsos algorithm :raw-latex:`\cite{br-refine}` to create rectangular
grids that contain the tagged cells.

The relevant runtime parameters are:

-  : name of file from which to read the
   grids (text; default: no file)

   If set to a filename, e.g. fixed_girds, then list of grids
   at each fine level are read in from this file during the gridding
   procedure. These grids must not violate the
   amr.max_grid_size criterion. The rest of the gridding procedure
   described below will not occur if amr.regrid_file is set.

-  : radius of additional tagging
   around already tagged cells (Integer :math:`\geq 0`; default: 1)

-  : maximum size of a grid in any
   direction (Integer :math:`> 0`; default: 128 (2-d), 32 (3-d))

   Note: amr.max_grid_size must be even, and a multiple of
   amr.blocking_factor at every level.

-  : grid size must be a multiple of this
   (Integer :math:`> 0`; default: 2)

   Note: amr.blocking_factor at every level must be a power of
   2 and the domain size must be a multiple of
   amr.blocking_factor at level 0.

   This can be very important for elliptic problems with
   multigrid. A higher blocking factor allows the
   multigrid algorithm to coarsen more at the lowest level, reducing
   the amount of work required by the bottom solver.

-  : grid efficiency (Real :math:`>0` and :math:`<1`;
   default: 0.7)

   When creating a refined grid, do we make boxes that only include
   the coarse cells that were explicitly tagged for refinement? or
   do we allow ourselves to encompass nearby, untagged cells in order
   to make larger and more regular boxes? This is the grid efficiency.

   When blocking_factor = 1, *grid efficiency* is exactly the
   fraction of refined cells in the fine BoxArray which correspond to
   coarse cells which were tagged. For other blocking factors,
   we actually apply grid_eff at the level which has been coarsened
   by blocking_factor, so it is no longer strictly this fraction,
   but the idea is still the same.

-  | : refine grids more if # of
     processors :math:`>` # of grids (0 if false, 1 if true; default: 1)

Note also that amr.n_error_buf, amr.max_grid_size and
amr.blocking_factor can be read in as a single value which is
assigned to every level, or as multiple values, one for each level.

As an example, consider:

::

    amr.grid_eff = 0.9
    amr.max_grid_size = 64 
    amr.blocking_factor} = 32

The grid efficiency, amr.grid_eff, means that during the grid
creation process, at least 90% of the cells in each grid at the level
at which the grid creation occurs must be tagged cells. A higher
grid efficiency means fewer cells at higher levels, but may result
in the production of lots of small grids, which have inefficient cache
and OpenMP performance and higher communication costs.

The amr.max_grid_size parameter means that the final grids
will be no longer than 64 cells on a side at every level.
Alternately, we could specify a value for each level of refinement as:
amr.max_grid_size = 64 32 16, in which case our final grids
will be no longer than 64 cells on a side at level 0, 32 cells on a
side at level 1, and 16 cells on a side at level 2. The amr.blocking_factor
means that all of the final grids will be multiples of 32 at all levels.
Again, this can be specified on a level-by-level basis, like
amr.blocking_factor = 32 16 8, in which case the
dimensions of all the final grids will be multiples of 32
at level 0, multiples of 16 at level 1, and multiples of 8 at level 2.

Getting good performance
^^^^^^^^^^^^^^^^^^^^^^^^

These parameters can have a large impact on the performance
of Castro, so taking the time to experiment with is worth the effort.
Having grids that are large enough to coarsen multiple levels in a
V-cycle is essential for good multigrid performance in simulations
that use self-gravity.

 Need more experience here

How grids are created
^^^^^^^^^^^^^^^^^^^^^

The gridding algorithm proceeds in this order:

#. Grids are created using the Berger-Rigoutsos clustering algorithm
   modified to ensure that all new fine grids are divisible by
   amr.blocking_factor.

#. Next, the grid list is chopped up if any grids are larger than max_grid_size.
   Note that because amr.max_grid_size is a multiple of
   amr.blocking_factor the amr.blocking_factor criterion is
   still satisfied.

#. Next, if amr.refine_grid_layout = 1 and there are more processors than grids, and
   if amr.max_grid_size / 2 is a multiple of amr.blocking_factor,
   then the grids will be redefined, at each level independently, so that
   the maximum length of a grid at level :math:`\ell`, in any dimension, is
   amr.max_grid_size[:math:`\ell`] / 2.

#. Finally, if amr.refine_grid_layout = 1, and there are still more processors
   than grids, and if amr.max_grid_size / 4 is a multiple of
   amr.blocking_factor, then the grids will be redefined, at each level
   independently, so that the maximum length of a grid at level :math:`\ell`,
   in any dimension, is amr.max_grid_size[:math:`\ell`] / 4.

Simulation Time
~~~~~~~~~~~~~~~

There are two paramters that can define when a simulation ends:

-  : maximum number of level 0 time steps (Integer
   :math:`\geq 0`; default: -1)

-  : final simulation time (Real :math:`\geq 0`; default:
   -1.0)

To control the number of time steps, you can limit by the maximum
number of level 0 time steps (max_step) or by the final
simulation time (stop_time), or both. The code will stop at
whichever criterion comes first.

Note that if the code reaches stop_time then the final time
step will be shortened so as to end exactly at stop_time, not
past it.

As an example:

::

    max_step  = 1000
    stop_time  = 1.0

will end the calculation when either the simulation time reaches 1.0 or
the number of level 0 steps taken equals 1000, whichever comes first.

Time Step
~~~~~~~~~

If castro.do_hydro = 1, then typically
the code chooses a time step based on the CFL number:

.. math::

   \Delta t = \mathtt{CFL}\, \cdot\, \min_{i,j,k}\left[\min\left\{\frac{\Delta x}{|u|_{i,j,k}+c_{i,j,k}},
                                                                  \frac{\Delta y}{|v|_{i,j,k}+c_{i,j,k}},
                                                                  \frac{\Delta z}{|w|_{i,j,k}+c_{i,j,k}}\right\}\right]
   \label{eq:cfl}

If method-of-lines integration is used instead, then we have

.. math::

   \Delta t = \mathtt{CFL}\, \cdot\, \min_{i,j,k}\left[\left(\frac{\Delta x}{|u|_{i,j,k}+c_{i,j,k}}\right)^{-1} +
                                                       \left(\frac{\Delta y}{|v|_{i,j,k}+c_{i,j,k}}\right)^{-1} +
                                                       \left(\frac{\Delta z}{|w|_{i,j,k}+c_{i,j,k}}\right)^{-1}\right]^{-1}

(If we are simulating in 1D or 2D, the extraneous parts related to :math:`v` and/or :math:`w` are removed.)

The following parameters affect the timestep choice:

-  : CFL number (Real :math:`> 0` and :math:`\leq 1`;
   default: 0.8)

-  : factor by which to shrink the initial
   time step (Real :math:`> 0` and :math:`\leq 1`; default: 1.0)

-  : factor by which the time step can
   grow in subsequent steps (Real :math:`\geq 1`; default: 1.1)

-  : level 0 time step regardless of cfl
   or other settings (Real :math:`> 0`; unused if not set)

-  : initial level 0 time
   step regardless of other settings (Real :math:`> 0`; unused if not set)

-  : time step below which calculation
   will abort (Real :math:`> 0`; default: 0.0)

-  : whether or not to abort the
   simulation if the hydrodynamics update creates velocities that
   violate the CFL criterion (Integer; default: 1)

As an example, consider:

::

    castro.cfl = 0.9 
    castro.init_shrink = 0.01 
    castro.change_max = 1.1
    castro.dt_cutoff = 1.e-20

This defines the :math:`\mathtt{cfl}` parameter in Eq. \ `[eq:cfl] <#eq:cfl>`__ to be
0.9, but sets (via init_shrink) the first timestep we take to
be 1% of what it would be otherwise. This allows us to ramp up to
the hydrodynamic timestep at the start of a simulation. The
change_max parameter restricts the timestep from increasing by
more than 10% over a coarse timestep. Note that the time step can
shrink by any factor; this only controls the extent to which it can
grow. The dt_cutoff parameter will force the code to abort if
the timestep ever drops below :math:`10^{-20}`. This is a safety
feature—if the code hits such a small value, then something likely
went wrong in the simulation, and by aborting, you won’t burn through
your entire allocation before noticing that there is an issue.

If we know what we are doing, then we can force a particular timestep:

::

    castro.fixed_dt = 1.e-4

This sets the level 0 time step to be 1.e-4 for the entire simulation,
ignoring the other timestep controls. Note that if
castro.init_shrink :math:`\neq 1` then the first time step will in fact
be castro.init_shrink :math:`\cdot` castro.fixed_dt.

::

    castro.initial_dt = 1.e-4

sets the *initial* level 0 time step to be :math:`10^{-4}` regardless of
castro.cfl or castro.fixed_dt. The time step can
grow in subsequent steps by a factor of castro.change_max each step.

[] If diffusion is enabled, the timestep will also
be limited by:

.. math::

   \Delta t = \frac{1}{2}\min_{i,j,k}\left[\min\left\{\frac{\Delta x^2}{D_{i,j,k}},
                                                      \frac{\Delta y^2}{D_{i,j,k}},
                                                      \frac{\Delta z^2}{D_{i,j,k}}\right\}\right]

where :math:`D \equiv k / (\rho c_V)` if we are diffusing temperature, and
:math:`D \equiv k / (\rho c_P)` if we are diffusing enthalpy. No input parameter
is necessary to enable this constraint. See Chapter `[ch:diffusion] <#ch:diffusion>`__ for more details.

[] If reactions are enabled, the timestep will also
be limited by two constraints:

.. math:: \Delta t = \mathtt{dtnuc\_e}\, \min_{i,j,k} \left\{\frac{e_{i,j,k}}{\dot{e}_{i,j,k}}\right\}

.. math:: \Delta t = \mathtt{dtnuc\_X}\, \min_{i,j,k} \left\{\min_n\frac{X^n_{i,j,k}}{\dot{X}^n_{i,j,k}}\right\}

where :math:`e` is the internal energy, and :math:`X^n` is the mass fraction of
the :math:`n`\ th species. The safety factors correspond to the runtime parameters
and . These limiters
say that the timestep must be small enough so that no zone can change
its internal energy by more than the fraction in one
step, and so that no zone can change the abundance of any isotope by
more than the fraction in one step. The time derivatives
:math:`\dot{e}` and :math:`\dot{X}^n` are estimated by calling the right-hand-side
of the nuclear network given the state at the time the timestep limiter
is being calculated. (We use a small number floor to prevent division by zero.)
To prevent the timestep from being dominated by trace species, there is
an additional option which is the
mass fraction threshold below which a species will not be considered in
the timestep constraint. and are set to
a large number by default, effectively disabling them. Typical choices
for these values in the literature are :math:`\sim 0.1`.

Subcycling
~~~~~~~~~~

Castro supports a number of different modes for subcycling in time,
set via .

-  amr.subcycling_mode = Auto (default): the code will run
   with equal refinement in space and time. In other words, if level
   :math:`n+1` is a factor of 2 refinement above level :math:`n`, then :math:`n+1` will
   take 2 steps of half the duration for every level :math:`n` step.

-  If amr.subcycling_mode = None: the code will not refine
   in time. All levels will advance together with a timestep dictated
   by the level with the strictest :math:`dt`. Note that this is identical to
   the deprecated command amr.nosub = 1.

-  If amr.subcycling_mode = Manual: the code will subcycle
   according to the values supplied by .

In the case of amr.subcycling_mode = Manual, we subcycle in
manual mode with largest allowable timestep. The number of iterations
at each level is then specified as:

::

    amr.subcycling_iterations = 1 2 1 2

Here, we take 1 level-0 timestep at a time (required). Take 2 level-1
timesteps for each level-0 step, 1 timestep at level-2 for each
level-1 step, and take 2 timesteps at level-3 for each level-2 step.

Alternately, we could do:

::

    amr.subcycling_iterations = 2

which will subcycle twice at every level (except level 0).

Restart Capability
~~~~~~~~~~~~~~~~~~

Castro has a standard sort of checkpointing and restarting capability.
In the inputs file, the following options control the generation of
checkpoint files (which are really directories):

-  : prefix for restart files (text;
   default: chk)

-  : how often (by level 0 time steps) to
   write restart files (Integer :math:`> 0`; default: -1)

-  : how often (by simulation time) to
   write restart files (Real :math:`> 0`; default: -1.0)

   Note that amr.check_per will write a checkpoint at the first
   timestep whose ending time is past an integer multiple of this interval.
   In particular, the timestep is not modified to match this interval, so
   you won’t get a checkpoint at exactly the time you requested.

-  : name of the file (directory) from
   which to restart (Text; not used if not set)

-  : should we write
   checkpoint files? (0 or 1; default: 1)

   If you are doing a scaling study then set
   amr.checkpoint_files_output = 0 so you can test scaling of the
   algorithm without I/O.

-  : how parallel is the writing of
   the checkpoint files? (Integer :math:`\geq 1`; default: 64)

   See the § \ `9 <#software:io>`__ for more details on parallel I/O and the
   amr.check_nfiles parameter.

-  : should we write a
   checkpoint immediately after restarting? (0 or 1; default: 0)

-  : factor by which domain has been
   grown (Integer :math:`\geq 1`; default: 1)

Note:

-  You can specify both amr.check_int or amr.check_per,
   if you so desire; the code will print a warning in case you did this
   unintentionally. It will work as you would expect – you will get checkpoints
   at integer multiples of amr.check_int timesteps and at integer
   multiples of amr.check_per simulation time intervals.

-  amr.plotfile_on_restart and
   amr.checkpoint_on_restart require amr.regrid_on_restart
   to be in effect.

As an example,

::

    amr.check_file = chk_run
    amr.check_int = 10

means that restart files (really directories) starting with the prefix
“chk_run” will be generated every 10 level-0 time steps. The
directory names will be chk_run00000, chk_run00010,
chk_run00020, etc.

If instead you specify

::

    amr.check_file = chk_run
    amr.check_per = 0.5

then restart files (really directories) starting with the prefix
“chk_run” will be generated every 0.1 units of
simulation time. The directory names will be chk_run00000,
chk_run00043, chk_run00061, etc, where :math:`t = 0.1` after
43 level-0 steps, :math:`t = 0.2` after 61 level-0 steps, etc.

To restart from chk_run00061, for example, then set

::

    amr.restart = chk_run00061

.. _sec:PlotFiles:

Controlling Plotfile Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main output from Castro is in the form of plotfiles (which are
really directories). The following options in the inputs file control
the generation of plotfiles:

-  : prefix for plotfiles (text; default:
   “plt”)

-  : how often (by level-0 time steps) to
   write plot files (Integer :math:`> 0`; default: -1)

-  : how often (by simulation time) to write
   plot files (Real :math:`> 0`; default: -1.0)

   Note that amr.plot_per will write a plotfile at the first
   timestep whose ending time is past an integer multiple of this interval.
   In particular, the timestep is not modified to match this interval, so
   you won’t get a checkpoint at exactly the time you requested.

-  : name of state variables to include in
   plotfiles (valid options: ALL, NONE or a list; default:
   ALL)

-  : name of derived variables to
   include in plotfiles (valid options: ALL, NONE or a
   list; default: NONE

-  : should we write plot files?
   (0 or 1; default: 1)

   If you are doing a scaling study then set
   amr.plot_files_output = 0 so you can test scaling of the
   algorithm without I/O.

-  : should we write a plotfile
   immediately after restarting? (0 or 1; default: 0)

-  : how parallel is the writing of the
   plotfiles? (Integer :math:`\geq 1`; default: 64)

   See the Software Section for more details on parallel I/O and the
   amr.plot_nfiles parameter.

-  : include all the species mass
   fractions in the plotfile (0 or 1; default: 0)

All the options for amr.derive_plot_vars are kept in
``derive_lst`` in Castro_setup.cpp. Feel free to look at
it and see what’s there.

Some notes:

-  You can specify both amr.plot_int or amr.plot_per,
   if you so desire; the code will print a warning in case you did this
   unintentionally. It will work as you would expect – you will get plotfiles
   at integer multiples of amr.plot_int timesteps and at integer
   multiples of amr.plot_per simulation time intervals.

As an example:

::

    amr.plot_file = plt_run
    amr.plot_int = 10

means that plot files (really directories) starting with the prefix
“plt_run” will be generated every 10 level-0 time steps. The
directory names will be plt_run00000, plt_run00010,
plt_run00020, etc.

If instead you specify

::

    amr.plot_file = plt_run
    amr.plot_per = 0.5

then restart files (really directories) starting with the prefix
“plt_run” will be generated every 0.1 units of simulation time. The
directory names will be plt_run00000, plt_run00043,
plt_run00061, etc, where :math:`t = 0.1` after 43 level-0 steps, :math:`t =
0.2` after 61 level-0 steps, etc.

Screen Output
~~~~~~~~~~~~~

There are several options that set how much output is written to the
screen as Castro runs:

-  : verbosity of Amr.cpp (0 or 1; default: 0)

-  : verbosity of Castro.cpp (0 or 1; default: 0)

-  : verbosity of Gravity.cpp (0 or 1; default: 0)

-  : verbosity of Diffusion.cpp (0 or 1;
   default: 0)

-  : verbosity of multigrid solver (for gravity) (allow
   values: 0,1,2,3,4; default: 0)

-  : name of the file to which the grids are
   written (text; not used if not set)

-  : name of the file to which certain output is
   written (text; not used if not set)

-  : name of the file to which certain
   (terser) output is written (text; not used if not set)

-  : if :math:`> 0`, how often (in level-0 time
   steps) to compute and print integral quantities (Integer; default: -1)

   The integral quantities include total mass, momentum and energy in
   the domain every castro.sum_interval level-0 steps.
   The print statements have the form

   ::

           TIME= 1.91717746 MASS= 1.792410279e+34
         

   for example. If this line is commented out then
   it will not compute and print these quanitities.

-  : allows the user to set a
   special flag based on user-specified criteria (0 or 1; default: 1)

   castro.do_special_tagging = 1 can be used, for example, to
   calculate the bounce time in a core collapse simulation; the bounce
   time is defined as the first time at which the maximum density in
   the domain exceeds a user-specified value. This time can then be
   printed into a special file as a useful diagnostic.

As an example:

::

    amr.grid_log = grdlog
    amr.run_log = runlog 

Every time the code regrids it prints a list of grids at all relevant
levels. Here the code will write these grids lists into the file
grdlog. Additionally, every time step the code prints certain
statements to the screen (if amr.v = 1), such as:

::

    STEP = 1 TIME = 1.91717746 DT = 1.91717746 
    PLOTFILE: file = plt00001 

The run_log option will output these statements into
*runlog* as well.

Terser output can be obtained via:

::

    amr.run_log_terse} = runlogterse

This file, runlogterse differs from runlog, in that it
only contains lines of the form

::

    10  0.2  0.005

in which “10” is the number of steps taken, “0.2” is the
simulation time, and “0.005” is the level-0 time step. This file
can be plotted very easily to monitor the time step.

Other parameters
~~~~~~~~~~~~~~~~

There are a large number of solver-specific runtime parameters. We describe these
together with the discussion of the physics solvers in later chapters.

Software Framework
==================

Code structure
--------------

Castro is built upon the AMReX C framework. This provides
high-level classes for managing an adaptive mesh refinement
simulation, including the core data structures we will deal with. A
key design pattern in AMReX is that the overall memory management
and parallelization is done in the C layer, while the heavy
computational work is done in Fortran kernels. AMReX provides
convenient data structures that allow for this workflow—high level
objects in C that communicate with Fortran through pointers to
data regions that appear as multidimensional arrays.

Castro uses a structured-grid approach to hydrodynamics. We work
with square/cubic zones that hold state variables (density, momentum,
etc.) and compute the fluxes of these quantities through the
interfaces of the zones (this is a finite-volume approach).
Parallelization is achieved by domain decomposition. We divide our
domain into many smaller boxes, and distributed these across
processors. When information is needed at the boundaries of the
boxes, messages are exchanged and this data is held in a perimeter of
*ghost cells*. AMReX manages this decompostion and
communication for us. Additionally, AMReX implements adaptive mesh
refinement. In addition to the coarse decomposition of our domain
into zones and boxes, we can refine rectangular regions by adding
finer-gridded boxes on top of the coarser grid. We call the
collection of boxes at the same resolution a *level*.

Castro uses a hybrid MPI + OpenMP approach to parallelism. MPI is
at used to communicate across nodes on a computer and OpenMP is used
within a node, to loop over subregions of a box with different
threads.

The code structure in the Castro/ directory reflects the
division between C and Fortran.

-  constants/: contains a file of useful constants in CGS units

-  Docs/: you’re reading this now!

-  Exec/: various problem implementations, sorted by category:

   -  gravity_tests/: test problems that primarily exercise the gravity solver

   -  hydro_tests/: test problems of the hydrodynamics (with or without reactions)

   -  radiation_tests/: test problems that primarily exercise the radiation hydrodynamics solver

   -  science/: problem setups that were used for scientific investigations

   -  unit_tests/: test problems that exercise primarily a single module

-  Microphysics/: contains directories for different default
   microphysics (these are all implemented in Fortran)

   -  conductivity/: the thermal conductivity

   -  EOS/: the equation of state

   -  networks/: the nuclear reaction networks

   -  opacity/: the radiative opacity (used with radiation)

   -  viscosity/: the viscous transport coefficient

-  Source/: source code. In this main directory is all of
   the code. Sources are mixed C and Fortran and are organized by topic as:

   -  diffusion/ : thermal diffusion code

   -  driver/ : the main driver, I/O, runtime parameter support

   -  gravity/ : self-gravity code

   -  hydro/ : the compressible hydrodynamics code

   -  particles/ : support for particles

   -  problems/ : template code for implementing a problem

   -  radiation/ : the implicit radiation solve code

   -  reactions/ : nuclear reaction code

   -  rotation/ : rotating code

   -  sources/ : hydrodynamics source terms support

-  Util/: a catch-all for additional things you may need

   -  ConvertCheckpoint/: a tool to convert a checkpoint file to
      a larger domain

   -  :math:`\ldots`

Major data structures
---------------------

The following data structures are the most commonly encountered when
working in the C portions of Castro. This are all
AMReX data-structures / classes.

Amr
~~~

This is the main class that drives the whole simulation. This is
the highest level in Castro.

AmrLevel and Castro classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

An is a virtual base class provided by AMReX that
stores all the state data on a single level in the AMR hierarchy and
understands how to advance that data in time.

The most important data managed by the AmrLevel is an array of
StateData, which holds the fluid quantities, etc., in the boxes
that together make up the level.

The Castro class is derived from the AmrLevel. It provides
the Castro-specific routines to evolve our system of equations. Like
the AmrLevel, there is one Castro object for each level in the
AMR hierarchry.

A lot of the member data in the Castro class are static member
variables—this means that they are shared across all instances of
the class. So, in this case, every level will have the same data.
This is done, in particular, for the values of the runtime parameters,
but also for the Gravity, Diffusion, and Radiation
objects. This means that those objects cover all levels and are the
same object in each instantiation of the Castro class.

Floating point data
~~~~~~~~~~~~~~~~~~~

Floating point data in the C AMReX frame work is declared as
Real. This is typedef to either float or
double depending on the make variable .

The corresponding type for Fortran is provided by the
as . We typically rename
this to when using it. An example of a declaration of a
parameter is:

::

      use amrex_fort_module, only : rt => amrex_real                                       

      real(rt) :: tol = 1.0e-10_rt

The bl_constants_module provides common constants that can
be used in the code, like ZERO, THIRD, ONE, etc.

Note: single precision support in Castro is not yet complete. In
particular, a lot of the supporting microphysics has not been updated.

Box and FArrayBox
~~~~~~~~~~~~~~~~~

A is simply a rectangular region in space. It does not hold
data. In AMReX, an AMR level has a global index space, with
:math:`(0,0,0)` being the lower left corner of the domain at that level, and
:math:`(N_x-1, N_y-1, N_z-1)` being the upper right corner of the domain
(for a domain of :math:`N_x \times N_y \times N_z` zones). The location of
any Box at a level can be uniquely specified with respect to this
global index space by giving the index of its lower-left and
upper-right corners. Figure \ `[fig:soft:indexspace] <#fig:soft:indexspace>`__ shows an
example of three boxes at the same level of refinement.

AMReX provides other data structures that collect Boxes together,
most importantly the . We generally do not use these
directly, with the exception of the BoxArray ,
which is defined as part of the AmrLevel class that Castro
inherits. grids is used when building new MultiFabs to give
the layout of the boxes at the current level.

.. raw:: latex

   \centering

.. figure:: index_grid2
   :alt: [fig:soft:indexspace] Three boxes that comprise a single level. At this
   resolution, the domain is 20\ :math:`\times`\ 18 zones. Note that the
   indexing in AMReX starts with :math:`0`.
   :width: 4in

   [fig:soft:indexspace] Three boxes that comprise a single level. At this
   resolution, the domain is 20\ :math:`\times`\ 18 zones. Note that the
   indexing in AMReX starts with :math:`0`.

A or *FAB*, for *Fortran array box* is a data
structure that contains a Box locating it in space, as well as a
pointer to a data buffer. The real floating point data are stored as
one-dimensional arrays in FArrayBoxes. The associated Boxcan be
used to reshape the 1D array into multi-dimensional arrays to be used
by Fortran subroutines. The key part of the C AMReX data
structures is that this data buffer can be sent to Fortran, where it
will appear as a DIM+1 dimensional array (DIM space + 1
component).

Note: Castro is complied for a specific dimensionality.

MultiFab
~~~~~~~~

At the highest abstraction level, we have the (mulitple
FArrayBoxes). A MultiFab contains an array of Boxes, including
Boxes owned by other processors for the purpose of communication,
an array of MPI ranks specifying which MPI processor owns each Box,
and an array of pointers to FArrayBoxes owned by this MPI
processor. Note: a
MultiFab is a collection of the boxes that together make up a single
level of data in the AMR hierarchy.

A MultiFab can have multiple components (like density, temperature,
...) as well as a perimeter of ghost cells to exchange data with
neighbors or implement boundary conditions (this is all reflected in
the underlying FArrayBox).

Parallelization in AMReX is done by distributing the FABs across
processors. Each processor knows which FABs are local to it. To loop
over all the boxes local to a processor, an MFIter is used (more
on this below).

High-level operations exist on MultiFabs to add, subtract, multiply,
etc., them together or with scalars, so you don’t need to write out
loops over the data directly.

In Castro, MultiFabs are one of the main data structures you will
interact with in the C portions of the code.

.. _soft:sec:statedata:

StateData
~~~~~~~~~

is a class that essentially holds a pair of MultiFabs: one
at the old time and one at the new time. AMReX knows how to
interpolate in time between these states to get data at any
intermediate point in time. The main data that we care about in
Castro (the fluid state, gravitational potential, etc.) will be
stored as StateData. Essentially, data is made StateData in
Castro if we need it to be stored in checkpoints / plotfiles, and/or
we want it to be automatically interpolated when we refine.

An AmrLevel stores an array of StateData (in a C array
called state). We index this array using integer keys (defined
via an enum in Castro.H). The state data is registered
with AMReX in .

Note that each of the different StateData carried in the state
array can have different numbers of components, ghost cells, boundary
conditions, etc. This is the main reason we separate all this data
into separate StateData objects collected together in an indexable
array.

The current StateData names Castro carries are:

-  : this is the NUM_STATE hydrodynamics
   components that make up the conserved hydrodynamics state (usually
   referred to as :math:`{\bf U}` in these notes. But note that this does
   not include the radiation energy density.

   In Fortran, the components of a FAB derived from State_Type
   is indexed using the integer keys defined in
   and stored in , e.g., URHO, UMX,
   UMY, ...

   Note: regardless of dimensionality, we always carry around all
   three velocity components. The “out-of-plane” components
   will simply be advected, but we will allow rotation (in particular,
   the Coriolis force) to affect them.

   State_Type MultiFabs have no ghost cells by default for
   pure hydro and a single ghost cell by default when RADIATION
   is enabled. There is an option to force them to have ghost cells by
   setting the parameter at runtime.

   Note that the prediction of the hydrodynamic state to the interface
   will require 4 ghost cells. This accomodated by creating a separate
   MultiFab, that lives at the old-time level and
   has the necessary ghost cells. We will describe this more later.

-  : this stores the radiation energy density,
   commonly denoted :math:`E_r` in these notes. It has
   components—the number of energy groups used in the multigroup
   radiation hydrodynamics approximation.

-  : this is simply the gravitational
   potential, usually denoted :math:`\Phi` in these notes.

-  : this is the gravitational
   acceleration. There are always 3 components, regardless of the
   dimensionality (consistent with our choice of always carrying all 3
   velocity components).

-  : this is the rotational potential.
   When rotation is enabled, this will store the effective potential
   corresponding to the centrifugal force.

-  : this is the rotational acceleration.
   There are always 3 components, regardless of the dimensionality
   (consistent with our choice of always carrying all 3 velocity
   components). This includes the terms corresponding to the Coriolis
   force, the centrifugal force, as well as optional terms due to the
   change in rotation rate, :math:`\Omega`.

-  : this holds the time-rate of change of
   the source terms, :math:`d{\bf S}/dt`, for each of the NUM_STATE
   State_Type variables.

   .. raw:: latex

      \marginpar{\vskip-\baselineskip\raggedright\tiny\sffamily
      \hrule\smallskip{\color{red}SDC does differently}\par\smallskip\hrule}

   Note: we do not make use of the old-time quantity here. In fact, we
   never allocate the FArrayBoxs for the old-time in the Source_Type
   StateData, so there is not wasted memory.

-  : this holds the data for the nuclear
   reactions. It has NumSpec+2 components: the species
   creation rates (usually denoted :math:`\dot\omega_k` in these notes),
   the specific energy generation rate (:math:`\dot{e}_\mathrm{nuc}`),
   and its density (:math:`\rho \dot{e}_\mathrm{nuc}`).

   These are stored as StateData so we have access to the reaction terms
   outside of advance, both for diagnostics (like flame speed estimation)
   and for reaction timestep limiting (this in particular needs the
   data stored in checkpoints for continuity of timestepping upon restart).

   .. raw:: latex

      \marginpar{\vskip-\baselineskip\raggedright\tiny\sffamily
      \hrule\smallskip{\color{red}why do we need rho edot and edot separately?}\par\smallskip\hrule}

-  : this is used with the SDC
   time-advancement algorithm. This stores the QVAR terms
   that describe how the primitive variables change over the timestep
   due only to reactions. These are used when predicting the interface
   states of the primitive variables for the hydrodynamics portion of the
   algorithm.

We access the multifabs that carry the data of interest by interacting
with the StateData using one of these keys. For instance:

::

    MultiFab& S_new = get_new_data(State_Type);

gets a pointer to the multifab containing the hydrodynamics state data
at the new time.

Various source MultiFabs
~~~~~~~~~~~~~~~~~~~~~~~~

There are a number of different MultiFabs (and arrays of MultiFabs)
that hold source term information.

-  : this is a MultiFab that holds the
   update to the hydrodynamics (basically the divergence of the
   fluxes). This is filled in the conservative update routine of the
   hydrodynamics.

   As this is expressed as a source term, what is actually stored is

   .. math:: {\bf S}_\mathrm{flux} = -\nabla \cdot {\bf F}

   So the update of the conserved state appears as:

   .. math:: \frac{\partial {\bf U}}{\partial t} = {\bf S}_\mathrm{flux}

-  : a single MultiFab that stores
   the sum of sources over each physical process.

MFIter and interacting with Fortran
-----------------------------------

The process of looping over boxes at a given level of refinement and
operating on their data in Fortran is linked to how Castro achieves
thread-level parallelism. The OpenMP approach in Castro has evolved
considerably since the original paper was written, with the modern
approach, called *tiling*, gearing up to meet the demands of
many-core processors in the next-generation of supercomputers. We
discuss the original and new approach together here.

In both cases, the key construct is the —this is a
C iterator that knows how to loop over the FArrayBoxes in the
MultiFab that are local to the processor (in this way, a lot of the
parallelism is hidden from view).

Non-Tiling MFIter
~~~~~~~~~~~~~~~~~

The non-tiling way to iterate over the FArrayBoxs is
 [4]_:

.. code:: c++

      for (MFIter mfi(mf); mfi.isValid(); ++mfi) // Loop over boxes
      {
        // Get the index space of this iteration
        const Box& box = mfi.validbox();

        // Get a reference to the FAB, which contains data and box
        FArrayBox& fab = mf[mfi];

        // Get the index space for the data region in th FAB.
        // Note "abox" may have ghost cells, and is thus larger than
        // or equal to "box" obtained using mfi.validbox().
        const Box& abox = fab.box();

        // We can now pass the information to a Fortran routine,
        // fab.dataPtr() gives a double*, which is reshaped into
        // a multi-dimensional array with dimensions specified by
        // the information in "abox". We will also pass "box",
        // which specifies our "work" region.
        do_work(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
                fab.dataPtr(), fab.nComp(),
                ARLIM_3D(abox.loVect()), ARLIM_3D(abox.hiVect())

      }

A few comments about this code

-  In this example, we are working off of a MultiFab named mf.
   This could, for example, come from state data as:

   ::

        MultiFab& mf = get_old_data(State_Type);

-  We are passing the data in mf one box at a time to the
   Fortran function do_work.

-  Here the MFIter iterator, mfi, will perform the loop
   only over the boxes that are local to the MPI task. If there are 3
   boxes on the processor, then this loop has 3 iterations.

   ++mfi iterates to the next FArrayBox owned by the
   MultiFab mf, and mfi.isValid() returns false
   after we’ve reached the last box contained in the MultiFab,
   terminating the loop.

-  box as returned from mfi.validbox() does not include
   ghost cells. This is the valid data region only.
   We can get the indices of the valid zones as box.loVect() and
   box.hiVect().

   In passing to the Fortran function, we use the macro
   , defined in to pass the lo
   and hi vectors as pointers to an int array. This array
   is defined to always be 3D, with 0s substituted for the
   higher dimension values if we are running in 1- or 2D.

   Passing the data in this 3D fashion is a newer approach in Castro.
   This enables writing *dimension agnostic code*. There are many
   other approaches that will pass only the DIM values of
   lo and hi using alternate macros in ArrayLim.H.

-  fab.dataPtr() returns a double \*—a pointer to the
   data region. This is what is passed to Fortran.

-  fab.nComp() gives an int—the number of components
   in the MultiFab. This will be used for dimensioning in Fortran.

-  To properly dimension the array in Fortran, we need the actual
   bounds of the data region, including any ghost cells. This is the
   Box abox, obtained as fab.box(). We pass the
   lo and hi of the full data region as well.

To properly compile, we need a prototype for the Fortran
function. These are placed in the \_F.H files in the
Castro Source/ directory. Here’s the prototype for
our function:

.. code:: c++

      void do_work
        (const int* lo, const int* hi,
         Real* state, const Real& ncomp
         const int* s_lo, const int* s_hi)

A few comments on the prototype:

-  we use the const qualifier on the many of the arguments.
   This indicates that the data that is pointed to cannot be
   modified [5]_
   means that the pointers themselves are to be unmodified. But the
   contents of the memory space that they point to can be modified.

-  For ncomp, we in the calling sequence, we just did
   fab.nComp(). This returns a int. But Fortran is a
   pass-by-reference language, so we make the argument in the prototype
   a reference. This ensures that it is passed by reference.

In our Fortran example, we want to loop over all of the data,
including 1 ghost cell all around. The corresponding Fortran function
will look like:

.. code:: fortran

      subroutine do_work(lo, hi, &
                         state, ncomp, &
                         s_lo, s_hi) bind(C, name="do_work")

        use prob_params_module, only : dg

        integer, intent(in) :: lo(3), hi(3)
        integer, intent(in) :: s_lo(3), s_hi(3), ncomp

        real (kind=dp_t), intent(inout) :: state(s_lo(1):s_hi(1), &
                                                 s_lo(2):s_hi(2), &
                                                 s_lo(3):s_hi(3), ncomp)

        ! loop over the data
        do k = lo(3)-1*dg(3), hi(3)+1*dg(3)
           do j = lo(2)-1*dg(2), hi(2)+1*dg(2)
              do i = lo(1)-1*dg(1), hi(1)+1*dg(1)

                 ! work on state(i,j,k,:), where the last index
                 ! is the component of the multifab

              enddo
           enddo
        enddo

      end subroutine do_work

Finally, comments on the Fortran routine;

-  We use the Fortran 2003 bind keyword to specify
   that we want this to be interoperable with C. Ordinarily
   we would not need to specify the optional argument name
   in the binding, but the PGI compiler requires this if our
   Fortran subroutine is part of a module.

-  We dimension state using s_lo and s_hi—these are
   the bounds we got from the FArrayBox, and are for the entire data
   region, including ghost cells.

   Note, in Fortran, the spatial indices of state don’t
   necessarily start at 1—they reflect the global index space
   for the entire domain at this level of refinement. This means that
   we know where the box is located.

   Later we’ll see how to compute the spatial coordinates using this
   information.

-  Our loop uses lo and hi—these are the indices
   of the valid data region (no ghost cells). Since we want a single
   ghost cell all around, we subtract 1 from lo and add 1
   to hi.

   Finally, since this is dimension-agnostic code (it should work
   correctly in 1-, 2-, and 3D), we need to ensure the loops over the
   higher dimensions do nothing when we compile for a lower
   dimensionality. This is the role of dg—dg is 1
   if our simulation includes that spatial dimension and 0
   otherwise.

   If we were not looping over ghost cells too, then we would not need
   to invoke dg, since lo and hi are both set to
   0 for any dimensions not represented in our simulation.

Up to this point, we have not said anything about threading. In this
style of using the MFIter, we implement the OpenMP in Fortran, for
instance by putting a pragma around the outer loop in this example.

.. _sec:boxlib1:

AMReX’s Current Tiling Approach In C++
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two types of tiling that people discuss. In *logical
tiling*, the data storage in memory is unchanged from how we do things
now in pure MPI. In a given box, the data region is stored
contiguously). But when we loop in OpenMP over a box, the tiling
changes how we loop over the data. The alternative is called
*separate tiling*—here the data storage in memory itself is changed
to reflect how the tiling will be performed. This is not considered
in AMReX.

We have recently introduced logical tiling into parts of AMReXİt
is off by default, to make the transition smooth and because not
everything should be tiled. It can be enabled on a loop-by-loop basis
by setting an optional argument to MFIter. We demonstrate this
below. Further examples can be found at Tutorials/Tiling_C,
and Src/LinearSolvers/C_CellMG/.

In our logical tiling approach, a box is logically split into tiles,
and a MFIter loops over each tile in each box. Note that the
non-tiling iteration approach can be considered as a special case of
tiling with the tile size equal to the box size.

Let us consider an example. Suppose there are four boxes—see
Figure \ `[fig:domain-tiling] <#fig:domain-tiling>`__.

.. raw:: latex

   \centering

.. figure:: domain-tile
   :alt: [fig:domain-tiling] A simple domain showing 4
   Boxes labeled 0–3, and their tiling regions (dotted lines)

   [fig:domain-tiling] A simple domain showing 4
   Boxes labeled 0–3, and their tiling regions (dotted lines)

The first box is divided into 4 logical tiles, the second and third
are divided into 2 tiles each (because they are small), and the fourth
into 4 tiles. So there are 12 tiles in total. The difference between
the tiling and non-tiling version are then:

-  In the tiling version, the loop body will be run 12 times. Note
   that tilebox is different for each tile, whereas fab
   might be referencing the same object if the tiles belong to the same
   box.

-  In the non-tiling version (by constructing MFIter without
   the optional second argument or setting to false), the loop
   body will be run 4 times because there are four boxes, and a call to
   mfi.tilebox() will return the traditional validbox. The
   non-tiling case is essentially having one tile per box.

The tiling implementation of the same call to our the Fortran
do_work routine is show below:

.. code:: c++

      bool tiling = true;
      for (MFIter mfi(mf, tiling); mfi.isValid(); ++mfi) // Loop over tiles
      {
        // Get the index space of this iteration.
        const Box& box = mfi.growntilebox(1);

        // Get a reference to the FAB, which contains data and box
        FArrayBox& fab = mf[mfi];

        // Get the index space for the data pointed by the double*.
        const Box& abox = fab.box();

        // We can now pass the information to a Fortran routine.
        do_work(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
                fab.dataPtr(), fab.nComp(),
                ARLIM_3D(abox.loVect()), ARLIM_3D(abox.hiVect())

      }

Note that the code is almost identical to the one in § \ `[sec:boxlib0] <#sec:boxlib0>`__.
Some comments:

-  The iterator now takes an extra argument to turn on tiling (set
   to true).

   There is another interface fo MFIter that can take an
   IntVect that explicitly gives the tile size in each coordinate
   direction. If we don’t explictly specify the tile size at the loop,
   then the runtime parameter
   can be used to set it globally.

-  .validBox() has the same meaning as in the non-tile
   approach, so we don’t use it.
   Since in this example, we want to include a single ghost cell in our
   loop over the data, we use .growntilebox(1) (where the 1
   here indicates a single ghost cells) to get the Box (and
   corresponding lo and hi) for the *current tile*, not
   the entire data region. If instead, we just wanted the valid
   region in Fortran, without any ghost cells, we would use
   .tilebox().

-  When passing into the Fortran routine, we still use the index
   space of the entire FArrayBox (including ghost cells), as seen in
   the abox construction. This is needed to properly dimension
   the array in Fortran.

   The Fortran routine will declare a multidimensional array that is of
   the same size as the entire box, but only work on the index space
   identified by the tile-box (box).

The Fortran code is almost the same as before, but now our loop
simply uses lo and hi, since, by construction with
.growntilebox(1), this already includes the single ghost cell
all around:

.. code:: fortran

      subroutine do_work(lo, hi, &
                         state, ncomp, &
                         s_lo, s_hi) bind(C, name="do_work")

        integer, intent(in) :: lo(3), hi(3)
        integer, intent(in) :: s_lo(3), s_hi(3), ncomp

        real (kind=dp_t), intent(inout) :: state(s_lo(1):s_hi(1), &
                                                 s_lo(2):s_hi(2), &
                                                 s_lo(3):s_hi(3), ncomp)

        ! loop over the data
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)

                 ! work on state(i,j,k,:), where the last index
                 ! is the component of the multifab

              enddo
           enddo
        enddo

      end subroutine do_work

The function prototype is unchanged.

Tiling provides us the opportunity of a coarse-grained approach for
OpenMP. Threading can be turned on by inserting the following line
above the for (MFIter...) line.

::

      #pragma omp parallel

Note that the OpenMP pragma does not have a for—this is not
used when working with an iterator.

Assuming four threads are used in the above example, thread 0 will
work on 3 tiles from the first box, thread 1 on 1 tile from the first
box and 2 tiles from the second box, and so forth. Note that
OpenMP can be used even when tiling is turned off. In that case, the
OpenMP granularity is at the box level (and good performance would need
many boxes per MPI task).

The tile size for the three spatial dimensions can be set by a
parameter, e.g., fabarray.mfiter_tile_size = 1024000 8 8. A
huge number like 1024000 will turn off tiling in that direction.
As noted above, the MFIter constructor can also take an explicit
tile size: MFIter(mfi(mf,IntVect(128,16,32))).

Note that tiling can naturally transition from all threads working
on a single box to each thread working on a separate box as the boxes
coarsen (e.g., in multigrid).

The MFIter class provides some other useful functions:

-  mfi.validbox() : The same meaning as before independent of tiling.

-  mfi.tilebox() : The standard way of getting the bounds of the
   current tile box. This will tile over the valid data region only.

-  mfi.growntilebox(int) : A grown tile box that includes
   ghost cells at box boundaries only. Thus the returned boxes for a
   FArrayBox are non-overlapping.

-  mfi.nodaltilebox(int) : Returns non-overlapping
   edge-type boxes for tiles. The argument is for direction.

-  mfi.fabbox() : Same as mf[mfi].box().

Finally we note that tiling is not always desired or better. The
traditional fine-grained approach coupled with dynamic scheduling is
more appropriate for work with unbalanced loads, such as chemistry
burning in cells by an implicit solver. Tiling can also create extra
work in the ghost cells of tiles.

Practical Details in Working with Tiling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With tiling, the OpenMP is now all in C, and not in Fortran for all
modules except reactions and initdata.

It is the responsibility of the coder to make sure that the routines
within a tiled region are safe to use with OpenMP. In particular,
note that:

-  tile boxes are non-overlapping

-  the union of tile boxes completely cover the valid region of the
   fab

-  Consider working with a node-centered MultiFab, ugdnv, and
   a cell-centered MultiFab, s:

   -  with mfi(s), the tiles are based on the cell-centered
      index space. If you have an :math:`8\times 8` box, then and 4 tiles,
      then your tiling boxes will range from :math:`0\rightarrow 3`,
      :math:`4\rightarrow 7`.

   -  with mfiugdnv, the tiles are based on nodal indices,
      so your tiling boxes will range from :math:`0\rightarrow 3`,
      :math:`4\rightarrow 8`.

-  When updating routines to work with tiling, we need to
   understand the distinction between the index-space of the entire box
   (which corresponds to the memory layout) and the index-space of the
   tile.

   -  In the C end, we pass (sometimes via the
      BL_TO_FORTRAN() macro) the loVect and hiVect of the
      entire box (including ghost cells). These are then used to
      allocate the array in Fortran as:

      ::

            double precision :: a(a_l1:a_h1, a_l2:a_h2, ...)

      When tiling is used, we do not want to loop as do a_l1,
      a_h1, but instead we need to loop over the tiling region. The
      indices of the tiling region need to be passed into the Fortran
      routine separately, and they come from the mfi.tilebox()
      or mfi.growntilebox() statement.

   -  In Fortran, when initializing an array to 0, do so only
      over the tile region, not for the entire box. For a Fortran array
      a, this means we cannot do:

      ::

            a = 0.0
            a(:,:,:,:) = 0.0

      but instead must do:

      ::

            a(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.0

      where lo() and hi() are the index-space for the tile box
      returned from mfi.tilebox() in C and passed into the Fortran
      routine.

   -  Look at r_old_s in Exec/gravity_tests/DustCollapse/probdata.f90 as an
      example of how to declare a threadprivate variable—this is then used
      in sponge_nd.f90.

Boundaries: FillPatch and FillPatchIterator
-------------------------------------------

AMReX calls the act of filling ghost cells a *fillpatch*
operation. Boundaries between grids are of two types. The first we
call “fine-fine”, which is two grids at the same level. The second
type is "coarse-fine", which needs interpolation from the coarse grid
to fill the fine grid ghost cells. Both of these are part of the
fillpatch operation. Fine-fine fills are just a straight copy from
“valid regions” to ghost cells. Coarse-fine fills are enabled
because the StateData is not just arrays, they’re “State Data”,
which means that the data knows how to interpolate itself (in an
anthropomorphical sense). The type of interpolation to use is defined
in Castro_setup.cpp—search for
cell_cons_interp, for example—that’s “cell conservative
interpolation”, i.e., the data is cell-based (as opposed to
node-based or edge-based) and the interpolation is such that the
average of the fine values created is equal to the coarse value from
which they came. (This wouldn’t be the case with straight linear
interpolation, for example.)

Additionally, since StateData has an old and new timelevel,
the fill patch operation can interpolate to an intermediate time.

Examples
~~~~~~~~

To illustrate the various ways we fill ghost cells and use the data,
let’s consider the following scenarios:

-  *You have state data that was defined with no ghost cells. You
   want to create a new MultiFab containing a copy of that data with
   NGROW ghost cells.*

   This is the case with —the MultiFab of the
   hydrodynamic state that we use to kick-off the hydrodynamics
   advance.

   Sborder is declared in Castro.H simply as:

   .. code:: c++

         Multifab Sborder;

   It is then allocated in

   .. code:: c++

         Sborder.define(grids, NUM_STATE, NUM_GROW, Fab_allocate);                   
         const Real prev_time = state[State_Type].prevTime();                        
         expand_state(Sborder, prev_time, NUM_GROW);      

   Note in the call to .define(), we tell AMReX to already
   allocate the data regions for the FArrayBoxs that are part of
   Sborder.

   The actually filling of the ghost cells is done by
   :

   .. code:: c++

         AmrLevel::FillPatch(*this, Sborder, NUM_GROW, 
                             prev_time, State_Type, 0, NUM_STATE);                

   Here, we are filling the ng ghost cells of MultiFab
   Sborder at time prev_time. We are using the
   StateData that is part of the current Castro object that we
   are part of. Note: FillPatch takes an object reference as its
   first argument, which is the object that contains the relevant
   StateData—that is what the this pointer indicates.
   Finally, we are copying the State_Type data components 0 to
   NUM_STATE [6]_.

   The result of this operation is that Sborder will now have
   NUM_GROW ghost cells consistent with the State_Type
   data at the old time-level.

-  *You have state data that was defined with NGROW ghost
   cells. You want to ensure that the ghost cells are filled
   (including any physical boundaries) with valid data.*

   This is very similar to the procedure shown above. The main
   difference is that for the MultiFab that will be the target
   of the ghost cell filling, we pass in a reference to the StateData itself.

   The main thing you need to be careful of here, is that you
   need to ensure that the the time you are at is consistent with
   the StateData’s time. Here’s an example from the radiation
   portion of the code MGFLDRadSolver.cpp:

   .. code:: c++

         Real time = castro->get_state_data(Rad_Type).curTime();
         MultiFab& S_new = castro->get_new_data(State_Type);

         AmrLevel::FillPatch(*castro, S_new, ngrow, time, State_Type,
                             0, S_new.nComp(), 0); 

   In this example, S_new is a pointer to the new-time-level
   State_Type MultiFab. So this operation will use the
   State_Type data to fill its own ghost cells. we fill the
   ngrow ghost cells of the new-time-level State_Type data,
   for all the components.

   Note that in this example, because the StateData lives in the
   Castro object and we are working from the Radiation object,
   we need to make reference to the current castro object
   pointer. If this were all done within the Castro object, then
   the pointer will simply be this, as we saw above.

-  *You have a MultiFab with some derived quantity. You want to
   fill its ghost cells.*

   MultiFabs have a FillBoundary() method that will fill all
   the ghost cells between boxes at the same level. It will not fill
   ghost cells at coarse-fine boundaries or at physical boundaries.

-  *You want to loop over the FABs in state data, filling ghost cells
   along the way*

   This is the job of the —this iterator is used to
   loop over the grids and fill ghostcells. A key thing to keep in
   mind about the FillPatchIterator is that you operate on a copy
   of the data—the data is disconnected from the original source. If
   you want to update the data in the source, you need to explicitly
   copy it back. Also note: FillPatchIterator takes a multifab,
   but this is not filled—this is only used to get the grid
   layout. Finally, the way the FillPatchIterator is implemented
   is that all the communication is done first, and then the iterating
   over boxes commences.

   For example, the loop that calls CA_UMDRV (all the
   hydrodynamics integration stuff) starts with

   ::

          for (FillPatchIterator fpi(*this, S_new, NUM_GROW,
                                     time, State_Type, strtComp, NUM_STATE);
                fpi.isValid(); ++fpi)
          {
            FArrayBox &state = fpi();
            Box bx(fpi.validbox());

            // work on the state FAB.  The interior (valid) cells will 
            // live between bx.loVect() and bx.hiVect()
          }

   Here the FillPatchIterator is the thing that distributes the
   grids over processors and makes parallel “just work”. This fills the
   single patch “fpi” , which has NUM_GROW ghost cells,
   with data of type “State_Type” at time “time”,
   starting with component strtComp and including a total of
   NUM_STATE components.

In general, one should never assume that ghostcells are valid, and
instead do a fill patch operation when in doubt. Sometimes we will
use a FillPatchIterator to fill the ghost cells into a multifab
without an explict look. This is done as:

::

      FillPatchIterator fpi(*this,S_old,1,time,State_Type,0,NUM_STATE);
      MultiFab& state_old = fpi.get_mf();     

In this operation, state_old points to the internal
MultiFab in the FillPatchIterator, by getting a reference to it as
fpi.get_mf(). This avoids a local copy.

Note that in the examples above, we see that only StateData can fill
physical boundaries (because these register how to fill the boundaries
when they are defined). There are some advanced operations in
AMReX that can get around this, but we do not use them in Castro.

.. _soft:phys_bcs:

Physical Boundaries
~~~~~~~~~~~~~~~~~~~

Physical boundary conditions are specified by an integer
index [7]_ in
the inputs file, using the and
runtime parameters. The generally
supported boundary conditions are, their corresponding integer key,
and the action they take for the normal velocity, transverse
velocity, and generic scalar are shown in Table \ `[table:castro:bcs] <#table:castro:bcs>`__

The definition of the specific actions are:

-  INT_DIR: data taken from other grids or interpolated

-  EXT_DIR: data specified on EDGE (FACE) of bndry

-  HOEXTRAP: higher order extrapolation to EDGE of bndry

-  FOEXTRAP: first order extrapolation from last cell in interior

-  REFLECT_EVEN: :math:`F(-n) = F(n)` true reflection from interior cells

-  REFLECT_ODD: :math:`F(-n) = -F(n)` true reflection from interior cells

The actual registration of a boundary condition action to a particular
variable is done in Castro_setup.cpp. At the top we define
arrays such as “scalar_bc”, “norm_vel_bc”, etc,
which say which kind of bc to use on which kind of physical boundary.
Boundary conditions are set in functions like “
set_scalar_bc”, which uses the scalar_bc pre-defined
arrays. We also specify the name of the Fortran routine that
is responsible for filling the data there (e.g., ).
These routines are discussed more below.

If you want to specify a value at a function (like at an inflow
boundary), then you choose an *inflow* boundary at that face of
the domain. You then need to write the implementation code for this.
An example is the problem which implements a
hydrostatic lower boundary (through its custom
routines.

.. raw:: latex

   \centering

.. table:: [table:castro:bcs] Physical boundary conditions supported in Castro. why does slipwall and noslipwall do the same thing?

   +-------------+-------------+-------------+-------------+-------------+
   | **name**    | **integer** | **normal    | **transvers | **scalars** |
   |             |             | velocity**  | e           |             |
   |             |             |             | velocity**  |             |
   +=============+=============+=============+=============+=============+
   | interior    | 0           | INT_DIR     | INT_DIR     | INT_DIR     |
   +-------------+-------------+-------------+-------------+-------------+
   | inflow      | 1           | EXT_DIR     | EXT_DIR     | EXT_DIR     |
   +-------------+-------------+-------------+-------------+-------------+
   | outflow     | 2           | FOEXTRAP    | FOEXTRAP    | FOEXTRAP    |
   +-------------+-------------+-------------+-------------+-------------+
   | symmetry    | 3           | REFLECT_ODD | REFLECT_EVE | REFLECT_EVE |
   |             |             |             | N           | N           |
   +-------------+-------------+-------------+-------------+-------------+
   | slipwall    | 4           | REFLECT_ODD | REFLECT_EVE | REFLECT_EVE |
   |             |             |             | N           | N           |
   +-------------+-------------+-------------+-------------+-------------+
   | noslipwall  | 5           | REFLECT_ODD | REFLECT_EVE | REFLECT_EVE |
   |             |             |             | N           | N           |
   +-------------+-------------+-------------+-------------+-------------+

FluxRegister
~~~~~~~~~~~~

A FluxRegister holds face-centered data at the boundaries of a box.
It is composed of a set of MultiFabs (one for each face, so 6 for
3D). A FluxRegister stores fluxes at coarse-fine interfaces,
and isused for the flux-correction step.

Other AMReX Concepts
--------------------

There are a large number of classes that help define the structure of
the grids, metadata associate with the variables, etc. A good way to
get a sense of these is to look at the .H files in the
amrex/Src/ directory.

Geometry class
~~~~~~~~~~~~~~

There is a Geometry object, for each level as part of
the Castro object (this is inhereted through AmrLevel).

ParmParse class
~~~~~~~~~~~~~~~

Error Estimators
~~~~~~~~~~~~~~~~

Gravity class
-------------

There is a single Gravity object, gravity, that is a
static class member of the Castro object. This means that all
levels refer to the same Gravity object.

Within the Gravity object, there are pointers to the Amr
object (as parent), and all of the AmrLevels (as a PArray,
LevelData). The gravity object gets the geometry
information at each level through the parent Amr class.

The main job of the gravity object is to provide the potential
and gravitation acceleration for use in the hydrodynamic sources.
Depending on the approximation used for gravity, this could mean
calling the AMReX multigrid solvers to solve the Poisson equation.

Fortran Helper Modules
----------------------

There are a number of modules that make data available to the Fortran
side of Castro or perform other useful tasks.

-  :

   This provides double precision constants as Fortran parameters, like
   ZERO, HALF, and ONE.

-  :

   This provides a double precision type, dp_t for use in
   Fortran. This should be identical to double precision on most
   architectures.

-  :

   This module provides access to the runtime parameters for the
   microphysics routines (EOS, reaction network, etc.). The source
   for this module is generated at compile type via a make rule
   that invokes a python script. This will search for all of the
   files in the external sources, parse them
   for runtime parameters, and build the module.

-  fundamental_constants_module:

   This provides the CGS values of many physical constants.

-  math_module:

   This provides simple mathematical functions. At the moment, a cross
   product routine.

-  meth_params_module:

   This module provides the integer keys used to access the state
   arrays for both the conserved variables (URHO, UMX, :math:`\ldots`)
   and primitive variables (QRHO, QU, :math:`\ldots`), as well
   as the number of scalar variables.

   It also provides the values of most of the castro.\ *xxxx*
   runtime parameters.

-  model_parser_module:

   This module is built if USE_MODELPARSER = TRUE is set in the
   problem’s GNUmakefile. It then provides storage for the an
   initial model and routines to read it in and interpolate onto the
   Castro grid.

-  prob_params_module:

   [soft:prob_params]

   This module stores information about the domain and current level,
   and is periodically synced up with the C driver. The information
   available here is:

   -  , : these are the boundary
      condition types at the low and high ends of the domain, for each
      coordinate direction. Integer keys, Interior, Inflow,
      Outflow, Symmetry, SlipWall, and
      NoSlipWall allow you to interpret the values.

   -  is the center of the problem. Note—this is up
      to the problem setup to define (in the probinit subroutine).
      Alternately, it can be set at runtime via
      .

      Usually center will be the physical center of the domain,
      but not always. For instance, for axisymmetric problems,
      center may be on the symmetry axis.

      center is used in the multipole gravity, hybrid advection
      algorithm, rotation sources, for the point mass gravity, in
      defining the center of the sponge, and in deriving the radial
      velocity.

   -  
   -  
   -  
   -  *refining information*

Setting Up Your Own Problem
---------------------------

To define a new problem, we create a new directory in one
of the subdirectories of Exec/,
and place in it a Prob_2d.f90 file (or 1d/3d,
depending on the dimensionality of the problem), a probdata.f90
file, the inputs and probin files, and a
Make.package file that tells the build system what problem-specific
routines exist. Finally, if you need custom boundary conditions, a
bc_fill_2d.F90 (or 1d/3d) file is needed. The
simplest way to get started is to copy these files from an existing
problem. Here we describe how to customize your problem.

The purpose of these files is:

-  : this holds the probdata_module Fortran module
   that allocates storage for all the problem-specific runtime parameters that
   are used by the problem (including those that are read from the probin
   file.

-  : this holds the main routines to
   initialize the problem and grid and perform problem-specific boundary
   conditions:

   -  probinit():

      This routine is primarily responsible for reading in the
      probin file (by defining the &fortin namelist and
      reading in an initial model (usually through the
      model_parser_module—see the toy_convect problem
      setup for an example). The parameters that are initialized
      here are those stored in the probdata_module.

   -  :

      This routine will initialize the state data for a single grid.
      The inputs to this routine are:

      -  level: the level of refinement of the grid we are filling

      -  time: the simulation time

      -  lo(), hi(): the integer indices of the box’s
         *valid data region* lower left and upper right corners. These
         integers refer to a global index space for the level and
         identify where in the computational domain the box lives.

      -  nscal: the number of scalar quantities—this is not typically
         used in Castro.

      -  state_l1, state_l2, (state_l3): the
         integer indices of the lower left corner of the box in each
         coordinate direction. These are for the box as allocated in memory,
         so they include any ghost cells as well as the valid data regions.

      -  state_h1, state_h2, (state_h3): the
         integer indices of the upper right corner of the box in each
         coordinate direction. These are for the box as allocated in memory,
         so they include any ghost cells as well as the valid data regions.

      -  state(): the main state array. This is dimensioned as:

         ::

             double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)

         (in 2-d), where NVAR comes from the meth_params_module.

         When accessing this array, we use the index keys provided by
         meth_params_module (e.g., URHO) to refer to specific
         quantities

      -  delta(): this is an array containing the zone width (:math:`\Delta x`)
         in each coordinate direction: :math:`\mathtt{delta(1)} = \Delta x`,
         :math:`\mathtt{delta(2)} = \Delta y`, :math:`\ldots`.

      -  xlo(), xhi(): these are the physical coordinates of the
         lower left and upper right corners of the *valid region*
         of the box. These can be used to compute the coordinates of the
         cell-centers of a zone as:

         ::

               do j = lo(2), hi(2)
                  y = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)
                  ...

         (Note: this method works fine for the problem initialization
         stuff, but for routines that implement tiling, as discussed below,
         lo and xlo may not refer to the same corner, and instead
         coordinates should be computed using problo() from the
         prob_params_module.)

-  :

   These routines handle how Castro fills ghostcells
   *at physical boundaries* for specific data. Most problem
   setups won’t need to do anything special here, and inclusion
   of this file is optional – only use it if you need to set
   specific boundary conditions.

   These routines are registered in Castro_setup.cpp, and
   called as needed. By default, they just
   pass the arguments through to filcc, which handles all of
   the generic boundary conditions (like reflecting, extrapolation,
   etc.). The specific ‘fill’ routines can then supply the
   problem-specific boundary conditions, which are typically just
   Dirichlet boundary conditions (usually this means looking to see
   if the bc() flag at a boundary is EXT_DIR. The
   problem-specific code implementing these specific conditions
   should *follow* the filcc call.

   -  ca_hypfill:
      This handles the boundary filling for the hyperbolic system.

   -  ca_denfill: At times, we need to fill just the density
      (always assumed to be the first element in the hyperbolic state)
      instead of the entire state. When the fill patch routine is called
      with first_comp = Density and num_comp = 1, then we
      use ca_denfill instead of ca_hypfill.

      (Note: it seems that this may be used for more than just
      density, but it is only used for tagging and the plotfile)

   -  ca_grav?fill: These routines fill will the ghostcells
      of the gravitational acceleration grids with the gravitational
      acceleration.

      Note: for constant gravity, these routines will never be called.
      For one of the Poisson-type gravities, you only need to do
      something special here if you are implementing an Interior
      boundary type (which you can test for by comparing
      bc(:,:,:) to EXT_DIR.

      For the other standard physical boundary types, the ghost cell
      filling will be handled automatically by the default filcc
      call in these routines.

      The gravitational acceleration in the ghost cells is used during
      the hydrodynamics portion of the code in predicting the
      interface states.

   -  ca_reactfill: This handles boundary filling for
      any Reactions_Type MultiFABs, which are sometimes used to interface
      with the nuclear burning module. It stores the normal state data
      in addition to components for the energy release and species change.

   These routines take the following arguments:

   -  adv_l1, adv_l2, (adv_l3): the indicies of
      the lower left corner of the box holding the data we are working on.
      These indices refer to the entire box, including ghost cells.

   -  adv_h1, adv_h2, (adv_h3): the indicies of
      the upper right corner of the box holding the data we are working on.
      These indices refer to the entire box, including ghost cells.

   -  adv(): the array of data whose ghost cells we are filling.
      Depending on the routine, this may have an additional index refering
      to the variable.

      This is dimensioned as:

      ::

            double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

   -  domlo(), domhi(): the integer indices of the lower
      left and upper right corners of the valid region of the *entire
      domain*. These are used to test against to see if we are filling
      physical boundary ghost cells.

      This changes according to refinement level: level-0 will
      range from 0 to castro.max_grid_size,
      and level-n will range from 0 to
      :math:`\mathtt{castro.max\_grid\_size} \cdot \prod_n \mathtt{castro.ref\_ratio(n)}`.

   -  delta(): is the zone width in each coordinate direction,
      as in initdata() above.

   -  xlo(): this is the physical coordinate of the lower
      left corner of the box we are filling—including the ghost cells.

      Note: this is different than how xlo() was defined in
      initdata() above.

   -  time: the simulation time

   -  bc(): an array that holds the type of boundary conditions
      to enforce at the physical boundaries for adv.

      Sometimes it appears of the form bc(:,:) and sometimes
      bc(:,:,:)—the last index of the latter holds the variable
      index, i.e., density, pressure, species, etc.

      The first index is the coordinate direction and the second index
      is the domain face (1 is low, 2 is hi), so
      bc(1,1) is the lower :math:`x` boundary type, bc(1,2) is
      the upper :math:`x` boundary type, bc(2,1) is the lower
      :math:`y` boundary type, etc.

      To interpret the array values, we test against the quantities
      defined in bc_types.fi included in each subroutine,
      for example, EXT_DIR, FOEXTRAP, :math:`\ldots`. The
      meaning of these are explained below.

Optional Files
~~~~~~~~~~~~~~

The follow problem-specific files are optional. There are stubs for
each of these in the main source tree.

-  :

   This provides two routines, and
   that can be used to add information to the
   checkpoint files and read it in upon restart. This is useful for
   some global problem-specific quantities. For instance, the
    [8]_ problem uses this
   to store center of mass position and velocity information in the
   checkpoint files that are used for runtime diagnostics.

   The name of the checkpoint directory is passed in as an argument.
   provides the C interfaces for these routines.

-  ,

   This implements problem-specific tagging for refinement, through a
   subroutine . The full hydrodynamic state
   (State_Type) is passed in, and the problem can mark zones for
   refinement by setting the variable for a zone to
   . An example is provided by the
   problem which refines a rectangular region (fuel layer) based on
   a density parameter and the H mass fraction.

-  , ,

   Together, these provide a mechanism to create derived quantities
   that can be stored in the plotfile. Problem_Derives.H
   provides the C code that defines these new plot variables. It
   does this by adding them to the —a list of
   derived variables that Castro knows about. When adding new
   variables, a descriptive name, Fortran routine that does the
   deriving, and component of StateData are specified.

   The Fortran routine that does the deriving is put in the
   problem-specific problem_derive_nd.f90 (and a prototype for
   C is put in Problem_Derives.H). A example is provided by
   the problem, which derives several new
   quantities (perturbations against a background one-dimensional
   model, in this case).

-  , ,

   These files provide problem-specific routines for computing global
   diagnostic information through the
   functionality that is part of the Castro class.

   An example is provided by , where an estimate
   of the flame speed is computed by integrating the mass of fuel on
   the grid.

Dimension Agnostic Problem Initialization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most of the problem setups have separate implementations for 1-, 2-,
and 3D. A new method exists that allows you to write just a single
set of files for any dimensionality (this is called the *dimension
agnostic* format). To use this mode, set
in your GNUmakefile.
Then write you problem initialization in .
Analogous routines exist for tagging and boundary conditions. See the
and problem setups for an
example.

.. _software:io:

Parallel I/O
------------

Both checkpoint files and plotfiles are really directories containing
subdirectories: one subdirectory for each level of the AMR hierarchy.
The fundamental data structure we read/write to disk is a MultiFab,
which is made up of multiple FAB’s, one FAB per grid. Multiple
MultiFabs may be written to each directory in a checkpoint file.
MultiFabs of course are shared across CPUs; a single MultiFab may be
shared across thousands of CPUs. Each CPU writes the part of the
MultiFab that it owns to disk, but they don’t each write to their own
distinct file. Instead each MultiFab is written to a runtime
configurable number of files N (N can be set in the inputs file as the
parameter and ; the
default is 64). That is to say, each MultiFab is written to disk
across at most N files, plus a small amount of data that gets written
to a header file describing how the file is laid out in those N files.

What happens is :math:`N` CPUs each opens a unique one of the :math:`N` files into
which the MultiFab is being written, seeks to the end, and writes
their data. The other CPUs are waiting at a barrier for those :math:`N`
writing CPUs to finish. This repeats for another :math:`N` CPUs until all the
data in the MultiFab is written to disk. All CPUs then pass some data
to CPU 0 which writes a header file describing how the MultiFab is
laid out on disk.

We also read MultiFabs from disk in a “chunky” manner, opening only :math:`N`
files for reading at a time. The number :math:`N`, when the MultiFabs were
written, does not have to match the number :math:`N` when the MultiFabs are
being read from disk. Nor does the number of CPUs running while
reading in the MultiFab need to match the number of CPUs running when
the MultiFab was written to disk.

Think of the number :math:`N` as the number of independent I/O pathways in
your underlying parallel filesystem. Of course a “real” parallel
filesytem should be able to handle any reasonable value of :math:`N`. The
value -1 forces :math:`N` to the number of CPUs on which you’re
running, which means that each CPU writes to a unique file, which can
create a very large number of files, which can lead to inode issues.

Single-Level Flow Chart
=======================

.. _introduction-1:

Introduction
------------

There are several different time-evolution methods currently
implemented in Castro. As best as possible, they share the same
driver routines and use preprocessor or runtime variables to separate
the different code paths.

-  Strang-splitting: the Strang evolution does the burning on the
   state for :math:`\Delta t/2`, then updates the hydrodynamics using the
   burned state, and then does the final :math:`\Delta t/2` burning. No
   explicit coupling of the burning and hydro is done. Within the
   Strang code path, there are two methods for doing the hydrodynamics,
   controlled by .

   -  Corner-transport upwind (CTU): this implements the unsplit,
      characteristic tracing method of :raw-latex:`\cite{colella:1990}`.

   -  Method of lines (MOL): this discretizes the space part of
      our system without any characteristic tracing and uses an
      ODE integrator to advance the state. Multiple stages can be done,
      each requiring reconstruction, Riemann solve, etc., and the final
      solution is pieced together from the intermediate stages.

-  SDC: the SDC path is enabled by the preprocessor
   variable. This iteratively couples the reactions and hydrodynamics together.

Several helper functions are used throughout:

-  :
   There are many ways that the hydrodynamics state may become
   unphysical in the evolution. The routine
   enforces some checks on the state. In particular, it

   #. enforces that the density is above

   #. normalizes the species so that the mass fractions sum to 1

   #. resets the internal energy if necessary (too small or negative)
      and computes the temperature for all zones to be thermodynamically
      consistent with the state.

.. _flow:sec:nosdc:

Overview of a single step (no SDC)
----------------------------------

The main evolution for a single step is contained in
, as . This does
the following advancement. Note, some parts of this are only done
depending on which preprocessor directives are defined at
compile-time—the relevant directive is noted in the [ ] at the start
of each step.

#. *Initialization* ()

   This sets up the current level for advancement. The following
   actions are performend (note, we omit the actions taken for a retry,
   which we will describe later):

   -  Sync up the level information to the Fortran-side of Castro

   -  Do any radiation initialization

   -  Initialize all of the intermediate storage arrays (like those
      that hold source terms, etc.).

   -  Swap the StateData from the new to old (e.g., ensures that
      the next evolution starts with the result from the previous step).

   -  Do a

   -  Create the MultiFabs that hold the primitive variable information
      for the hydro solve.

   -  For method of lines integration: allocate the storage for the
      intermediate stage updates, , and the  that holds the post burn state.

   -  Zero out all of the fluxes

#. *Advancement*

   The update strategy differs for CTU vs MOL:

   -  CTU: Calls to take a single step,
      incorporating hydrodynamics, reactions, and source terms.

   -  MOL: Call times
      (i.e., once for each of the intermediate stages in the ODE
      integration). Within do_advance we will use the stage
      number, , to do an pre- or post-hydro
      sources (e.g., burning).

   In either case, for radiation-hydrodynamics, this step does the
   advective (hyperbolic) portion of the radiation update only.
   Source terms, including gravity, rotation, and diffusion are
   included in this step, and are time-centered to achieve second-order
   accuracy.

   If is set, then we subcycle the current
   step if we violated any stability criteria to reach the desired
   :math:`\Delta t`. The idea is the following: if the timestep that you
   took had a timestep that was not sufficient to enforce the stability
   criteria that you would like to achieve, such as the CFL criterion
   for hydrodynamics or the burning stability criterion for reactions,
   you can retry the timestep by setting castro.use_retry = 1 in
   your inputs file. This will save the current state data at the
   beginning of the level advance, and then if the criteria are not
   satisfied, will reject that advance and start over from the old
   data, with a series of subcycled timesteps that should be small
   enough to satisfy the criteria. Note that this will effectively
   double the memory footprint on each level if you choose to use it.

#. [] *Auxiliary quantitiy evolution*

   Auxiliary variables in Castro are those that obey a continuity
   equation (with optional sources) that are passed into the EOS, but
   not subjected to the constraint on mass fractions (summing to one).

   The advection and source terms are already dealt with in the
   main hydrodynamics advance (above step). A user-supplied routine
   can be provided here to further update these
   quantities.

#. *Radial data and [] point mass*

   If is set, then we average the state data
   over angles here to create a radial profile. This is then used in the
   boundary filling routines to properly set Dirichlet BCs when our domain
   is smaller than the star, so the profile on the boundaries will not
   be uniform.

   If is set, then we
   change the mass of the point mass that optionally contributes to the
   gravitational potential by taking mass from the surrounding zones
   (keeping the density in those zones constant).

#. [] *Radiation implicit update*

   The do_advance() routine only handled the hyperbolic
   portion of the radiation update. This step does the implicit solve
   (either gray or multigroup) to advance the radiation energies to the
   new time level. Note that at the moment, this is backward-difference
   implicit (first-order in time) for stability.

   This is handled by .

#. [] *Particles*

   If we are including passively-advected particles, they are
   advanced in this step.

#. *Finalize*

   This cleans up at the end of a step:

   -  Update the flux registers to account for mismatches at
      coarse-fine interfaces. This cleans up the memory used during
      the step.

   -  If is set, then we
      also add up the mass that left through the boundary over this
      step. [9]_

   -  Free any memory allocated for the level advance.

Main Hydro, Reaction, and Gravity Advancement (CTU w/ Strang-splitting)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The explicit portion of the system advancement (reactions,
hydrodynamics, and gravity) is done by . Consider
our system of equations as:

.. math:: \frac{\partial{\bf U}}{\partial t} = -{\bf A}({\bf U}) + {\bf R}({\bf U}) + {\bf S},

where :math:`{\bf A}({\bf U}) = \nabla \cdot {\bf F}({\bf U})`, with :math:`{\bf F}` the flux vector, :math:`{\bf R}` are the reaction
source terms, and :math:`{\bf S}` are the non-reaction source terms, which
includes any user-defined external sources, :math:`{\bf S}_{\rm ext}`. We use
Strang splitting to discretize the advection-reaction equations. In
summary, for each time step, we update the conservative variables,
:math:`{\bf U}`, by reacting for half a time step, advecting for a full time
step (ignoring the reaction terms), and reacting for half a time step.
The treatment of source terms complicates this a little. The actual
update, in sequence, looks like:

.. math::

   \begin{aligned}
   {\bf U}^\star &= {\bf U}^n + \frac{\Delta t}{2}{\bf R}({\bf U}^n) \\
   {\bf U}^{n+1,(a)} &= {\bf U}^\star + \Delta t\, {\bf S}({\bf U}^\star) \\
   {\bf U}^{n+1,(b)} &= {\bf U}^{n+1,(a)} - \Delta t\, {\bf A}({\bf U}^\star) \\
   {\bf U}^{n+1,(c)} &= {\bf U}^{n+1,(b)} + \frac{\Delta t}{2}\, [{\bf S}({\bf U}^{n+1,(b)}) - {\bf S}({\bf U}^\star)] \label{eq:source_correct}\\
   {\bf U}^{n+1}     &= {\bf U}^{n+1,(c)} + \frac{\Delta t}{2} {\bf R}({\bf U}^{n+1,(c)})\end{aligned}

Note that in the first step, we add a full :math:`\Delta t` of the old-time
source to the state. This prediction ensures consistency when it
comes time to predicting the new-time source at the end of the update.
The construction of the advective terms, :math:`{\bf A({\bf U})}` is purely
explicit, and based on an unsplit second-order Godunov method. We
predict the standard primitive variables, as well as :math:`\rho e`, at
time-centered edges and use an approximate Riemann solver construct
fluxes.

At the beginning of the time step, we assume that :math:`{\bf U}` and :math:`\phi` are
defined consistently, i.e., :math:`\rho^n` and :math:`\phi^n` satisfy equation
(`[eq:Self Gravity] <#eq:Self Gravity>`__). Note that in
Eq. \ `[eq:source_correct] <#eq:source_correct>`__, we actually can actually do some
sources implicitly by updating density first, and then momentum,
and then energy. This is done for rotating and gravity, and can
make the update more akin to:

.. math:: {\bf U}^{n+1,(c)} = {\bf U}^{n+1,(b)} + \frac{\Delta t}{2} [{\bf S}({\bf U}^{n+1,(c)}) - {\bf S}({\bf U}^n)]

Castro also supports radiation. This part of the update algorithm
only deals with the advective / hyperbolic terms in the radiation update.

Here is the single-level algorithm. The goal here is to update the
 from the old to new time (see
§ \ `2.6 <#soft:sec:statedata>`__). We will use the following notation
here, consistent with the names used in the code:

-  is a multifab reference to the old-time-level
   State_Type data.

-  is a multifab that has ghost cells and is
   initialized from S_old. This is what the hydrodynamic
   reconstruction will work from.

-  is a multifab reference to the new-time-level
   State_Type data.

In the code, the objective is to evolve the state from the old time,
S_old, to the new time, S_new.

#. [strang:init] *Initialize*

   #. In :

      #. Create , initialized from S_old

   #. Check for NaNs in the initial state, S_old.

#. *React :math:`\Delta t/2`.* []

   Update the solution due to the effect of reactions over half a time
   step. The integration method and system of equations used here is
   determined by a host of runtime parameters that are part of the
   Microphysics package. But the basic idea is to evolve the energy
   release from the reactions, the species mass fractions, and
   temperature through :math:`\Delta t/2`.

   Using the notation above, we begin with the time-level :math:`n` state,
   :math:`{\bf U}^n`, and produce a state that has evolved only due to reactions,
   :math:`{\bf U}^\star`.

   .. math::

      \begin{aligned}
          (\rho e)^\star &= (\rho e)^\star - \frac{\Delta t}{2} \rho H_\mathrm{nuc} \\
          (\rho E)^\star &= (\rho E)^\star - \frac{\Delta t}{2} \rho H_\mathrm{nuc} \\
          (\rho X_k)^\star &= (\rho X_k)^\star + \frac{\Delta t}{2}(\rho\dot\omega_k)^n.
        \end{aligned}

   Here, :math:`H_\mathrm{nuc}` is the energy release (erg/g/s) over the
   burn, and :math:`\dot\omega_k` is the creation rate for species :math:`k`.

   After exiting the burner, we call the EOS with :math:`\rho^\star`,
   :math:`e^\star`, and :math:`X_k^\star` to get the new temperature, :math:`T^\star`.

   Note that the density, :math:`\rho`, does not change via reactions in the
   Strang-split formulation.

   The reaction data needs to be valid in the ghost cells. The logic
   in this routine (accomplished throuh the use of a mask) will burn
   only in the valid interior cells or in any ghost cells that are on a
   coarse-fine interface or physical boundary. This allows us to just
   use a level FillBoundary() call to fill all of the ghost cells
   on the same level with valid data.

   An experimental option (enabled via
   ) will create a custom
   distribution map based on the work needed in burning a zone and
   redistribute the boxes across processors before burning, to better
   load balance..

   After reactions, is called.

   At the end of this step, sees the effects of the
   reactions.

#. [strang:oldsource] *Construct time-level :math:`n` sources and apply*
   [, ]

   The time level :math:`n` sources are computed, and added to the
   StateData . The sources are then applied
   to the state after the burn, :math:`{\bf U}^\star` with a full :math:`\Delta t`
   weighting (this will be corrected later). This produces the
   intermediate state, :math:`{\bf U}^{n+1,(a)}`.

   The sources that we deal with here are:

   #. sponge : the sponge is a damping term added to
      the momentum equation that is designed to drive the velocities to
      zero over some timescale. Our implementation of the sponge
      follows that of Maestro :raw-latex:`\cite{maestro:III}`

   #. external sources : users can define problem-specific sources
      in the file. Sources for the different
      equations in the conservative state vector, :math:`{\bf U}`, are indexed
      using the integer keys defined in meth_params_module
      (e.g., URHO).

      This is most commonly used for external heat sources (see the
      problem setup) for an example. But most
      problems will not use this.

   #. [] diffusion : thermal diffusion can be
      added in an explicit formulation. Second-order accuracy is
      achieved by averaging the time-level :math:`n` and :math:`n+1` terms, using
      the same predictor-corrector strategy described here.

      Note: thermal diffusion is distinct from radiation hydrodynamics.

      Also note that incorporating diffusion brings in an additional
      timestep constraint, since the treatment is explicit. See
      Chapter \ `[ch:diffusion] <#ch:diffusion>`__ for more details.

   #. [] angular momentum

      .. raw:: latex

         \marginpar{\vskip-\baselineskip\raggedright\tiny\sffamily
         \hrule\smallskip{\color{red}need to write this up}\par\smallskip\hrule}

   #. [] gravity:

      For full Poisson gravity, we solve for for gravity using:

      .. math::

         {\bf g}^n = -\nabla\phi^n, \qquad
               \Delta\phi^n = 4\pi G\rho^n,

      The construction of the form of the gravity source for the
      momentum and energy equation is dependent on the parameter
      . Full details of the gravity
      solver are given in Chapter \ `[ch:gravity] <#ch:gravity>`__.

      .. raw:: latex

         \marginpar{\vskip-\baselineskip\raggedright\tiny\sffamily
         \hrule\smallskip{\color{red}we should add a description of whether we do a level solve or a composite solve}\par\smallskip\hrule}

      .. raw:: latex

         \marginpar{\vskip-\baselineskip\raggedright\tiny\sffamily
         \hrule\smallskip{\color{red}what do we store? phi and g? source?}\par\smallskip\hrule}

   #. [] rotation

      We compute the rotational potential (for use in the energy update)
      and the rotational acceleration (for use in the momentum
      equation). This includes the Coriolis and centrifugal terms in a
      constant-angular-velocity co-rotating frame. The form of the
      rotational source that is constructed then depends on the
      parameter . More details are
      given in Chapter \ `[ch:rotation] <#ch:rotation>`__.

   The source terms here are evaluated using the post-burn state,
   :math:`{\bf U}^\star` (), and later corrected by using the
   new state just before the burn, :math:`{\bf U}^{n+1,(b)}`. This is compatible
   with Strang-splitting, since the hydro and sources takes place
   completely inside of the surrounding burn operations.

   Note that the source terms are already applied to
   in this step, with a full :math:`\Delta t`—this will be corrected later.

#. [strang:hydro] *Construct the hydro update* []

   The goal is to advance our system considering only the advective
   terms (which in Cartesian coordinates can be written as the
   divergence of a flux).

   We do the hydro update in two parts—first we construct the
   advective update and store it in the , then we do the conservative update in a separate step. This
   separation allows us to use the advective update separately in more
   complex time-integration schemes.

   In the Strang-split formulation, we start the reconstruction using
   the state after burning, :math:`{\bf U}^\star` (). There
   are two approaches we use, the corner transport upwind (CTU) method
   that uses characteristic tracing as described in
   :raw-latex:`\cite{colella:1990}`, and a method-of-lines approach. The choice is
   determined by the parameter .

   #. CTU method:

      For the CTU method, we predict to the half-time (:math:`n+1/2`) to get a
      second-order accurate method. Note: does not
      know of any sources except for reactions. The advection step is
      complicated, and more detail is given in Section
      `6 <#Sec:Advection Step>`__. Here is the summarized version:

      #. Compute primitive variables.

      #. Convert the source terms to those acting on primitive variables

      #. Predict primitive variables to time-centered edges.

      #. Solve the Riemann problem.

      #. Compute fluxes and update.

      To start the hydrodynamics, we need to know the hydrodynamics source
      terms at time-level :math:`n`, since this enters into the prediction to
      the interface states. This is essentially the same vector that was
      computed in the previous step, with a few modifications. The most
      important is that if we set
      , then we extrapolate the
      source terms from :math:`n` to :math:`n+1/2`, using the change from the previous
      step.

      Note: we neglect the reaction source terms, since those are already
      accounted for in the state directly, due to the Strang-splitting
      nature of this method.

      The update computed here is then immediately applied to
      .

   #. method of lines

#. [strang:clean] *Clean State* []

   .. raw:: latex

      \marginpar{\vskip-\baselineskip\raggedright\tiny\sffamily
      \hrule\smallskip{\color{red}we only seem to do this for the MOL integration}\par\smallskip\hrule}

   This is done on .

   After these checks, we check the state for NaNs.

#. [strang:radial] *Update radial data and center of mass for monopole gravity*

   These quantities are computed using .

#. [strang:newsource] *Correct the source terms with the :math:`n+1` contribution*
   [, ]

   Previously we added :math:`\Delta t\, {\bf S}({\bf U}^\star)` to the state, when
   we really want a time-centered approach, :math:`(\Delta t/2)[{\bf S}({\bf U}^\star
       + {\bf S}({\bf U}^{n+1,(b)})]`. We fix that here.

   We start by computing the source term vector :math:`{\bf S}({\bf U}^{n+1,(b)})`
   using the updated state, :math:`{\bf U}^{n+1,(b)}`. We then compute the
   correction, :math:`(\Delta t/2)[{\bf S}({\bf U}^{n+1,(b)}) - {\bf S}({\bf U}^\star)]` to
   add to :math:`{\bf U}^{n+1,(b)}` to give us the properly time-centered source,
   and the fully updated state, :math:`{\bf U}^{n+1,(c)}`. This correction is stored
   in the  [10]_.

   In the process of updating the sources, we update the temperature to
   make it consistent with the new state.

#. *React :math:`\Delta t/2`.* []

   We do the final :math:`\Delta t/2` reacting on the state, begining with :math:`{\bf U}^{n+1,(c)}` to
   give us the final state on this level, :math:`{\bf U}^{n+1}`.

   This is largely the same as , but
   it does not currently fill the reactions in the ghost cells.

#. [strang:finalize] *Finalize* []

   Finalize does the following:

   #. for the momentum sources, we compute :math:`d{\bf S}/dt`, to use in the
      source term prediction/extrapolation for the hydrodynamic
      interface states during the next step.

   #. If we are doing the hybrid momentum algorithm, then we sync up
      the hybrid and linear momenta

A summary of which state is the input and which is updated for each of
these processes is presented below:

+--------------------+-------+-----------------+-----------------+
| *step*             | S_old | Sborder         | S_new           |
+====================+=======+=================+=================+
| 1. init            | input | updated         |                 |
+--------------------+-------+-----------------+-----------------+
| 2. react           |       | input / updated |                 |
+--------------------+-------+-----------------+-----------------+
| 3. old sources     |       | input           | updated         |
+--------------------+-------+-----------------+-----------------+
| 4. hydro           |       | input           | updated         |
+--------------------+-------+-----------------+-----------------+
| 5. clean           |       |                 | input / updated |
+--------------------+-------+-----------------+-----------------+
| 6. radial / center |       |                 | input           |
+--------------------+-------+-----------------+-----------------+
| 7. correct sources |       |                 | input / updated |
+--------------------+-------+-----------------+-----------------+
| 8. react           |       |                 | input / updated |
+--------------------+-------+-----------------+-----------------+

Main Hydro, Reaction, and Gravity Advancement (MOL w/ Strang-splitting)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The handling of sources differs in the MOL integration, as compared to CTU.
Again, consider our system as:

.. math:: \frac{\partial{\bf U}}{\partial t} = -{\bf A}({\bf U}) + {\bf R}({\bf U}) + {\bf S}\, .

We will again use Strang splitting to discretize the
advection-reaction equations, but the hydro update will consist of :math:`s`
stages. The update first does the reactions, as with CTU:

.. math:: {\bf U}^\star = {\bf U}^n + \frac{\Delta t}{2}{\bf R}({\bf U}^n)

We then consider the hydro update discretized in space, but not time, written
as:

.. math:: \frac{\partial {\bf U}}{\partial t} = -{\bf A}({\bf U}) + {\bf S}({\bf U})

Using a Runge-Kutta (or similar) integrator, we write the update as:

.. math:: {\bf U}^{n+1,\star} = {\bf U}^\star + \Delta t\sum_{l=1}^s b_i {\bf k}_l

where :math:`b_i` is the weight for stage :math:`i` and :math:`k_i` is the stage update:

.. math:: {\bf k}_l = -{\bf A}({\bf U}_l) + {\bf S}({\bf U}_l)

with

.. math:: {\bf U}_l = {\bf U}^\star  + \Delta t\sum_{m=1}^{l-1} a_{lm} {\bf k}_m

Finally, there is the last part of the reactions:

.. math:: {\bf U}^{n+1} = {\bf U}^{n+1,\star} + \frac{\Delta t}{2} {\bf R}({\bf U}^{n+1,\star})

In contrast to the CTU method, the sources are treated together
with the advection here.

The time at the intermediate stages is evaluated as:

.. math:: t_l = c_l \Delta t

The integration coefficients are stored in the vectors
, , and , and the
stage updates are stored in the MultiFab .

Here is the single-level algorithm. We use the same notation
as in the CTU flowchart.

In the code, the objective is to evolve the state from the old time,
S_old, to the new time, S_new.

#. [strang:init] *Initialize*

   In , set the starting point for the stage’s integration:

   #. if : initialize
      from

   #. if : we need to create
      the starting point for the current stage. We store this,
      temporarily in the new-time slot (what we normally refer to as
      ):

      .. math:: \mathtt{S\_new}_\mathrm{iter} = \mathtt{Sburn} + \Delta t\sum_{l=0}^{\mathrm{iter}-1} a_{\mathrm{iter},l} \mathtt{k\_mol}_l

      Then initialize from .

   Check for NaNs in the initial state, S_old.

#. *React :math:`\Delta t/2`.* []

   This step is unchanged from the CTU version. At the end of this
   step, sees the effects of the reactions.

   Each stage needs to build its starting point from this point, so we
   store the effect of the burn in a new MultiFab, ,
   for use in the stage initialization.

#. [strang:oldsource] *Construct sources from the current
   stage’s state*
   [, ]

   .. raw:: latex

      \marginpar{\vskip-\baselineskip\raggedright\tiny\sffamily
      \hrule\smallskip{\color{red}fix: gravity is still using{\tt S\_old}}\par\smallskip\hrule}

   The time level :math:`n` sources are computed, and added to the
   StateData . The sources are then applied
   to the state after the burn, :math:`{\bf U}^\star` with a full :math:`\Delta t`
   weighting (this will be corrected later). This produces the
   intermediate state, :math:`{\bf U}^{n+1,(a)}`.

   For full Poisson gravity, we solve for for gravity using:

   .. math::

      {\bf g}^n = -\nabla\phi^n, \qquad
          \Delta\phi^n = 4\pi G\rho^n,

#. [strang:hydro] *Construct the hydro update* []

   The hydro update in the MOL branch will include both the advective
   and source terms. In each stage, store in the righthand side for the current stage.

   In constructing the stage update, we use the source evaluated earlier,
   and compute:

   .. math:: \mathtt{k\_mol}_l = - {\bf A}({\bf U}_l) + {\bf S}({\bf U}_l)

   Each call to do_advance_mol only computes this update for
   a single stage. On the last stage, we compute the final update
   as:

   .. math:: \mathtt{S\_new} = \mathtt{Sburn} + \Delta t\sum_{l=0}^{\mathrm{n\_stages}-1} b_l \, \mathrm{k\_mol}_l

#. [strang:clean] *Clean State* []

   .. raw:: latex

      \marginpar{\vskip-\baselineskip\raggedright\tiny\sffamily
      \hrule\smallskip{\color{red}we only seem to do this for the MOL integration}\par\smallskip\hrule}

   This is done on .

   After these checks, we check the state for NaNs.

#. *React :math:`\Delta t/2`.* []

   We do the final :math:`\Delta t/2` reacting on the state, begining with :math:`{\bf U}^{n+1,(c)}` to
   give us the final state on this level, :math:`{\bf U}^{n+1}`.

   This is largely the same as , but
   it does not currently fill the reactions in the ghost cells.

#. [strang:finalize] *Finalize* []

   Finalize does the following:

   #. for the momentum sources, we compute :math:`d{\bf S}/dt`, to use in the
      source term prediction/extrapolation for the hydrodynamic
      interface states during the next step.

   #. If we are doing the hybrid momentum algorithm, then we sync up
      the hybrid and linear momenta

A summary of which state is the input and which is updated for each of
these processes is presented below:

+--------------------+-------+-----------------+-----------------+
| *step*             | S_old | Sborder         | S_new           |
+====================+=======+=================+=================+
| 1. init            | input | updated         |                 |
+--------------------+-------+-----------------+-----------------+
| 2. react           |       | input / updated |                 |
+--------------------+-------+-----------------+-----------------+
| 3. old sources     |       | input           | updated         |
+--------------------+-------+-----------------+-----------------+
| 4. hydro           |       | input           | updated         |
+--------------------+-------+-----------------+-----------------+
| 5. clean           |       |                 | input / updated |
+--------------------+-------+-----------------+-----------------+
| 6. radial / center |       |                 | input           |
+--------------------+-------+-----------------+-----------------+
| 7. correct sources |       |                 | input / updated |
+--------------------+-------+-----------------+-----------------+
| 8. react           |       |                 | input / updated |
+--------------------+-------+-----------------+-----------------+

Overview of a single step (with SDC)
------------------------------------

We express our system as:

.. math:: {\bf U}_t = \mathcal{A}({\bf U}) + {\bf R}({\bf U})

here :math:`\mathcal{A}` is the advective source, which includes both the
flux divergence and the hydrodynamic source terms (e.g. gravity):

.. math:: \mathcal{A}({\bf U}) = -\nabla \cdot {\bf F}({\bf U}) + {\bf S}

The SDC version of the main advance loop looks similar to the no-SDC
version, but includes an iteration loop over the hydro, gravity, and
reaction update. So the only difference happens in step 2 of the
flowchart outlined in § \ `2 <#flow:sec:nosdc>`__. In particular this
step now proceeds as:

2. *Advancement*

   Loop :math:`k` from 0 to sdc_iters, doing:

   #. *Hydrodynamics advance*: This is done through
      do_advance—in SDC mode, this only updates the hydrodynamics,
      including the non-reacting sources. However, in predicting the
      interface states, we use an iteratively-lagged approximation to the
      reaction source on the primitive variables, :math:`\mathcal{I}_q^{k-1}`.

      The result of this is an approximation to :math:`\mathcal{A}({\bf U})`,
      stored in (the flux divergence)
      and and .

   #. *React*: Reactions are integrated with the advective
      update as a source—this way the reactions see the
      time-evolution due to advection as we integrate:

      .. math:: \frac{d{\bf U}}{dt} = \left [ \mathcal{A}({\bf U}) \right ]^{n+1/2} + {\bf R}({\bf U})

      The advective source includes both the divergence of the fluxes
      as well as the time-centered source terms. This is computed by
      by summing over all source components
      , , and
      .

   #. *Clean state*: This ensures that the thermodynamic state is
      valid and consistent.

   #. *Construct reaction source terms*: Construct the change
      in the primitive variables due only to reactions over the
      timestep, :math:`\mathcal{I}_q^{k}`. This will be used in the next
      iteration.

Note that is it likely that some of the other updates (like any
non-advective auxiliary quantity updates) should be inside the SDC
loop, but presently they are only done at the end. Also note that the
radiation implicit update is not done as part of the SDC iterations.

Main Hydro and Gravity Advancement (SDC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The evolution in do_advance is substantially different than the
Strang case. In particular, reactions are not evolved. Here we
summarize those differences.

#. *Initialize* []

   This is unchanged from step `[strang:init] <#strang:init>`__ in the Strang algorithm.

#. *Construct time-level :math:`n` sources and apply*
   [, ]

   This corresponds to step `[strang:oldsource] <#strang:oldsource>`__ in the Strang
   algorithm. There are not differences compared to the Strang
   algorithm, although we note, this only needs to be done for the first
   SDC iteration in the advancement, since the old state does not change.

#. *Construct the hydro update* []

   This corresponds to step \ `[strang:hydro] <#strang:hydro>`__ in the Strang
   algorithm. There are a few major differences with the Strang case:

   -  There is no need to extrapolate source terms to the half-time
      for the prediction (the
      parameter), since SDC provides a natural way to approximate the
      time-centered source—we simply use the iteratively-lagged new-time
      source.

   -  The primitive variable source terms that are used for the
      prediction include the contribution due to reactions (from the last
      SDC iteration). This addition is done in
      after the source terms are
      converted to primitive variables.

#. *Update radial data and center of mass for monopole gravity*

   This is the same as the Strang step \ `[strang:radial] <#strang:radial>`__

#. *Clean State* []

   This is the same as the Strang step \ `[strang:clean] <#strang:clean>`__

#. [strang:newsource] *Correct the source terms with the :math:`n+1` contribution*
   [, ]

   This is the same as the Strang step \ `[strang:newsource] <#strang:newsource>`__

#. *Finalize* []

   This differs from Strang step \ `[strang:finalize] <#strang:finalize>`__ in that we do not
   construct :math:`d{\bf S}/dt`, but instead store the total hydrodynamical source
   term at the new time. As discussed above, this will be used in the
   next iteration to approximate the time-centered source term.

Runtime Parameters
==================

[chapter:parameters]

Introduction to Runtime Parameters
----------------------------------

Castro has 2 sets of runtime parameters—those controlled by
C and those controlled by Fortran. The C parameters are set
in the inputs file and managed by the AMReX ParmParse
class. For Castro-specific parameters, we list the runtime
parameters in a file and generate the
C code and headers using the script—note
this script needs to be run every time the \_cpp_parameters
file is updated.

The behavior of the network, EOS, and other microphysics routines are
controlled by a different set of runtime parameters. These parameters are defined
in plain-text files located in the different
directories that hold the microphysics code. At compile time, a
script in the AMReX bulid system, , locates all
of the \_parameters files that are needed for the given choice
of network, integrator, and EOS, and assembles all of the runtime
parameters into a module named (using the
script). The parameters are set in your
probin file in the &extern namelist.

C++ parameter format
~~~~~~~~~~~~~~~~~~~~

The C parameters take the form of:

::

    # comment describing the parameter
    name   type   default   need in Fortran?   ifdef    fortran name    fortran type

Here, name is the name of the parameter that will be looked for
in the inputs file, type is one of int, Real,
or string, and default is the default value of the
parameter. The next columns are optional, but you need to fill in all
of the information up to and including any of the optional columns you
need (e.g., if you are going to provide the fortran name, you
also need to provide need in Fortran? and ifdef. The
need in Fortran? column is y if the runtime parameter should
be made available in Fortran (through ).
The ifdef field provides the name of a preprocessor name that
should wrap this parameter definition—it will only be compiled in if
that name is defined to the preprocessor. The fortran name is
the name that the parameter should use in Fortran—by default it will
be the same as name. The fortran type is the data type of
the parameter in Fortran—by default it will be the
Fortran-equivalent to type. Finally, any comment immediately
before the parameter definition will be used to generate the LaTeX documentation
describing the parameters.

Microphysics/extern parameter format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The microphysics/extern parameter definitions take the form of:

::

    # comment describing the parameter
    name              data-type       default-value      priority

Here, the priority is simply an integer. When two directories
define the same parameter, but with different defaults, the version of
the parameter with the highest priority takes precedence. This allows
specific implementations to override the general parameter defaults.

The documentation below for the Castro C parameters is
automatically generated, using the comments in the \_cpp_parameters
file.

Removed Runtime Parameters
--------------------------

The following runtime parameters have been removed for Castro.

-  castro.ppm_flatten_before_integrals : this parameter
   controlled whether we applied the flattening of the parabolic
   profiles before we integrated under their profiles or afterwards.
   The default was switched to flattening before the integration,
   which is more consistent with the original PPM methodology. This
   parameter was removed since the variation enabled by this parameter
   was not that great.

   (removed in commit: 9cab697268997714919de16db1ca7e77a95c4f98)

-  castro.ppm_reference and
   castro.ppm_reference_edge_limit : these parameters controlled
   whether we use the integral under the parabola for the fastest wave
   moving toward the interface for the reference state and whether in
   the case that the wave is moving away from the interface we use the
   cell-center value or the limit of the parabola on the interface.
   These were removed because there is little reason to not use the
   reference state.

   (removed in commit: 213f4ffc53463141084551c7be4b37a2720229aa)

.. _ch:parameters:

 castro  Namespace
-----------------

.. raw:: latex

   \small

.. table:: castro : AMR
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | do we do the          | 1                     |
   |                       | hyperbolic reflux at  |                       |
   |    \endfoot           | coarse-fine           |                       |
   |                       | interfaces?           |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | how to do limiting of | 0                     |
   |                       | the state data when   |                       |
   |                       | interpolating 0: only |                       |
   |                       | prevent new extrema   |                       |
   |                       | 1: preserve linear    |                       |
   |                       | combinations of state |                       |
   |                       | variables             |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | highest order used in | 1                     |
   |                       | interpolation         |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Number of ghost zones | 0                     |
   |                       | for state data to     |                       |
   |                       | have. Note that if    |                       |
   |                       | you are using         |                       |
   |                       | radiation, choosing   |                       |
   |                       | this to be zero will  |                       |
   |                       | be overridden since   |                       |
   |                       | radiation needs at    |                       |
   |                       | least one ghost zone. |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | whether to re-compute | 1                     |
   |                       | new-time source terms |                       |
   |    \rowcolor{tableSha | after a reflux        |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | should we have state  | 0                     |
   |                       | data for custom       |                       |
   |                       | load-balancing        |                       |
   |                       | weighting?            |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : diagnostics, I/O
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | abort if we exceed    | 1                     |
   |                       | CFL = 1 over the      |                       |
   |    \endfoot           | cource of a timestep  |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | a string describing   | ""                    |
   |                       | the simulation that   |                       |
   |                       | will be copied into   |                       |
   |                       | the plotfile’s        |                       |
   |                       | job_info file         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | write a final         | 1                     |
   |                       | plotfile and          |                       |
   |    \rowcolor{tableSha | checkpoint upon       |                       |
   | de}                   | completion            |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | display warnings in   | (0, 1)                |
   |                       | Fortran90 routines    |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | display information   | (0, 1)                |
   |                       | about updates to the  |                       |
   |    \rowcolor{tableSha | state (how much mass, |                       |
   | de}                   | momentum, energy      |                       |
   |                       | added)                |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Do we want to reset   | -1                    |
   |                       | the number of steps   |                       |
   |                       | in the checkpoint?    |                       |
   |                       | This ONLY takes       |                       |
   |                       | effect if             |                       |
   |                       | amr.regrid_on_restart |                       |
   |                       | = 1 and               |                       |
   |                       | amr.checkpoint_on_res |                       |
   |                       | tart                  |                       |
   |                       | = 1, (which require   |                       |
   |                       | that max_step and     |                       |
   |                       | stop_time be less     |                       |
   |                       | than the value in the |                       |
   |                       | checkpoint) and you   |                       |
   |                       | set it to value       |                       |
   |                       | greater than this     |                       |
   |                       | default value.        |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Do we want to reset   | -1.e200               |
   |                       | the time in the       |                       |
   |    \rowcolor{tableSha | checkpoint? This ONLY |                       |
   | de}                   | takes effect if       |                       |
   |                       | amr.regrid_on_restart |                       |
   |                       | = 1 and               |                       |
   |                       | amr.checkpoint_on_res |                       |
   |                       | tart                  |                       |
   |                       | = 1, (which require   |                       |
   |                       | that max_step and     |                       |
   |                       | stop_time be less     |                       |
   |                       | than the value in the |                       |
   |                       | checkpoint) and you   |                       |
   |                       | set it to value       |                       |
   |                       | greater than this     |                       |
   |                       | default value.        |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | display center of     | 0                     |
   |                       | mass diagnostics      |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | how often (number of  | -1                    |
   |                       | coarse timesteps) to  |                       |
   |    \rowcolor{tableSha | compute integral sums |                       |
   | de}                   | (for runtime          |                       |
   |                       | diagnostics)          |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | how often (simulation | -1.0e0                |
   |                       | time) to compute      |                       |
   |                       | integral sums (for    |                       |
   |                       | runtime diagnostics)  |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | calculate losses of   | 0                     |
   |                       | material through      |                       |
   |    \rowcolor{tableSha | physical grid         |                       |
   | de}                   | boundaries            |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : diffusion
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | scaling factor for    | 1.0                   |
   |                       | conductivity          |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | set a cutoff density  | -1.e200               |
   |                       | for diffusion – we    |                       |
   |                       | zero the term out     |                       |
   |                       | below this density    |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | enable enthalpy       | 0                     |
   |                       | diffusion             |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | enable species        | 0                     |
   |                       | diffusion             |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | enable thermal        | 0                     |
   |                       | diffusion             |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | enable velocity       | 0                     |
   |                       | diffusion             |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : embiggening
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the factor by which   | 1                     |
   |                       | to extend the domain  |                       |
   |    \endfoot           | upon restart for      |                       |
   |                       | embiggening           |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | used with the         | -1                    |
   |                       | embiggening routines  |                       |
   |                       | to determine how to   |                       |
   |                       | extend the domain     |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : gravity and rotation
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | permits gravity       | -1                    |
   |                       | calculation to be     |                       |
   |    \endfoot           | turned on and off     |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | permits rotation      | -1                    |
   |                       | calculation to be     |                       |
   |                       | turned on and off     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | determines how the    | 4                     |
   |                       | gravitational source  |                       |
   |    \rowcolor{tableSha | term is added to the  |                       |
   | de}                   | momentum and energy   |                       |
   |                       | state variables.      |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | we can do a implicit  | 1                     |
   |                       | solution of the       |                       |
   |                       | rotation update to    |                       |
   |                       | allow for better      |                       |
   |                       | coupling of the       |                       |
   |                       | Coriolis terms        |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | to we recompute the   | 0                     |
   |                       | center used for the   |                       |
   |    \rowcolor{tableSha | multipole gravity     |                       |
   | de}                   | solve each step?      |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | mass of the point     | 0.0                   |
   |                       | mass                  |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if we have a central  | 0                     |
   |                       | point mass, we can    |                       |
   |    \rowcolor{tableSha | prevent mass from     |                       |
   | de}                   | building up in the    |                       |
   |                       | zones adjacent to it  |                       |
   |                       | by keeping their      |                       |
   |                       | density constant and  |                       |
   |                       | adding their mass to  |                       |
   |                       | the point mass object |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the coordinate axis   | 3                     |
   |                       | (:math:`x=1`,         |                       |
   |                       | :math:`y=2`,          |                       |
   |                       | :math:`z=3`) for the  |                       |
   |                       | rotation vector       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | determines how the    | 4                     |
   |                       | rotation source terms |                       |
   |    \rowcolor{tableSha | are added to the      |                       |
   | de}                   | momentum and energy   |                       |
   |                       | equations             |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | permits the           | 1                     |
   |                       | centrifugal terms in  |                       |
   |                       | the rotation to be    |                       |
   |                       | turned on and off     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | permits the Coriolis  | 1                     |
   |                       | terms in the rotation |                       |
   |    \rowcolor{tableSha | to be turned on and   |                       |
   | de}                   | off                   |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | permits the           | 1                     |
   |                       | d(omega)/dt terms in  |                       |
   |                       | the rotation to be    |                       |
   |                       | turned on and off     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the rotation periods  | 0.0                   |
   |                       | time evolution—this   |                       |
   |    \rowcolor{tableSha | allows the rotation   |                       |
   | de}                   | rate to change        |                       |
   |                       | durning the           |                       |
   |                       | simulation time       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the rotation period   | -1.e200               |
   |                       | for the corotating    |                       |
   |                       | frame                 |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Which reference frame | 1                     |
   |                       | to measure the state  |                       |
   |    \rowcolor{tableSha | variables with        |                       |
   | de}                   | respect to. The       |                       |
   |                       | standard in the       |                       |
   |                       | literature when using |                       |
   |                       | a rotating reference  |                       |
   |                       | frame is to measure   |                       |
   |                       | the state variables   |                       |
   |                       | with respect to an    |                       |
   |                       | observer fixed in     |                       |
   |                       | that rotating frame.  |                       |
   |                       | If this option is     |                       |
   |                       | disabled by setting   |                       |
   |                       | it to 0, the state    |                       |
   |                       | variables will be     |                       |
   |                       | measured with respect |                       |
   |                       | to an observer fixed  |                       |
   |                       | in the inertial frame |                       |
   |                       | (but the frame will   |                       |
   |                       | still rotate).        |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | include a central     | 1                     |
   |                       | point mass            |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : hydrodynamics
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if true, define an    | 0                     |
   |                       | additional source     |                       |
   |    \endfoot           | term                  |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Whether or not to     | 1                     |
   |                       | allow the internal    |                       |
   |                       | energy to be less     |                       |
   |                       | than the internal     |                       |
   |                       | energy corresponding  |                       |
   |                       | to small_temp         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | for the Colella &     | 2                     |
   |                       | Glaz Riemann solver,  |                       |
   |    \rowcolor{tableSha | what to do if we do   |                       |
   | de}                   | not converge to a     |                       |
   |                       | solution for the star |                       |
   |                       | state. 0 = do         |                       |
   |                       | nothing; print        |                       |
   |                       | iterations and exit 1 |                       |
   |                       | = revert to the       |                       |
   |                       | original guess for    |                       |
   |                       | p-star 2 = do a       |                       |
   |                       | bisection search for  |                       |
   |                       | another 2 \*          |                       |
   |                       | cg_maxiter            |                       |
   |                       | iterations.           |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | for the Colella &     | 12                    |
   |                       | Glaz Riemann solver,  |                       |
   |                       | the maximum number of |                       |
   |                       | iterations to take    |                       |
   |                       | when solving for the  |                       |
   |                       | star state            |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | for the Colella &     | 1.0e-5                |
   |                       | Glaz Riemann solver,  |                       |
   |    \rowcolor{tableSha | the tolerance to      |                       |
   | de}                   | demand in finding the |                       |
   |                       | star state            |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Which method to use   | 1                     |
   |                       | when resetting a      |                       |
   |                       | negative/small        |                       |
   |                       | density 1 = Reset to  |                       |
   |                       | characteristics of    |                       |
   |                       | adjacent zone with    |                       |
   |                       | largest density 2 =   |                       |
   |                       | Use average of all    |                       |
   |                       | adjacent zones for    |                       |
   |                       | all state variables 3 |                       |
   |                       | = Reset to the        |                       |
   |                       | original zone state   |                       |
   |                       | before the hydro      |                       |
   |                       | update                |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the coefficient of    | 0.1                   |
   |                       | the artificial        |                       |
   |    \rowcolor{tableSha | viscosity             |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | do we do the CTU      | 1                     |
   |                       | unsplit method or a   |                       |
   |                       | method-of-lines       |                       |
   |                       | approach?             |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | permits hydro to be   | -1                    |
   |                       | turned on and off for |                       |
   |    \rowcolor{tableSha | running pure rad      |                       |
   | de}                   | problems              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | permits sponge to be  | 0                     |
   |                       | turned on and off     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Threshold value of (E | 1.0e0                 |
   |                       | - K) / E such that    |                       |
   |    \rowcolor{tableSha | above eta1, the       |                       |
   | de}                   | hydrodynamic pressure |                       |
   |                       | is derived from E -   |                       |
   |                       | K; otherwise, we use  |                       |
   |                       | the internal energy   |                       |
   |                       | variable UEINT.       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Threshold value of (E | 1.0e-4                |
   |                       | - K) / E such that    |                       |
   |                       | above eta2, we update |                       |
   |                       | the internal energy   |                       |
   |                       | variable UEINT to     |                       |
   |                       | match E - K. Below    |                       |
   |                       | this, UEINT remains   |                       |
   |                       | unchanged.            |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | set the flattening    | 0                     |
   |                       | parameter to zero to  |                       |
   |    \rowcolor{tableSha | force the             |                       |
   | de}                   | reconstructed         |                       |
   |                       | profiles to be flat,  |                       |
   |                       | resulting in a        |                       |
   |                       | first-order method    |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       | 0                     |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | do we do fourth-order | 0                     |
   |                       | accurate MOL hydro?   |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if we are doing HSE   | 0                     |
   |                       | boundary conditions,  |                       |
   |                       | should we get the     |                       |
   |                       | temperature via       |                       |
   |                       | interpolation (using  |                       |
   |                       | model_parser) or hold |                       |
   |                       | it constant?          |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if we are doing HSE   | 0                     |
   |                       | boundary conditions,  |                       |
   |    \rowcolor{tableSha | how do we treat the   |                       |
   | de}                   | velocity? reflect? or |                       |
   |                       | outflow?              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if we are doing HSE   | 0                     |
   |                       | boundary conditions,  |                       |
   |                       | do we zero the        |                       |
   |                       | velocity?             |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | whether to use the    | 0                     |
   |                       | hybrid advection      |                       |
   |    \rowcolor{tableSha | scheme that updates   |                       |
   | de}                   | z-angular momentum,   |                       |
   |                       | cylindrical momentum, |                       |
   |                       | and azimuthal         |                       |
   |                       | momentum (3D only)    |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | do we drop from our   | 0                     |
   |                       | regular Riemann       |                       |
   |                       | solver to HLL when we |                       |
   |                       | are in shocks to      |                       |
   |                       | avoid the odd-even    |                       |
   |                       | decoupling            |                       |
   |                       | instability?          |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Should we limit the   | 0                     |
   |                       | density fluxes so     |                       |
   |    \rowcolor{tableSha | that we do not create |                       |
   | de}                   | small densities?      |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | integration order for | 2                     |
   |                       | MOL integration 1 =   |                       |
   |                       | first order, 2 =      |                       |
   |                       | second order TVD, 3 = |                       |
   |                       | 3rd order TVD, 4 =    |                       |
   |                       | 4th order RK          |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | for piecewise linear, | 2                     |
   |                       | reconstruction order  |                       |
   |    \rowcolor{tableSha | to use                |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | do we construct       | 0                     |
   |                       | :math:`\gamma_e = p/( |                       |
   |                       | \rho e) + 1`          |                       |
   |                       | and bring it to the   |                       |
   |                       | interfaces for        |                       |
   |                       | additional            |                       |
   |                       | thermodynamic         |                       |
   |                       | information (this is  |                       |
   |                       | the Colella & Glaz    |                       |
   |                       | technique) or do we   |                       |
   |                       | use :math:`(\rho e)`  |                       |
   |                       | (the classic          |                       |
   |                       | Castro behavior).     |                       |
   |                       | Note this also uses   |                       |
   |                       | :math:`\tau = 1/\rho` |                       |
   |                       | instead of            |                       |
   |                       | :math:`\rho`.         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | do we use the         | 0                     |
   |                       | reference state in    |                       |
   |    \rowcolor{tableSha | evaluating the        |                       |
   | de}                   | eigenvectors?         |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | various methods of    | 0                     |
   |                       | giving temperature a  |                       |
   |                       | larger role in the    |                       |
   |                       | reconstruction—see    |                       |
   |                       | Zingale & Katz 2015   |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | reconstruction type:  | 1                     |
   |                       | 0: piecewise linear;  |                       |
   |    \rowcolor{tableSha | 1: classic Colella &  |                       |
   | de}                   | Woodward ppm; 2:      |                       |
   |                       | extrema-preserving    |                       |
   |                       | ppm                   |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | which Riemann solver  | 0                     |
   |                       | do we use: 0:         |                       |
   |                       | Colella, Glaz, &      |                       |
   |                       | Ferguson (a two-shock |                       |
   |                       | solver); 1: Colella & |                       |
   |                       | Glaz (a two-shock     |                       |
   |                       | solver) 2: HLLC       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the small density     | -1.e200               |
   |                       | cutoff. Densities     |                       |
   |    \rowcolor{tableSha | below this value will |                       |
   | de}                   | be reset              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the small specific    | -1.e200               |
   |                       | internal energy       |                       |
   |                       | cutoff. Internal      |                       |
   |                       | energies below this   |                       |
   |                       | value will be reset   |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the small pressure    | -1.e200               |
   |                       | cutoff. Pressures     |                       |
   |    \rowcolor{tableSha | below this value will |                       |
   | de}                   | be reset              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the small temperature | -1.e200               |
   |                       | cutoff. Temperatures  |                       |
   |                       | below this value will |                       |
   |                       | be reset              |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | extrapolate the       | 0                     |
   |                       | source terms (gravity |                       |
   |    \rowcolor{tableSha | and rotation) to      |                       |
   | de}                   | :math:`n+1/2`         |                       |
   |                       | timelevel for use in  |                       |
   |                       | the interface state   |                       |
   |                       | prediction            |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if we are using the   | 1                     |
   |                       | sponge, whether to    |                       |
   |                       | use the implicit      |                       |
   |                       | solve for it          |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if the transverse     | 1                     |
   |                       | interface state       |                       |
   |    \rowcolor{tableSha | correction, if the    |                       |
   | de}                   | new density is        |                       |
   |                       | negative, then        |                       |
   |                       | replace all of the    |                       |
   |                       | interface quantities  |                       |
   |                       | with their values     |                       |
   |                       | without the           |                       |
   |                       | transverse            |                       |
   |                       | correction.           |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if the interface      | 0                     |
   |                       | state for             |                       |
   |                       | :math:`(\rho e)` is   |                       |
   |                       | negative after we add |                       |
   |                       | the transverse terms, |                       |
   |                       | then replace the      |                       |
   |                       | interface value of    |                       |
   |                       | :math:`(\rho e)` with |                       |
   |                       | a value constructed   |                       |
   |                       | from the              |                       |
   |                       | :math:`(\rho e)`      |                       |
   |                       | evolution equation    |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | after we add the      | 0                     |
   |                       | transverse correction |                       |
   |    \rowcolor{tableSha | to the interface      |                       |
   | de}                   | states, replace the   |                       |
   |                       | predicted pressure    |                       |
   |                       | with an EOS call      |                       |
   |                       | (using :math:`e` and  |                       |
   |                       | :math:`\rho`).        |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | should we use the EOS | 0                     |
   |                       | in the Riemann solver |                       |
   |                       | to ensure             |                       |
   |                       | thermodynamic         |                       |
   |                       | consistency?          |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | flatten the           | 1                     |
   |                       | reconstructed         |                       |
   |    \rowcolor{tableSha | profiles around       |                       |
   | de}                   | shocks to prevent     |                       |
   |                       | them from becoming    |                       |
   |                       | too thin              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | for the piecewise     | 1                     |
   |                       | linear                |                       |
   |                       | reconstruction, do we |                       |
   |                       | subtract off          |                       |
   |                       | :math:`(\rho g)` from |                       |
   |                       | the pressure before   |                       |
   |                       | limiting?             |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if we are doing an    | ""                    |
   |                       | external -x boundary  |                       |
   |    \rowcolor{tableSha | condition, who do we  |                       |
   | de}                   | interpret it?         |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if we are doing an    | ""                    |
   |                       | external +x boundary  |                       |
   |                       | condition, who do we  |                       |
   |                       | interpret it?         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if we are doing an    | ""                    |
   |                       | external -y boundary  |                       |
   |    \rowcolor{tableSha | condition, who do we  |                       |
   | de}                   | interpret it?         |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if we are doing an    | ""                    |
   |                       | external +y boundary  |                       |
   |                       | condition, who do we  |                       |
   |                       | interpret it?         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if we are doing an    | ""                    |
   |                       | external -z boundary  |                       |
   |    \rowcolor{tableSha | condition, who do we  |                       |
   | de}                   | interpret it?         |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if we are doing an    | ""                    |
   |                       | external +z boundary  |                       |
   |                       | condition, who do we  |                       |
   |                       | interpret it?         |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : parallelization
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1                     |
   |                       |                       |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | determines whether we | -1                    |
   |                       | use accelerators for  |                       |
   |                       | specific loops        |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : particles
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | permits tracer        | 0                     |
   |                       | particle calculation  |                       |
   |    \endfoot           | to be turned on and   |                       |
   |                       | off                   |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : reactions
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | disable burning       | 0                     |
   |                       | inside hydrodynamic   |                       |
   |    \endfoot           | shock regions         |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | permits reactions to  | -1                    |
   |                       | be turned on and off  |                       |
   |                       | – mostly for          |                       |
   |                       | efficiency’s sake     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Limit the timestep    | 1.e200                |
   |                       | based on how much the |                       |
   |    \rowcolor{tableSha | burning can change    |                       |
   | de}                   | the species mass      |                       |
   |                       | fractions of a zone.  |                       |
   |                       | The timestep is equal |                       |
   |                       | to dtnuc              |                       |
   |                       | :math:`\cdot\,(X / \d |                       |
   |                       | ot{X})`.              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | If we are using the   | 1.e-3                 |
   |                       | timestep limiter      |                       |
   |                       | based on changes in   |                       |
   |                       | :math:`X`, set a      |                       |
   |                       | threshold on the      |                       |
   |                       | species abundance     |                       |
   |                       | below which the       |                       |
   |                       | limiter is not        |                       |
   |                       | applied. This helps   |                       |
   |                       | prevent the timestep  |                       |
   |                       | from becoming very    |                       |
   |                       | small due to changes  |                       |
   |                       | in trace species.     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Limit the timestep    | 1.e200                |
   |                       | based on how much the |                       |
   |    \rowcolor{tableSha | burning can change    |                       |
   | de}                   | the internal energy   |                       |
   |                       | of a zone. The        |                       |
   |                       | timestep is equal to  |                       |
   |                       | dtnuc                 |                       |
   |                       | :math:`\cdot\,(e / \d |                       |
   |                       | ot{e})`.              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | limit the zone size   | 1.e200                |
   |                       | based on how much the |                       |
   |                       | burning can change    |                       |
   |                       | the internal energy   |                       |
   |                       | of a zone. The zone   |                       |
   |                       | size on the finest    |                       |
   |                       | level must be smaller |                       |
   |                       | than dxnuc            |                       |
   |                       | :math:`\cdot\, c_s\cd |                       |
   |                       | ot (e / \dot{e})`,    |                       |
   |                       | where :math:`c_s` is  |                       |
   |                       | the sound speed. This |                       |
   |                       | ensures that the      |                       |
   |                       | sound-crossing time   |                       |
   |                       | is smaller than the   |                       |
   |                       | nuclear energy        |                       |
   |                       | injection timescale.  |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Disable limiting      | 1.e200                |
   |                       | based on dxnuc above  |                       |
   |    \rowcolor{tableSha | this threshold. This  |                       |
   | de}                   | allows zones that     |                       |
   |                       | have already ignited  |                       |
   |                       | or are about to       |                       |
   |                       | ignite to be          |                       |
   |                       | de-refined.           |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Disable limiting      | -1                    |
   |                       | based on dxnuc above  |                       |
   |                       | this AMR level.       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | maximum temperature   | 1.e200                |
   |                       | for allowing          |                       |
   |    \rowcolor{tableSha | reactions to occur in |                       |
   | de}                   | a zone                |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | minimum temperature   | 0.0                   |
   |                       | for allowing          |                       |
   |                       | reactions to occur in |                       |
   |                       | a zone                |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | maximum density for   | 1.e200                |
   |                       | allowing reactions to |                       |
   |    \rowcolor{tableSha | occur in a zone       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | minimum density for   | 0.0                   |
   |                       | allowing reactions to |                       |
   |                       | occur in a zone       |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : refinement
parameters

   +--------------------------+--+---+
   |                          |  |   |
   +--------------------------+--+---+
   | Table —continued         |  |   |
   +--------------------------+--+---+
   |                          |  |   |
   +--------------------------+--+---+
   |                          |  |   |
   +--------------------------+--+---+
   | .. raw:: latex           |  | 0 |
   |                          |  |   |
   |    \endfoot              |  |   |
   |                          |  |   |
   | .. raw:: latex           |  |   |
   |                          |  |   |
   |    \hline                |  |   |
   |                          |  |   |
   | .. raw:: latex           |  |   |
   |                          |  |   |
   |    \endlastfoot          |  |   |
   |                          |  |   |
   | .. raw:: latex           |  |   |
   |                          |  |   |
   |    \rowcolor{tableShade} |  |   |
   +--------------------------+--+---+
   |                          |  | 0 |
   +--------------------------+--+---+

.. raw:: latex

   \small

.. table:: castro : timestep control
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the effective Courant | 0.8                   |
   |                       | number to use—we will |                       |
   |    \endfoot           | not allow the         |                       |
   |                       | hydrodynamic waves to |                       |
   | .. raw:: latex        | cross more than this  |                       |
   |                       | fraction of a zone    |                       |
   |    \hline             | over a single         |                       |
   |                       | timestep              |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the maximum factor by | 1.1                   |
   |                       | which the timestep    |                       |
   |                       | can increase from one |                       |
   |                       | step to the next.     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | If we do request more | 1                     |
   |                       | than the maximum      |                       |
   |    \rowcolor{tableSha | number of subcycles,  |                       |
   | de}                   | should we fail, or    |                       |
   |                       | should we clamp to    |                       |
   |                       | that maximum number   |                       |
   |                       | and perform that      |                       |
   |                       | many?                 |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the smallest valid    | 0.0                   |
   |                       | timestep—if we go     |                       |
   |                       | below this, we abort  |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | a fixed timestep to   | -1.0                  |
   |                       | use for all steps     |                       |
   |    \rowcolor{tableSha | (negative turns it    |                       |
   | de}                   | off)                  |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | a factor by which to  | 1.0                   |
   |                       | reduce the first      |                       |
   |                       | timestep from that    |                       |
   |                       | requested by the      |                       |
   |                       | timestep estimators   |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the initial timestep  | -1.0                  |
   |                       | (negative uses the    |                       |
   |    \rowcolor{tableSha | step returned from    |                       |
   | de}                   | the timestep          |                       |
   |                       | constraints)          |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the largest valid     | 1.e200                |
   |                       | timestep—limit all    |                       |
   |                       | timesteps to be no    |                       |
   |                       | larger than this      |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Do not permit more    | 10                    |
   |                       | subcycled timesteps   |                       |
   |    \rowcolor{tableSha | than this parameter.  |                       |
   | de}                   | Set to a negative     |                       |
   |                       | value to disable this |                       |
   |                       | criterion.            |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | enforce that the AMR  | 0                     |
   |                       | plot interval must be |                       |
   |                       | hit exactly           |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | If we’re doing        | 1.e-1                 |
   |                       | retries, set the      |                       |
   |    \rowcolor{tableSha | target threshold for  |                       |
   | de}                   | changes in density if |                       |
   |                       | a retry is triggered  |                       |
   |                       | by a negative         |                       |
   |                       | density. If this is   |                       |
   |                       | set to a negative     |                       |
   |                       | number then it will   |                       |
   |                       | disable retries using |                       |
   |                       | this criterion.       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | When performing a     | 0.5                   |
   |                       | retry, the factor to  |                       |
   |                       | multiply the current  |                       |
   |                       | timestep by when      |                       |
   |                       | trying again.         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Tolerance to use when | 0.02                  |
   |                       | evaluating whether to |                       |
   |    \rowcolor{tableSha | do a retry. The       |                       |
   | de}                   | timestep suggested by |                       |
   |                       | the retry will be     |                       |
   |                       | multiplied by (1 +    |                       |
   |                       | this factor) before   |                       |
   |                       | comparing the actual  |                       |
   |                       | timestep to it. If    |                       |
   |                       | set to some number    |                       |
   |                       | slightly larger than  |                       |
   |                       | zero, then this       |                       |
   |                       | prevents retries that |                       |
   |                       | are caused by small   |                       |
   |                       | numerical             |                       |
   |                       | differences.          |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Number of iterations  | 2                     |
   |                       | for the SDC advance.  |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | enforce that the AMR  | 0                     |
   |                       | small plot interval   |                       |
   |    \rowcolor{tableSha | must be hit exactly   |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Check for a possible  | 0                     |
   |                       | post-timestep regrid  |                       |
   |                       | if certain stability  |                       |
   |                       | criteria were         |                       |
   |                       | violated.             |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Retry a timestep if   | 0                     |
   |                       | it violated the       |                       |
   |    \rowcolor{tableSha | timestep-limiting     |                       |
   | de}                   | criteria over the     |                       |
   |                       | course of an advance. |                       |
   |                       | The criteria will     |                       |
   |                       | suggest a new         |                       |
   |                       | timestep that         |                       |
   |                       | satisfies the         |                       |
   |                       | criteria, and we will |                       |
   |                       | do subcycled          |                       |
   |                       | timesteps on the same |                       |
   |                       | level until we reach  |                       |
   |                       | the original target   |                       |
   |                       | time.                 |                       |
   +-----------------------+-----------------------+-----------------------+

.. _ch:parameters:

 diffusion  Namespace
--------------------

.. raw:: latex

   \small

.. table:: diffusion parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Use MLMG as the       | 4                     |
   |                       | operator              |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the level of          | 0                     |
   |                       | verbosity for the     |                       |
   |                       | diffusion solve       |                       |
   |                       | (higher number means  |                       |
   |                       | more output)          |                       |
   +-----------------------+-----------------------+-----------------------+

.. _ch:parameters:

 gravity  Namespace
------------------

.. raw:: latex

   \small

.. table:: gravity parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if doing constant     | 0.0                   |
   |                       | gravity, what is the  |                       |
   |    \endfoot           | acceleration          |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Check if the user     | 0                     |
   |                       | wants to compute the  |                       |
   |                       | boundary conditions   |                       |
   |                       | using the brute force |                       |
   |                       | method. Default is    |                       |
   |                       | false, since this     |                       |
   |                       | method is slow.       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | should we apply a     | 1                     |
   |                       | lagged correction to  |                       |
   |    \rowcolor{tableSha | the potential that    |                       |
   | de}                   | gets us closer to the |                       |
   |                       | composite solution?   |                       |
   |                       | This makes the        |                       |
   |                       | resulting fine grid   |                       |
   |                       | calculation slightly  |                       |
   |                       | more accurate, at the |                       |
   |                       | cost of an additional |                       |
   |                       | Poisson solve per     |                       |
   |                       | timestep.             |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | ratio of dr for       | 1                     |
   |                       | monopole gravity      |                       |
   |                       | binning to grid       |                       |
   |                       | resolution            |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | For non-Poisson       | 0                     |
   |                       | gravity, do we want   |                       |
   |    \rowcolor{tableSha | to construct the      |                       |
   | de}                   | gravitational         |                       |
   |                       | acceleration by       |                       |
   |                       | taking the gradient   |                       |
   |                       | of the potential,     |                       |
   |                       | rather than           |                       |
   |                       | constructing it       |                       |
   |                       | directly?             |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | what type             | "fillme"              |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the maximum mulitpole | 0                     |
   |                       | order to use for      |                       |
   |    \rowcolor{tableSha | multipole BCs when    |                       |
   | de}                   | doing Poisson gravity |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | For all gravity       | MAX_LEV-1             |
   |                       | types, we can choose  |                       |
   |                       | a maximum level for   |                       |
   |                       | explicitly            |                       |
   |                       | calculating the       |                       |
   |                       | gravity and           |                       |
   |                       | associated potential. |                       |
   |                       | Above that level, we  |                       |
   |                       | interpolate from      |                       |
   |                       | coarser levels.       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Do agglomeration?     | 1                     |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       | 1                     |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | how many FMG cycles?  | 0                     |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Do N-Solve?           | 0                     |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | do we do a composite  | 0                     |
   |                       | solve?                |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | do we perform the     | 0                     |
   |                       | synchronization at    |                       |
   |                       | coarse-fine           |                       |
   |                       | interfaces?           |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the level of          | 0                     |
   |                       | verbosity for the     |                       |
   |    \rowcolor{tableSha | gravity solve (higher |                       |
   | de}                   | number means more     |                       |
   |                       | output on the status  |                       |
   |                       | of the solve /        |                       |
   |                       | multigrid             |                       |
   +-----------------------+-----------------------+-----------------------+

.. _ch:parameters:

 particles  Namespace
--------------------

.. raw:: latex

   \small

.. table:: particles parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the name of an input  | ""                    |
   |                       | file containing the   |                       |
   |    \endfoot           | total particle number |                       |
   |                       | and the initial       |                       |
   | .. raw:: latex        | position of each      |                       |
   |                       | particle.             |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the name of timestamp | ""                    |
   |                       | files.                |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the name of a file    | ""                    |
   |                       | with new particles at |                       |
   |    \rowcolor{tableSha | restart               |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | to restart from a     | 0                     |
   |                       | checkpoint that was   |                       |
   |                       | written with          |                       |
   |                       | USE_PARTICLES=FALSE   |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | whether the local     | 1                     |
   |                       | densities at given    |                       |
   |    \rowcolor{tableSha | positions of          |                       |
   | de}                   | particles are stored  |                       |
   |                       | in output files       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the name of a         | ""                    |
   |                       | directory in which    |                       |
   |                       | timestamp files are   |                       |
   |                       | stored.               |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | whether the local     | 0                     |
   |                       | temperatures at given |                       |
   |    \rowcolor{tableSha | positions of          |                       |
   | de}                   | particles are stored  |                       |
   |                       | in output files       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the level of          | 0                     |
   |                       | verbosity for the     |                       |
   |                       | tracer particle (0 or |                       |
   |                       | 1)                    |                       |
   +-----------------------+-----------------------+-----------------------+

Compilers, Managing Jobs, and Scaling
=====================================

Castro requires both a C compiler (with support for C 11) and
a Fortran compiler with 2003+ support.

General Compilers
-----------------

GCC
~~~

The GCC compilers are the preferred compilers on Linux machines.
Castro runs without issue with GCC 7.x and GCC 8.x.

Intel
~~~~~

Intel compilers 17.x, 18.x, and 19.x produce internal compiler errors
and should not be used.

PGI
~~~

The PGI compilers (16.4–16.10) are also tested with Castro and have
no known issues. Since PGI uses the GCC header files, a compatible
version of GCC is needed. This is particularly problematic with C 11
support. At the moment, with PGI 16.10 compilers the latest GCC compiler
that has compatible C 11 support is GCC 4.9.4.

Working at OLCF (ORNL)
----------------------

The build system needs python 2.7+ or 3.5+. The default python at OLCF
is too old, so you will need to manually load a newer python as:

::

    module load python

Cray
~~~~

Note, you will need to swap the default PGI environment to use Cray:

::

    module swap PrgEnv-pgi PrgEnv-cray

cce/8.5.7 through 8.6.3 fail to compile AMReX. The default version at OLCF as of May 2018,
cce/8.6.4, resolves this issue.

.. _pgi-1:

PGI
~~~

If you are using OpenACC on the GPUs, then you need to use the PGI
compilers instead. It is best to use the latest available PGI
compilers, which are 17.7 at this writing.

::

    module swap pgi pgi/17.7

These seem to build the code without problems.

This is perhaps out of date:
Finally, to use OpenACC, you need to make the CUDA libraries available to the PGI compiler:

::

    module load craype-accel-nvidia35

Automatic Restarting and Archiving of Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See the Maestro User’s Guide

Working at NERSC
----------------

General issues
~~~~~~~~~~~~~~

Edison compilers
~~~~~~~~~~~~~~~~

-  Intel compilers

   The default compilers on edison are the Intel compilers.

   The Intel 17.0.2 compilers produce code that crashes at runtime when self-gravity
   is used (the error is in the F_MG directory. This appears to be an issue since
   16.0.2.

   At the moment, avoid Intel compilers on Edison.

-  Cray compilers

   -  The Cray compilers from versions 8.5.7 to 8.6.3 have issues compiling AMReX. This
      can be resolved by loading version 8.6.4 or later. As of May 2018, you can do

      ::

                module swap cce cce/8.6.5
              

-  GNU

   The GNU compilers (6.3) seem to work fine.

Cori (KNL) compilers
~~~~~~~~~~~~~~~~~~~~

To compile with Cray compilers on Cori, we need to swap the compiler
wrappers to use the AVX-512 instruction set supported on the Intel Phi
processors instead of the AVX-2 extensions used by the Intel Haswell
architecture. This is done as:

::

    module swap craype-{haswell,mic-knl}

It could happen that even when the various verbosities are set to 0, when using several nodes (more than 64) in a run compiled with Intel, the follwing error shows:

::

    "forrtl: severe (40): recursive I/O operation, unit -1, file unknown"

Seems like the error is due to all threads printing to stdout. Adding the following to the inputs file, prevents this error to occur:

::

    castro.print_fortran_warnings = 0

Hypre and radiation
~~~~~~~~~~~~~~~~~~~

On edison, the Cray *Third Party Scientific Libraries* provide
hypre in a form that works directly with the compiler wrappers
used on that machine (CC, ftn, :math:`\ldots`). To use this,
simply do:

::

    module load cray-tpsl

There is no need to set HYPRE_DIR, but note however that the
dependency checker script (BoxLib/Tools/C_scripts/mkdep) will
complain about:

::

    /path/to/Hypre--with-openmp/include does not exist

This can be ignored an compilation will finish. If you do wish to
silence it, you can set HYPRE_DIR to the path shown by

::

    module show cray-tpsl

as

::

    export HYPRE_DIR=${CRAY_TPSL_PREFIX_DIR}

This path will change dynamically to reflect which compiler programming
environment you have loaded. (You can also see that this is the path
sent to the compilation by doing ftn -craype-verbose).

Running jobs
~~~~~~~~~~~~

edison is configured with 24 cores per node split between two Intel
IvyBridge 12-core processors. Each processor connects to 1/2 of the
node’s memory and is called a NUMA node, so there are 2 NUMA nodes per
edison node. Best performance is seen when running with 6 or 12 threads.

Jobs should be run in your $SCRATCH or $CSCATCH directory.
By default, SLURM will change directory into the submission directory.

A sample job submission script, edison.MPI.OMP.slurm is in
Castro/Util/job_scripts/edison/, and includes logic to
automatically add the correct restart options to the run to continue a
simulation from the last checkpoint file in the submission directory.

To chain jobs, such that one queues up after the previous job finished,
use the chainslurm.sh script in that same directory:

::

    chainslurm.sh jobid number script

where jobid is the existing job you want to start you chain
from, number is the number of new jobs to chain from this
starting job, and script is the job submission script to use
(the same one you used originally most likely). You can view the job
dependency using:

::

    squeue -l -j job-id                                                             

where job-id is the number of the job.

Jobs are submitted with sbatch. A job can be canceled using
scancel, and the status can be checked using squeue -u
*username*.

Archiving data to HPSS
~~~~~~~~~~~~~~~~~~~~~~

The script edison.xfer.slurm in
Castro/Util/job_scripts/edison/ can be used to archive data to
HPSS automatically. This is submitted to the xfer queue and
runs the script process.xrb which continually looks for output
and stores it to HPSS.

To use the scripts, first create a directory in HPSS that has the same
name as the directory on lustre you are running in (just the directory
name, not the full path). E.g. if you are running in a directory
call wdconvect_run, then do:

::

    hsi                                                                             
    mkdir wdconvect_run                                                             

(Note: if the hsi command prompts you for your password, you will need
to talk to the NERSC help desk to ask for password-less access to
HPSS).

The script process.xrb is called from the xfer job and will
run in the background and continually wait until checkpoint or
plotfiles are created (actually, it always leaves the most recent one
alone, since data may still be written to it, so it waits until there
are more than 1 in the directory).

Then the script will use htar to archive the plotfiles and
checkpoints to HPSS. If the htar command was successful, then
the plotfiles are copied into a plotfile/ subdirectory. This is
actually important, since you don’t want to try archiving the data a
second time and overwriting the stored copy, especially if a purge
took place. The same is done with checkpoint files.

Additionally, if the ftime executable is in your path (
ftime.f90 lives in BoxLib/Tools/Postprocessing/F_src/), then
the script will create a file called ftime.out that lists the
name of the plotfile and the corresponding simulation time.

Finally, right when the job is submitted, the script will tar up all
of the diagnostic files, ftime.out, submission script, inputs
and probin, and archive them on HPSS. The .tar file is given a
name that contains the date-string to allow multiple archives to
co-exist.

When process.xrb is running, it creates a lockfile (called
process.pid) that ensures that only one instance of the script
is running at any one time. Sometimes if the machine crashes, the
process.pid file will be left behind, in which case, the script
aborts. Just delete that if you know the script is not running.

Jobs in the xfer queue start up quickly. The best approach is
to start one as you start your main job (or make it dependent on the
main job). The sample process.xrb script will wait for output
and then archive it as it is produced, using the techniques described
for titan above.

To check the status of a job in the xfer queue, use:

::

    squeue -u username -M all                                                       

Working at LANL
---------------

For the following LANL systems, which all have access to a joint file system,

-  Cielito

-  Conejo

-  Lightshow

-  Moonlight

-  Pinto

-  Wolf

-  Mustang

-  Trinitite

the following steps are needed to get Castro compiling (reported by Platon Karpov, 6/13/2016):

| 1) Compile on the login node (thus host is lanl.gov)
| 2) Do *not* set env variable ``BOXLIB_USE_MPI_WRAPPERS`` at all
| 3) Execute ``module load python-epd`` at the command line, but do *not* load anaconda

Scaling
-------

Data from scaling studies is archived in Castro/Docs/ManagingJobs/scaling/

Needs to be updated

Gotyas
------

#. 3/4/16: The default version loaded of Python on Mira is not
   recent enough to support the Python scripts in our build system. Add
   +python to your .soft to fix this.

#. 2/18/16: The default version loaded of Python on Titan (2.6.9)
   is not recent enough to support the Python scripts in our build
   system. At the terminal, do module load python to fix this.

GPUs
----

Bender
------

Compile as:

::

    make CUDA_VERSION=cc60 COMPILE_CUDA_PATH=/usr/local/cuda-9.2 USE_CUDA=TRUE COMP=PGI -j 4

To run the CUDA code path without device launching, do:

::

    make -j4 COMP=PGI USE_CUDA=TRUE USE_MPI=FALSE DEBUG=TRUE NO_DEVICE_LAUNCH=TRUE CUDA_VERSION=cc60 COMPILE_CUDA_PATH=/usr/local/cuda-9.2

Frequently Asked Questions
==========================

Compiling
---------

#. *Compiling fails giving me a cryptic message about a module not
   being found, usually bl_types or bl_error_module, like:*

   ::

       mpif90 -fno-range-check -fno-second-underscore
        -Jo/3d.Linux.gcc.gfortran.MPI.EXE -I o/3d.Linux.gcc.gfortran.MPI.EXE
        -ffixed-line-length-0 -g -O3
        -I. -I/home/zingale/development/Microphysics/util
        -I../../Microphysics/EOS -I../../Microphysics/EOS/gamma_law
        -I../../Microphysics/networks
        -I../../Microphysics/networks/general_null
        -I. -I/home/zingale/development/BoxLib//Src/C_BaseLib
        -I/home/zingale/development/BoxLib//Src/C_AMRLib
        -I/home/zingale/development/BoxLib//Src/C_BoundaryLib -I../../Source
        -I../../Source/Src_3d -I../../Source/Src_nd -I../../constants
        -I/home/zingale/development/BoxLib//Src/F_BaseLib
        -I/home/zingale/development/BoxLib//Tools/C_scripts -c
        ../../Microphysics/EOS/eos.f90 -o
        o/3d.Linux.gcc.gfortran.MPI.EXE/eos.o
       ../../Microphysics/EOS/eos.f90:3:6:

          use bl_types
             1
       Fatal Error: Can't open module file ‘bl_types.mod’ for reading at (1):
        No such file or directory compilation terminated.

       /home/zingale/development/BoxLib//Tools/C_mk/Make.rules:122: recipe
        for target 'o/3d.Linux.gcc.gfortran.MPI.EXE/eos.o' failed

   This usually indicates that the build system cannot find a source file
   (note: the problem above is not bl_types, that just seems to be
   the way the error manifests itself). The source files are specified
   in the various Make.package files throughout the
   Castro directory hierarchy. make will look through the
   directories in the VPATH_LOCATIONS to find the files.

   There are 2 things you can do to check what’s happening. First, inspect
   the directories in VPATH_LOCATIONS. This can be done via:

   ::

       make print-VPATH_LOCATIONS

   Next, ask make to tell you where it is finding each of the source
   files. This is done through a script find_files_vpath.py
   that is hooked into Castro’s build system. You can run this as:

   ::

       make file_locations

   At the end of the report, it will list any files it cannot find in
   the vpath. Some of these are to be expected (like extern.f90
   and buildInfo.cpp—these are written at compiletime. But any
   other missing files need to be investigated.

#. *I’m still having trouble compiling. How can I find out what
   all of the make variables are set to?*

   Use:

   ::

       make help

   This will tell you the value of all the compilers and their options.

#. *How do I use a system’s BLAS library instead of compiling and
   linking the one that comes with the StarKiller microphysics?*

   To use a system’s BLAS library, set the Make variable
   to TRUE. This will then look at
   the Make variable for the library to link
   (defaults to -lopenblas).

#. *How can I check to make sure the function signatures defined
   in C are consistent with their implementations in Fortran?*

   Use:

   ::

       make typecheck

   This will compile the code and report on any mismatched function signatures.

Debugging
---------

#. *Castro crashes with a floating point exception—how can
   I get more information?*

   The best thing to do is to recompile the code with TEST=TRUE
   set in the GNUmakefile. This will have AMReX catch the
   signals raised in both C and Fortran functions. Behind the
   scenes, this defines the BL_TESTING preprocessor flag, which
   will initialize memory allocated in fabs or multifabs to
   signaling NaNs (sNaN), and use the BLBackTrace::handler()
   function to handle various signals raised in both C and Fortran
   functions. This is a Linux/UNIX capability. This gives us a chance
   to print out backtrace information. The signals include seg fault,
   floating point exceptions (NaNs, divided by zero and overflow), and
   interruption by the user and system. What signals are handed to
   AMReX are controlled by AMReX(e.g., using interruption by the
   user, this was once used to find an MPI deadlock.) It also includes
   the BL_ASSERTION statements if USE_ASSERTION=TRUE or
   DEBUG=TRUE.

   The AMReX parameters that affect the behavior are:

   -  amrex.fpe_trap_invalid

   -  amrex.fpe_trap_zero

   -  amrex.fpe_trap_overflow

   For further capabilities, defining BACKTRACE=TRUE enables you
   to get more information than the backtrace of the call stack info by
   instrumenting the code. (This is in C code only). Here is an
   example. You know the line “Real rho = state(cell,0);” is
   causing a segfault. You could add a print statement before that.
   But it might print out thousands (or even millions) of line before
   it hits the segfault. With BACKTRACE, you could do

   ::

             #ifdef BL_BACKTRACING
                std::ostringstream ss;
                ss << ``state.box() = `` << state.box() << `` cell = `` << cell;
                BL_BACKTRACE_PUSH(ss.str()); // PUSH takes std::string
             #endif

             Real rho = state(cell,0);  // state is a Fab, and cell is an IntVect.

   The “print” prints to a stack of string, not stdout. When it
   hits the segfault, you will only see the last print out.

#. *How can I monitor the state in a zone from the C side
   at various points in the evolution?*

   Given a MultiFab mf, you can dump out the state as:

   ::

           print_state(mf, IntVect(D_DECL(10, 20, 30)));

   Here, the IntVect has the dimension that we were compiled with
   (and this is handled through the preprocessor D_DECL). In
   this case, we are inspecting zone (10, 20, 30), in the global index
   space. Note that since a multifab exists only on a single level, the
   integer indices here refer to the global index space on that level.

#. *What if I want to see all the data in a FArrayBox?*

   You can simply output a FAB to std::cout. Imagine that you
   are in an MFIter loop, with a MultiFab mf:

   ::

           S = FArrayBox& mf[mfi];
           std::cout << S << std::endl;

   This will output the contents on the FAB, one zone per line.

Profiling
---------

#. *How can I get line-by-line profiling information?*

   With the GNU compliers, you can enabling profiling with gprof
   by compiling with

   ::

         USE_GPROF=TRUE

   in your GNUmakefile.

   When you run, a file named gmon.out will be produced. This can
   be processed with gprof by running:

   ::

         gprof exec-name

   where *exec-name* is the name of the executable. More detailed
   line-by-line information can be obtained by passing the -l
   argument to gprof.

Managing Runs
-------------

#. *How can I force the running code to output, even it the plot or
   checkpoint interval parameters don’t require it?*

   Create a file called , e.g., as:

   ::

       touch dump_and_continue

   This will force the code to output a checkpoint file that can be used
   to restart. Other options are to output
   a plotfile, to output a checkpoint file
   and halt the code, and to simply stop the code.
   Note that the parameter controls how often
   the existence of these files is checked; by default it is 10, so the
   check will be done at the end of every timestep that is a multiple of 10.
   Set that to 1 in your inputs file if you’d like it to check every timestep.

#. *How can I output plotfiles in single precision?*

   The AMReX runtime parameter:

   ::

       fab.format = NATIVE_32

   controls this (put this in your inputs file). Note: checkpoint files are unaffected
   by this and will always be written out in the native precision (the ‘fab.format‘ parameter
   is overridden in the checkpoint code in AMReX).

.. _ch:faq:vis:

Visualization
-------------

#. *When I try to use Amrvis with the Nvidia driver, all I see is
   black—no data. How do I fix this?*

   You need to edit your xorg.conf file (usually found in /etc/X11/
   to enable the Dac8Bit option. The section will look like:

   ::

       Section "Device"
           Identifier     "Device0"
           Driver         "nvidia"
           VendorName     "NVIDIA Corporation"
           Option         "Dac8bit" "True"
       EndSection

   If you don’t already have an xorg.conf then you can create one
   by running nvidia-xconfig first.

Hydrodynamics
=============

.. _introduction-2:

Introduction
------------

The hydrodynamics scheme in Castro implements an unsplit
second-order Godunov method. Characteristic tracing is used to
time-center the input states to the Riemann solver. The same
hydrodynamics routines are used for pure hydro and radiation
hydrodynamics.

Some general notes:

-  Regardless of the dimensionality, we always carry around all 3
   components of velocity/momentum—this allows for rotation sources easily.

-  When radiation is enabled (via ), we discuss
   the gas and radiation quantities separately. This generally applies
   to the temperature, pressure, internal energy, various adiabatic
   indices, and sound speed. When we refer to the “total” value of
   one of these, it means that both gas and radiation contributions
   are included. When we refer to the “gas” quantity, this is what
   the equation of state would return.

   For continuity, we continue to use the “gas” qualifier even if we
   are not solving the radiation hydrodynamics equations. In this
   case, it still means that it comes through the equation of state,
   but note some of our equations of state (like the helmeos) include a
   radiation pressure contribution when we are running without
   radiation hydrodynamics enabled. In this case, we still refer to
   this as the “gas”.

Hydrodynamics Data Structures
-----------------------------

Within the Fortran routines that implement the hydrodynamics, there are
several main data structures that hold the state.

-  conserved state: these arrays generally begin with u,
   e.g., , . The
   components for the state data in the array are accessed using
   integer keys defined in `[table:consints] <#table:consints>`__.

   .. raw:: latex

      \centering

   .. table:: [table:consints] The integer variables to index the conservative state array

      +-----------------------+-----------------------+-----------------------+
      | **variable**          | **quantity**          | **note**              |
      +=======================+=======================+=======================+
      | URHO                  | :math:`\rho`          |                       |
      +-----------------------+-----------------------+-----------------------+
      | UMX                   | :math:`\rho u`        |                       |
      +-----------------------+-----------------------+-----------------------+
      | UMY                   | :math:`\rho v`        |                       |
      +-----------------------+-----------------------+-----------------------+
      | UMZ                   | :math:`\rho w`        |                       |
      +-----------------------+-----------------------+-----------------------+
      | UEDEN                 | :math:`\rho E`        |                       |
      +-----------------------+-----------------------+-----------------------+
      | UEINT                 | :math:`\rho e`        | this is computed from |
      |                       |                       | the other quantities  |
      |                       |                       | using                 |
      |                       |                       | :math:`\rho e = \rho  |
      |                       |                       | E - \rho {\bf u}\cdot |
      |                       |                       |  {\bf u}/ 2`          |
      +-----------------------+-----------------------+-----------------------+
      | UTEMP                 | :math:`T`             | this is computed from |
      |                       |                       | the other quantities  |
      |                       |                       | using the EOS         |
      +-----------------------+-----------------------+-----------------------+
      | UFA                   | :math:`\rho A_1`      | the first advected    |
      |                       |                       | quantity              |
      +-----------------------+-----------------------+-----------------------+
      | UFS                   | :math:`\rho X_1`      | the first species     |
      +-----------------------+-----------------------+-----------------------+
      | UFX                   | :math:`\rho Y_1`      | the first auxiliary   |
      |                       |                       | variable              |
      +-----------------------+-----------------------+-----------------------+
      | USHK                  | a shock flag          | (used for shock       |
      |                       |                       | detection)            |
      +-----------------------+-----------------------+-----------------------+
      | UMR                   | radial momentum       | (if is defined)       |
      +-----------------------+-----------------------+-----------------------+
      | UML                   | angular momentum      | (if is defined)       |
      +-----------------------+-----------------------+-----------------------+
      | UMP                   | vertical momentum     | (if is defined)       |
      +-----------------------+-----------------------+-----------------------+

-  primitive variable state: these arrays generally simply called
   q, and has components. Note: if
   is defined, then there are
   components that are pure hydro out of the total NQ components,
   and the pure hydro components always come first in the state array.

   Table \ `[table:primlist] <#table:primlist>`__ gives the names of the primitive variable integer
   keys for accessing these arrays.

   .. raw:: latex

      \centering

   .. table:: [table:primlist] The integer variable keys for
   accessing the primitive state vector. Note, unless otherwise
   specified the quantities without a subscript are “gas” only
   and those with the “tot” subscript are “gas + radiation”.

      +-----------------------+-----------------------+-----------------------+
      | **variable**          | **quantity**          | **note**              |
      +=======================+=======================+=======================+
      | QRHO                  | :math:`\rho`          |                       |
      +-----------------------+-----------------------+-----------------------+
      | QU                    | :math:`u`             |                       |
      +-----------------------+-----------------------+-----------------------+
      | QV                    | :math:`v`             |                       |
      +-----------------------+-----------------------+-----------------------+
      | QW                    | :math:`w`             |                       |
      +-----------------------+-----------------------+-----------------------+
      | QPRES                 | :math:`p`             |                       |
      +-----------------------+-----------------------+-----------------------+
      | QREINT                | :math:`\rho e`        |                       |
      +-----------------------+-----------------------+-----------------------+
      | QTEMP                 | :math:`T`             |                       |
      +-----------------------+-----------------------+-----------------------+
      | QGAME                 | :math:`p/(\rho e) + 1 |                       |
      |                       | `                     |                       |
      +-----------------------+-----------------------+-----------------------+
      | QFA                   | :math:`A_1`           | the first advected    |
      |                       |                       | quantity              |
      +-----------------------+-----------------------+-----------------------+
      | QFS                   | :math:`X_1`           | the first species     |
      +-----------------------+-----------------------+-----------------------+
      | QFX                   | :math:`Y_1`           | the first auxiliary   |
      |                       |                       | variable              |
      +-----------------------+-----------------------+-----------------------+
      | QPTOT                 | :math:`p_\mathrm{tot} | the total pressure,   |
      |                       | `                     | gas + radiation       |
      +-----------------------+-----------------------+-----------------------+
      | QREITOT               | :math:`e_\mathrm{tot} | the total specific    |
      |                       | `                     | internal energy, gas  |
      |                       |                       | + radiation           |
      +-----------------------+-----------------------+-----------------------+
      | QRAD                  | :math:`E_r`           | the radiation energy  |
      |                       |                       | (there are of these)  |
      +-----------------------+-----------------------+-----------------------+

-  auxiliary primitive variables: these arrays are generally called
   . The main difference between these and the regular
   primitive variables is that we do not attempt to do any
   reconstruction on their profiles. There are quantities, indexed
   by the integer keys listed in table \ `[table:qauxlist] <#table:qauxlist>`__.

   .. raw:: latex

      \centering

   .. table:: [table:qauxlist] The integer variable keys for
   accessing the auxiliary primitive state vector, quax.
   Note, unless otherwise specified the quantities without a
   subscript are “gas” only and those with the “tot” subscript
   are “gas + radiation”.

      +-----------------------+-----------------------+-----------------------+
      | **variable**          | **quantity**          | **note**              |
      +=======================+=======================+=======================+
      | QGAMC                 | :math:`\gamma_1`      | the first adiabatic   |
      |                       |                       | exponent, as returned |
      |                       |                       | from the EOS          |
      +-----------------------+-----------------------+-----------------------+
      | QC                    | :math:`c_s`           | the sound speed, as   |
      |                       |                       | returned from the EOS |
      +-----------------------+-----------------------+-----------------------+
      | QCSML                 |                       | a small sound speed   |
      |                       |                       | used for cutoffs      |
      +-----------------------+-----------------------+-----------------------+
      | QDPDR                 | :math:`\partial p/\pa | computed via the EOS  |
      |                       | rtial \rho |_e`       |                       |
      +-----------------------+-----------------------+-----------------------+
      | QDPDE                 | :math:`\partial p/\pa | computed via the EOS  |
      |                       | rtial e|_\rho`        |                       |
      +-----------------------+-----------------------+-----------------------+
      | QGAMCG                | :math:`{\Gamma_1}_\ma | includes radiation    |
      |                       | thrm{tot}`            | components (defined   |
      |                       |                       | only if is defined)   |
      +-----------------------+-----------------------+-----------------------+
      | QCG                   | :math:`{c_s}_\mathrm{ | total sound speed     |
      |                       | tot}`                 | including radiation   |
      |                       |                       | (defined only if is   |
      |                       |                       | defined)              |
      +-----------------------+-----------------------+-----------------------+
      | QLAMS                 | :math:`\lambda_f`     | the flux limiters     |
      |                       |                       | (defined only if is   |
      |                       |                       | defined)              |
      +-----------------------+-----------------------+-----------------------+

-  interface variable: these are the time-centered interface states
   returned by the Riemann solver. They are used to discretize some
   non-conservative terms in the equations. These arrays are generally
   called qx, qy, and qz for the x, y, and z
   interfaces respectively (in some places the numbers 1, 2, and 3 are
   used instead). There are components accessed with
   the integer keys defined in table \ `[table:gdlist] <#table:gdlist>`__

   .. raw:: latex

      \centering

   .. table:: [table:gdlist] The integer variable keys for
   accessing the Godunov interface state vectors.
   Note, unless otherwise specified the quantities without a
   subscript are “gas” only and those with the “tot” subscript
   are “gas + radiation”.

      +-----------------------+-----------------------+-----------------------+
      | **variable**          | **quantity**          | **note**              |
      +=======================+=======================+=======================+
      | QGDRHO                | :math:`\rho`          |                       |
      +-----------------------+-----------------------+-----------------------+
      | QDU                   | :math:`u`             |                       |
      +-----------------------+-----------------------+-----------------------+
      | QDV                   | :math:`v`             |                       |
      +-----------------------+-----------------------+-----------------------+
      | QDW                   | :math:`w`             |                       |
      +-----------------------+-----------------------+-----------------------+
      | QDPRES                | :math:`p`             | regardless of whether |
      |                       |                       | is defined, this is   |
      |                       |                       | always just the gas   |
      |                       |                       | pressure              |
      +-----------------------+-----------------------+-----------------------+
      | QDGAME                | :math:`\gamma_e = p/( | regardless of whether |
      |                       | \rho e) + 1`          | is defined, this is   |
      |                       |                       | always just the gas   |
      |                       |                       | contribution          |
      +-----------------------+-----------------------+-----------------------+
      | QDLAMS                | :math:`{\lambda_f}`   | the starting index    |
      |                       |                       | for the flux          |
      |                       |                       | limiter—there are     |
      |                       |                       | components (defined   |
      |                       |                       | only if is defined)   |
      +-----------------------+-----------------------+-----------------------+
      | QDERADS               | :math:`E_r`           | the starting index    |
      |                       |                       | for the radiation     |
      |                       |                       | energy—there are      |
      |                       |                       | components (defined   |
      |                       |                       | only if is defined)   |
      +-----------------------+-----------------------+-----------------------+

Conservation Forms
------------------

We begin with the fully compressible equations for the conserved state vector,
:math:`{\bf U}= (\rho, \rho {\bf u}, \rho E, \rho A_k, \rho X_k, \rho Y_k):`

.. math::

   \begin{aligned}
   \frac{\partial \rho}{\partial t} &=& - \nabla \cdot (\rho {\bf u}) + S_{{\rm ext},\rho}, \\
   \frac{\partial (\rho {\bf u})}{\partial t} &=& - \nabla \cdot (\rho {\bf u}{\bf u}) - \nabla p +\rho {\bf g}+ {\bf S}_{{\rm ext},\rho{\bf u}}, \\
   \frac{\partial (\rho E)}{\partial t} &=& - \nabla \cdot (\rho {\bf u}E + p {\bf u}) + \rho {\bf u}\cdot {\bf g}- \sum_k {\rho q_k \dot\omega_k} + \nabla\cdot{k_\mathrm{th}}\nabla T + S_{{\rm ext},\rho E}, \\
   \frac{\partial (\rho A_k)}{\partial t} &=& - \nabla \cdot (\rho {\bf u}A_k) + S_{{\rm ext},\rho A_k}, \\
   \frac{\partial (\rho X_k)}{\partial t} &=& - \nabla \cdot (\rho {\bf u}X_k) + \rho \dot\omega_k + S_{{\rm ext},\rho X_k}, \\
   \frac{\partial (\rho Y_k)}{\partial t} &=& - \nabla \cdot (\rho {\bf u}Y_k) + S_{{\rm ext},\rho Y_k}.\label{eq:compressible-equations}\end{aligned}

Here :math:`\rho, {\bf u}, T, p`, and :math:`{k_\mathrm{th}}` are the density, velocity,
temperature, pressure, and thermal conductivity, respectively, and :math:`E
= e + {\bf u}\cdot {\bf u}/ 2` is the total energy with :math:`e` representing the
internal energy. In addition, :math:`X_k` is the abundance of the :math:`k^{\rm
  th}` isotope, with associated production rate, :math:`\dot\omega_k`, and
energy release, :math:`q_k`. Here :math:`{\bf g}` is the gravitational vector, and
:math:`S_{{\rm ext},\rho}, {\bf S}_{{\rm ext}\rho{\bf u}}`, etc., are user-specified
source terms. :math:`A_k` is an advected quantity, i.e., a tracer. We also
carry around auxiliary variables, :math:`Y_k`, which have a user-defined
evolution equation, but by default are treated as advected quantities.

In the code we also carry around :math:`T` and :math:`\rho e` in the conservative
state vector even though they are derived from the other conserved
quantities. The ordering of the elements within :math:`{\bf U}` is defined
by integer variables into the array—see
Table \ `[table:consints] <#table:consints>`__

Some notes:

-  Regardless of the dimensionality of the problem, we always carry
   all 3 components of the velocity. This allows for, e.g., 2.5-d
   rotation (advecting the component of velocity out of the plane in
   axisymmetric coordinates).

   You should always initialize all velocity components to zero, and
   always construct the kinetic energy with all three velocity components.

-  There are advected quantities, which range from
   UFA: UFA+nadv-1. The advected quantities have no effect at all on
   the rest of the solution but can be useful as tracer quantities.

-  There are species (defined in the network
   directory), which range from UFS: UFS+nspec-1.

-  There are auxiliary variables, from
   UFX:UFX+naux-1 The auxiliary variables are passed into the equation
   of state routines along with the species; An example of an auxiliary
   variable is the electron fraction, :math:`Y_e`, in core collapse simulations.

.. raw:: latex

   \marginpar{\vskip-\baselineskip\raggedright\tiny\sffamily
   \hrule\smallskip{\color{red}note about qpass\_map here}\par\smallskip\hrule}

Source Terms
------------

We now compute explicit source terms for each variable in :math:`{\bf Q}` and
:math:`{\bf U}`. The primitive variable source terms will be used to construct
time-centered fluxes. The conserved variable source will be used to
advance the solution. We neglect reaction source terms since they are
accounted for in **Steps 1** and **6**. The source terms are:

.. math::

   {\bf S}_{{\bf Q}}^n =
   \left(\begin{array}{c}
   S_\rho \\
   {\bf S}_{{\bf u}} \\
   S_p \\
   S_{\rho e} \\
   S_{A_k} \\
   S_{X_k} \\
   S_{Y_k}
   \end{array}\right)^n
   =
   \left(\begin{array}{c}
   S_{{\rm ext},\rho} \\
   {\bf g}+ \frac{1}{\rho}{\bf S}_{{\rm ext},\rho{\bf u}} \\
   \frac{1}{\rho}\frac{\partial p}{\partial e}S_{{\rm ext},\rho E} + \frac{\partial p}{\partial\rho}S_{{\rm ext}\rho} \\
   \nabla\cdot{k_\mathrm{th}}\nabla T + S_{{\rm ext},\rho E} \\
   \frac{1}{\rho}S_{{\rm ext},\rho A_k} \\
   \frac{1}{\rho}S_{{\rm ext},\rho X_k} \\
   \frac{1}{\rho}S_{{\rm ext},\rho Y_k}
   \end{array}\right)^n,

.. math::

   {\bf S}_{{\bf U}}^n =
   \left(\begin{array}{c}
   {\bf S}_{\rho{\bf u}} \\
   S_{\rho E} \\
   S_{\rho A_k} \\
   S_{\rho X_k} \\
   S_{\rho Y_k}
   \end{array}\right)^n
   =
   \left(\begin{array}{c}
   \rho {\bf g}+ {\bf S}_{{\rm ext},\rho{\bf u}} \\
   \rho {\bf u}\cdot {\bf g}+ \nabla\cdot{k_\mathrm{th}}\nabla T + S_{{\rm ext},\rho E} \\
   S_{{\rm ext},\rho A_k} \\
   S_{{\rm ext},\rho X_k} \\
   S_{{\rm ext},\rho Y_k}
   \end{array}\right)^n.

Primitive Forms
---------------

Castro uses the primitive form of the fluid equations, defined in terms of
the state :math:`{\bf Q}= (\rho, {\bf u}, p, \rho e, A_k, X_k, Y_k)`, to construct the
interface states that are input to the Riemann problem.

The primitive variable equations for density, velocity, and pressure are:

.. math::

   \begin{aligned}
     \frac{\partial\rho}{\partial t} &=& -{\bf u}\cdot\nabla\rho - \rho\nabla\cdot{\bf u}+ S_{{\rm ext},\rho} \\
   %
     \frac{\partial{\bf u}}{\partial t} &=& -{\bf u}\cdot\nabla{\bf u}- \frac{1}{\rho}\nabla p + {\bf g}+ 
   \frac{1}{\rho} ({\bf S}_{{\rm ext},\rho{\bf u}} - {\bf u}\; S_{{\rm ext},\rho}) \\
   \frac{\partial p}{\partial t} &=& -{\bf u}\cdot\nabla p - \rho c^2\nabla\cdot{\bf u}+
   \left(\frac{\partial p}{\partial \rho}\right)_{e,X}S_{{\rm ext},\rho}\nonumber\\
   &&+\  \frac{1}{\rho}\sum_k\left(\frac{\partial p}{\partial X_k}\right)_{\rho,e,X_j,j\neq k}\left(\rho\dot\omega_k + S_{{\rm ext},\rho X_k} - X_kS_{{\rm ext},\rho}\right)\nonumber\\
   && +\  \frac{1}{\rho}\left(\frac{\partial p}{\partial e}\right)_{\rho,X}\left[-eS_{{\rm ext},\rho} - \sum_k\rho q_k\dot\omega_k + \nabla\cdot{k_\mathrm{th}}\nabla T \right.\nonumber\\
   && \quad\qquad\qquad\qquad+\ S_{{\rm ext},\rho E} - {\bf u}\cdot\left({\bf S}_{{\rm ext},\rho{\bf u}} - \frac{{\bf u}}{2}S_{{\rm ext},\rho}\right)\Biggr] \end{aligned}

The advected quantities appear as:

.. math::

   \begin{aligned}
   \frac{\partial A_k}{\partial t} &=& -{\bf u}\cdot\nabla A_k + \frac{1}{\rho}
                                        ( S_{{\rm ext},\rho A_k} - A_k S_{{\rm ext},\rho} ), \\
   \frac{\partial X_k}{\partial t} &=& -{\bf u}\cdot\nabla X_k + \dot\omega_k + \frac{1}{\rho}
                                        ( S_{{\rm ext},\rho X_k}  - X_k S_{{\rm ext},\rho} ), \\
   \frac{\partial Y_k}{\partial t} &=& -{\bf u}\cdot\nabla Y_k + \frac{1}{\rho} 
                                        ( S_{{\rm ext},\rho Y_k}  - Y_k S_{{\rm ext},\rho} ).\end{aligned}

All of the primitive variables are derived from the conservative state
vector, as described in Section `6.1 <#Sec:Compute Primitive Variables>`__.
When accessing the primitive variable state vector, the integer variable
keys for the different quantities are listed in Table \ `[table:primlist] <#table:primlist>`__.

Internal energy and temperature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We augment the above system with an internal energy equation:

.. math::

   \begin{aligned}
   \frac{\partial(\rho e)}{\partial t} &=& - {\bf u}\cdot\nabla(\rho e) - (\rho e+p)\nabla\cdot{\bf u}- \sum_k \rho q_k\dot\omega_k 
                                           + \nabla\cdot{k_\mathrm{th}}\nabla T + S_{{\rm ext},\rho E} \nonumber\\
   && -\  {\bf u}\cdot\left({\bf S}_{{\rm ext},\rho{\bf u}}-\frac{1}{2}S_{{\rm ext},\rho}{\bf u}\right), \end{aligned}

This has two benefits. First, for a general equation of state,
carrying around an additional thermodynamic quantity allows us to
avoid equation of state calls (in particular, in the Riemann solver,
see e.g. :raw-latex:`\cite{colglaz}`). Second, it is sometimes the case that the
internal energy calculated as

.. math:: e_T \equiv E - \frac{1}{2} \mathbf{v}^2

is
unreliable. This has two usual causes: one, for high Mach number
flows, the kinetic energy can dominate the total gas energy, making
the subtraction numerically unreliable; two, if you use gravity or
other source terms, these can indirectly alter the value of the
internal energy if obtained from the total energy.

To provide a more reasonable internal energy for defining the
thermodynamic state, we have implemented the dual energy formalism
from ENZO :raw-latex:`\cite{bryan:1995,bryan:2014}`, where we switch between :math:`(\rho
e)` and :math:`(\rho e_T)` depending on the local state of the fluid. To do
so, we define parameters :math:`\eta_1`, :math:`\eta_2`, and :math:`\eta_3`,
corresponding to the code parameters
,
, and
. We then consider the ratio :math:`e_T
/ E`, the ratio of the internal energy (derived from the total energy)
to the total energy. These parameters are used as follows:

-  :math:`\eta_1`: If :math:`e_T > \eta_1 E`, then we use :math:`e_T` for the purpose
   of calculating the pressure in the hydrodynamics update. Otherwise,
   we use the :math:`e` from the internal energy equation in our EOS call to
   get the pressure.

-  :math:`\eta_2`: At the end of each hydro advance, we examine whether
   :math:`e_T > \eta_2 E`. If so, we reset :math:`e` to be equal to :math:`e_T`,
   discarding the results of the internal energy equation. Otherwise,
   we keep :math:`e` as it is.

   Optionally we can also update :math:`E` so that it gains the difference of
   the old and and new :math:`e`, by setting
   to 1.

-  :math:`\eta_3`: Similar to :math:`\eta_1`, if :math:`e_T > \eta_3 E`, we use
   :math:`e_T` for the purposes of our nuclear reactions, otherwise, we use
   :math:`e`.

Note that our version of the internal energy equation does not require
an artificial viscosity, as used in some other hydrodynamics
codes. The update for :math:`(\rho e)` uses information from the Riemann
solve to calculate the fluxes, which contains the information
intrinsic to the shock-capturing part of the scheme.

In the code we also carry around :math:`T` in the primitive state vector.

Primitive Variable System
~~~~~~~~~~~~~~~~~~~~~~~~~

The full primitive variable form (without the advected or auxiliary
quantities) is

.. math:: \frac{\partial{\bf Q}}{\partial t} + \sum_d {\bf A}_d\frac{\partial{\bf Q}}{\partial x_d} = {\bf S}_{{\bf Q}}.

For example, in 2D:

.. math::

   \left(\begin{array}{c}
   \rho \\
   u \\
   v \\
   p \\
   \rho e \\
   X_k
   \end{array}\right)_t
   +
   \left(\begin{array}{cccccc}
   u & \rho & 0 & 0 & 0 & 0 \\
   0 & u & 0 & \frac{1}{\rho} & 0 & 0 \\
   0 & 0 & u & 0 & 0 & 0 \\
   0 & \rho c^2 & 0 & u & 0 & 0 \\
   0 & \rho e + p & 0 & 0 & u & 0 \\
   0 & 0 & 0 & 0 & 0 & u
   \end{array}\right)
   \left(\begin{array}{c}
   \rho \\
   u \\
   v \\
   p \\
   \rho e \\
   X_k
   \end{array}\right)_x
   +
   \left(\begin{array}{cccccc}
   v & 0 & \rho & 0 & 0 & 0 \\
   0 & v & 0 & 0 & 0 & 0 \\
   0 & 0 & v & \frac{1}{\rho} & 0 & 0 \\
   0 & 0 & \rho c^2 & v & 0 & 0 \\
   0 & 0 & \rho e + p & 0 & v & 0 \\
   0 & 0 & 0 & 0 & 0 & v
   \end{array}\right)
   \left(\begin{array}{c}
   \rho \\
   u \\
   v \\
   p \\
   \rho e \\
   X_k
   \end{array}\right)_y
   =
   {\bf S}_{\bf Q}

The eigenvalues are:

.. math:: {\bf \Lambda}({\bf A}_x) = \{u-c,u,u,u,u,u+c\}, \qquad {\bf \Lambda}({\bf A}_y) = \{v-c,v,v,v,v,v+c\} .

The right column eigenvectors are:

.. math::

   {\bf R}({\bf A}_x) =
   \left(\begin{array}{cccccc}
   1 & 1 & 0 & 0 & 0 & 1 \\
   -\frac{c}{\rho} & 0 & 0 & 0 & 0 & \frac{c}{\rho} \\
   0 & 0 & 1 & 0 & 0 & 0 \\
   c^2 & 0 & 0 & 0 & 0 & c^2 \\
   h & 0 & 0 & 1 & 0 & h \\
   0 & 0 & 0 & 0 & 1 & 0 \\
   \end{array}\right),
   \qquad
   {\bf R}({\bf A}_y) =
   \left(\begin{array}{cccccc}
   1 & 1 & 0 & 0 & 0 & 1 \\
   0 & 0 & 1 & 0 & 0 & 0 \\
   -\frac{c}{\rho} & 0 & 0 & 0 & 0 & \frac{c}{\rho} \\
   c^2 & 0 & 0 & 0 & 0 & c^2 \\
   h & 0 & 0 & 1 & 0 & h \\
   0 & 0 & 0 & 0 & 1 & 0 \\
   \end{array}\right).

The left row eigenvectors, normalized so that :math:`{\bf R}_d\cdot{\bf L}_d = {\bf I}` are:

.. math::

   {\bf L}_x =
   \left(\begin{array}{cccccc}
   0 & -\frac{\rho}{2c} & 0 & \frac{1}{2c^2} & 0 & 0 \\
   1 & 0 & 0 & -\frac{1}{c^2} & 0 & 0 \\
   0 & 0 & 1 & 0 & 0 & 0 \\
   0 & 0 & 0 & -\frac{h}{c^2} & 1 & 0 \\
   0 & 0 & 0 & 0 & 0 & 1 \\
   0 & \frac{\rho}{2c} & 0 & \frac{1}{2c^2} & 0 & 0
   \end{array}\right),
   \qquad
   {\bf L}_y =
   \left(\begin{array}{cccccc}
   0 & 0 & -\frac{\rho}{2c} & \frac{1}{2c^2} & 0 & 0 \\
   1 & 0 & 0 & -\frac{1}{c^2} & 0 & 0 \\
   0 & 1 & 0 & 0 & 0 & 0 \\
   0 & 0 & 0 & -\frac{h}{c^2} & 1 & 0 \\
   0 & 0 & 0 & 0 & 0 & 1 \\
   0 & 0 & \frac{\rho}{2c} & \frac{1}{2c^2} & 0 & 0
   \end{array}\right).

.. _Sec:Advection Step:

Hydrodynamics Update
--------------------

There are four major steps in the hydrodynamics update:

#. Converting to primitive variables

#. Construction the edge states

#. Solving the Riemann problem

#. Doing the conservative update

Each of these steps has a variety of runtime parameters that
affect their behavior. Additionally, there are some general
runtime parameters for hydrodynamics:

-  : time-advance the fluid dynamical
   equations (0 or 1; must be set)

-  : include additional user-specified
   source term (0 or 1; default 0)

-  : call the sponge routine
   after the solution update (0 or 1; default: 0)

   The purpose of the sponge is to damp velocities outside of a star, to
   prevent them from dominating the timestep constraint. The sponge parameters
   are set in your probin file, in the &sponge namelist. You can sponge either
   on radius from the center (using and
   ) or on density (using
   and ). The timescale of the damping is
   set through .

-  : enforce that :math:`\sum_i X_i = 1`
   (0 or 1; default: 0)

-  : enforce constant mass flux at
   domain boundary (0 or 1; default: 1)

-  : is internal energy allowed to be
   negative? (0 or 1; default: 1)

-  : this is used to set the boundary
   conditions by assuming the star is spherically symmetric in
   the outer regions (0 or 1; default: 0)

   When used, Castro averages the values at a given radius over the
   cells that are inside the domain to define a radial function. This
   function is then used to set the values outside the domain in
   implementing the boundary conditions.

-  : (0 or 1; default: 0)

Several floors are imposed on the thermodynamic quantities to prevet unphysical
behavior:

-  : (Real; default: -1.e20)

-  : (Real; default: -1.e20)

-  : (Real; default: -1.e20)

.. _Sec:Compute Primitive Variables:

Compute Primitive Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We compute the primtive variables from the conserved variables.

-  :math:`\rho, \rho e`: directly copy these from the conserved state
   vector

-  :math:`{\bf u}, A_k, X_k, Y_k`: copy these from the conserved state
   vector, dividing by :math:`\rho`

-  :math:`p,T`: use the EOS.

   First, if is 0 (it defaults to
   1) and :math:`e < 0`, we do the following:

   #. Use the EOS to set :math:`e = e(\rho,T_{\rm small},X_k)`.

   #. If :math:`e < 0`, abort the program with an error message.

   Now, use the EOS to compute :math:`p,T = p,T(\rho,e,X_k)`.

We also compute the flattening coefficient, :math:`\chi\in[0,1]`, used in
the edge state prediction to further limit slopes near strong shocks.
We use the same flattening procedure described in the the the original
PPM paper :raw-latex:`\cite{ppm}` and the Flash paper :raw-latex:`\cite{flash}`.
A flattening coefficient of 1 indicates that no additional limiting
takes place; a flattening coefficient of 0 means we effectively drop
order to a first-order Godunov scheme (this convention is opposite of
that used in the Flash paper). For each cell, we compute the
flattening coefficient for each spatial direction, and choose the
minimum value over all directions. As an example, to compute the
flattening for the x-direction, here are the steps:

#. Define :math:`\zeta`

   .. math:: \zeta_i = \frac{p_{i+1}-p_{i-1}}{\max\left(p_{\rm small},|p_{i+2}-p_{i-2}|\right)}.

#. Define :math:`\tilde\chi`

   .. math:: \tilde\chi_i = \min\left\{1,\max[0,a(\zeta_i - b)]\right\},

   where :math:`a=10` and :math:`b=0.75` are tunable parameters. We are essentially
   setting :math:`\tilde\chi_i=a(\zeta_i-b)`, and then constraining
   :math:`\tilde\chi_i` to lie in the range :math:`[0,1]`. Then, if either
   :math:`u_{i+1}-u_{i-1}<0` or

   .. math:: \frac{p_{i+1}-p_{i-1}}{\min(p_{i+1},p_{i-1})} \le c,

   where :math:`c=1/3` is a tunable parameter, then set :math:`\tilde\chi_i=0`.

#. Define :math:`\chi`

   .. math::

      \chi_i =
      \begin{cases}
      1 - \max(\tilde\chi_i,\tilde\chi_{i-1}) & p_{i+1}-p_{i-1} > 0 \\
      1 - \max(\tilde\chi_i,\tilde\chi_{i+1}) & \text{otherwise}
      \end{cases}.

The following runtime parameters affect the behavior here:

-  turns on/off the flattening of parabola
   near shocks (0 or 1; default 1)

Edge State Prediction
~~~~~~~~~~~~~~~~~~~~~

We wish to compute a left and right state of primitive variables at
each edge to be used as inputs to the Riemann problem. There
are several reconstruction techniques, a piecewise
linear method that follows the description in Colella (1990) :raw-latex:`\cite{colella:1990}`,
the classic PPM limiters :raw-latex:`\cite{ppm}`, and the new PPM limiters introduced
in Colella & Sekora (2008) :raw-latex:`\cite{colellasekora}`. The choice of
limiters is determined by .

For the new PPM limiters, we have further modified the method
of :raw-latex:`\cite{colellasekora}` to eliminate sensitivity due to roundoff error
(modifications via personal communication with Colella).

We also use characteristic tracing with corner coupling in 3D, as
described in Miller & Colella (2002) :raw-latex:`\cite{millercolella:2002}`. We
give full details of the new PPM algorithm, as it has not appeared before
in the literature, and summarize the developments from Miller &
Colella.

The PPM algorithm is used to compute time-centered edge states by
extrapolating the base-time data in space and time. The edge states
are dual-valued, i.e., at each face, there is a left state and a right
state estimate. The spatial extrapolation is one-dimensional, i.e.,
transverse derivatives are ignored. We also use a flattening
procedure to further limit the edge state values. The Miller &
Colella algorithm, which we describe later, incorporates the
transverse terms, and also describes the modifications required for
equations with additional characteristics besides the fluid velocity.
There are four steps to compute these dual-valued edge states (here,
we use :math:`s` to denote an arbitrary scalar from :math:`{\bf Q}`, and we write the
equations in 1D, for simplicity):

-  **Step 1**: Compute :math:`s_{i,+}` and :math:`s_{i,-}`, which are spatial
   interpolations of :math:`s` to the hi and lo side of the face with special
   limiters, respectively. Begin by interpolating :math:`s` to edges using a
   4th-order interpolation in space:

   .. math::

      s_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
         \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
         \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
         \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}} = \frac{7}{12}\left(s_{i+1}+s_i\right) - \frac{1}{12}\left(s_{i+2}+s_{i-1}\right).

   Then, if :math:`(s_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
      \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
      \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
      \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}}-s_i)(s_{i+1}-s_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
      \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
      \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
      \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}}) < 0`, we limit
   :math:`s_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
      \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
      \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
      \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}}` a nonlinear combination of approximations to the
   second derivative. The steps are as follows:

   #. Define:

      .. math::

         \begin{aligned}
         (D^2s)_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}} &=& 3\left(s_{i}-2s_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}}+s_{i+1}\right) \\
         (D^2s)_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2},L} &=& s_{i-1}-2s_{i}+s_{i+1} \\
         (D^2s)_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2},R} &=& s_{i}-2s_{i+1}+s_{i+2}\end{aligned}

   #. Define

      .. math::

         s = \text{sign}\left[(D^2s)_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}}\right],

      .. math::

         (D^2s)_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2},\text{lim}} = s\max\left\{\min\left[Cs\left|(D^2s)_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2},L}\right|,Cs\left|(D^2s)_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2},R}\right|,s\left|(D^2s)_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}}\right|\right],0\right\},

      where :math:`C=1.25` as used in Colella and Sekora 2009. The limited value
      of :math:`s_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
         \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
         \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
         \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}}` is

      .. math::

         s_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}} = \frac{1}{2}\left(s_{i}+s_{i+1}\right) - \frac{1}{6}(D^2s)_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2},\text{lim}}.

   Now we implement an updated implementation of the Colella & Sekora
   algorithm which eliminates sensitivity to roundoff. First we
   need to detect whether a particular cell corresponds to an
   “extremum”. There are two tests.

   -  For the first test, define

      .. math::

         \alpha_{i,\pm} = s_{i\pm\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}} - s_i.

      If :math:`\alpha_{i,+}\alpha_{i,-} \ge 0`, then we are at an extremum.

   -  We only apply the second test if either :math:`|\alpha_{i,\pm}| >
        2|\alpha_{i,\mp}|`. If so, we define:

      .. math::

         \begin{aligned}
         (Ds)_{i,{\rm face},-} &=& s_{i-\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}} - s_{i-\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 3}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 3}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 3}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{3\!/2}} \\
         (Ds)_{i,{\rm face},+} &=& s_{i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 3}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 3}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 3}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{3\!/2}} - s_{i-\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}}\end{aligned}

      .. math:: (Ds)_{i,{\rm face,min}} = \min\left[\left|(Ds)_{i,{\rm face},-}\right|,\left|(Ds)_{i,{\rm face},+}\right|\right].

      .. math::

         \begin{aligned}
         (Ds)_{i,{\rm cc},-} &=& s_{i} - s_{i-1} \\
         (Ds)_{i,{\rm cc},+} &=& s_{i+1} - s_{i}\end{aligned}

      .. math:: (Ds)_{i,{\rm cc,min}} = \min\left[\left|(Ds)_{i,{\rm cc},-}\right|,\left|(Ds)_{i,{\rm cc},+}\right|\right].

      If :math:`(Ds)_{i,{\rm face,min}} \ge (Ds)_{i,{\rm cc,min}}`, set
      :math:`(Ds)_{i,\pm} = (Ds)_{i,{\rm face},\pm}`. Otherwise, set
      :math:`(Ds)_{i,\pm} = (Ds)_{i,{\rm cc},\pm}`. Finally, we are at an extreumum if
      :math:`(Ds)_{i,+}(Ds)_{i,-} \le 0`.

   Thus concludes the extremum tests. The remaining limiters depend on
   whether we are at an extremum.

   -  If we are at an extremum, we modify :math:`\alpha_{i,\pm}`. First, we
      define

      .. math::

         \begin{aligned}
         (D^2s)_{i} &=& 6(\alpha_{i,+}+\alpha_{i,-}) \\
         (D^2s)_{i,L} &=& s_{i-2}-2s_{i-1}+s_{i} \\
         (D^2s)_{i,R} &=& s_{i}-2s_{i+1}+s_{i+2} \\
         (D^2s)_{i,C} &=& s_{i-1}-2s_{i}+s_{i+1}\end{aligned}

      Then, define

      .. math:: s = \text{sign}\left[(D^2s)_{i}\right],

      .. math:: (D^2s)_{i,\text{lim}} = \max\left\{\min\left[s(D^2s)_{i},Cs\left|(D^2s)_{i,L}\right|,Cs\left|(D^2s)_{i,R}\right|,Cs\left|(D^2s)_{i,C}\right|\right],0\right\}.

      Then,

      .. math:: \alpha_{i,\pm} = \frac{\alpha_{i,\pm}(D^2s)_{i,\text{lim}}}{\max\left[(D^2s)_{i},1\times 10^{-10}\right]}

   -  If we are not at an extremum and :math:`|\alpha_{i,\pm}| >
        2|\alpha_{i,\mp}|`, then define

      .. math:: s = \text{sign}(\alpha_{i,\mp})

      .. math:: \delta\mathcal{I}_{\text{ext}} = \frac{-\alpha_{i,\pm}^2}{4\left(\alpha_{j,+}+\alpha_{j,-}\right)},

      .. math:: \delta s = s_{i\mp 1} - s_i,

      If :math:`s\delta\mathcal{I}_{\text{ext}} \ge s\delta s`, then we perform
      the following test. If :math:`s\delta s - \alpha_{i,\mp} \ge 1\times
      10^{-10}`, then

      .. math::

         \alpha_{i,\pm} =  -2\delta s - 2s\left[(\delta s)^2 - \delta s \alpha_{i,\mp}\right]^{\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
            \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
            \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}}

      otherwise,

      .. math:: \alpha_{i,\pm} =  -2\alpha_{i,\mp}

   Finally, :math:`s_{i,\pm} = s_i + \alpha_{i,\pm}`.

-  **Step 2**: Construct a quadratic profile using :math:`s_{i,-},s_i`,
   and :math:`s_{i,+}`.

   .. math:: s_i^I(x) = s_{i,-} + \xi\left[s_{i,+} - s_{i,-} + s_{6,i}(1-\xi)\right],\label{Quadratic Interp}

   .. math:: s_6 = 6s_{i} - 3\left(s_{i,-}+s_{i,+}\right),

   .. math:: \xi = \frac{x - ih}{h}, ~ 0 \le \xi \le 1.

-  | **Step 3:** Integrate quadratic profiles. We are essentially
     computing the average value swept out by the quadratic profile
     across the face assuming the profile is moving at a speed
     :math:`\lambda_k`.
   | Define the following integrals, where :math:`\sigma_k =
       |\lambda_k|\Delta t/h`:

     .. math::

        \begin{aligned}
        \mathcal{I}^{(k)}_{+}(s_i) &=& \frac{1}{\sigma_k h}\int_{(i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
           \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
           \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
           \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2})h-\sigma_k h}^{(i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
           \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
           \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
           \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2})h}s_i^I(x)dx \\
        \mathcal{I}^{(k)}_{-}(s_i) &=& \frac{1}{\sigma_k h}\int_{(i-\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
           \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
           \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
           \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2})h}^{(i-\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
           \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
           \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
           \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2})h+\sigma_k h}s_i^I(x)dx\end{aligned}

     Plugging in (`[Quadratic Interp] <#Quadratic Interp>`__) gives:

     .. math::

        \begin{aligned}
        \mathcal{I}^{(k)}_{+}(s_i) &=& s_{i,+} - \frac{\sigma_k}{2}\left[s_{i,+}-s_{i,-}-\left(1-\frac{2}{3}\sigma_k\right)s_{6,i}\right], \\
        \mathcal{I}^{(k)}_{-}(s_i) &=& s_{i,-} + \frac{\sigma_k}{2}\left[s_{i,+}-s_{i,-}+\left(1-\frac{2}{3}\sigma_k\right)s_{6,i}\right].\end{aligned}

-  **Step 4:** Obtain 1D edge states by performing a 1D
   extrapolation to get left and right edge states. Note that we
   include an explicit source term contribution.

   .. math::

      \begin{aligned}
      s_{L,i+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
         \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
         \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
         \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}} &=& s_i - \chi_i\sum_{k:\lambda_k \ge 0}{\bf l}_k\cdot\left[s_i-\mathcal{I}^{(k)}_{+}(s_i)\right]{\bf r}_k + \frac{\Delta t}{2}S_i^n, \\
      s_{R,i-\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
         \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
         \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
         \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}} &=& s_i - \chi_i\sum_{k:\lambda_k < 0}{\bf l}_k\cdot\left[s_i-\mathcal{I}^{(k)}_{-}(s_i)\right]{\bf r}_k + \frac{\Delta t}{2}S_i^n.\end{aligned}

   Here, :math:`{\bf r}_k` is the :math:`k^{\rm th}` right column eigenvector of
   :math:`{\bf R}({\bf A}_d)` and :math:`{\bf l}_k` is the :math:`k^{\rm th}` left row eigenvector lf
   :math:`{\bf L}({\bf A}_d)`. The flattening coefficient is :math:`\chi_i`.

In order to add the transverse terms in an spatial operator unsplit
framework, the details follow exactly as given in Section 4.2.1 in
Miller & Colella, except for the details of the Riemann solver,
which are given below.

For the reconstruction of the interface states, the following apply:

-  : use piecewise linear vs PPM algorithm
   (0, 1, 2; default: 1)

   Values of 1 and 2 are both piecewise parabolic reconstruction, with
   2 using updated limiters that better preserve extrema.

-  does various attempts to use the
   temperature in the reconstruction of the interface states. This
   is experimental.

-  reconstructs :math:`\gamma_e = p/(\rho e) + 1`
   to the interfaces and does the necessary transverse terms to aid in
   the conversion between the conserved and primitive interface states
   in the transverse flux routines (0 or 1; default 0)

-  uses the reference states in
   the evaluation of the eigenvectors for the characteristic projection
   (0 or 1; default 0)

The interface states are corrected with information from the
transverse directions to make this a second-order update. These
transverse directions involve separate Riemann solves. Sometimes, the
update to the interface state from the transverse directions can make
the state ill-posed. There are several parameters that help fix this:

-  : If this is 1, then we call
   the equation of state on the interface, using :math:`\rho`, :math:`e`, and
   :math:`X_k`, to get the interface pressure. This should result in a
   thermodynamically consistent interface state.

-  : If the transverse
   corrections result in a negative density on the interface, then we
   reset all of the interface states to their values before the
   transverse corrections.

-  : The transverse updates operate
   on the conserved state. Usually, we construct the interface
   :math:`(\rho e)` in the transverse update from total energy and the
   kinetic energy, however, if the interface :math:`(rho e)` is negative,
   and transverse_reset_rhoe = 1, then we explicitly
   discretize an equation for the evolution of :math:`(\rho e)`, including
   its transverse update.

Riemann Problem
~~~~~~~~~~~~~~~

Castro has three main options for the Riemann solver—the
Colella & Glaz solver :raw-latex:`\cite{colglaz}` (the same solver used
by Flash), a simpler solver described in an unpublished
manuscript by Colella, Glaz, & Ferguson, and an HLLC
solver. The first two are both
two-shock approximate solvers, but differ in how they approximate
the thermodynamics in the “star” region.

Inputs from the edge state prediction are :math:`\rho_{L/R}, u_{L/R},
v_{L/R}, p_{L/R}`, and :math:`(\rho e)_{L/R}` (:math:`v` represents all of the
transverse velocity components). We also compute :math:`\Gamma \equiv d\log
p / d\log \rho |_s` at cell centers and copy these to edges directly
to get the left and right states, :math:`\Gamma_{L/R}`. We also define
:math:`c_{\rm avg}` as a face-centered value that is the average of the
neighboring cell-centered values of :math:`c`. We have also computed
:math:`\rho_{\rm small}, p_{\rm small}`, and :math:`c_{\rm small}` using
cell-centered data.

Here are the steps. First, define :math:`(\rho c)_{\rm small} = \rho_{\rm
  small}c_{\rm small}`. Then, define:

.. math:: (\rho c)_{L/R} = \max\left[(\rho c)_{\rm small},\left|\Gamma_{L/R},p_{L/R},\rho_{L/R}\right|\right].

Define star states:

.. math:: p^* = \max\left[p_{\rm small},\frac{\left[(\rho c)_L p_R + (\rho c)_R p_L\right] + (\rho c)_L(\rho c)_R(u_L-u_R)}{(\rho c)_L + (\rho c)_R}\right],

.. math:: u^* = \frac{\left[(\rho c)_L u_L + (\rho c)_R u_R\right]+ (p_L - p_R)}{(\rho c)_L + (\rho c)_R}.

If :math:`u^* \ge 0` then define :math:`\rho_0, u_0, p_0, (\rho e)_0` and :math:`\Gamma_0` to be the left state. Otherwise, define them to be the right state. Then, set

.. math:: \rho_0 = \max(\rho_{\rm small},\rho_0),

and define

.. math:: c_0 = \max\left(c_{\rm small},\sqrt{\frac{\Gamma_0 p_0}{\rho_0}}\right),

.. math:: \rho^* = \rho_0 + \frac{p^* - p_0}{c_0^2},

.. math:: (\rho e)^* = (\rho e)_0 + (p^* - p_0)\frac{(\rho e)_0 + p_0}{\rho_0 c_0^2},

.. math:: c^* = \max\left(c_{\rm small},\sqrt{\left|\frac{\Gamma_0 p^*}{\rho^*}\right|}\right)

Then,

.. math::

   \begin{aligned}
   c_{\rm out} &=& c_0 - {\rm sign}(u^*)u_0, \\
   c_{\rm in} &=& c^* - {\rm sign}(u^*)u^*, \\
   c_{\rm shock} &=& \frac{c_{\rm in} + c_{\rm out}}{2}.\end{aligned}

If :math:`p^* - p_0 \ge 0`, then :math:`c_{\rm in} = c_{\rm out} = c_{\rm shock}`.
Then, if :math:`c_{\rm out} = c_{\rm in}`, we define :math:`c_{\rm temp} =
\epsilon c_{\rm avg}`. Otherwise, :math:`c_{\rm temp} = c_{\rm out} -
c_{\rm in}`. We define the fraction

.. math:: f = \frac{1}{2}\left[1 + \frac{c_{\rm out} + c_{\rm in}}{c_{\rm temp}}\right],

and constrain :math:`f` to lie in the range :math:`f\in[0,1]`.

To get the final “Godunov” state, for the transverse velocity, we
upwind based on :math:`u^*`.

.. math::

   v_{\rm gdnv} =
   \begin{cases}
   v_L, & u^* \ge 0 \\
   v_R, & {\rm otherwise}
   \end{cases}.

Then, define

.. math::

   \begin{aligned}
   \rho_{\rm gdnv} &=& f\rho^* + (1-f)\rho_0, \\
   u_{\rm gdnv} &=& f u^* + (1-f)u_0, \\
   p_{\rm gdnv} &=& f p^* + (1-f)p_0, \\
   (\rho e)_{\rm gdnv} &=& f(\rho e)^* + (1-f)(\rho e)_0.\end{aligned}

Finally, if :math:`c_{\rm out} < 0`, set :math:`\rho_{\rm gdnv}=\rho_0, u_{\rm
  gdnv}=u_0, p_{\rm gdnv}=p_0`, and :math:`(\rho e)_{\rm gdnv}=(\rho e)_0`.
If :math:`c_{\rm in}\ge 0`, set :math:`\rho_{\rm gdnv}=\rho^*, u_{\rm gdnv}=u^*,
p_{\rm gdnv}=p^*`, and :math:`(\rho e)_{\rm gdnv}=(\rho e)^*`.

If instead the Colella & Glaz solver is used, then we define

.. math:: \gamma \equiv \frac{p}{\rho e} + 1

on each side of the interface and follow the rest of the algorithm as
described in the original paper.

For the construction of the fluxes in the Riemann solver, the following
parameters apply:

-  : this can be one of the following values:

   -  0: the Colella, Glaz, & Ferguson solver.

   -  1: the Colella & Glaz solver

   -  2: the HLLC solver. Note: this should only be used with Cartesian
      geometries because it relies on the pressure term being part of the flux
      in the momentum equation.

   The default is to use the solver based on an unpublished Colella,
   Glaz, & Ferguson manuscript (it also appears in :raw-latex:`\cite{pember:1996}`),
   as described in the original Castro paper :raw-latex:`\cite{castro_I}`.

   The Colella & Glaz solver is iterative, and two runtime parameters are used
   to control its behavior:

   -  : number of iterations for CG algorithm
      (Integer; default: 12)

   -  : tolerance for CG solver when solving
      for the “star” state (Real; default: 1.0e-5)

   -  : this controls what happens if the root
      finding in the CG solver fails. There is a nonlinear equation to find
      the pressure in the *star* region from the jump conditions for a
      shock (this is the two-shock approximation—the left and right states
      are linked to the star region each by a shock). The default root
      finding algorithm is a secant method, but this can sometimes fail.

      The options here are:

      -  0 : do nothing. The pressure from each iteration is
         printed and the code aborts with a failure

      -  1 : revert to the original guess for p-star and carry
         through on the remainder of the Riemann solve. This is almost like
         dropping down to the CGF solver. The p-star used is very approximate.

      -  2 : switch to bisection and do an additional cg_maxiter
         iterations to find the root. Sometimes this can work where the
         secant method fails.

-  castro.hybrid_riemann: switch to an HLL Riemann solver when we are
   in a zone with a shock (0 or 1; default 0)

   This eliminates an odd-even decoupling issue (see the oddeven
   problem). Note, this cannot be used with the HLLC solver.

Compute Fluxes and Update
~~~~~~~~~~~~~~~~~~~~~~~~~

Compute the fluxes as a function of the primitive variables, and then
advance the solution:

.. math::

   {\bf U}^{n+1} = {\bf U}^n - \Delta t\nabla\cdot{\bf F}^{n+\mathchoice{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
      \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptfont 0 1}\kern-.15em/
      \kern-.15em\lower.25ex\hbox{\the\scriptfont 0 2}}{\kern 0em\raise.5ex\hbox{\the\scriptscriptfont 0 1}\kern-.2em/
      \kern-.15em\lower.25ex\hbox{\the\scriptscriptfont 0 2}}{1\!/2}}+ \Delta t{\bf S}^n.

Again, note that since the source term is not time centered, this is
not a second-order method. After the advective update, we correct the
solution, effectively time-centering the source term.

Temperature Fixes
-----------------

There are a number of experimental options for improving the behavior
of the temperature in the reconstruction and interface state
prediction. The options are controlled by ,
which takes values:

-  0: the default method—temperature is not considered

-  1: do parabolic reconstruction on :math:`T`, giving
   :math:`\mathcal{I}_{+}^{(k)}(T_i)`. We then derive the pressure and
   internal energy (gas portion) via the equation of state as:

   .. math::

      \begin{aligned}
            \mathcal{I}_{+}^{(k)}(p_i) &= p(\mathcal{I}_{+}^{(k)}(\rho_i), \mathcal{I}_{+}^{(k)}(T_i)) \\
            \mathcal{I}_{+}^{(k)}((\rho e)_i) &= (\rho e)(\mathcal{I}_{+}^{(k)}(\rho_i), \mathcal{I}_{+}^{(k)}(T_i))
          \end{aligned}

   The remainder of the hydrodynamics algorithm then proceeds unchanged.

-  2: on entering the Riemann solver, we recompute the
   thermodynamics on the interfaces to ensure that they are all
   consistent. This is done by taking the interface values of
   :math:`\rho`, :math:`e`, :math:`X_k`, and computing the corresponding pressure, :math:`p`
   from this.

-  3: This does the characteristic tracing using the
   :math:`(\tau, u, T)` eigensystem. Note: this is not widely
   implemented—see the for an
   implementation.

Resets
------

Density Resets
~~~~~~~~~~~~~~

Need to document density_reset_method

Energy
~~~~~~

Need to document allow_negative_energy and allow_small_energy

.. _app:hydro:flux_limiting:

Flux Limiting
~~~~~~~~~~~~~

Multi-dimensional hydrodynamic simulations often have numerical
artifacts that result from the sharp density gradients. A somewhat
common issue, especially at low resolution, is negative densities that
occur as a result of a hydro update. Castro contains a prescription
for dealing with negative densities, that resets the negative density
to be similar to nearby zones. Various choices exist for how to do
this, such as resetting it to the original zone density before the
update or resetting it to some linear combination of the density of
nearby zones. The reset is problematic because the strategy is not
unique and no choice is clearly better than the rest in all
cases. Additionally, it is not specified at all how to reset momenta
in such a case. Consequently, we desired to improve the situation by
limiting fluxes such that negative densities could not occur, so that
such a reset would in practice always be avoided. Our solution
implements the positivity-preserving method of :raw-latex:`\cite{hu:2013}`. This
behavior is controlled by
.

A hydrodynamical update to a zone can be broken down into an update
over every face of the zone where a flux crosses the face over the
timestep. The central insight of the positivity-preserving method is
that if the update over every face is positivity-preserving, then the
total update must be positivity-preserving as well. To guarantee
positivity preservation at the zone edge :math:`{\rm i}+1/2`, the flux
:math:`\mathbf{F}^{n+1/2}_{{\rm i}+1/2}` at that face is modified to become:

.. math:: \mathbf{F}^{n+1/2}_{{\rm i}+1/2} \rightarrow \theta_{{\rm i}+1/2} \mathbf{F}^{n+1/2}_{{\rm i}+1/2} + (1 - \theta_{{\rm i}+1/2}) \mathbf{F}^{LF}_{{\rm i}+1/2}, \label{eq:limited_flux}

where :math:`0 \leq \theta_{{\rm i}+1/2} \leq 1` is a scalar, and :math:`\mathbf{F}^{LF}_{{\rm i}+1/2}` is the Lax-Friedrichs flux,

.. math:: \mathbf{F}^{LF}_{{\rm i}+1/2} = \frac{1}{2}\left[\mathbf{F}^{n}_{{\rm i}} + \mathbf{F}^{n}_{{\rm i}+1} + \text{CFL}\frac{\Delta x}{\Delta t} \frac{1}{\alpha}\left(\mathbf{U}^{n}_{{\rm i}} - \mathbf{U}^{n}_{{\rm i}+1}\right)\right],

where :math:`0 < \text{CFL} < 1` is the CFL safety factor (the method is
guaranteed to preserve positivity as long as :math:`\text{CFL} < 1/2`), and
:math:`\alpha` is a scalar that ensures multi-dimensional correctness
(:math:`\alpha = 1` in 1D, :math:`1/2` in 2D, :math:`1/3` in 3D). :math:`\mathbf{F}_{{\rm
    i}}` is the flux of material evaluated at the zone center :math:`{\rm
  i}` using the cell-centered quantities :math:`\mathbf{U}`. The scalar
:math:`\theta_{{\rm i}+1/2}` is chosen at every interface by calculating the
update that would be obtained from , setting
the density component equal to a value just larger than the density floor,
, and solving
for the value of :math:`\theta` at the interface that makes the equality
hold. In regions where the density is not at risk of going negative,
:math:`\theta \approx 1` and the original hydrodynamic update is recovered.
Further discussion, including a proof of the method, a description of
multi-dimensional effects, and test verification problems, can be
found in :raw-latex:`\cite{hu:2013}`.

Gravity
=======

[ch:gravity]

.. _introduction-3:

Introduction
------------

Integration Strategy
~~~~~~~~~~~~~~~~~~~~

Castro uses subcycling to integrate levels at different timesteps.
The gravity algorithm needs to respect this. Self-gravity is computed
via multigrid. At coarse-fine interfaces, the stencil used in the
Laplacian understands the coarse-fine interface and is different than
the stencil used in the interior.

There are two types of
solves that we discuss with AMR:

-  *composite solve* : This is a multilevel solve, starting at
   a coarse level (usually level 0) and solving for the potential on
   all levels up to the finest level.

-  *level solve* : This solves for the potential only on
   a particular level. Finer levels are ignored. At coarse-fine
   interfaces, the data from the coarse levels acts as Dirichlet
   boundary conditions for the current-level-solve.

The overall integration strategy is unchanged from the discussion in
:raw-latex:`\cite{castro_I}`. Briefly:

-  At the beginning of a simulation, we do a multilevel composite
   solve (if ).

   We also do a multilevel composite solve after each regrid.

-  The old-time gravity on the coarse level is defined based on
   this composite solve, but we also do a level solve on the coarse
   level, and use it to determine the difference between the composite
   solve and the level solve, and store that in a MultiFab.

-  After the hydro advance on the coarse level, we do another level
   solve, and use the (level solve :math:`-` compositive solve) as a lagged
   predictor of how much we need to add back to that level solve to get
   an effective new-time value for phi on the coarse level, and that’s
   what defines the phi used for the new-time gravity

-  Then we do the fine grid timestep(s), each using the same
   strategy

-  At an AMR synchronization step across levels (see Section `2 <#sec:amr_synchronization>`__
   for a description of when these synchronizations occur), if we’re choosing
   to synchronize the gravitational field across levels ()
   we then do a solve starting from
   the coarse grid that adjusts for the mismatch between the fine-grid
   phi and the coarse-grid phi, as well as the mismatch between the
   fine-grid density fluxes and the coarse-grid density fluxes, and add
   the resulting sync solve phi to both the coarse and the fine level

   Thus, to within the gravity error tolerance, you get the same final
   result as if you had done a full composite solve at the end of the
   timestep (assuming ).

If you do , then you never do a full
multilevel solve, and the gravity on any level is defined only by the
solve on that level. The only time this would be appropriate is if
the fine level(s) cover essentially all of the mass on the grid for
all time.

Controls
~~~~~~~~

Castro can incorporate gravity as a constant, monopole approximation,
or a full Poisson solve. To enable gravity in the code, set:

::

    USE_GRAV = TRUE

in the GNUmakefile, and then turn it on in the inputs file
via castro.do_grav = 1. If you want to incorporate a point mass
(through castro.point_mass), you must have

::

    USE_POINTMASS = TRUE

in the GNUmakefile.

There are currently four options for how gravity is calculated,
controlled by setting . The options are
ConstantGrav, PoissonGrav, Monopole Grav or
PrescribedGrav. Again, these are only relevant if USE_GRAV =
TRUE in the GNUmakefile and castro.do_grav = 1 in the
inputs file. If both of these are set then the user is required
to specify the gravity type in the inputs file or the program will
abort.

Some additional notes:

-  For the full Poisson solver
   (), the behavior
   of the full Poisson solve / multigrid solver is controlled by
   and .

-  For isolated boundary conditions, and when
   , the parameters
   gravity.max_multipole_order and
   control the accuracy of
   the Dirichlet boundary conditions. These are described in
   Section `2.3.2 <#sec-poisson-3d-bcs>`__.

-  For MonopoleGrav, in 1D we must have coord_sys = 2, and in
   2D we must have coord_sys = 1.

The following parameters apply to gravity
solves:

-  : how should we calculate gravity?
   Can be ConstantGrav, PoissonGrav, MonopoleGrav, or
   PrescribedGrav

-  : if gravity.gravity_type =
   ConstantGrav, set the value of constant gravity (default: 0.0)

-  : gravity.gravity_type =
   PoissonGrav, do we perform the “sync solve"? (0 or 1; default: 0)

-  : if gravity.gravity_type
   = PoissonGrav, whether to perform a composite solve (0 or 1;
   default: 0)

-  : maximum level to solve
   for :math:`\phi` and :math:`\mathbf{g}`; above this level, interpolate from
   below (default: :math:`{\tt MAX\_LEV} - 1`)

-  : if gravity.gravity_type =
   PoissonGrav, this is the absolute tolerance for the Poisson
   solve. You can specify a single value for this tolerane (or do
   nothing, and get a reasonable default value), and then the absolute
   tolerance used by the multigrid solve is :math:`\text{abs\_tol} \times
     4\pi G\, \rho_{\text{max}}` where :math:`\rho_{\text{max}}` is the maximum
   value of the density on the domain. On fine levels, this absolute
   tolerance is multiplied by :math:`\text{ref\_ratio}^2` to account for the
   change in the absolute scale of the Laplacian operator. You can
   also specify an array of values for , one for each
   possible level in the simulation, and then the scaling by
   :math:`\text{ref\_ratio}^2` is not applied.

-  : if gravity.gravity_type
   = PoissonGrav, this is the relative tolerance for the Poisson
   solve. By default it is zero. You can specify a single value for
   this tolerance and it will apply on every level, or you can specify
   an array of values for , one for each possible level
   in the simulation. This replaces the old parameter
   .

-  : if
   gravity.gravity_type = PoissonGrav, this is the max :math:`\ell` value
   to use for multipole BCs (must be :math:`\geq 0`; default: 0)

-  : if
   gravity.gravity_type = PoissonGrav, evaluate BCs using exact sum
   (0 or 1; default: 0)

-  : ratio of dr for monopole gravity
   binning to grid resolution

The follow parameters affect the coupling of hydro and gravity:

-  : turn on/off gravity

-  : do we recompute the center
   used for the multipole gravity solver each step?

-  : point mass at the center of the star
   (must be :math:`\geq 0`; default: 0.0)

Note that in the following, MAX_LEV is a hard-coded parameter
in Source/Gravity.cpp which is currently set to 15. It
determines how many levels can be tracked by the Gravity object.

Types of Approximations
-----------------------

ConstantGrav
~~~~~~~~~~~~

Gravity can be defined as constant in direction and magnitude,
defined in the inputs file by

for example, to set the gravity to have magnitude :math:`9.8` in the
negative :math:`y`-direction if in 2D, negative :math:`z`-direction if in 3-D.
The actual setting is done in Gravity.cpp as:

::

     grav.setVal(const_grav, BL_SPACEDIM-1, 1, ng);

Note that at present we do not fill the gravitational potential :math:`\phi` in
this mode; it will be set to zero.

Note: ConstantGrav can only be used along a Cartesian direction
(vertical for 2D axisymmetric).

.. _sec-monopole-grav:

MonopoleGrav
~~~~~~~~~~~~

MonopoleGrav integrates the mass distribution on the grid
in spherical shells, defining an enclosed mass and uses this
to compute the gravitational potential and acceleration in a
spherically-symmetric fashion.

-  In 1D spherical coordinates we compute

   .. math:: g(r) = -\frac{G M_{\rm enclosed}}{ r^2}

   where :math:`M_{\rm enclosed}` is calculated from the density at the time
   of the call.

   For levels above the coarsest level we define the extent of that
   level’s radial arrays as ranging from the center of the star (:math:`r=0`)
   to the cell at that level farthest away from the origin. If there
   are gaps between fine grids in that range then we interpolate the
   density from a coarser level in order to construct a continuous
   density profile. We note that the location of values in the density
   profile and in the gravitational field exactly match the location of
   data at that level so there is no need to interpolate between points
   when mapping the 1D radial profile of :math:`g` back onto the original
   grid.

-  In 2D or 3D we compute a 1D radial average of density and use
   this to compute gravity as a one-dimensional integral, then
   interpolate the gravity vector back onto the Cartesian grid
   cells. At the coarsest level we define the extent of the 1D arrays
   as ranging from the center of the star to the farthest possible
   point in the grid (plus a few extra cells so that we can fill ghost
   cell values of gravity). At finer levels we first define a single
   box that contains all boxes on that fine level, then we interpolate
   density from coarser levels as needed to fill the value of density
   at every fine cell in that box. The extent of the radial array is
   from the center of the star to the *nearest* cell on one of the
   faces of the single box. This ensures that all cells at that
   maximum radius of the array are contained in this box.

   We then average the density onto a 1D radial array. We note that
   there is a mapping from the Cartesian cells to the radial array and
   back; unlike the 1D case this requires interpolation. We use quadratic
   interpolation with limiting so that the interpolation does not create
   new maxima or minima.

   The default resolution of the radial arrays at a level is the grid
   cell spacing at that level, i.e., :math:`\Delta r = \Delta x`. O For
   increased accuracy, one can define as a number
   greater than :math:`1` (:math:`2` or :math:`4` are recommended) and the spacing of the
   radial array will then satisfy :math:`\Delta x / \Delta r =` drdxfac.
   Individual Cartesian grid cells are subdivided by drdxfac in
   each coordinate direction for the purposing of averaging the density,
   and the integration that creates :math:`g` is done at the finer resolution
   of the new :math:`\Delta r`.

   Note that the center of the star is defined in the subroutine PROBINIT
   and the radius is computed as the distance from that center.

   .. raw:: latex

      \marginpar{\vskip-\baselineskip\raggedright\tiny\sffamily
      \hrule\smallskip{\color{red}there is an additional correction at the corners in{\tt
          make\_radial\_grav} that accounts for the volume in a shell that
        is not part of the grid}\par\smallskip\hrule}

 What about the potential in this case? when does
make_radial_phi come into play?

PoissonGrav
~~~~~~~~~~~

The most general case is a self-induced gravitational field,

.. math:: \mathbf{g}(\mathbf{x},t) = \nabla \phi

where :math:`\phi` is defined by solving

.. math:: \mathbf{\Delta} \phi = 4 \pi G \rho .\label{eq:Self Gravity}

We only allow PoissonGrav in 2D or 3D because in 1D, computing
the monopole approximation in spherical coordinates is faster and more
accurate than solving the Poisson equation.

Poisson Boundary Conditions: 2D
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In 2D, if boundary conditions are not periodic in both directions, we
use a monopole approximation at the coarsest level. This involves
computing an effective 1D radial density profile (on level =
0 only), integrating it outwards from the center to get the
gravitational acceleration :math:`\mathbf{g}`, and then integrating :math:`g`
outwards from the center to get :math:`\phi` (using :math:`\phi(0) = 0` as a
boundary condition, since no mass is enclosed at :math:`r = 0`). For more
details, see Section `2.2 <#sec-monopole-grav>`__.

.. _sec-poisson-3d-bcs:

Poisson Boundary Conditions: 3D
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following describes methods for doing isolated boundary
conditions. The best reference for Castro’s implementation of this
is :raw-latex:`\cite{katz:2016}`.

-  **Multipole Expansion**

   In 3D, by default, we use a multipole expansion to estimate the value
   of the boundary conditions. According to, for example, Jackson’s
   *Classical Electrodynamics* (with the corresponding change to
   Poisson’s equation for electric charges and gravitational
   ”charges”), an expansion in spherical harmonics for :math:`\phi` is

   .. math:: \phi(\mathbf{x}) = -G\sum_{l=0}^{\infty}\sum_{m=-l}^{l} \frac{4\pi}{2l + 1} q_{lm} \frac{Y_{lm}(\theta,\phi)}{r^{l+1}}, \label{spherical_harmonic_expansion}

   The origin of the coordinate system is taken to be the ``center``
   variable, that must be declared and stored in the ``probdata``
   module in your project directory. The validity of the expansion used
   here is based on the assumption that a sphere centered on
   ``center``, of radius approximately equal to the size of half the
   domain, would enclose all of the mass. Furthermore, the lowest order
   terms in the expansion capture further and further departures from
   spherical symmetry. Therefore, it is crucial that ``center`` be
   near the center of mass of the system, for this approach to achieve
   good results.

   The multipole moments :math:`q_{lm}` can be calculated by expanding the
   Green’s function for the Poisson equation as a series of spherical
   harmonics, which yields

   .. math:: q_{lm} = \int Y^*_{lm}(\theta^\prime, \phi^\prime)\, {r^\prime}^l \rho(\mathbf{x}^\prime)\, d^3x^\prime. \label{multipole_moments_original}

   Some simplification of Equation `[spherical_harmonic_expansion] <#spherical_harmonic_expansion>`__ can
   be achieved by using the addition theorem for spherical harmonics:

   .. math::

      \begin{aligned}
        &\frac{4\pi}{2l+1} \sum_{m=-l}^{l} Y^*_{lm}(\theta^\prime,\phi^\prime)\, Y_{lm}(\theta, \phi) = P_l(\text{cos}\, \theta) P_l(\text{cos}\, \theta^\prime) \notag \\
        &\ \ + 2 \sum_{m=1}^{l} \frac{(l-m)!}{(l+m)!} P_{l}^{m}(\text{cos}\, \theta)\, P_{l}^{m}(\text{cos}\, \theta^\prime)\, \left[\text{cos}(m\phi)\, \text{cos}(m\phi^\prime) + \text{sin}(m\phi)\, \text{sin}(m\phi^\prime)\right].\end{aligned}

   Here the :math:`P_{l}^{m}` are the associated Legendre polynomials and the
   :math:`P_l` are the Legendre polynomials. After some algebraic
   simplification, the potential outside of the mass distribution can be
   written in the following way:

   .. math:: \phi(\mathbf{x}) \approx -G\sum_{l=0}^{l_{\text{max}}} \left[Q_l^{(0)} \frac{P_l(\text{cos}\, \theta)}{r^{l+1}} + \sum_{m = 1}^{l}\left[ Q_{lm}^{(C)}\, \text{cos}(m\phi) + Q_{lm}^{(S)}\, \text{sin}(m\phi)\right] \frac{P_{l}^{m}(\text{cos}\, \theta)}{r^{l+1}} \right].

   The modified multipole moments are:

   .. math::

      \begin{aligned}
        Q_l^{(0)}   &= \int P_l(\text{cos}\, \theta^\prime)\, {r^{\prime}}^l \rho(\mathbf{x}^\prime)\, d^3 x^\prime \\
        Q_{lm}^{(C)} &= 2\frac{(l-m)!}{(l+m)!} \int P_{l}^{m}(\text{cos}\, \theta^\prime)\, \text{cos}(m\phi^\prime)\, {r^\prime}^l \rho(\mathbf{x}^\prime)\, d^3 x^\prime \\
        Q_{lm}^{(S)} &= 2\frac{(l-m)!}{(l+m)!} \int P_{l}^{m}(\text{cos}\, \theta^\prime)\, \text{sin}(m\phi^\prime)\, {r^\prime}^l \rho(\mathbf{x}^\prime)\, d^3 x^\prime.\end{aligned}

   Our strategy for the multipole boundary conditions, then, is to pick
   some value :math:`l_{\text{max}}` that is of sufficiently high order to
   capture the distribution of mass on the grid, evaluate the discretized
   analog of the modified multipole moments for :math:`0 \leq l \leq
   l_{\text{max}}` and :math:`1 \leq m \leq l`, and then directly compute the
   value of the potential on all of the boundary zones. This is
   ultimately an :math:`\mathcal{O}(N^3)` operation, the same order as the
   monopole approximation, and the wall time required to calculate the
   boundary conditions will depend on the chosen value of
   :math:`l_{\text{max}}`.

   The number of :math:`l` values calculated is controlled by
   in your inputs file. By
   default, it is set to ``0``, which means that a monopole
   approximation is used. There is currently a hard-coded limit of
   :math:`l_{\text{max}} = 50`. This is because the method used to generate the
   Legendre polynomials is not numerically stable for arbitrary :math:`l`
   (because the polynomials get very large, for large enough :math:`l`).

-  **Direct Sum**

   Up to truncation error caused by the discretization itself, the
   boundary values for the potential can be computed exactly by a direct
   sum over all cells in the grid. Suppose I consider some ghost cell
   outside of the grid, at location :math:`\mathbf{r}^\prime \equiv (x^\prime,
   y^\prime, z^\prime)`. By the principle of linear superposition as
   applied to the gravitational potential,

   .. math:: \phi(\mathbf{r}^\prime) = \sum_{\text{ijk}} \frac{-G \rho_{\text{ijk}}\, \Delta V_{\text{ijk}}}{\left[(x - x^\prime)^2 + (y - y^\prime)^2 + (z - z^\prime)^2\right]^{1/2}},

   where :math:`x = x(i)`, :math:`y = y(j)` and :math:`z = z(k)` are constructed in the
   usual sense from the zone indices. The sum here runs over every cell
   in the physical domain (that is, the calculation is :math:`\mathcal{O}(N^3)`
   for each boundary cell). There are :math:`6N^2` ghost cells needed for the
   Poisson solve (since there are six physical faces of the domain), so
   the total cost of this operation is :math:`\mathcal{O}(N^5)` (this only
   operates on the coarse grid, at present). In practice, we use the
   domain decomposition inherent in the code to implement this solve: for
   the grids living on any MPI task, we create six :math:`N^2` arrays
   representing each of those faces, and then iterate over every cell on
   each of those grids, and compute their respective contributions to all
   of the faces. Then, we do a global reduce to add up the contributions
   from all cells together. Finally, we place the boundary condition
   terms appropriate for each grid onto its respective cells.

   This is quite expensive even for reasonable sized domains, so this
   option is recommended only for analysis purposes, to check if the
   other methods are producing accurate results. It can be enabled by
   setting in your inputs file.

PrescribedGrav
~~~~~~~~~~~~~~

With PrescribedGrav [11]_, gravity can be defined as a function that
is specified by the user. The option is allowed in 2D and 3D. To
define the gravity vector, copy prescribe_grav_nd.f90 from
Source/gravity/ to your run directory. The makefile system will always
choose this local copy of the file over the one in another directory.
Then define the components of gravity inside a loop over the grid
inside the file. If your problem uses a radial gravity in the form
:math:`g(r)`, you can simply adapt
ca_prescribe_grav_gravityprofile, otherwise you will have to
adapt **ca_prescribe_grav**, both are located in
prescribed_grav_nd.90.

Point Mass
~~~~~~~~~~

Pointmass gravity works with all other forms of gravity, it is not a
separate option. Since the Poisson equation is linear in potential
(and its derivative, the acceleration, is also linear), the point mass
option works by adding the gravitational acceleration of the point
mass onto the acceleration from whatever other gravity type is under
in the simulation.

Note that point mass can be :math:`< 0`.

A useful option is point_mass_fix_solution. If set to
1, then it takes all zones that are adjacent to the location of the
center variable and keeps their density constant. Any changes
in density that occur after a hydro update in those zones are reset,
and the mass deleted is added to the pointmass. (If there is
expansion, and the density lowers, then the point mass is reduced and
the mass is added back to the grid). This calculation is done in
pm_compute_delta_mass() in
Source/gravity/pointmass_nd.f90.

GR correction
-------------

In the cases of compact objects or very massive stars, the general
relativity (GR) effect starts to play a role [12]_. First, we consider the hydrostatic equilibrium due to
effects of GR then derive GR-correction term for Newtonian gravity.
The correction term is applied to the monopole approximation only when
USE_GR = TRUE is set in the GNUmakefile.

The formulae of GR-correction here are based on :raw-latex:`\cite{grbk1}`. For
detailed physics, please refer to :raw-latex:`\cite{grbk2}`. For describing very
strong gravitational field, we need to use Einstein field equations

.. math::

   \label{field}
   R_{ik}-\frac{1}{2}g_{ik}R=\frac{\kappa}{c^{2}}T_{ik} \quad , \quad
   \kappa=\frac{8\pi G}{c^{2}}\quad ,

where :math:`R_{ik}` is the Ricci tensor, :math:`g_{ik}` is the metric tensor, :math:`R`
is the Riemann curvature, :math:`c` is the speed of light and :math:`G` is
gravitational constant. :math:`T_{ik}` is the energy momentum tensor, which
for ideal gas has only the non-vanishing components :math:`T_{00}` =
:math:`\varrho c^2` , :math:`T_{11}` = :math:`T_{22}` = :math:`T_{33}` = :math:`P` ( contains rest
mass and energy density, :math:`P` is pressure). We are interested in
spherically symmetric mass distribution. Then the line element :math:`ds`
for given spherical coordinate :math:`(r, \vartheta, \varphi)` has the
general form

.. math::

   \label{metric}
     ds^{2} = e^{\nu}c^{2}dt^{2}-e^{\lambda}dr^{2}-r^{2}(d\vartheta^{2}+\sin^{2}
     \vartheta d\varphi) \quad ,

with :math:`\nu = \nu(r)`, :math:`\lambda = \lambda(r)`. Now we can put the
expression of :math:`T_{ik}` and :math:`ds` into (`[field] <#field>`__), then field
equations can be reduced to 3 ordinary differential equations:

.. math::

   \label{diff1}
      \frac{\kappa P}{c^{2}} =
      e^{-\lambda}\left (\frac{\nu^{\prime}}{r}+\frac{1}{r^{2}} \right )-\frac{1}{r^{2}}
      \quad ,

.. math::

   \label{diff2}
     \frac{\kappa P}{c^{2}} =
     \frac{1}{2}e^{-\lambda}\left (\nu^{\prime\prime}+\frac{1}{2}{\nu^{\prime}}^{2}+\frac{\nu^
       {\prime}-\lambda^{\prime}}{r}
      -\frac{\nu^{\prime}\lambda^{\prime}}{2} \right ) \quad ,

.. math::

   \label{diff3}
     \kappa \varrho =
     e^{-\lambda}\left (\frac{\lambda^{\prime}}{r}-\frac{1}{r^{2}}\right )+\frac{1}{r^{2}} \quad ,

where primes means the derivatives with respect to :math:`r`. After
multiplying with :math:`4\pi r^2`, (`[diff3] <#diff3>`__) can be integrated and
yields

.. math::

   \label{gmass1}
     \kappa m = 4\pi r (1-e^{-\lambda}) \quad ,

the :math:`m` is called “gravitational mass” inside r defined as

.. math::

   \label{gmass2}
     m = \int_{0}^{r}4\pi r^{2}  \varrho dr\quad .

For the :math:`r = R`, :math:`m` becomes the mass :math:`M` of the star. :math:`M` contains
not only the rest mass but the whole energy (divided by :math:`c^2`), that
includes the internal and gravitational energy. So the :math:`\varrho =
\varrho_0 +U/c^2` contains the whole energy density :math:`U` and rest-mass
density :math:`\varrho_0`. Differentiation of (`[diff1] <#diff1>`__) with respect to
:math:`r` gives :math:`P = P^{\prime}(\lambda,\lambda^{\prime},
\nu,\nu^{\prime},r)`, where
:math:`\lambda,\lambda^{\prime},\nu,\nu^{\prime}` can be eliminated by
(`[diff1] <#diff1>`__), (`[diff2] <#diff2>`__), (`[diff3] <#diff3>`__). Finally we reach
*Tolman-Oppenheinmer-Volkoff(TOV)* equation for hydrostatic
equilibrium in general relativity:

.. math::

   \label{tov}
     \frac{dP}{dr} = -\frac{Gm}{r^{2}}\varrho \left (1+\frac{P}{\varrho
       c^{2}}\right )\left (1+\frac{4\pi r^3 P}{m c^{2}}\right ) \left (1-\frac{2Gm}{r c^{2}} \right)^{-1} \quad .

For Newtonian case :math:`c^2 \rightarrow  \infty`, it reverts to usual form

.. math::

   \label{newton}
     \frac{dP}{dr} = -\frac{Gm}{r^{2}}\varrho \quad .

Now we take effective monopole gravity as

.. math::

   \label{tov2}
   \tilde{g} = -\frac{Gm}{r^{2}} (1+\frac{P}{\varrho
     c^{2}})(1+\frac{4\pi r^3 P}{m c^{2}}) (1-\frac{2Gm}{r c^{2}})^{-1}  \quad .

For general situations, we neglect the :math:`U/c^2` and potential energy in
m because they are usually much smaller than :math:`\varrho_0`. Only when
:math:`T` reaches :math:`10^{13} K` (:math:`KT \approx m_{p} c^2`, :math:`m_p` is proton mass)
before it really makes a difference. So (`[tov2] <#tov2>`__) can be expressed
as

.. math::

   \label{tov3}
     \tilde{g} = -\frac{GM_{\rm enclosed}}{r^{2}} \left (1+\frac{P}{\varrho
       c^{2}} \right )\left (1+\frac{4\pi r^3 P}{M_{\rm enclosed} c^{2}} \right ) \left (1-\frac{2GM_{\rm enclosed}}{r c^{2}} \right )^{-1} \quad ,

where :math:`M_{enclosed}` has the same meaning as with the
MonopoleGrav approximation.

Hydrodynamics Source Terms
--------------------------

There are several options to incorporate the effects of gravity into
the hydrodynamics system. The main parameter here is
castro.grav_source_type.

-  castro.grav_source_type = 1 : we use a
   standard predictor-corrector formalism for updating the momentum and
   energy. Specifically, our first update is equal to :math:`\Delta t \times
     \mathbf{S}^n` , where :math:`\mathbf{S}^n` is the value of the source
   terms at the old-time (which is usually called time-level :math:`n`). At
   the end of the timestep, we do a corrector step where we subtract
   off :math:`\Delta t / 2 \times \mathbf{S}^n` and add on :math:`\Delta t / 2
     \times \mathbf{S}^{n+1}`, so that at the end of the timestep the
   source term is properly time centered.

-  castro.grav_source_type = 2 : we do something very
   similar to 1. The major difference is that when evaluating the
   energy source term at the new time (which is equal to :math:`\mathbf{u}
     \cdot \mathbf{S}^{n+1}_{\rho \mathbf{u}}`, where the latter is the
   momentum source term evaluated at the new time), we first update the
   momentum, rather than using the value of :math:`\mathbf{u}` before
   entering the gravity source terms. This permits a tighter coupling
   between the momentum and energy update and we have seen that it
   usually results in a more accurate evolution.

-  castro.grav_source_type = 3 : we do the same momentum
   update as the previous two, but for the energy update, we put all of
   the work into updating the kinetic energy alone. In particular, we
   explicitly ensure that :math:`(rho e)` maintains the same, and update
   :math:`(rho K)` with the work due to gravity, adding the new kinetic
   energy to the old internal energy to determine the final total gas
   energy. The physical motivation is that work should be done on the
   velocity, and should not directly update the temperature—only
   indirectly through things like shocks.

-  castro.grav_source_type = 4 : the energy update is done
   in a “conservative” fashion. The previous methods all evaluate
   the value of the source term at the cell center, but this method
   evaluates the change in energy at cell edges, using the
   hydrodynamical mass fluxes, permitting total energy to be conserved
   (excluding possible losses at open domain boundaries). See
   :raw-latex:`\cite{katzthesis}` for some more details.

Diffusion
=========

[ch:diffusion]

Thermal Diffusion
-----------------

Castro incorporates explicit thermal diffusion into the energy equation.
In terms of the specific internal energy, :math:`e`, this appears as:

.. math:: \rho \frac{De}{Dt} + p \nabla \cdot {\bf u}= \nabla \cdot {k_\mathrm{th}}\nabla T

where :math:`{k_\mathrm{th}}` is the thermal conductivity, with units
:math:`\mathrm{erg~cm^{-1}~s^{-1}~K^{-1}}`.

To see the similarity to the thermal diffusion equation, consider the special
case of constant conductivity, :math:`{k_\mathrm{th}}`, and density, and assume an
ideal gas, so :math:`e = c_v T`, where :math:`c_v` is the specific heat at constant volume.
Finally, ignore hydrodynamics, so :math:`{\bf u}= 0`. This gives:

.. math:: \frac{\partial T}{\partial t} = D \nabla^2 T

where :math:`D \equiv {k_\mathrm{th}}/(\rho c_v)`. Solving this equation
explicitly requires a timestep limiter of

.. math:: \Delta t_\mathrm{diff} \le \frac{1}{2} \frac{\Delta x^2}{D}

(this is implemented in in
Castro/Source/driver/timestep.F90).

Support for diffusion must be compiled into the code by setting
USE_DIFFUSION = TRUE in your GNUmakefile. It is treated
explicitly, by constructing the contribution to the evolution as a
source term. This is time-centered to achieve second-order accuracy
in time.

The following parameter affects diffusion:

-  : enable thermal diffusion (0 or 1; default 0)

A pure diffusion problem (with no hydrodynamics) can be run by setting

::

    castro.diffuse_temp = 1
    castro.do_hydro = 0

To complete the setup, a thermal conductivity must be specified. The
interface for the conductivity is:

::

      subroutine thermal_conductivity(eos_state, therm_cond)
        
        use extern_probin_module, only: const_conductivity

        type (eos_t), intent(in) :: eos_state
        real (kind=dp_t), intent(inout) :: therm_cond

The density, temperature, and mass fractions come in through the
eos_state type. An EOS call is done in Castro just before the
call to , so you can assume that the entire
state is consistent.

There are two conductivity routines provided with Castro by default:

-  constant : A simple constant thermal conductivity. This can be
   selected by setting

   ::

       Conductivity_dir := constant

   in your GNUmakefile. To set the value of the conductivity (e.g., to
   :math:`100`), you add to your probin file’s &extern namelist:

   ::

       const_conductivity = 100.0

-  constant_opacity : A simple constant opacity. This is
   converted to an opacity as:

   .. math:: {k_\mathrm{th}}= \frac{16 \sigma_B T^3}{3 \kappa_\mathrm{const} \rho}

   where :math:`\kappa_\mathrm{const}` is the opacity, with units :math:`\mathrm{cm^2~g^{-1}}`.
   This is selected by setting

   ::

       Conductivity_dir := constant_opacity

   in your GNUmakefile. To set the value of the opacity, e.g., to
   0.2 (e.g., for electron scattering), set:

   ::

       const_opacity = 0.2

   in the &extern namelist of your probin.

The diffusion approximation breaks down at the surface of stars,
where the density rapidly drops and the mean free path becomes
large. In those instances, you should use the flux limited diffusion
module in Castro to evolve a radiation field. However, if your
interest is only on the diffusion in the interior, you can use
the parameter to specify a density,
below which, diffusion is not modeled. This is implemented in the
code by zeroing out the conductivity and skipping the estimation
of the timestep limit in these zones.

A simple test problem that sets up a Gaussian temperature profile
and does pure diffusion is provided as diffusion_test.

Enthalpy Diffusion
------------------

Castro can also diffuse enthalpy

Note this uses the same interface for the transport coefficients as
thermal diffusion, so the two cannot be used at the same time.

Species Diffusion
-----------------

Castro can also diffuse species.

Note this uses the same interface for the transport coefficients as
thermal diffusion, so the two cannot be used at the same time.

Viscosity
---------

Rotation
========

[ch:rotation]

.. _introduction-4:

Introduction
------------

Currently, Castro supports contant, solid-body rotation about a fixed
(in space and time) axis in 2D and 3D by transforming the evolution
equations to the rotating frame of reference.

To include rotation you must set

::

    USE_ROTATION = TRUE

in the GNUMakefile. Rotation can then be enabled via

::

    castro.do_rotation = 1

in the inputs file. The rotational period must then be set via
. The rotational period is internally
converted to an angular frequency for use in the source term
equations.

The axis of rotation currently depends on the dimensionality of the
problem and the value of coord_sys; in all cases, however, the
default axis of rotation points from center, which is typically
defined in a Prob_$(DIM)d.f90 routine, to the typical “vertical
direction.” The vertical direction is defined as follows:

-  2D

   -  coord_sys = 0, (x,y): out of the (x,y)-plane along the “z”-axis

   -  coord_sys = 1, (r,z): along the z-axis

-  3D

   -  coord_sys = 0, (x,y,z): along the z-axis

To change these defaults, modify the omega vector in the
ca_rotate routine found in the Rotate_$(DIM)d.f90 file.

The main parameters that affect rotation are:

-  : include rotation as a forcing
   term (0 or 1; default: 0)

-  : period (s) of rotation
   (default: 0.0)

-  : d(period) / dt for rotation
   (default: 0.0)

-  : whether to
   include the centrifugal forcing (default: 1)

-  : whether to
   include the Coriolis forcing (default: 1)

-  : whether to
   include the forcing from the time derivative of the rotation
   frequency (default: 1)

-  : whether state
   variables are measured in the rotating frame (default: 1)

-  : method of updating the
   energy during a rotation update (default: 4)

-  : for the Coriolis
   term, which mixes momenta in the source term, whether we should
   solve for the update implicitly (default: 1)

-  : rotation axis (default: 3
   (Cartesian); 2 (cylindrical))

For completeness, we show below a derivation of the source terms that
appear in the momentum and total energy evolution equations upon
switching to a rotating reference frame.

Coordinate transformation to rotating frame
-------------------------------------------

.. raw:: latex

   \centering

.. figure:: tframes
   :alt: [fig:sec:rot:frames] Inertial frame :math:`C` and
   non-inertial frame :math:`\tilde{C}`. We consider a fluid element
   :math:`P`, whose distance in the two frames is related by
   :math:`{\bf r} = \tilde{\bf{r}} + {\bf l}`

   [fig:sec:rot:frames] Inertial frame :math:`C` and
   non-inertial frame :math:`\tilde{C}`. We consider a fluid element
   :math:`P`, whose distance in the two frames is related by
   :math:`{\bf r} = \tilde{\bf{r}} + {\bf l}`

Consider an inertial reference frame :math:`C` and a non-inertial
reference frame :math:`\widetilde{C}` whose origins are separated by the
vector :math:`\boldsymbol{l}` (see Figure \ `[fig:sec:rot:frames] <#fig:sec:rot:frames>`__). The non-inertial frame is rotating about the axis
:math:`\boldsymbol{\omega}` with a *constant* angular velocity :math:`\omega`;
furthermore, we assume the *direction* of the rotational axis is
fixed. Consider a fluid element at the point :math:`P` whose location is
given by :math:`{\bf r}` in :math:`C` and by :math:`\widetilde{{\bf r}}` in
:math:`\widetilde{C}`:

.. math:: {\bf r}= \widetilde{{\bf r}}+ \boldsymbol{l},

or in component notation

.. math::

   \label{eq:r}
       r_i\boldsymbol{e_i} = \widetilde{r_i}\widetilde{\boldsymbol{e_i}} + l_i\boldsymbol{e_i},

where :math:`\boldsymbol{e_i}` and :math:`\widetilde{\boldsymbol{e_i}}` are the :math:`i`\ th unit
vectors in the :math:`C` and :math:`\widetilde{C}` coordinate systems,
respectively. The total time rate of change of `[eq:r] <#eq:r>`__ is given by

.. math::

   \label{eq:vcomp}
       \frac{Dr_i}{Dt}\boldsymbol{e_i} = \frac{D\widetilde{r_i}}{Dt}\widetilde{\boldsymbol{e_i}} + \widetilde{r_i}\frac{D\widetilde{\boldsymbol{e_i}}}{Dt} + \frac{Dl_i}{Dt}\boldsymbol{e_i},

where we have used the fact that the unit vectors of the inertial
frame :math:`C` are not moving (or at least can be considered stationary,
and the change in :math:`\boldsymbol{l}` gives the relative motion of the two
coordinate systems). By definition, a unit vector can not change its
length, and therefore the only change of :math:`\widetilde{\boldsymbol{e_i}}` with
time can come from changing direction. This change is carried out by
a rotation about the :math:`\boldsymbol{\omega}` axis, and the tip of the unit
vector moves circumferentially, that is

.. math::

   \label{eq:etilde-rot}
       \frac{D\widetilde{\boldsymbol{e_i}}}{Dt} = \boldsymbol{\omega}\times\widetilde{\boldsymbol{e_i}}.

Plugging `[eq:etilde-rot] <#eq:etilde-rot>`__ into `[eq:vcomp] <#eq:vcomp>`__ and switching back to
vector notation, we have

.. math::

   \label{eq:r-dot}
       \frac{D{\bf r}}{Dt} = \frac{D\widetilde{{\bf r}}}{Dt} + \boldsymbol{\omega}\times\widetilde{{\bf r}}+ \frac{D\boldsymbol{l}}{Dt}.

The left hand side of `[eq:r-dot] <#eq:r-dot>`__ is interpretted as the velocity
of the fluid element as seen in the inertial frame; the first term on the
right hand side is the velocity of the fluid element as seen by a
stationary observer in the rotating frame :math:`\widetilde{C}`. The second
and third terms on the right hand side of `[eq:r-dot] <#eq:r-dot>`__ describe the
additional velocity due to rotation and translation of the frame
:math:`\widetilde{C}` as seen in :math:`C`. In other words,

.. math::

   \label{eq:v}
       \boldsymbol{v}= \widetilde{\boldsymbol{v}}+ \boldsymbol{\omega}\times\widetilde{{\bf r}}+ \boldsymbol{v_l},

where we use :math:`\boldsymbol{v_l}` to represent the translational velocity.

Similarly, by taking a second time derivative of `[eq:v] <#eq:v>`__ we have

.. math::

   \label{eq:a}
       \frac{D\boldsymbol{v}}{Dt} = \frac{D\widetilde{\boldsymbol{v}}}{Dt} + 2\boldsymbol{\omega}\times\widetilde{\boldsymbol{v}}+ \boldsymbol{\omega}\times\left[\boldsymbol{\omega}\times\widetilde{{\bf r}}\right] + \frac{D\boldsymbol{v_l}}{Dt}.

Henceforth we will assume the two coordinate systems are not
translating relative to one another, :math:`\boldsymbol{v_l} = 0`. It is
also worth mentioning that derivatives with respect to spatial
coordinates do not involve additional terms due to rotation,
i.e. :math:`\boldsymbol{\nabla}\cdot\boldsymbol{v}= \boldsymbol{\nabla}\cdot\widetilde{\boldsymbol{v}}`.
Because of this, the continuity equation remains unchanged in the
rotating frame:

.. math::

   \label{eq:cont-rot}
       \frac{\partial \rho}{\partial t} = -\boldsymbol{\nabla}\cdot\left(\rho\widetilde{\boldsymbol{v}}\right),

or

.. math::

   \label{eq:cont-rot-total}
       \frac{D\rho}{Dt} = -\rho\boldsymbol{\nabla}\cdot\widetilde{\boldsymbol{v}}.

Momentum equation in rotating frame
-----------------------------------

The usual momentum equation applies in an inertial frame:

.. math::

   \label{eq:mom1}
       \frac{D\left(\rho\boldsymbol{v}\right)}{Dt} = -\rho\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}- \boldsymbol{\nabla}p + \rho{\bf g}.

Using the continuity equation, `[eq:cont-rot-total] <#eq:cont-rot-total>`__, and substituting for
the terms in the rotating frame from `[eq:a] <#eq:a>`__, we have from `[eq:mom1] <#eq:mom1>`__:

.. math::

   \begin{aligned}
       \rho\left(\frac{D\widetilde{\boldsymbol{v}}}{Dt} + 2\boldsymbol{\omega}\times\widetilde{\boldsymbol{v}}+ \boldsymbol{\omega}\times\left[\boldsymbol{\omega}\times\widetilde{{\bf r}}\right]\right) - \rho\boldsymbol{v}\boldsymbol{\nabla}\cdot\boldsymbol{v}&=& -\rho\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}- \boldsymbol{\nabla}p + \rho{\bf g}\nonumber \\
       \rho\left(\frac{\partial\widetilde{\boldsymbol{v}}}{\partial t} + \widetilde{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\widetilde{\boldsymbol{v}}\right) &=& -\boldsymbol{\nabla}p + \rho{\bf g}- 2\rho\boldsymbol{\omega}\times\widetilde{\boldsymbol{v}}- \rho\boldsymbol{\omega}\times\left[\boldsymbol{\omega}\times\widetilde{{\bf r}}\right] \nonumber \\
     \frac{\partial\left(\rho\widetilde{\boldsymbol{v}}\right)}{\partial t} &=& -\boldsymbol{\nabla}\cdot\left(\rho\widetilde{\boldsymbol{v}}\widetilde{\boldsymbol{v}}\right) - \boldsymbol{\nabla}p + \rho{\bf g}- 2\rho\boldsymbol{\omega}\times\widetilde{\boldsymbol{v}}\nonumber \\
     & & -\ \rho\boldsymbol{\omega}\times\left[\boldsymbol{\omega}\times\widetilde{{\bf r}}\right]\label{eq:mom-rot}
     \end{aligned}

or

.. math::

   \label{eq:mom-rot-tot}
       \frac{D\left(\rho\widetilde{\boldsymbol{v}}\right)}{Dt} = -\rho\widetilde{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\widetilde{\boldsymbol{v}}- \boldsymbol{\nabla}p + \rho{\bf g}- 2\rho\boldsymbol{\omega}\times\widetilde{\boldsymbol{v}}- \rho\boldsymbol{\omega}\times\left[\boldsymbol{\omega}\times\widetilde{{\bf r}}\right].

Energy equations in rotating frame
----------------------------------

From `[eq:mom-rot-tot] <#eq:mom-rot-tot>`__, we have the velocity evolution equation in
a rotating frame

.. math::

   \label{eq:v-rot}
       \frac{D\widetilde{\boldsymbol{v}}}{Dt} = -\frac{1}{\rho}\boldsymbol{\nabla}p + {\bf g}- 2\boldsymbol{\omega}\times\widetilde{\boldsymbol{v}}- \boldsymbol{\omega}\times\left[\boldsymbol{\omega}\times\widetilde{{\bf r}}\right].

The kinetic energy equation can be obtained from `[eq:v-rot] <#eq:v-rot>`__ by
mulitplying by :math:`\rho\widetilde{\boldsymbol{v}}`:

.. math::

   \begin{aligned}
       \rho\widetilde{\boldsymbol{v}}\cdot\frac{D\widetilde{\boldsymbol{v}}}{Dt} &=& -\widetilde{\boldsymbol{v}}\cdot\boldsymbol{\nabla}p + \rho\widetilde{\boldsymbol{v}}\cdot{\bf g}- 2\rho\widetilde{\boldsymbol{v}}\cdot\left[\boldsymbol{\omega}\times\widetilde{\boldsymbol{v}}\right] - \rho\widetilde{\boldsymbol{v}}\cdot\left\{\boldsymbol{\omega}\times\left[\boldsymbol{\omega}\times\widetilde{{\bf r}}\right]\right\} \nonumber \\
       \frac{1}{2}\frac{D\left(\rho\widetilde{\boldsymbol{v}}\cdot\widetilde{\boldsymbol{v}}\right)}{Dt} - \frac{1}{2}\widetilde{\boldsymbol{v}}\cdot\widetilde{\boldsymbol{v}}\frac{D\rho}{Dt} &=& -\widetilde{\boldsymbol{v}}\cdot\boldsymbol{\nabla}p + \rho\widetilde{\boldsymbol{v}}\cdot{\bf g}- \rho\widetilde{\boldsymbol{v}}\cdot\left[\left(\boldsymbol{\omega}\cdot\widetilde{{\bf r}}\right)\boldsymbol{\omega}- \rho\omega^2\widetilde{{\bf r}}\right] \nonumber \\
       \frac{1}{2}\frac{D\left(\rho\widetilde{\boldsymbol{v}}\cdot\widetilde{\boldsymbol{v}}\right)}{Dt} &=& -\frac{1}{2}\rho\widetilde{\boldsymbol{v}}\cdot\widetilde{\boldsymbol{v}}\boldsymbol{\nabla}\cdot\widetilde{\boldsymbol{v}}- \widetilde{\boldsymbol{v}}\cdot\boldsymbol{\nabla}p + \rho\widetilde{\boldsymbol{v}}\cdot{\bf g}- \rho\widetilde{\boldsymbol{v}}\cdot\left[\left(\boldsymbol{\omega}\cdot\widetilde{{\bf r}}\right)\boldsymbol{\omega}- \rho\omega^2\widetilde{{\bf r}}\right]. \label{eq:ekin-rot-total}
     \end{aligned}

The internal energy is simply advected, and, from the first law of
thermodynamics, can change due to :math:`pdV` work:

.. math::

   \label{eq:eint-rot-total}
       \frac{D\left(\rho e\right)}{Dt} = -\left(\rho e + p\right)\boldsymbol{\nabla}\cdot\widetilde{\boldsymbol{v}}.

Combining `[eq:ekin-rot-total] <#eq:ekin-rot-total>`__ and `[eq:eint-rot-total] <#eq:eint-rot-total>`__ we can
get the evolution of the total specific energy in the rotating frame,
:math:`\rho \widetilde{E} = \rho e + \frac{1}{2}\rho\widetilde{\boldsymbol{v}}\cdot\widetilde{\boldsymbol{v}}`:

.. math::

   \begin{aligned}
       \frac{D\left(\rho e\right)}{Dt} + \frac{1}{2}\frac{D\left(\rho\widetilde{\boldsymbol{v}}\cdot\widetilde{\boldsymbol{v}}\right)}{Dt} &=& -\left(\rho e + p + \frac{1}{2}\rho\widetilde{\boldsymbol{v}}\cdot\widetilde{\boldsymbol{v}}\right)\boldsymbol{\nabla}\cdot\widetilde{\boldsymbol{v}}- \widetilde{\boldsymbol{v}}\cdot\boldsymbol{\nabla}p + \rho\widetilde{\boldsymbol{v}}\cdot{\bf g}-\rho\widetilde{\boldsymbol{v}}\cdot\left[\left(\boldsymbol{\omega}\cdot\widetilde{{\bf r}}\right)\boldsymbol{\omega}- \rho\omega^2\widetilde{{\bf r}}\right]\nonumber \\
       \frac{D\left(\rho \widetilde{E}\right)}{Dt} &=& -\rho\widetilde{E}\boldsymbol{\nabla}\cdot\widetilde{\boldsymbol{v}}- \boldsymbol{\nabla}\cdot\left(p\widetilde{\boldsymbol{v}}\right) + \rho\widetilde{\boldsymbol{v}}\cdot{\bf g}- \rho\widetilde{\boldsymbol{v}}\cdot\left[\left(\boldsymbol{\omega}\cdot\widetilde{{\bf r}}\right)\boldsymbol{\omega}- \rho\omega^2\widetilde{{\bf r}}\right] \label{eq:etot-rot-total}
     \end{aligned}

or

.. math::

   \label{eq:etot-rot}
       \frac{\partial\left(\rho\widetilde{E}\right)}{\partial t} = -\boldsymbol{\nabla}\cdot\left(\rho\widetilde{E}\widetilde{\boldsymbol{v}}+ p\widetilde{\boldsymbol{v}}\right) + \rho\widetilde{\boldsymbol{v}}\cdot{\bf g}- \rho\widetilde{\boldsymbol{v}}\cdot\left[\left(\boldsymbol{\omega}\cdot\widetilde{{\bf r}}\right)\boldsymbol{\omega}- \rho\omega^2\widetilde{{\bf r}}\right].

Switching to the rotating reference frame
-----------------------------------------

If we choose to be a stationary observer in the rotating reference
frame, we can drop all of the tildes, which indicated terms in the
non-inertial frame :math:`\widetilde{C}`. Doing so, and making sure we
account for the offset, :math:`\boldsymbol{l}`, between the two coordinate systems, we obtain
the following equations for hydrodynamics in a rotating frame of
reference:

.. math::

   \begin{aligned}
       \frac{\partial\rho}{\partial t} &=& -\boldsymbol{\nabla}\cdot\left(\rho\boldsymbol{v}\right) \label{eq:cont-rot-switch} \\
       \frac{\partial \left(\rho\boldsymbol{v}\right)}{\partial t} &=& -\boldsymbol{\nabla}\cdot\left(\rho\boldsymbol{v}\boldsymbol{v}\right) - \boldsymbol{\nabla}p + \rho{\bf g}- 2\rho\boldsymbol{\omega}\times\boldsymbol{v}- \rho\left(\boldsymbol{\omega}\cdot{\bf r}\right)\boldsymbol{\omega}+ \rho\omega^2{\bf r}\label{eq:mom-rot-switch} \\
       \frac{\partial\left(\rho E\right)}{\partial t} &=& -\boldsymbol{\nabla}\cdot\left(\rho E\boldsymbol{v}+ p\boldsymbol{v}\right) + \rho\boldsymbol{v}\cdot{\bf g}- \rho\left(\boldsymbol{\omega}\cdot{\bf r}\right)\left(\boldsymbol{\omega}\cdot\boldsymbol{v}\right) + \rho\omega^2\left(\boldsymbol{v}\cdot{\bf r}\right). \label{eq:etot-rot-switch}
     \end{aligned}

Adding the forcing to the hydrodynamics
---------------------------------------

There are several ways to incorporate the effect of the rotation
forcing on the hydrodynamical evolution. We control this through the
use of the runtime parameter . This
is an integer with values currently ranging from 1 through 4, and
these values are all analogous to the way that gravity is used to
update the momentum and energy. For the most part, the differences are
in how the energy update is done:

-  castro.rot_source_type = 1 : we use a
   standard predictor-corrector formalism for updating the momentum and
   energy. Specifically, our first update is equal to :math:`\Delta t \times
     \mathbf{S}^n` , where :math:`\mathbf{S}^n` is the value of the source
   terms at the old-time (which is usually called time-level :math:`n`). At
   the end of the timestep, we do a corrector step where we subtract
   off :math:`\Delta t / 2 \times \mathbf{S}^n` and add on :math:`\Delta t / 2
     \times \mathbf{S}^{n+1}`, so that at the end of the timestep the
   source term is properly time centered.

-  castro.rot_source_type = 2 : we do something very
   similar to 1. The major difference is that when evaluating the
   energy source term at the new time (which is equal to :math:`\mathbf{u}
     \cdot \mathbf{S}^{n+1}_{\rho \mathbf{u}}`, where the latter is the
   momentum source term evaluated at the new time), we first update the
   momentum, rather than using the value of :math:`\mathbf{u}` before
   entering the rotation source terms. This permits a tighter coupling
   between the momentum and energy update and we have seen that it
   usually results in a more accurate evolution.

-  castro.rot_source_type = 3 : we do the same momentum
   update as the previous two, but for the energy update, we put all of
   the work into updating the kinetic energy alone. In particular, we
   explicitly ensure that :math:`(rho e)` maintains the same, and update
   :math:`(rho K)` with the work due to rotation, adding the new kinetic
   energy to the old internal energy to determine the final total gas
   energy. The physical motivation is that work should be done on the
   velocity, and should not directly update the temperature – only
   indirectly through things like shocks.

-  castro.rot_source_type = 4 : the energy update is done
   in a “conservative” fashion. The previous methods all evaluate
   the value of the source term at the cell center, but this method
   evaluates the change in energy at cell edges, using the
   hydrodynamical mass fluxes, permitting total energy to be conserved
   (excluding possible losses at open domain boundaries). Additionally,
   the velocity update is slightly different—for the corrector step,
   we note that there is an implicit coupling between the velocity
   components, and we directly solve this coupled equation, which
   results in a slightly better coupling and a more accurate evolution.

The other major option is castro.implicit_rotation_update.
This does the update of the Coriolis term in the momentum equation
implicitly (e.g., the velocity in the Coriolis force for the zone
depends on the updated momentum). The energy update is unchanged.

A detailed discussion of these options and some verification
tests is presented in :raw-latex:`\cite{katz:2016}`.

Radiation
=========

.. _introduction-5:

Introduction
------------

Castro has three radiation solvers:

-  SingleGroupSolver: this solver does not have radiation
   pressre. It is pure hydro plus radiation diffusion. This is only
   applicable when the medium is optically thick and the pressure is small.

-  SGFLDSolver: this is the gray flux-limited diffusion
   radiation hydrodymamics solver. Here the radiation pressure is
   separate from the gas pressure, and both pressures participate in
   the Riemann solver.

-  MGFLDSolver: this is the multigroup flux-limited diffusion
   radiation hydrodynamics solver. As with the gray solver, radiation
   pressure contributes to the pressure in the Riemann solver. Here a
   number of energy groups are used to represent the radiation field,
   and the opacity can be frequency-dependent.

The gray solver has a comoving frame mode and a mixed frame mode,
whereas the MG solver uses the comoving frame approach. More details
about the formulation and algorithm can be found in the series of
Castro papers.

.. _getting-started-1:

Getting Started
---------------

Getting the Code
~~~~~~~~~~~~~~~~

The Castro radiation solver is part of the main Castro git repo,
so you already have all the Castro code and problem setups
to exercise radiation. The only other requirement is a copy
of the Hypre library. Hypre provides the algebraic multigrid
solvers used by the implicit radiation update. You can get
a copy at https://computation.llnl.gov/casc/linear_solvers/sls_hypre.html. You will need to follow their installation instructions.

In addition to the environment variables you set for the main
Castro hydrodynamics problems, you also need to tell the code
where to find Hypre. This is done via one of two variables:

-  the environment variable HYPRE_DIR should
   point to the location of your Hypre installation
   (e.g., ) or
   can be set directly in
   .
   This applies is you build with USE_OMP=FALSE

-  the variable HYPRE_OMP_DIR should be set (either as an
   environment variable or in
   ) to the directory
   for openmp enabled Hypre (e.g.,
   ) if you build with
   USE_OMP=TRUE

Now go to a “run” directory, say
,
edit the file GNUmakefile, and set

-  COMP = your favorite compiler suite (e.g., gnu, pgi, intel)

-  DIM = 1 or 2 or 3

-  —this is important. This tells the build system to
   compile in, and link the the radiation code.

Then type make to generate an executable file.

Microphysics: EOS, Network, and Opacity
---------------------------------------

EOS
~~~

Castro provides several types of equation of state (EOS), including
gamma-law and Helmholtz. To use the gamma-law EOS, set

::

    EOS_DIR := gamma_law_general

in the GNUmakefile.

The original Helmholtz EOS for stellar interiors includes a radiation
contribution. However, for radiation hydrodynamics calculations, the
radition contribution should be taken out of EOS because radiation has
been treated in other places. To use Helmholtz EOS, we will use the
version in Microphysics, as with the pure hydrodynamics code, but
this will interpret the preprocessor variable and
disable the radiation portion of the EOS [13]_ If you have your own EOS, you
can put it in Microphysics.

EOS Parameters
^^^^^^^^^^^^^^

The following parameters affect how the radiation solver used the EOS:

-  radiation.do_real_eos = 1

   Usually you do not want to change this from the default. Setting
   this to 0 is only for contrived tests that assume the
   specific heat is in the form of a power-law,

   .. math:: c_v = \mathrm{const}\ \rho^m T^{-n}

Network
~~~~~~~

The radiation solver uses the same networks as we saw for pure hydro,
so nothing needs to change here. Again, if you are not modeling
reactions, then the general_null network can be used to define
the appropriate composition for your problem.

Opacity
~~~~~~~

By default, we assume that

.. math:: \kappa = \mathrm{const}\ \rho^{m} T^{-n} \nu^{p} , \label{eq:kappa}

where :math:`\kappa` is either Planck or Rosseland mean absorption
coefficients, :math:`\rho` is density, :math:`T` is temperature, :math:`\nu` is
frequency, and :math:`m`, :math:`n` and :math:`p` are constants. For the gray solver,
:math:`p = 0`. If Equation (\ `[eq:kappa] <#eq:kappa>`__) is sufficient, set

::

    Opacity_dir := null

in GNUmakefile. Otherwise, put your own opacity in
and set
the input parameter, radiation.use_opacity_table_module = 1 (see
§ \ `3.3.1 <#sec:opacpars>`__).

Some notes:

-  Here, :math:`\kappa` has units of :math:`\mathrm{cm}^{-1}`. Some papers or
   texts may instead have an implicit density factor in :math:`\kappa`,
   yielding units :math:`\mathrm{cm}^2~\mathrm{g}^{-1}`.

-  Castro allows for two temperatures (different radiation and gas
   temperature, so :math:`E_\mathrm{r} \ne a T_\mathrm{gas}^4`).
   Correspondingly, Castro cares about both the Planck mean,
   :math:`\kappa_P`, and Rosseland mean, :math:`\kappa_R`, opacities—these have
   different weightings.

   If we set :math:`\kappa_P \Delta x \gg 1` (:math:`\kappa_P` is really large),
   then the two temperatures become the same.

   If we set :math:`\kappa_P = \kappa_R`, then we can see how different the
   two temperature are.

   In an optically thick medium, we would not expect the two temperatures
   to be very different.

.. _sec:opacpars:

Opacity Parameters
^^^^^^^^^^^^^^^^^^

The parameters describing the opacity include:

-  radiation.use_opacity_table_module = 0

   For neutrino problems, this parameter is not ignored. For photon
   problems, this determines whether the opacity module at
   Opacity_dir (which is set in GNUmakefile) will be used to
   compute opacities. If this is set to 1, the following parameters
   for opacities will be ignored.

-  For the Planck mean opacity of the form in Eq. (\ `[eq:kappa] <#eq:kappa>`__),
   the following parameters set the coefficient and exponents:

   -  radiation.const_kappa_p = -1.0

   -  radiation.kappa_p_exp_m = 0.0

   -  radiation.kappa_p_exp_n = 0.0

   -  radiation.kappa_p_exp_p = 0.0

-  For the Rosseland mean opacity of the form in Eq. (\ `[eq:kappa] <#eq:kappa>`__),
   the following parameters set the coefficient and exponents:

   -  radiation.const_kappa_r = -1.0

   -  radiation.kappa_r_exp_m = 0.0

   -  radiation.kappa_r_exp_n = 0.0

   -  radiation.kappa_r_exp_p = 0.0

-  For the scattering coefficient of the form in Eq. (\ `[eq:kappa] <#eq:kappa>`__),
   the following parameters set the coefficient and exponents:

   -  radiation.const_scattering = 0.0

   -  radiation.scattering_exp_m = 0.0

   -  radiation.scattering_exp_n = 0.0

   -  radiation.scattering_exp_p = 0.0

-  radiation.kappa_r_floor = 0.0

   Floor for Rosseland mean.

-  radiation.do_kappa_stm_emission = 0

   If it is 1, correction for stimulated emission is applied to Planck mean as
   follows

   .. math::

      \kappa = \mathrm{const}\ \rho^{m} T^{-n} \nu^{p}
          \left [1-\exp{\left (-\frac{h\nu}{k T} \right )} \right ].

-  radiation.surface_average = 2

   How the averaging of opacity is done from faces to center for
   the radiation solver. 0 is arithmetic averaging, 1
   is harmonic averaging, and 2 is a combination of the two.
   This is implemented in RAD_?D.F in kavg.

Note that the unit for opacities is :math:`\mathrm{cm}^{-1}`. For
the gray solver, the total opacity in the diffusion coefficient is the sum
of kappa_r and scattering, whereas for the MG solver,
there are two possibilities. If const_kappa_r is greater than
0, then the total opacity is set by kappa_r alone, otherwise
the total opacity is the sum of kappa_p and scattering.

Radiation Solver Physics
------------------------

In this section, we list some radiation related parameters that you
can set in an inputs file. Here are some important parameters:

-  radiation.SolverType:

   Set it to 5 for the gray solver, and 6 for the MG solver.

-  castro.do_hydro

   Usually you want to set it to 1. If it is set to 0,
   hydro will be turned off, and the calculation will only solve
   radiation diffusion equation.

-  castro.do_radiation

   If it is 0, the calculation will be pure hydro.

Below are more parameters. For each parameter, the default value is
on the right-hand side of the equal sign.

.. _sec:bothpar:

Verbosity and I/O
~~~~~~~~~~~~~~~~~

-  radiation.v = 0

   Verbosity

-  radiation.verbose = 0

   Verbosity

-  radiation.plot_lambda = 0

   If 1, save flux limiter in plotfiles.

-  radiation.plot_kappa_p = 0

   If 1, save Planck mean opacity in plotfiles.

-  radiation.plot_kappa_r = 0

   If 1, save Rosseland mean opacity in plotfiles.

-  radiation.plot_lab_Er = 0

   If 1, save lab frame radiation energy density in plotfiles.
   This flag is ignored when the mixed-frame gray solver is used.

-  radiation.plot_com_flux = 0

   If 1, save comoving frame radiation flux in plotfiles.

-  radiation.plot_lab_flux = 0

   If 1, save lab frame radiation flux in plotfiles.

.. _sec:fluxlimiter:

Flux Limiter and Closure
~~~~~~~~~~~~~~~~~~~~~~~~

-  radiation.limiter = 2

   Possible values are:

   -   0: No flux limiter

   -   2: Approximate limiter of Levermore & Pomraning

   -  12: Bruenn’s limiter

   -  22: Larsen’s square root limiter

   -  32: Minerbo’s limiter

-  radiation.closure = 3

   Possible values are:

   -  0: :math:`f = \lambda`, where :math:`f` is the scalar Eddington factor
      and :math:`\lambda` is the flux limiter.

   -  1: :math:`f = \frac{1}{3}`

   -  2: :math:`f = 1 - 2 \lambda`

   -  3: :math:`f = \lambda + (\lambda R)^2`, where :math:`R` is the radiation
      Knudsen number.

   -  4: :math:`f = \frac{1}{3} + \frac{2}{3} (\frac{F}{cE})^2`, where
      :math:`F` is the radiation flux, :math:`E` is the radiation energy density,
      and :math:`c` is the speed of light.

Note the behavior of the radiative flux in the optically thin and
optically thick limits. The flux limiter, :math:`\lambda = \lambda(R)`,
where

.. math:: R = \frac{|\nabla E_r^{(0)}|}{\chi_R E_r^{(0)}}

Regardless of the limiter chosen, when we are optically thick,
:math:`\chi_R \rightarrow \infty`, :math:`R \rightarrow 0`, and :math:`\lambda \rightarrow 1/3`.
The radiative flux then becomes

.. math::

   F_r^{(0)} = -\frac{c\lambda}{\chi_R} \nabla E_r^{(0)} \rightarrow
     \frac{1}{3} \frac{c}{\chi_R} \nabla E_r^{(0)}

And when we are optically thin, :math:`\chi_R \rightarrow 0`, :math:`R \rightarrow \infty`,
and :math:`\lambda \rightarrow 1/R = \chi_R E_r^{(0)}/{|\nabla E_r^{0}|}`, and
the radiative flux then becomes

.. math::

   F_r^{(0)} = -\frac{c\lambda}{\chi_R} \nabla E_r^{(0)} \rightarrow
     -\frac{c}{\chi_R}\frac{\chi_R E_r^{(0)}}{|\nabla E_r^{0}|}
       \nabla E_r^{(0)} = -c E_r^{0}

See Krumholz et al. 2007 for some discussion on this.

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

Castro needs to know about the boundary conditions for both
the hydrodynamics and radiation portions of the evolution.

Hydrodynamics Evolution
^^^^^^^^^^^^^^^^^^^^^^^

For the hydrodynamics portion of the solve, the boundary conditions
for the normal hydrodynamic state values will be set by the problem’s
hypfill routine (which typically just calls filcc to handle
the usual hydrodynamics boundary types: outflow, symmetry, etc.).

A corresponding radfill routine needs to be written to fill the
ghost cells for the radiation energy density during the hydrodynamics
evolution. Again, this usually will just default to calling
filcc.

Note: if any of the hydrodynamic boundary conditions types are set
to Inflow, then you will need to ensure that the radfill
routine explicitly handles the boundary condition implementation
for the radiation energy density in that case—the filcc
routine will not do a hydrodynamic Inflow boundary.

Radiation Evolution
^^^^^^^^^^^^^^^^^^^

The following parameters are for radiation boundary in the diffusion
equation. They do not affect hydrodynamic boundaries.

-  radiation.lo_bc

   This sets the action to take at the lower edge of the domain in
   each coordinate direction. Possible values are:

   -  101 *Dirchlet*:

      Specify the radiation energy density on the boundary.
      For gray radiation, this could be :math:`E_r = a T^4`.

      For multigroup radiation, Castro stores the energy density as
      :math:`\mathrm{erg}~\mathrm{cm}^{-3}`, so the total radiation energy
      can be found by simply summing over the groups. So if you want
      to set the radiation BCs using the Planck function, you simply
      multiply by the group width—see Exec/radiation_tests/RadSphere/Tools/radbc.f90
      for an example.

   -  102 *Neumann*:

      Here, you specify the radiation flux on the boundary. For gray
      radiation, this is the expression given in the gray Castro paper
      (Eq. 7, 8),

      .. math:: F_r = - \frac{c\lambda}{\kappa_R} \nabla E_r

      where :math:`\lambda` is the flux limiter.

      Note that if your boundary represents an incoming flux through
      a vacuum (like stellar irradiation), then :math:`\kappa \rightarrow 0`, leaving

      .. math:: F_r = -c E_r

      (see § \ `4.2 <#sec:fluxlimiter>`__) in that case.

   -  104 *Marshak* (vacuum):

      Here, you specify the incident flux and the outside is a vacuum.
      This differs from the Neumann condition because there is also a
      flux coming from inside, for the net flux across the boundary is
      different than the incident flux.

   -  105 *Sanchez-Pomraning*:

      This is a modified form of the Marshak boundary condition that works with FLD.
      This is like the Marshak condition, but :math:`\lambda = 1/3` is not assumed inside
      the boundary (optical thickness).

-  radiation.hi_bc

   See radiation.lo_bc.

-  radiation.lo_bcflag = 0 0 0

   If it is 0, bcval is used for that dimension, otherwise
   subroutine rbndry in RadBndry_1d.f90 is called to set
   boundary conditions.

-  radiation.hi_bcflag = 0 0 0

   See radiation.lo_bcflag

-  radiation.lo_bcval = 0.0 0.0 0.0

   The actual value to impose for the boundary condition type set by
   radiation.lo_bc. This parameter is interpreted differently
   depending on the boundary condition:

   -  Dirchlet: Dirichlet value of rad energy density

   -  Neumann: inward flux of rad energy

   -  Marshak: incident flux

   -  Sanchez-Pomraning: incident flux

-  radiation.hi_bcval = 0.0 0.0 0.0

   See radiation.lo_bcval

Convergence
~~~~~~~~~~~

For the gray solver, there is only one iteration in the scheme,
whereas for the MG solver, there are two iterations with an inner
iteration embedded inside an outer iteration. In the following, the
iteration in the gray solver will also be referred as the outer
iteration for convenience. The parameters for the inner iteration are
irrelevant to the gray solver.

radiation.maxiter = 50
    | 
    | Maximal number of outer iteration steps.

radiation.miniter = 1
    | 
    | Minimal number of outer iteration steps.

radiation.reltol = 1.e-6
    | 
    | Relative tolerance for the outer iteration.

radiation.abstol = 0.0
    | 
    | Absolute tolerance for the outer iteration.

radiation.maxInIter = 30
    | 
    | Maximal number of inner iteration steps.

radiation.minInIter = 1
    | 
    | Minimal number of inner iteration steps.

radiation.relInTol = 1.e-4
    | 
    | Relative tolerance for the inner iteration.

radiation.absInTol = 0.0
    | 
    | Absolute tolerance for the inner iteration.

radiation.convergence_check_type = 0
    | 
    | For the MG solver only. This specifiy the way of checking the
      convergence of an outer iteration. Possible values are

    -  0: Check :math:`T`, :math:`Y_e`, and the residues of the equations for
       :math:`\rho e` and :math:`\rho Y_e`

    -  1: Check :math:`\rho e`

    -  2: Check the residues of the equations for :math:`\rho e` and :math:`\rho Y_e`

    -  3: Check :math:`T` and :math:`Y_e`

.. _sec:graypar:

Parameters for Gray Solver
~~~~~~~~~~~~~~~~~~~~~~~~~~

radiation.comoving = 1
    | 
    | Do we use the comoving frame approach?

radiation.Er_Lorentz_term = 1
    | 
    | If the mixed-frame approach is taken, this parameter decides whether
      Lorentz transformation terms are retained.

radiation.delta_temp = 1.0
    | 
    | This is used in computing numerical derivativas with respect to :math:`T`.
      So it should be a small number compared with :math:`T`, but not too small.

radiation.update_limiter = 1000
    | 
    | Stop updating flux limiter after update_limiter iteration steps.

radiation.update_planck = 1000
    | 
    | Stop updating Planck mean opacity after update_planck iteration steps.

radiation.update_rosseland = 1000
    | 
    | Stop updating Rosseland mean opacity after update_rosseland iteration steps.

Grouping in the MG Solver
~~~~~~~~~~~~~~~~~~~~~~~~~

We provide two methods of setting up groups based upon logarithmic
spacing. In both methods, you must provide:

radiation.nGroups
    | 
    | Number of groups.

radiation.lowestGroupHz
    | 
    | Frequency of the lower bound for the first group.

In addition, if the parameter groupGrowFactor is provided, then
the first method will be used, otherwise the second method will be
used. In the first way, you must also provide firstGroupWidthHz
(the width of the first group). The width of other groups is set to
be groupGrowFactor times the width of its immediately preceding
group. In the second way, you must provide highestGroupHz as
the upper bound of the last group. It should be noted that
lowestGroupHz can be 0 in the first method, but not the second
method. However, when we compute the group-integrated Planck
function, the lower bound for the first group and the upper bound for
the last group are assumed to be 0 and :math:`\infty`, respectively.

.. _sec:mgpar:

Parameters for MG Solver
~~~~~~~~~~~~~~~~~~~~~~~~

radiation.delta_e_rat_dt_tol = 100.0
    | 
    | Maximally allowed relative change in :math:`e` during one time step.

radiation.delta_T_rat_dt_tol = 100.0
    | 
    | Maximally allowed relative change in :math:`T` during one time step.

radiation.delta_Ye_dt_tol = 100.0
    | 
    | Maximally allowed absolute change in :math:`Y_e` during one tim estep.

radiation.fspace_advection_type = 2
    | 
    | Possible value is 1 or 2. The latter is better.

radiation.integrate_Planck = 1
    | 
    | If 1, integrate Planck function for each group. For the first
      group, the lower bound in the integration is assumed to be 0 no
      matter what the grouping is. For the last group, the upper bound in
      the integration is assumed to be :math:`\infty`.

radiation.matter_update_type = 0
    | 
    | How to update matter. 0 is proabaly the best.

radiation.accelerate = 2
    | 
    | The inner iteration of the MG solver usually requires an
      acceleration scheme. Choices are

    -  0: No acceleration

    -  1: Local acceleration

    -  2: Gray acceleration

radiation.skipAccelAllowed = 0
    | 
    | If it is set to 1, skip acceleration if it does not help.

radiation.n_bisect = 1000
    | 
    | Do bisection for the outer iteration after n_bisec iteration steps.

radiation.use_dkdT = 1
    | 
    | If it is 1, :math:`\frac{\partial \kappa}{\partial T}` is retained in the
      Jacobi matrix for the outer (Newton) iteration.

radiation.update_opacity = 1000
    | 
    | Stop updating opacities after update_opacity outer iteration steps.

radiation.inner_update_limiter = 0
    | 
    | Stop updating flux limiter after inner_update_limiter inner
      iteration steps. If it is 0, the limiter is lagged by one outer
      iteration. If it is -1, the limiter is lagged by one time step. If
      the inner iteration has difficulty in converging, setting this
      parameter it to -1 can help. Since the flux limiter is only a
      kludge, it is justified to lag it.

.. _sec:hypre:

Linear System Solver
~~~~~~~~~~~~~~~~~~~~

There are a number of choices for the linear system solver. The
performance of the solvers usually depends on problems and the
computer. So it is worth trying a few solvers to find out which one
is best for your problem and computer.

: the linear solver
in Hypre to use. The available choices are:

-  0: SMG

-  1: PFMG (:math:`\ge` 2-d only)

-  100: AMG using ParCSR ObjectType

-  102: GMRES using ParCSR ObjectType

-  103: GMRES using SStruct ObjectType

-  104: GMRES using AMG as preconditioner

-  109: GMRES using Struct SMG/PFMG as preconditioner

-  150: AMG using ParCSR ObjectType

-  1002: PCG using ParCSR ObjectType

-  1003: PCG using SStruct ObjectType

As a general rule, the SMG is the most stable solver, but is usually
the slowest. The asymmetry in the linear system comes from the
adaptive mesh, so the PFMG should be your first choice. Note: in
you cannot use PFMG.

Setting this to 109 (GMRES using Struct SMG/PFMG as preconditioner)
should work reasonably well for most problems.

(default: 40):
Maximal number of iteration in Hypre.

(default: 1.e-10):
Relative tolerance in Hypre

(default: 0):
Absolute tolerance in Hypre

(default: 0):
Verbosity

(default: 0):
Verbosity

(default: 0):
Verbosity for level_solver_flag :math:`<` 100

(default: 0):
Verbosity for level_solver_flag :math:`>=` 100

Output
------

Gray Solver
~~~~~~~~~~~

For the gray radiation solver, the radiation energy density is stored in plotfiles
as rad. Note that this quantity has units of :math:`\mathrm{erg~cm^{-3}}`, which
is different that the specify internal energy of the gas :math:`\mathrm{erg~g^{-1}}`.

Particles
=========

Tracer particles
----------------

Tracer particles are to track the Lagrangian evolution of a model fluid using discrete particles. In hydrodynamical simulations based on an Eulerian grid (including CASTRO), thermodynamic variables at a given time are derived by solving the equations of motion of a fluid between cells. Therefore, in this scheme, the physical quantities that we can access to are not discretized quantities at any given position, but rather average values over each cell. However, employing discrete particles, passively advected with the fluid flow, allows us to obtain local instantaneous thermodynamic variables, such as the temperature and the density, at well-defined positions, independent of the spatial resolution, i.e., the spatial cell size. This means that we can follow the evolution of the fluid at any given position and time.

CASTRO provides a tracer particle scheme with useful options. In this scheme, particles are advanced using the midpoint method either with the cell-centered velocities or the face-centered velocities (Marker-And-Cell method) [14]_. The number and the initial positions of particles are flexibly determined according to the purpose of a given model.

Initializing the Particles
--------------------------

One must include the tracer particles in the GNUmakefile by setting

.. raw:: latex

   \vspace{0.1in}

.

And the particles can be initialized via

.. raw:: latex

   \vspace{0.1in}

= 1

in the file.

If one wants to investigate the evolution of fluid motions starting from specific positions (or a certain range of area or volume), one should manually specify the positions of particles by providing an input file containing the total number and the initial positions of the particles.
The input file should be in the same directory where your inputs file is located. The name of the input file is determined via :

.. raw:: latex

   \vspace{0.1in}

| Here *particle_file* is the user-specified name of the file. The first line in this file is
  assumed to contain the number of particles. Each line after that contains the positions in a coordinate system adopted for your model. For 3-D cartesian coordinates,
| :math:`x ~y ~z`
| For example, an input file for a model fluid with 6 particles in 2-D Cartesian coordinates may look like,

::

    6
    3.28125e+08 9.9198e+08 
    5.46875e+08 9.9198e+08 
    7.65625e+08 9.9198e+08 
    9.84375e+08 9.9198e+08 
    1.20312e+09 9.9198e+08 
    1.42188e+09 9.9198e+08 

According to this input file, the 6 particles will be positioned at the same height (same :math:`y` coordinate in the second column), equally spaced in :math:`x` direction (the first column except for the particle number on the first line) from :math:`3.28\times10^{8} {\rm ~cm}` to :math:`1.42\times 10^{9} {\rm ~cm}`.

.. _particles:output_file:

Output file
-----------

The output files are stored in a directory whose name is determined by a variable particles.timestamp_dir. For example, if the variable is set as follows,

.. raw:: latex

   \vspace{0.1in}
   {\tt  {\bf particles.timestamp\_dir=}} {\em particle\_dir}

,

A directory *particle_dir* is automatically made with the directories for the main CASTRO output file (plt***\*) once a simulation starts and the particle output files are stored inside that directory.

.. raw:: latex

   \vspace{0.05in}

The name of the output file consists of Timestamp\_ along with a number at the end. The number increases (typically from 00) as more processors are involved in following the trajectories of particles. In parallel computing, a computational domain is divided according to the number of processors requested. Then each processor only follows the particles on the domain assigned to that processor and records their positions and velocities at any given time in a different output file. Since it is possible for particles to move from one domain to another during the evolution, its history can be stored in different files. More output files (with larger numbers at the end of the file name) can be produced as more processors track the particles.

.. raw:: latex

   \vspace{0.05in}

By default, the output file contains the positions and velocities of all particles at a given time, meaning [:math:`3+ 2\times`\ dimensionality] columns. For example, for particles in a 3-D domain, the columns in the output file are,

.. raw:: latex

   \vspace{0.1in}

:math:`{\rm index1}~~{\rm index2}~~x~~ y~~ z~~ t~~ v_{\rm x} ~~v_{\rm y}~~ v_{\rm z}~~ [\rho ~~ T]`

The first two integers correspond to the particle index and the processor number.
One should use the two numbers in order to identify a particle and extract its history (i.e., the trajectory in Figure `[fig:particletrajectory] <#fig:particletrajectory>`__).

.. raw:: latex

   \centering

|A model atmosphere (*left* panel) and the trajectories of 500 particles (*right* panel) following the fluid motion on the atmosphere. The particles are initially positioned at five different heights, :math:`y=13000\mathrm{~km},~11000\mathrm{~km},~ 8000\mathrm{~km},~ 6000\mathrm{~km}, ~38000\mathrm{~km}` (100 particles at each height). In the *left* panel, the arrows roughly show the fluid motion. In the *right* panel, the solid lines represent the trajectories of the particles. |
|A model atmosphere (*left* panel) and the trajectories of 500 particles (*right* panel) following the fluid motion on the atmosphere. The particles are initially positioned at five different heights, :math:`y=13000\mathrm{~km},~11000\mathrm{~km},~ 8000\mathrm{~km},~ 6000\mathrm{~km}, ~38000\mathrm{~km}` (100 particles at each height). In the *left* panel, the arrows roughly show the fluid motion. In the *right* panel, the solid lines represent the trajectories of the particles. |

One can also add the last two columns :math:`[\rho ~~ T]`, i.e., the local density and local temperature of fluid at the position of each particle by setting the following,

.. raw:: latex

   \vspace{0.1in}

| = 1,
| = 1.

For example, let’s consider 10 particles on a domain. If 4 out 10 particles are initially on a processor and the rest are on another processor, this means two processors are tracking the particles and two output files are produced. In the output file written by the processor with 4 particles, one can find that four lines are stored at the same time and each line corresponds to each particle info. while in the other output file for the other 6 particles, 6 lines are stored at the same time.

.. raw:: latex

   \vspace{0.05in}

If , the particle data are stored in a binary file along with the main CASTRO output plotfile in directories plt*****/Tracer/.

Run-time Screen Output
~~~~~~~~~~~~~~~~~~~~~~

The verbosity written to the screen at run-time is constrolled by setting:

| = 0 or 1 (default: 0)

Equation of State and Burning Network
=====================================

Equation of State
-----------------

Standard Castro EOSes
~~~~~~~~~~~~~~~~~~~~~

Castro is written in a modular fashion so that the EOS and network
burning routines can be supplied by the user. However, for the
examples presented later we use several EOS and network routines
that come with the Microphysics distribution.

Castro relies on routines to calculate the equation of state (EOS)
of a fluid, as well as a species network to define the components of
the fluid. The network optionally has the ability to do nuclear burning,
but for this section its main purpose is in defining the species so that
the EOS can calculate fluid properties that depend on composition, such
as electron fraction.

By default, Castro comes with the gamma_law
EOS. This represents a gamma law gas, with equation of state:

.. math:: P = (\gamma - 1) \rho e.

The gas is currently assumed to be monatomic and ideal. (Only a
restricted set of thermodynamic variables are actually calculated,
the minimum necessary for the hydrodynamics. A fuller set of
thermodynamic variables, for example the entropy from the
Sackur-Tetrode equation, are calculated in the gamma_law_general
EOS inside the Microphysics repository.)

EOS Interfaces and Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each EOS should have two main routines by which it interfaces to the
rest of Castro. At the beginning of the simulation, eos_init
will perform any initialization steps and save EOS variables (mainly
``smallt``, the temperature floor, and ``smalld``, the
density floor). Then, whenever you want to call the EOS, use

.. math:: {\tt call\ eos (eos\_input, eos\_state)}.

The first argument specifies the inputs to the EOS. The options
that are currently available are stored in
EOS/eos_data.F90, and are always a combination of two
thermodynamic quantities. For example, eos_input_rt means
that we call the EOS with rho (density) and T (temperature)
and we expect the EOS to return the associated thermodynamic
quantities such as internal energy e and entropy s.

We note that for real (non-analytic) equations of state
in which :math:`\rho`, :math:`T` and species are the independent variables, such
as the Helmholtz EOS, eos_input_rt directly calls the EOS
and obtains the other thermodynamic variables. But for other inputs,
e.g. eos_input_re, a Newton-Raphson iteration is performed
to find the density or temperature that corresponds to the given
input.

The eos_state variable is a Fortran derived type (similar to
a C struct). It stores a complete set of thermodynamic
variables. When calling the EOS, you should first fill the variables
that are the inputs, for example with

::

      use eos_type_module
      ...
      type (eos_t) :: eos_state
      ...
      eos_state % rho = state(i,j,k,URHO)
      eos_state % T   = state(i,j,k,UTEMP)
      eos_state % e   = state(i,j,k,UEINT) / state(i,j,k,URHO)
      eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)
      eos_state % aux = state(i,j,k,UFX:UFX+naux-1) / state(i,j,k,URHO)

Whenever the ``eos_state`` type is initialized, the
thermodynamic state variables are filled with unphysical numbers. If
you do not input the correct arguments to match your input quantities,
the EOS will call an error. This means that it is good
practice to fill the quantities that will be iterated over with an
initial guess. Indeed, this initial guess is typically required for
equations of state that iterate over this variable, as the values
they are initialized with will likely not
converge. Usually a prior value of the temperature or density suffices
if it’s available, but if not then use ``small_temp`` or
``small_dens``.

If you are interested in using more realistic and sophisticated equations of
state, you should download the `Microphysics <https://github.com/starkiller-astro/Microphysics>`__
repository. This is a collection of microphysics routines that are compatible with the
BoxLib codes. We refer you to the documentation in that repository for how to set it up
and for information on the equations of state provided. That documentation
also goes into more detail about the details of the EOS code, in case you are interested in
how it works (and in case you want to develop your own EOS).

Nuclear Network
---------------

The nuclear network serves two purposes: it defines the fluid components used
in both the equation of state and the hydrodynamics, and it evolves those
components through a nuclear burning step. Castro comes with a general_null
network (which lives in the Networks/ directory). This is a bare interface for a
nuclear reaction network. No reactions are enabled, and no auxiliary variables
are accepted. It contains several sets of isotopes; for example,
Networks/general_null/triple_alpha_plus_o.net would describe the
isotopes needed to represent the triple-\ :math:`\alpha` reaction converting
helium into carbon, as well as oxygen and iron.

The main interface file, network.f90, is a wrapper function. The
actual network details are defined in actual_network.f90, a
file which is automatically generated in your work directory when you compile.
It supplies the number and names of species and auxiliary variables, as
well as other initializing data, such as their mass numbers, proton numbers,
and the binding energies.

The burning front-end interface, Networks/burner.f90, accepts a different
derived type called the burn_t type. Like the eos_t, it has entries
for the basic thermodynamic quantities:

::

      use burn_type_module
      ...
      type (burn_t) :: burn_state
      ...
      burn_state % rho = state(i,j,k,URHO)
      burn_state % T   = state(i,j,k,UTEMP)
      burn_state % e   = state(i,j,k,UEINT) / state(i,j,k,URHO)
      burn_state % xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)

It takes in an input burn_t and returns an output burn_t after
the burning has completed. The nuclear energy release can be computed by
taking the difference of burn_state_out % e and
burn_state_in % e. The species change can be computed analogously.
In normal operation in Castro  the integration occurs over a time interval
of :math:`\Delta t/2`, where :math:`\Delta t` is the hydrodynamics timestep.

If you are interested in using actual nuclear burning networks,
you should download the `Microphysics <https://github.com/starkiller-astro/Microphysics>`__
repository. This is a collection of microphysics routines that are compatible with the
BoxLib codes. We refer you to the documentation in that repository for how to set it up
and for information on the networks provided. That documentation
also goes into more detail about the details of the network code, in case you are interested in
how it works (and in case you want to develop your own network).

Required Thermodynamics Quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three input modes are required of any EOS:

-  eos_input_re: :math:`\rho`, :math:`e`, and :math:`X_k` are input

-  eos_input_rt: :math:`\rho`, :math:`T`, and :math:`X_k` are input

-  eos_input_rp: :math:`\rho`, :math:`P`, and :math:`X_k` are input

The eos_t derived type holds a large number of thermodynamics
quantities, but not all of these are needed for basic
Castro operation. The main quantities that any EOS in any mode needs to
supply, if they are not input, are:

-  eos_state % T: the temperature

-  eos_state % P: total pressure

-  eos_state % e: the specific energy

-  eos_state % gam1: the first adiabatic index,
   :math:`\Gamma_1 = d\log P / d\log \rho |_s`

Additionally the eos_input_re mode also needs to supply:

-  eos_state % cs: the adiabatic sound speed

-  eos_state % dpdr_e: the derivative, :math:`\partial
       p/\partial \rho |_e`—note that the specific internal energy, :math:`e`
   is held constant here.

-  eos_state % dpde: the derivative, :math:`\partial p /
       \partial e |_\rho`

For radiation hydro, the eos_input_rt model needs to supply:

-  eos_state % cv: the specific heat capacity.

Other quantities (e.g., entropy) might be needed for the derived
variables that are optional output into the plotfiles.

For burning, when the temperature equation is evolved, the EOS
needs to supply:

-  eos_state % dhdX(nspec): the derivative of the
   specific enthalpy with respect to mass fraction at constant
   :math:`T` and :math:`p`. This is commonly computed as:

   .. math:: \xi_k = e_{X_k} + \frac{1}{p_\rho} \left (\frac{p}{\rho^2} - e_\rho \right ) p_{X_k}\enskip .

   with

   .. math::

      \begin{aligned}
      p_{X_k} &=& \left .\frac{\partial p}{\partial \bar{A}} \right |_{\rho, T, \bar{Z}}
                \frac{\partial \bar{A}}{\partial X_k} +
                \left . \frac{\partial p}{\partial \bar{Z}} \right |_{\rho, T, \bar{A}}
                \frac{\partial \bar{Z}}{\partial X_k} \nonumber \\
              &=& -\frac{\bar{A}^2}{A_k}
                \left .\frac{\partial p}{\partial \bar{A}} \right |_{\rho, T, \bar{Z}} +
                \frac{\bar{A}}{A_k} \left (Z_k - \bar{Z} \right )
                \left . \frac{\partial p}{\partial \bar{Z}} \right |_{\rho, T, \bar{A}}\enskip,\end{aligned}

   .. math::

      \begin{aligned}
      e_{X_k} &=& \left . \frac{\partial e }{\partial \bar{A}} \right |_{\rho, T, \bar{Z}}
              \frac{\partial \bar{A}}{\partial X_k} +
              \left .\frac{\partial e}{\partial \bar{Z}} \right |_{\rho, T, \bar{A}}
              \frac{\partial \bar{Z}}{\partial X_k} \nonumber \\
              &=& -\frac{\bar{A}^2}{A_k}
              \left . \frac{\partial e }{\partial \bar{A}} \right |_{\rho, T, \bar{Z}} +
              \frac{\bar{A}}{A_k} \left (Z_k - \bar{Z}\right )
              \left .\frac{\partial e}{\partial \bar{Z}} \right |_{\rho, T, \bar{A}}\enskip.\end{aligned}

(see :raw-latex:`\cite{maestro:III}`, Appendix A).

.. _amr-1:

AMR
===

Our approach to adaptive refinement in Castro uses a nested hierarchy
of logically-rectangular grids with simultaneous refinement of the
grids in both space and time. The integration algorithm on the grid
hierarchy is a recursive procedure in which coarse grids are advanced
in time, fine grids are advanced multiple steps to reach the same time
as the coarse grids and the data at different levels are then
synchronized.

During the regridding step, increasingly finer grids
are recursively embedded in coarse grids until the solution is
sufficiently resolved. An error estimation procedure based on
user-specified criteria (described in Section `1 <#sec:tagging>`__)
evaluates where additional refinement is needed
and grid generation procedures dynamically create or
remove rectangular fine grid patches as resolution requirements change.

A good introduction to the style of AMR used here is in Lecture 1
of the Adaptive Mesh Refinement Short Course at
https://ccse.lbl.gov/people/jbb/shortcourse/lecture1.pdf.

.. _sec:tagging:

Tagging for Refinement
----------------------

Castro determines what zones should be tagged for refinement at the
next regridding step by using a set of built-in routines that test on
quantities such as the density and pressure and determining whether
the quantities themselves or their gradients pass a user-specified
threshold. This may then be extended if amr.n_error_buf :math:`> 0`
to a certain number of zones beyond these tagged zones. This section
describes the process by which zones are tagged, and describes how to
add customized tagging criteria.

The routines for tagging cells are located in the
file in the Source/driver/ directory. (These are
dimension-agnostic routines that loop over all three dimensional
indices even for 1D or 2D problems.) The main routines are
, , ,
, and . They refine based on
density, temperature, pressure, velocity, and radiation energy density
(if enabled), respectively. The same approach is used for all of
them. As an example, we consider the density tagging routine. There
are four parameters that control tagging. If the density in a zone is
greater than the user-specified parameter , then that
zone will be tagged for refinement, but only if the current AMR level
is less than the user-specified parameter .
Similarly, if the absolute density gradient between a zone and any
adjacent zone is greater than the user-specified parameter
, that zone will be tagged for refinement, but only
if we are currently on a level below
. Note that setting denerr alone
will not do anything; you’ll need to set max_dengrad_lev :math:`>=
1` for this to have any effect.

All four of these parameters are set in the namelist
in your probin file. If left unmodified, they
default to a value that means we will never tag. The complete set of
parameters that can be controlled this way is the following:

-  density:

   -  value: denerr, max_denerr_lev

   -  gradient: dengrad, max_dengrad_lev

-  temperature:

   -  value: temperr, max_temperr_lev

   -  gradient: tempgrad, max_tempgrad_lev

-  velocity (magnitude):

   -  value: velerr, max_velerr_lev

   -  gradient: velgrad, max_velgrad_lev

-  pressure:

   -  value: presserr, max_presserr_lev

   -  gradient: pressgrad, max_pressgrad_lev

-  radiation energy density:

   -  value: raderr, max_raderr_lev

   -  gradient: radgrad, max_radgrad_lev

Since there are multiple algorithms for determining
whether a zone is tagged or not, it is worthwhile to specify
in detail what is happening to a zone in the code during this step.
We show this in the following pseudocode section. A zone
is tagged if the variable itag = SET, and is not tagged
if itag = CLEAR (these are mapped to 1 and 0, respectively).

::

    itag = CLEAR

    for errfunc[k] from k = 1 ... N
        // Three possibilities for itag: SET or CLEAR or remaining unchanged
        call errfunc[k](itag)  
    end for

In particular, notice that there is an order dependence of this operation; if errfunc[2]
CLEARs a zone and then errfunc[3] SETs that zone, the final operation will
be to tag that zone (and vice versa). In practice by default this does not matter, because the
built-in tagging routines never explicitly perform a ``CLEAR``. However,
it is possible to overwrite the Tagging_nd.f90 file if you want to change how
ca_denerror, ca_temperror, etc. operate. This is not recommended, and if you do so
be aware that CLEARing a zone this way may not have the desired effect.

We provide also the ability for the user to define their own tagging criteria.
This is done through the Fortran function in the
files. This function is provided the entire
state (including density, temperature, velocity, etc.) and the array
of tagging status for every zone. As an example of how to use this, suppose we
have a 3D Cartesian simulation where we want to tag any zone that has a
density gradient greater than 10, but we don’t care about any regions
outside a radius :math:`r > 75` from the problem origin; we leave them always unrefined.
We also want to ensure that the region :math:`r \leq 10` is always refined.
In our probin file we would set denerr = 10 and max_denerr_lev = 1
in the &tagging namelist. We would also make a copy of
problem_tagging_3d.f90 to our work directory and set it up as follows:

::

    subroutine set_problem_tags(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                                state,state_l1,state_l2,state_l3, &
                                state_h1,state_h2,state_h3,&
                                set,clear,&
                                lo,hi,&
                                dx,problo,time,level)

      use bl_constants_module, only: ZERO, HALF
      use prob_params_module, only: center
      use meth_params_module, only: URHO, UMX, UMY, UMZ, UEDEN, NVAR
     
      implicit none
      
      integer         ,intent(in   ) :: lo(3),hi(3)
      integer         ,intent(in   ) :: state_l1,state_l2,state_l3, &
                                        state_h1,state_h2,state_h3
      integer         ,intent(in   ) :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      double precision,intent(in   ) :: state(state_l1:state_h1, &
                                              state_l2:state_h2, &
                                              state_l3:state_h3,NVAR)
      integer         ,intent(inout) :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision,intent(in   ) :: problo(3),dx(3),time
      integer         ,intent(in   ) :: level,set,clear

      double precision :: x, y, z, r

      do k = lo(3), hi(3)
         z = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
         do j = lo(2), hi(2)
            y = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
            do i = lo(1), hi(1)
               x = problo(1) + (dble(i) + HALF) * dx(1) - center(2)

               r = (x**2 + y**2 + z**2)**(HALF)

               if (r > 75.0) then
                 tag(i,j,k) = clear
               elseif (r <= 10.0) then
                 tag(i,j,k) = set
               endif
            enddo
         enddo
      enddo
      
    end subroutine set_problem_tags

.. _sec:amr_synchronization:

Synchronization Algorithm
-------------------------

Here we present the AMR algorithm for the compressible equations with
self-gravity. The gravity component of the algorithm is closely
related to (but not identical to) that in Miniati and Colella, JCP,
2007. The content here is largely based on the content in the original
Castro paper (:raw-latex:`\cite{castro_I}`). The most significant difference is the
addition of a different strategy for when to employ the synchronization;
but regardless of whether the original or new strategy is used, the fundamental
synchronization step is identical.

.. _sec:synchronization_methodology:

Synchronization Methodology
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Over a coarse grid time step we collect flux register information for
the hyperbolic part of the synchronization:

.. math:: \delta{\bf F}= -\Delta t_c A^c F^c + \sum \Delta t_f A^f F^f

Analogously, at the end of a coarse grid time step we store the
mismatch in normal gradients of :math:`\phi` at the coarse-fine interface:

.. math::

   \delta F_\phi =  - A^c \frac{\partial \phi^c}{\partial n}
   + \sum A^f \frac{\partial \phi^f}{\partial n}

We want the composite :math:`\phi^{c-f}` to satisfy the multilevel
version of (`[eq:Self Gravity] <#eq:Self Gravity>`__) at the synchronization time, just
as we want the coarse and fine fluxes at that time to match. So the goal
is to synchronize :math:`\phi` across levels at that time and then zero out
this mismatch register.

At the end of a coarse grid time step we can define
:math:`{\overline{{\bf U}}}^{c-f}` and :math:`\overline{\phi}^{c-f}` as the composite
of the data from coarse and fine grids as a provisional solution at
time :math:`n+1`. (Assume :math:`\overline{{\bf U}}` has been averaged down so that
the data on coarse cells underlying fine cells is the average of the
fine cell data above it.)

The synchronization consists of two parts:

-  Step 1: Hyperbolic reflux

   In the hyperbolic reflux step, we update the conserved variables with
   the flux synchronization and adjust the gravitational terms to reflect
   the changes in :math:`\rho` and :math:`{\bf u}`.

   .. math:: {{\bf U}}^{c, \star} = \overline{{\bf U}}^{c} + \frac{\delta{\bf F}}{V},

   where :math:`V` is the volume of the cell and the correction from
   :math:`\delta{\bf F}` is supported only on coarse cells adjacent to fine grids.

   Note: this can be enabled/disabled via castro.do_reflux. Generally,
   it should be enabled (1).

   Also note that for axisymmetric or 1D spherical coordinates, the
   reflux of the pressure gradient is different, since it cannot be
   expressed as a divergence in those geometries. We use a separate
   flux register in the hydro code to store the pressure term in these
   cases.

-  Step 2: Gravitational synchronization

   In this step we correct for the mismatch in normal derivative in
   :math:`\phi^{c-f}` at the coarse-fine interface, as well as accounting for
   the changes in source terms for :math:`(\rho {\bf u})` and :math:`(\rho E)` due to the
   change in :math:`\rho.`

   On the coarse grid only, we define

   .. math:: (\delta \rho)^{c} =  \rho^{c, \star} - {\overline{\rho}}^{c}  .

   We then form the composite residual, which is composed of two
   contributions. The first is the degree to which the current :math:`\overline{\phi}^{c-f}` does not satisfy the original equation on a
   composite grid (since we have solved for :math:`\overline{\phi}^{c-f}`
   separately on the coarse and fine levels). The second is the response
   of :math:`\phi` to the change in :math:`\rho.` We define

   .. math::

      R \equiv  4 \pi G \rho^{\star,c-f} - \Delta^{c-f} \; \overline{\phi}^{c-f} 
      = - 4 \pi G (\delta \rho)^c - (\nabla \cdot \delta F_\phi ) |_c   .

   Then we solve

   .. math::

      \Delta^{c-f} \; \delta \phi^{c-f} = R
      \label{eq:gravsync}

   as a two level solve at the coarse and fine levels.
   We define the update to gravity,

   .. math:: \delta {\bf g}^{c-f} = \nabla (\delta \phi^{c-f})  .

   Finally, we need to

   -  add :math:`\delta \phi^{c-f}` directly to
      to :math:`\phi^{c}` and :math:`\phi^{f}` and interpolate :math:`\delta \phi^{c-f}` to any finer
      levels and add to the current :math:`\phi` at those levels.

   -  if level :math:`c` is not the coarsest level in the calculation, then we must transmit the
      effect of this change in :math:`\phi` to the coarser levels by updating the flux register between
      level :math:`c` and the next coarser level, :math:`cc.` In particular, we set

      .. math::

         \delta {F_\phi}^{cc-c} = \delta F_\phi^{cc-c} 
         + \sum A^c \frac{\partial (\delta \phi)^{c-f}}{\partial n}  .

   The gravity synchronization algorithm can be disabled with
   gravity.no_sync = 1. This should be done with care. Generally,
   it is okay only if he refluxing happens in regions of low density that
   don’t affect the gravity substantially.

.. _sec:synchronization_sources:

Source Terms
~~~~~~~~~~~~

After a synchronization has been applied, the state on the coarse grid
has changed, due to the change in fluxes at the coarse-fine boundary as
well as the change in the gravitational field. This poses a problem
regarding the source terms, all of which generally rely either on the
state itself, or on the global variables affected by the synchronization
such as the gravitational field. The new-time sources constructed on the
coarse grid all depended on what the state was after the coarse-grid
hydrodynamic update, but the synchronization and associated flux
correction step retroactively changed that hydrodynamic update. So one
can imagine that in a perfect world, we would have calculated the
hydrodynamic update first, including the coarse-fine mismatch
corrections, and only then computed the source terms at the new time.
Indeed, an algorithm that did not subcycle, but marched every zone along
at the same timestep, could do so – and some codes, like FLASH,
actually do this, where no new-time source terms are computed on any
level until the hydrodynamic update has been fully completed and the
coarse-fine mismatches corrected. But in Castro we cannot do this; in
general we assume the ability to subcycle, so the architecture is set up
to always calculate the new-time source terms on a given level
immediately after the hydrodynamic update on that level. Hence on the
coarse level we calculate the new-time source terms before any fine grid
timesteps occur.

One way to fix this, as suggested by Miniati and Colella for the case
of gravity, is to explicitly compute what the difference in the source
term is as a result of any flux corrections across coarse-fine
boundaries. They work out the form of this update for the case of a
cell-centered gravitational force, which has contributions from both
the density advected across the coarse-fine interfaces
(i.e. :math:`\delta \rho \mathbf{g}`, where :math:`\delta \rho` is the density
change due to the coarse-fine synchronization on the coarse rid), as
well as the global change in the gravitational field due to the
collective mass motion (see Miniati and Colella for the explicit form
of the source term). This has a couple of severe limitations. First,
it means that when the form of the source term is changed, the form of
the corrector term is changed too. For example, it is less easy to
write down the form of this corrector term for the flux-based
gravitational energy source term that is now standard in Castro.
Second, gravity is a relatively easy case due to its linearity in the
density and the gravitational acceleration; other source terms
representing more complicated physics might not have an easily
expressible representation in terms of the reflux contribution. For
example, for a general nuclear reaction network (that does not have an
analytic solution), it is not possible to write down an analytic
expression for the nuclear reactions that occur because of
:math:`\delta \rho`.

Instead we choose a more general approach. On the coarse level, we save
the new-time source terms that were applied until the end of the fine
timesteps. We also save the fine level new-time source terms. Then, when
we do the AMR synchronization after a fine timestep, we first subtract
the previously applied new-time source terms to both the coarse and the
fine level, then do the flux correction and associated gravitational
sync solve, and then re-compute the new-time source terms on both the
coarse and the fine level [15]_. In this way, we get almost
the ideal behavior – if we aren’t subcycling, then we get essentially
the same state at the end of the fine timestep as we would in a code
that explicitly had no subcycling. The cost is re-computing the new-time
source terms that second time on each level. For most common source
terms such as gravity, this is not a serious problem – the cost of
re-computing :math:`\rho \mathbf{g}` (for example, once you already know
:math:`\mathbf{g}`) is negligible compared to the cost of actually computing
:math:`\mathbf{g}` itself (say, for self-gravity). If you believe that the
error in not recomputing the source terms is sufficiently low, or the
computational cost of computing them too high, you can disable this
behavior [16]_ using the
code parameter castro.update_sources_after_reflux.

Note that at present nuclear reactions are not enabled as part of this
scheme, and at present are not automatically updated after an AMR
synchronization. This will be amended in a future release of Castro.

.. _sec:synchronization_timing:

Synchronization Timing
~~~~~~~~~~~~~~~~~~~~~~

The goal of the synchronization step is for the coarse and fine grid to
match at the end of a coarse timesteps, after all subcycled fine grid
timesteps have been completed and the two levels have reached the same
simulation time. If subcycling is disabled, so that the coarse and fine
grid take the same timestep, then this is sufficient. However, in the
general subcycling case, the situation is more complicated. Consider the
discussion about source terms in `2.2 <#sec:synchronization_sources>`__. If
we have a coarse level and one fine level with a refinement ratio of
two, then for normal subcycling the fine grid takes two timesteps for
every one timestep taken by the coarse level. The strategy advocated by
the original Castro paper (and Miniati and Colella) is to only do the
AMR synchronization at the actual synchronization time between coarse
and fine levels, that is, at the end of the second fine timestep.
Consequently, we actually only update the source terms after that second
fine timestep. Thus note that on the fine grid, only the *new-time*
source terms in the *second* fine timestep are updated. But a
moment’s thought should reveal a limitation of this. The first fine grid
timestep was also responsible for modifying the fluxes on the coarse
grid, but the algorithm as presented above didn’t take full account of
this information. So, the gravitational field at the old time in
the second fine timestep is actually missing information that would have
been present if we had updated the coarse grid already. Is there a way
to use this information? For the assumptions we make in Castro, the
answer is actually yes. If we apply the effect of the synchronization
not at the synchronization time but at the end of every fine
timestep, then every fine timestep always has the most up-to-date
information possible about the state of the gravitational field. Now, of
course, in fine timesteps before the last one, we have not actually
reached the synchronization time. But we already know at the end of the
first fine timestep what the synchronization correction will be from
that fine timestep: it will be equal to 1/2 of the coarse contribution
to the flux register and the normal contribution to the flux register
for just that timestep. This is true because in Castro, we assume that
the fluxes provided by the hydrodynamic solver are piecewise-constant
over the timestep, which is all that is needed to be second-order
accurate in time if the fluxes are time centered [17]_. So it is fair to say
that halfway through the coarse timestep, half of the coarse flux has
been advected, and we can mathematically split the flux register into
two contributions that have equal weighting from the coarse flux. (In
general, of course, the coarse flux contribution at each fine timestep
is weighted by :math:`1/R` where :math:`R` is the refinement ratio between the
coarse and fine levels.) So, there is nothing preventing us from
updating the coarse solution at the synchronization time :math:`t^{n+1}_c`
after this first fine timestep; we already know at that point how the
coarse solution will change, so why not use that information? We can
then update the gravitational potential at :math:`t^{n+1/2}_c` that is used to
construct the boundary conditions for the gravitational potential solve
on the fine grid at the beginning of the second fine timestep.

In practice, this just means calling the synchronization routine
described in `2.1 <#sec:synchronization_methodology>`__, with the only
modification being that the flux register contribution from the coarse
grid is appropriately weighted by the fine grid timestep instead of
the coarse grid timestep, and we only include the current fine step:

.. math:: \delta{\bf F}= -\Delta t_f A^c F^c + \Delta t_f A^f F^f

The form of the :math:`\phi` flux register remains unchanged, because the
intent of the gravity sync solve is to simply instantaneously correct
the mismatch between the fine and coarse grid. The only difference,
then, between the old strategy and this new method is that we call the
synchronization at the end of every fine timestep instead of only the
last subcycled one, and we change the weighting appropriately. This
new method is more expensive as currently implemented because we have
to do :math:`R` gravitational sync solves, refluxes, and source term
recalculations instead of only one. However, it results in maximal
possible accuracy, especially in cases where significant amounts of
material are crossing refinement boundaries. The reflux strategy is
controlled by the parameter castro.reflux_strategy. At present
the old method is still the default.

Note that one does not need to be using self-gravity for this to be
beneficial. Even in pure hydrodynamics this can matter. If a regrid
occurs on the fine level, new zones on the boundaries of the current
fine level are filled by interpolation from the coarse level. In the
old method, that interpolation is not using the most up-to-date data
that accounts for the synchronization.

For multiple levels of refinement, the scheme extends naturally. In
the old method, we always call the synchronization at the
synchronization time between any two levels. So for example with two
jumps in refinement by a factor of two, there is a synchronization at
the end of the first two timesteps on level 2 (between level 1 and
level 2), a synchronization after the next two timesteps on level 2
(again between level 1 and level 2), and then a synchronization
between level 0 and level 1. In the new method, we always call the
synchronization at the end of every timestep *on the finest level
only*, and we simultaneously do the synchronization *on every
level*. The timestep :math:`\Delta t_f` in the flux register is just the
timestep on the finest level. (If this is unclear, give it a sanity
check: when the sum of all flux register totals is added up, the level
0 contribution will have a factor of :math:`\Delta t` equal to the coarse
grid timestep since the sum of the timesteps on the finest level over
the entire advance must equal the level 0 timestep. So, the final
contribution from the flux register is the same as if we had saved up
the flux mismatch until the end of the level 0 timestep.) The
synchronization no longer needs to be called at the end of any coarser
level’s timestep because it will already be up to date as a result of
the synchronizations applied at the end of the fine level timesteps.

ConvertCheckpoint
=================

Within the CASTRO distribution, there is the capability to “grow” a
checkpoint file so that a calculation can be restarted in a larger
domain covered by grid cells a factor of two or four coarser than the
existing coarsest level. Instructions for how to do so are in the
Castro/Util/ConvertCheckpoint/README file and are included here.
Upon restart the existing data in the checkpoint file will be used to
fill the region of the previous computational domain, and the new
regions will be filled by some other means, typically interpolation
from a 1D model file.

Star in Corner (**star_at_center = 0**) 
----------------------------------------

In this section we consider the case where the star (or feature of interest)
is centered at the lower left corner of the domain, e.g. you are modeling only one
quarter of the star in 2D, or an octant of the star in 3D. Then you only want
to grow the domain in the “high side” directions (e.g., to the upper right).

Converting the Checkpoint File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let’s say you have a checkpoint file, *chk00100*, say, with 5 levels of refinement
and a (real) problem domain size :math:`P` and (integer) domain size :math:`D` at level 0.

The inputs file that created this might have contained:

-  **max_step** = 100

-  **amr.max_level** = 5

-  **amr.n_cell** = D D

-  **geometry.prob_lo** = 0 0

-  **geometry.prob_hi** = P P

-  **amr.ref_ratio** = 4 4 4 4

Now let’s suppose that you want to grow the domain by a factor of 8 and cover that
new larger domain with a level that is a factor of 2 coarser than the existing level 0 grids.

#. First, set DIM = in the GNUmakefile, and type “
   make” in the Util/ConvertCheckpoint/ directory. This will
   make an executable from the Embiggen.cpp code.

#. Run the embiggening code as follows:

   ::

       Embiggen2d.Linux.Intel.Intel.ex checkin=chk00100 checkout=newchk00050 ref_ratio=2 grown_factor=8 star_at_center=0

   (Your executable may have a slightly different name depending on the compilers you
   built it with.)

   This will create a new checkpoint directory, called *newchk00050*, that represents a simulation
   with *one* additional level of refinement *coarser* than the previous level 0 grids by
   a factor of ref_ratio (in this case, 2).
   The new domain will be a factor of grown_factor (in this case, 8) larger than the previous domain.

   Note that ref_ratio must be 2 or 4, because those are the only acceptable values of ref_ratio
   in Castro.

   grown_factor can be any reasonable integer; but it’s only been
   tested with 2, 3, 4 and 8. It does not need to be a multiple of 2.

Restarting from a Grown Checkpoint File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You should now be able to restart your calculation using *newchk00050*.

Your inputs file should now contain lines like:

-  **max_step** = 51

-  **amr.restart** = newchk00050

-  **amr.max_level** = 6

-  **amr.n_cell** = 4D 4D

-  **geometry.prob_lo** = 0 0

-  **geometry.prob_hi** = 8P 8P

-  **castro.grown_factor** = 8

-  **castro.star_at_center** = 0

-  **amr.ref_ratio** = 2 4 4 4 4

IMPORTANT:

#. Unlike earlier, you may now set **amr.max_level** to be at most one greater than before,
   but you need not set it that high. For example, you could set **amr.max_level** the same as before
   and you would lose data at the finest refinement level. You may not set **amr.max_level** = 0,
   however, because we have no data at the new level 0 until we average down from the new level 1 after
   the restart.

#. You must set **amr.n_cell** = (**grown_factor** / **ref_ratio**) times the previous
   value of **amr.n_cell**. In this case **amr.n_cell** = (8/2)*D = 4D.

#. You must set **amr.prob_hi** to be a factor of **grown_factor** greater than the previous
   value of **amr.prob_hi**.

#. You must insert the value of **ref_ratio** used in the Embiggen call as the first
   value in the list of **amr.ref_ratio**, since that will now be the refinement ratio between
   the new level 0 and the new level 1.

#. You must set **castro.grown_factor** in your inputs file equal to the value of
   **grown_factor** you used when you called Embiggen*ex so that the CASTRO code knows
   how big the original domain was.

#. Note that if you have run 100 steps at the original level 0, that would be equivalent
   to 50 steps at the new level 0 because you coarsened by a factor of 2.
   Thus once you re-start from the new checkpoint directory,
   the next step will be 51, not 101. Make sure to keep track of your plotfiles accordingly.

#. Don’t forget to adjust max_denerr_lev and comparable variables to control
   the number of fine levels you now want. If you want to have 6 levels of refinement
   after restart, then make sure max_denerr_lev, etc, are set high enough. If you
   only want to have 5 levels of refinement (where the new level 5 would now be
   a factor of ref_ratio coarser than the previous level 5), make sure to adjust
   max_denerr_lev accordingly as well.

Star at Center of Domain (**star_at_center = 1**) 
--------------------------------------------------

Now let’s assume that the star (or feature of interest) is centered at the center of the
domain in 2D or 3D Cartesian coordinates. We will later consider the case of 2D cylindrical (r-z)
coordinates in which the star is centered at the left midpoint.

.. _converting-the-checkpoint-file-1:

Converting the Checkpoint File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose that you want to grow the domain by a factor of 2 and cover that
new larger domain with a level that is a factor of 2 coarser than the existing level 0 grids.

After you build the Embiggen executable, you type:

-  | Embiggen2d.Linux.Intel.Intel.ex **checkin**\ =chk00100 **checkout**\ =newchk00050 **ref_ratio**\ =2
   | **grown_factor**\ =2 **star_at_center**\ =1

Note that

-  **ref_ratio** must still be 2 or 4

-  **grown_factor** can only be 2 or 3 in this case.

.. _restarting-from-a-grown-checkpoint-file-1:

Restarting from a Grown Checkpoint File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Your inputs file for restarting would now look like

-  **max_step** = 51

-  **amr.restart** = newchk00050

-  **amr.max_level** = 6

-  **amr.n_cell** = D D

-  **geometry.prob_lo** = -P/2 -P/2

-  **geometry.prob_hi** = 3P/2 3P/2

-  **castro.grown_factor** = 2

-  **castro.star_at_center** = 1

-  **amr.ref_ratio** = 2 4 4 4 4

Cylindrical Coordinates
~~~~~~~~~~~~~~~~~~~~~~~

In the case of 2D cylindrical (r-z) coordinates in which the star is centered at the left edge
but vertical midpoint of the domain, the embiggening procedure is the same as above
(with **star_at_center = 1**) but the inputs file for restart is slightly different in that
**geometry.prob_lo** is modified in the z- but not the r-direction. If we consider the original
inputs file to look like:

-  **max_step** = 100

-  **amr.max_level** = 6

-  **amr.n_cell** = D 2D

-  **geometry.prob_lo** = 0 0

-  **geometry.prob_hi** = P 2P

-  **amr.ref_ratio** = 4 4 4 4

then an inputs file for restart would look like:

-  **amr.restart** = newchk00050

-  **amr.max_level** = 6

-  **amr.n_cell** = D 2D

-  **geometry.prob_lo** = 0 -P

-  **geometry.prob_hi** = 2P 3P

-  **castro.grown_factor** = 2

-  **castro.star_at_center** = 1

-  **amr.ref_ratio** = 2 4 4 4 4

.. raw:: latex

   \centering

|Data from checkpoint file before and after the domain has been coarsened and grown. This case
uses **star_at_center = 0** and **ref_ratio**\ =2. The first grown example has
**grown_factor**\ =2, the second has **grown_factor**\ =3. In all figures the level 0 grids
are shown in white, the level 1 grids in red, the level 2 grids in yellow, and in the grown figures,
the level 3 grids are in pink.|
|Data from checkpoint file before and after the domain has been coarsened and grown. This case
uses **star_at_center = 0** and **ref_ratio**\ =2. The first grown example has
**grown_factor**\ =2, the second has **grown_factor**\ =3. In all figures the level 0 grids
are shown in white, the level 1 grids in red, the level 2 grids in yellow, and in the grown figures,
the level 3 grids are in pink.|
|Data from checkpoint file before and after the domain has been coarsened and grown. This case
uses **star_at_center = 0** and **ref_ratio**\ =2. The first grown example has
**grown_factor**\ =2, the second has **grown_factor**\ =3. In all figures the level 0 grids
are shown in white, the level 1 grids in red, the level 2 grids in yellow, and in the grown figures,
the level 3 grids are in pink.|

.. raw:: latex

   \centering

|Data from checkpoint file before and after the domain has been coarsened and grown. This case
uses **star_at_center = 0** and **ref_ratio**\ =2. The first grown example has
**grown_factor**\ =2, the second has **grown_factor**\ =3. In all figures the level 0 grids
are shown in white, the level 1 grids in red, the level 2 grids in yellow, and in the grown figure,
the level 3 grids are in pink. |
|Data from checkpoint file before and after the domain has been coarsened and grown. This case
uses **star_at_center = 0** and **ref_ratio**\ =2. The first grown example has
**grown_factor**\ =2, the second has **grown_factor**\ =3. In all figures the level 0 grids
are shown in white, the level 1 grids in red, the level 2 grids in yellow, and in the grown figure,
the level 3 grids are in pink. |
|Data from checkpoint file before and after the domain has been coarsened and grown. This case
uses **star_at_center = 0** and **ref_ratio**\ =2. The first grown example has
**grown_factor**\ =2, the second has **grown_factor**\ =3. In all figures the level 0 grids
are shown in white, the level 1 grids in red, the level 2 grids in yellow, and in the grown figure,
the level 3 grids are in pink. |

Initializing Castro with Maestro Data
=====================================

Overview
--------

We can now initialize a Castro simulation using data from a Maestro plotfile. This should not be thought of as a restart mode, but rather
a new simulation with a special initialization. In order to use this
feature, you must make sure the Maestro plotfile has the proper
variables, add some new parameters to your inputs file, and add a few
subroutines to Prob_Xd.f90. You need to build a special executable
with “USE_MAESTRO_INIT=TRUE”, which will add “.MAESTRO” to the
executable string. For multilevel problems, there are a few extra
steps relating to the fact that you have to supply a grids file
consistent with the Maestro grid structure.

MAESTRO Plotfile Requirements
-----------------------------

The Maestro plotfile needs to have the following variables:

-  “x_vel”, “y_vel”, (and “z_vel”, depending on
   dimensionality of the problem)

-  “density” (**castro.MAESTRO_init_type** = 1 and 2 only)

-  Optional species (such as “X(C12)”) - there is an option to
   not read any species from the Maestro plotfile. In this case, you
   must make sure your code manually defines the species cell-by-cell
   in the initial Castro data

-  “tfromp”

-  “pi” (**castro.MAESTRO_init_type** = 2, 3, and 4 only)

-  “entropy” (**castro.MAESTRO_init_type** = 4 only)

Also, model_cc_XXXXX needs to list variables in the following order,
which is the default order found in MAESTRO/Source/base_io.f90: r,
base_r, rho0, p0, gamma1bar, rhoh0, div_coeff, psi, tempbar,
etarho_cc, tempbar_init.

List of Parameters
------------------

Here are the additional parameters you must add to your inputs file.

+-----------------+-----------------+-----------------+-----------------+
| Parameter       | Definition      | Type            | Default         |
+=================+=================+=================+=================+
| **castro.MAESTR | name of the     | std::string     | must be set     |
| O_plotfile**    | Maestro plotfil |                 |                 |
|                 | e               |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| **castro.MAESTR | name of the     | std::string     | must be set     |
| O_modelfile**   | Maestro “model_ |                 |                 |
|                 | cc”             |                 |                 |
|                 | file            |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| **castro.MAESTR | number of       | int             | must be set     |
| O_npts_model**  | points in the   |                 |                 |
|                 | Maestro model_c |                 |                 |
|                 | c               |                 |                 |
|                 | file            |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| **castro.MAESTR | name of the     | std::string     | must be set or  |
| O_first_species | first species   |                 | else nothing    |
| **              |                 |                 | will be read in |
+-----------------+-----------------+-----------------+-----------------+
| **castro.MAESTR | number of       | std::string     | NumSpec in      |
| O_nspec**       | species in the  |                 | Castro          |
|                 | Maestro plotfil |                 |                 |
|                 | e               |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| **castro.MAESTR | controls how we | Real            | must be set     |
| O_cutoff_densit | overwrite data  |                 |                 |
| y**             | at the edge of  |                 |                 |
|                 | the star        |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| **castro.MAESTR | determines how  | int             | must be set     |
| O_init_type**   | we initialize   |                 |                 |
|                 | the             |                 |                 |
|                 | Castro state    |                 |                 |
+-----------------+-----------------+-----------------+-----------------+
| **castro.MAESTR | specifies       | int             | must be set     |
| O_spherical**   | planar or       |                 |                 |
|                 | spherical       |                 |                 |
|                 | problem         |                 |                 |
+-----------------+-----------------+-----------------+-----------------+

Examples of Usage
~~~~~~~~~~~~~~~~~

-  **castro.MAESTRO_plotfile** = "wd_384_6.25e8K_norotate_plt120578"

-  **castro.MAESTRO_modelfile** = "./wd_384_6.25e8K_norotate_plt120578/model_cc_120578"

-  | **castro.MAESTRO_npts_model** = 1663
   | This is the number of
     points in **castro.MAESTRO_modelfile**. Note that this is not
     the same thing as “npts_model”, which is the number of points in
     the initial model file used for standard simulations where we do not
     initialize from a Maestro plotfile.

-  **castro.MAESTRO_first_species** = “X(C12)” If you do not
   specify this, no species will be read in. You can always manually
   specify or overwrite the species cell-by-cell later.

-  | **castro.MAESTRO_nspec** = 3
   | If you do not specify this, it
     will default to the number of species in the Castro network,
     “NumSpec”. We have this here because sometimes Maestro and Castro will use different networks with different number of species.

-  | **castro.MAESTRO_cutoff_density** = 1.e6
   | The code will use
     this density to figure out the radial coordinate, r_model_start,
     which is the last radial coordinate before rho0 falls below
     **castro.MAESTRO_cutoff_density**. It is possible to set
     **castro.MAESTRO_cutoff_density** to a tiny value, such that rho0
     never falls below this value, in which case we set r_model_start
     to :math:`\infty`. In INITDATA_MAKEMODEL, we create a new 1D model
     integrating outward starting from r_model_start. Then, in
     INITDATA_OVERWRITE, we overwrite newly initialized Castro data in
     any cell that maps into a radial coordinate greater than
     r_model_start by interpolating from the new 1D model.

-  | **castro.MAESTRO_init_type** = 2
   | Castro will read in data
     from the Maestro plotfile, and then call the EOS to make sure that
     :math:`\rho`, :math:`e`, :math:`T`, and :math:`X_k` are consistent. The inputs to the EOS
     are based on the value of **castro.MAESTRO_init_type**:

   #. :math:`e = e(\rho,T,X_k)`

   #. :math:`e,T = e,T(\rho,p_0+\pi,X_k)`

   #. :math:`\rho,e = \rho,e(p_0+\pi,T,X_k)`

   #. :math:`\rho,T,e = \rho,T,e(p_0+\pi,s,X_k)`

-  | **castro.MAESTRO_spherical** = 1
   | 0 = planar; 1 = spherical.

New Subroutines in Prob_Xd.f90
------------------------------

There are three routines that need to be added to your local copy of
Prob_Xd.f90. See Castro/Exec/wdconvect/Prob_3d.f90 for
a standard spherical Maestro initialization.

#. | INITDATA_MAESTRO
   | This fills in the Castro state by taking
     the Maestro data, calling the EOS, and making the proper variables
     conserved quantities. Specifically, we need a thermodynamically
     consistent :math:`\rho`, :math:`T`, :math:`e`, and :math:`X_k`, and then algebraically
     compute :math:`\rho{\bf u}`, :math:`\rho e`, :math:`\rho E`, and :math:`\rho X_k`,

#. | INITDATA_MAKEMODEL
   | This creates a user-defined 1D initial model starting from r_model_start.

#. | INITDATA_OVERWRITE
   | This overwrites the initialized Castro data using the new 1D initial model for all cells that map into
     radial coordinates greater than r_model_start.

Additional Notes
----------------

Note that for both single-level and multilevel Maestro to Castro initialization, the Castro base grid structure does not have to match
the Maestro base grid structure, as long as the problem domain is the
same. For example, if the coarsest level in a Maestro plotfile
contains :math:`64^3` cells divided into 8-\ :math:`32^3` grids, it is ok to use a
Castro base grid structure with 1-\ :math:`64^3` grid, 64-\ :math:`16^3` grids, or
anything else you can imagine - the grids don’t even have to be the
same size. As is normally the case, the Castro base grid structure is
created based on the parameters in the Castro inputs file, such as
**amr.max_grid_size**, **amr.blocking_factor**, etc.

Multilevel Restart
~~~~~~~~~~~~~~~~~~

When initialing from a multilevel Maestro plotfile, there are some
extra steps. First, you need to create a Castro-compatible grids file
from the Maestro plotfile. This can be done with the
BoxLib/Tools/Postprocessing/F_Src/fboxinfo.f90 utility. Compile
and run this using the “``–``\ castro” option, e.g.,
“fboxinfo.Linux.gfortran.exe ``–``\ castro pltxxxxx ``|``
tee gr0.maestro”, to generate the Castro-compatible grids file. Note
that the base grid structure is still controlled by
**amr.max_grid_size**, **amr.blocking_factor**, etc., since in C BoxLib, the grids file only indicates the refined grid structure,
whereas in Fortran BoxLib the grids file contains the base grid and
refined grid structures.

Now, when you initialize the Castro simulation, you need to specify
the grid file using **amr.regrid_file = "gr0_3d.128_2levels"**,
for example. You can happily run this now, but note that the
regridding algorithm will never be called (since Castro thinks it’s
started a new simulation from scratch with a grids file, thus
disabling the regridding). If you wish for the grid structure to be
changed, you must do a traditional Castro restart from the
Castro-generated checkpoint file (you can still use the same
“.MAESTRO” executable or an executable built with
USE_MAESTRO_INIT=FALSE), making sure that you **do not** specity
**amr.regrid_file** (or else the grids will stay fixed). You are
free to specify **amr.regrid_on_restart**,
**amr.compute_new_dt_on_regrid**, and
**amr.plotfile_on_restart**.

Sometimes a Maestro plotfile will only have 1 or 2 total levels, but
you ultimately want to run a Castro simulation with many more levels
of refinement. My recommended strategy is the following:

#. Initialize a Castro simulation from the Maestro plotfile
   while preserving the exact same grid structure and run for 10 time
   steps.

#. Do a traditional Castro restart from chk00010, but do not
   increase **amr.max_level**, and run for 10 more time steps. This
   allows a new grid structure with the same effective resolution as
   before settle in using the C BoxLib regridding algorithm.

#. Do a traditional Castro restart from chk00020, but increase
   **amr.max_level** by 1, and run for 10 time steps.

#. Repeat the procedure from the previous step (using the most
   updated checkpoint of course) as many times as desired.

Visualization
=============

There are a large number of tools that can be used to read in Castro or BoxLib data and make plots. Here we give a brief overview of some
of the tools as well as some examples.

Controlling What’s in the PlotFile
----------------------------------

There are a few options that can be set at runtime to control what
variables appear in the plotfile.

-  : this controls which of the main
   state variables appear in the plotfile. The default is for all of
   them to be stored. But you can specify a subset by name, e.g.

   ::

             amr.plot_vars = density
           

   to only store that subset.

-  : this controls which of the
   derived variables to be stored in the plotfile. Derived variables
   are created only when the plotfile is being created, using the
   infrastructure provided by BoxLib to register variables and the
   associated Fortran routine to do the deriving
   ().

   By default, no derived variables are stored. You can store all
   derived variables that Castro knows about by doing:

   ::

             amr.derive_plot_vars = ALL
           

   or a subset by explicitly listing them, e.g.,

   ::

             amr.derive_plot_vars = entropy pressure
           

amrvis
------

Our favorite visualization tool is amrvis. We heartily encourage you
to build the amrvis2d and amrvis3d executables, and to try using them
to visualize your data. A very useful feature is View/Dataset, which
allows you to actually view the numbers – this can be handy for
debugging. You can modify how many levels of data you want to see,
whether you want to see the grid boxes or not, what palette you use,
etc.

If you like to have amrvis display a certain variable, at a certain
scale, when you first bring up each plotfile (you can always change it
once the amrvis window is open), you can modify the amrvis.defaults
file in your directory to have amrvis default to these settings every
time you run it. The directories CoreCollapse, HSE_test, Sod and
Sedov have amrvis.defaults files in them. If you are working in a new
run directory, simply copy one of these and modify it.

VisIt
-----

VisIt is also a great visualization tool, and it directly handles our
plotfile format (which it calls Boxlib). For more information check
out visit.llnl.gov.

[Useful tip:] To use the Boxlib3D plugin, select it from File
:math:`\rightarrow` Open file :math:`\rightarrow` Open file as type Boxlib, and
then the key is to read the Header file, plt00000/Header, for example,
rather than telling to to read plt00000.

yt
--

yt is a free and open-source software that provides data analysis and
publication-level visualization tools for astrophysical simulation
results such as those CASTRO produces. As yt is script-based, it’s not
as easy to use as VisIt, and certainly not as easy as amrvis, but the
images can be worth it! Here we do not flesh out yt, but give an
overview intended to get a person started. Full documentation and
explanations from which this section was adapted can be found at
http://yt-project.org/doc/index.html.

yt can be installed by the following commands:

$ wget https://raw.githubusercontent.com/yt-project/yt/master/doc/install_script.sh

$ bash install_script.sh

This installs yt in your current directory. To update ytin the
future, simply do

$ conda update yt

assuming you have conda installed.

Castro-Specific Data
~~~~~~~~~~~~~~~~~~~~

yt was originally created for simple analysis and visualization of
data from the Enzo code. Since, it has grown to include support for a
variety of codes, including Castro. However, ytwill still sometimes
make assumptions, especially about data field names, that favor Enzo
and cause errors with Castro data. These problems can usually be
avoided by taking care to specify the data fields desired in
visualization. For example, Enzo’s density field is called
“Density,” and is the default for many plotting mechanisms when the
user does not specify the field. However, Castro does not have a field
called “Density”; instead, the density field is called “density.”
If a user does not specify a field while plotting with Castro data,
chances are that yt will try (and fail) to find “Density” and return
an error. As you will see in the examples, however, there is a way to
create your own fields from existing ones. You can use these derived
fields as you would use any other field.

There are also a few imperatives when it comes to reading in your
Castro simulation data and associated information. First and foremost
is that the inputs file for the simulation **must** exist in the
same directory as where the plotfile directory is located, and it
**must** be named “**inputs**.” yt reads information from the
inputs file such as the number of levels in the simulation run, the
number of cells, the domain dimensions, and the simulation time. yt will also optionally parse the probin file for pertinent information
if it is similarly included with the name “**probin**” in the same
directory as the plotfile of interest. When specifying a plotfile as
the data source for plots, you may simply call it by its directory
name, rather than using the Header file as in VisIt. As a final
caveat, the existence of the job_info file within the plotfile
directory is what currently distinguishes Castro data from MAESTRO
data in yt; unless you like surprises, we suggest you ensure your
plotfile has one.

Interacting with yt: Command Line and Scripting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ytis written completely in python (if you don’t have python, yt will
install it for you) and there are a number of different ways to
interact with it, including a web-based gui. Here we will cover
command-line yt and scripts/the python interactive prompt, but other
methods are outlined on the yt webpage at
http://yt-project.org/doc/interacting/index.html.

The first step in starting up yt is to activate the yt environment:

$ source $YT_DEST/bin/activate

From the command line you can create simple plots, perform simple
volume renderings, print the statistics of a field for your data set,
and do a few other things. Try $ yt to see a list of commands,
and $ yt :math:`<`\ command\ :math:`>` --help
to see the details of a command. The command line is the easiest way
to get quick, preliminary plots – but the simplicity comes at a
price, as yt will make certain assumptions for you. We could plot a
projection of density along the x-axis for the plotfile (yt calls it a
parameter file) plt_def_00020 by doing the following:

$ yt plot -p -a 0 -f density plt_def_00020

Or a temperature-based volume rendering with 14 contours:

$ yt render -f Temp --contours 14 plt_def_00020

Any plots created from the command line will be saved into a
subdirectory called “frames.” The command line is nice for fast
visualization without immersing yourself too much in the data, but
usually you’ll want to specify and control more details about your
plots. This can be done either through scripts or the python
interactive prompt. You can combine the two by running scripts within
the interactive prompt by the command

:math:`>>>` execfile(‘script.py’)

which will leave you in the interactive prompt, allowing you to
explore the data objects you’ve created in your script and debug
errors you may encounter. While in the yt environment, you can access
the interactive prompt by $ *python* or the shortcut

$ pyyt

Once you’re in the yt environment and in a .py script or the
interactive prompt, there are a couple of points to know about the
general layout of yt scripting. Usually there are five sections to a
yt script:

#. Import modules

#. Load parameter files and saved objects

#. Define variables

#. Create and modify data objects, image arrays, plots,
   etc. :math:`\rightarrow` this is the meat of the script

#. Save images and objects

Note that neither saving nor loading objects is necessary, but can be
useful when the creation of these objects is time-consuming, which is
often the case during identification of clumps or contours.

yt Basics
~~~~~~~~~

The first thing you will always want to do is to import yt:

:math:`>>>` from yt.mods import \*

Under certain circumstances you will be required to import more, as we
will see in some of the examples, but this covers most of it,
including all of the primary functions and data objects provided by
yt. Next, you’ll need yt to access the plotfile you’re interested in
analyzing. Remember, you must have the “inputs” file in the same
directory:

:math:`>>>` ds = load(‘plt_def_00020’)

When this line is executed, it will print out some key parameters from
the simulation. However, in order to access information about all of
the fluid quantities in the simulation, we must use the “index”
object. (Note that for yt versions earlier than 3.0, this information
was contained in the “hierarchy” object; for these versions, replace
pf.index with pf. h in the following examples. The “hierarchy” object
was removed in yt-3.0 and its associated functionality for accessing data
was moved directly to the datasets themselves.) It contains the geometry
of the grid zones, their parentage relationships, and the fluid states
within each one. It is easily created:

:math:`>>>` ds.index

Upon execution, yt may print out a number of lines saying it’s adding
unknown fields to the list of fields. This is because Castro has
different names for fields than what yt expects. We can see what
fields exist through the commands

:math:`>>>` print ds.index.field_list

:math:`>>>` print ds.index.derived_field_list

There may not be any derived fields for Castro data. We can find out
the number of grids and cells at each level, the simulation time, and
information about the finest resolution cells:

:math:`>>>` ds.index.print_stats()

The dataset itself also stores a number of associated methods; for example,
you can find the value and location of the maximum of a field in the domain:

:math:`>>>` value, location = ds.find_max(‘density’)

(Note that in yt versions before 3.0, this type of method was primarily
associated with the hierarchy object and was accessed with ds.h.find_max.)

The list goes on. A full list of methods and attributes associated
with the index object (and most any yt object or function) can be
accessed by the help function:

:math:`>>>` help(pf.index)

You can also use :math:`>>>` *dir()* on an object or
function to find out which names it defines. Don’t be shy about
searching the yt documentation for help. Note that creating the
index object in its own line is not always needed before calling
functions like find_max; yt will construct it automatically if it
does not already exist.

Data Containers and Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes, you’ll want to select, analyze, or plot only portions of
your simulation data. To that end, yt includes a way to create data
“containers” that select data based on geometric bounds or fluid
quantity values. There are many, including rays, cylinders, and clumps
(some in the examples, all described in the documentation), but the
easiest to create is a sphere, centered on the location of the maximum
density cell we found above:

:math:`>>>` my_data_container = ds.sphere(location, (5.0e4, ‘km’))

Here, specify that the radius is in units of kilometers using a dimensionful
quantity. When specifying distances in yt, the default is to use the
simulation-native unit named “code_length”, which for Castro is “cm”, and
if you just put in 5.0e4 instead of (5.0e4, ‘km’), you will get a 50,000 cm radius.
The pf.index.print_stats() command lists available units. We can access the data
within the container:

:math:`>>>` print my_data_container[‘density’]

:math:`>>>` print my_data_container.quantities[‘Extrema’]([‘density’, ‘pressure’])

When the creation of objects is time-consuming, it can be convenient
to save objects so they can be used in another session. To save an
object as part of the .yt file affiliated with the index:

:math:`>>>` pf.index.save_object(my_data_container, ‘sphere_to_analyze_later’)

Once it has been saved, it can be easily loaded later:

:math:`>>>` sphere_to_analyze = pf.index.load_object(‘sphere_to_analyze_later’)

Grid Inspection
~~~~~~~~~~~~~~~

yt also allows for detailed grid inspection. The index object
possesses an array of grids, from which we can select and examine
specific ones:

:math:`>>>` print pf.index.grids

:math:`>>>` my_grid = pf.index.grids[4]

Each grid is a data object that carries information about its
location, parentage relationships (grids within which it resides, and
grids that reside within it, at least in part), fluid quantities, and
more. Here are some of the commands:

:math:`>>>` print my_grid.Level

:math:`>>>` print my_grid_ActiveDimensions

:math:`>>>` print my_grid.LeftEdge

:math:`>>>` print my_grid.RightEdge

:math:`>>>` print my_grid.dds

(dds is the size of each cell within the grid).

:math:`>>>` print my_grid.Parent

:math:`>>>` print my_grid.Children[2].LeftEdge

:math:`>>>` print my_grid[‘Density’]

You can examine which cells within the grid have been refined with the
child_mask attribute, a representative array set to zero everywhere
there is finer resolution.To find the fraction of your grid that isn’t
further refined:

:math:`>>>`\ print my_grid.child_mask.sum()/float(my_grid.ActiveDimensions.prod())

Rather than go into detail about the many possibilities for plotting
in yt, we’ll provide some examples.

Example Scripts
~~~~~~~~~~~~~~~

In these examples, we investigate 3-D simulation data of two stars
orbiting in the center of the domain, which is a box of sides
:math:`10^{10}\:cm`.

*# Pressure Contours*

.. raw:: latex

   \setlength{\parskip}{0pt}

from yt.mods import \*

pf = load(‘plt00020’)

field = ‘pressure’

pf.index

*# Most Castro fields have no inherent units, so we add them in,
in the form of a raw string*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# with some LaTeX-style formatting.*

pf.field_info[field]._units = r‘\\rm{Ba}’

*# SlicePlot parameters include: parameter file, axis, field, window width (effectively the*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# x and y zoom), and fontsize. We can also create projections with ProjectionPlot().*

p = SlicePlot(pf, ‘z’, field, width=((5.0e9, ‘cm’), (3.0e9, ‘cm’)),

fontsize=13)

*# Zlim is the range of the colorbar. In other words, the range of the data we want to display.*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# Names for many colormaps can be found at wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps.*

p.set_zlim(field, 2.85e13, 2.95e13)

p.set_cmap(field, ‘jet’)

*# Here we add 5 density contour lines within certain limits on top of the image. We overlay*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# our finest grids with a transparency of 0.2 (lower is more transparent). We add a quiver*

*# plot with arrows every 16 pixels with x_velocity in the x-direction and y_velocity in*

*# the y-direction. We also mark the center with an ‘x’ and label one of our stars.*

p.annotate_contour(‘density’, clim=(1.05e-4, 1.16e-4), ncont=5, label=False)

p.annotate_grids(alpha=0.2, min_level=2)

p.annotate_quiver(‘x_velocity’, ‘y_velocity’, factor=16)

p.annotate_marker([5.0e9, 5.0e9], marker=‘x’)

p.annotate_point([5.95e9, 5.1e9], ‘Star!’)

*# This saves the plot to a file with the given prefix. We can alternatively specify*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# the entire filename.*

p.save(‘contours.press_den\_’)

.. raw:: latex

   \centering

.. figure:: Slice_z_pressure
   :alt: Pressure slice with annotations
   :width: 6in

   Pressure slice with annotations

*#————————*

*# Volume Rendering*

.. raw:: latex

   \setlength{\parskip}{0pt}

from yt.mods import \*

pf = load(‘plt00020’)

field = ‘pressure’
dd = pf.all_data()

*# We take the log of the extrema of the pressure field, as well as a couple other interesting*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# value ranges we’d like to visualize.*

h_mi, h_ma = dd.quantities[‘Extrema’](field)[0]

h_mi, h_ma = np.log10(h_mi), np.log10(h_ma)

s_mi, s_ma = np.log10(2.90e13), np.log10(3.10e13)

pf.index

*# We deal in terms of logarithms here because we have such a large range of values.*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# It can make things easier, but is not necessary.*

pf.field_info[field].take_log=True

*# This is what we use to visualize volumes. There are a couple of other, more complex*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# ways. We set the range of values we’re interested in and the number of bins in the*

*# function. Make sure to have a lot of bins if your data spans many orders of magnitude!*

*# Our raw data ranges from about :math:`10^{13}` to :math:`10^{22}`.*

tf = ColorTransferFunction((h_mi-1, h_ma+1), nbins=1.0e6)

*# Here we add several layers to our function, either one at a time or in groups. We*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# specify the value-center and width of the layer. We can manipulate the color by*

*# individually setting the colormaps and ranges to spread them over. We can also*

*# change the transparency, which will usually take some time to get perfect.*

tf.sample_colormap(np.log10(2.0e21), 0.006, col_bounds=[h_mi,h_ma],

alpha=[27.0], colormap=‘RdBu_r’)

tf.sample_colormap(np.log10(2.0e19), 0.001, col_bounds=[h_mi,h_ma],

alpha=[5.5], colormap=‘RdBu_r’)

tf.add_layers(6, mi=np.log10(2.95e13), ma=s_ma,

col_bounds=[s_mi,s_ma],

alpha=19*na.ones(6,dtype=‘float64’), colormap=‘RdBu_r’)

tf.sample_colormap(np.log10(2.95e13), 0.000005, col_bounds=[s_mi,s_ma],

alpha=[13.0], colormap=‘RdBu_r’)

tf.sample_colormap(np.log10(2.90e13), 0.000007, col_bounds=[s_mi,s_ma],

alpha=[11.5], colormap=‘RdBu_r’)

tf.sample_colormap(np.log10(2.85e13), 0.000008, col_bounds=[s_mi,s_ma],

alpha=[9.5], colormap=‘RdBu_r’)

*# By default each color channel is only opaque to itself. If we set grey_opacity=True,*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# this is no longer the case. This is good to use if we want to obscure the inner*

*# portions of our rendering. Here it only makes a minor change, as we must set our*

*# alpha values for the outer layers higher to see a strong effect.*

tf.grey_opacity=True

*# Volume rendering uses a camera object which centers the view at the coordinates we’ve*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# called ‘c.’ ‘L’ is the normal vector (automatically normalized) between the camera*

*# position and ‘c,’ and ‘W’ determines the width of the image—again, like a zoom.*

*# ‘Nvec’ is the number of pixels in the x and y directions, so it determines the actual*

*# size of the image.*

c = [5.0e9, 5.0e9, 5.0e9]

L = [0.15, 1.0, 0.40]

W = (pf.domain_right_edge - pf.domain_left_edge)*0.5

Nvec = 768

*# ‘no_ghost’ is an optimization option that can speed up calculations greatly, but can*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# also create artifacts at grid edges and affect smoothness. For our data, there is no*

*# speed difference, so we opt for a better-looking image.*

cam = pf.camera(c, L, W, (Nvec,Nvec), transfer_function = tf,

fields=[field], pf=pf, no_ghost=False)

*# Obtain an image! However, we’ll want to annotate it with some other things before*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# saving it.*

im = cam.snapshot()

*# Here we draw a box around our stars, and visualize the gridding of the top two levels.*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# Note that draw_grids returns a new image while draw_box does not. Also, add\_*

*# background_color in front of draw_box is necessary to make the box appear over*

*# blank space (draw_grids calls this internally). For draw_box we specify the left*

*# (lower) and right(upper) bounds as well its color and transparency.*

im.add_background_color(‘black’, inline=True)

cam.draw_box(im, np.array([3.0e9, 4.0e9, 4.0e9]),

np.array([7.0e9, 6.0e9, 6.0e9]), np.array([1.0, 1.0, 1.0, 0.14]))

im = cam.draw_grids(im, alpha=0.12, min_level=2)

im = cam.draw_grids(im, alpha=0.03, min_level=1, max_level=1)

*# ‘im’ is an image array rather than a plot object, so we save it using a different*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# function. There are others, such as ‘write_bitmap.’*

im.write_png(‘pressure_shell_volume.png’)

.. raw:: latex

   \centering

.. figure:: volume
   :alt: Volume rendering
   :width: 3.5in

   Volume rendering

*#————————*

*# Isocontour Rendering*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# Here we extract isocontours using some extra modules and plot them using matplotlib.*

from mpl_toolkits.mplot3d import Axes3D

from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import matplotlib.pyplot as plt

from yt.mods import \*

pf = load(‘plt00020’)

field = ‘pressure’

field_weight = ‘magvel’

contour_value = 2.83e13

domain = pf.all_data()

*# This object identifies isocontours at a given value for a given field. It returns*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# the vertices of the triangles in that isocontour. It requires a data source, which*

*# can be an object—but here we just give it all of our data. Here we find a pressure*

*# isocontour and color it the magnitude of velocity over the same contour.*

surface = pf.surface(domain, field, contour_value)

colors = apply_colormap(np.log10(surface[field_weight]), cmap_name=‘RdBu’)

fig = plt.figure()

ax = fig.gca(projection=‘3d’)

p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)

p3dc.set_facecolors(colors[0,:,:]/255.)

ax.add_collection(p3dc)

*# By setting the scaling on the plot to be the same in all directions (using the x scale),*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# we ensure that no warping or stretching of the data occurs.*

ax.auto_scale_xyz(surface.vertices[0,:], surface.vertices[0,:],

surface.vertices[0,:])

ax.set_aspect(1.0)

plt.savefig(‘pres_magvel_isocontours.png’)

.. raw:: latex

   \centering

.. figure:: isocontours
   :alt: Pressure isocontour rendering colored with velocity magnitude
   :width: 4in

   Pressure isocontour rendering colored with velocity magnitude

*#————————*

*#1-D and 2-D Profiles*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# Line plots and phase plots can be useful for analyzing data in detail.*

from yt.mods import \*

pf = load(‘plt00020’)

pf.index

*# Just like with the pressure_contours script, we can set the units for fields that*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# have none.*

pf.field_info[‘magvel’]._units = r‘\\rm{cm}/\rm{s}’

pf.field_info[‘kineng’]._units = r‘\\rm{ergs}’

*# We can create new fields from existing ones. ytassumes all units are in cgs, and*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# does not do any unit conversions on its own (but we can make it). Creating new fields*

*# requires us to define a function that acts on our data and returns the new data,*

*# then call add_field while supplying the field name, the function the data comes from,*

*# and the units. Here, we create new fields simply to rename our data to make the plot*

*# look prettier.*

def \_newT(field, data):

return data[‘t’]

add_field(‘X’, function=_newT, units=r‘\\rm{domain} \rm{fraction}’)

def \_newDen(field, data):

return data[‘density’]

add_field(‘Density’, function=_newDen, units=r‘\\rm{g}/\rm{cm}^{3}’)

*# PlotCollections are one of the most commonly used tools in yt, alongside SlicePlots and*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# ProjectionPlots. They are useful when we want to create multiple plots from the same*

*# parameter file, linked by common characteristics such as the colormap, its bounds, and*

*# the image width. It is easy to create 1-D line plots and 2-D phase plots through a*

*# PlotCollection, but we can also create thin projections and so on. When we create a*

*# PlotCollection, it is empty, and only requires the parameter file and the ’center’ that*

*# will be supplied to plots like slices and sphere plots.*

pc = PlotCollection(pf, ‘c’)

*# Now we add a ray—a sample of our data field along a line between two points we define*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# in the function call.*

ray = pc.add_ray([0.0, 5.0e9, 5.0e9],[1.e10, 5.0e9, 5.0e9], ‘magvel’)

*# This is where our derived fields come in handy. Our ray is drawn along the x-axis*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# through the center of the domain, but by default the fraction of the ray we have gone*

*# along is called ‘t.’ We now have the same data in another field we called ‘X,’ whose*

*# name makes more sense, so we’ll reassign the ray’s first field to be that. If we wanted,*

(*# we could also reassign names to ‘magvel’ and ‘kineng.’*

ray.fields = [‘X’, ‘magvel’]

*# Next, we’ll create a phase plot. The function requires a data source, and we can’t*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# just hand it our parameter file, but as a substitute we can quickly create an object*

*# that spans our entire domain (or use the method in the isocontour example). The*

*# specifications of the region (a box) are the center, left bound, and right bound.*

region = pf.region([5.0e9, 5.0e9, 5.0e9], [0.0, 0.0, 0.0],

[1.0e10, 1.0e10, 1.0e10])

*# The phase object accepts a data source, fields, a weight, a number of bins along both*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# axes, and several other things, including its own colormap, logarithm options,*

*# normalization options, and an accumulation option. The first field is binned onto*

*# the x-axis, the second field is binned onto the y-axis, and the third field is*

*# binned with the colormap onto the other two. Subsequent fields go into an underlying*

*# profile and do not appear on the image.*

phase = pc.add_phase_object(region, [‘Density’, ‘magvel’,‘kineng’], weight=None,

x_bins=288, y_bins=288)

pc.save(‘profile’)

.. raw:: latex

   \centering

.. figure:: LineQueryPlot_0_t_magvel
   :alt: Density/velocity magnitude/kinetic energy phase plot
   :width: 4in

   Density/velocity magnitude/kinetic energy phase plot

.. figure:: Profile2D_1_Density_magvel_kineng
   :alt: Density/velocity magnitude/kinetic energy phase plot
   :width: 4in

   Density/velocity magnitude/kinetic energy phase plot

.. raw:: latex

   \quad

*#————————*

*#Off-Axis Projection*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# If we don’t want to take a projection (this can be done for a slice as well) along*

*# one of the coordinate axes, we can take one from any direction using an*

*# OffAxisProjectionPlot. To accomplish the task of setting the view up, the plot*

*# requires some of the same parameters as the camera object: a normal vector, center,*

*# width, and field, and optionally we can set no_ghost (default is False). The normal*

*# vector is automatically normalized as in the case of the camera. The plot also*

*# requires a depth—that is, how much data we want to sample along the line of sight,*

*# centered around the center. In this case ‘c’ is a shortcut for the domain center.*

pf = load(‘plt00020’)

field = ‘density’

L = [0.25, 0.9, 0.40]

plot = OffAxisProjectionPlot(pf, L, field, center=‘c’,

width=(5.0e9, 4.0e9), depth=3.0e9)

*# Here we customize our newly created plot, dictating the font, colormap, and title.*

.. raw:: latex

   \setlength{\parskip}{0pt}

*# Logarithmic data is used by default for this plot, so we turn it off.*

plot.set_font({‘family’:‘Bitstream Vera Sans’, ‘style’:‘italic’,

‘weight’:‘normal’, ‘size’:14, ‘color’:‘red’})

plot.set_log(field, False)

plot.set_cmap(field, ‘jet’)

plot.annotate_title(‘Off-Axis Density Projection’)

*# The actual size of the image can also be set. Note that the units are in inches.*

.. raw:: latex

   \setlength{\parskip}{0pt}

plot.set_window_size(8.0)

plot.save(‘off_axis_density’)

.. raw:: latex

   \centering

.. figure:: OffAxisProjection_density
   :alt: Off-axis density projection
   :width: 4in

   Off-axis density projection

Verification Test Problems
==========================

Hydrodynamics Test Problems
---------------------------

Sod’s Problem (and Other Shock Tube Problems)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Exec/hydro_tests/Sod problem directory sets up a general one-dimensional
shock tube. The left and right primitive-variable states are specified
and the solution evolves until a user-specified end time. For a simple
discontinuity, the exact solution can be found from an exact Riemann
solver. For this problem, the exact solutions were computed with the
exact Riemann solver from Toro :raw-latex:`\cite{toro:1997}`, Chapter 4.

Sod’s Problem
^^^^^^^^^^^^^

The problem :raw-latex:`\cite{sod:1978}` is a simple shock tube problem that
exhibits a shock, contact discontinuity, and a rarefaction wave.
The initial conditions are:

.. math::

   \begin{array}{l}
   \rho_L = 1 \\
   u_L = 0 \\
   p_L = 1
   \end{array}
   \qquad
   \begin{array}{l}
   \rho_R = 0.125 \\
   u_R = 0 \\
   p_R = 0.1
   \end{array}

The gamma_law equation of state is used with :math:`\gamma = 1.4`.
The system is evolved until :math:`t = 0.2` s. Setups for 1-, 2-, and 3-d
are provided. The following inputs files and probin files setup the
Sod’s problem:

.. raw:: latex

   \centering

+--------------+--------------+-----------------------------------------+
| inputs-sod-x | probin-sod-x | Sod’s problem along :math:`x`-direction |
+--------------+--------------+-----------------------------------------+
| inputs-sod-y | probin-sod-y | Sod’s problem along :math:`y`-direction |
+--------------+--------------+-----------------------------------------+
| inputs-sod-z | probin-sod-z | Sod’s problem along :math:`z`-direction |
+--------------+--------------+-----------------------------------------+

[Table:Sod]

For multi-dimensional runs, the directions transverse to the jump are
kept constant. We use a CFL number of 0.9, an initial timestep shrink
(castro.init_shrink) of 0.1, and the maximum factor by which
the timestep can increase (castro.change_max) of 1.05.

.. raw:: latex

   \centering

.. figure:: sod_3d
   :alt: [fig:sod] Castro solution for Sod’s problem run in 3-d,
   with the newest ppm limiters,
   along the :math:`x`, :math:`y`, and :math:`z` axes. A coarse grid of 32 zones in the
   direction of propagation, with 2 levels of refinement was used. The
   analytic solution appears as the red line.
   :width: 4.75in

   [fig:sod] Castro solution for Sod’s problem run in 3-d,
   with the newest ppm limiters,
   along the :math:`x`, :math:`y`, and :math:`z` axes. A coarse grid of 32 zones in the
   direction of propagation, with 2 levels of refinement was used. The
   analytic solution appears as the red line.

.. raw:: latex

   \centering

.. figure:: sod_3d_ppm0
   :alt: [fig:sod_ppm0] Castro solution for Sod’s problem run in 3-d,
   with the piecewise-linear Godunov method with limiters,
   along the :math:`x`, :math:`y`, and :math:`z` axes. A coarse grid of 32 zones in the
   direction of propagation, with 2 levels of refinement was used. The
   analytic solution appears as the red line.
   :width: 4.75in

   [fig:sod_ppm0] Castro solution for Sod’s problem run in 3-d,
   with the piecewise-linear Godunov method with limiters,
   along the :math:`x`, :math:`y`, and :math:`z` axes. A coarse grid of 32 zones in the
   direction of propagation, with 2 levels of refinement was used. The
   analytic solution appears as the red line.

Figure \ `[fig:sod] <#fig:sod>`__ shows the Castro solution using the newest PPM limiters
compared to the analytic
solution, showing the density, velocity, pressure, and internal energy.
Figure \ `[fig:sod_ppm0] <#fig:sod_ppm0>`__ is the same as Figure \ `[fig:sod] <#fig:sod>`__,
but with the piecewise-linear Godunov method with limiters,
shown for comparison.

The Verification subdirectory includes the analytic solution for
the Sod problem sod-exact.out, with :math:`\gamma = 1.4`. 1-d slices
can be extracted from the Castro plotfile using the fextract tool
from BoxLib/Tools/Postprocessing/F_Src/.
The steps to generate this verification plot with Castro are:

#. in Exec/hydro_tests/Sod, build the Castro executable in 3-d

#. | run the Sod problem with Castro in the :math:`x`, :math:`y`, and :math:`z` directions:
   | ./Castro3d.Linux.Intel.Intel.ex inputs-sod-x
   | ./Castro3d.Linux.Intel.Intel.ex inputs-sod-y
   | ./Castro3d.Linux.Intel.Intel.ex inputs-sod-z

#. build the fextract tool in BoxLib/Tools/Postprocessing/F_Src/.

#. | run fextract on the Castro output to generate 1-d slices
     through the output:
   | fextract3d.Linux.Intel.exe -d 1 -s sodx.out -p sod_x_plt00034
   | fextract3d.Linux.Intel.exe -d 2 -s sody.out -p sod_y_plt00034
   | fextract3d.Linux.Intel.exe -d 3 -s sodz.out -p sod_z_plt00034

#. copy the sodx/y/z.out files into the Verification directory.

#. | in Verification run the gnuplot script sod_3d.gp as:
   | gnuplot sod_3d.gp
   | This will produce the figure sod_3d.eps.

Double Rarefaction
^^^^^^^^^^^^^^^^^^

The double rarefaction is the “Test 2” problem described by Toro
:raw-latex:`\cite{toro:1997}`, Chapter 6. In this test, the center of the domain
is evacuated as two rarefaction waves propagate in each direction, outward
from the center. It is difficult to get the internal energy to
behave at the center of the domain because we are creating a vacuum.
The initial conditions are:

.. math::

   \begin{array}{l}
   \rho_L = 1 \\
   u_L = -2 \\
   p_L = 0.4
   \end{array}
   \qquad
   \begin{array}{l}
   \rho_R = 1 \\
   u_R = 2 \\
   p_R = 0.4
   \end{array}

The gamma_law equation of state is used with :math:`\gamma = 1.4`.
The system is evolved until :math:`t = 0.15` s. Setups for 1-, 2-, and 3-d
are provided. The following inputs files and probin files setup the
Sod’s problem:

.. raw:: latex

   \centering

+-----------------------+-----------------------+-----------------------+
| inputs-test2-x        | probin-test2-x        | Double rarefaction    |
|                       |                       | problem along         |
|                       |                       | :math:`x`-direction   |
+-----------------------+-----------------------+-----------------------+
| inputs-test2-y        | probin-test2-y        | Double rarefaction    |
|                       |                       | problem along         |
|                       |                       | :math:`y`-direction   |
+-----------------------+-----------------------+-----------------------+
| inputs-test2-z        | probin-test2-z        | Double rarefaction    |
|                       |                       | problem along         |
|                       |                       | :math:`z`-direction   |
+-----------------------+-----------------------+-----------------------+

[Table:Sod]

We use a CFL number of 0.8, an initial
timestep shrink (castro.init_shrink) of 0.1, and the maximum factor by which
the timestep can increase (castro.change_max) of 1.05. The PPM
solver with the new limiters are used.

.. raw:: latex

   \centering

.. figure:: test2_3d
   :alt: [fig:test2] Castro solution for the double rarefaction
   problem run in 3-d, along the :math:`x`, :math:`y`, and :math:`z` axes. A coarse grid
   of 32 zones in the direction of propagation, with 2 levels of
   refinement was used. The analytic solution appears as the red
   line.
   :width: 5in

   [fig:test2] Castro solution for the double rarefaction
   problem run in 3-d, along the :math:`x`, :math:`y`, and :math:`z` axes. A coarse grid
   of 32 zones in the direction of propagation, with 2 levels of
   refinement was used. The analytic solution appears as the red
   line.

Figure \ `[fig:test2] <#fig:test2>`__ shows the Castro output, run along all 3
coordinate axes in 3-d, compared to the analytic solution.

The comparison to the analytic solution follows the same procedure as
described for the Sod’s problem above. The gnuplot script
test2_3d.gp will generate the figure, from the 1-d slices created by
fextract named test2x.out, test2y.out, and test2z.out.

Strong Shock
^^^^^^^^^^^^

The strong shock test is the “Test 3” problem described by Toro
:raw-latex:`\cite{toro:1997}`, Chapter 6. In this test, a large pressure jump
at the initial interface creates a very strong rightward moving
shock, followed very closely by a contact discontinuity.
The initial conditions are:

.. math::

   \begin{array}{l}
   \rho_L = 1 \\
   u_L = 0 \\
   p_L = 1000
   \end{array}
   \qquad
   \begin{array}{l}
   \rho_R = 1 \\
   u_R = 0 \\
   p_R = 0.01
   \end{array}

The gamma_law equation of state is used with :math:`\gamma = 1.4`.
The system is evolved until :math:`t = 0.012` s. Setups for 1-, 2-, and 3-d
are provided. The following inputs files and probin files setup the
Sod’s problem:

.. raw:: latex

   \centering

+-----------------------+-----------------------+-----------------------+
| inputs-test3-x        | probin-test3-x        | Strong shock problem  |
|                       |                       | along                 |
|                       |                       | :math:`x`-direction   |
+-----------------------+-----------------------+-----------------------+
| inputs-test3-y        | probin-test3-y        | Strong shock problem  |
|                       |                       | along                 |
|                       |                       | :math:`y`-direction   |
+-----------------------+-----------------------+-----------------------+
| inputs-test3-z        | probin-test3-z        | Strong shock problem  |
|                       |                       | along                 |
|                       |                       | :math:`z`-direction   |
+-----------------------+-----------------------+-----------------------+

[Table:Sod]

We use a CFL number of 0.9, an initial
timestep shrink (castro.init_shrink) of 0.1, and the maximum factor by which
the timestep can increase (castro.change_max) of 1.05. The PPM
solver with the new limiters are used.

.. raw:: latex

   \centering

.. figure:: test3_3d
   :alt: [fig:test3] Castro solution for the strong shock
   problem run in 3-d, along the :math:`x`, :math:`y`, and :math:`z` axes. A coarse grid
   of 32 zones in the direction of propagation, with 2 levels of
   refinement was used. The analytic solution appears as the red
   line.
   :width: 5in

   [fig:test3] Castro solution for the strong shock
   problem run in 3-d, along the :math:`x`, :math:`y`, and :math:`z` axes. A coarse grid
   of 32 zones in the direction of propagation, with 2 levels of
   refinement was used. The analytic solution appears as the red
   line.

Figure \ `[fig:test3] <#fig:test3>`__ shows the Castro output, run along all 3
coordinate axes in 3-d, compared to the analytic solution.

The comparison to the analytic solution follows the same procedure as
described for the Sod’s problem above. The gnuplot script
test3_3d.gp will generate the figure, from the 1-d slices created by
fextract named test3x.out, test3y.out, and test3z.out.

Sedov Problem
~~~~~~~~~~~~~

The (or Sedov-Taylor) blast wave is a standard hydrodynamics
test problem. A large amount of energy is placed into a very small
volume, driving a spherical (or cylindrical in 2-d Cartesian
coordinates) blast wave. Analytic solutions were found by Sedov
:raw-latex:`\cite{sedov:1959}`.

A cylindrical blast wave (e.g. a point explosion in a 2-d plane) can
be modeled in 2-d Cartesian coordinates. A spherical blast wave can
be modeled in 1-d spherical, 2-d axisymmetric (cylindrical :math:`r`-:math:`z`), or 3-d
Cartesian coordinates. This provides a good test on the geometric
factors in the hydrodynamics solver.
We use a publically available code, sedov3.f
:raw-latex:`\cite{timmes_sedov_code}`, to generate the analytic solutions.

The Castro implementation of the Sedov problem is in Exec/hydro_tests/Sedov.
A number of different inputs/probin files are provided, corresponding
to different Sedov/Castro geometries. The main ones are:

[Table:Sod]

In the Sedov problem, the explosion energy, :math:`E_\mathrm{exp}` (in units
of energy, not energy/mass or energy/volume)
is to be deposited into a single point, in a medium of uniform ambient
density, :math:`\rho_\mathrm{ambient}`, and pressure, :math:`p_\mathrm{ambient}`.
Initializing the problem can be difficult because the small volume is
typically only a cell in extent. This can lead to grid imprinting in
the solution. A standard solution (see for example :raw-latex:`\cite{omang:2006}`
and the references therein)
is to convert the explosion energy into a pressure contained within a
certain volume, :math:`V_\mathrm{init}`, of radius :math:`r_\mathrm{init}` as

.. math:: p_\mathrm{init} = \frac{(\gamma - 1) E_\mathrm{exp}}{V_\mathrm{init}} \enskip .

This pressure is then deposited in all of the cells where :math:`r <
r_\mathrm{init}`.

To further minimize any grid effects, we do subsampling
in each zone: each zone is divided it into :math:`N_\mathrm{sub}` subzones in each
coordinate direction, each subzone is initialized independently, and
then the subzones are averaged together (using a volume weighting for
spherical or cylindrical/axisymmetric Castro grids) to determine the
initial state of the full zone.

For these runs, we use :math:`\rho_\mathrm{ambient} = 1`,
:math:`p_\mathrm{ambient} = 10^{-5}`, :math:`E_\mathrm{exp} = 1`, :math:`r_\mathrm{init}
 = 0.01`, and :math:`N_\mathrm{sub} = 10`. A base grid with 32 zones in each
coordinate direction plus 3 levels of refinement is used (the finest
mesh would coorespond to 256 zones in a coordinate direction). The
domain runs from 0 to 1 in each coordinate direction.

Analysis routines for the Sedov problem are provided in
Castro/Diagnostics/Sedov/. These routines will
average the Castro solution over angles, using the proper geometric
weighting, to produce an average profile as a function of radius.
The following routines correspond to the inputs files described above:

Spherical Blast Wave
^^^^^^^^^^^^^^^^^^^^

A spherical Sedov explosion can be modeled in 1-d spherical, 2-d
cylindrical (axisymmetric), or 3-d Cartesian coordinates, using the
inputs files described in Table \ `[table:sedov_inputs] <#table:sedov_inputs>`__. A 1-d radial
profile can be extracted using the appropriate fsedov routine,
as listed in Table \ `[table:fsedov] <#table:fsedov>`__. For example, to run and process
the 2-d cylindrical Sedov explosion, one would do:

#. in Exec/hydro_tests/Sedov, build the Castro executable in 2-d

#. | run the spherical Sedov problem with Castro in 2-d cylindrical coordinates:
   | ./Castro2d.Linux.Intel.Intel.ex inputs.2d.sph_in_cylcoords

#. build the fsedov2d_sph_in_cylcoords tool in
   Castro/Diagnostics/Sedov.

#. | run fsedov2d_sph_in_cylcoords on the Castro output to generate 1-d radial
     profiles:
   | fsedov2d_sph_in_cylcoords.Linux.Intel.exe -s sedov_2d_sph_in_cyl.out :math:`\mathtt{\backslash}` 
   | :math:`~~~~~`\ -p sedov_2d_sph_in_cyl_plt00246

A similar procedure can be used for the 1-d and 3-d spherical Sedov
explosions (with the output named sedov_1d_sph.out and
sedov_3d_sph.out respectively). Once this is done, the
sedov_sph.gp gnuplot script can be used to make a plot comparing
the 3 solutions to the analytic solution, spherical_sedov.dat.

Figure \ `[fig:sedov_sph] <#fig:sedov_sph>`__ shows the comparison of the 3 Castro spherical Sedov explosion simulations to the analytic solution.

.. raw:: latex

   \centering

.. figure:: sedov_sph
   :alt: [fig:sedov_sph] Castro solution for the Sedov blast wave problem
   run in 1-d spherical, 2-d axisymmetric, and 3-d Cartesian coordinates.
   Each of these geometries produces a spherical Sedov explosion.
   :width: 5in

   [fig:sedov_sph] Castro solution for the Sedov blast wave problem
   run in 1-d spherical, 2-d axisymmetric, and 3-d Cartesian coordinates.
   Each of these geometries produces a spherical Sedov explosion.

Cylindrical Blast Wave
^^^^^^^^^^^^^^^^^^^^^^

.. raw:: latex

   \centering

.. figure:: sedov_cyl
   :alt: [fig:sedov_cyl] Castro solution for the Sedov blast wave problem
   run in 2-d Cartesian coordinates. This corresponds to a cylindrical
   Sedov explosion.
   :width: 5in

   [fig:sedov_cyl] Castro solution for the Sedov blast wave problem
   run in 2-d Cartesian coordinates. This corresponds to a cylindrical
   Sedov explosion.

Rayleigh-Taylor
~~~~~~~~~~~~~~~

2D. Domain size 0.5 by 1.0. 256 by 512 cells, single level
calculation. Periodic in x, solid walls on top and bottom in y.
Gamma law gas with :math:`\gamma=1.4`, no reactions. Zero initial velocity.
Constant :math:`|{\bf g}|=1`. The density profile is essentially :math:`\rho=1` on
bottom, :math:`\rho=2` on top, but with a perturbation. A single-mode
perturbation is constructed as:

.. math:: \tilde y(x) = 0.5 + 0.01 \frac{\cos(4\pi x) + \cos(4\pi(L_x - x))}{2}

We note that the symmetric form of the cosine is done to ensure that
roundoff error does not introduce a left-right asymmetry in the problem.
Without this construction, the R-T instability will lose its symmetry
as it evolves. This then applied to the interface with a tanh profile
to smooth the transition between the high and low density material:

.. math:: \rho(x,y) = 1 + 0.5\left[1+\tanh\left(\frac{y-\tilde y(x)}{0.005}\right)\right]

Hydrostatic pressure with :math:`p=5.0` at bottom of domain, assuming
:math:`\rho=1` on the lower half of the domain, and :math:`\rho=2` on the upper
half and no density perturbation. We run to :math:`t=2.5` with piecewise
linear, old PPM, and new PPM. CFL=0.9. See Figure `[fig:RT] <#fig:RT>`__.

.. raw:: latex

   \centering

.. figure:: RT_ppm_type
   :alt: [fig:RT]Rayleigh-Taylor with different PPM types.
   :width: 6.5in

   [fig:RT]Rayleigh-Taylor with different PPM types.

Gravity Test Problems
---------------------

Radiation Test Problems
-----------------------

There are two photon radiation solvers in Castro—a gray solver and a
multigroup solver. The gray solver follows the algorithm outlined
in :raw-latex:`\cite{howellgreenough:2003}`. We use the notation described in that
paper. In particular, the radiation energy equation takes the form
of:

.. math::

   \frac{\partial E_R}{\partial t} = 
    \nabla \cdot \left ( \frac{c \lambda(E_R)}{\kappa_R} \nabla E_R \right ) +
    \kappa_P (4 \sigma T^4 - c E_R )

Here, :math:`E_R` is the radiation energy density, :math:`\kappa_R` is the
Roseland-mean opacity, :math:`\kappa_P` is the Planck-mean opaciy, and
:math:`\lambda` is a quantity :math:`\le 1/3` that is subjected to limiting to
keep the radiation field causal. Castro allows for :math:`\kappa_R`
and :math:`\kappa_P` to be set independently as power-laws.

Light Front
~~~~~~~~~~~

The light front problem tests the ability of the radiation solver to
operate in the free-streaming limit. A radiation front is
estabilished by initializing one end of the computational domain with
a finite radiation field, and zero radiation field everywhere else.
The speed of propagation of the radiation front is keep in check by
the flux-limiters, to prevent it from exceeding :math:`c`.

Diffusion of a Gaussian Pulse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The diffusion of a Gaussian pulse problem tests the diffusion term in
the radiation energy equation. The radiation energy density is
initialized at time :math:`t = t_0` to a Gaussian distribution:

.. math:: E_R = (E_R)_0 \exp \left \{ - \frac{1}{4 D t_0} |r - r_0|^2 \right \} \enskip .

As the radiation diffuses, the overall distribution will remain
Gaussian, with the time-dependent solution of:

.. math:: E_R = (E_R)_0 \frac{t_0}{t_0 + t} \exp \left \{ -\frac{1}{4 D (t_0 + t)} |r - r_0|^2 \right \}

Radiation Source Problem
~~~~~~~~~~~~~~~~~~~~~~~~

The radiation source problem tests the coupling between the radiation
field and the gas energy through the radiation source term. The
problem begins with the radiation field and gas temperature out of
equilibrium. If the gas is too cool, then the radiation field will
heat it. If the gas is too hot, then it will radiate and cool. In
each case, the gas energy and radiation field will evolve until
thermal equilibrium is achieved.

Our implementation of this problem follows that of
:raw-latex:`\cite{swestymyra:2009}`.

.. raw:: latex

   \centering

.. figure:: radiating_source
   :alt: [fig:radsource] Castro solution for radiating source
   test problem. Heating and cooling solutions are shown as a function
   of time, compared to the analytic solution. The gray photon solver
   was used.
   :width: 5in

   [fig:radsource] Castro solution for radiating source
   test problem. Heating and cooling solutions are shown as a function
   of time, compared to the analytic solution. The gray photon solver
   was used.

Radiating Sphere
~~~~~~~~~~~~~~~~

The radiating sphere () is a multigroup radiation
test problem. A hot sphere is centered at the origin in a spherical
geometry. The spectrum from this sphere follows a Planck
distribution. The ambient medium is at a much lower temperature. A
frequency-dependent opacity makes the domain optically thin for high
frequecies and optically thick for low frequency. At long times, the
solution will be a combination of the blackbody radiation from the
ambient medium plus the radiation that propagated from the hot sphere.
An analytic solution exists :raw-latex:`\cite{graziani:2008}` which gives the
radiation energy as a function of energy group at a specified time and
distance from the radiating sphere.

Our implementation of this problem is in Exec/radiation_tests/RadSphere and
follows that of :raw-latex:`\cite{swestymyra:2009}`. The routine that computes
the analytic solution is provided as analytic.f90.

.. raw:: latex

   \centering

.. figure:: radiating_sphere
   :alt: [fig:radsphere] Castro solution for radiating sphere problem,
   showing the radiation energy density as a function of energy group.
   This test was run with 64 photon energy groups.
   :width: 5in

   [fig:radsphere] Castro solution for radiating sphere problem,
   showing the radiation energy density as a function of energy group.
   This test was run with 64 photon energy groups.

Regression Testing
------------------

An automated regression test suite for Castro (or any BoxLib-based
code) written in Python exists in BoxLib/Tools/RegressionTesting.
Details of its use are provided in the BoxLib User’s Guide.

.. raw:: latex

   \backmatter

.. raw:: latex

   \addcontentsline{toc}{chapter}{References}

.. raw:: latex

   \bibliographystyle{plain}

.. raw:: latex

   \clearpage

.. raw:: latex

   \if@twoside

.. raw:: latex

   \ifodd

.. raw:: latex

   \c@page

.. raw:: latex

   \else

.. raw:: latex

   \hbox{}

.. raw:: latex

   \thispagestyle{empty}

.. raw:: latex

   \newpage

.. raw:: latex

   \if@twocolumn

.. raw:: latex

   \hbox{}

.. raw:: latex

   \newpage

.. raw:: latex

   \fi

.. raw:: latex

   \fi

.. raw:: latex

   \fi

.. raw:: latex

   \phantomsection

.. raw:: latex

   \addcontentsline{toc}{chapter}{Index}

.. raw:: latex

   \printindex

.. [1]
   earlier versions of Castro used the
   BoxLib library

.. [2]
   Note: previously the radiation
   solver was distributed separately as CastroRadiation.git,
   but this has been merged into the main Castro respository

.. [3]
   Each of these will recognize it as the
   BoxLib format.

.. [4]
   Note: some older code will use a special AMReX preprocessor macro,
   , defined in , that converts
   the C multifab into a Fortran array and its lo and hi indices.
   Additionally, some older code will wrap the Fortran subroutine name
   in an additional preprocessor macro,
   to handle the name mangling between Fortran and C. This later
   macro is generally not needed any more because of Fortran 2003
   interoperability with C (through the Fortran bind keyword).

.. [5]
   the way to read these complicated
   C declarations is right-to-left. So ‘const int\* lo‘ means
   ‘lo‘ is a integer pointer to a memory space that is constant. See
   https://isocpp.org/wiki/faq/const-correctness#ptr-to-const

.. [6]
   for clarity and continuity in this
   documentation, some of the variable names have been changed
   compared to the actual code

.. [7]
   the integer values are defined in

.. [8]
   available separately at
   https://github.com/BoxLib-Codes/wdmerger

.. [9]
   Note: this functionality assumes that only the
   coarse grid touches the physical boundary. It does not use
   any use masks to prevent double counting if multiple levels
   touch the boundary.

.. [10]
   The correction for gravity is slightly different since we directly compute the time-centered gravitational source term using the hydrodynamic fluxes.

.. [11]
   Note: The PrescribedGrav
   option and text here were contributed by Jan Frederik Engels of
   University of Gottingen.

.. [12]
   Note: The GR
   code and text here were contributed by Ken Chen of Univ. of
   Minnesota.

.. [13]
   at the moment, we
   don’t have a way to allow for the EOS to provide radiation pressure
   if the Castro radiation is used solely for neutrinos, but this is
   something that could be added easily.

.. [14]
   One can simplify interpolation with the cell-centered velocity. However, this can lead to decoupling of the pressure and the velocity components, possibly resulting in instability. This can be avoided with the face-centered velocity

.. [15]
   In the absence of a global field like
   the gravitational potential, this would only need to be done on the
   coarse level, as we always assume that the solution on the fine grid is
   correct and average it down to the coarse grid. In Castro we do it by
   default on the fine level too in anticipation of the fact that gravity
   is a common component of many of our production science
   simulations. This could be generalized so that if you aren’t using any
   global force fields, you don’t bother updating the fine level. If this
   is important to the science you want to do, please let the Castro developers know and we can look into it.

.. [16]
   in general it may be desirable for this to be a
   source-term specific setting, so that some source terms that are cheap
   or physically important are re-computed after a synchronization can be
   set to update, while others can be disabled. If this is important for
   your science application, please let the developers know, as this would
   be a straightforward extension of the current architecture.

.. [17]
   If this scheme
   is generalized to higher-order methods, in principle all one would need
   to do is integrate the fluxes until :math:`\Delta t / 2`, which is what we are
   doing here for the constant-in-time flux case.

.. |A model atmosphere (*left* panel) and the trajectories of 500 particles (*right* panel) following the fluid motion on the atmosphere. The particles are initially positioned at five different heights, :math:`y=13000\mathrm{~km},~11000\mathrm{~km},~ 8000\mathrm{~km},~ 6000\mathrm{~km}, ~38000\mathrm{~km}` (100 particles at each height). In the *left* panel, the arrows roughly show the fluid motion. In the *right* panel, the solid lines represent the trajectories of the particles. | image:: fluid_motion
   :width: 2.5in
.. |A model atmosphere (*left* panel) and the trajectories of 500 particles (*right* panel) following the fluid motion on the atmosphere. The particles are initially positioned at five different heights, :math:`y=13000\mathrm{~km},~11000\mathrm{~km},~ 8000\mathrm{~km},~ 6000\mathrm{~km}, ~38000\mathrm{~km}` (100 particles at each height). In the *left* panel, the arrows roughly show the fluid motion. In the *right* panel, the solid lines represent the trajectories of the particles. | image:: tracer_trajectory
   :width: 2.39in
.. |Data from checkpoint file before and after the domain has been coarsened and grown. This case
uses **star_at_center = 0** and **ref_ratio**\ =2. The first grown example has
**grown_factor**\ =2, the second has **grown_factor**\ =3. In all figures the level 0 grids
are shown in white, the level 1 grids in red, the level 2 grids in yellow, and in the grown figures,
the level 3 grids are in pink.| image:: orig_corner
   :width: 3in
.. |Data from checkpoint file before and after the domain has been coarsened and grown. This case
uses **star_at_center = 0** and **ref_ratio**\ =2. The first grown example has
**grown_factor**\ =2, the second has **grown_factor**\ =3. In all figures the level 0 grids
are shown in white, the level 1 grids in red, the level 2 grids in yellow, and in the grown figures,
the level 3 grids are in pink.| image:: grown_corner_2
   :width: 3in
.. |Data from checkpoint file before and after the domain has been coarsened and grown. This case
uses **star_at_center = 0** and **ref_ratio**\ =2. The first grown example has
**grown_factor**\ =2, the second has **grown_factor**\ =3. In all figures the level 0 grids
are shown in white, the level 1 grids in red, the level 2 grids in yellow, and in the grown figures,
the level 3 grids are in pink.| image:: grown_corner_3
   :width: 3in
.. |Data from checkpoint file before and after the domain has been coarsened and grown. This case
uses **star_at_center = 0** and **ref_ratio**\ =2. The first grown example has
**grown_factor**\ =2, the second has **grown_factor**\ =3. In all figures the level 0 grids
are shown in white, the level 1 grids in red, the level 2 grids in yellow, and in the grown figure,
the level 3 grids are in pink. | image:: orig_center
   :width: 3in
.. |Data from checkpoint file before and after the domain has been coarsened and grown. This case
uses **star_at_center = 0** and **ref_ratio**\ =2. The first grown example has
**grown_factor**\ =2, the second has **grown_factor**\ =3. In all figures the level 0 grids
are shown in white, the level 1 grids in red, the level 2 grids in yellow, and in the grown figure,
the level 3 grids are in pink. | image:: grown_center_2
   :width: 3in
.. |Data from checkpoint file before and after the domain has been coarsened and grown. This case
uses **star_at_center = 0** and **ref_ratio**\ =2. The first grown example has
**grown_factor**\ =2, the second has **grown_factor**\ =3. In all figures the level 0 grids
are shown in white, the level 1 grids in red, the level 2 grids in yellow, and in the grown figure,
the level 3 grids are in pink. | image:: grown_center_3
   :width: 3in
