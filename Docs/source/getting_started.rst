***************
Getting Started
***************

.. note::

   Castro has two source dependencies: `AMReX <https://github.com/AMReX-Codes/amrex>`_, the adaptive mesh
   library, and `Microphysics <https://github.com/AMReX-Astro/Microphysics>`_, the collection of equations
   of state, reaction networks, and other microphysics.  The
   instructions below describe how to get these dependencies automatically
   with Castro.


The compilation process is managed by AMReX and its build system.  The
general requirements to build Castro are:

 * A C++17 (or later) compiler (for GCC, we need >= 9.0 for CUDA compilation)

 * python (>= 3.10)

 * GNU make (>= 3.82)

GCC is the main compiler suite used by the developers.

For running in parallel, an MPI library is required.  For running on GPUs:

* CUDA 11 or later is required for NVIDIA GPUs

* ROCM 6.3.1 or later is required for AMD GPUs (earlier versions have a register allocation bug)

More information on parallel builds is given in section
:ref:`ch:mpiplusx`.

Downloading the Code
====================

Castro is maintained as a repository on GitHub, and can be obtained
via standard git clone commands. First, make sure that git
is installed on your machine—we recommend version 1.7.x or higher.


#. Clone/fork the Castro repository from the AMReX-Astro GitHub
   organization, using either HTTP access:

   .. prompt:: bash

      git clone --recursive https://github.com/AMReX-Astro/Castro.git

   or SSH access if you have an SSH key enabled with GitHub:

   .. prompt:: bash

      git clone --recursive git@github.com:AMReX-Astro/Castro.git

   The ``--recursive`` option to ``git clone`` is used to ensure
   that all of Castro's dependencies are downloaded. Currently this
   requirement is for the AMReX mesh refinement framework, which is
   maintained in the AMReX-Codes organization on GitHub, and the
   Microphysics repository from the AMReX-Astro organization.
   AMReX adds the necessary code for the driver code for the simulation,
   while Microphysics adds the equations of state, reaction
   networks, and other microphysics needed to run Castro.

   If you forget to do a recursive clone, you can rectify the
   situation by running the following from the top-level of the Castro
   directory:

   .. prompt:: bash

      git submodule update --init --recursive

   .. note::

      By default, you will be on the ``main`` branch of the source.
      Development on Castro (and its primary dependencies, AMReX and
      Microphysics) is done in the ``development`` branch, so you
      should work there if you want the latest source:

      .. prompt:: bash

         git checkout development

      The Castro team runs nightly regression testing on the
      ``development`` branch, so bugs are usually found and fixed
      relatively quickly, but it is generally less stable than staying
      on the ``main`` branch.

#. We recommend setting the ``CASTRO_HOME`` environment
   variable to point to the path name where you have put Castro.
   Add the following to your ``.bashrc``:

   .. code:: bash

      export CASTRO_HOME="/path/to/Castro/"

   (or use the analogous form for a different shell).

#. You can keep the code up to date with:

   .. prompt:: bash

      git pull --recurse-submodules

   The recommended frequency for doing this is monthly, if you are on the
   stable ``main`` branch of the code; we issue a new release of the code
   at the beginning of each month.

#. *optional, for developers*: If you prefer, you can maintain AMReX and
   Microphysics as standalone repositories rather than as git submodules.
   To do so, you can clone them from GitHub using:

   .. prompt:: bash

      git clone https://github.com/AMReX-Codes/amrex.git
      git clone https://github.com/AMReX-Astro/Microphysics.git

   or via SSH as:

   .. prompt:: bash

      git clone git@github.com:/AMReX-Codes/amrex.git
      git clone git@github.com:/AMReX-Astro/Microphysics.git

   Then, set the ``AMREX_HOME`` environment variable to point to the
   ``amrex/`` directory, and the ``MICROPHYSICS_HOME`` environment
   variable to point to the ``Microphysics/`` directory. Castro will
   look there instead of in its local ``external/`` subdirectory.

Building the Code
=================

In Castro each different problem setup is stored in its own
sub-directory under ``Castro/Exec/``. You build the
Castro executable in the problem sub-directory. Here we’ll
build the Sedov problem:

#. From the directory in which you checked out the Castro git repo,
   type:

   .. prompt:: bash

      cd Castro/Exec/hydro_tests/Sedov

   This will put you into a directory in which you can run the Sedov
   problem in 1-d, 2-d or 3-d.

#. In ``Sedov/``, edit the ``GNUmakefile``, and set

   * ``DIM = 2``

     This is the dimensionality—here we pick 2-d.

   * ``COMP = gnu``

     This is the set of compilers. GNU are a good default choice (this
     will use g++). You can also choose ``intel`` for example.

     If you want to try other compilers than the GNU suite and they
     don’t work, please let us know.

   * ``DEBUG = FALSE``

     This disables debugging checks and results in a more optimized
     executable.

   * ``USE_MPI = FALSE``

     This turns off parallelization via MPI. Set it to ``TRUE`` to build
     with MPI—this requires that you have the MPI library installed on
     your machine. In this case, the build system will need to know
     about your MPI installation. This can be done by editing the
     makefiles in the AMReX tree, but the default fallback is to look
     for the standard MPI wrappers (e.g. ``mpic++`` and ``mpif90``) to do
     the build.

#. Now type ``make``.

   The resulting executable will look something like
   ``Castro2d.gnu.ex``, which means this is a 2-d version
   of the code compiled with ``COMP = gnu``.

More information on the various build options is given in :ref:`ch:buildsystem`.

Running the Code
================

#. Castro takes an input file that overrides the runtime parameter defaults.
   The code is run as:

   .. prompt:: bash

      ./Castro2d.gnu.ex inputs.2d.cyl_in_cartcoords

   This will run the 2-d cylindrical Sedov problem in Cartesian
   (:math:`x`-:math:`y` coordinates). You can see other possible
   options, which should be clear by the names of the inputs files.

#. You will notice that running the code generates directories that
   look like ``plt00000/``, ``plt00020/``, etc, and ``chk00000/``,
   ``chk00020/``, etc. These are “plotfiles” and “checkpoint”
   files. The plotfiles are used for visualization, the checkpoint
   files are used for restarting the code.

Visualization of the Results
============================

There are several options for visualizing the data. The popular
packages yt and VisIt both support the AMReX file format
natively [1]_. The standard tool used within the AMReX-community is
Amrvis, which we demonstrate here. Amrvis is available on github.


.. _sec:gettingstartedyt:

yt
^^

yt is the primary visualization and analysis tool used by the
developers.  Install yt following their instructions: `Getting yt
<https://yt-project.org/#getyt>`_ .

You should be able to read in your plotfiles using ``yt.load()`` and
do any of the plots described in the `yt Cookbook
<https://yt-project.org/doc/cookbook/index.html>`_ .

Here we do a sample visualization and analysis of the
plotfiles generated.  This section was generated from a
Jupyter notebook which can be found in
``Docs/source/yt_example.ipynb`` in the Castro repo.

.. include:: yt_example.rst


Amrvis
^^^^^^

Amrvis is a tool developed at LBNL to visualize AMReX data.  It
provides a simple GUI that allows you to quickly visualize slices and
the grid structure.

#. Get Amrvis:

   .. prompt:: bash

      git clone https://github.com/AMReX-Codes/Amrvis

   Then cd into ``Amrvis/``, edit the ``GNUmakefile`` there
   to set ``DIM = 2``, and again set ``COMP`` to compilers that
   you have. Leave ``DEBUG = FALSE``.

   Type ``make`` to build, resulting in an executable that
   looks like ``amrvis2d...ex``.

   If you want to build amrvis with ``DIM = 3``, you must first
   download and build volpack:

   .. prompt:: bash

      git clone https://ccse.lbl.gov/pub/Downloads/volpack.git

   Then cd into ``volpack/`` and type ``make``.

   Note: Amrvis requires the OSF/Motif libraries and headers. If you
   don’t have these you will need to install the development version
   of motif through your package manager.  On most Linux
   distributions, the motif library is provided by the openmotif
   package, and its header files (like ``Xm.h``) are provided by
   openmotif-devel. If those packages are not installed, then use the
   package management tool to install them, which varies from
   distribution to distribution, but is straightforward.  ``lesstif``
   gives some functionality and will allow you to build the Amrvis
   executable, but Amrvis may not run properly.

   You may then want to create an alias to amrvis2d, for example:

   .. code:: bash

      alias amrvis2d=/tmp/Amrvis/amrvis2d...ex

   where ``/tmp/Amrvis/amrvis2d...ex`` is the full path and name of
   the Amrvis executable.

#. Configure Amrvis:

   Copy the ``amrvis.defaults`` file to your home directory (you can
   rename it to ``.amrvis.defaults`` if you wish). Then edit the
   file, and change the palette line to point to the full
   path/filename of the ``Palette`` file that comes with Amrvis.

#. Visualize:

   Return to the ``Castro/Exec/hydro_tests/Sedov`` directory. You should
   have a number of output files, including some in the form ``pltXXXXX``,
   where XXXXX is a number corresponding to the timestep the file
   was output.

   ``amrvis2d filename`` to see a single plotfile, or ``amrvis2d -a
   plt*``, which will animate the sequence of plotfiles.

   Try playing around with this—you can change which variable you are
   looking at, select a region and click “Dataset” (under View) in
   order to look at the actual numbers, etc. You can also export the
   pictures in several different formats under "File/Export".

   Some users have found that Amrvis does not work properly under X
   with the proprietary Nvidia graphics driver. A fix for this is
   provided in the FAQ (§ :ref:`ch:faq:vis`)—this is due
   to the default behavior of the DAC in mappuing colors.


.. [1]
   Each of these will recognize it as the
   BoxLib format.
