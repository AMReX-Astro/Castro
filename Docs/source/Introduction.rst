**********************
Introduction to Castro
**********************

Castro is a adaptive mesh, radiation/MHD hydrodynamics code that is
designed to model astrophysical reacting flows on massively parallel
computers.

Castro's major capabilities:

  * 1-, 2-, and 3-dimensional unsplit, 2nd-order finite-volume
    hydrodynamics; 4th order hydro for uniform grids.
    (see :ref:`ch:hydro`)

  * 3-dimension constrained transport ideal MHD (single level only currently)
    (see :ref:`ch:mhd`)

  * multigroup flux-limited diffusion radiation hydrodynamics
    (see :ref:`ch:radiation`)

  * generalized retry mechanism for recovering from physical
    violations over a timestep (see :ref:`ch:retry`)

  * adaptive mesh refinement with subcycling; jumps of 2x and 4x
    between levels (see :ref:`ch:amr`)

  * arbitrary equation of state (provided by the companion StarKiller
    Microphysics project)

  * general nuclear reaction networks

  * explicit thermal diffusion (see :ref:`ch:diffusion`)

  * full Poisson gravity (with isolated boundary conditions)
    and a conservative energy formulation (see :ref:`ch:gravity`)

  * rotation (in the co-rotating frame) in 2-d axisymmetric and 3-d
    (see :ref:`ch:rotation`)

  * spectral deferred corrections time integration for coupling hydro
    and reactions (see :ref:`ch:sdc`)

  * parallelization via MPI + OpenMP (CPUs), MPI + CUDA (NVIDIA GPUs), or MPI + HIP (AMD GPUs)


Development Model
=================

Castro is developed on github (https://github.com/amrex-astro/Castro
). The ``main`` branch is stable and can be used for day-to-day
science.  New changes are made via pull requests to the
``development`` branch.  This is where the ongoing regression testing
is done (both on CPU and GPU).

At the start of each month, we merge ``development`` → ``main`` and
apply a tag of the form ``YY.MM`` (e.g. ``20.02`` for Feb. 2020).  We
also create a github release and mint a Zenodo DOI using the
information in the ``.zenodo.json`` file at the root level.

Castro "core developers" are those who have made substantial code
contributions (details are in the main ``README.md``).  These
developers are coauthors on the Zenodo DOI and of any papers
describing Castro generally (science papers coauthors are decided by
the science paper lead).


Units and Conventions
=====================

Castro works in CGS units unless otherwise specified.
:numref:`table:units` shows some of the common symbols / names used
throughout the code documentation and papers.

.. _table:units:

.. table:: Common quantities and units.

   +-----------------------+-----------------------+-----------------------+
   | name                  | units                 | description           |
   +=======================+=======================+=======================+
   | :math:`t`             | s                     | time                  |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\rho`          | :math:`\gcc`          | mass density          |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\ub`           | :math:`\cms`          | velocity vector       |
   +-----------------------+-----------------------+-----------------------+
   | :math:`p`             | :math:`\presunit`     | pressure              |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\gb`           | :math:`\accelunit`    | gravitational         |
   |                       |                       | acceleration          |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\Sb`           | varies                | source term           |
   +-----------------------+-----------------------+-----------------------+
   | :math:`E`             | :math:`\ergg`         | specific total energy |
   +-----------------------+-----------------------+-----------------------+
   | :math:`e`             | :math:`\ergg`         | specific internal     |
   |                       |                       | energy                |
   +-----------------------+-----------------------+-----------------------+
   | :math:`T`             | :math:`K`             | temperature           |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\kth`          | :math:`\mathrm{erg~cm | thermal conductivity  |
   |                       | ^{-1}~s^{-1}~K~{-1}}` |                       |
   +-----------------------+-----------------------+-----------------------+
   | :math:`X_k`           | –                     | mass fraction of      |
   |                       |                       | species :math:`k`     |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\omegadot_k`   | :math:`\mathrm{s^{-1} | species creation rate |
   |                       | }`                    | (from reactions)      |
   +-----------------------+-----------------------+-----------------------+

Physical constants, again using the CGS system are available
in ``Microphysics/constants/``.


