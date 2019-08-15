**********************
Introduction to Castro
**********************

Castro is a adaptive mesh, radiation hydrodynamics code that is
designed to model astrophysical reacting flows on massively parallel
computers.

Castro's major capabilities:

  * 1-, 2-, and 3-dimensional unsplit, 2nd-order hydrodynamics

  * multigroup flux-limited diffusion radiation hydrodynamics

  * adaptive mesh refinement with subcycling; jumps of 2x and 4x between levels

  * arbitrary equation of state (gamma-law and stellar EOSes are bundled)

  * general nuclear reaction networks

  * explicit thermal diffusion

  * full Poisson gravity (with isolated boundary conditions)

  * rotation (in the co-rotating frame) in 2-d axisymmetric and 3-d.

  * parallelization via MPI + OpenMP or MPI + CUDA

Units and Conventions
=====================

Castro works in CGS units unless otherwise specified.
Table \ `[table:units] <#table:units>`__ shows some of the common symbols / names used
throughout the code documentation and papers.

.. table:: [table:units] Common quantities and units.

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
in ``Castro/constants/constants_cgs.f90``

