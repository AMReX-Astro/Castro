Introduction to Castro
======================

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
=====================

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
   | :math:`\ub`           | :math:`\mathrm{cm~s^{ | velocity vector       |
   |                       | -1}}`                 |                       |
   +-----------------------+-----------------------+-----------------------+
   | :math:`p`             | :math:`\mathrm{dyn~cm | pressure              |
   |                       | ^{-2}}`               |                       |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\gb`           | :math:`\mathrm{cm~s^{ | gravitational         |
   |                       | -2}}`                 | acceleration          |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\Sb`           | varies                | source term           |
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
   | :math:`\omegadot_k`   | :math:`\mathrm{s^{-1} | species creation rate |
   |                       | }`                    | (from reactions)      |
   +-----------------------+-----------------------+-----------------------+

Physical constants, again using the CGS system are available
in Castro/constants/constants_cgs.f90
