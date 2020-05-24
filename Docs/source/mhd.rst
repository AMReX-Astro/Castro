.. _ch:mhd:

***
MHD
***

Introduction
============

Castro implements a constrained transport (CT) corner transport upwind
(CTU) ideal MHD scheme based on the work of Miniati & Martin
:cite:`miniati_martin`.  MHD is enabled by compiling with ``USE_MHD =
TRUE``.  This replaces the pure hydrodynamics solver, but uses the
same driver as the CTU hydrodynamics solver.  This means that all of
the source terms supported by the hydrodynamics solver are also
supported by MHD.

.. note::

   The MHD solver supports 3-d only.

   Currently the MHD solver is single-level only.  AMR support is forthcoming.

Equations and Data Structures
=============================

The ideal MHD equations we solve appear as:

.. math::

   \frac{\partial \rho X_k}{\partial t} + \frac{\partial}{\partial x_j} ( \rho U_j X_k) = \rho \omegadot_k


