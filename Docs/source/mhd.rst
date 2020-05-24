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

   \begin{align}
   \frac{\partial \rho X_k}{\partial t} + \frac{\partial}{\partial x_j} ( \rho U_j X_k) &= \rho \omegadot_k \\
   \frac{\partial \rho U_j}{\partial t} + \frac{\partial}{\partial x_j} (\rho U_i U_j + p \delta_{ij} - B_j B_i) &= S_{\rho U_j} \\
   \frac{\partial \rho E}{\partial t} + \frac{\partial}{\partial x_j} \left [ U_j (\rho E + p) - B_j B_i u_i \right ] &= S_{\rho E} + \rho H_\mathrm{nuc} \\
   \frac{\partial B_i}{\partial t} &= -\frac{\partial}{\partial x_j} (U_j B_i - B_j U_i)
   \end{align}

where

.. math::

   p = p_g + \frac{1}{2} B^2

and the :math:`S` sources represent the hydrodynamical sources and
the remainder are reaction sources.  The MHD solver uses the same
time-advancement driver as the hydrodynamic CTU driver
(:ref:`sec:strangctu`), using Strang splitting for the reactions (by
default).

The constrained transport algorithm algebraically ensures that
:math:`\nabla \cdot {\bf B} = 0`.  Throughout the algorithm, the
magnetic fields are face-centered, the electric fields are
edge-centered, and the conserved state is cell-centered (or averages),
unless specified otherwise.  We use the following indexing
conventions:

  * ``U(i,j,k)`` is the cell-centered state, :math:`U_{i,j,k}`

  * ``qleft(i,j,k)`` is the interface value of a state, for example,
    in the x-direction this would be :math:`q_{i-1/2,j,k}`

  * ``Bx(i,j,k)`` is the face-centered x-component of the magnetic field,
    :math:`B_{x,i-1/2,j,k}`

  * ``Ex(i,j,k)`` is the edge-centered x-component of the electric field,
    :math:`E_{x,i,j-1/2,k-1/2}`

Hydrodynamics Update
====================

We use piecewise linear reconstruction with a characteristic projection
and the full 12-Riemann solve corner transport upwind method.  The HLLD
Riemann solver is used.


Electric Update
===============
