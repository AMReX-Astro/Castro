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

   The MHD solver is still under development and should not be used
   for science simulations.  See the issues in the github repo.

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

.. note::

   Note: we are following the convention in the MHD community of setting the permeability to unity.  This is
   why the magnetic pressure has the form :math:`\frac{1}{2} B^2` instead of :math:`\frac{1}{8\pi} B^2`.  In
   effect, we are carrying our magnetic field as :math:`{\bf B}^\prime = {\bf B}/\sqrt{4\pi}`.

   If you wish to express the magnetic field in Gauss, then you will need to multiply the value Castro carries
   by :math:`\sqrt{4\pi}`.


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


Problem Initialization
======================

.. index:: ca_initmag

There is an additional initialization routine for MHD,
``problem_initialize_mhd_data()``,
that is used to initialize the face-centered magnetic field
components.  This is done separately from the main conserved fluid
state.

The conserved fluid state is initialized in ``problem_initialize_state_data()`` just as
with pure hydrodynamics problems. Note that you do not need to include
the magnetic energy contribution to the total energy density, ``UEDEN``.
After this initialization, the driver handles the addition of the magnetic
contribution.


Hydrodynamics Update
====================

We use piecewise linear or piecewise parabolic reconstruction with a
characteristic projection and the full 12-Riemann solve corner
transport upwind method.  The HLLD Riemann solver is used.

Within the solver, we use the same indexing into the primitive state
as defined in :ref:`table:primlist`, with the additions of ``QMAGX``,
``QMAGY``, and ``QMAGZ`` for the cell-centered magnetic field
components and ``QPTOT`` for the total pressure (gas + magnetic).

Just like with pure hydrodynamics, the reconstruction type is
controlled by ``castro.ppm_type``, with ``0`` selecting piecewise
linear and ``1`` selecting piecewise parabolic.

For the piecewise linear method, the slope limiting is controlled by
the same ``plm_iorder`` runtime parameter as for hydrodynamics.
Additionally, we have an option ``mhd_limit_characteristic`` that
controls whether you want to do the slope limiting on the
characteristic variables (the default) or the primitive variables.

Electric Update
===============

Coupled to the hydrodynamics update is the electric field update.
Here we update the components of E using the contact upwind scheme
first proposed in :cite:`GS2005`.  The updated electric field then
gives the magnetic field via Faraday's law and the discretization ensures
that :math:`\nabla \cdot {\bf B} = 0`.
