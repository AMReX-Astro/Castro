.. _ch:diffusion:

*****************
Thermal Diffusion
*****************

Castro incorporates explicit thermal diffusion into the energy equations.
In terms of the specific internal energy, :math:`e`, this appears as:

.. math:: \rho \frac{De}{Dt} + p \nabla \cdot \ub = \nabla \cdot \kth \nabla T

where :math:`\kth` is the thermal conductivity, with units
:math:`\mathrm{erg~cm^{-1}~s^{-1}~K^{-1}}`.

.. note::

   To enable diffusion, you need to compile with:

   ::

     USE_DIFFUSION=TRUE

Thermal Diffusion related source codes are contained in the ``diffusion`` directory.
Thermal Diffusion is treated explicitly, by constructing the contribution to the
evolution as a source term. This is time-centered to achieve second-order accuracy
in time.

Overall Procedure
=================

Computing Thermal Conductivity
------------------------------
The main function that computes the diffusion term is ``getTempDiffusionTerm()``.
Within ``getTempDiffusionTerm()``, it first calculates the cell centered
thermal conductivity, :math:`\kth` contained in the variable ``coeff_cc``
using the function ``fill_temp_cond()`` located in ``diffusion_util.cpp``.
``fill_temp_cond()`` fills an ``eos_state`` using the
input conserved variables, which is used to calculate :math:`\kth` via
``conductivity(eos_state)``. ``conductivity()`` routine is supplied via
the ``Microphysics`` package. See :ref:`sec:conductivities` to see the
specific choices of conductivity routines available.

.. note::
   The diffusion approximation breaks down at the surface of stars,
   where the density rapidly drops and the mean free path becomes
   large. In those instances, you should use the flux limited diffusion
   module in Castro to evolve a radiation field.

Now :math:`\kth` is reset to 0 unless
:math:`\rho \gt \mathrm{castro::diffuse\_cutoff\_density}`.
And if :math:`\rho \lt \mathrm{castro::diffuse\_cutoff\_density\_hi}`,
a linear scaling of :math:`\kth` is done as:

.. math::

   \kth = \kth \cdot \frac{\rho - \mathtt{castro.diffuse\_cutoff\_density}}{\mathtt{castro.diffuse\_cutoff\_density\_hi} - \mathtt{castro.diffuse\_cutoff\_density}}

Lastly, :math:`\kth` is scaled with ``castro::diffuse_cond_scale_fac``,
a runtime parameter controlled by the user.

After obtaining cell-centered :math:`\kth`, we do an average along
i, j, and k depending on the direction to obtain face-centered MultiFabs.
This is stored in ``coeffs``, a vector of MultiFabs, and the number of
MultiFabs corresponds to geometry dimension, since a :math:`\nabla` operator
will be applied to it later.
These Multifabs have 1 ghost cells due to the nature of MLMG solvers.

.. _sec:thermal_diffusion:

Computing Thermal Diffusion
---------------------------
We are now ready to compute :math:`\nabla \cdot \kth \nabla T`
after obtaining :math:`\kth`. This is done in the ``applyop_mlmg()`` function
in ``Diffusion.cpp``. It defines ``mlabec`` an instance of class
``MLABecLaplacian`` which defines the Laplacian of the form:

.. math::
   (A\alpha - B\nabla \cdot \beta \nabla) \phi = f

where A and B are constant scalars, and :math:`\alpha` and :math:`\beta`
are scalar fields. In order to make it correctly represents our diffusion term,
we set A = 0 and B = -1, which is done via ``mlabec.setScalars(0.0, -1.0)``.
Now we recognize :math:`\beta = \kth`, which needs to be an array of MultiFab,
corresponding to dimension. This is done via ``mlabec.setBCoeffs()``.

One of the important flag that we need to pass in is to set ``setMetricTerm(true)``.
This enables modifications due to curvilinear coordinates.

Finally we create an instance of ``MLMG`` using ``mlabec``, and call
``mlmg.apply()``, which simply evaluates the LHS but do not solve it.
See more information in the amrex documentation:
https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html


Timestep Limiter
================

Castro integrates diffusion explicitly in time; this means that
there is a diffusion timestep limiter.

To see the similarity to the thermal diffusion equation, consider the
special case of constant conductivity, :math:`\kth`, and density, and
assume an ideal gas, so :math:`e = c_v T`, where :math:`c_v` is the
specific heat at constant volume.  Finally, ignore hydrodynamics, so
:math:`\ub = 0`. This gives:

.. math:: \frac{\partial T}{\partial t} = D \nabla^2 T

where :math:`D \equiv \kth/(\rho c_v)`.

The timestep limiter for this is:

.. math:: \Delta t_\mathrm{diff} \le \frac{1}{2} \frac{\Delta x^2}{D}

This is implemented in ``estdt_temp_diffusion``.


Runtime Parameters
==================

The following parameter affects diffusion:

*  ``castro.diffuse_temp``: enable thermal diffusion (0 or 1; default 0)

   A pure diffusion problem (with no hydrodynamics) can be run by setting::

      castro.diffuse_temp = 1
      castro.do_hydro = 0
* ``castro.diffuse_cond_scale_fac``: a linear scaling to :math:`\kth`. (default 0).

* ``castro.diffuse_cutoff_density``: density under which :math:`\kth` is set to 0.
  (Default: -1e200)

* ``castro.diffuse_cutoff_density_hi``: density under which a linear scaling is
  applied to :math:`\kth`, see section :ref:`sec:thermal_diffusion` for details.
  (Default: -1e200)

.. _sec:conductivities:

Conductivities
==============

To complete the setup, a thermal conductivity must be specified. These
are supplied by Microphysics, and use an interface similar to the
equation of state interface.

.. index:: CONDUCTIVITY_DIR

.. note::

   The choice of conductivity must be specified at compile-time via
   the ``CONDUCTIVITY_DIR`` option.

The current choices of conductivity are:

* ``constant`` : A simple constant thermal conductivity. This can be
  selected by setting::

       CONDUCTIVITY_DIR := constant

  in your ``GNUmakefile``. To set the value of the conductivity (e.g., to
  :math:`100`), you add to your input file::

       conductivity.const_conductivity = 100.0

* ``constant_opacity`` : A simple constant opacity. This is
  converted to an opacity as:

  .. math::

      \kth = \frac{16 \sigma_B T^3}{3 \kappa_\mathrm{const} \rho}

  where :math:`\kappa_\mathrm{const}` is the opacity, with units :math:`\mathrm{cm^2~g^{-1}}`.
  This is selected by setting::

       CONDUCTIVITY_DIR := constant_opacity

  in your ``GNUmakefile``. To set the value of the opacity, e.g., to
  0.2 (for electron scattering), set::

       conductivity.const_opacity = 0.2

  in the inputs file.

* ``stellar`` : This is the set of conductivities and radiative opacities
  appropriate for stellar interiors described in :cite:`timmes_he_flames`.


Unit Tests
==========

A simple test problem that sets up a Gaussian temperature profile
and does pure diffusion is provided as ``diffusion_test``.
