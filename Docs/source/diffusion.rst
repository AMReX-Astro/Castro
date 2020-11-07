.. _ch:diffusion:

*********
Diffusion
*********


Thermal Diffusion
=================

Castro incorporates explicit thermal diffusion into the energy equation.
In terms of the specific internal energy, :math:`e`, this appears as:

.. math:: \rho \frac{De}{Dt} + p \nabla \cdot \ub = \nabla \cdot \kth \nabla T

where :math:`\kth` is the thermal conductivity, with units
:math:`\mathrm{erg~cm^{-1}~s^{-1}~K^{-1}}`.

To see the similarity to the thermal diffusion equation, consider the
special case of constant conductivity, :math:`\kth`, and density, and
assume an ideal gas, so :math:`e = c_v T`, where :math:`c_v` is the
specific heat at constant volume.  Finally, ignore hydrodynamics, so
:math:`\ub = 0`. This gives:

.. math:: \frac{\partial T}{\partial t} = D \nabla^2 T

where :math:`D \equiv \kth/(\rho c_v)`. Solving this equation
explicitly requires a timestep limiter of

.. math:: \Delta t_\mathrm{diff} \le \frac{1}{2} \frac{\Delta x^2}{D}

(this is implemented in ``ca_estdt_temp_diffusion`` in
``Castro/Source/driver/timestep.F90``).

Support for diffusion must be compiled into the code by setting
``USE_DIFFUSION = TRUE`` in your ``GNUmakefile``. It is treated
explicitly, by constructing the contribution to the evolution as a
source term. This is time-centered to achieve second-order accuracy
in time.

The following parameter affects diffusion:

-  ``castro.diffuse_temp``: enable thermal diffusion (0 or 1; default 0)

A pure diffusion problem (with no hydrodynamics) can be run by setting::

    castro.diffuse_temp = 1
    castro.do_hydro = 0

To complete the setup, a thermal conductivity must be specified. The
interface for the conductivity is::

      subroutine actual_conductivity(eos_state)

        type (eos_t), intent(inout) :: eos_state

The density, temperature, and mass fractions come in through the
``eos_state`` type. An EOS call is done in Castro just before the call to
``thermal_conductivity``, so you can assume that the entire state is
consistent.  The conductivity is filled in ``eos_state % conductivity``.

.. index:: CONDUCTIVITY_DIR

There are two conductivity routines provided with Castro by default:

-  ``constant`` : A simple constant thermal conductivity. This can be
   selected by setting::

       CONDUCTIVITY_DIR := constant

   in your ``GNUmakefile``. To set the value of the conductivity (e.g., to
   :math:`100`), you add to your probin file’s ``&extern`` namelist::

       const_conductivity = 100.0

-  ``constant_opacity`` : A simple constant opacity. This is
   converted to an opacity as:

   .. math:: \kth = \frac{16 \sigma_B T^3}{3 \kappa_\mathrm{const} \rho}

   where :math:`\kappa_\mathrm{const}` is the opacity, with units :math:`\mathrm{cm^2~g^{-1}}`.
   This is selected by setting::

       CONDUCTIVITY_DIR := constant_opacity

   in your ``GNUmakefile``. To set the value of the opacity, e.g., to
   0.2 (for electron scattering), set::

       const_opacity = 0.2

   in the ``&extern`` namelist of your probin.

.. index:: castro.diffusion_cutoff_density, castro.diffusion_cutoff_density_hi

The diffusion approximation breaks down at the surface of stars,
where the density rapidly drops and the mean free path becomes
large. In those instances, you should use the flux limited diffusion
module in Castro to evolve a radiation field. However, if your
interest is only on the diffusion in the interior, you can use
the parameters:

 * ``castro.diffuse_cutoff_density``

 * ``castro.diffuse_cutoff_density_hi``

to specify a density,
below which, diffusion is not modeled. This is implemented in the
code by linearly scaling the conductivity to zero between these limits, e.g.,

.. math::

   \kth = \kth \cdot \frac{\rho - \mathtt{castro.diffuse\_cutoff\_density}}{\mathtt{castro.diffuse\_cutoff\_density\_hi} - \mathtt{castro.diffuse\_cutoff\_density}}


A simple test problem that sets up a Gaussian temperature profile
and does pure diffusion is provided as ``diffusion_test``.
