***********************
Hydrostatic Equilibrium
***********************

There are a few places where we take care to maintain HSE.  These differ
by integration method:

* CTU: this does characteristic tracing of the interface states and is
  used by Strang and simplified-SDC.  CTU requires that we predict the
  interface states at the midpoint in time -- this can make enforcing
  HSE much more difficult.

* true-SDC: this uses a method of lines time reconstruction of the interface
  states -- they are only spatially reconstructed, not extrapolated in time.


Reconstruction
==============

Piecewise linear
----------------

.. index:: castro.use_pslope, castro.pslope_cutoff_density


Piecewise linear reconstruction can be used with both CTU and
true-SDC.

For the 4th-order MC limiter, if we set ``castro.use_pslope = 1``,
then we subtract off the hydrostatic pressure before limiting.
The idea here is that the hydrostatic pressure balances gravity
so it is not available to generate waves or cause over/undershoots.
We define the excess pressure as :math:`\tilde{p}`:

.. math::

   \tilde{p}_i = p_i - p_{\mathrm{hse},i}

This is designed to work with HSE discretized as:

.. math::

   p_{i+1} = p_i + \frac{1}{4} \Delta x (\rho_i + \rho_{i+1}) (g_i + g_{i+1})

The 4th order MC limiter uses information in zones
:math:`\{i-2,i-1,i,i+1,i+2\}`.  We pick zone ``i`` as the reference
and integrate from ``i`` to the 4 other zone centers, and construct
:math:`\{\tilde{p}_{i-2}, \tilde{p}_{i-1}, \tilde{p}_{i}, \tilde{p}_{i+1}, \tilde{p}_{i+2}\}`.  We then limit on this, giving us the change in the excess
pressure over the zone, :math:`\Delta \tilde{p}_i`.  Finally, we
construct the total slope (including the hydrostatic part) as:

.. math::

   \Delta p_i = \Delta \tilde{p}_i + \rho_i g_i \Delta x

.. note::

   This can be disabled at low densities by setting ``castro.pslope_cutoff_density``.

The remainder of the PLM algorithm is unchanged.  Since this was just
the reconstruction part, this step is identical between CTU and true-SDC.



PPM
---

The PPM version of ``use_pslope`` is essentially the same as the PLM
version described above.  We first compute the excess pressure,
:math:`\tilde{p}`, and then fit a cubic to the 4 cells surrounding an
interface to get the initial interface values.  These become the left
and right values of the parabola in each zone.  The PPM limiting is
then done on the parabola, again working with :math:`\tilde{p}`.
Finally, the parabola values are updated to include the hydrostatic
pressure.

.. index:: castro.ppm_well_balanced

We can do better with PPM, and only use the perturbational pressure,
$\tilde{p}$, in the characteristic tracing and then add back the
hydrostatic pressure to the interface afterwards.  This is done via
``castro.ppm_well_balanced=1``.

Fully fourth-order method
-------------------------

The 4th order accurate true SDC solver does not appear to need any
well-balancing.  It maintains HSE to very high precision, as shown in
:cite:`castro-sdc`.


Boundary conditions
===================

reflecting
----------

For piecewise linear reconstruction, in ``slope.H``, when using the
4th order MC slopes (``castro.plm_limiter = 2``), we impose special
ghost cell values on the normal velocity at reflecting boundaries to
ensure that the velocity goes to 0 at the boundary. This follows
:cite:`saltzman:1994` , page 162 (but note that they have a sign
error).  Only the normal velocity is treated specially.




HSE
---

.. index:: castro.hse_zero_vels, castro.hse_reflect_vels, castro.hse_interp_temp, castro.hse_fixed_temp

For hydrostatic boundary conditions, we follow the method from
:cite:`ppm-hse`.  Essentially, this starts with the last
zone inside of the boundary and integrates HSE into the ghost cells,
keeping the density either constant or extrapolating it.  Together
the discretized HSE equation and the EOS yield the density and pressure
in the ghost cells.

.. warning::

   The HSE boundary condition only works with constant gravity at the moment.

To enable this, we set the appropriate boundary's ``xl_ext_bc_type``, ``xr_ext_bc_type``,
``yl_ext_bc_type``, ``yr_ext_bc_type``, ``zl_ext_bc_type``,
``zr_ext_bc_type`` to ``1``.

We then control the behavior via the following options.

For the velocity, we have:

* ``hse_zero_vels`` : all 3 components of the velocity in the ghost
  cell are set to ``0``.

* ``hse_reflect_vels`` : the normal velocity is reflected and the transverse
  velocity components are given a zero gradient.

If neither of these are set, then all components of the velocity are
simply given a zero gradient.

The temperature in the ghost cells is controlled by:

* ``hse_interp_temp`` : if this is set to ``1``, then we fill the
  temperatures in the ghost cells via linear extrapolation, using the
  2 interior zones just inside the domain.  Otherwise, we take the
  temperature in the ghost cells to be constant.

* ``hse_fixed_temp`` : if this is positive, then we set the
  temperature in the ghost cells to the value specified.  This
  requires ``hse_interp_temp = 0``.



Interface states at reflecting boundary
=======================================

For all methods, we enforce the reflecting condition on the interface
states directly by reflecting the state just inside the domain to
overwrite the state on the reflecting boundary just outside of the
domain.  This is done for all variables (flipping the sign on the
normal velocity state).  This is especially important for
reconstruction that used a one-sided stencil (like the 4th order
method).


Test problems
=============

``Castro/Exec/gravity_tests/hse_convergence_general`` can be used to
test the different HSE approaches.  This sets up a 1-d X-ray burst
atmosphere (based on the ``flame_wave`` setup).  Richardson
extrapolation can be used to measure the convergence rate (or just
look at how the peak velocity changes).
