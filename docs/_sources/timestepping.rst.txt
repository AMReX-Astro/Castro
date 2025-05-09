************************
Timestepping and Retries
************************

Simulation Time
---------------

There are two parameters that can define when a simulation ends:

  * ``max_step``: maximum number of level 0 time steps (integer
    :math:`\geq 0`; default: -1)

  * ``stop_time``: final simulation time (Real :math:`\geq 0`; default:
    -1.0)

To control the number of time steps, you can limit by the maximum
number of level 0 time steps (``max_step``) or by the final
simulation time (``stop_time``), or both. The code will stop at
whichever criterion comes first.

Note that if the code reaches ``stop_time`` then the final time
step will be shortened so as to end exactly at ``stop_time``, not
past it.

As an example::

    max_step  = 1000
    stop_time  = 1.0

will end the calculation when either the simulation time reaches 1.0 or
the number of level 0 steps taken equals 1000, whichever comes first.

Time Step Constraints
---------------------

Hydrodynamics
^^^^^^^^^^^^^

If ``castro.do_hydro = 1``, then typically
the code chooses a time step based on the CFL number:

.. math::
   \Delta t = \mathtt{CFL}\, \cdot\, \min_{i,j,k}\left[\min\left\{\frac{\Delta x}{|u|_{i,j,k}+c_{i,j,k}},
                                                                  \frac{\Delta y}{|v|_{i,j,k}+c_{i,j,k}},
                                                                  \frac{\Delta z}{|w|_{i,j,k}+c_{i,j,k}}\right\}\right]
   :label: eq:cfl

If SDC integration is used instead, then we have

.. math::

   \Delta t = \mathtt{CFL}\, \cdot\, \min_{i,j,k}\left[\left(\frac{\Delta x}{|u|_{i,j,k}+c_{i,j,k}}\right)^{-1} +
                                                       \left(\frac{\Delta y}{|v|_{i,j,k}+c_{i,j,k}}\right)^{-1} +
                                                       \left(\frac{\Delta z}{|w|_{i,j,k}+c_{i,j,k}}\right)^{-1}\right]^{-1}

(If we are simulating in 1D or 2D, the extraneous parts related to :math:`v` and/or :math:`w` are removed.)

Additional Controls
^^^^^^^^^^^^^^^^^^^

The following parameters affect the timestep choice:

  * ``castro.cfl``: CFL number (Real :math:`> 0` and :math:`\leq 1`;
    default: 0.8)

  * ``castro.init_shrink``: factor by which to shrink the initial
    time step (Real :math:`> 0` and :math:`\leq 1`; default: 1.0)

  * ``castro.change_max``: factor by which the time step can
    grow in subsequent steps (Real :math:`\geq 1`; default: 1.1)

  * ``castro.fixed_dt``: level 0 time step regardless of cfl
    or other settings (Real :math:`> 0`; unused if not set)

  * ``castro.initial_dt``: initial level 0 time
    step regardless of other settings (Real :math:`> 0`; unused if not set)

  * ``castro.dt_cutoff``: as a fraction of the current simulation time,
    the time step below which the calculation will abort (Real
    :math:`> 0`; default: 1.e-12); typically not user-defined

As an example, consider::

    castro.cfl = 0.9
    castro.init_shrink = 0.01
    castro.change_max = 1.1

This defines the :math:`\mathtt{cfl}` parameter in :eq:`eq:cfl` to be
0.9, but sets (via ``init_shrink``) the first timestep we take to
be 1% of what it would be otherwise. This allows us to ramp up to
the hydrodynamic timestep at the start of a simulation. The
``change_max`` parameter restricts the timestep from increasing by
more than 10% over a coarse timestep. Note that the time step can
shrink by any factor; this only controls the extent to which it can
grow. The ``dt_cutoff`` parameter will force the code to abort if
the timestep ever drops below :math:`10^{-12}` of the current time. This is a safety
feature—if the code hits such a small value, then something likely
went wrong in the simulation, and by aborting, you won’t burn through
your entire allocation before noticing that there is an issue.

If we know what we are doing, then we can force a particular timestep::

    castro.fixed_dt = 1.e-4

This sets the level 0 time step to be 1.e-4 for the entire simulation,
ignoring the other timestep controls. Note that if
``castro.init_shrink`` :math:`\neq 1` then the first time step will in fact
be ``castro.init_shrink`` :math:`\cdot` ``castro.fixed_dt``.

::

    castro.initial_dt = 1.e-4

sets the *initial* level 0 time step to be :math:`10^{-4}` regardless of
``castro.cfl`` or ``castro.fixed_dt``. The time step can
grow in subsequent steps by a factor of castro.change_max each step.


Diffusion
^^^^^^^^^

If diffusion is enabled, the timestep will also be limited by:

.. math::

   \Delta t = \frac{1}{2}\min_{i,j,k}\left[\min\left\{\frac{\Delta x^2}{D_{i,j,k}},
                                                      \frac{\Delta y^2}{D_{i,j,k}},
                                                      \frac{\Delta z^2}{D_{i,j,k}}\right\}\right]

where :math:`D \equiv k / (\rho c_V)` if we are diffusing temperature,
and :math:`D \equiv k / (\rho c_P)` if we are diffusing enthalpy. No
input parameter is necessary to enable this constraint. See Chapter
:ref:`ch:diffusion` for more details.

Reactions
^^^^^^^^^

If reactions are enabled, the timestep will also
be limited by two constraints:

.. math:: \Delta t = \mathtt{dtnuc\_e}\, \min_{i,j,k} \left\{\frac{e_{i,j,k}}{\dot{e}_{i,j,k}}\right\}

.. math:: \Delta t = \mathtt{dtnuc\_X}\, \min_{i,j,k} \left\{\min_n\frac{X^n_{i,j,k}}{\dot{X}^n_{i,j,k}}\right\}

where :math:`e` is the internal energy, and :math:`X^n` is the mass fraction of
the :math:`n`\ th species. The safety factors correspond to the runtime parameters
``castro.dtnuc_e`` and ``castro.dtnuc_X``. These limiters
say that the timestep must be small enough so that no zone can change
its internal energy by more than the fraction in one
step, and so that no zone can change the abundance of any isotope by
more than the fraction in one step. The time derivatives
:math:`\dot{e}` and :math:`\dot{X}^n` are estimated by calling the right-hand-side
of the nuclear network given the state at the time the timestep limiter
is being calculated. (We use a small number floor to prevent division by zero.)
To prevent the timestep from being dominated by trace species, there is
an additional option ``castro.dtnuc_X_threshold`` which is the
mass fraction threshold below which a species will not be considered in
the timestep constraint. and are set to
a large number by default, effectively disabling them. Typical choices
for these values in the literature are :math:`\sim 0.1`.

Subcycling
----------

Subcycling with AMR means that coarser grids can take a larger timestep
than finer grids.
Castro supports a number of different modes for subcycling in time,
set via ``amr.subcycling_mode``.

  * ``amr.subcycling_mode`` = ``Auto`` (default): the code will run with
    equal refinement in space and time. In other words, if level
    :math:`n+1` is a factor of 2 refinement above level :math:`n`,
    then :math:`n+1` will take 2 steps of half the duration for every
    level :math:`n` step.

  * If ``amr.subcycling_mode`` = ``None``: the code will not refine in
    time. All levels will advance together with a timestep dictated by
    the level with the strictest :math:`dt`. Note that this is
    identical to the deprecated command ``amr.nosub = 1``.

  * If ``amr.subcycling_mode`` = ``Manual``: the code will subcycle
    according to the values supplied by ``amr.subcycling_iterations``.

In the case of ``amr.subcycling_mode`` = Manual, we subcycle in
manual mode with largest allowable timestep. The number of iterations
at each level is then specified as::

    amr.subcycling_iterations = 1 2 1 2

Here, we take 1 level-0 timestep at a time (required). Take 2 level-1
timesteps for each level-0 step, 1 timestep at level-2 for each
level-1 step, and take 2 timesteps at level-3 for each level-2 step.

Alternately, we could do::

    amr.subcycling_iterations = 2

which will subcycle twice at every level (except level 0).


.. index:: retry

.. _ch:retry:

Retry Mechanism
---------------

.. index:: castro.use_retry, castro.abundance_failure_tolerance, castro.abundance_failure_rho_cutoff, castro.retry_small_density_cutoff, castro.small_dens

Castro's Strang CTU solver has a retry mechanism that can discard a
time step on a level and restart with a smaller timestep, subcycling
within the level to make up the full time step needed for that level.
It is enabled by setting::

   castro.use_retry = 1

.. note::

   The Castro retry mechanism is enabled by default for CTU + Strang
   and Simplified SDC integration.

The number of subcycles to try in the level is controlled via the
``castro.max_subcycles`` parameter.  It is not really suggested to go
beyond ``16``---any more is usually an indication of a bigger problem.

A retry can be triggered by a number of conditions:

  * Exceeding the CFL condition for a level

  * A negative density is encountered.  This check can be disabled
    in low density regions by setting ``castro.retry_small_density_cutoff`` to the density below which we silently reset the density to
    ``castro.small_dens``.

  * The mass fractions fall outside of :math:`[0, 1]` -- we use
    ``castro.abundance_failure_tolerance`` with a default value of
    ``0.01`` to trigger the retry.  This check can be disabled at low
    densities by setting ``castro.abundance_failure_rho_cutoff`` to
    the density below which we want to silently renormalize the species.

  * Integration failure in the burner

    By construction, the integration routines in Microphysics will not
    abort if the integration fails, but instead return control to the
    calling function and set ``burn_t burn_state.success=false``.  This
    allows Castro to handle the failure.
