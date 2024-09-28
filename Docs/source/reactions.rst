*********
Reactions
*********


.. index:: burn_t

The nuclear network serves two purposes: it defines the fluid
components used in both the equation of state and the hydrodynamics,
and it evolves those components through a nuclear burning step.  All
of the reaction networks that Castro uses are provided by the
`Microphysics repository <https://github.com/amrex-astro/Microphysics>`_.

.. index:: USE_REACT

.. note::

   To enable reactions in a simulation, you must compile with

   ::

      USE_REACT=TRUE

   This can be set in your ``GNUmakefile``.

Microphysics comes with a ``general_null``
network. This is a bare interface for a
nuclear reaction network. No reactions are enabled, and no auxiliary variables
are accepted.  It contains several sets of isotopes; for example,
``Microphysics/networks/general_null/triple_alpha_plus_o.net`` would describe the
isotopes needed to represent the triple-\ :math:`\alpha` reaction converting
helium into carbon, as well as oxygen and iron.

Other reaction networks can be found in ``Microphysics/networks``.  To select
a network, you set the make variable ``NETWORK_DIR`` to the name of the network
directory.

.. note::

   An arbitrary reaction network can be created for Castro via the
   `pynucastro library <https://pynucastro.github.io/pynucastro/>`_.


The main interface for burning is in ``Microphysics/interfaces/burner.H``:

.. code:: c++

   void burner (burn_t& state, Real dt)

Here the ``burn_t`` type contains all of the information needed for the reaction
network.  It is similar to the equation of state ``eos_t``.

.. tip::

   The equation of state routines can be called directly with a ``burn_t`` in place
   of an ``eos_t``.

Castro has several different modes of coupling reactions and
hydrodynamics, selected through the parameter
``castro.time_integration_method``.  See :ref:`sec:flowchart` for the
details.

Controlling burning
===================

.. index:: castro.react_T_min, castro.react_T_max, castro.react_rho_min, castro.react_rho_max, castro.do_react

There are a number of reactions-related parameters that can be set at runtime
in the inputs file. Reactions are enabled by setting::

    castro.do_react = 1

(Note: turning reactions off for problems where they're not required can help improve
the efficiency).

It is possible to set the maximum and minimum temperature and density for allowing
reactions to occur in a zone using the parameters:

* ``castro.react_T_min`` and ``castro.react_T_max`` for temperature

* ``castro.react_rho_min`` and ``castro.react_rho_max`` for density


Burning in Shocks
-----------------

.. index:: USE_SHOCK_VAR, castro.diable_shock_burning, castro.shock_detection_threshold, castro.shock_detection_include_sources

Burning can also be disabled inside shocks.  This requires that the code be
compiled with::

  USE_SHOCK_VAR = TRUE

in the ``GNUmakefile``.  This will allocate storage for a shock flag in the conserved
state array.  This flag is computed via a multidimensional shock detection algorithm
described in :cite:`doubledet2024`.  A zone is tagged as a shock if the following
conditions are true:

.. math::

   \begin{align*}
   \nabla \cdot \ub &< 0 \\
   \frac{|(\nabla p - \rho {\bf g}) \cdot \ub|}{p |\ub_\mathrm{cell}|} &> f_\mathrm{shock}
   \end{align*}

This requires that there is compression and that the pressure jump (excluding
the part of the pressure that balances gravity) is large.  The runtime parameter

::

   castro.disable_shock_burning = 1

will skip reactions in a zone where we've detected a shock.  The runtime parameters
``castro.shock_detection_threshold`` and ``castro.shock_detection_include_sources``
will set the value of $f_\mathrm{shock}$ and whether to subtract $\rho {\bf g}$
from the pressure gradient.

.. note::

   Both the compilation with ``USE_SHOCK_VAR = TRUE`` and the runtime parameter
   ``castro.disable_shock_burning = 1`` are needed to turn off burning in shocks.

Reactions Flowchart
===================

Here we describe how the ``burn_t`` is setup before the burn and how we update the
castro state afterwards for both Strang and simplified-SDC.

Strang
------

In ``Castro_react.cpp``, the flow is:

* create ``burn_t burn_state``

* if ``NSE_NET`` is defined, initialize the chemical potentials that
  will be used as an initial guess for the NSE solve

  * ``burn_state.mu_p`` $= U(\mu_p)$

  * ``burn_state.mu_n`` $= U(\mu_n)$

  * ``burn_state.y_e`` $= 0$ (this will be filled if needed by the NSE routines)

* initialize ``burn_state.dx`` -- this is used for some NSE conditions.

* set ``burn_state.success = true`` : we assume that the burn was successful.  The
  integrator will set this to ``false`` is a problem occurred.

* fill the thermodynamic quantities for input to the burner:

  * ``burn_state.rho`` $= U(\rho)$

  * ``burn_state.e`` $= U(\rho e) / U(\rho)$

  * ``burn_state.T`` $= U(T)$

    .. note::

       It is assumed here that the temperature is thermodynamically
       consistent with the energy.  For most networks, the temperature
       passed in will be used to set the thermodynamics in the burner.

  * ``burn_state.xn[]`` $= U(\rho X_k) / U(\rho)$

  * if ``NAUX_NET > 0``: ``burn_state.aux[]`` $= U(\rho \alpha_k) / U(\rho)$

* If we are doing ``castro.drive_initial_convection`` then we set
  ``burn_state.T_fixed`` by interpolating from the initial model.

* Initialize the metadata that is used for diagnostics

* Call the burner:

  * We check to make sure that $T$ and $\rho$ are within the limits given
    by ``castro.react_T_min``, ``castro.react_T_max``, ``castro_react_rho_min``,
    and ``castro.react_rho_max``.

  * The burner will set ``burn_state.success = false`` if it failed.  This can happen
    for a number of reasons and is integrator-dependent.

    .. note::

       Castro will not abort by default here if the burn failed.
       Instead we leave it to the :ref:`ch:retry` mechanism to attempt
       the step again with a smaller timestep.

* Store the burning sources for plotting

  .. index:: Reactions_Type

  We use the ``Reactions_Type`` ``StateData`` to hold the reactive
  sources that are output to the plotfile and the ``burn_weights``
  ``MultiFab`` to hold the number of righthand side evaluations for
  diagnostics.

  We fill these as:

  .. index:: castro.store_omega_dot

  * energy generation rate:

    $\mathtt{reactions}(\rho e) = \dfrac{U(\rho) \, \cdot\, \mathtt{burn\_state.e}\, -\, U(\rho e)}{\Delta t}$

  * species and auxiliary creation rates (only if ``castro.store_omegadot = 1``):

    * $\mathtt{reactions}(\rho X_k) = U(\rho) \dfrac{\mathtt{burn\_state.xn[k]}\, -\, U(\rho X_k) / U(\rho)}{\Delta t}$

    * $\mathtt{reactions}(\rho \alpha_k) = U(\rho) \dfrac{\mathtt{burn\_state.aux[k]}\, -\, U(\rho \alpha_k) / U(\rho)}{\Delta t}$

  * NSE flag (only if ``NSE`` is defined).  This simply stores the value of ``burn_state.nse``.

* Update the conserved state:

  .. note::

     $\rho$ and $\rho \ub$ are unchanged by reactions so those variables are not
     updated here.  They are already the "new" state.

  * $U^\mathrm{new}(\rho e) = U^\mathrm{new}(\rho) \cdot \mathtt{burn\_state.e}$

  * $U^\mathrm{new}(\rho E) = U^\mathrm{old}(\rho E) + (U^\mathrm{new}(\rho e) - U^\mathrm{old}(\rho e))$

  * $U^\mathrm{new}(\rho X_k) = U^\mathrm{new}(\rho) \cdot \mathtt{burn\_state.xn[k]}$

  * if ``NAUX_NET > 0``: $U^\mathrm{new}(\rho \alpha_k) = U^\mathrm{new}(\rho) \cdot \mathtt{burn\_state.aux[k]}$

  * if ``NSE_NET`` :

    * $U(\mu_p) = \mathtt{burn\_state.mu\_p}$

    * $U(\mu_n) = \mathtt{burn\_state.mu\_n}$



Simplified-SDC
--------------

In ``Castro_react.cpp``, the flow is:

* create ``burn_t burn_state``

* if ``NSE_NET`` is defined, initialize the chemical potentials that
  will be used as an initial guess for the NSE solve

  * ``burn_state.mu_p`` $= U(\mu_p)$

  * ``burn_state.mu_n`` $= U(\mu_n)$

  * ``burn_state.y_e`` $= 0$ (this will be filled if needed by the NSE routines)

* initialize ``burn_state.dx`` -- this is used for some NSE conditions.

* set ``burn_state.success = true`` : we assume that the burn was successful.  The
  integrator will set this to ``false`` is a problem occurred.

* fill the conserved state -- this is stored in the ``burn_t`` only when
  we are using simplified-SDC.

  * ``burn_state.y[SRHO]`` $= U(\rho)$

  * ``burn_state.y[SMX]`` $= U(\rho u)$

  * ``burn_state.y[SMY]`` $= U(\rho v)$

  * ``burn_state.y[SMZ]`` $= U(\rho w)$

  * ``burn_state.y[SEDEN]`` $= U(\rho E)$

  * ``burn_state.y[SEINT]`` $= U(\rho e)$

  * ``burn_state.y[SFS+k]`` $= U(\rho X_k)$ for $k = 0 \ldots N_{\mathrm{spec}} - 1$

  * if ``NAUX_NET > 0`` : ``burn_state.y[SFX+k]`` $= U(\rho \alpha_k)$ for $k = 0 \ldots N_{\mathrm{aux}} - 1$


* fill the thermodynamic quantities in the ``burn_t`` :

  * ``burn_state.rho`` $= U(\rho)$

  * ``burn_state.T`` $= U(T)$ -- this is mainly going to be used as an initial guess

  .. note::

     We don't initialize ``burn_state.xn[]`` or ``burn_state.aux[]``

  * if ``NAUX_NET > 0``: ``burn_state.aux[]`` $= U(\rho \alpha_k) / U(\rho)$

* If we are doing ``castro.drive_initial_convection`` then we set
  ``burn_state.T_fixed`` by interpolating from the initial model.

* Store the advective update that will be used during the SDC integration.

* Compute

* Initialize the metadata that is used for diagnostics

* Call the burner:

  * We check to make sure that $T$ and $\rho$ are within the limits given
    by ``castro.react_T_min``, ``castro.react_T_max``, ``castro_react_rho_min``,
    and ``castro.react_rho_max``.

  * The burner will set ``burn_state.success = false`` if it failed.  This can happen
    for a number of reasons and is integrator-dependent.

    .. note::

       Castro will not abort by default here if the burn failed.
       Instead we leave it to the :ref:`ch:retry` mechanism to attempt
       the step again with a smaller timestep.

* Store the burning sources for plotting

  .. index:: Reactions_Type

  We use the ``Reactions_Type`` ``StateData`` to hold the reactive
  sources that are output to the plotfile and the ``burn_weights``
  ``MultiFab`` to hold the number of righthand side evaluations for
  diagnostics.

  We fill these as:

  .. index:: castro.store_omega_dot

  * energy generation rate:

    $\mathtt{reactions}(\rho e) = \dfrac{U(\rho) \, \cdot\, \mathtt{burn\_state.e}\, -\, U(\rho e)}{\Delta t}$

  * species and auxiliary creation rates (only if ``castro.store_omegadot = 1``):

    * $\mathtt{reactions}(\rho X_k) = U(\rho) \dfrac{\mathtt{burn\_state.xn[k]}\, -\, U(\rho X_k) / U(\rho)}{\Delta t}$

    * $\mathtt{reactions}(\rho \alpha_k) = U(\rho) \dfrac{\mathtt{burn\_state.aux[k]}\, -\, U(\rho \alpha_k) / U(\rho)}{\Delta t}$

  * NSE flag (only if ``NSE`` is defined).  This simply stores the value of ``burn_state.nse``.

* Update the conserved state:

  .. note::

     $\rho$ and $\rho \ub$ are unchanged by reactions so those variables are not
     updated here.  They are already the "new" state.

  * $U^\mathrm{new}(\rho e) = U^\mathrm{new}(\rho) \cdot \mathtt{burn\_state.e}$

  * $U^\mathrm{new}(\rho E) = U^\mathrm{old}(\rho E) + (U^\mathrm{new}(\rho e) - U^\mathrm{old}(\rho e))$

  * $U^\mathrm{new}(\rho X_k) = U^\mathrm{new}(\rho) \cdot \mathtt{burn\_state.xn[k]}$

  * if ``NAUX_NET > 0``: $U^\mathrm{new}(\rho \alpha_k) = U^\mathrm{new}(\rho) \cdot \mathtt{burn\_state.aux[k]}$

  * if ``NSE_NET`` :

    * $U(\mu_p) = \mathtt{burn\_state.mu\_p}$

    * $U(\mu_n) = \mathtt{burn\_state.mu\_n}$


