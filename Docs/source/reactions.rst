*********
Reactions
*********


.. index:: burn_t

The nuclear network serves two purposes: it defines the fluid
components used in both the equation of state and the hydrodynamics,
and it evolves those components through a nuclear burning step.  All
of the reaction networks that Castro uses are provided by the
[Microphysics
repository](https://github.com/amrex-astro/Microphysics).

.. note::

   An arbitrary reaction network can be created for Castro via the
   [pynucastro library](https://pynucastro.github.io/pynucastro/).

Microphysics comes with a ``general_null``
network. This is a bare interface for a
nuclear reaction network. No reactions are enabled, and no auxiliary variables
are accepted.  It contains several sets of isotopes; for example,
``Microphysics/networks/general_null/triple_alpha_plus_o.net`` would describe the
isotopes needed to represent the triple-\ :math:`\alpha` reaction converting
helium into carbon, as well as oxygen and iron.

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

.. index:: castro.react_T_min, castro.react_T_max, castro.react_rho_min, castro.react_rho_max

There are a number of reactions-related parameters that can be set at runtime
in the inputs file. Reactions are enabled by setting::

    castro.do_react = 1

(Note: turning reactions off for problems where they're not required can help improve
the efficiency).

It is possible to set the maximum and minimum temperature and density for allowing
reactions to occur in a zone using the parameters ``castro.react_T_min``,
``castro.react_T_max``, ``castro.react_rho_min`` and ``castro.react_rho_max``.


