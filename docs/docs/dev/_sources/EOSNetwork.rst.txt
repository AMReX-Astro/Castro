************
Microphysics
************

Equation of State
=================

Standard Castro EOSes
---------------------

Castro is written in a modular fashion so that the EOS and network
burning routines can be supplied by the user. However, for the
examples presented later we use several EOS and network routines
that come with the Microphysics distribution.

Castro relies on routines to calculate the equation of state (EOS)
of a fluid, as well as a species network to define the components of
the fluid. The network optionally has the ability to do nuclear burning,
but for this section its main purpose is in defining the species so that
the EOS can calculate fluid properties that depend on composition, such
as electron fraction.

By default, Castro comes with the ``gamma_law``
EOS. This represents a gamma law gas, with equation of state:

.. math:: p = (\gamma - 1) \rho e.

The gas is currently assumed to be monatomic and ideal. (Only a
restricted set of thermodynamic variables are actually calculated,
the minimum necessary for the hydrodynamics. A fuller set of
thermodynamic variables, for example the entropy from the
Sackur-Tetrode equation, are calculated in the ``gamma_law_general``
EOS inside the Microphysics repository.)

Runtime Parameters
------------------

When inverting the EOS (e.g. by calling `eos_input_re`), an initial guess for
the temperature is required. This guess is provided by the runtime parameter
`castro.T_guess`, and should be set to a sensible value for each problem
(it will vary depending on which EOS is used).

EOS Interfaces and Parameters
-----------------------------

.. index:: eos_t

Each EOS should have two main routines by which it interfaces to the
rest of Castro. At the beginning of the simulation, ``eos_init``
will perform any initialization steps and save EOS variables (mainly
``smallt``, the temperature floor, and ``smalld``, the
density floor). Then, whenever you want to call the EOS, use::

 call eos (eos_input, eos_state)

The first argument specifies the inputs to the EOS. The options
that are currently available are stored in
``EOS/eos_data.F90``, and are always a combination of two
thermodynamic quantities. For example, ``eos_input_rt`` means
that we call the EOS with :math:`\rho` (density) and :math:`T` (temperature)
and we expect the EOS to return the associated thermodynamic
quantities such as internal energy :math:`e` and entropy :math:`s`.

We note that for real (non-analytic) equations of state
in which :math:`\rho`, :math:`T` and species are the independent variables, such
as the Helmholtz EOS, ``eos_input_rt`` directly calls the EOS
and obtains the other thermodynamic variables. But for other inputs,
e.g. ``eos_input_re``, a Newton-Raphson iteration is performed
to find the density or temperature that corresponds to the given
input.

The eos_state variable is a Fortran derived type (``eos_t``, similar to
a C struct). It stores a complete set of thermodynamic
variables. When calling the EOS, you should first fill the variables
that are the inputs, for example with

::

      use eos_type_module
      ...
      type (eos_t) :: eos_state
      ...
      eos_state % rho = state(i,j,k,URHO)
      eos_state % T   = state(i,j,k,UTEMP)
      eos_state % e   = state(i,j,k,UEINT) / state(i,j,k,URHO)
      eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)
      eos_state % aux = state(i,j,k,UFX:UFX+naux-1) / state(i,j,k,URHO)

Whenever the ``eos_state`` type is initialized, the thermodynamic
state variables are filled with unphysical numbers. If you do not
input the correct arguments to match your input quantities, the EOS
will call an error. This means that it is good practice to fill the
quantities that will be iterated over with an initial guess. Indeed,
this initial guess is typically required for equations of state that
iterate over this variable, as the values they are initialized with
will likely not converge. Usually a prior value of the temperature or
density suffices if it’s available, but if not then use ``T_guess`` or
``small_dens``.

If you are interested in using more realistic and sophisticated equations of
state, you should download the `Microphysics <https://github.com/starkiller-astro/Microphysics>`__
repository. This is a collection of microphysics routines that are compatible with the
BoxLib codes. We refer you to the documentation in that repository for how to set it up
and for information on the equations of state provided. That documentation
also goes into more detail about the details of the EOS code, in case you are interested in
how it works (and in case you want to develop your own EOS).

Required Thermodynamics Quantities
----------------------------------

Three input qwuantities are required of any EOS:

-  ``eos_input_re``: :math:`\rho`, :math:`e`, and :math:`X_k` are input

-  ``eos_input_rt``: :math:`\rho`, :math:`T`, and :math:`X_k` are input

-  ``eos_input_rp``: :math:`\rho`, :math:`P`, and :math:`X_k` are input

The ``eos_t`` derived type holds a large number of thermodynamics
quantities, but not all of these are needed for basic
Castro operation. The main quantities that any EOS in any mode needs to
supply, if they are not input, are:

-  ``eos_state % T``: the temperature

-  ``eos_state % P``: total pressure

-  ``eos_state % e``: the specific energy

-  ``eos_state % gam1``: the first adiabatic index,
   :math:`\Gamma_1 = d\log P / d\log \rho |_s`

Additionally the ``eos_input_re`` mode also needs to supply:

-  ``eos_state % cs``: the adiabatic sound speed

-  ``eos_state % dpdr_e``: the derivative, :math:`\partial p/\partial \rho |_e`
   — note that the specific internal energy, :math:`e`
   is held constant here.

-  ``eos_state % dpde``: the derivative, :math:`\partial p / \partial e |_\rho`

For radiation hydro, the ``eos_input_rt`` model needs to supply:

-  ``eos_state % cv``: the specific heat capacity.

Other quantities (e.g., entropy) might be needed for the derived
variables that are optional output into the plotfiles.


Composition derivatives
-----------------------

.. index:: eos_xderivs_t

A separate type, ``eos_xderivs_t`` provides access to derivatives with respect to mass fraction.

-  ``eos_xderivs % dhdX(nspec)``: the derivative of the
   specific enthalpy with respect to mass fraction at constant
   :math:`T` and :math:`p`:

   .. math:: \xi_k = e_{X_k} + \frac{1}{p_\rho} \left (\frac{p}{\rho^2} - e_\rho \right ) p_{X_k}

-  ``eos_xderivs % dpdx(nspec)``: the derivative of the pressure with respect to mass fraction:

   .. math::

      \begin{align}
      p_{X_k} &= \left .\frac{\partial p}{\partial \bar{A}} \right |_{\rho, T, \bar{Z}}
                \frac{\partial \bar{A}}{\partial X_k} +
                \left . \frac{\partial p}{\partial \bar{Z}} \right |_{\rho, T, \bar{A}}
                \frac{\partial \bar{Z}}{\partial X_k} \nonumber \\
              &= -\frac{\bar{A}^2}{A_k}
                \left .\frac{\partial p}{\partial \bar{A}} \right |_{\rho, T, \bar{Z}} +
                \frac{\bar{A}}{A_k} \left (Z_k - \bar{Z} \right )
                \left . \frac{\partial p}{\partial \bar{Z}} \right |_{\rho, T, \bar{A}}
      \end{align}

-  ``eos_xderivs % dedx(nspec)``: the derivative of the specific internal energy with respect to mass fraction:

   .. math::

      \begin{align}
      e_{X_k} &= \left . \frac{\partial e }{\partial \bar{A}} \right |_{\rho, T, \bar{Z}}
              \frac{\partial \bar{A}}{\partial X_k} +
              \left .\frac{\partial e}{\partial \bar{Z}} \right |_{\rho, T, \bar{A}}
              \frac{\partial \bar{Z}}{\partial X_k} \nonumber \\
              &= -\frac{\bar{A}^2}{A_k}
              \left . \frac{\partial e }{\partial \bar{A}} \right |_{\rho, T, \bar{Z}} +
              \frac{\bar{A}}{A_k} \left (Z_k - \bar{Z}\right )
              \left .\frac{\partial e}{\partial \bar{Z}} \right |_{\rho, T, \bar{A}}
      \end{align}

(see :cite:`maestro:III`, Appendix A).


Nuclear Network
===============

.. index:: burn_t

The nuclear network serves two purposes: it defines the fluid components used
in both the equation of state and the hydrodynamics, and it evolves those
components through a nuclear burning step. Castro comes with a ``general_null``
network (which lives in the ``networks/`` directory). This is a bare interface for a
nuclear reaction network. No reactions are enabled, and no auxiliary variables
are accepted.  It contains several sets of isotopes; for example,
``networks/general_null/triple_alpha_plus_o.net`` would describe the
isotopes needed to represent the triple-\ :math:`\alpha` reaction converting
helium into carbon, as well as oxygen and iron.

The main interface file, ``network.f90``, is a wrapper function. The
actual network details are defined in ``actual_network.f90``, a
file which is automatically generated in your work directory when you compile.
It supplies the number and names of species and auxiliary variables, as
well as other initializing data, such as their mass numbers, proton numbers,
and the binding energies.

The burning front-end interface, ``networks/burner.f90``, accepts a different
derived type called the ``burn_t`` type. Like the ``eos_t``, it has entries
for the basic thermodynamic quantities:

::

      use burn_type_module
      ...
      type (burn_t) :: burn_state
      ...
      burn_state % rho = state(i,j,k,URHO)
      burn_state % T   = state(i,j,k,UTEMP)
      burn_state % e   = state(i,j,k,UEINT) / state(i,j,k,URHO)
      burn_state % xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)

It takes in an input ``burn_t`` and returns an output ``burn_t`` after
the burning has completed. The nuclear energy release can be computed by
taking the difference of ``burn_state_out % e`` and
``burn_state_in % e``. The species change can be computed analogously.
In normal operation in Castro  the integration occurs over a time interval
of :math:`\Delta t/2`, where :math:`\Delta t` is the hydrodynamics timestep.

If you are interested in using actual nuclear burning networks,
you should download the `Microphysics <https://github.com/starkiller-astro/Microphysics>`__
repository. This is a collection of microphysics routines that are compatible with the
AMReX Astro codes. We refer you to the documentation in that repository for how to set it up
and for information on the networks provided. That documentation
also goes into more detail about the details of the network code, in case you are interested in
how it works (and in case you want to develop your own network).


Controlling burning
-------------------

There are a number of reactions-related parameters that can be set at runtime
in the inputs file. Reactions are enabled by setting::

    castro.do_react = 1

(Note: turning reactions off for problems where they're not required can help improve
the efficiency).

It is possible to set the maximum and minimum temperature and density for allowing
reactions to occur in a zone using the parameters ``castro.react_T_min``,
``castro.react_T_max``, ``castro.react_rho_min`` and ``castro.react_rho_max``.

