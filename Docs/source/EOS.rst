*****************
Equation of State
*****************

Castro is written in a modular fashion so that the EOS
can be supplied by the user.   No equations of state
are distributed with Castro, instead they are part
of the separate `Microphysics repository <https://github.com/amrex-astro/Microphysics>`_.

Most equations of state are written to take :math:`(\rho, T, X_k)` as
input and return the needed thermodynamic quantities.  For other
inputs, a Newton-Raphson root find is done.

Most of the standard problem setups in Castro (such as the Sedov blast wave)
use the ``gamma_law`` EOS. This represents a gamma law gas, with equation of state:

.. math:: p = (\gamma - 1) \rho e.

The gas is currently assumed to be monatomic and ideal.

Runtime Parameters
==================

When inverting the EOS (e.g. by using ``eos_input_re``), an initial guess for
the temperature is required. This guess is provided by the runtime parameter
``castro.T_guess``, and should be set to a sensible value for each problem
(it will vary depending on which EOS is used).

EOS Interfaces and Parameters
=============================

.. index:: eos_t

Each EOS should have two main routines by which it interfaces to the
rest of Castro. At the beginning of the simulation, ``eos_init``
will perform any initialization steps and save EOS variables (mainly
``smallt``, the temperature floor, and ``smalld``, the
density floor). Then, whenever you want to call the EOS, use:

.. code:: c++

   eos (eos_input, eos_state)

The first argument specifies the inputs to the EOS. The options
that are currently available are stored in Microphysics in
``interfaces/eos_type.H``, and are always a combination of two
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

The eos_state variable is a C struct, ``eos_t``. It stores a complete
set of thermodynamic
variables. When calling the EOS, you should first fill the variables
that are the inputs, for example with

.. code:: c++

   eos_t eos_state;
   ...
   eos_state.rho = state(i,j,k,URHO);
   eos_state.T   = state(i,j,k,UTEMP);
   eos_state.e   = state(i,j,k,UEINT) / state(i,j,k,URHO);
   for (int n = 0; n < NumSpec; ++n) {
       eos_state.xn[n] = state(i,j,k,UFS+n) / state(i,j,k,URHO);
   }
   for (int n = 0; n < NumAux; ++n) {
       eos_state.aux[n] = state(i,j,k,UFX+n) / state(i,j,k,URHO);
   }

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


Required Thermodynamics Quantities
==================================

Three input quantities are required of any EOS:

-  ``eos_input_re``: :math:`\rho`, :math:`e`, and :math:`X_k` are input

-  ``eos_input_rt``: :math:`\rho`, :math:`T`, and :math:`X_k` are input

-  ``eos_input_rp``: :math:`\rho`, :math:`P`, and :math:`X_k` are input

The ``eos_t`` derived type holds a large number of thermodynamics
quantities, but not all of these are needed for basic
Castro operation. The main quantities that any EOS in any mode needs to
supply, if they are not input, are:

-  ``eos_state.T``: the temperature

-  ``eos_state.p``: total pressure

-  ``eos_state.e``: the specific energy

-  ``eos_state.gam1``: the first adiabatic index,
   :math:`\Gamma_1 = d\log P / d\log \rho |_s`

Additionally the ``eos_input_re`` mode also needs to supply:

-  ``eos_state.cs``: the adiabatic sound speed

-  ``eos_state.dpdr_e``: the derivative, :math:`\partial p/\partial \rho |_e`
   — note that the specific internal energy, :math:`e`
   is held constant here.

-  ``eos_state.dpde``: the derivative, :math:`\partial p / \partial e |_\rho`

For radiation hydro, the ``eos_input_rt`` model needs to supply:

-  ``eos_state.cv``: the specific heat capacity.

Other quantities (e.g., entropy) might be needed for the derived
variables that are optional output into the plotfiles.


Composition derivatives
=======================

.. index:: eos_xderivs_t

A separate type, ``eos_xderivs_t`` provides access to derivatives with respect to mass fraction.

-  ``eos_xderivs.dhdX[NumSpec]``: the derivative of the
   specific enthalpy with respect to mass fraction at constant
   :math:`T` and :math:`p`:

   .. math:: \xi_k = e_{X_k} + \frac{1}{p_\rho} \left (\frac{p}{\rho^2} - e_\rho \right ) p_{X_k}

-  ``eos_xderivs.dpdX[NumSpec]``: the derivative of the pressure with respect to mass fraction:

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

-  ``eos_xderivs.dedX[NumSpec]``: the derivative of the specific internal energy with respect to mass fraction:

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
