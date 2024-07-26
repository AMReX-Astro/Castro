.. _ch:radiation:

*********
Radiation
*********

Introduction
============

Castro has three radiation solvers:

-  SingleGroupSolver: this solver does not have radiation
   pressure. It is pure hydro plus radiation diffusion. This is only
   applicable when the medium is optically thick and the pressure is small.

-  SGFLDSolver: this is the gray flux-limited diffusion
   radiation hydrodymamics solver. Here the radiation pressure is
   separate from the gas pressure, and both pressures participate in
   the Riemann solver.

-  MGFLDSolver: this is the multigroup flux-limited diffusion
   radiation hydrodynamics solver. As with the gray solver, radiation
   pressure contributes to the pressure in the Riemann solver. Here a
   number of energy groups are used to represent the radiation field,
   and the opacity can be frequency-dependent.

The gray solver has a comoving frame mode and a mixed frame mode,
whereas the MG solver uses the comoving frame approach. More details
about the formulation and algorithm can be found in the series of
Castro papers.

Getting Started
===============

Getting the Code
----------------

The Castro radiation solver is part of the main Castro git repo,
so you already have all the Castro code and problem setups
to exercise radiation. The only other requirement is a copy
of the Hypre library. Hypre provides the algebraic multigrid
solvers used by the implicit radiation update. You can get
a copy at https://github.com/hypre-space/hypre (the minimum
supported release version is 2.23.0). Their install
instructions describe what to do; we recommend using the autotools
and GNU Make build. On HPC clusters, you typically want to build
with the same compiler you're using to build Castro, and you also
want to make sure that the options you're using for Castro are
compatible with the Hypre options, in particular when it comes to
``USE_MPI``, ``USE_OMP``, and ``USE_CUDA``.

As an example, to build Hypre on Summit with MPI and CUDA, you
should load the ``gcc/7.4.0`` and ``spectrum-mpi`` modules and
then do the following from the Hypre ``src/`` directory,
replacing ``/path/to/Hypre/install`` with the target location
where you want the Hypre files to be installed.
::

   CUDA_HOME=$OLCF_CUDA_ROOT HYPRE_CUDA_SM=70 CXX=mpicxx CC=mpicc FC=mpifort ./configure --prefix=/path/to/Hypre/install --with-MPI --with-cuda --enable-unified-memory
   make install

Then, when you are building Castro, you would build with
``USE_MPI=TRUE`` and ``USE_CUDA=TRUE``.

Castro looks for Hypre in the environment variable ``HYPRE_DIR``,
which you should point to the install directory you chose above.
Other than that, the only difference for builds with radiation
is that you must set ``USE_RAD=TRUE``.

Microphysics: EOS, Network, and Opacity
=======================================

EOS
---

Castro provides several types of equation of state (EOS), including
gamma-law and Helmholtz. To use the gamma-law EOS, set

::

    EOS_DIR := gamma_law

in the GNUmakefile.

The original Helmholtz EOS for stellar interiors includes a radiation
contribution. However, for radiation hydrodynamics calculations, the
radition contribution should be taken out of EOS because radiation has
been treated in other places. To use Helmholtz EOS, we will use the
version in Microphysics, as with the pure hydrodynamics code, but
this will interpret the RADIATION preprocessor variable and
disable the radiation portion of the EOS. If you have your own EOS, you
can put it in Microphysics.

There is also an artificial EOS that is used for several test cases called

::

   EOS_DIR := rad_power_law

This EOS should only be used for pure radiation-diffusion tests (i.e.
``castro.do_hydro = 0``). It defines the specific heat as a power law,

   .. math:: c_v = \mathrm{const}\ \rho^m T^{-n}

Network
-------

The radiation solver uses the same networks as we saw for pure hydro,
so nothing needs to change here. Again, if you are not modeling
reactions, then the general_null network can be used to define
the appropriate composition for your problem.

Opacity
-------

The most commonly used opacity setup is

.. math::
   \kappa_\nu = \mathrm{const}\ \rho^{m} T^{-n} \nu^{p} ,
   :label: eq:kappa

where :math:`\kappa` is either the Planck or Rosseland absorption
coefficients, :math:`\rho` is density, :math:`T` is temperature, :math:`\nu` is
frequency, and :math:`m`, :math:`n` and :math:`p` are constants. For the gray solver,
:math:`p = 0`. If :eq:`eq:kappa` is sufficient, set

::

    Opacity_dir := rad_power_law

in your GNUmakefile. See § \ `3.3.1 <#sec:opacpars>`__ for instructions on how
to configure the parameters used for this opacity setup. If you would prefer a different
opacity mechanism, you will need to create your own opacity module by creating a new
directory in the Microphysics/opacity directory, creating the same set of subroutines
that the others have.

Some notes:

-  Here, :math:`\kappa` has units of :math:`\mathrm{cm}^{-1}`. Some papers or
   texts may instead have an implicit density factor in :math:`\kappa`,
   yielding units :math:`\mathrm{cm}^2~\mathrm{g}^{-1}`.

-  Castro allows for two temperatures (different radiation and gas
   temperature, so :math:`E_\mathrm{r} \ne a T_\mathrm{gas}^4`).
   Correspondingly, Castro cares about both the Planck mean,
   :math:`\kappa_P`, and Rosseland mean, :math:`\kappa_R`, opacities—these have
   different weightings.

   If we set :math:`\kappa_P \Delta x \gg 1` (:math:`\kappa_P` is really large),
   then the two temperatures become the same.

   If we set :math:`\kappa_P = \kappa_R`, then we can see how different the
   two temperatures are.

   In an optically thick medium, we would not expect the two temperatures
   to be very different.

.. _sec:opacpars:

Opacity Parameters
~~~~~~~~~~~~~~~~~~

The parameters describing the opacity include:

-  For the Planck opacity of the form in :eq:`eq:kappa`,
   the following parameters set the coefficient and exponents.
   These are set via the ``opacity`` namespace in the inputs file.
   ``opacity.const_kappa_p`` must be set positive to be used.

   -  ``opacity.const_kappa_p = -1.0``

   -  ``opacity.kappa_p_exp_m = 0.0``

   -  ``opacity.kappa_p_exp_n = 0.0``

   -  ``opacity.kappa_p_exp_p = 0.0``

-  For the Rosseland opacity of the form in :eq:`eq:kappa`,
   the following parameters set the coefficient and exponents.
   These are set via the ``opacity`` namespace in the inputs file.
   ``opacity.const_kappa_r`` must be set positive to be used.

   -  ``opacity.const_kappa_r = -1.0``

   -  ``opacity.kappa_r_exp_m = 0.0``

   -  ``opacity.kappa_r_exp_n = 0.0``

   -  ``opacity.kappa_r_exp_p = 0.0``

-  For the scattering coefficient of the form in :eq:`eq:kappa`,
   the following parameters set the coefficient and exponents.
   These are set via the ``opacity`` namespace in the inputs file.

   -  ``opacity.const_scatter = 0.0``

   -  ``opacity.scatter_exp_m = 0.0``

   -  ``opacity.scatter_exp_n = 0.0``

   -  ``opacity.scatter_exp_p = 0.0``

-  Since the formula above, :eq:`eq:kappa`, is non-physical and
   singular, we must set some floors in practice to prevent
   numerical issues. We have one floor for the opacity, which is
   applied to both the Planck and Rosseland opacities, and we
   also have a temperature floor.

   -  ``opacity.kappa_floor = 1.d-50``

   -  ``opacity.rad_temp_floor = 0.0``

-  ``radiation.do_kappa_stm_emission = 0``

   If it is 1, correction for stimulated emission is applied to Planck mean as
   follows

   .. math::

      \kappa = \mathrm{const}\ \rho^{m} T^{-n} \nu^{p}
          \left [1-\exp{\left (-\frac{h\nu}{k T} \right )} \right ].

Note that the unit for opacities is :math:`\mathrm{cm}^{-1}`. For
the gray solver, the total opacity in the diffusion coefficient is the sum
of kappa_r and scattering, whereas for the MG solver,
there are two possibilities. If const_kappa_r is greater than
0, then the total opacity is set by kappa_r alone, otherwise
the total opacity is the sum of kappa_p and scattering.

Radiation Solver Physics
========================

In this section, we list some radiation related parameters that you
can set in an inputs file. Here are some important parameters:

-  radiation.SolverType:

   Set it to 5 for the gray solver, and 6 for the MG solver.

-  castro.do_hydro

   Usually you want to set it to 1. If it is set to 0,
   hydro will be turned off, and the calculation will only solve
   radiation diffusion equation.

-  castro.do_radiation

   If it is 0, the calculation will be pure hydro.

Below are more parameters. For each parameter, the default value is
on the right-hand side of the equal sign.

.. _sec:bothpar:

Verbosity and I/O
-----------------

-  radiation.v = 0

   Verbosity

-  radiation.verbose = 0

   Verbosity

-  radiation.plot_lambda = 0

   If 1, save flux limiter in plotfiles.

-  radiation.plot_kappa_p = 0

   If 1, save Planck mean opacity in plotfiles.

-  radiation.plot_kappa_r = 0

   If 1, save Rosseland mean opacity in plotfiles.

-  radiation.plot_lab_Er = 0

   If 1, save lab frame radiation energy density in plotfiles.
   This flag is ignored when the mixed-frame gray solver is used.

-  radiation.plot_com_flux = 0

   If 1, save comoving frame radiation flux in plotfiles.

-  radiation.plot_lab_flux = 0

   If 1, save lab frame radiation flux in plotfiles.

.. _sec:fluxlimiter:

Flux Limiter and Closure
------------------------

-  radiation.limiter = 2

   Possible values are:

   -   0: No flux limiter

   -   2: Approximate limiter of Levermore & Pomraning

   -  12: Bruenn’s limiter

   -  22: Larsen’s square root limiter

   -  32: Minerbo’s limiter

-  radiation.closure = 3

   Possible values are:

   -  0: :math:`f = \lambda`, where :math:`f` is the scalar Eddington factor
      and :math:`\lambda` is the flux limiter.

   -  1: :math:`f = \frac{1}{3}`

   -  2: :math:`f = 1 - 2 \lambda`

   -  3: :math:`f = \lambda + (\lambda R)^2`, where :math:`R` is the radiation
      Knudsen number.

   -  4: :math:`f = \frac{1}{3} + \frac{2}{3} (\frac{F}{cE})^2`, where
      :math:`F` is the radiation flux, :math:`E` is the radiation energy density,
      and :math:`c` is the speed of light.

Note the behavior of the radiative flux in the optically thin and
optically thick limits. The flux limiter, :math:`\lambda = \lambda(R)`,
where

.. math:: R = \frac{|\nabla E_r^{(0)}|}{\chi_R E_r^{(0)}}

Regardless of the limiter chosen, when we are optically thick,
:math:`\chi_R \rightarrow \infty`, :math:`R \rightarrow 0`, and :math:`\lambda \rightarrow 1/3`.
The radiative flux then becomes

.. math::

   F_r^{(0)} = -\frac{c\lambda}{\chi_R} \nabla E_r^{(0)} \rightarrow
     \frac{1}{3} \frac{c}{\chi_R} \nabla E_r^{(0)}

And when we are optically thin, :math:`\chi_R \rightarrow 0`, :math:`R \rightarrow \infty`,
and :math:`\lambda \rightarrow 1/R = \chi_R E_r^{(0)}/{|\nabla E_r^{0}|}`, and
the radiative flux then becomes

.. math::

   F_r^{(0)} = -\frac{c\lambda}{\chi_R} \nabla E_r^{(0)} \rightarrow
     -\frac{c}{\chi_R}\frac{\chi_R E_r^{(0)}}{|\nabla E_r^{0}|}
       \nabla E_r^{(0)} = -c E_r^{0}

See Krumholz et al. 2007 for some discussion on this.

Boundary Conditions
-------------------

The following parameters are for the radiation boundary in the diffusion
equation. They do not affect hydrodynamic boundaries.

-  radiation.lo_bc

   This sets the action to take at the lower edge of the domain in
   each coordinate direction. Possible values are:

   -  101 *Dirichlet*:

      Specify the radiation energy density on the boundary.
      For gray radiation, this could be :math:`E_r = a T^4`.

      For multigroup radiation, Castro stores the energy density as
      :math:`\mathrm{erg}~\mathrm{cm}^{-3}`, so the total radiation energy
      can be found by simply summing over the groups. So if you want
      to set the radiation BCs using the Planck function, you simply
      multiply by the group width—see Exec/radiation_tests/RadSphere/Tools/radbc.f90
      for an example.

   -  102 *Neumann*:

      Here, you specify the radiation flux on the boundary. For gray
      radiation, this is the expression given in the gray Castro paper
      (Eq. 7, 8),

      .. math:: F_r = - \frac{c\lambda}{\kappa_R} \nabla E_r

      where :math:`\lambda` is the flux limiter.

      Note that if your boundary represents an incoming flux through
      a vacuum (like stellar irradiation), then :math:`\kappa \rightarrow 0`, leaving

      .. math:: F_r = -c E_r

      (see § \ `4.2 <#sec:fluxlimiter>`__) in that case.

   -  104 *Marshak* (vacuum):

      Here, you specify the incident flux and the outside is a vacuum.
      This differs from the Neumann condition because there is also a
      flux coming from inside, for the net flux across the boundary is
      different than the incident flux.

   -  105 *Sanchez-Pomraning*:

      This is a modified form of the Marshak boundary condition that works with FLD.
      This is like the Marshak condition, but :math:`\lambda = 1/3` is not assumed inside
      the boundary (optical thickness).

-  radiation.hi_bc

   See radiation.lo_bc.

-  radiation.lo_bcflag = 0 0 0

   If it is 0, bcval is used for that dimension, otherwise
   subroutine rbndry in RadBndry_1d.f90 is called to set
   boundary conditions.

-  radiation.hi_bcflag = 0 0 0

   See radiation.lo_bcflag

-  radiation.lo_bcval = 0.0 0.0 0.0

   The actual value to impose for the boundary condition type set by
   radiation.lo_bc. This parameter is interpreted differently
   depending on the boundary condition:

   -  Dirchlet: Dirichlet value of rad energy density

   -  Neumann: inward flux of rad energy

   -  Marshak: incident flux

   -  Sanchez-Pomraning: incident flux

-  radiation.hi_bcval = 0.0 0.0 0.0

   See radiation.lo_bcval

Convergence
-----------

For the gray solver, there is only one iteration in the scheme,
whereas for the MG solver, there are two iterations with an inner
iteration embedded inside an outer iteration. In the following, the
iteration in the gray solver will also be referred as the outer
iteration for convenience. The parameters for the inner iteration are
irrelevant to the gray solver.

radiation.maxiter = 50
    |
    | Maximal number of outer iteration steps.

radiation.miniter = 1
    |
    | Minimal number of outer iteration steps.

radiation.reltol = 1.e-6
    |
    | Relative tolerance for the outer iteration.

radiation.abstol = 0.0
    |
    | Absolute tolerance for the outer iteration.

radiation.maxInIter = 30
    |
    | Maximal number of inner iteration steps.

radiation.minInIter = 1
    |
    | Minimal number of inner iteration steps.

radiation.relInTol = 1.e-4
    |
    | Relative tolerance for the inner iteration.

radiation.absInTol = 0.0
    |
    | Absolute tolerance for the inner iteration.

radiation.convergence_check_type = 0
    |
    | For the MG solver only. This specify the way of checking the
      convergence of an outer iteration. Possible values are

    -  0: Check :math:`T`, :math:`Y_e`, and the residues of the equations for
       :math:`\rho e` and :math:`\rho Y_e`

    -  1: Check :math:`\rho e`

    -  2: Check the residues of the equations for :math:`\rho e` and :math:`\rho Y_e`

    -  3: Check :math:`T` and :math:`Y_e`

.. _sec:graypar:

Parameters for Gray Solver
--------------------------

radiation.comoving = 1
    |
    | Do we use the comoving frame approach?

radiation.Er_Lorentz_term = 1
    |
    | If the mixed-frame approach is taken, this parameter decides whether
      Lorentz transformation terms are retained.

radiation.delta_temp = 1.0
    |
    | This is used in computing numerical derivativas with respect to :math:`T`.
      So it should be a small number compared with :math:`T`, but not too small.

radiation.update_limiter = 1000
    |
    | Stop updating flux limiter after update_limiter iteration steps.

radiation.update_planck = 1000
    |
    | Stop updating Planck mean opacity after update_planck iteration steps.

radiation.update_rosseland = 1000
    |
    | Stop updating Rosseland mean opacity after update_rosseland iteration steps.

Grouping in the MG Solver
-------------------------

We provide two methods of setting up groups based upon logarithmic
spacing. In both methods, you must provide:

radiation.nGroups
    |
    | Number of groups.

radiation.lowestGroupHz
    |
    | Frequency of the lower bound for the first group.

In addition, if the parameter groupGrowFactor is provided, then
the first method will be used, otherwise the second method will be
used. In the first way, you must also provide firstGroupWidthHz
(the width of the first group). The width of other groups is set to
be groupGrowFactor times the width of its immediately preceding
group. In the second way, you must provide highestGroupHz as
the upper bound of the last group. It should be noted that
lowestGroupHz can be 0 in the first method, but not the second
method. However, when we compute the group-integrated Planck
function, the lower bound for the first group and the upper bound for
the last group are assumed to be 0 and :math:`\infty`, respectively.

.. _sec:mgpar:

Parameters for MG Solver
------------------------

radiation.delta_e_rat_dt_tol = 100.0
    |
    | Maximally allowed relative change in :math:`e` during one time step.

radiation.delta_T_rat_dt_tol = 100.0
    |
    | Maximally allowed relative change in :math:`T` during one time step.

radiation.delta_Ye_dt_tol = 100.0
    |
    | Maximally allowed absolute change in :math:`Y_e` during one tim estep.

radiation.fspace_advection_type = 2
    |
    | Possible value is 1 or 2. The latter is better.

radiation.integrate_Planck = 1
    |
    | If 1, integrate Planck function for each group. For the first
      group, the lower bound in the integration is assumed to be 0 no
      matter what the grouping is. For the last group, the upper bound in
      the integration is assumed to be :math:`\infty`.

radiation.matter_update_type = 0
    |
    | How to update matter. 0 is proabaly the best.

radiation.accelerate = 2
    |
    | The inner iteration of the MG solver usually requires an
      acceleration scheme. Choices are

    -  0: No acceleration

    -  1: Local acceleration

    -  2: Gray acceleration

radiation.skipAccelAllowed = 0
    |
    | If it is set to 1, skip acceleration if it does not help.

radiation.n_bisect = 1000
    |
    | Do bisection for the outer iteration after n_bisec iteration steps.

radiation.use_dkdT = 1
    |
    | If it is 1, :math:`\frac{\partial \kappa}{\partial T}` is retained in the
      Jacobi matrix for the outer (Newton) iteration.

radiation.update_opacity = 1000
    |
    | Stop updating opacities after update_opacity outer iteration steps.

radiation.inner_update_limiter = 0
    |
    | Stop updating flux limiter after inner_update_limiter inner
      iteration steps. If it is 0, the limiter is lagged by one outer
      iteration. If it is -1, the limiter is lagged by one time step. If
      the inner iteration has difficulty in converging, setting this
      parameter it to -1 can help. Since the flux limiter is only a
      kludge, it is justified to lag it.

.. _sec:hypre:

Linear System Solver
--------------------

There are a number of choices for the linear system solver. The
performance of the solvers usually depends on problems and the
computer. So it is worth trying a few solvers to find out which one
is best for your problem and computer.

radsolve.level_solver_flag: the linear solver
in Hypre to use. The available choices are:

-  0: SMG

-  1: PFMG (:math:`\ge` 2-d only)

-  100: AMG using ParCSR ObjectType

-  102: GMRES using ParCSR ObjectType

-  103: GMRES using SStruct ObjectType

-  104: GMRES using AMG as preconditioner

-  109: GMRES using Struct SMG/PFMG as preconditioner

-  150: AMG using ParCSR ObjectType

-  1002: PCG using ParCSR ObjectType

-  1003: PCG using SStruct ObjectType

As a general rule, the SMG is the most stable solver, but is usually
the slowest. The asymmetry in the linear system comes from the
adaptive mesh, so the PFMG should be your first choice. Note: in
you cannot use PFMG.

Setting this to 109 (GMRES using Struct SMG/PFMG as preconditioner)
should work reasonably well for most problems.

radsolve.maxiter (default: 40):
Maximal number of iteration in Hypre.

radsolve.reltol (default: 1.e-10):
Relative tolerance in Hypre

radsolve.abstol (default: 0):
Absolute tolerance in Hypre

radsolve.v (default: 0):
Verbosity

radsolve.verbos (default: 0):
Verbosity

habec.verbose (default: 0):
Verbosity for level_solver_flag :math:`<` 100

hmabec.verbose (default: 0):
Verbosity for level_solver_flag :math:`>=` 100

Output
======

Gray Solver
-----------

For the gray radiation solver, the radiation energy density is stored in plotfiles
as rad. Note that this quantity has units of :math:`\mathrm{erg~cm^{-3}}`, which
is different that the specify internal energy of the gas :math:`\mathrm{erg~g^{-1}}`.
