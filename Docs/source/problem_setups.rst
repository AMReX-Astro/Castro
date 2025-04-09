Distributed Problem Setups
==========================

There are a number of standard problem setups that come with Castro.
These can be used as a starting point toward writing your own setup.
We organize these into subdirectories by broad type (radiation, hydro,
gravity, etc.): The standard categories and *some* of the included
problems are:

* ``gravity_tests``:

   * ``DustCollapse``: A pressureless cloud collapse that is a
     standard test problem for gravity. An analytic solution that
     describes the radius of the sphere as a function of time is found
     in Colgate and White :cite:`colgwhite`. This problem is also
     found in the FLASH User’s Guide.

   * ``evrard_collapse``: This is the collapse of an isothermal
     spherical gas cloud.  This problem was originally discussed in
     :cite:`Evrard1988`.
     This implementation of the test comes from section 9.1 of
     :cite:`springel:2010`.


   * ``hydrostatic_adjust``: Model a 1-d stellar atmosphere (plane-parallel or
     spherical/self-gravitating) and dump energy in via an analytic
     heat source and watch the atmosphere’s hydrostatic state adjust
     in response. This is the counterpart to the Maestro
     ``test_basestate`` unit test.

   * ``hse_convergence``: This is meant to be a simple 1-d test for assessing the convergence of
     hydro + gravity in maintaining HSE.  Convergence can be measured either
     by looking at the max :math:`|U|` in the plotfiles.

   * ``hydrostatic_adjust``: This is a problem that explores the
     change in a hydrostatic structure due to heating.  This was used
     originally in :cite:`maestro:II`.

   * ``StarGrav``: This problem sets up a single spherical star in
     hydrostatic equilibrium and is used to assess the ability to
     maintain HSE.

   * ``uniform_cube_sphere``: This is used to compute the
     gravitational potential of a perfect cube, for which there is an
     analytic solution.  It tests our isolated boundary conditions.
     This was demonstrated in :cite:`katz:2016`.


* ``hydro_tests``:

   * ``acoustic_pulse``: The acoustic pulse problem from
     :cite:`mccorquodalecolella` used to measure convergence of pure
     hydrodynamics problems (as used for Castro in
     :cite:`castro-sdc`).

   * ``acoustic_pulse_general``: a general equation of state version
     of ``acoustic_pulse`` used for measuring convergence in
     :cite:`castro-sdc`.

   * ``double_bubble``: Initialize 1 or 2 bubbles in a stratified
     atmosphere (isothermal or isentropic) and allow for the bubbles
     to have the same or a different :math:`\gamma` from one another /
     the background atmosphere.  This uses the multigamma EOS.
     An analogous problem is implemented in Maestro.

   * ``KH``: A Kelvin-Helmholtz shear instability problem.

   * ``oddeven``: A grid-aligned shock hitting a very small density
     perturbation.  This demonstrates the odd-even decoupling problem
     discussed in :cite:`quirk1997`. This setup serves to test the
     castro.hybrid_riemann option to hydrodynamics.

   * ``RT``: A single-model Rayleigh-Taylor instability problem.

   * ``Sedov``: The standard Sedov-Taylor blast wave problem. This
     setup was used in the first Castro paper :cite:`castro_I`.

   * ``Sod``: A one-dimensional shock tube setup, including the
     classic Sod problem. This setup was used in the original Castro
     paper :cite:`castro_I`.

   * ``Sod_stellar``: A version of the Sod shock tube for the general
     stellar equation of state. This setup and the included inputs
     files was used in :cite:`zingalekatz`.

   * ``toy_convect``: A simple nova-like convection problem with an
     external heating source. This problem shows how to use the model
     parser to initialize a 1-d atmosphere on the Castro grid,
     incorporate a custom tagging routine, sponge the fluid above the
     atmosphere, and write a custom diagnostics routine.
     A MAESTROeX version of this problem setup also exists.

* ``mhd_tests``:

   * ``Alfven``: a linearized MHD wave test problem from :cite:`crockett:2005` and :cite:`miniati_martin`.

   * ``BrioWu``: the Brio Wu shock tube problem as described in :cite:`briowu`.  This is a standard
     test problem used in many MHD code papers (e.g. :cite:`athena`).

   * ``DaiWoodward``: a shock tube problem described in :cite:`Dai_1998`

   * ``FastRarefaction``: a shock tube problem dominated by kinetic energy, as described in :cite:`miniati_martin`

   * ``MagnetosonicWaves``: the fast and slow magnetosonic wave problem from :cite:`crockett:2005`

   * ``OrszagTang``: a two-dimensional magnetized vortex problem, following :cite:`athena`

   * ``RT``: a magnetized Rayleigh-Taylor instability problem

   * ``species``: a simple test problem to ensure that species are accurately advected.


* ``radiation_tests``:

   * ``Rad2Tshock``: This sets up a radiating shock that can be
     compared to a semi-analytic solution described in :cite:`lowrieedwards`.

   * ``RadFront``: This is the optically-thin streaming of a radiation front problem
     demonstrated originally in Castro in :cite:`CastroII`.

   * ``RadShestakovBolstad``: This is a linear multigroup diffusion test problem first described
     by :cite:`SHESTAKOV2005` and demonstrated in Castro in :cite:`CastroIII`.

   * ``RadSourceTest``: Test the implementation of the source terms in the gray radiation
     solver.  This does the "relaxation to thermal equilibrium" test as
     described in :cite:`swestymyra:2009`  (originally described in :cite:`turnerstone2001`).

   * ``RadSphere``: This is a multigroup radiating sphere test problem with an analytic solution,
     described in :cite:`graziani:2008` and :cite:`swestymyra:2009` and shown in Castro in :cite:`CastroIII`.

   * ``RadSuOlson``: This is a non-equlibrium Marshak wave test described in :cite:`suolson:1996` and shown
     in Castro in :cite:`CastroII`.

   * ``RadSuOlsonMG``: This is a multigroup version of ``RadSuOlson`` described in :cite:`suolson:1999`
     and shown in Castro in :cite:`CastroIII`.

   * ``RadThermalWave``: A thermal wave test adapted from :cite:`howellgreenough:2003` and shown in Castro
     in :cite:`CastroII`.

* ``reacting_tests``:

   * ``bubble_convergence``: a reacting bubble problem designed for measuring the convergence of
     the reactive hydro algorithms in Castro.  This was used in :cite:`castro-sdc`.

   * ``reacting_bubble``: A reacting bubble in a stratified white
     dwarf atmosphere. This problem was featured in the
     Maestro reaction paper :cite:`maestro:III`.

   * ``reacting_convergence``: a simple reacting hydrodynamics problem for measuring convergence,
     used in :cite:`castro-sdc` and :cite:`strang_rnaas`

* ``science``:

  The problems in the science directory are science problems that have
  appeared in papers (or will shortly).  Many of these are being actively used and are shared
  here for reproducibility.

   * ``Detonation``: this sets up a 1-d detonation that propagates through the domain.

   * ``flame``: this sets up a 1-d deflagration that propagates through the domain.  This setup
     was used for the testing in :cite:`eiden:2020`.

   * ``flame_wave``: this is a model of a flame propagating across a neutron star as a model for
     an X-ray burst.  This was presented in :cite:`eiden:2020` and :cite:`harpole2021dynamics`.

   * ``nova``: this models convection at the base of an accreted layer
     on a white dwarf as a model of a nova.

   * ``planet``: this is the problem setup from :cite:`ryu:2018` that models shear and turbulence in a
     hot Jupiter atmosphere.

   * ``subchandra``: a model of sub-Chandra Type Ia supernova that initializes a hot spot in a helium
     layer on a low mass carbon-oxygen white dwarf.

   * ``wdmerger``: a problem setup for modeling white dwarf mergers.  This was used in :cite:`katz:2016`.

   * ``xrb_mixed``: a compressible version of the X-ray burst convection problem from :cite:`zingale:2015`.

* ``unit_tests``:

   * ``diffusion_test``: a test of thermal diffusion (without hydro).  This was used to demonstrate convergence
     in both :cite:`castro-sdc` and :cite:`eiden:2020`.

   * ``particles_test``: a test of passive particles.

