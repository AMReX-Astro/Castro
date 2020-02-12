***************************
Plotfiles and Visualization
***************************

There are a large number of tools that can be used to read in Castro
or AMReX data and make plots.  These tools all work from Castro
plotfiles.  Here we give an overview of the variables in plotfiles and
controlling their output, as well as some of the tools that can be used
for visualization.


Plotfile Outputting
===================

Castro has two levels of plotfiles, `regular` plotfiles and `small`
plotfiles.  The idea behind this distinction is that we can output a
small number of variables very frequently in the small plotfiles and
output a large number (or all variables) less frequently.  This helps
keep the data sizes down while allowing for fine-grained temporal
analysis of important quantities.

.. index:: amr.plot_files_output, amr.plotfile_on_restart, amr.write_plotfile_with_checkpoint

A few general controls determines whether we want to output plotfiles and when:

  * ``amr.plot_files_output`` : this is set to 1 to output plotfiles

  * ``amr.plotfile_on_restart`` : set this to 1 to dump out a plotfile
    immediately when we restart.

  * ``amr.write_plotfile_with_checkpoint`` : always output a plotfile
    when we dump a checkpoint file.

.. index:: amr.plot_file, amr.plot_per, amr.plot_int

The frequency of outputting and naming of regular plotfiles is
controlled by:

  * ``amr.plot_file`` : this is the base name for the plotfile,
    e.g. ``plt``.

  * ``amr.plot_per`` : this is the amount of simulation time between
    plotfile output

  * ``amr.plot_int`` this is the number of timesteps between plotfiles.
    Set this to -1 to rely on the simulation-time-based outputting.

.. index:: amr.small_plot_file, amr.small_plot_per, amr.small_plot_int

Similarly, the frequency of outputting and naming of small plotfiles
is controlled by:

  * ``amr.small_plot_file`` : this is the base name for the small plotfile,
    e.g. ``smallplt``.

  * ``amr.small_plot_per`` : this is the amount of simulation time between
    small plotfile output

  * ``amr.small_plot_int`` this is the number of timesteps between small plotfiles.
    Set this to -1 to rely on the simulation-time-based outputting.


Controlling What’s in the PlotFile
==================================

.. index:: amr.plot_vars, amr.derive_plot_vars

There are a few options that can be set at runtime to control what
variables appear in the regular plotfile.

  * ``amr.plot_vars``: this controls which of the main
    state variables appear in the plotfile. The default is for all of
    them to be stored. But you can specify a subset by name, e.g.::

        amr.plot_vars = density

    to only store that subset.

  * ``amr.derive_plot_vars``: this controls which of the derived
    variables to be stored in the plotfile. Derived variables are
    created only when the plotfile is being created, using the
    infrastructure provided by AMReX to register variables and the
    associated Fortran routine to do the deriving (``Derive_nd.F90``).

    By default, no derived variables are stored. You can store all
    derived variables that Castro knows about by doing::

       amr.derive_plot_vars = ALL

   or a subset by explicitly listing them, e.g.::

      amr.derive_plot_vars = entropy pressure

   To not output any derived variable,s this is set to ``NONE``.

.. index:: amr.small_plot_vars

For small plotfiles, the controls that lists the variables is:

  * ``amr.small_plot_vars`` : this is a list of which variables
    to include in the small plotfile.

  * ``amr.derive_small_plot_vars`` : this is a list of which derived
    variables to include in the small plotfile.


Plotfile Variables
==================

Native variables
----------------

These variables come directly from the ``StateData``, either the
``State_Type`` (for the hydrodynamic variables), ``Reactions_Type``
(for the nuclear energy generation quantities). ``PhiGrav_Type`` and
``Gravity_Type`` (for the gravity quantities), ``PhiRot_Type`` and
``Rotation_Type`` (for the rotation quantities) and ``Rad_Type`` (for
radiation quantities).


+-----------------------------------+---------------------------------------------------+--------------------------------------+
| variable name                     | description                                       | units                                |
+===================================+===================================================+======================================+
| ``density``                       | Mass density, :math:`\rho`                        | :math:`\gcc`                         |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``xmom``                          | x-momentum, :math:`(\rho u)`                      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``ymom``                          | y-momentum, :math:`(\rho v)`                      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``zmom``                          | z-momentum, :math:`(\rho w)`                      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_E``                         | Total energy density                              | :math:`{\rm erg~cm^{-3}}`            |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_e``                         | Internal energy density                           | :math:`{\rm erg~cm^{-3}}`            |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Temp``                          | Temperature                                       | :math:`{\rm K}`                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_X``                         | Mass density of species X                         | :math:`\gcc`                         |
| (where X is any of the species    |                                                   |                                      |
| defined in the network)           |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``omegadot_X``                    | Creation rate of species X                        | :math:`{\rm s^{-1}}`                 |
| (where X is any of the species    | :math:`\omegadot_k = DX_k/Dt`                     |                                      |
| defined in the network)           |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``enuc``                          | Nuclear energy generation rate / gram             | :math:`{\rm erg~g^{-1}~s^{-1}}`      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_enuc``                      | Nuclear energy generation rate density            | :math:`{\rm erg~cm^{-3}~s^{-1}}`     |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiGrav``                       | Gravitational potential                           | :math:`{\rm erg~g^{-1}}`             |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``grav_x``, ``grav_y``,           | Gravitational acceleration                        | :math:`{\rm cm~s^{-2}}`              |
| ``grav_z``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiRot``                        | Effective centrifugal potential                   | :math:`{\rm erg~g^{-1}}`             |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rot_x``. ``rot_y``, ``rot_z``   | Rotational acceleration                           | :math:`{\rm cm~s^{-2}}`              |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rmom``                          | Radial momentum (defined for                      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
|                                   | ``HYBRID_MOMENTUM``)                              |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``lmom``                          | Angular momentum (:math:`\theta`; defined for     | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
|                                   | ``HYBRID_MOMENTUM``)                              |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``pmom``                          | z-momentum (defined for ``HYBRID_MOMENTUM``)      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Shock``                         | Shock flag (= 1 if a zone has a shock;            | --                                   |
|                                   | defined for ``SHOCK``)                            |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rad``, ``rad0``, ``rad1``,      | Radiation energy density                          |                                      |
| ...                               | (for multigroup radiation, each group has its     |                                      |
|                                   | own variable)                                     |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+



Derived variables
-----------------

+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| variable name                     | description                                       | derive routine              | units                                   |
+===================================+===================================================+=============================+=========================================+
| ``angular_momentum_x``,           | Angular momentum / volume in the x, y, or z dir   | ``derangmomx``,             | :math:`{\rm g~cm^{-1}~s^{-1}`           |
| ``angular_momentum_y``,           | computed as :math:`[(\rho \ub) \times {\bf r}]_n` | ``derangmomy``,             |                                         |
| ``angular_momentum_z``            | where :math:`{\bf r}` is the distance from        | ``derangmomz``              |                                         |
|                                   | ``center`` and :math:`n` is either x, y, or z     |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``diff_coeff``                    | Thermal diffusion coefficient,                    | ``derdiffcoeff``            | :math:`{\rm cm^2~s^{-1}}`               |
|                                   | :math:`\kth/(\rho c_v)`                           |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``diff_term``                     | :math:`\nabla\cdot(\kth\nabla T)`                 | ``derdiffterm``             | :math:`{\rm erg~cm^{-3}~s^{-1}}`        |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``divu``                          | :math:`\nabla \cdot \ub`                          | ``derdivu``                 | :math:`{\rm s^{-1}}`                    |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``eint_e``                        | Specific internal energy computed from the        | ``dereint2``                | :math:`{\rm erg~g^{-1}}`                |
|                                   | conserved :math:`(\rho e)` state variable as      |                             |                                         |
|                                   | :math:`e = (\rho e)/\rho`                         |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``eint_E``                        | Specific internal energy computed from the        | ``dereint1``                | :math:`{\rm erg~g^{-1}}`                |
|                                   | total energy and momentum conserved state as      |                             |                                         |
|                                   | :math:`e=[(\rho E)-\frac{1}{2}(\rho \ub^2)]/\rho` |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``entropy``                       | Specific entropy, :math:`s`, computed as          | ``derentropy``              | :math:`{\rm erg~g^{-1}~K^{-1}}`         |
|                                   | :math:`s = s(\rho, e, X_k)`, where `e` is         |                             |                                         |
|                                   | computed from :math:`(\rho e)`                    |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``Ertot``                         | Total radiation energy density                    | ``derertot``                |                                         |
|                                   | (for multigroup radiation problems)               |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``Frcomx``, ``Frcomy``,           | Comoving radiation flux                           | ``Radiation.cpp``           |                                         |
| ``Frcomz``                        |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``Frlabx``, ``Frlaby``,           | Lab-frame radiation flux                          | ``Radiation.cpp``           |                                         |
| ``Frlabz``                        |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``Gamma_1``                       | Adiabatic index,                                  | ``dergamma1``               | --                                      |
|                                   | :math:`d\log p/d\log \rho|_s`                     |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``kineng``                        | Kinetic energy density,                           | ``derkineng``               | :math:`{\rm erg~cm^{-3}}`               |
|                                   | :math:`K = \frac{1}{2} |(\rho \ub)|^2`            |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``lambda``                        | Radiation flux limiter                            |                             | --                                      |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``logden``                        | :math:`\log_{10} \rho`                            | ``derlogten``               | dimensionless, assuming :math:`\rho`    |
|                                   |                                                   |                             | is in CGS                               |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``MachNumber``                    | Fluid Mach number, :math:`|\ub|/c_s`              | ``dermachnumber``           | --                                      |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``maggrav``                       | Gravitational acceleration magnitude              | ``dermaggrav``              | :math:`{\rm cm~s^{-2}}`                 |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magmom``                        | Momentum density magnitude,                       | ``dermagmom``               | :math:`{\rm g~cm^{-2}~s^{-1}}`          |
|                                   | :math:`|\rho \ub|`                                |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magvel``                        | Velocity magnitude, :math:`|\ub|`                 | ``dermagvel``               | :math:`\cms`                            |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magvort``                       | Vorticity magnitude, :math:`|\nabla\times\ub|`    | ``dermagvort``              | :math:`{\rm s^{-1}}`                    |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``pressure``                      | Total pressure, including ions, electrons,        | ``derpres``                 | :math:`{\rm dyn~cm^{-2}}`               |
|                                   | and radiation (for non radhydro problems)         |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``radvel``                        | Radial velocity (measured with respect to         | ``derradialvel``            | :math:`\cms`                            |
|                                   | `center`),                                        |                             |                                         |
|                                   | :math:`(xu + yv + zw)/r`                          |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``soundspeed``                    | Sound speed                                       | ``dersoundspeed``           | :math:`\cms`                            |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``StateErr``                      |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``thermal_cond``                  | Thermal conductivity, :math:`\kth`                | ``dercond``                 | :math:`{\rm erg~cm^{-1}~s^{-1}~K^{-1}}` |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``t_sound_t_enuc``                |                                                   | ``derenuctimescale``        | --                                      |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``uminusc``                       | (only for 1D) x-velocity :math:`-` sound          | ``deruminusc``              | :math:`\cms`                            |
|                                   | speed                                             |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``uplusc``                        | (only for 1D) x-velocity + sound speed            | ``deruplusc``               | :math:`\cms`                            |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``X(q)``                          | Mass fraction of species q                        | ``derspec``                 | --                                      |
|                                   | :math:`X_k = (\rho X_k)/\rho`                     |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``x_velocity``,                   | Fluid velocity,                                   | ``dervel``                  | :math:`\cms`                            |
| ``y_velocity``,                   | :math:`\ub = (\rho \ub)/\rho`                     |                             |                                         |
| ``z_velocity``                    |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+


problem-specific plotfile variables
-----------------------------------

+-----------------------------------+---------------------------------------------------+--------------------------------------+
| variable name                     | description                                       | units                                |
+===================================+===================================================+======================================+
| ``analytic``                      |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``pi``                            |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``pioverp0``                      |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``primarymask``                   |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``secondarymask``                 |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Terror``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Texact``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``inertial_angular_momentum_x``,  |                                                   |                                      |
| ``inertial_angular_momentum_y``,  |                                                   |                                      |
| ``inertial_angular_momentum_z``   |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``inertial_momentum_x``,          |                                                   |                                      |
| ``inertial_momentum_y``,          |                                                   |                                      |
| ``inertial_momentum_z``           |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``inertial_radial_momentum_x``,   |                                                   |                                      |
| ``inertial_radial_momentum_y``,   |                                                   |                                      |
| ``inertial_radial_momentum_z``    |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiEff``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiEffPM_P``                    |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiEffPM_S``                    |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``tpert``                         |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+



Visualization Tools
===================

.. toctree::

   vis_tools


