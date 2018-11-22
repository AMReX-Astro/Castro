namespace: ``castro``
---------------------

**diffusion**

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``diffuse_temp``                       | enable thermal diffusion                                | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``diffuse_enth``                       | enable enthalpy diffusion                               | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``diffuse_spec``                       | enable species diffusion                                | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``diffuse_vel``                        | enable velocity diffusion                               | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``diffuse_cutoff_density``             | set a cutoff density for diffusion -- we zero the term  | -1.e200       |
|                                        | out below this density                                  |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``diffuse_cond_scale_fac``             | scaling factor for conductivity                         | 1.0           |
+----------------------------------------+---------------------------------------------------------+---------------+



**particles**

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``do_tracer_particles``                | permits tracer particle calculation to be turned on and | 0             |
|                                        | off                                                     |               |
+----------------------------------------+---------------------------------------------------------+---------------+



**embiggening**

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``grown_factor``                       | the factor by which to extend the domain upon restart   | 1             |
|                                        | for embiggening                                         |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``star_at_center``                     | used with the embiggening routines to determine how to  | -1            |
|                                        | extend the domain                                       |               |
+----------------------------------------+---------------------------------------------------------+---------------+



**reactions**

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``dtnuc_e``                            | Limit the timestep based on how much the burning can    | 1.e200        |
|                                        | change the internal energy of a zone. The timestep is   |               |
|                                        | equal to {\tt dtnuc}  $\cdot\,(e / \dot{e})$.           |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``dtnuc_X``                            | Limit the timestep based on how much the burning can    | 1.e200        |
|                                        | change the species mass fractions of a zone. The        |               |
|                                        | timestep is equal to {\tt dtnuc}  $\cdot\,(X /          |               |
|                                        | \dot{X})$.                                              |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``dtnuc_X_threshold``                  | If we are using the timestep limiter based on changes   | 1.e-3         |
|                                        | in $X$, set a threshold on the species abundance below  |               |
|                                        | which the limiter is not applied. This helps prevent    |               |
|                                        | the timestep from becoming very small due to changes in |               |
|                                        | trace species.                                          |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``dxnuc``                              | limit the zone size based on how much the burning can   | 1.e200        |
|                                        | change the internal energy of a zone. The zone size on  |               |
|                                        | the finest level must be smaller than {\tt dxnuc}       |               |
|                                        | $\cdot\, c_s\cdot (e / \dot{e})$, where $c_s$ is the    |               |
|                                        | sound speed. This ensures that the sound-crossing time  |               |
|                                        | is smaller than the nuclear energy injection timescale. |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``dxnuc_max``                          | Disable limiting based on dxnuc above this threshold.   | 1.e200        |
|                                        | This allows zones that have already ignited or are      |               |
|                                        | about to ignite to be de-refined.                       |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``max_dxnuc_lev``                      | Disable limiting based on dxnuc above this AMR level.   | -1            |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``do_react``                           | permits reactions to be turned on and off -- mostly for | -1            |
|                                        | efficiency's sake                                       |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``react_T_min``                        | minimum temperature for allowing reactions to occur in  | 0.0           |
|                                        | a zone                                                  |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``react_T_max``                        | maximum temperature for allowing reactions to occur in  | 1.e200        |
|                                        | a zone                                                  |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``react_rho_min``                      | minimum density for allowing reactions to occur in a    | 0.0           |
|                                        | zone                                                    |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``react_rho_max``                      | maximum density for allowing reactions to occur in a    | 1.e200        |
|                                        | zone                                                    |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``disable_shock_burning``              | disable burning inside hydrodynamic shock regions       | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+



**diagnostics, I/O**

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``print_fortran_warnings``             | display warnings in Fortran90 routines                  | (0, 1)        |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``print_update_diagnostics``           | display information about updates to the state (how     | (0, 1)        |
|                                        | much mass, momentum, energy added)                      |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``track_grid_losses``                  | calculate losses of material through physical grid      | 0             |
|                                        | boundaries                                              |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``sum_interval``                       | how often (number of coarse timesteps) to compute       | -1            |
|                                        | integral sums (for runtime diagnostics)                 |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``sum_per``                            | how often (simulation time) to compute integral sums    | -1.0e0        |
|                                        | (for runtime diagnostics)                               |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``show_center_of_mass``                | display center of mass diagnostics                      | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``hard_cfl_limit``                     | abort if we exceed CFL = 1 over the cource of a         | 1             |
|                                        | timestep                                                |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``job_name``                           | a string describing the simulation that will be copied  | ""            |
|                                        | into the plotfile's {\tt job\_info} file                |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``output_at_completion``               | write a final plotfile and checkpoint upon completion   | 1             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``reset_checkpoint_time``              | Do we want to reset the time in the checkpoint? This    | -1.e200       |
|                                        | ONLY takes effect if amr.regrid\_on\_restart = 1 and    |               |
|                                        | amr.checkpoint\_on\_restart = 1, (which require that    |               |
|                                        | max\_step and stop\_time be less than the value in the  |               |
|                                        | checkpoint) and you set it to value greater than this   |               |
|                                        | default value.                                          |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``reset_checkpoint_step``              | Do we want to reset the number of steps in the          | -1            |
|                                        | checkpoint? This ONLY takes effect if                   |               |
|                                        | amr.regrid\_on\_restart = 1 and                         |               |
|                                        | amr.checkpoint\_on\_restart = 1, (which require that    |               |
|                                        | max\_step and stop\_time be less than the value in the  |               |
|                                        | checkpoint) and you set it to value greater than this   |               |
|                                        | default value.                                          |               |
+----------------------------------------+---------------------------------------------------------+---------------+



**gravity and rotation**

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``do_grav``                            | permits gravity calculation to be turned on and off     | -1            |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``moving_center``                      | to we recompute the center used for the multipole       | 0             |
|                                        | gravity solve each step?                                |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``grav_source_type``                   | determines how the gravitational source term is added   | 4             |
|                                        | to the momentum and energy state variables.             |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``do_rotation``                        | permits rotation calculation to be turned on and off    | -1            |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``rotational_period``                  | the rotation period for the corotating frame            | -1.e200       |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``rotational_dPdt``                    | the rotation periods time evolution---this allows the   | 0.0           |
|                                        | rotation rate to change durning the simulation time     |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``rotation_include_centrifugal``       | permits the centrifugal terms in the rotation to be     | 1             |
|                                        | turned on and off                                       |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``rotation_include_coriolis``          | permits the Coriolis terms in the rotation to be turned | 1             |
|                                        | on and off                                              |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``rotation_include_domegadt``          | permits the d(omega)/dt terms in the rotation to be     | 1             |
|                                        | turned on and off                                       |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``state_in_rotating_frame``            | Which reference frame to measure the state variables    | 1             |
|                                        | with respect to. The standard in the literature when    |               |
|                                        | using a rotating reference frame is to measure the      |               |
|                                        | state variables with respect to an observer fixed in    |               |
|                                        | that rotating frame. If this option is disabled by      |               |
|                                        | setting it to 0, the state variables will be measured   |               |
|                                        | with respect to an observer fixed in the inertial frame |               |
|                                        | (but the frame will still rotate).                      |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``rot_source_type``                    | determines how the rotation source terms are added to   | 4             |
|                                        | the momentum and energy equations                       |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``implicit_rotation_update``           | we can do a implicit solution of the rotation update to | 1             |
|                                        | allow for better coupling of the Coriolis terms         |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``rot_axis``                           | the coordinate axis ($x=1$, $y=2$, $z=3$) for the       | 3             |
|                                        | rotation vector                                         |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``use_point_mass``                     | include a central point mass                            | 1             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``point_mass``                         | mass of the point mass                                  | 0.0           |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``point_mass_fix_solution``            | if we have a central point mass, we can prevent mass    | 0             |
|                                        | from building up in the zones adjacent to it by keeping |               |
|                                        | their density constant and adding their mass to the     |               |
|                                        | point mass object                                       |               |
+----------------------------------------+---------------------------------------------------------+---------------+



**AMR**

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``state_interp_order``                 | highest order used in interpolation                     | 1             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``lin_limit_state_interp``             | how to do limiting of the state data when interpolating | 0             |
|                                        | 0: only prevent new extrema 1: preserve linear          |               |
|                                        | combinations of state variables                         |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``state_nghost``                       | Number of ghost zones for state data to have. Note that | 0             |
|                                        | if you are using radiation, choosing this to be zero    |               |
|                                        | will be overridden since radiation needs at least one   |               |
|                                        | ghost zone.                                             |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``do_reflux``                          | do we do the hyperbolic reflux at coarse-fine           | 1             |
|                                        | interfaces?                                             |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``update_sources_after_reflux``        | whether to re-compute new-time source terms after a     | 1             |
|                                        | reflux                                                  |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``use_custom_knapsack_weights``        | should we have state data for custom load-balancing     | 0             |
|                                        | weighting?                                              |               |
+----------------------------------------+---------------------------------------------------------+---------------+



**refinement**

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``do_special_tagging``                 |                                                         | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``spherical_star``                     |                                                         | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+



**timestep control**

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``fixed_dt``                           | a fixed timestep to use for all steps (negative turns   | -1.0          |
|                                        | it off)                                                 |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``initial_dt``                         | the initial timestep (negative uses the step returned   | -1.0          |
|                                        | from the timestep constraints)                          |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``dt_cutoff``                          | the smallest valid timestep---if we go below this, we   | 0.0           |
|                                        | abort                                                   |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``max_dt``                             | the largest valid timestep---limit all timesteps to be  | 1.e200        |
|                                        | no larger than this                                     |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``cfl``                                | the effective Courant number to use---we will not allow | 0.8           |
|                                        | the hydrodynamic waves to cross more than this fraction |               |
|                                        | of a zone over a single timestep                        |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``init_shrink``                        | a factor by which to reduce the first timestep from     | 1.0           |
|                                        | that requested by the timestep estimators               |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``change_max``                         | the maximum factor by which the timestep can increase   | 1.1           |
|                                        | from one step to the next.                              |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``plot_per_is_exact``                  | enforce that the AMR plot interval must be hit exactly  | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``small_plot_per_is_exact``            | enforce that the AMR small plot interval must be hit    | 0             |
|                                        | exactly                                                 |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``use_retry``                          | Retry a timestep if it violated the timestep-limiting   | 0             |
|                                        | criteria over the course of an advance. The criteria    |               |
|                                        | will suggest a new timestep that satisfies the          |               |
|                                        | criteria, and we will do subcycled timesteps on the     |               |
|                                        | same level until we reach the original target time.     |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``retry_tolerance``                    | Tolerance to use when evaluating whether to do a retry. | 0.02          |
|                                        | The timestep suggested by the retry will be multiplied  |               |
|                                        | by (1 + this factor) before comparing the actual        |               |
|                                        | timestep to it. If set to some number slightly larger   |               |
|                                        | than zero, then this prevents retries that are caused   |               |
|                                        | by small numerical differences.                         |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``retry_neg_dens_factor``              | If we're doing retries, set the target threshold for    | 1.e-1         |
|                                        | changes in density if a retry is triggered by a         |               |
|                                        | negative density. If this is set to a negative number   |               |
|                                        | then it will disable retries using this criterion.      |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``retry_subcycle_factor``              | When performing a retry, the factor to multiply the     | 0.5           |
|                                        | current timestep by when trying again.                  |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``use_post_step_regrid``               | Check for a possible post-timestep regrid if certain    | 0             |
|                                        | stability criteria were violated.                       |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``max_subcycles``                      | Do not permit more subcycled timesteps than this        | 10            |
|                                        | parameter. Set to a negative value to disable this      |               |
|                                        | criterion.                                              |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``clamp_subcycles``                    | If we do request more than the maximum number of        | 1             |
|                                        | subcycles, should we fail, or should we clamp to that   |               |
|                                        | maximum number and perform that many?                   |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``sdc_iters``                          | Number of iterations for the SDC advance.               | 2             |
+----------------------------------------+---------------------------------------------------------+---------------+



**hydrodynamics**

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``difmag``                             | the coefficient of the artificial viscosity             | 0.1           |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``small_dens``                         | the small density cutoff.  Densities below this value   | -1.e200       |
|                                        | will be reset                                           |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``small_temp``                         | the small temperature cutoff.  Temperatures below this  | -1.e200       |
|                                        | value will be reset                                     |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``small_pres``                         | the small pressure cutoff.  Pressures below this value  | -1.e200       |
|                                        | will be reset                                           |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``small_ener``                         | the small specific internal energy cutoff.  Internal    | -1.e200       |
|                                        | energies below this value will be reset                 |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``do_hydro``                           | permits hydro to be turned on and off for running pure  | -1            |
|                                        | rad problems                                            |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``do_ctu``                             | do we do the CTU unsplit method or a method-of-lines    | 1             |
|                                        | approach?                                               |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``fourth_order``                       | do we do fourth-order accurate MOL hydro?               | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``add_ext_src``                        | if true, define an additional source term               | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``hybrid_hydro``                       | whether to use the hybrid advection scheme that updates | 0             |
|                                        | z-angular momentum, cylindrical momentum, and azimuthal |               |
|                                        | momentum (3D only)                                      |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``ppm_type``                           | reconstruction type: 0: piecewise linear; 1: classic    | 1             |
|                                        | Colella \& Woodward ppm; 2: extrema-preserving ppm      |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``ppm_temp_fix``                       | various methods of giving temperature a larger role in  | 0             |
|                                        | the reconstruction---see Zingale \& Katz 2015           |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``ppm_predict_gammae``                 | do we construct $\gamma_e = p/(\rho e) + 1$ and bring   | 0             |
|                                        | it to the interfaces for additional thermodynamic       |               |
|                                        | information (this is the Colella \& Glaz technique) or  |               |
|                                        | do we use $(\rho e)$ (the classic \castro\ behavior).   |               |
|                                        | Note this also uses $\tau = 1/\rho$ instead of $\rho$.  |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``ppm_reference_eigenvectors``         | do we use the reference state in evaluating the         | 0             |
|                                        | eigenvectors?                                           |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``plm_iorder``                         | for piecewise linear, reconstruction order to use       | 2             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``hybrid_riemann``                     | do we drop from our regular Riemann solver to HLL when  | 0             |
|                                        | we are in shocks to avoid the odd-even decoupling       |               |
|                                        | instability?                                            |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``riemann_solver``                     | which Riemann solver do we use: 0: Colella, Glaz, \&    | 0             |
|                                        | Ferguson (a two-shock solver); 1: Colella \& Glaz (a    |               |
|                                        | two-shock solver) 2: HLLC                               |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``cg_maxiter``                         | for the Colella \& Glaz Riemann solver, the maximum     | 12            |
|                                        | number of iterations to take when solving for the star  |               |
|                                        | state                                                   |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``cg_tol``                             | for the Colella \& Glaz Riemann solver, the tolerance   | 1.0e-5        |
|                                        | to demand in finding the star state                     |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``cg_blend``                           | for the Colella \& Glaz Riemann solver, what to do if   | 2             |
|                                        | we do not converge to a solution for the star state. 0  |               |
|                                        | = do nothing; print iterations and exit 1 = revert to   |               |
|                                        | the original guess for p-star 2 = do a bisection search |               |
|                                        | for another 2 * cg\_maxiter iterations.                 |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``use_eos_in_riemann``                 | should we use the EOS in the Riemann solver to ensure   | 0             |
|                                        | thermodynamic consistency?                              |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``use_flattening``                     | flatten the reconstructed profiles around shocks to     | 1             |
|                                        | prevent them from becoming too thin                     |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``transverse_use_eos``                 | after we add the transverse correction to the interface | 0             |
|                                        | states, replace the predicted pressure with an EOS call |               |
|                                        | (using $e$ and $\rho$).                                 |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``transverse_reset_density``           | if the transverse interface state correction, if the    | 1             |
|                                        | new density is negative, then replace all of the        |               |
|                                        | interface quantities with their values without the      |               |
|                                        | transverse correction.                                  |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``transverse_reset_rhoe``              | if the interface state for $(\rho e)$ is negative after | 0             |
|                                        | we add the transverse terms, then replace the interface |               |
|                                        | value of $(\rho e)$ with a value constructed from the   |               |
|                                        | $(\rho e)$ evolution equation                           |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``dual_energy_eta1``                   | Threshold value of (E - K) / E such that above eta1,    | 1.0e0         |
|                                        | the hydrodynamic pressure is derived from E - K;        |               |
|                                        | otherwise, we use the internal energy variable UEINT.   |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``dual_energy_eta2``                   | Threshold value of (E - K) / E such that above eta2, we | 1.0e-4        |
|                                        | update the internal energy variable UEINT to match E -  |               |
|                                        | K. Below this, UEINT remains unchanged.                 |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``use_pslope``                         | for the piecewise linear reconstruction, do we subtract | 1             |
|                                        | off $(\rho g)$ from the pressure before limiting?       |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``fix_mass_flux``                      |                                                         | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``limit_fluxes_on_small_dens``         | Should we limit the density fluxes so that we do not    | 0             |
|                                        | create small densities?                                 |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``density_reset_method``               | Which method to use when resetting a negative/small     | 1             |
|                                        | density 1 = Reset to characteristics of adjacent zone   |               |
|                                        | with largest density 2 = Use average of all adjacent    |               |
|                                        | zones for all state variables 3 = Reset to the original |               |
|                                        | zone state before the hydro update                      |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``allow_small_energy``                 | Whether or not to allow the internal energy to be less  | 1             |
|                                        | than the internal energy corresponding to small\_temp   |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``do_sponge``                          | permits sponge to be turned on and off                  | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``sponge_implicit``                    | if we are using the sponge, whether to use the implicit | 1             |
|                                        | solve for it                                            |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``source_term_predictor``              | extrapolate the source terms (gravity and rotation) to  | 0             |
|                                        | $n+1/2$ timelevel for use in the interface state        |               |
|                                        | prediction                                              |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``first_order_hydro``                  | set the flattening parameter to zero to force the       | 0             |
|                                        | reconstructed profiles to be flat, resulting in a       |               |
|                                        | first-order method                                      |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``xl_ext_bc_type``                     | if we are doing an external -x boundary condition, who  | ""            |
|                                        | do we interpret it?                                     |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``xr_ext_bc_type``                     | if we are doing an external +x boundary condition, who  | ""            |
|                                        | do we interpret it?                                     |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``yl_ext_bc_type``                     | if we are doing an external -y boundary condition, who  | ""            |
|                                        | do we interpret it?                                     |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``yr_ext_bc_type``                     | if we are doing an external +y boundary condition, who  | ""            |
|                                        | do we interpret it?                                     |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``zl_ext_bc_type``                     | if we are doing an external -z boundary condition, who  | ""            |
|                                        | do we interpret it?                                     |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``zr_ext_bc_type``                     | if we are doing an external +z boundary condition, who  | ""            |
|                                        | do we interpret it?                                     |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``hse_zero_vels``                      | if we are doing HSE boundary conditions, do we zero the | 0             |
|                                        | velocity?                                               |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``hse_interp_temp``                    | if we are doing HSE boundary conditions, should we get  | 0             |
|                                        | the temperature via interpolation (using model\_parser) |               |
|                                        | or hold it constant?                                    |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``hse_reflect_vels``                   | if we are doing HSE boundary conditions, how do we      | 0             |
|                                        | treat the velocity? reflect? or outflow?                |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``mol_order``                          | integration order for MOL integration 1 = first order,  | 2             |
|                                        | 2 = second order TVD, 3 = 3rd order TVD, 4 = 4th order  |               |
|                                        | RK                                                      |               |
+----------------------------------------+---------------------------------------------------------+---------------+



**parallelization**

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``do_acc``                             | determines whether we use accelerators for specific     | -1            |
|                                        | loops                                                   |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``bndry_func_thread_safe``             |                                                         | 1             |
+----------------------------------------+---------------------------------------------------------+---------------+



namespace: ``diffusion``
------------------------

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``v``                                  | the level of verbosity for the diffusion solve (higher  | 0             |
|                                        | number means more output)                               |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``mlmg_maxorder``                      | Use MLMG as the operator                                | 4             |
+----------------------------------------+---------------------------------------------------------+---------------+



namespace: ``gravity``
----------------------

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``gravity_type``                       | what type                                               | "fillme"      |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``const_grav``                         | if doing constant gravity, what is the acceleration     | 0.0           |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``direct_sum_bcs``                     | Check if the user wants to compute the boundary         | 0             |
|                                        | conditions using the brute force method.  Default is    |               |
|                                        | false, since this method is slow.                       |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``drdxfac``                            | ratio of dr for monopole gravity binning to grid        | 1             |
|                                        | resolution                                              |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``max_multipole_order``                | the maximum mulitpole order to use for multipole BCs    | 0             |
|                                        | when doing Poisson gravity                              |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``v``                                  | the level of verbosity for the gravity solve (higher    | 0             |
|                                        | number means more output on the status of the solve /   |               |
|                                        | multigrid                                               |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``no_sync``                            | do we perform the synchronization at coarse-fine        | 0             |
|                                        | interfaces?                                             |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``no_composite``                       | do we do a composite solve?                             | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``do_composite_phi_correction``        | should we apply a lagged correction to the potential    | 1             |
|                                        | that gets us closer to the composite solution? This     |               |
|                                        | makes the resulting fine grid calculation slightly more |               |
|                                        | accurate, at the cost of an additional Poisson solve    |               |
|                                        | per timestep.                                           |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``max_solve_level``                    | For all gravity types, we can choose a maximum level    | MAX\_LEV-1    |
|                                        | for explicitly calculating the gravity and associated   |               |
|                                        | potential. Above that level, we interpolate from        |               |
|                                        | coarser levels.                                         |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``get_g_from_phi``                     | For non-Poisson gravity, do we want to construct the    | 0             |
|                                        | gravitational acceleration by taking the gradient of    |               |
|                                        | the potential, rather than constructing it directly?    |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``mlmg_max_fmg_iter``                  | how many FMG cycles?                                    | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``mlmg_agglomeration``                 | Do agglomeration?                                       | 1             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``mlmg_consolidation``                 |                                                         | 1             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``mlmg_nsolve``                        | Do N-Solve?                                             | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+



namespace: ``particles``
------------------------

+----------------------------------------+---------------------------------------------------------+---------------+
| parameter                              | description                                             | default value |
+========================================+=========================================================+===============+
| ``v``                                  | the level of verbosity for the tracer particle (0 or 1) | 0             |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``particle_init_file``                 | the name of an input file containing the total particle | ""            |
|                                        | number and the initial position of each particle.       |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``particle_restart_file``              | the name of a file with new particles at restart        | ""            |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``restart_from_nonparticle_chkfile``   | to restart from a checkpoint that was written with {\tt | 0             |
|                                        | USE\_PARTICLES}=FALSE                                   |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``particle_output_file``               | the name of timestamp files.                            | ""            |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``timestamp_dir``                      | the name of a directory in which timestamp files are    | ""            |
|                                        | stored.                                                 |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``timestamp_density``                  | whether the local densities at given positions of       | 1             |
|                                        | particles are stored in output files                    |               |
+----------------------------------------+---------------------------------------------------------+---------------+
| ``timestamp_temperature``              | whether the local temperatures at given positions of    | 0             |
|                                        | particles are stored in output files                    |               |
+----------------------------------------+---------------------------------------------------------+---------------+



