
#ifdef TRUE_SDC
///
/// Calculate primitive variables from conserved variables (uses StateData)
///
/// @param time     current time
///
    void cons_to_prim(const amrex::Real time);
#endif

///
/// Calculate primitive variables from given conserved variables
///
/// @param u        MultiFab of conserved variables
/// @param q        MultiFab to save primitive variables to
/// @param qaux     MultiFab of auxiliary variables
/// @param time     Current time
///
    void cons_to_prim(amrex::MultiFab& u, amrex::MultiFab& q, amrex::MultiFab& qaux, const amrex::Real time);

///
/// convert the conservative state cell averages to primitive cell
/// averages with 4th order accuracy
///
/// @param time     current time
///
    void cons_to_prim_fourth(const amrex::Real time);

///
/// Check to see if the CFL condition has been violated
///
/// @param State MultiFab of conserved variables
/// @param dt timestep
///
    void check_for_cfl_violation(const amrex::MultiFab& State, const amrex::Real dt);

///
/// this constructs the hydrodynamic source (essentially the flux
/// divergence) using the CTU framework for unsplit hydrodynamics
///
/// @param time     current time
/// @param dt       timestep
///
    void construct_ctu_hydro_source(amrex::Real time, amrex::Real dt);

///
/// this constructs the hydrodynamic source (essentially the flux
/// divergence) using method of lines integration.  The output, is the
/// divergence of the fluxes, A = -div{F(U)}
///
/// @param time      current time
/// @param dt        timestep
/// @param A_update  divergence of the fluxes
///
    void construct_mol_hydro_source(amrex::Real time, amrex::Real dt, amrex::MultiFab& A_update);

///
/// convert the conserved form of the source terms, S(U), into source terms
/// for the primitive variable equations, S(q).
///
/// @param bx      the box to operate over
/// @param q_arr   primitive variable state array
/// @param src     conserved source array
/// @param srcQ    primitive variable source array
///
    void src_to_prim(const amrex::Box& bx,
                     amrex::Array4<amrex::Real const> const& q_arr,
                     amrex::Array4<amrex::Real const> const& src,
                     amrex::Array4<amrex::Real> const& srcQ);

///
/// the actual work routine that does the conversion of conserved to primtiive
/// variables.
///
/// @param bx        the box to operate over
/// @param time      current time
/// @param uin       input conserved state
/// @param Bx        x-component of magentic field (if USE_MHD=TRUE)
/// @param By        y-component of magentic field (if USE_MHD=TRUE)
/// @param Bz        z-component of magentic field (if USE_MHD=TRUE)
/// @param Erin      radiation energy (if USE_RAD=TRUE)
/// @param lam       radiation flux limiter (if USE_RAD=TRUE)
/// @param q_arr     output primitive state
/// @param qaux_arr  output auxillary quantities
///
    void ctoprim(const amrex::Box& bx,
                 const amrex::Real time,
                 amrex::Array4<amrex::Real const> const& uin,
#ifdef MHD
                 amrex::Array4<amrex::Real const> const& Bx,
                 amrex::Array4<amrex::Real const> const& By,
                 amrex::Array4<amrex::Real const> const& Bz,
#endif
#ifdef RADIATION
                 amrex::Array4<amrex::Real const> const& Erin,
                 amrex::Array4<amrex::Real const> const& lam,
#endif
                 amrex::Array4<amrex::Real> const& q_arr,
                 amrex::Array4<amrex::Real> const& qaux_arr);

///
/// compute the flattening coefficient.  This is 0 if we are in a shock and
/// 1 if we are not applying any flattening
///
/// @param bx         the box to operate over
/// @param q_arr      the primitive variable state
/// @param flatn      the flattening coefficient
/// @param pres_comp  index into q_arr of the pressure component
///
    void uflatten(const amrex::Box& bx,
                  amrex::Array4<amrex::Real const> const& q_arr,
                  amrex::Array4<amrex::Real> const& flatn,
                  const int pres_comp);

///
/// A multidimensional shock detection algorithm
///
/// @param bx     the box to operate over
/// @param q_arr  the primitive variable state
/// @param shk    the shock flag (1 = shock, 0 = no shock)
///
    void shock(const amrex::Box& bx,
               amrex::Array4<amrex::Real const> const& q_arr,
               amrex::Array4<amrex::Real> const& shk);

///
/// Compute the node-centered velocity divergence (DU)_{i-1/2.j-1/2,k-1/2}
///
/// @param bx     the box to operate over
/// @param q_arr  input primitive variable state
/// @param div    output velocity divergence
///
    void divu(const amrex::Box& bx,
              amrex::Array4<amrex::Real const> const& q_arr,
              amrex::Array4<amrex::Real> const& div);

///
/// Update the flux with the artifical viscosity.  This is a 2nd order
/// accurate implementation.
///
/// @param bx    the box to operate over
/// @param idir  the coordinate direction (0 = x, 1 = y, 2 = z)
/// @param div   the node-centered velocity divergence
/// @param uin   the conserved state
/// @param flux  the flux in direction idir
///
    void apply_av(const amrex::Box& bx,
                  const int idir,
                  amrex::Array4<amrex::Real const> const& div,
                  amrex::Array4<amrex::Real const> const& uin,
                  amrex::Array4<amrex::Real> const& flux);

#ifdef RADIATION
///
/// Update the radiation flux with the artifical viscosity.  This is a 2nd order
/// accurate implementation.
///
/// @param bx       the box to operate over
/// @param idir     the coordinate direction (0 = x, 1 = y, 2 = z)
/// @param div      the node-centered velocity divergence
/// @param Erin     the radiation energy
/// @param radflux  the radiation energy flux in direction idir
///
    void apply_av_rad(const amrex::Box& bx,
                      const int idir,
                      amrex::Array4<amrex::Real const> const& div,
                      amrex::Array4<amrex::Real const> const& Erin,
                      amrex::Array4<amrex::Real> const& radflux);


    void
    ctu_rad_consup(const amrex::Box& bx,
                   amrex::Array4<amrex::Real> const& update,
                   amrex::Array4<amrex::Real const> const& Erin,
                   amrex::Array4<amrex::Real> const& uout,
                   amrex::Array4<amrex::Real> const& Erout,
                   amrex::Array4<amrex::Real const> const& radflux1,
                   amrex::Array4<amrex::Real const> const& qx,
                   amrex::Array4<amrex::Real const> const& area1,
#if AMREX_SPACEDIM >= 2
                   amrex::Array4<amrex::Real const> const& radflux2,
                   amrex::Array4<amrex::Real const> const& qy,
                   amrex::Array4<amrex::Real const> const& area2,
#endif
#if AMREX_SPACEDIM == 3
                   amrex::Array4<amrex::Real const> const& radflux3,
                   amrex::Array4<amrex::Real const> const& qz,
                   amrex::Array4<amrex::Real const> const& area3,
#endif
                   int& nstep_fsp,
                   amrex::Array4<amrex::Real const> const& vol,
                   const amrex::Real dt);

#endif

///
/// Compute the normal left and right primitive variable interface
/// states at each interface using the piecewise parabolic method with
/// characteristic tracing.  This implementation is for pure
/// hydrodynamics
///
/// @param bx       the box to operate over
/// @param vbx      the valid region box (excludes ghost cells)
/// @param q_arr    the primitive variable state
/// @param flatn    flattening coefficient
/// @param qaux_arr the auxillary state
/// @param srcQ     the source terms for the primitive variable equations
/// @param qxm      left interface state in x, q_{i-1/2,j,k,L}
/// @param qxp      right interface state in x, q_{i-1/2,j,k,R}
/// @param qym      left interface state in y, q_{i,j-1/2,k,L}
/// @param qyp      right interface state in y, q_{i,j-1/2,k,R}
/// @param qzm      left interface state in z, q_{i,j,k-1/2,L}
/// @param qzp      right interface state in z, q_{i,j,k-1/2,R}
/// @param dloga    the geometric factor d(log area)
/// @param dt       timestep
///
    void ctu_ppm_states(const amrex::Box& bx, const amrex::Box& vbx,
                        amrex::Array4<amrex::Real const> const& q_arr,
                        amrex::Array4<amrex::Real const> const& flatn,
                        amrex::Array4<amrex::Real const> const& qaux_arr,
                        amrex::Array4<amrex::Real const> const& srcQ,
                        amrex::Array4<amrex::Real> const& qxm,
                        amrex::Array4<amrex::Real> const& qxp,
#if AMREX_SPACEDIM >= 2
                        amrex::Array4<amrex::Real> const& qym,
                        amrex::Array4<amrex::Real> const& qyp,
#endif
#if AMREX_SPACEDIM == 3
                        amrex::Array4<amrex::Real> const& qzm,
                        amrex::Array4<amrex::Real> const& qzp,
#endif
#if AMREX_SPACEDIM < 3
                        amrex::Array4<amrex::Real const> const& dloga,
#endif
                        const amrex::Real dt);

///
/// Compute the normal left and right primitive variable interface
/// states at each interface using piecewise linear reconstruction with
/// characteristic tracing.  This implementation is for pure
/// hydrodynamics
///
/// @param bx       the box to operate over
/// @param vbx      the valid region box (excludes ghost cells)
/// @param q_arr    the primitive variable state
/// @param flatn    flattening coefficient
/// @param qaux_arr the auxillary state
/// @param srcQ     the source terms for the primitive variable equations
/// @param qxm      left interface state in x, q_{i-1/2,j,k,L}
/// @param qxp      right interface state in x, q_{i-1/2,j,k,R}
/// @param qym      left interface state in y, q_{i,j-1/2,k,L}
/// @param qyp      right interface state in y, q_{i,j-1/2,k,R}
/// @param qzm      left interface state in z, q_{i,j,k-1/2,L}
/// @param qzp      right interface state in z, q_{i,j,k-1/2,R}
/// @param dloga    the geometric factor d(log area)
/// @param dt       timestep
///
    void ctu_plm_states(const amrex::Box& bx, const amrex::Box& vbx,
                        amrex::Array4<amrex::Real const> const& q_arr,
                        amrex::Array4<amrex::Real const> const& flatn,
                        amrex::Array4<amrex::Real const> const& qaux_arr,
                        amrex::Array4<amrex::Real const> const& srcQ,
                        amrex::Array4<amrex::Real> const& qxm,
                        amrex::Array4<amrex::Real> const& qxp,
#if AMREX_SPACEDIM >= 2
                        amrex::Array4<amrex::Real> const& qym,
                        amrex::Array4<amrex::Real> const& qyp,
#endif
#if AMREX_SPACEDIM == 3
                        amrex::Array4<amrex::Real> const& qzm,
                        amrex::Array4<amrex::Real> const& qzp,
#endif
#if AMREX_SPACEDIM < 3
                        amrex::Array4<amrex::Real const> const& dloga,
#endif
                        const amrex::Real dt);

#ifdef RADIATION
///
/// Compute the normal left and right primitive variable interface
/// states at each interface using the piecewise parabolic method with
/// characteristic tracing.  This implementation is for the radiation
/// hydrodynamics.
///
/// @param bx       the box to operate over
/// @param vbx      the valid region box (excludes ghost cells)
/// @param q_arr    the primitive variable state
/// @param flatn    flattening coefficient
/// @param qaux_arr the auxillary state
/// @param srcQ     the source terms for the primitive variable equations
/// @param qxm      left interface state in x, q_{i-1/2,j,k,L}
/// @param qxp      right interface state in x, q_{i-1/2,j,k,R}
/// @param qym      left interface state in y, q_{i,j-1/2,k,L}
/// @param qyp      right interface state in y, q_{i,j-1/2,k,R}
/// @param qzm      left interface state in z, q_{i,j,k-1/2,L}
/// @param qzp      right interface state in z, q_{i,j,k-1/2,R}
/// @param dloga    the geometric factor d(log area)
/// @param dt       timestep
///
    void ctu_ppm_rad_states(const amrex::Box& bx, const amrex::Box& vbx,
                            amrex::Array4<amrex::Real const> const& q_arr,
                            amrex::Array4<amrex::Real const> const& flatn,
                            amrex::Array4<amrex::Real const> const& qaux_arr,
                            amrex::Array4<amrex::Real const> const& srcQ,
                            amrex::Array4<amrex::Real> const& qxm,
                            amrex::Array4<amrex::Real> const& qxp,
#if AMREX_SPACEDIM >= 2
                            amrex::Array4<amrex::Real> const& qym,
                            amrex::Array4<amrex::Real> const& qyp,
#endif
#if AMREX_SPACEDIM == 3
                            amrex::Array4<amrex::Real> const& qzm,
                            amrex::Array4<amrex::Real> const& qzp,
#endif
#if AMREX_SPACEDIM < 3
                            amrex::Array4<amrex::Real const> const& dloga,
#endif
                            const amrex::Real dt);
#endif

///
/// Compute the left and right primitive variable interface states at
/// each interface using piecewise linear reconstruction for a
/// method-of-lines type of integration (i.e., no prediction to n+1/2
/// or characteristic tracing).  This implementation is for pure
/// hydrodynamics
///
/// @param bx          the box to operate over
/// @param idir        coordinate direction (0 = x, 1 = y, 2 = z)
/// @param q_arr       the primitive variable state
/// @param flatn_arr   flattening coefficient
/// @param dq          slope of the primtiive variables
/// @param qm          left interface state, e.g., q_{i-1/2,j,k,L}
/// @param qp          right interface state, e.g., q_{i-1/2,j,k,R}
///
    void mol_plm_reconstruct(const amrex::Box& bx,
                             const int idir,
                             amrex::Array4<amrex::Real const> const& q_arr,
                             amrex::Array4<amrex::Real const> const& flatn_arr,
                             amrex::Array4<amrex::Real const> const& src_q_arr,
                             amrex::Array4<amrex::Real> const& dq,
                             amrex::Array4<amrex::Real> const& qm,
                             amrex::Array4<amrex::Real> const& qp);
///
/// Compute the left and right primitive variable interface states at
/// each interface using piecewise parabolic reconstruction for a
/// method-of-lines type of integration (i.e., no prediction to n+1/2
/// or characteristic tracing).  This implementation is for pure
/// hydrodynamics
///
/// @param bx          the box to operate over
/// @param idir        coordinate direction (0 = x, 1 = y, 2 = z)
/// @param q_arr       the primitive variable state
/// @param flatn_arr   flattening coefficient
/// @param qm          left interface state, e.g., q_{i-1/2,j,k,L}
/// @param qp          right interface state, e.g., q_{i-1/2,j,k,R}
///
    void mol_ppm_reconstruct(const amrex::Box& bx,
                             const int idir,
                             amrex::Array4<amrex::Real const> const& q_arr,
                             amrex::Array4<amrex::Real const> const& flatn_arr,
                             amrex::Array4<amrex::Real> const& qm,
                             amrex::Array4<amrex::Real> const& qp);

///
/// Compute the conservative update of the fluid state given the fluxes.
/// This returns a source term, update, of the form -div{F} + S(U) which
/// can then be added to the old state to update in time.  This version
/// is written for the method-of-lines / true SDC framework.
///
/// @param bx      the box to operate over
/// @param shk     the shock flag
/// @param srcU    source terms for the conservative equations
/// @param update  the final conserved state update, -div{F} + S(U)
/// @param flux0   flux in the x-direction
/// @param flux1   flux in the y-direction
/// @param flux2   flux in the z-direction
/// @param area0   area of x faces
/// @param area1   area of y faces
/// @param area2   area of z faces
/// @param q0      Godunuv state on x faces
/// @param vol     cell volume
///
    void mol_consup(const amrex::Box& bx,
#ifdef SHOCK_VAR
                    amrex::Array4<amrex::Real const> const& shk,
#endif
                    amrex::Array4<amrex::Real const> const& srcU,
                    amrex::Array4<amrex::Real> const& update,
                    const amrex::Real dt,
                    amrex::Array4<amrex::Real const> const& flux0,
#if AMREX_SPACEDIM >= 2
                    amrex::Array4<amrex::Real const> const& flux1,
#endif
#if AMREX_SPACEDIM == 3
                    amrex::Array4<amrex::Real const> const& flux2,
#endif
                    amrex::Array4<amrex::Real const> const& area0,
#if AMREX_SPACEDIM >= 2
                    amrex::Array4<amrex::Real const> const& area1,
#endif
#if AMREX_SPACEDIM == 3
                    amrex::Array4<amrex::Real const> const& area2,
#endif
#if AMREX_SPACEDIM <= 2
                    amrex::Array4<amrex::Real const> const& q0,
#endif
                    amrex::Array4<amrex::Real const> const& vol);

///
/// Compute the diffusive flux term and add it to the energy fluxes.
/// This implementation is second-order accurate.
///
/// @param bx    the box to operate over
/// @param idir  the coordinate direction (0 = x, 1 = y, 2 = z)
/// @param U     the conserved state
/// @param cond  thermal conductivity
/// @param flux  hydrodynamic flux in direction idir
///
     void
     mol_diffusive_flux(const amrex::Box& bx,
                        const int idir,
                        amrex::Array4<amrex::Real const> const& U,
                        amrex::Array4<amrex::Real const> const& cond,
                        amrex::Array4<amrex::Real> const& flux);

///
/// Add the transverse flux difference in idir_t to the interface states in
/// direction idir_n as part of the CTU hydrodynamics update.  This is the main driver.
///
/// @param bx        the box to operate over
/// @param idir_t    direction for the transverse flux difference (0 = x, 1 = y, 2 = z)
/// @param idir_n    direction of the interface states normal (0 = x, 1 = y, 2 = z)
/// @param qm        input left interface state
/// @param qmo       updated left interface state
/// @param qp        input right interface state
/// @param qpo       updated right interface state
/// @param qaux      auxillary state
/// @param flux_t    flux in the idir_t direction
/// @param rflux_t   radiation flux in the idir_t direction
/// @param q_t       Godunov state in the idir_t direction
/// @param area_t    face area in the idir_t direction
/// @param vol       cell volume
/// @param hdt       1/2 * timestep
/// @param cdt       weight * timestep, where the weight comes from the CTU algorithm
///
     void trans_single(const amrex::Box& bx,
                       int idir_t, int idir_n,
                       amrex::Array4<amrex::Real const> const& qm,
                       amrex::Array4<amrex::Real> const& qmo,
                       amrex::Array4<amrex::Real const> const& qp,
                       amrex::Array4<amrex::Real> const& qpo,
                       amrex::Array4<amrex::Real const> const& qaux,
                       amrex::Array4<amrex::Real const> const& flux_t,
#ifdef RADIATION
                       amrex::Array4<amrex::Real const> const& rflux_t,
#endif
                       amrex::Array4<amrex::Real const> const& q_t,
#if AMREX_SPACEDIM == 2
                       amrex::Array4<amrex::Real const> const& area_t,
                       amrex::Array4<amrex::Real const> const& vol,
#endif
                       amrex::Real hdt, amrex::Real cdtdx);

///
/// Add the transverse flux difference in idir_t to the interface states in
/// direction idir_n as part of the CTU hydrodynamics update.
///
/// @param bx        the box to operate over
/// @param idir_t    direction for the transverse flux difference (0 = x, 1 = y, 2 = z)
/// @param idir_n    direction of the interface states normal (0 = x, 1 = y, 2 = z)
/// @param d         a flag set by the calling driver that determines which interface
///                  we are updating
/// @param q_arr     input interface state
/// @param qo_arr    updated interface state
/// @param qaux      auxillary state
/// @param flux_t    flux in the idir_t direction
/// @param rflux_t   radiation flux in the idir_t direction
/// @param q_t       Godunov state in the idir_t direction
/// @param area_t    face area in the idir_t direction
/// @param vol       cell volume
/// @param hdt       1/2 * timestep
/// @param cdt       weight * timestep, where the weight comes from the CTU algorithm
///
     void actual_trans_single(const amrex::Box& bx,
                              int idir_t, int idir_n, int d,
                              amrex::Array4<amrex::Real const> const& q_arr,
                              amrex::Array4<amrex::Real> const& qo_arr,
                              amrex::Array4<amrex::Real const> const& qaux,
                              amrex::Array4<amrex::Real const> const& flux_t,
#ifdef RADIATION
                              amrex::Array4<amrex::Real const> const& rflux_t,
#endif
                              amrex::Array4<amrex::Real const> const& q_t,
#if AMREX_SPACEDIM == 2
                              amrex::Array4<amrex::Real const> const& area_t,
                              amrex::Array4<amrex::Real const> const& vol,
#endif
                              amrex::Real hdt, amrex::Real cdtdx);

///
/// The final transverse update where the flux difference in both
/// perpendicular directions are added to the normal flux.  This is
/// the main driver.
///
/// @param bx        the box to operate over
/// @param idir_n    direction of the interface states normal (0 = x, 1 = y, 2 = z)
/// @param idir_t1   direction for the first transverse flux difference (0 = x, 1 = y, 2 = z)
/// @param idir_t2   direction for the second transverse flux difference (0 = x, 1 = y, 2 = z)
/// @param qm        input left interface state
/// @param qmo       updated left interface state
/// @param qp        input right interface state
/// @param qpo       updated right interface state
/// @param qaux      auxillary state
/// @param flux_t1   flux in the idir_t1 direction
/// @param rflux_t1  radiation flux in the idir_t1 direction
/// @param flux_t2   flux in the idir_t2 direction
/// @param rflux_t2  radiation flux in the idir_t2 direction
/// @param q_t1      Godunov state in the idir_t1 direction
/// @param q_t2      Godunov state in the idir_t2 direction
/// @param hdt       1/2 * timestep
/// @param cdtdx_t1  weight * timestep in idir_t1 direction
/// @param cdtdx_t2  weight * timestep in idir_t2 direction
///
     void trans_final(const amrex::Box& bx,
                      int idir_n, int idir_t1, int idir_t2,
                      amrex::Array4<amrex::Real const> const& qm,
                      amrex::Array4<amrex::Real> const& qmo,
                      amrex::Array4<amrex::Real const> const& qp,
                      amrex::Array4<amrex::Real> const& qpo,
                      amrex::Array4<amrex::Real const> const& qaux,
                      amrex::Array4<amrex::Real const> const& flux_t1,
#ifdef RADIATION
                      amrex::Array4<amrex::Real const> const& rflux_t1,
#endif
                      amrex::Array4<amrex::Real const> const& flux_t2,
#ifdef RADIATION
                      amrex::Array4<amrex::Real const> const& rflux_t2,
#endif
                      amrex::Array4<amrex::Real const> const& q_t1,
                      amrex::Array4<amrex::Real const> const& q_t2,
                      amrex::Real cdtdx_t1, amrex::Real cdtdx_t2);

///
/// The final transverse update where the flux difference in both
/// perpendicular directions are added to the normal flux.
///
/// @param bx        the box to operate over
/// @param idir_n    direction of the interface states normal (0 = x, 1 = y, 2 = z)
/// @param idir_t1   direction for the first transverse flux difference (0 = x, 1 = y, 2 = z)
/// @param idir_t2   direction for the second transverse flux difference (0 = x, 1 = y, 2 = z)
/// @param d         a flag set by the driver to indicate which state we are operating on
/// @param q_arr     input interface state
/// @param qo_arr    updated interface state
/// @param qaux      auxillary state
/// @param flux_t1   flux in the idir_t1 direction
/// @param rflux_t1  radiation flux in the idir_t1 direction
/// @param flux_t2   flux in the idir_t2 direction
/// @param rflux_t2  radiation flux in the idir_t2 direction
/// @param q_t1      Godunov state in the idir_t1 direction
/// @param q_t2      Godunov state in the idir_t2 direction
/// @param hdt       1/2 * timestep
/// @param cdtdx_n   weight * timestep in idir_n direction
/// @param cdtdx_t1  weight * timestep in idir_t1 direction
/// @param cdtdx_t2  weight * timestep in idir_t2 direction
///
     void actual_trans_final(const amrex::Box& bx,
                             int idir_n, int idir_t1, int idir_t2, int d,
                             amrex::Array4<amrex::Real const> const& q_arr,
                             amrex::Array4<amrex::Real> const& qo_arr,
                             amrex::Array4<amrex::Real const> const& qaux,
                             amrex::Array4<amrex::Real const> const& flux_t1,
#ifdef RADIATION
                             amrex::Array4<amrex::Real const> const& rflux_t1,
#endif
                             amrex::Array4<amrex::Real const> const& flux_t2,
#ifdef RADIATION
                             amrex::Array4<amrex::Real const> const& rflux_t2,
#endif
                             amrex::Array4<amrex::Real const> const& q_t1,
                             amrex::Array4<amrex::Real const> const& q_t2,
                             amrex::Real cdtdx_t1, amrex::Real cdtdx_t2);

///
/// Reconstruct the primtive state as parabola, integrate under them,
/// and perform the characteristic tracing to get the interface states.
/// This is for the CTU hydrodynamics scheme.
///
/// @param bx        the box to operate over
/// @param idir      coordinate direction of the interface (0 = x, 1 = y, 2 = z)
/// @param q_arr     primitive variable state
/// @param qaux_arr  the auxillary state
/// @param srcQ      primitive variable equation source terms
/// @param flatn     flattening coefficient
/// @param qm        left interface state (e.g., q_{i-1/2,j,k,L})
/// @param qp        right interface state (e.g., q_{i-1/2,j,k,R})
/// @param dloga     geometric factor dlog(area)
/// @param vbx       the valid region box (excluding ghost cells)
/// @param dt        timestep
///
    void trace_ppm(const amrex::Box& bx,
                   const int idir,
                   amrex::Array4<amrex::Real const> const& q_arr,
                   amrex::Array4<amrex::Real const> const& qaux_arr,
                   amrex::Array4<amrex::Real const> const& srcQ,
                   amrex::Array4<amrex::Real const> const& flatn,
                   amrex::Array4<amrex::Real> const& qm,
                   amrex::Array4<amrex::Real> const& qp,
#if (AMREX_SPACEDIM < 3)
                   amrex::Array4<amrex::Real const> const& dloga,
#endif
                   const amrex::Box& vbx,
                   const amrex::Real dt);

///
/// Reconstruct the primtive state as pieceeise linear, integrate under them,
/// and perform the characteristic tracing to get the interface states.
/// This is for the CTU hydrodynamics scheme.
///
/// @param bx           the box to operate over
/// @param idir         coordinate direction of the interface (0 = x, 1 = y, 2 = z)
/// @param q_arr        primitive variable state
/// @param qaux_arr     the auxillary state
/// @param flatn_arr    flattening coefficient
/// @param qm           left interface state (e.g., q_{i-1/2,j,k,L})
/// @param qp           right interface state (e.g., q_{i-1/2,j,k,R})
/// @param dloga        geometric factor dlog(area)
/// @param srcQ         primitive variable equation source terms
/// @param vbx          the valid region box (excluding ghost cells)
/// @param dt           timestep
///
    void trace_plm(const amrex::Box& bx,
                   const int idir,
                   amrex::Array4<amrex::Real const> const& q_arr,
                   amrex::Array4<amrex::Real const> const& qaux_arr,
                   amrex::Array4<amrex::Real const> const& flatn_arr,
                   amrex::Array4<amrex::Real> const& qm,
                   amrex::Array4<amrex::Real> const& qp,
#if (AMREX_SPACEDIM < 3)
                   amrex::Array4<amrex::Real const> const& dloga,
#endif
                   amrex::Array4<amrex::Real const> const& SrcQ,
                   const amrex::Box& vbx,
                   const amrex::Real dt);

///
/// Reconstruct the primtive state as parabola, integrate under them,
/// and perform the characteristic tracing to get the interface states.
/// This is for the CTU radiation hydrodynamics scheme.
///
/// @param bx        the box to operate over
/// @param idir      coordinate direction of the interface (0 = x, 1 = y, 2 = z)
/// @param q_arr     primitive variable state
/// @param qaux_arr  the auxillary state
/// @param srcQ      primitive variable equation source terms
/// @param flatn     flattening coefficient
/// @param qm        left interface state (e.g., q_{i-1/2,j,k,L})
/// @param qp        right interface state (e.g., q_{i-1/2,j,k,R})
/// @param dloga     geometric factor dlog(area)
/// @param vbx       the valid region box (excluding ghost cells)
/// @param dt        timestep
///
    void trace_ppm_rad(const amrex::Box& bx,
                       const int idir,
                       amrex::Array4<amrex::Real const> const& q_arr,
                       amrex::Array4<amrex::Real const> const& qaux_arr,
                       amrex::Array4<amrex::Real const> const& srcQ,
                       amrex::Array4<amrex::Real const> const& flatn,
                       amrex::Array4<amrex::Real> const& qm,
                       amrex::Array4<amrex::Real> const& qp,
#if (AMREX_SPACEDIM < 3)
                       amrex::Array4<amrex::Real const> const& dloga,
#endif
                       const amrex::Box& vbx,
                       const amrex::Real dt);

///
/// Compute the conservative update of the hydrodynamics state and store it
/// in update so it can be applied later to advance the state
///
/// @param bx      the box to operate over
/// @param shk     the shock flag
/// @param update  the hydrodynamics update term, -div{F} + S(U)
/// @param flux0   flux in the x direction
/// @param qx      Godunov state in the x direction
/// @param area0   area of x faces
/// @param flux1   flux in the y direction
/// @param qy      Godunov state in the y direction
/// @param area1   area of y faces
/// @param flux2   flux in the z direction
/// @param qz      Godunov state in the z direction
/// @param area2   area of z faces
/// @param vol     volume of the cell
/// @param dt      timestep
///
    void consup_hydro(const amrex::Box& bx,
#ifdef SHOCK_VAR
                      amrex::Array4<amrex::Real const> const& shk,
#endif
                      amrex::Array4<amrex::Real> const& update,
                      amrex::Array4<amrex::Real> const& flux0,
                      amrex::Array4<amrex::Real const> const& qx,
                      amrex::Array4<amrex::Real const> const& area0,
#if AMREX_SPACEDIM >= 2
                      amrex::Array4<amrex::Real> const& flux1,
                      amrex::Array4<amrex::Real const> const& qy,
                      amrex::Array4<amrex::Real const> const& area1,
#endif
#if AMREX_SPACEDIM == 3
                      amrex::Array4<amrex::Real> const& flux2,
                      amrex::Array4<amrex::Real const> const& qz,
                      amrex::Array4<amrex::Real const> const& area2,
#endif
                      amrex::Array4<amrex::Real const> const& vol,
                      const amrex::Real dt);


///
/// A general interface to the pure hydrodynamics and radiation
/// Riemann solvers which takes the left and right interface states,
/// qm and qp, and computes the flux through the interface and the
/// godunov state (qint) on the interface
///
/// @param bx          the box to operate over
/// @param qm          left state on the interface
/// @param qp          right state on the interface
/// @param flx         flux through the interface
/// @param qint        tempoary storage of the full interface state
/// @param rflux       radiation fluxes through the interface
/// @param lambda_int  radiation flux limiter on the interface
/// @param qgdnv       Godunov state on the interface
/// @param qaux        auxillary state
/// @param shk         shock flag
/// @param idir        coordinate direction of the solve (0 = x, 1 = y, 2 = z)
///
    void cmpflx_plus_godunov(const amrex::Box& bx,
                             amrex::Array4<amrex::Real> const& qm,
                             amrex::Array4<amrex::Real> const& qp,
                             amrex::Array4<amrex::Real> const& flx,
                             amrex::Array4<amrex::Real> const& qint,
#ifdef RADIATION
                             amrex::Array4<amrex::Real> const& rflx,
                             amrex::Array4<amrex::Real> const& lambda_int,
#endif
                             amrex::Array4<amrex::Real> const& qgdnv,
                             amrex::Array4<amrex::Real const> const& qaux,
                             amrex::Array4<amrex::Real const> const& shk,
                             const int idir);

///
/// Call the hydrodynamic Riemann solvers and return the interface state, without
/// computing the fluxes
///
/// @param bx              the box to operate over
/// @param qm              left state on the interface
/// @param qp              right state on the interface
/// @param qint            the full hydrodynamic interface state
/// @param lambda_int      radiation flux limiter on the interface
/// @param qaux            auxillary state
/// @param idir            coordinate direction of the solve (0 = x, 1 = y, 2 = z)
/// @param compute_gammas  do we call the EOS using the interface states to get Gamma_1 on interfaces
///
    void riemann_state(const amrex::Box& bx,
                       amrex::Array4<amrex::Real> const& qm,
                       amrex::Array4<amrex::Real> const& qp,
                       amrex::Array4<amrex::Real> const& qint,
#ifdef RADIATION
                       amrex::Array4<amrex::Real> const& lambda_int,
#endif
                       amrex::Array4<amrex::Real const> const& qaux,
                       const int idir, const int compute_gammas);

///
/// The Colella-Glaz Riemann solver for pure hydrodynamics.  This is a
/// two shock approximate state Riemann solver.
///
/// @param bx         the box to operate over
/// @param ql         the left interface state
/// @param qr         the right interface state
/// @param qaux_arr   the auxillary state
/// @param qint       the full Godunov state on the interface
/// @param idir       coordinate direction for the solve (0 = x, 1 = y, 2 = z)
///
    void riemanncg(const amrex::Box& bx,
                   amrex::Array4<amrex::Real> const& ql,
                   amrex::Array4<amrex::Real> const& qr,
                   amrex::Array4<amrex::Real const> const& qaux_arr,
                   amrex::Array4<amrex::Real> const& qint,
                   const int idir);

///
/// The Colella-Glaz-Ferguson Riemann solver for hydrodynamics and
/// radiation hydrodynamics.  This is a two shock approximate state
/// Riemann solver.
///
/// @param bx          the box to operate over
/// @param ql          the left interface state
/// @param qr          the right interface state
/// @param qaux_arr    the auxillary state
/// @param qint        the full Godunov state on the interface
/// @param lambda_int  the radiation flux limiter on the interface
/// @param idir        coordinate direction for the solve (0 = x, 1 = y, 2 = z)
///
    void riemannus(const amrex::Box& bx,
                   amrex::Array4<amrex::Real> const& ql,
                   amrex::Array4<amrex::Real> const& qr,
                   amrex::Array4<amrex::Real const> const& qaux_arr,
                   amrex::Array4<amrex::Real> const& qint,
#ifdef RADIATION
                   amrex::Array4<amrex::Real> const& lambda_int,
#endif
                   const int idir, const int compute_gammas);

///
/// A HLLC Riemann solver for pure hydrodynamics
///
/// @param bx          the box to operate over
/// @param ql          the left interface state
/// @param qr          the right interface state
/// @param qaux_arr    the auxillary state
/// @param uflx        the flux through the interface
/// @param qint        an approximate Godunov state on the interface
/// @param idir        coordinate direction for the solve (0 = x, 1 = y, 2 = z)
///
    void HLLC(const amrex::Box& bx,
              amrex::Array4<amrex::Real const> const& ql,
              amrex::Array4<amrex::Real const> const& qr,
              amrex::Array4<amrex::Real const> const& qaux_arr,
              amrex::Array4<amrex::Real> const& uflx,
              amrex::Array4<amrex::Real> const& qint,
              const int idir);

    void
    compute_flux_q(const amrex::Box& bx,
                   amrex::Array4<amrex::Real const> const& qint,
                   amrex::Array4<amrex::Real> const& F,
#ifdef RADIATION
                   amrex::Array4<amrex::Real const> const& lambda,
                   amrex::Array4<amrex::Real> const& rF,
#endif
                   const int idir);

    void
    store_godunov_state(const amrex::Box& bx,
                        amrex::Array4<amrex::Real const> const& qint,
#ifdef RADIATION
                        amrex::Array4<amrex::Real const> const& lambda,
#endif
                        amrex::Array4<amrex::Real> const& qgdnv);

    void normalize_species_fluxes(const amrex::Box& bx, amrex::Array4<amrex::Real> const& flux);

    void
    limit_hydro_fluxes_on_small_dens(const amrex::Box& bx,
                                     int idir,
                                     amrex::Array4<amrex::Real const> const& u,
                                     amrex::Array4<amrex::Real const> const& q,
                                     amrex::Array4<amrex::Real const> const& vol,
                                     amrex::Array4<amrex::Real> const& flux,
                                     amrex::Array4<amrex::Real const> const& area,
                                     amrex::Real dt);

    void
    limit_hydro_fluxes_on_large_vel(const amrex::Box& bx,
                                    int idir,
                                    amrex::Array4<amrex::Real const> const& u,
                                    amrex::Array4<amrex::Real const> const& q,
                                    amrex::Array4<amrex::Real const> const& vol,
                                    amrex::Array4<amrex::Real> const& flux,
                                    amrex::Array4<amrex::Real const> const& area,
                                    amrex::Real dt);

    void scale_flux(const amrex::Box& bx,
#if AMREX_SPACEDIM == 1
                    amrex::Array4<amrex::Real const> const& qint,
#endif
                    amrex::Array4<amrex::Real> const& flux,
                    amrex::Array4<amrex::Real const> const& area,
                    const amrex::Real dt);

#ifdef RADIATION
    void scale_rad_flux(const amrex::Box& bx,
                        amrex::Array4<amrex::Real> const& rflux,
                        amrex::Array4<amrex::Real const> const& area,
                        const amrex::Real dt);
#endif

    void reset_edge_state_thermo(const amrex::Box& bx,
                                 amrex::Array4<amrex::Real> const& qedge);

    void edge_state_temp_to_pres(const amrex::Box& bx,
                                 amrex::Array4<amrex::Real> const& qm,
                                 amrex::Array4<amrex::Real> const& qp);


    void do_enforce_minimum_density(const amrex::Box& bx,
                                    amrex::Array4<amrex::Real> const& state_arr,
                                    const int verbose);

#ifdef HYBRID_MOMENTUM

///
/// Construct hybrid source terms at old time
///
/// @param source   MultiFab to save source to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_hybrid_source(amrex::MultiFab& source, amrex::MultiFab& state, amrex::Real time, amrex::Real dt);


///
/// Construct hybrid source terms at new time
///
/// @param source       MultiFab to save source to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_hybrid_source(amrex::MultiFab& source, amrex::MultiFab& state_old, amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);


///
/// Fill ``source`` with hybrid source terms
///
/// @param state        Current state
/// @param source       MultiFab to save source to
/// @param mult_factor
///
    void fill_hybrid_hydro_source(amrex::MultiFab& state, amrex::MultiFab& source, const amrex::Real mult_factor);


///
/// Synchronize linear momentum with hybrid momentum
///
/// @param state    Current state
/// @param ng       Number of ghost cells
///
    void hybrid_to_linear_momentum(amrex::MultiFab& state, int ng = 0);


///
/// Synchronize hybrid momentum with linear momentum
///
/// @param state    Current state
/// @param ng       Number of ghost cells
///
    void linear_to_hybrid_momentum(amrex::MultiFab& state, int ng = 0);

#endif

#ifdef TRUE_SDC

    void
    fourth_interfaces(const amrex::Box& bx,
                      const int idir, const int ncomp,
                      amrex::Array4<amrex::Real const> const& a,
                      amrex::Array4<amrex::Real> const& a_int);

    void
    states(const amrex::Box& bx,
           const int idir, const int ncomp,
           amrex::Array4<amrex::Real const> const& a,
           amrex::Array4<amrex::Real const> const& a_int,
           amrex::Array4<amrex::Real const> const& flatn,
           amrex::Array4<amrex::Real> const& al,
           amrex::Array4<amrex::Real> const& ar);

    void
    fourth_avisc(const amrex::Box& bx,
                 amrex::Array4<amrex::Real const> const& q,
                 amrex::Array4<amrex::Real const> const& qaux,
                 amrex::Array4<amrex::Real> const& avis,
                 const int idir);

    void
    fourth_add_diffusive_flux(const amrex::Box& bx,
                              amrex::Array4<amrex::Real const> const& q,
                              const int temp_comp,
                              amrex::Array4<amrex::Real const> const& qint,
                              amrex::Array4<amrex::Real> const& F,
                              const int idir, const bool is_avg);

    void make_cell_center(const Box& bx,
                          Array4<Real const> const& U,
                          Array4<Real> const& U_cc,
                          GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi);

    void make_cell_center_in_place(const Box& bx,
                                   Array4<Real> const& U,
                                   Array4<Real> const& tmp,
                                   GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi);

    void compute_lap_term(const Box& bx,
                          Array4<Real const> const& U,
                          Array4<Real> const& lap, const int ncomp,
                          GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi);

    void make_fourth_average(const Box& bx,
                             Array4<Real> const& q,
                             Array4<Real const> const& q_bar,
                             GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi);

    void make_fourth_in_place(const Box& bx,
                              Array4<Real> const& q,
                              Array4<Real> const& tmp,
                              GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi);

    void make_fourth_in_place_n(const Box& bx,
                                Array4<Real> const& q, const int ncomp,
                                Array4<Real> const& tmp,
                                GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi);

#endif
