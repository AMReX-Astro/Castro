
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
/// @param time     current time
/// @param dt   timestep
/// @param A_update   divergence of the fluxes
///
    void construct_mol_hydro_source(amrex::Real time, amrex::Real dt, amrex::MultiFab& A_update);

    void src_to_prim(const amrex::Box& bx,
                     amrex::Array4<amrex::Real const> const& q,
                     amrex::Array4<amrex::Real const> const& src,
                     amrex::Array4<amrex::Real> const& srcQ);

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

    void uflatten(const amrex::Box& bx,
                  amrex::Array4<amrex::Real const> const& q,
                  amrex::Array4<amrex::Real> const& flatn,
                  const int pres_comp);

    void shock(const amrex::Box& bx,
               amrex::Array4<amrex::Real const> const& q,
               amrex::Array4<amrex::Real> const& shk);

    void divu(const amrex::Box& bx,
              amrex::Array4<amrex::Real const> const& q,
              amrex::Array4<amrex::Real> const& div);

    void apply_av(const amrex::Box& bx,
                  const int idir,
                  amrex::Array4<amrex::Real const> const& div,
                  amrex::Array4<amrex::Real const> const& uin,
                  amrex::Array4<amrex::Real> const& flux);

#ifdef RADIATION
    void apply_av_rad(const amrex::Box& bx,
                      const int idir,
                      amrex::Array4<amrex::Real const> const& div,
                      amrex::Array4<amrex::Real const> const& Erin,
                      amrex::Array4<amrex::Real> const& radflux);
#endif

    void ctu_ppm_states(const amrex::Box& bx, const amrex::Box& vbx,
                        amrex::Array4<amrex::Real const> const& q,
                        amrex::Array4<amrex::Real const> const& flatn,
                        amrex::Array4<amrex::Real const> const& qaux,
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
    void ctu_ppm_rad_states(const amrex::Box& bx, const amrex::Box& vbx,
                            amrex::Array4<amrex::Real const> const& q,
                            amrex::Array4<amrex::Real const> const& flatn,
                            amrex::Array4<amrex::Real const> const& qaux,
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

    void mol_plm_reconstruct(const amrex::Box& bx,
                             const int idir,
                             amrex::Array4<amrex::Real const> const& q_arr,
                             amrex::Array4<amrex::Real const> const& flatn_arr,
                             amrex::Array4<amrex::Real const> const& src_q_arr,
                             amrex::Array4<amrex::Real> const& dq,
                             amrex::Array4<amrex::Real> const& qm,
                             amrex::Array4<amrex::Real> const& qp);

     void mol_ppm_reconstruct(const amrex::Box& bx,
                              const int idir,
                              amrex::Array4<amrex::Real const> const& q_arr,
                              amrex::Array4<amrex::Real const> const& flatn_arr,
                              amrex::Array4<amrex::Real> const& qm,
                              amrex::Array4<amrex::Real> const& qp);

     void mol_consup(const amrex::Box& bx,
                     amrex::Array4<amrex::Real const> const& shk,
                     amrex::Array4<amrex::Real const> const& uin,
                     amrex::Array4<amrex::Real const> const& srcU,
                     amrex::Array4<amrex::Real> const& update,
                     const amrex::Real dt,
                     amrex::Array4<amrex::Real const> const& flux1,
#if AMREX_SPACEDIM >= 2
                     amrex::Array4<amrex::Real const> const& flux2,
#endif
#if AMREX_SPACEDIM == 3
                     amrex::Array4<amrex::Real const> const& flux3,
#endif
                     amrex::Array4<amrex::Real const> const& area1,
#if AMREX_SPACEDIM >= 2
                     amrex::Array4<amrex::Real const> const& area2,
#endif
#if AMREX_SPACEDIM == 3
                     amrex::Array4<amrex::Real const> const& area3,
#endif
                     amrex::Array4<amrex::Real const> const& q1,
#if AMREX_SPACEDIM >= 2
                     amrex::Array4<amrex::Real const> const& q2,
#endif
#if AMREX_SPACEDIM == 3
                     amrex::Array4<amrex::Real const> const& q3,
#endif
                     amrex::Array4<amrex::Real const> const& vol);

     void
     mol_diffusive_flux(const amrex::Box& bx,
                        const int idir,
                        amrex::Array4<amrex::Real const> const& U,
                        amrex::Array4<amrex::Real const> const& cond,
                        amrex::Array4<amrex::Real> const& flux);

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
                      amrex::Real hdt, amrex::Real cdtdx_n,
                      amrex::Real cdtdx_t1, amrex::Real cdtdx_t2);

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
                             amrex::Real hdt, amrex::Real cdtdx_n,
                             amrex::Real cdtdx_t1, amrex::Real cdtdx_t2);

    void trace_ppm(const amrex::Box& bx,
                   const int idir,
                   amrex::Array4<amrex::Real const> const& q,
                   amrex::Array4<amrex::Real const> const& qaux,
                   amrex::Array4<amrex::Real const> const& srcQ,
                   amrex::Array4<amrex::Real const> const& flatn,
                   amrex::Array4<amrex::Real> const& qm,
                   amrex::Array4<amrex::Real> const& qp,
#if (AMREX_SPACEDIM < 3)
                   amrex::Array4<amrex::Real const> const& dloga,
#endif
                   const amrex::Box& vbx,
                   const amrex::Real dt);

    void ctu_plm_states(const amrex::Box& bx, const amrex::Box& vbx,
                        amrex::Array4<amrex::Real const> const& q_arr,
                        amrex::Array4<amrex::Real const> const& flatn,
                        amrex::Array4<amrex::Real const> const& qaux__arr,
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

    void trace_plm(const amrex::Box& bx,
                   const int idir,
                   amrex::Array4<amrex::Real const> const& q,
                   amrex::Array4<amrex::Real const> const& qaux,
                   amrex::Array4<amrex::Real const> const& flatn_arr,
                   amrex::Array4<amrex::Real> const& qm,
                   amrex::Array4<amrex::Real> const& qp,
#if (AMREX_SPACEDIM < 3)
                   amrex::Array4<amrex::Real const> const& dloga,
#endif
                   amrex::Array4<amrex::Real const> const& SrcQ,
                   const amrex::Box& vbx,
                   const amrex::Real dt);

    void trace_ppm_rad(const amrex::Box& bx,
                       const int idir,
                       amrex::Array4<amrex::Real const> const& q,
                       amrex::Array4<amrex::Real const> const& qaux,
                       amrex::Array4<amrex::Real const> const& srcQ,
                       amrex::Array4<amrex::Real const> const& flatn,
                       amrex::Array4<amrex::Real> const& qm,
                       amrex::Array4<amrex::Real> const& qp,
#if (AMREX_SPACEDIM < 3)
                       amrex::Array4<amrex::Real const> const& dloga,
#endif
                       const amrex::Box& vbx,
                       const amrex::Real dt);

    void consup_hydro(const amrex::Box& bx,
                      amrex::Array4<amrex::Real const> const& shk,
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

    void riemann_state(const amrex::Box& bx,
                       amrex::Array4<amrex::Real> const& qm,
                       amrex::Array4<amrex::Real> const& qp,
                       amrex::Array4<amrex::Real> const& qint,
#ifdef RADIATION
                       amrex::Array4<amrex::Real> const& lambda_int,
#endif
                       amrex::Array4<amrex::Real const> const& qaux,
                       const int idir, const int compute_gammas);

    void riemanncg(const amrex::Box& bx,
                   amrex::Array4<amrex::Real> const& ql,
                   amrex::Array4<amrex::Real> const& qr,
                   amrex::Array4<amrex::Real const> const& qaux,
                   amrex::Array4<amrex::Real> const& qint,
                   const int idir);

    void riemannus(const amrex::Box& bx,
                   amrex::Array4<amrex::Real> const& ql,
                   amrex::Array4<amrex::Real> const& qr,
                   amrex::Array4<amrex::Real const> const& qaux,
                   amrex::Array4<amrex::Real> const& qint,
#ifdef RADIATION
                   amrex::Array4<amrex::Real> const& lambda_int,
#endif
                   const int idir, const int compute_gammas);


    void HLLC(const amrex::Box& bx,
              amrex::Array4<amrex::Real const> const& ql,
              amrex::Array4<amrex::Real const> const& qr,
              amrex::Array4<amrex::Real const> const& qaux,
              amrex::Array4<amrex::Real> const& uflx,
              amrex::Array4<amrex::Real> const& qint,
              const int idir);

    AMREX_GPU_HOST_DEVICE
    void HLL(const amrex::Real* ql, const amrex::Real* qr,
             const amrex::Real cl, const amrex::Real cr,
             const int idir, const int coord,
             amrex::Real* f);

    void
    compute_flux_q(const amrex::Box& bx,
                   amrex::Array4<amrex::Real const> const& qint,
                   amrex::Array4<amrex::Real> const& F,
#ifdef RADIATION
                   amrex::Array4<amrex::Real const> const& lambda,
                   amrex::Array4<amrex::Real> const& rF,
#endif
                   const int idir, const int enforce_eos);

    void
    store_godunov_state(const amrex::Box& bx,
                        amrex::Array4<amrex::Real const> const& qint,
#ifdef RADIATION
                        amrex::Array4<amrex::Real const> const& lambda,
#endif
                        amrex::Array4<amrex::Real> const& qgdnv);


    AMREX_GPU_HOST_DEVICE
    void
    dflux(const GpuArray<Real, NUM_STATE>& u,
          const GpuArray<Real, NQ>& q,
          int dir, int coord,
          const GeometryData& geomdata,
          const GpuArray<Real, 3>& center,
          const GpuArray<int, 3>& idx,
          GpuArray<Real, NUM_STATE>& flux);

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
