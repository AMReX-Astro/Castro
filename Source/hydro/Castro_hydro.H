
///
/// Calculate primitive variables from conserved variables (uses StateData)
///
/// @param time     current time
///
    void cons_to_prim(const amrex::Real time);

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
/// @param dt timestep
///
    void check_for_cfl_violation(const amrex::Real dt);

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
                     amrex::Array4<amrex::Real const> const q,
                     amrex::Array4<amrex::Real const> const qaux,
                     amrex::Array4<amrex::Real const> const src,
                     amrex::Array4<amrex::Real> const srcQ);

    void uflatten(const amrex::Box& bx,
                  amrex::Array4<amrex::Real const> const q,
                  amrex::Array4<amrex::Real> const flatn,
                  const int pres_comp);

    void shock(const amrex::Box& bx,
               amrex::Array4<amrex::Real const> const q,
               amrex::Array4<amrex::Real> const shk);

    void divu(const amrex::Box& bx,
              amrex::Array4<amrex::Real const> const q,
              amrex::Array4<amrex::Real> const div);

    void apply_av(const amrex::Box& bx,
                  const int idir,
                  amrex::Array4<amrex::Real const> const div,
                  amrex::Array4<amrex::Real const> const uin,
                  amrex::Array4<amrex::Real> const flux);

#ifdef RADIATION
    void apply_av_rad(const amrex::Box& bx,
                      const int idir,
                      amrex::Array4<amrex::Real const> const div,
                      amrex::Array4<amrex::Real const> const Erin,
                      amrex::Array4<amrex::Real> const radflux);
#endif

    void ctu_ppm_states(const amrex::Box& bx, const amrex::Box& vbx,
                        amrex::Array4<amrex::Real const> const q,
                        amrex::Array4<amrex::Real const> const flatn,
                        amrex::Array4<amrex::Real const> const qaux,
                        amrex::Array4<amrex::Real const> const srcQ,
                        amrex::Array4<amrex::Real> const qxm,
                        amrex::Array4<amrex::Real> const qxp,
#if AMREX_SPACEDIM >= 2
                        amrex::Array4<amrex::Real> const qym,
                        amrex::Array4<amrex::Real> const qyp,
#endif
#if AMREX_SPACEDIM == 3
                        amrex::Array4<amrex::Real> const qzm,
                        amrex::Array4<amrex::Real> const qzp,
#endif
#if AMREX_SPACEDIM < 3
                        amrex::Array4<amrex::Real const> const dloga,
#endif
                        const amrex::Real dt);

     void mol_ppm_reconstruct(const amrex::Box& bx,
                              const int idir,
                              amrex::Array4<amrex::Real const> const q_arr,
                              amrex::Array4<amrex::Real const> const flatn_arr,
                              amrex::Array4<amrex::Real> const qm,
                              amrex::Array4<amrex::Real> const qp);


    void trace_ppm(const amrex::Box& bx,
                   const int idir,
                   amrex::Array4<amrex::Real const> const q,
                   amrex::Array4<amrex::Real const> const qaux,
                   amrex::Array4<amrex::Real const> const srcQ,
                   amrex::Array4<amrex::Real const> const flatn,
                   amrex::Array4<amrex::Real> const qm,
                   amrex::Array4<amrex::Real> const qp,
#if (AMREX_SPACEDIM < 3)
                   amrex::Array4<amrex::Real const> const dloga,
#endif
                   const amrex::Box& vbx,
                   const amrex::Real dt);

    void consup_hydro(const amrex::Box& bx,
                      amrex::Array4<amrex::Real const> const shk,
                      amrex::Array4<amrex::Real> const update,
                      amrex::Array4<amrex::Real> const flux0,
                      amrex::Array4<amrex::Real const> const qx,
                      amrex::Array4<amrex::Real const> const area0,
#if AMREX_SPACEDIM >= 2
                      amrex::Array4<amrex::Real> const flux1,
                      amrex::Array4<amrex::Real const> const qy,
                      amrex::Array4<amrex::Real const> const area1,
#endif
#if AMREX_SPACEDIM == 3
                      amrex::Array4<amrex::Real> const flux2,
                      amrex::Array4<amrex::Real const> const qz,
                      amrex::Array4<amrex::Real const> const area2,
#endif
                      amrex::Array4<amrex::Real const> const vol,
                      const amrex::Real dt);



    void cmpflx_plus_godunov(const amrex::Box& bx,
                             amrex::Array4<amrex::Real> const qm,
                             amrex::Array4<amrex::Real> const qp,
                             amrex::Array4<amrex::Real> const flx,
                             amrex::Array4<amrex::Real> const qint,
#ifdef RADIATION
                             amrex::Array4<amrex::Real> const rflx,
                             amrex::Array4<amrex::Real> const lambda_int,
#endif
                             amrex::Array4<amrex::Real> const qgdnv,
                             amrex::Array4<amrex::Real const> const qaux,
                             amrex::Array4<amrex::Real const> const shk,
                             const int idir);

    void riemann_state(const amrex::Box& bx,
                       amrex::Array4<amrex::Real> const qm,
                       amrex::Array4<amrex::Real> const qp,
                       amrex::Array4<amrex::Real> const qint,
#ifdef RADIATION
                       amrex::Array4<amrex::Real> lambda_int,
#endif
                       amrex::Array4<amrex::Real const> const qaux,
                       const int idir, const int compute_gammas);

    void riemanncg(const amrex::Box& bx,
                   amrex::Array4<amrex::Real> const ql,
                   amrex::Array4<amrex::Real> const qr,
                   amrex::Array4<amrex::Real const> const qaux,
                   amrex::Array4<amrex::Real> const qint,
                   const int idir);

    void riemannus(const amrex::Box& bx,
                   amrex::Array4<amrex::Real> const ql,
                   amrex::Array4<amrex::Real> const qr,
                   amrex::Array4<amrex::Real const> const qaux,
                   amrex::Array4<amrex::Real> const qint,
#ifdef RADIATION
                   amrex::Array4<amrex::Real> const lambda_int,
#endif
                   const int idir, const int compute_gammas);


    void HLLC(const amrex::Box& bx,
              amrex::Array4<amrex::Real const> const ql,
              amrex::Array4<amrex::Real const> const qr,
              amrex::Array4<amrex::Real const> const qaux,
              amrex::Array4<amrex::Real> const uflx,
              amrex::Array4<amrex::Real> const qint,
              const int idir);

    void HLL(const amrex::Real* ql, const amrex::Real* qr,
             const amrex::Real cl, const amrex::Real cr,
             const int idir,
             amrex::Real* f);

    AMREX_GPU_HOST_DEVICE
    void
    wsqge(const amrex::Real p, const amrex::Real v,
          const amrex::Real gam, const amrex::Real gdot, amrex::Real& gstar,
          const amrex::Real gmin, const amrex::Real gmax, const amrex::Real csq,
          amrex::Real& pstar, amrex::Real& wsq);

    AMREX_GPU_HOST_DEVICE
    void
    pstar_bisection(amrex::Real& pstar_lo, amrex::Real& pstar_hi,
                    const amrex::Real ul, const amrex::Real pl, const amrex::Real taul,
                    const amrex::Real gamel, const amrex::Real clsql,
                    const amrex::Real ur, const amrex::Real pr, const amrex::Real taur,
                    const amrex::Real gamer, const amrex::Real clsqr,
                    const amrex::Real gdot, const amrex::Real gmin, const amrex::Real gmax,
                    amrex::Real& pstar, amrex::Real& gamstar,
                    bool& converged, amrex::GpuArray<amrex::Real, 2*HISTORY_SIZE>& pstar_hist_extra);

    AMREX_GPU_HOST_DEVICE
    void
    cons_state(const amrex::Real* q, amrex::Real* U);

    AMREX_GPU_HOST_DEVICE
    void
    HLLC_state(const int idir, const amrex::Real S_k, const amrex::Real S_c,
               const amrex::Real* q, amrex::Real* U);

    AMREX_GPU_HOST_DEVICE
    void
    compute_flux_q(const amrex::Box& bx,
                   amrex::Array4<amrex::Real const> const qint,
                   amrex::Array4<amrex::Real> const F,
#ifdef RADIATION
                   amrex::Array4<amrex::Real const> const lambda,
                   amrex::Array4<amrex::Real> const rF,
#endif
                   const int idir, const int enforce_eos);

    AMREX_GPU_HOST_DEVICE
    void
    store_godunov_state(const amrex::Box& bx,
                        amrex::Array4<amrex::Real const> const qint,
#ifdef RADIATION
                        amrex::Array4<amrex::Real const> const lambda,
#endif
                        amrex::Array4<amrex::Real> const qgdnv);

    AMREX_GPU_HOST_DEVICE
    void
    compute_flux(const int idir, const amrex::Real bnd_fac,
                 const amrex::Real* U, const amrex::Real p,
                 amrex::Real* F);



    void normalize_species_fluxes(const amrex::Box& bx, amrex::Array4<amrex::Real> const flux);

    void scale_flux(const amrex::Box& bx,
#if AMREX_SPACEDIM == 1
                    amrex::Array4<amrex::Real const> const qint,
#endif
                    amrex::Array4<amrex::Real> const flux,
                    amrex::Array4<amrex::Real const> const area,
                    const amrex::Real dt);

#ifdef RADIATION
    void scale_rad_flux(const amrex::Box& bx,
                        amrex::Array4<amrex::Real> const rflux,
                        amrex::Array4<amrex::Real const> area,
                        const amrex::Real dt);
#endif

    void reset_edge_state_thermo(const amrex::Box& bx,
                                 amrex::Array4<amrex::Real> const qedge);

    void edge_state_temp_to_pres(const Box& bx,
                                 amrex::Array4<amrex::Real> const qm,
                                 amrex::Array4<amrex::Real> const qp);

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
