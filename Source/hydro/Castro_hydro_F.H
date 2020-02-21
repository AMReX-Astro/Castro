#ifndef _Castro_hydro_F_H_
#define _Castro_hydro_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ca_compute_cfl(const int* lo, const int* hi,
                      const BL_FORT_FAB_ARG_ANYD(q),
                      const BL_FORT_FAB_ARG_ANYD(qaux),
                      const amrex::Real dt, const amrex::Real* dx,
                      amrex::Real* courno, const int verbose);

#ifdef RADIATION
  void ctu_rad_consup
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(hydro_source),
     const BL_FORT_FAB_ARG_3D(Erborder),
     BL_FORT_FAB_ARG_3D(S_new),
     BL_FORT_FAB_ARG_3D(Er_out),
     const BL_FORT_FAB_ARG_3D(rad_flux0),
     const BL_FORT_FAB_ARG_3D(qe0),
     const BL_FORT_FAB_ARG_3D(area0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(rad_flux1),
     const BL_FORT_FAB_ARG_3D(qe1),
     const BL_FORT_FAB_ARG_3D(area1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(rad_flux2),
     const BL_FORT_FAB_ARG_3D(qe2),
     const BL_FORT_FAB_ARG_3D(area2),
#endif
     int* nstep_fsp,
     const BL_FORT_FAB_ARG_3D(volume),
     const amrex::Real* dx, const amrex::Real dt);
#endif

  void ca_track_grid_losses
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(flux0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(flux1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(flux2),
#endif
     amrex::Real* mass_lost, amrex::Real* xmom_lost, amrex::Real* ymom_lost, amrex::Real* zmom_lost,
     amrex::Real* eden_lost, amrex::Real* xang_lost, amrex::Real* yang_lost, amrex::Real* zang_lost);

  void cmpflx_plus_godunov
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(ql),
     const BL_FORT_FAB_ARG_3D(qr),
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(q_int),
#ifdef RADIATION
     BL_FORT_FAB_ARG_3D(radflux),
     BL_FORT_FAB_ARG_3D(lambda_int),
#endif
     BL_FORT_FAB_ARG_3D(qe),
     const BL_FORT_FAB_ARG_3D(qaux),
     const BL_FORT_FAB_ARG_3D(shk),
     const int idir, const int* domlo, const int* domhi);

  void limit_hydro_fluxes_on_small_dens
    (const int* lo, const int* hi,
     const int idir,
     BL_FORT_FAB_ARG_3D(Sborder),
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(volume),
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(area),
     const amrex::Real dt, const amrex::Real* dx);

  void limit_hydro_fluxes_on_large_vel
    (const int* lo, const int* hi,
     const int idir,
     BL_FORT_FAB_ARG_3D(Sborder),
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(volume),
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(area),
     const amrex::Real dt, const amrex::Real* dx);

  void ctu_plm_states
     (const int* lo, const int* hi,
      const int* vlo, const int* vhi,
      const BL_FORT_FAB_ARG_3D(q),
      const BL_FORT_FAB_ARG_3D(flatn),
      const BL_FORT_FAB_ARG_3D(qaux),
      const BL_FORT_FAB_ARG_3D(srcQ),
      BL_FORT_FAB_ARG_3D(dq),
      BL_FORT_FAB_ARG_3D(qxm),
      BL_FORT_FAB_ARG_3D(qxp),
#if AMREX_SPACEDIM >= 2
      BL_FORT_FAB_ARG_3D(qym),
      BL_FORT_FAB_ARG_3D(qyp),
#endif
#if AMREX_SPACEDIM == 3
      BL_FORT_FAB_ARG_3D(qzm),
      BL_FORT_FAB_ARG_3D(qzp),
#endif
      const amrex::Real* dx, const amrex::Real dt,
#if AMREX_SPACEDIM < 3
      const BL_FORT_FAB_ARG_3D(dloga),
#endif
      const int* domlo, const int* domhi);

  void ca_mol_plm_reconstruct
    (const int* lo, const int* hi,
     const int idir,
     BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(dq),
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp),
     const amrex::Real* dx,
     const int* domlo, const int* domhi);

  void edge_state_temp_to_pres
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp));

  void ca_mol_consup
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(shk),
     const BL_FORT_FAB_ARG_3D(statein),
     const BL_FORT_FAB_ARG_3D(srcU),
     BL_FORT_FAB_ARG_3D(update),
     const amrex::Real* dx, const amrex::Real dt,
     BL_FORT_FAB_ARG_3D(flux1),
#if AMREX_SPACEDIM >= 2
     BL_FORT_FAB_ARG_3D(flux2),
#endif
#if AMREX_SPACEDIM == 3
     BL_FORT_FAB_ARG_3D(flux3),
#endif
     const BL_FORT_FAB_ARG_3D(area1),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(area2),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(area3),
#endif
     const BL_FORT_FAB_ARG_3D(q1),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(q2),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(q3),
#endif
     const BL_FORT_FAB_ARG_3D(volume));

  void ca_fourth_interfaces
    (const int* lo, const int* hi,
     const int idir, const int ncomp,
     const BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(q_int),
     const int* domlo, const int* domhi);

  void ca_states
    (const int* lo, const int* hi,
     const int idir, const int ncomp,
     const BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(q_int),
     const BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp),
     const int* domlo, const int* domhi);

  void ca_mol_diffusive_flux
    (const int* lo, const int* hi,
     const int idir,
     const BL_FORT_FAB_ARG_3D(U),
     const BL_FORT_FAB_ARG_3D(cond),
     BL_FORT_FAB_ARG_3D(flux),
     const amrex::Real* dx);

  void trans_laplacian
    (const int i, const int j, const int k, const int n,
     const int idir,
     const BL_FORT_FAB_ARG_3D(q_avg), const int nq,
     amrex::Real* lap,
     const int* domlo, const int* domhi);

  void avisc
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(qaux),
     const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(avis),
     const int idir);

  void add_diffusive_flux
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(q), const int nq, const int comp,
     const BL_FORT_FAB_ARG_3D(q_avg),
     BL_FORT_FAB_ARG_3D(f_avg),
     const amrex::Real* dx, const int idir, const int is_avg);

#ifndef AMREX_USE_CUDA

#ifdef REACTIONS
  void ca_sdc_compute_C2_lobatto
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const amrex::Real* dt,
     const BL_FORT_FAB_ARG_3D(A_m),
     const BL_FORT_FAB_ARG_3D(A_0_old),
     const BL_FORT_FAB_ARG_3D(A_1_old),
     const BL_FORT_FAB_ARG_3D(R_0_old),
     const BL_FORT_FAB_ARG_3D(R_1_old),
     BL_FORT_FAB_ARG_3D(C2),
     const int* m_start);

  void ca_sdc_compute_C2_radau
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const amrex::Real* dt,
     const BL_FORT_FAB_ARG_3D(A_m),
     const BL_FORT_FAB_ARG_3D(A_0_old),
     const BL_FORT_FAB_ARG_3D(A_1_old),
     const BL_FORT_FAB_ARG_3D(A_2_old),
     const BL_FORT_FAB_ARG_3D(R_0_old),
     const BL_FORT_FAB_ARG_3D(R_1_old),
     const BL_FORT_FAB_ARG_3D(R_2_old),
     BL_FORT_FAB_ARG_3D(C2),
     const int* m_start);

  void ca_sdc_update_o2
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     const BL_FORT_FAB_ARG_3D(A_m),
     const BL_FORT_FAB_ARG_3D(R_m_old),
     const BL_FORT_FAB_ARG_3D(C2),
     const int* sdc_iteration,
     const int* m_start);

  void ca_sdc_update_centers_o4
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(U_center),
     BL_FORT_FAB_ARG_3D(U_new_center),
     const BL_FORT_FAB_ARG_3D(C_center),
     const int* sdc_iteration);

  void ca_sdc_conservative_update
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(U_old),
     BL_FORT_FAB_ARG_3D(U_new),
     const BL_FORT_FAB_ARG_3D(C),
     const BL_FORT_FAB_ARG_3D(R_new));

  void ca_sdc_compute_C4_lobatto(const int* lo, const int* hi,
                                 const amrex::Real* dt_m,
                                 const amrex::Real* dt,
                                 const BL_FORT_FAB_ARG_3D(A_m),
                                 const BL_FORT_FAB_ARG_3D(A_0),
                                 const BL_FORT_FAB_ARG_3D(A_1),
                                 const BL_FORT_FAB_ARG_3D(A_2),
                                 const BL_FORT_FAB_ARG_3D(R_0),
                                 const BL_FORT_FAB_ARG_3D(R_1),
                                 const BL_FORT_FAB_ARG_3D(R_2),
                                 BL_FORT_FAB_ARG_3D(C),
                                 const int* m_start);

  void ca_sdc_compute_C4_radau(const int* lo, const int* hi,
                               const amrex::Real* dt_m,
                               const amrex::Real* dt,
                               const BL_FORT_FAB_ARG_3D(A_m),
                               const BL_FORT_FAB_ARG_3D(A_0),
                               const BL_FORT_FAB_ARG_3D(A_1),
                               const BL_FORT_FAB_ARG_3D(A_2),
                               const BL_FORT_FAB_ARG_3D(A_3),
                               const BL_FORT_FAB_ARG_3D(R_0),
                               const BL_FORT_FAB_ARG_3D(R_1),
                               const BL_FORT_FAB_ARG_3D(R_2),
                               const BL_FORT_FAB_ARG_3D(R_3),
                               BL_FORT_FAB_ARG_3D(C),
                               const int* m_start);

  void ca_sdc_compute_initial_guess(const int* lo, const int* hi,
                                    const BL_FORT_FAB_ARG_3D(U_old),
                                    const BL_FORT_FAB_ARG_3D(U_new),
                                    const BL_FORT_FAB_ARG_3D(A_old),
                                    const BL_FORT_FAB_ARG_3D(R_old),
                                    BL_FORT_FAB_ARG_3D(U_guess),
                                    const amrex::Real* dt_m, const int* sdc_iteration);
#endif
#endif


#ifdef __cplusplus
}
#endif

#endif
