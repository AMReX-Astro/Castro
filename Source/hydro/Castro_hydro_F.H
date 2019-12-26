#ifndef _Castro_hydro_F_H_
#define _Castro_hydro_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ca_ctoprim(const int* lo, const int* hi,
	          const amrex::Real* u, const int* u_lo, const int* u_hi,
#ifdef RADIATION
		  const amrex::Real* Er, const int* Er_lo, const int* Er_hi,
		  const amrex::Real* lam, const int* lam_lo, const int* lam_hi,
#endif
		  amrex::Real* q, const int* q_lo, const int* q_hi,
		  amrex::Real* qaux, const int* qa_lo, const int* qa_hi);

  void ca_srctoprim(const int* lo,   const int* hi,
		    const amrex::Real* q,   const int* q_lo,   const int* q_hi,
		    const amrex::Real* qaux, const int* qa_lo, const int* qa_hi,
		    const amrex::Real* src, const int* src_lo, const int* src_hi,
		          amrex::Real* srQ, const int* srQ_lo, const int* srQ_hi);

  void ca_compute_cfl(const int* lo, const int* hi,
                      const BL_FORT_FAB_ARG_ANYD(q),
                      const BL_FORT_FAB_ARG_ANYD(qaux),
                      const amrex::Real dt, const amrex::Real* dx,
                      amrex::Real* courno, const int verbose);

  void divu
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q),
     const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(div));

  void ca_uflatten
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn), int pres_comp);

  void ca_rad_flatten
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(flatg));

  void ctu_consup
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(shk),
     BL_FORT_FAB_ARG_3D(hydro_source),
     const BL_FORT_FAB_ARG_3D(flux0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(flux1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(flux2),
#endif
#ifdef RADIATION
     const BL_FORT_FAB_ARG_3D(Erborder),
     BL_FORT_FAB_ARG_3D(S_new),
     BL_FORT_FAB_ARG_3D(Er_out),
     const BL_FORT_FAB_ARG_3D(rad_flux0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(rad_flux1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(rad_flux2),
#endif
     int* nstep_fsp,
#endif
     const BL_FORT_FAB_ARG_3D(qe0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(qe1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(qe2),
#endif
     const BL_FORT_FAB_ARG_3D(area0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(area1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(area2),
#endif
     const BL_FORT_FAB_ARG_3D(volume),
     const amrex::Real* dx, const amrex::Real dt);

  void scale_flux
    (const int* lo, const int* hi,
#if AMREX_SPACEDIM == 1
     const BL_FORT_FAB_ARG_3D(qe0),
#endif
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(area),
     const amrex::Real dt);

  void scale_rad_flux
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(rad_flux),
     BL_FORT_FAB_ARG_3D(area),
     const amrex::Real dt);

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
     const BL_FORT_FAB_ARG_3D(qr), const int nc, const int comp,
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

  void apply_av
    (const int* lo, const int* hi,
     const int idir, const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(div),
     BL_FORT_FAB_ARG_3D(Sborder),
     BL_FORT_FAB_ARG_3D(flux));

#ifdef RADIATION
  void apply_av_rad
    (const int* lo, const int* hi,
     int idir, const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(div),
     BL_FORT_FAB_ARG_3D(Er),
     BL_FORT_FAB_ARG_3D(radflux));
#endif

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

  void normalize_species_fluxes
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(flux));

  void ca_shock
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(shk),
     const amrex::Real* dx);

  void ca_ctu_update
    (const int* lo, const int* hi,
     const int* is_finest_level,
     const amrex::Real time,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(stateout),
#ifdef RADIATION
     BL_FORT_FAB_ARG_3D(Er),
     BL_FORT_FAB_ARG_3D(Erout),
#endif
     BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(qaux),
     const BL_FORT_FAB_ARG_3D(src_q),
     BL_FORT_FAB_ARG_3D(update),
     const amrex::Real* dx, const amrex::Real dt,
     BL_FORT_FAB_ARG_3D(xflux),
#if AMREX_SPACEDIM >= 2
     BL_FORT_FAB_ARG_3D(yflux),
#endif
#if AMREX_SPACEDIM == 3
     BL_FORT_FAB_ARG_3D(zflux),
#endif
#ifdef RADIATION
     BL_FORT_FAB_ARG_3D(rxflux),
#if AMREX_SPACEDIM >= 2
     BL_FORT_FAB_ARG_3D(ryflux),
#endif
#if AMREX_SPACEDIM == 3
     BL_FORT_FAB_ARG_3D(rzflux),
#endif
#endif
     const BL_FORT_FAB_ARG_3D(xarea),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(yarea),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(zarea),
#endif
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG_3D(pradial),
     const BL_FORT_FAB_ARG_3D(dloga),
#endif
     const BL_FORT_FAB_ARG_3D(volume),
#ifdef RADIATION
     const int* priv_nstep_fsp,
#endif
     amrex::Real& mass_lost,
     amrex::Real& xmom_lost,
     amrex::Real& ymom_lost,
     amrex::Real& zmom_lost,
     amrex::Real& eden_lost,
     amrex::Real& xang_lost,
     amrex::Real& yang_lost,
     amrex::Real& zang_lost);

  void ctu_ppm_states
     (const int* lo, const int* hi,
      const int* vlo, const int* vhi,
      const BL_FORT_FAB_ARG_3D(q),
      const BL_FORT_FAB_ARG_3D(flatn),
      const BL_FORT_FAB_ARG_3D(qaux),
      const BL_FORT_FAB_ARG_3D(srcQ),
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

  void transx_on_ystates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qym),
                         BL_FORT_FAB_ARG_3D(qymo),
                         const BL_FORT_FAB_ARG_3D(qyp),
                         BL_FORT_FAB_ARG_3D(qypo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fx),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfx),
#endif
                         const BL_FORT_FAB_ARG_3D(qx),
#if AMREX_SPACEDIM == 2
                         const BL_FORT_FAB_ARG_3D(area1),
                         const BL_FORT_FAB_ARG_3D(vol),
#endif
                         const amrex::Real hdt, const amrex::Real cdtdx);

  void transx_on_zstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qzm),
                         BL_FORT_FAB_ARG_3D(qzmo),
                         const BL_FORT_FAB_ARG_3D(qzp),
                         BL_FORT_FAB_ARG_3D(qzpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fx),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfx),
#endif
                         const BL_FORT_FAB_ARG_3D(qx),
                         const amrex::Real hdt, const amrex::Real cdtdx);

  void transy_on_xstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qxm),
                         BL_FORT_FAB_ARG_3D(qxmo),
                         const BL_FORT_FAB_ARG_3D(qxp),
                         BL_FORT_FAB_ARG_3D(qxpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fy),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfy),
#endif
                         const BL_FORT_FAB_ARG_3D(qy),
                         const amrex::Real cdtdy);

  void transy_on_zstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qzm),
                         BL_FORT_FAB_ARG_3D(qzmo),
                         const BL_FORT_FAB_ARG_3D(qzp),
                         BL_FORT_FAB_ARG_3D(qzpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fy),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfy),
#endif
                         const BL_FORT_FAB_ARG_3D(qy),
                         const amrex::Real cdtdy);

  void transz_on_xstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qxm),
                         BL_FORT_FAB_ARG_3D(qxmo),
                         const BL_FORT_FAB_ARG_3D(qxp),
                         BL_FORT_FAB_ARG_3D(qxpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fz),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfz),
#endif
                         const BL_FORT_FAB_ARG_3D(qz),
                         const amrex::Real cdtdz);

  void transz_on_ystates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qym),
                         BL_FORT_FAB_ARG_3D(qymo),
                         const BL_FORT_FAB_ARG_3D(qyp),
                         BL_FORT_FAB_ARG_3D(qypo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fz),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfz),
#endif
                         const BL_FORT_FAB_ARG_3D(qz),
                         const amrex::Real cdtdz);

  void transyz(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_3D(qm),
               BL_FORT_FAB_ARG_3D(qmo),
               const BL_FORT_FAB_ARG_3D(qp),
               BL_FORT_FAB_ARG_3D(qpo),
               const BL_FORT_FAB_ARG_3D(qaux),
               const BL_FORT_FAB_ARG_3D(fyz),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfyz),
#endif
               const BL_FORT_FAB_ARG_3D(fzy),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfzy),
#endif
               const BL_FORT_FAB_ARG_3D(qy),
               const BL_FORT_FAB_ARG_3D(qz),
               const amrex::Real hdt, const amrex::Real cdtdy, const amrex::Real cdtdz);

  void transxz(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_3D(qm),
               BL_FORT_FAB_ARG_3D(qmo),
               const BL_FORT_FAB_ARG_3D(qp),
               BL_FORT_FAB_ARG_3D(qpo),
               const BL_FORT_FAB_ARG_3D(qaux),
               const BL_FORT_FAB_ARG_3D(fxz),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfxz),
#endif
               const BL_FORT_FAB_ARG_3D(fzx),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfzx),
#endif
               const BL_FORT_FAB_ARG_3D(qx),
               const BL_FORT_FAB_ARG_3D(qz),
               const amrex::Real hdt, const amrex::Real cdtdx, const amrex::Real cdtdz);

  void transxy(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_3D(qm),
               BL_FORT_FAB_ARG_3D(qmo),
               const BL_FORT_FAB_ARG_3D(qp),
               BL_FORT_FAB_ARG_3D(qpo),
               const BL_FORT_FAB_ARG_3D(qaux),
               const BL_FORT_FAB_ARG_3D(fxy),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfxy),
#endif
               const BL_FORT_FAB_ARG_3D(fyx),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfyx),
#endif
               const BL_FORT_FAB_ARG_3D(qx),
               const BL_FORT_FAB_ARG_3D(qy),
               const amrex::Real hdt, const amrex::Real cdtdx, const amrex::Real cdtdy);

  void ca_mol_plm_reconstruct
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(dq),
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp),
     const amrex::Real* dx,
     const int* domlo, const int* domhi);

  void ca_mol_ppm_reconstruct
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp),
     const amrex::Real* dx);

  void ca_mol_consup
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(statein),
     BL_FORT_FAB_ARG_3D(stateout),
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

  void ca_fourth_single_stage
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(statein),
     BL_FORT_FAB_ARG_3D(stateout),
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(q_bar),
     BL_FORT_FAB_ARG_3D(qaux),
     BL_FORT_FAB_ARG_3D(qaux_bar),
#ifdef DIFFUSION
     BL_FORT_FAB_ARG_3D(T_Cc),
#endif
     const BL_FORT_FAB_ARG_3D(srcU),
     BL_FORT_FAB_ARG_3D(update),
     const amrex::Real* dx, const amrex::Real* dt,
     BL_FORT_FAB_ARG_3D(flux1),
#if AMREX_SPACEDIM >= 2
     BL_FORT_FAB_ARG_3D(flux2),
#endif
#if AMREX_SPACEDIM >= 3
     BL_FORT_FAB_ARG_3D(flux3),
#endif
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG_3D(pradial),
#endif
     const BL_FORT_FAB_ARG_3D(area1),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(area2),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(area3),
#endif
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG_3D(dloga),
#endif
     const BL_FORT_FAB_ARG_3D(volume),
     const int&  verbose);

  void ca_mol_diffusive_flux
    (const int* lo, const int* hi,
     const int idir,
     const BL_FORT_FAB_ARG_3D(U),
     const BL_FORT_FAB_ARG_3D(cond),
     BL_FORT_FAB_ARG_3D(flux),
     const amrex::Real* dx);

#ifndef AMREX_USE_CUDA
  void ca_sdc_update_advection_o2_lobatto
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const amrex::Real* dt,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     const int* m_start);

  void ca_sdc_update_advection_o2_radau
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const amrex::Real* dt,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     BL_FORT_FAB_ARG_3D(A_2_old),
     const int* m_start);

  void ca_sdc_update_advection_o4_lobatto
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const amrex::Real* dt,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     BL_FORT_FAB_ARG_3D(A_2_old),
     const int* m_start);

  void ca_sdc_update_advection_o4_radau
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const amrex::Real* dt,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     BL_FORT_FAB_ARG_3D(A_2_old),
     BL_FORT_FAB_ARG_3D(A_3_old),
     const int* m_start);

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
