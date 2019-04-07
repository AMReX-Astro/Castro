#ifndef _Castro_hydro_F_H_
#define _Castro_hydro_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  ///
  /// Convert the conservative state U to the primitive state, q_core,
  /// q_pass, q_rad
  ///
  void ca_ctoprim
    (const int* lo, const int* hi,
     const amrex::Real* U, const int* u_lo, const int* u_hi,
#ifdef RADIATION
     const amrex::Real* Er, const int* Er_lo, const int* Er_hi,
     const amrex::Real* lam, const int* lam_lo, const int* lam_hi,
#endif
     amrex::Real* q_core, const int* qc_lo, const int* qc_hi,
     amrex::Real* q_pass, const int* qp_lo, const int* qp_hi,
#ifdef RADIATION
     amrex::Real* q_rad, const int* qr_lo, const int* qr_hi,
#endif
     amrex::Real* qaux, const int* qa_lo, const int* qa_hi);

  ///
  /// Convert the conservative source terms to their primitive
  /// variable equivalents
  ///
  void ca_srctoprim
     (const int* lo, const int* hi,
      const amrex::Real* q_core, const int* qc_lo, const int* qc_hi,
      const amrex::Real* q_pass, const int* qp_lo, const int* qp_hi,
      const amrex::Real* qaux, const int* qa_lo, const int* qa_hi,
      const amrex::Real* src, const int* src_lo, const int* src_hi,
      amrex::Real* q_core_src, const int* qcs_lo, const int* qcs_hi
#ifdef PRIM_SPECIES_HAVE_SOURCES
     ,amrex::Real* q_pass_src, const int* qps_lo, const int* qps_hi
#endif
      );

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
     const BL_FORT_FAB_ARG_3D(Sborder),
     const BL_FORT_FAB_ARG_3D(shk),
     BL_FORT_FAB_ARG_3D(S_new),
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
     const BL_FORT_FAB_ARG_3D(pdivu),
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
     const BL_FORT_FAB_ARG_3D(ql_core),
     const BL_FORT_FAB_ARG_3D(ql_pass),
#ifdef RADIATION
     const BL_FORT_FAB_ARG_3D(ql_rad),
#endif
     const BL_FORT_FAB_ARG_3D(qr_core),
     const BL_FORT_FAB_ARG_3D(qr_pass),
#ifdef RADIATION
     const BL_FORT_FAB_ARG_3D(qr_pass),
#endif
     const int nc, const int comp,
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(q_int_core),
     BL_FORT_FAB_ARG_3D(q_int_pass),
#ifdef RADIATION
     BL_FORT_FAB_ARG_3D(q_int_rad),
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

  void normalize_species_fluxes
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(flux));

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
      const BL_FORT_FAB_ARG_3D(q_core),
      const BL_FORT_FAB_ARG_3D(q_pass),
      const BL_FORT_FAB_ARG_3D(flatn),
      const BL_FORT_FAB_ARG_3D(qaux),
      const BL_FORT_FAB_ARG_3D(q_core_src),
#ifdef PRIM_SPECIES_HAVE_SOURCES
      const BL_FORT_FAB_ARG_3D(q_pass_src),
#endif
      const BL_FORT_FAB_ARG_3D(shk),
      BL_FORT_FAB_ARG_3D(Ip_core),
      BL_FORT_FAB_ARG_3D(Ip_pass),
      BL_FORT_FAB_ARG_3D(Im_core),
      BL_FORT_FAB_ARG_3D(Im_pass),
      BL_FORT_FAB_ARG_3D(Ip_core_src),
#ifdef PRIM_SPECIES_HAVE_SOURCES
      BL_FORT_FAB_ARG_3D(Ip_pass_src),
#endif
      BL_FORT_FAB_ARG_3D(Im_core_src),
#ifdef PRIM_SPECIES_HAVE_SOURCES
      BL_FORT_FAB_ARG_3D(Im_pass_src),
#endif
      BL_FORT_FAB_ARG_3D(Ip_gc),
      BL_FORT_FAB_ARG_3D(Im_gc),
      BL_FORT_FAB_ARG_3D(sm),
      BL_FORT_FAB_ARG_3D(sp),
      BL_FORT_FAB_ARG_3D(qxm_core),
      BL_FORT_FAB_ARG_3D(qxm_pass),
      BL_FORT_FAB_ARG_3D(qxp_core),
      BL_FORT_FAB_ARG_3D(qxp_pass),
#if AMREX_SPACEDIM >= 2
      BL_FORT_FAB_ARG_3D(qym_core),
      BL_FORT_FAB_ARG_3D(qym_pass),
      BL_FORT_FAB_ARG_3D(qyp_core),
      BL_FORT_FAB_ARG_3D(qyp_pass),
#endif
#if AMREX_SPACEDIM == 3
      BL_FORT_FAB_ARG_3D(qzm_core),
      BL_FORT_FAB_ARG_3D(qzm_pass),
      BL_FORT_FAB_ARG_3D(qzp_core),
      BL_FORT_FAB_ARG_3D(qzp_pass),
#endif
      const amrex::Real* dx, const amrex::Real dt,
#if AMREX_SPACEDIM < 3
      const BL_FORT_FAB_ARG_3D(dloga),
#endif
      const int* domlo, const int* domhi);

  void ctu_plm_states
     (const int* lo, const int* hi,
      const int* vlo, const int* vhi,
      const BL_FORT_FAB_ARG_3D(q_core),
      const BL_FORT_FAB_ARG_3D(q_pass),
      const BL_FORT_FAB_ARG_3D(flatn),
      const BL_FORT_FAB_ARG_3D(qaux),
      const BL_FORT_FAB_ARG_3D(q_core_src),
#ifdef PRIM_SPECIES_HAVE_SOURCES
      const BL_FORT_FAB_ARG_3D(q_pass_src),
#endif
      const BL_FORT_FAB_ARG_3D(shk),
      BL_FORT_FAB_ARG_3D(dq_core),
      BL_FORT_FAB_ARG_3D(dq_pass),
      BL_FORT_FAB_ARG_3D(qxm_core),
      BL_FORT_FAB_ARG_3D(qxm_pass),
      BL_FORT_FAB_ARG_3D(qxp_core),
      BL_FORT_FAB_ARG_3D(qxp_pass),
#if AMREX_SPACEDIM >= 2
      BL_FORT_FAB_ARG_3D(qym_core),
      BL_FORT_FAB_ARG_3D(qym_pass),
      BL_FORT_FAB_ARG_3D(qyp_core),
      BL_FORT_FAB_ARG_3D(qyp_pass),
#endif
#if AMREX_SPACEDIM == 3
      BL_FORT_FAB_ARG_3D(qzm_core),
      BL_FORT_FAB_ARG_3D(qzm_pass),
      BL_FORT_FAB_ARG_3D(qzp_core),
      BL_FORT_FAB_ARG_3D(qzp_pass),
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
     BL_FORT_FAB_ARG_3D(q_core),
     BL_FORT_FAB_ARG_3D(q_pass),
     const BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(shk),
     BL_FORT_FAB_ARG_3D(dq_core),
     BL_FORT_FAB_ARG_3D(dq_pass),
     BL_FORT_FAB_ARG_3D(qm_core),
     BL_FORT_FAB_ARG_3D(qm_pass),
     BL_FORT_FAB_ARG_3D(qp_core),
     BL_FORT_FAB_ARG_3D(qp_pass),
     const amrex::Real* dx);

  void ca_mol_ppm_reconstruct
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q_core),
     BL_FORT_FAB_ARG_3D(q_pass),
     const BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(shk),
     BL_FORT_FAB_ARG_3D(qm_core),
     BL_FORT_FAB_ARG_3D(qm_pass),
     BL_FORT_FAB_ARG_3D(qp_core),
     BL_FORT_FAB_ARG_3D(qp_pass),
     const amrex::Real* dx);

  void ca_mol_thermo_states
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(qm_core),
     BL_FORT_FAB_ARG_3D(qm_pass),
     BL_FORT_FAB_ARG_3D(qp_core),
     BL_FORT_FAB_ARG_3D(qp_pass));

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
     BL_FORT_FAB_ARG_3D(q_core),
     BL_FORT_FAB_ARG_3D(q_pass),
     BL_FORT_FAB_ARG_3D(q_core_bar),
     BL_FORT_FAB_ARG_3D(q_pass_bar),
     BL_FORT_FAB_ARG_3D(qaux),
     BL_FORT_FAB_ARG_3D(qaux_bar),
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

#ifndef AMREX_USE_CUDA
  void ca_sdc_update_advection_o2
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     const int* m_start);

  void ca_sdc_update_advection_o4
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     BL_FORT_FAB_ARG_3D(A_2_old),
     const int* m_start);

#ifdef REACTIONS
  void ca_sdc_update_o2
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     BL_FORT_FAB_ARG_3D(R_0_old),
     BL_FORT_FAB_ARG_3D(R_1_old),
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

  void ca_sdc_compute_C4(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(A_m),
                         const BL_FORT_FAB_ARG_3D(A_0),
                         const BL_FORT_FAB_ARG_3D(A_1),
                         const BL_FORT_FAB_ARG_3D(A_2),
                         const BL_FORT_FAB_ARG_3D(R_0),
                         const BL_FORT_FAB_ARG_3D(R_1),
                         const BL_FORT_FAB_ARG_3D(R_2),
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
