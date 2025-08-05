#ifndef riemann_solvers_H
#define riemann_solvers_H

#include <Castro_util.H>
#include <riemann_type.H>
#include <riemann_2shock_solvers.H>
#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif
#ifdef RADIATION
#include <Radiation.H>
#include <rad_util.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
compute_flux_q(const int i, const int j, const int k, const int idir,
               const GeometryData& geomdata,
               const RiemannState& qint,
               Array4<Real> const& F,
#ifdef RADIATION
               Array4<Real> const& rF,
#endif
               Array4<Real> const& qgdnv, const bool store_full_state) {

    // given a primitive state, compute the flux in direction idir
    //

    int imx1, imx2, imx3;

    auto coord = geomdata.Coord();
    auto mom_check = mom_flux_has_p(idir, idir, coord);

    if (idir == 0) {
        imx1 = UMX;
        imx2 = UMY;
        imx3 = UMZ;

    } else if (idir == 1) {
        imx1 = UMY;
        imx2 = UMX;
        imx3 = UMZ;

    } else {
        imx1 = UMZ;
        imx2 = UMX;
        imx3 = UMY;
    }

    // Compute fluxes, order as conserved state (not q)
    F(i,j,k,URHO) = qint.rho * qint.un;

    F(i,j,k,imx1) = F(i,j,k,URHO) * qint.un;
    F(i,j,k,imx2) = F(i,j,k,URHO) * qint.ut;
    F(i,j,k,imx3) = F(i,j,k,URHO) * qint.utt;

#ifndef RADIATION
    if (mom_check) {
      F(i,j,k,imx1) += qint.p;
    }

    Real rhoetot = qint.rhoe + 0.5_rt * qint.rho *
        (qint.un * qint.un +
         qint.ut * qint.ut +
         qint.utt * qint.utt);

    F(i,j,k,UEDEN) = qint.un * (rhoetot + qint.p);
    F(i,j,k,UEINT) = qint.un * qint.rhoe;
#else
    if (mom_check) {
      F(i,j,k,imx1) += qint.p_g;
    }

    Real rhoetot = qint.rhoe_g + 0.5_rt * qint.rho *
        (qint.un * qint.un +
         qint.ut * qint.ut +
         qint.utt * qint.utt);

    F(i,j,k,UEDEN) = qint.un * (rhoetot + qint.p_g);
    F(i,j,k,UEINT) = qint.un * qint.rhoe_g;
#endif

    F(i,j,k,UTEMP) = 0.0;
#ifdef SHOCK_VAR
    F(i,j,k,USHK) = 0.0;
#endif

#ifdef NSE_NET
    F(i,j,k,UMUP) = 0.0;
    F(i,j,k,UMUN) = 0.0;
#endif

#ifdef RADIATION
    if (radiation::fspace_advection_type == 1) {
        for (int g = 0; g < NGROUPS; g++) {
            Real eddf = Edd_factor(qint.lam[g]);
            Real f1 = 0.5e0_rt*(1.0_rt-eddf);
            rF(i,j,k,g) = (1.0_rt + f1) * qint.er[g] * qint.un;
        }
    } else {
        // type 2
        for (int g = 0; g < NGROUPS; g++) {
            rF(i,j,k,g) = qint.er[g] * qint.un;
        }
    }
#endif

#ifdef HYBRID_MOMENTUM
    // the hybrid routine uses the Godunov indices, not the full NQ state
    GpuArray<Real, NGDNV> qgdnv_zone;
    qgdnv_zone[GDRHO] = qint.rho;
    if (idir == 0) {
        qgdnv_zone[GDU] = qint.un;
        qgdnv_zone[GDV] = qint.ut;
        qgdnv_zone[GDW] = qint.utt;
    } else if (idir == 1) {
        qgdnv_zone[GDV] = qint.un;
        qgdnv_zone[GDU] = qint.ut;
        qgdnv_zone[GDW] = qint.utt;
    } else {
        qgdnv_zone[GDW] = qint.un;
        qgdnv_zone[GDU] = qint.ut;
        qgdnv_zone[GDV] = qint.utt;
    }
    qgdnv_zone[GDPRES] = qint.p;
#ifdef RADIATION
    for (int g = 0; g < NGROUPS; g++) {
        qgdnv_zone[GDLAMS+g] = qint.lam[g];
        qgdnv_zone[GDERADS+g] = qint.er[g];
    }
#endif
    GpuArray<Real, NUM_STATE> F_zone;
    for (int n = 0; n < NUM_STATE; n++) {
        F_zone[n] = F(i,j,k,n);
    }
    compute_hybrid_flux(qgdnv_zone, geomdata, idir, i, j, k, F_zone);
    for (int n = 0; n < NUM_STATE; n++) {
        F(i,j,k,n) = F_zone[n];
    }
#endif

    // now we store the interface state in the full 4-d array.  We
    // either store just the NGDNV subset or the whole NQ hydro state,
    // depending on the value of store_full_state

    // note: we don't handle the passives here

    if (store_full_state) {

        // this is the primitive variable state on the interface

        qgdnv(i,j,k,QRHO) = qint.rho;
        if (idir == 0) {
            qgdnv(i,j,k,QU) = qint.un;
            qgdnv(i,j,k,QV) = qint.ut;
            qgdnv(i,j,k,QW) = qint.utt;
        } else if (idir == 1) {
            qgdnv(i,j,k,QV) = qint.un;
            qgdnv(i,j,k,QU) = qint.ut;
            qgdnv(i,j,k,QW) = qint.utt;
        } else {
            qgdnv(i,j,k,QW) = qint.un;
            qgdnv(i,j,k,QU) = qint.ut;
            qgdnv(i,j,k,QV) = qint.utt;
        }

        qgdnv(i,j,k,QTEMP) = 0.0_rt;

#ifdef RADIATION
        qgdnv(i,j,k,QPRES) = qint.p_g;
        qgdnv(i,j,k,QREINT) = qint.rhoe_g;
        qgdnv(i,j,k,QPTOT) = qint.p;
        qgdnv(i,j,k,QREITOT) = qint.rhoe;
        for (int g = 0; g < NGROUPS; g++) {
            qgdnv(i,j,k,QRAD) = qint.er[g];
        }
#else
        qgdnv(i,j,k,QPRES) = qint.p;
        qgdnv(i,j,k,QREINT) = qint.rhoe;
#endif

    } else {

        // this is the reduced "Godunov" state on the interface,
        // with NGDNV components

#ifdef HYBRID_MOMENTUM
        qgdnv(i,j,k,GDRHO) = qint.rho;
#endif
        if (idir == 0) {
            qgdnv(i,j,k,GDU) = qint.un;
            qgdnv(i,j,k,GDV) = qint.ut;
            qgdnv(i,j,k,GDW) = qint.utt;
        } else if (idir == 1) {
            qgdnv(i,j,k,GDV) = qint.un;
            qgdnv(i,j,k,GDU) = qint.ut;
            qgdnv(i,j,k,GDW) = qint.utt;
        } else {
            qgdnv(i,j,k,GDW) = qint.un;
            qgdnv(i,j,k,GDU) = qint.ut;
            qgdnv(i,j,k,GDV) = qint.utt;
        }
#ifndef RADIATION
        qgdnv(i,j,k,GDPRES) = qint.p;
#else
        qgdnv(i,j,k,GDPRES) = qint.p_g;
        for (int g = 0; g < NGROUPS; g++) {
            qgdnv(i,j,k,GDLAMS+g) = qint.lam[g];
            qgdnv(i,j,k,GDERADS+g) = qint.er[g];
        }
#endif

    }

}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
riemann_state(const int i, const int j, const int k, const int idir,
              Array4<Real> const& qm,
              Array4<Real> const& qp,
              Array4<Real const> const& qaux_arr,
              RiemannState& qint,
              const bool special_bnd_lo, const bool special_bnd_hi,
              GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi) {

  // just compute the hydrodynamic state on the interfaces
  // don't compute the fluxes

  // note: bx is not necessarily the limits of the valid (no ghost
  // cells) domain, but could be hi+1 in some dimensions.  We rely on
  // the caller to specify the interfaces over which to solve the
  // Riemann problems


  if (ppm_temp_fix == 2) {
      // recompute the thermodynamics on the interface to make it
      // all consistent

      // we want to take the edge states of rho, e, and X, and get
      // new values for p on the edges that are
      // thermodynamically consistent.

      eos_rep_t eos_state;

      // this is an initial guess for iterations, since we
      // can't be certain what temp is on interfaces
      eos_state.T = T_guess;

      // minus state
      eos_state.rho = qm(i,j,k,QRHO);
      eos_state.e = qm(i,j,k,QREINT)/qm(i,j,k,QRHO);
      for (int n = 0; n < NumSpec; n++) {
          eos_state.xn[n] = qm(i,j,k,QFS+n);
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
          eos_state.aux[n] = qm(i,j,k,QFX+n);
      }
#endif

      eos(eos_input_re, eos_state);

      qm(i,j,k,QREINT) = eos_state.e * eos_state.rho;
      qm(i,j,k,QPRES) = eos_state.p;

      // plus state
      eos_state.rho = qp(i,j,k,QRHO);
      eos_state.e = qp(i,j,k,QREINT)/qp(i,j,k,QRHO);
      for (int n = 0; n < NumSpec; n++) {
          eos_state.xn[n] = qp(i,j,k,QFS+n);
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
          eos_state.aux[n] = qp(i,j,k,QFX+n);
      }
#endif

      eos(eos_input_re, eos_state);

      qp(i,j,k,QREINT) = eos_state.e * eos_state.rho;
      qp(i,j,k,QPRES) = eos_state.p;
  }


  RiemannState ql;
  RiemannState qr;
  RiemannAux raux;


  load_input_states(i, j, k, idir,
                    qm, qp, qaux_arr,
                    ql, qr, raux);

  // deal with hard walls
  raux.bnd_fac = 1.0_rt;

  if (idir == 0) {
      if ((i == domlo[0] && special_bnd_lo) ||
          (i == domhi[0]+1 && special_bnd_hi)) {
          raux.bnd_fac = 0.0_rt;
      }

  } else if (idir == 1) {
      if ((j == domlo[1] && special_bnd_lo) ||
          (j == domhi[1]+1 && special_bnd_hi)) {
          raux.bnd_fac = 0.0_rt;
      }

  } else {
      if ((k == domlo[2] && special_bnd_lo) ||
          (k == domhi[2]+1 && special_bnd_hi)) {
          raux.bnd_fac = 0.0_rt;
      }
  }


  // Solve Riemann problem
  if (riemann_solver == 0) {
      // Colella, Glaz, & Ferguson solver

      TwoShock::riemannus(ql, qr, raux, qint);

  } else if (riemann_solver == 1) {
      // Colella & Glaz solver

#ifndef RADIATION
      TwoShock::riemanncg(ql, qr, raux, qint);
#endif

#ifndef AMREX_USE_GPU
  } else {
      amrex::Error("ERROR: invalid value of riemann_solver");
#endif
  }


}


#endif
