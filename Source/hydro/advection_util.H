#ifndef advection_util_H
#define advection_util_H

#include <Castro_util.H>

#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif

#ifdef RADIATION
#include <rad_util.H>
#endif

namespace hydro {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void src_to_prim(int i, int j, int k, const Real dt,
                 Array4<Real const> const& U_arr,
                 Array4<Real const> const& q_arr,
                 Array4<Real const> const& old_src,
#ifndef TRUE_SDC
                 Array4<Real const> const& src_corr,
#endif
                 Array4<Real> const& srcQ)
{

    amrex::ignore_unused(dt);

    for (int n = 0; n < NQSRC; ++n) {
      srcQ(i,j,k,n) = 0.0_rt;
    }

    // the conserved source may have a predictor that time-centers it

    Real srcU[NSRC] = {0.0_rt};

    for (int n = 0; n < NSRC; n++) {

#ifndef TRUE_SDC
        if (time_integration_method == CornerTransportUpwind && source_term_predictor == 1) {
            if (n == UMX || n == UMY || n == UMZ) {
                srcU[n] += 0.5 * dt * src_corr(i,j,k,n);
            }
        } else if (time_integration_method == SimplifiedSpectralDeferredCorrections  && source_term_predictor == 1) {
            srcU[n] += src_corr(i,j,k,n);
        }
#endif

        srcU[n] += old_src(i,j,k,n);
    }

    Real rhoinv = 1.0_rt / q_arr(i,j,k,QRHO);

    // get the needed derivatives
    eos_rep_t eos_state;
    eos_state.T = q_arr(i,j,k,QTEMP);
    eos_state.rho = q_arr(i,j,k,QRHO);
    eos_state.e = q_arr(i,j,k,QREINT) * rhoinv;
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n]  = U_arr(i,j,k,UFS+n) * rhoinv;
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = U_arr(i,j,k,UFX+n) * rhoinv;
    }
#endif

    eos(eos_input_re, eos_state);

    srcQ(i,j,k,QRHO) = srcU[URHO];
    srcQ(i,j,k,QU) = (srcU[UMX] - q_arr(i,j,k,QU) * srcQ(i,j,k,QRHO)) * rhoinv;
    srcQ(i,j,k,QV) = (srcU[UMY] - q_arr(i,j,k,QV) * srcQ(i,j,k,QRHO)) * rhoinv;
    srcQ(i,j,k,QW) = (srcU[UMZ] - q_arr(i,j,k,QW) * srcQ(i,j,k,QRHO)) * rhoinv;
    srcQ(i,j,k,QREINT) = srcU[UEINT];
    srcQ(i,j,k,QPRES ) = eos_state.dpde *
      (srcQ(i,j,k,QREINT) - q_arr(i,j,k,QREINT) * srcQ(i,j,k,QRHO)*rhoinv) *
      rhoinv + eos_state.dpdr_e * srcQ(i,j,k,QRHO);

#ifdef CONS_SPECIES_HAVE_SOURCES
    // Do the passive variables, i.e. advected, species, and auxiliary

    for (int ipassive = 0; ipassive < npassive; ipassive++) {
        const int nc = upassmap(ipassive);
        const int n = qpassmap(ipassive);
        srcQ(i,j,k,n) = (srcU[nc] - q_arr(i,j,k,n) * srcQ(i,j,k,QRHO)) * rhoinv;
    }
#endif
}

///
/// Convert conservative variables to primitive variables. The method is
/// templated so that q/qaux can either be either Array4 or Array1D.
///
/// @param i         x index
/// @param j         y index
/// @param k         z index
/// @param uin       input conserved state
/// @param Bx        x-component of magentic field (if USE_MHD=TRUE)
/// @param By        y-component of magentic field (if USE_MHD=TRUE)
/// @param Bz        z-component of magentic field (if USE_MHD=TRUE)
/// @param Erin      radiation energy (if USE_RAD=TRUE)
/// @param lam       radiation flux limiter (if USE_RAD=TRUE)
/// @param q         output primitive state
/// @param qaux      output auxiliary quantities
///

template<class T, class U>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void conservative_to_primitive(int i, int j, int k,
                               Array4<Real const> const& uin,
#ifdef MHD
                               Array4<Real const> const& Bx,
                               Array4<Real const> const& By,
                               Array4<Real const> const& Bz,
#endif
#ifdef RADIATION
                               Array4<Real const> const& Erin,
                               Array4<Real const> const& lam,
#endif
                               T& q, U& qaux, const bool fill_passives)
{
#ifndef AMREX_USE_GPU
    if (uin(i,j,k,URHO) <= 0.0_rt) {
        std::cout << std::endl;
        std::cout << ">>> Error: conservative_to_primitive " << i << " " << j << " " << k << std::endl;
        std::cout << ">>> ... negative density " << uin(i,j,k,URHO) << std::endl;
        amrex::Error("Error:: hydro::conservative_to_primitive()");
    } else if (uin(i,j,k,URHO) < castro::small_dens) {
        std::cout << std::endl;
        std::cout << ">>> Error: conservative_to_primitive " << i << " " << j << " " << k << std::endl;
        std::cout << ">>> ... small density " << uin(i,j,k,URHO) << std::endl;
        amrex::Error("Error:: hydro::conservative_to_primitive()");
    }
#endif

    q(QRHO) = uin(i,j,k,URHO);
    Real rhoinv = 1.0_rt / q(QRHO);

    q(QU) = uin(i,j,k,UMX) * rhoinv;
    q(QV) = uin(i,j,k,UMY) * rhoinv;
    q(QW) = uin(i,j,k,UMZ) * rhoinv;

#ifdef MHD
    q(QMAGX) = 0.5_rt * (Bx(i+1,j,k) + Bx(i,j,k));
    q(QMAGY) = 0.5_rt * (By(i,j+1,k) + By(i,j,k));
    q(QMAGZ) = 0.5_rt * (Bz(i,j,k+1) + Bz(i,j,k));
#endif

    // Get the internal energy, which we'll use for
    // determining the pressure.  We use a dual energy
    // formalism. If (E - K) < eta1 and eta1 is suitably
    // small, then we risk serious numerical truncation error
    // in the internal energy.  Therefore we'll use the result
    // of the separately updated internal energy equation.
    // Otherwise, we'll set e = E - K.

    Real kineng = 0.5_rt * q(QRHO) * (q(QU) * q(QU) + q(QV) * q(QV) + q(QW) * q(QW));

    if ((uin(i,j,k,UEDEN) - kineng) > castro::dual_energy_eta1*uin(i,j,k,UEDEN)) {
        q(QREINT) = (uin(i,j,k,UEDEN) - kineng) * rhoinv;
    } else {
        q(QREINT) = uin(i,j,k,UEINT) * rhoinv;
    }

    q(QTEMP) = uin(i,j,k,UTEMP);
#ifdef RADIATION
    for (int g = 0; g < NGROUPS; g++) {
        q(QRAD+g) = Erin(i,j,k,g);
    }
#endif

    if (fill_passives) {
        // Load passively advected quantities into q
        for (int ipassive = 0; ipassive < npassive; ipassive++) {
            int n  = upassmap(ipassive);
            int iq = qpassmap(ipassive);
            q(iq) = uin(i,j,k,n) * rhoinv;
        }
    }

    // get gamc, p, T, c, csml using q state
    eos_rep_t eos_state;
    eos_state.T = q(QTEMP);
    eos_state.rho = q(QRHO);
    eos_state.e = q(QREINT);
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n]  = uin(i,j,k,UFS+n) * rhoinv;
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = uin(i,j,k,UFX+n) * rhoinv;
    }
#endif

    eos(eos_input_re, eos_state);

    q(QTEMP) = eos_state.T;
    q(QREINT) = eos_state.e * q(QRHO);
    q(QPRES) = eos_state.p;
#ifdef TRUE_SDC
    q(QGC) = eos_state.gam1;
#endif

#ifdef MHD
    q(QPTOT) = q(QPRES) +
               0.5_rt * (q(QMAGX) * q(QMAGX) +
                         q(QMAGY) * q(QMAGY) +
                         q(QMAGZ) * q(QMAGZ));
#endif

#ifdef RADIATION
    qaux(QGAMCG) = eos_state.gam1;
    qaux(QCG) = eos_state.cs;

    Real lams[NGROUPS];
    for (int g = 0; g < NGROUPS; g++) {
        lams[g] = lam(i,j,k,g);
    }
    Real qs[NQ];
    for (int n = 0; n < NQ; n++) {
        qs[n] = q(n);
    }
    Real ptot;
    Real ctot;
    Real gamc_tot;
    compute_ptot_ctot(lams, qs,
                      qaux(QCG),
                      ptot, ctot, gamc_tot);

    q(QPTOT) = ptot;

    qaux(QC) = ctot;
    qaux(QGAMC) = gamc_tot;

    q(QREITOT) = q(QREINT);
    for (int g = 0; g < NGROUPS; g++) {
        qaux(QLAMS+g) = lam(i,j,k,g);
        q(QREITOT) += q(QRAD+g);
    }
#else
    qaux(QGAMC) = eos_state.gam1;
    qaux(QC) = eos_state.cs;
#endif
}

}

#endif
