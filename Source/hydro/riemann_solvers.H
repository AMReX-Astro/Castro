#ifndef riemann_solvers_H
#define riemann_solvers_H

#include "Castro_util.H"

///
/// A simple HLL Riemann solver for pure hydrodynamics.  This takes just a
/// single interface's data and returns the HLL flux
///
/// @param ql     the left interface state
/// @param qr     the right interface state
/// @param cl     sound speed on the left interface
/// @param cr     sound speed on the right interface
/// @param idir   coordinate direction for the solve (0 = x, 1 = y, 2 = z)
/// @param coord  geometry type (0 = Cartesian, 1 = axisymmetric, 2 = spherical)
/// @param f      the HLL fluxes
///
AMREX_GPU_HOST_DEVICE
void
HLL(const Real* ql, const Real* qr,
    const Real cl, const Real cr,
    const int idir, const int coord,
    Real* flux_hll) {

  // This is the HLLE solver.  We should apply it to zone averages
  // (not reconstructed states) at an interface in the presence of
  // shocks to avoid the odd-even decoupling / carbuncle phenomenon.
  //
  // See: Einfeldt, B.  et al. 1991, JCP, 92, 273
  //      Einfeldt, B. 1988, SIAM J NA, 25, 294


  constexpr Real small_hll = 1.e-10_rt;

  int ivel, ivelt, iveltt;
  int imom, imomt, imomtt;

  if (idir == 0) {
    ivel = QU;
    ivelt = QV;
    iveltt = QW;

    imom = UMX;
    imomt = UMY;
    imomtt = UMZ;

  } else if (idir == 1) {
    ivel = QV;
    ivelt = QU;
    iveltt = QW;

    imom = UMY;
    imomt = UMX;
    imomtt = UMZ;

  } else {
    ivel = QW;
    ivelt = QU;
    iveltt = QV;

    imom = UMZ;
    imomt = UMX;
    imomtt = UMY;
  }

  Real rhol_sqrt = std::sqrt(ql[QRHO]);
  Real rhor_sqrt = std::sqrt(qr[QRHO]);

  Real rhod = 1.0_rt/(rhol_sqrt + rhor_sqrt);


  // compute the average sound speed. This uses an approximation from
  // E88, eq. 5.6, 5.7 that assumes gamma falls between 1
  // and 5/3
  Real cavg = std::sqrt( (rhol_sqrt*cl*cl + rhor_sqrt*cr*cr)*rhod +
                         0.5_rt*rhol_sqrt*rhor_sqrt*rhod*rhod*std::pow(qr[ivel] - ql[ivel], 2));


  // Roe eigenvalues (E91, eq. 5.3b)
  Real uavg = (rhol_sqrt*ql[ivel] + rhor_sqrt*qr[ivel])*rhod;

  Real a1 = uavg - cavg;
  Real a4 = uavg + cavg;


  // signal speeds (E91, eq. 4.5)
  Real bl = amrex::min(a1, ql[ivel] - cl);
  Real br = amrex::max(a4, qr[ivel] + cr);

  Real bm = amrex::min(0.0_rt, bl);
  Real bp = amrex::max(0.0_rt, br);

  Real bd = bp - bm;

  if (std::abs(bd) < small_hll*amrex::max(std::abs(bm), std::abs(bp))) return;

  // we'll overwrite the passed in flux with the HLL flux

  bd = 1.0_rt/bd;

  // compute the fluxes according to E91, eq. 4.4b -- note that the
  // min/max above picks the correct flux if we are not in the star
  // region

  // density flux
  Real fl_tmp = ql[QRHO]*ql[ivel];
  Real fr_tmp = qr[QRHO]*qr[ivel];

  flux_hll[URHO] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO] - ql[QRHO]);

  // normal momentum flux.  Note for 1-d and 2-d non cartesian
  // r-coordinate, we leave off the pressure term and handle that
  // separately in the update, to accommodate different geometries
  fl_tmp = ql[QRHO]*ql[ivel]*ql[ivel];
  fr_tmp = qr[QRHO]*qr[ivel]*qr[ivel];
  if (mom_flux_has_p(idir, idir, coord)) {
    fl_tmp = fl_tmp + ql[QPRES];
    fr_tmp = fr_tmp + qr[QPRES];
  }

  flux_hll[imom] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[ivel] - ql[QRHO]*ql[ivel]);

  // transverse momentum flux
  fl_tmp = ql[QRHO]*ql[ivel]*ql[ivelt];
  fr_tmp = qr[QRHO]*qr[ivel]*qr[ivelt];

  flux_hll[imomt] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[ivelt] - ql[QRHO]*ql[ivelt]);


  fl_tmp = ql[QRHO]*ql[ivel]*ql[iveltt];
  fr_tmp = qr[QRHO]*qr[ivel]*qr[iveltt];

  flux_hll[imomtt] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[iveltt] - ql[QRHO]*ql[iveltt]);

  // total energy flux
  Real rhoEl = ql[QREINT] + 0.5_rt*ql[QRHO]*(ql[ivel]*ql[ivel] + ql[ivelt]*ql[ivelt] + ql[iveltt]*ql[iveltt]);
  fl_tmp = ql[ivel]*(rhoEl + ql[QPRES]);

  Real rhoEr = qr[QREINT] + 0.5_rt*qr[QRHO]*(qr[ivel]*qr[ivel] + qr[ivelt]*qr[ivelt] + qr[iveltt]*qr[iveltt]);
  fr_tmp = qr[ivel]*(rhoEr + qr[QPRES]);

  flux_hll[UEDEN] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(rhoEr - rhoEl);


  // eint flux
  fl_tmp = ql[QREINT]*ql[ivel];
  fr_tmp = qr[QREINT]*qr[ivel];

  flux_hll[UEINT] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QREINT] - ql[QREINT]);


  // passively-advected scalar fluxes
  for (int ipassive = 0; ipassive < npassive; ipassive++) {
    int n  = upassmap(ipassive);
    int nqs = qpassmap(ipassive);

    fl_tmp = ql[QRHO]*ql[nqs]*ql[ivel];
    fr_tmp = qr[QRHO]*qr[nqs]*qr[ivel];

    flux_hll[n] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[nqs] - ql[QRHO]*ql[nqs]);
  }
}

#endif
