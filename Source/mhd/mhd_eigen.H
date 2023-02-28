#ifndef CASTRO_MHD_EIGEN_H
#define CASTRO_MHD_EIGEN_H

#include <mhd_util.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
evals(Array1D<Real, 0, NEIGN-1>& lam,
      const Real as_in,
      Array1D<Real, 0, NQ-1>& Q,
      const int dir) {

  // The characteristic speeds of the system

  // Alfven

  Real ca  = (Q(QMAGX)*Q(QMAGX) + Q(QMAGY)*Q(QMAGY) + Q(QMAGZ)*Q(QMAGZ)) / Q(QRHO);

  int QUN;
  Real cad; // Alfven
  Real cs;  // slow
  Real cf;  // fast

  Real as = as_in * as_in;

  if (dir == 0) {
    QUN = QU;
    cad = (Q(QMAGX)*Q(QMAGX)) / Q(QRHO);
    cs = 0.5_rt * ((as + ca) - std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*cad));
    cf = 0.5_rt * ((as + ca) + std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*cad));

  } else if (dir == 1) {
    QUN = QV;
    cad = (Q(QMAGY)*Q(QMAGY)) / Q(QRHO);
    cs = 0.5_rt * ((as + ca) - std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*cad));
    cf = 0.5_rt * ((as + ca) + std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*cad));

  } else {
    QUN = QW;
    cad = (Q(QMAGZ)*Q(QMAGZ)) / Q(QRHO);
    cs = 0.5_rt * ((as + ca) - std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*cad));
    cf = 0.5_rt * ((as + ca) + std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*cad));

  }

  // eigenvalues
  lam(0) = Q(QUN) - std::sqrt(cf);
  lam(1) = Q(QUN) - std::sqrt(cad);
  lam(2) = Q(QUN) - std::sqrt(cs);
  lam(3) = Q(QUN);
  lam(4) = Q(QUN) + std::sqrt(cs);
  lam(5) = Q(QUN) + std::sqrt(cad);
  lam(6) = Q(QUN) + std::sqrt(cf);

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
evecx(Array2D<Real, 0, NEIGN-1, 0, NEIGN-1>& leig,
      Array2D<Real, 0, NEIGN-1, 0, NEIGN-1>& reig,
      const Real as_in,
      Array1D<Real, 0, NQ-1>& Q) {

  // returns Leig, where the rows are the left eigenvectors of the characteristic matrix Ax

  // The characteristic speeds of the system

  Real as = as_in * as_in;

  // Alfven

  Real ca = (Q(QMAGX)*Q(QMAGX) + Q(QMAGY)*Q(QMAGY) + Q(QMAGZ)*Q(QMAGZ)) / Q(QRHO);
  Real cax = (Q(QMAGX)*Q(QMAGX)) / Q(QRHO);

  // Slow

  Real csx = 0.5_rt*((as + ca) - std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*cax));

  // Fast

  Real cfx = 0.5_rt*((as + ca) + std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*cax));

  // useful constants

  Real alf;
  Real als;
 
  if (std::abs(cfx-csx) <= 1e-14){
    alf = 1.0_rt;
    als = 0.0_rt;
  }
  
  if (as - csx < 0.0) {
    alf = 0.0_rt;
  } else {
    alf = std::sqrt((as - csx)/(cfx - csx));
  }

  if (cfx - as < 0.0) {
    als = 0.0_rt;
  } else {
    als = std::sqrt((cfx - as)/(cfx - csx));
  }
  
  

  Real bety;
  Real betz;

  if (std::abs(Q(QMAGY)) <= 1.e-14_rt && std::abs(Q(QMAGZ)) <= 1.e-14_rt) {
    bety = 1.0_rt / std::sqrt(2.0_rt);
    betz = bety;

  } else {
    bety = Q(QMAGY) / (std::sqrt(Q(QMAGY)*Q(QMAGY) + Q(QMAGZ)*Q(QMAGZ)));
    betz = Q(QMAGZ) / (std::sqrt(Q(QMAGY)*Q(QMAGY) + Q(QMAGZ)*Q(QMAGZ)));
  }

  Real cff = std::sqrt(cfx) * alf;
  Real css = std::sqrt(csx) * als;

  Real S = std::copysign(1.0_rt, Q(QMAGX));

  Real Qf = std::sqrt(cfx) * alf * S;
  Real Qs = std::sqrt(csx) * als * S;

  Real N = 0.5_rt / as;

  Real AAf = std::sqrt(as) * alf * std::sqrt(Q(QRHO));
  Real AAs = std::sqrt(as) * als * std::sqrt(Q(QRHO));


  // u - cf wave
  leig(0, IEIGN_RHO) = 0.0;
  leig(0, IEIGN_U) = -N*cff;
  leig(0, IEIGN_V) = N*Qs*bety;
  leig(0, IEIGN_W) = N*Qs*betz;
  leig(0, IEIGN_P) = N*alf/Q(QRHO);
  leig(0, IEIGN_BT) = N*AAs*bety/Q(QRHO);
  leig(0, IEIGN_BTT) = N*AAs*betz/Q(QRHO);

  // u - cax
  leig(1, IEIGN_RHO) = 0.0;
  leig(1, IEIGN_U) = 0.0;
  leig(1, IEIGN_V) = -0.5_rt*betz;
  leig(1, IEIGN_W) = 0.5_rt*bety;
  leig(1, IEIGN_P) = 0.0;
  leig(1, IEIGN_BT) = -0.5_rt*betz*S/std::sqrt(Q(QRHO));
  leig(1, IEIGN_BTT) = 0.5_rt*bety*S/std::sqrt(Q(QRHO));

  // u - cs
  leig(2, IEIGN_RHO) = 0.0;
  leig(2, IEIGN_U) = -N*css;
  leig(2, IEIGN_V) = -N*Qf*bety;
  leig(2, IEIGN_W) = -N*Qf*betz;
  leig(2, IEIGN_P) = N*als/Q(QRHO);
  leig(2, IEIGN_BT) = -N*AAf*bety/Q(QRHO);
  leig(2, IEIGN_BTT) = -N*AAf*betz/Q(QRHO);

  // u
  leig(3, IEIGN_RHO) = 1.0;
  leig(3, IEIGN_U) = 0.0;
  leig(3, IEIGN_V) = 0.0;
  leig(3, IEIGN_W) = 0.0;
  leig(3, IEIGN_P) = -1.0_rt/as;
  leig(3, IEIGN_BT) = 0.0;
  leig(3, IEIGN_BTT) = 0.0;

  // u + cs
  leig(4, IEIGN_RHO) = 0.0;
  leig(4, IEIGN_U) = N*css;
  leig(4, IEIGN_V) = N*Qf*bety;
  leig(4, IEIGN_W) = N*Qf*betz;
  leig(4, IEIGN_P) = N*als/Q(QRHO);
  leig(4, IEIGN_BT) = -N*AAf*bety/Q(QRHO);
  leig(4, IEIGN_BTT) = -N*AAf*betz/Q(QRHO);

  // u + cax
  leig(5, IEIGN_RHO) = 0.0;
  leig(5, IEIGN_U) = 0.0;
  leig(5, IEIGN_V) = 0.5_rt*betz;
  leig(5, IEIGN_W) = -0.5_rt*bety;
  leig(5, IEIGN_P) = 0.0;
  leig(5, IEIGN_BT) = -0.5_rt*betz*S/std::sqrt(Q(QRHO));
  leig(5, IEIGN_BTT) = 0.5_rt*bety*S/std::sqrt(Q(QRHO));

  // u + cf
  leig(6, IEIGN_RHO) = 0.0;
  leig(6, IEIGN_U) = N*cff;
  leig(6, IEIGN_V) = -N*Qs*bety;
  leig(6, IEIGN_W) = -N*Qs*betz;
  leig(6, IEIGN_P) = N*alf/Q(QRHO);
  leig(6, IEIGN_BT) = N*AAs*bety/Q(QRHO);
  leig(6, IEIGN_BTT) = N*AAs*betz/Q(QRHO);

  // u - cf ; u - Cax ; u - cs ; u ; u + cs ; u + Cax ; u + cf

  reig(IEIGN_RHO, 0) = Q(QRHO)*alf;
  reig(IEIGN_RHO, 1) = 0.0;
  reig(IEIGN_RHO, 2) = Q(QRHO)*als;
  reig(IEIGN_RHO, 3) = 1.0_rt;
  reig(IEIGN_RHO, 4) = Q(QRHO)*als;
  reig(IEIGN_RHO, 5) = 0.0;
  reig(IEIGN_RHO, 6) = Q(QRHO)*alf;

  reig(IEIGN_U, 0) = -cff;
  reig(IEIGN_U, 1) = 0.0;
  reig(IEIGN_U, 2) = -css;
  reig(IEIGN_U, 3) = 0.0;
  reig(IEIGN_U, 4) = css;
  reig(IEIGN_U, 5) = 0.0;
  reig(IEIGN_U, 6) = cff;

  reig(IEIGN_V, 0) = Qs*bety;
  reig(IEIGN_V, 1) = -betz;
  reig(IEIGN_V, 2) = -Qf*bety;
  reig(IEIGN_V, 3) = 0.0;
  reig(IEIGN_V, 4) = Qf*bety;
  reig(IEIGN_V, 5) = betz;
  reig(IEIGN_V, 6) = -Qs*bety;

  reig(IEIGN_W, 0) = Qs*betz;
  reig(IEIGN_W, 1) = bety;
  reig(IEIGN_W, 2) = -Qf*betz;
  reig(IEIGN_W, 3) = 0.0;
  reig(IEIGN_W, 4) = Qf*betz;
  reig(IEIGN_W, 5) = -bety;
  reig(IEIGN_W, 6) = -Qs*betz;

  reig(IEIGN_P, 0) = Q(QRHO)*as*alf;
  reig(IEIGN_P, 1) = 0.0;
  reig(IEIGN_P, 2) = Q(QRHO)*as*als;
  reig(IEIGN_P, 3) = 0.0;
  reig(IEIGN_P, 4) = Q(QRHO)*as*als;
  reig(IEIGN_P, 5) = 0.0;
  reig(IEIGN_P, 6) = Q(QRHO)*as*alf;

  reig(IEIGN_BT, 0) = AAs*bety;
  reig(IEIGN_BT, 1) = -betz*S*std::sqrt(Q(QRHO));
  reig(IEIGN_BT, 2) = -AAf*bety;
  reig(IEIGN_BT, 3) = 0.0;
  reig(IEIGN_BT, 4) = -AAf*bety;
  reig(IEIGN_BT, 5) = -betz*S*std::sqrt(Q(QRHO));
  reig(IEIGN_BT, 6) = AAs*bety;

  reig(IEIGN_BTT, 0) = AAs*betz;
  reig(IEIGN_BTT, 1) = bety*S*std::sqrt(Q(QRHO));
  reig(IEIGN_BTT, 2) = -AAf*betz;
  reig(IEIGN_BTT, 3) = 0.0;
  reig(IEIGN_BTT, 4) = -AAf*betz;
  reig(IEIGN_BTT, 5) = bety*S*std::sqrt(Q(QRHO));
  reig(IEIGN_BTT, 6) = AAs*betz;

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
evecy(Array2D<Real, 0, NEIGN-1, 0, NEIGN-1>& leig,
      Array2D<Real, 0, NEIGN-1, 0, NEIGN-1>& reig,
      const Real as_in,
      Array1D<Real, 0, NQ-1>& Q) {


  // returns Leig, where the rows are the left eigenvectors of the characteristic matrix Ay

  // The characteristic speeds of the system

  Real as = as_in * as_in;

  // Alfven

  Real ca = (Q(QMAGX)*Q(QMAGX) + Q(QMAGY)*Q(QMAGY) + Q(QMAGZ)*Q(QMAGZ)) / Q(QRHO);
  Real cay = (Q(QMAGY)*Q(QMAGY)) / Q(QRHO);

  // Slow

  Real csy = 0.5_rt*((as + ca) - std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*cay));

  // Fast

  Real cfy = 0.5_rt*((as + ca) + std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*cay));

  // useful constants

  Real alf;
  Real als;
  
  if (std::abs(cfy-csy) <= 1e-14){
    alf = 1.0_rt;
    als = 0.0_rt;
  }

  if (as - csy < 0.0) {
    alf = 0.0_rt;
  } else {
    alf = std::sqrt((as - csy)/(cfy - csy));
  }

  if (cfy - as < 0.0) {
    als = 0.0_rt;
  } else {
    als = std::sqrt((cfy - as)/(cfy - csy));
  }

  Real betx;
  Real betz;

  if(std::abs(Q(QMAGX)) <= 1.e-14_rt && std::abs(Q(QMAGZ)) <= 1.e-14_rt) {
    betx = 1.0_rt / std::sqrt(2.0_rt);
    betz = betx;

  } else {
    betx = Q(QMAGX) / std::sqrt(Q(QMAGX)*Q(QMAGX) + Q(QMAGZ)*Q(QMAGZ));
    betz = Q(QMAGZ) / std::sqrt(Q(QMAGX)*Q(QMAGX) + Q(QMAGZ)*Q(QMAGZ));
  }

  Real cff = std::sqrt(cfy) * alf;
  Real css = std::sqrt(csy) * als;

  Real S = std::copysign(1.0_rt, Q(QMAGY));

  Real Qf = std::sqrt(cfy) * alf * S;
  Real Qs = std::sqrt(csy) * als * S;

  Real N = 0.5_rt / as;

  Real AAf = std::sqrt(as) * alf * std::sqrt(Q(QRHO));
  Real AAs = std::sqrt(as) * als * std::sqrt(Q(QRHO));


  // v - cf
  leig(0, IEIGN_RHO) = 0.0;
  leig(0, IEIGN_U) = N*Qs*betx;
  leig(0, IEIGN_V) = -N*cff;
  leig(0, IEIGN_W) = N*Qs*betz;
  leig(0, IEIGN_P) = N*alf/Q(QRHO);
  leig(0, IEIGN_BT) = N*AAs*betx/Q(QRHO);
  leig(0, IEIGN_BTT) = N*AAs*betz/Q(QRHO);

  // v - cay
  leig(1, IEIGN_RHO) = 0.0;
  leig(1, IEIGN_U) = -0.5_rt*betz;
  leig(1, IEIGN_V) = 0.0;
  leig(1, IEIGN_W) = 0.5_rt*betx;
  leig(1, IEIGN_P) = 0.0;
  leig(1, IEIGN_BT) = -0.5_rt*betz*S/std::sqrt(Q(QRHO));
  leig(1, IEIGN_BTT) = 0.5_rt*betx*S/std::sqrt(Q(QRHO));

  // v - cs
  leig(2, IEIGN_RHO) = 0.0;
  leig(2, IEIGN_U) = -N*Qf*betx;
  leig(2, IEIGN_V) = -N*css;
  leig(2, IEIGN_W) = -N*Qf*betz;
  leig(2, IEIGN_P) = N*als/Q(QRHO);
  leig(2, IEIGN_BT) = -N*AAf*betx/Q(QRHO);
  leig(2, IEIGN_BTT) = -N*AAf*betz/Q(QRHO);

  // v
  leig(3, IEIGN_RHO) = 1.0_rt;
  leig(3, IEIGN_U) = 0.0;
  leig(3, IEIGN_V) = 0.0;
  leig(3, IEIGN_W) = 0.0;
  leig(3, IEIGN_P) = -1.0_rt/as;
  leig(3, IEIGN_BT) = 0.0;
  leig(3, IEIGN_BTT) = 0.0;

  // v + cs
  leig(4, IEIGN_RHO) = 0.0;
  leig(4, IEIGN_U) = N*Qf*betx;
  leig(4, IEIGN_V) = N*css;
  leig(4, IEIGN_W) = N*Qf*betz;
  leig(4, IEIGN_P) = N*als/Q(QRHO);
  leig(4, IEIGN_BT) = -N*AAf*betx/Q(QRHO);
  leig(4, IEIGN_BTT) = -N*AAf*betz/Q(QRHO);

  // v + cay
  leig(5, IEIGN_RHO) = 0.0;
  leig(5, IEIGN_U) = 0.5_rt*betz;
  leig(5, IEIGN_V) = 0.0;
  leig(5, IEIGN_W) = -0.5_rt*betx;
  leig(5, IEIGN_P) = 0.0;
  leig(5, IEIGN_BT) = -0.5_rt*betz*S/std::sqrt(Q(QRHO));
  leig(5, IEIGN_BTT) = 0.5_rt*betx*S/std::sqrt(Q(QRHO));

  // v + cf
  leig(6, IEIGN_RHO) = 0.0;
  leig(6, IEIGN_U) = -N*Qs*betx;
  leig(6, IEIGN_V) = N*cff;
  leig(6, IEIGN_W) = -N*Qs*betz;
  leig(6, IEIGN_P) = N*alf/Q(QRHO);
  leig(6, IEIGN_BT) = N*AAs*betx/Q(QRHO);
  leig(6, IEIGN_BTT) = N*AAs*betz/Q(QRHO);

  // v - cf ; v - Cay ; v - cs ; v ; v + cs ; v + Cay ; v + cf

  reig(IEIGN_RHO, 0) = Q(QRHO)*alf;
  reig(IEIGN_RHO, 1) = 0.0;
  reig(IEIGN_RHO, 2) = Q(QRHO)*als;
  reig(IEIGN_RHO, 3) = 1.0_rt;
  reig(IEIGN_RHO, 4) = Q(QRHO)*als;
  reig(IEIGN_RHO, 5) = 0.0;
  reig(IEIGN_RHO, 6) = Q(QRHO)*alf;

  reig(IEIGN_V, 0) = -cff;
  reig(IEIGN_V, 1) = 0.0;
  reig(IEIGN_V, 2) = -css;
  reig(IEIGN_V, 3) = 0.0;
  reig(IEIGN_V, 4) = css;
  reig(IEIGN_V, 5) = 0.0;
  reig(IEIGN_V, 6) = cff;

  reig(IEIGN_U, 0) = Qs*betx;
  reig(IEIGN_U, 1) = -betz;
  reig(IEIGN_U, 2) = -Qf*betx;
  reig(IEIGN_U, 3) = 0.0;
  reig(IEIGN_U, 4) = Qf*betx;
  reig(IEIGN_U, 5) = betz;
  reig(IEIGN_U, 6) = -Qs*betx;

  reig(IEIGN_W, 0) = Qs*betz;
  reig(IEIGN_W, 1) = betx;
  reig(IEIGN_W, 2) = -Qf*betz;
  reig(IEIGN_W, 3) = 0.0;
  reig(IEIGN_W, 4) = Qf*betz;
  reig(IEIGN_W, 5) = -betx;
  reig(IEIGN_W, 6) = -Qs*betz;

  reig(IEIGN_P, 0) = Q(QRHO)*as*alf;
  reig(IEIGN_P, 1) = 0.0;
  reig(IEIGN_P, 2) = Q(QRHO)*as*als;
  reig(IEIGN_P, 3) = 0.0;
  reig(IEIGN_P, 4) = Q(QRHO)*as*als;
  reig(IEIGN_P, 5) = 0.0;
  reig(IEIGN_P, 6) = Q(QRHO)*as*alf;

  reig(IEIGN_BT, 0) = AAs*betx;
  reig(IEIGN_BT, 1) = -betz*S*std::sqrt(Q(QRHO));
  reig(IEIGN_BT, 2) = -AAf*betx;
  reig(IEIGN_BT, 3) = 0.0;
  reig(IEIGN_BT, 4) = -AAf*betx;
  reig(IEIGN_BT, 5) = -betz*S*std::sqrt(Q(QRHO));
  reig(IEIGN_BT, 6) = AAs*betx;

  reig(IEIGN_BTT, 0) = AAs*betz;
  reig(IEIGN_BTT, 1) = betx*S*std::sqrt(Q(QRHO));
  reig(IEIGN_BTT, 2) = -AAf*betz;
  reig(IEIGN_BTT, 3) = 0.0;
  reig(IEIGN_BTT, 4) = -AAf*betz;
  reig(IEIGN_BTT, 5) = betx*S*std::sqrt(Q(QRHO));
  reig(IEIGN_BTT, 6) = AAs*betz;

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
evecz(Array2D<Real, 0, NEIGN-1, 0, NEIGN-1>& leig,
      Array2D<Real, 0, NEIGN-1, 0, NEIGN-1>& reig,
      const Real as_in,
      Array1D<Real, 0, NQ-1>& Q) {

  // returns Leig, where the rows are the left eigenvectors of the characteristic matrix Az

  // The characteristic speeds of the system

  Real as = as_in * as_in;

  // Alfven

  Real ca = (Q(QMAGX)*Q(QMAGX) + Q(QMAGY)*Q(QMAGY) + Q(QMAGZ)*Q(QMAGZ)) / Q(QRHO);
  Real caz = (Q(QMAGZ)*Q(QMAGZ)) / Q(QRHO);

  // Slow

  Real csz = 0.5_rt*((as + ca) - std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*caz));

  // Fast

  Real cfz = 0.5_rt*((as + ca) + std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*caz));

  // useful constants

  Real alf;
  Real als;

  if (std::abs(cfz-csz) <= 1e-14){
    alf = 1.0_rt;
    als = 0.0_rt;
  }

  if (as - csz < 0.0) {
    alf = 0.0_rt;
  } else {
    alf = std::sqrt((as - csz)/(cfz - csz));
  }

  if (cfz - as < 0.0) {
    als = 0.0_rt;
  } else {
    als = std::sqrt((cfz - as)/(cfz - csz));
  }

  Real betx;
  Real bety;

  if (std::abs(Q(QMAGX)) <= 1.e-14_rt && std::abs(Q(QMAGY)) <= 1.e-14_rt) {
    betx = 1.0_rt/std::sqrt(2.0_rt);
    bety = betx;

  } else {
    betx = Q(QMAGX) / std::sqrt(Q(QMAGX)*Q(QMAGX) + Q(QMAGY)*Q(QMAGY));
    bety = Q(QMAGY) / std::sqrt(Q(QMAGX)*Q(QMAGX) + Q(QMAGY)*Q(QMAGY));
  }

  Real cff = std::sqrt(cfz) * alf;
  Real css = std::sqrt(csz) * als;

  Real S = std::copysign(1.0_rt, Q(QMAGZ));

  Real Qf = std::sqrt(cfz) * alf * S;
  Real Qs = std::sqrt(csz) * als * S;

  Real N = 0.5_rt / as;

  Real AAf = std::sqrt(as) * alf * std::sqrt(Q(QRHO));
  Real AAs = std::sqrt(as) * als * std::sqrt(Q(QRHO));


  // w - cf
  leig(0, IEIGN_RHO) = 0.0;
  leig(0, IEIGN_U) = N*Qs*betx;
  leig(0, IEIGN_V) = N*Qs*bety;
  leig(0, IEIGN_W) = -N*cff;
  leig(0, IEIGN_P) = N*alf/Q(QRHO);
  leig(0, IEIGN_BT) = N*AAs*betx/Q(QRHO);
  leig(0, IEIGN_BTT) = N*AAs*bety/Q(QRHO);

  // w - caz
  leig(1, IEIGN_RHO) = 0.0;
  leig(1, IEIGN_U) = -0.5_rt*bety;
  leig(1, IEIGN_V) = 0.5_rt*betx;
  leig(1, IEIGN_W) = 0.0;
  leig(1, IEIGN_P) = 0.0;
  leig(1, IEIGN_BT) = -0.5_rt*bety*S/std::sqrt(Q(QRHO));
  leig(1, IEIGN_BTT) = 0.5_rt*betx*S/std::sqrt(Q(QRHO));

  // w - cs
  leig(2, IEIGN_RHO) = 0.0;
  leig(2, IEIGN_U) = -N*Qf*betx;
  leig(2, IEIGN_V) = -N*Qf*bety;
  leig(2, IEIGN_W) = -N*css;
  leig(2, IEIGN_P) = N*als/Q(QRHO);
  leig(2, IEIGN_BT) = -N*AAf*betx/Q(QRHO);
  leig(2, IEIGN_BTT) = -N*AAf*bety/Q(QRHO);

  // w
  leig(3, IEIGN_RHO) = 1.0_rt;
  leig(3, IEIGN_U) = 0.0;
  leig(3, IEIGN_V) = 0.0;
  leig(3, IEIGN_W) = 0.0;
  leig(3, IEIGN_P) = -1.0_rt/as;
  leig(3, IEIGN_BT) = 0.0;
  leig(3, IEIGN_BTT) = 0.0;

  // w + cs
  leig(4, IEIGN_RHO) = 0.0;
  leig(4, IEIGN_U) = N*Qf*betx;
  leig(4, IEIGN_V) = N*Qf*bety;
  leig(4, IEIGN_W) = N*css;
  leig(4, IEIGN_P) = N*als/Q(QRHO);
  leig(4, IEIGN_BT) = -N*AAf*betx/Q(QRHO);
  leig(4, IEIGN_BTT) = -N*AAf*bety/Q(QRHO);

  // w + caz
  leig(5, IEIGN_RHO) = 0.0;
  leig(5, IEIGN_U) = 0.5_rt*bety;
  leig(5, IEIGN_V) = -0.5_rt*betx;
  leig(5, IEIGN_W) = 0.0;
  leig(5, IEIGN_P) = 0.0;
  leig(5, IEIGN_BT) = -0.5_rt*bety*S/std::sqrt(Q(QRHO));
  leig(5, IEIGN_BTT) = 0.5_rt*betx*S/std::sqrt(Q(QRHO));

  // w + cf
  leig(6, IEIGN_RHO) = 0.0;
  leig(6, IEIGN_U) = -N*Qs*betx;
  leig(6, IEIGN_V) = -N*Qs*bety;
  leig(6, IEIGN_W) = N*cff;
  leig(6, IEIGN_P) = N*alf/Q(QRHO);
  leig(6, IEIGN_BT) = N*AAs*betx/Q(QRHO);
  leig(6, IEIGN_BTT) = N*AAs*bety/Q(QRHO);

  //  w - cf ; w - Caz ; w - cs ; w ; w + cs ; w + Caz ; w + cf

  reig(IEIGN_RHO, 0) = Q(QRHO)*alf;
  reig(IEIGN_RHO, 1) = 0.0;
  reig(IEIGN_RHO, 2) = Q(QRHO)*als;
  reig(IEIGN_RHO, 3) = 1.0_rt;
  reig(IEIGN_RHO, 4) = Q(QRHO)*als;
  reig(IEIGN_RHO, 5) = 0.0;
  reig(IEIGN_RHO, 6) = Q(QRHO)*alf;

  reig(IEIGN_W, 0) = -cff;
  reig(IEIGN_W, 1) = 0.0;
  reig(IEIGN_W, 2) = -css;
  reig(IEIGN_W, 3) = 0.0;
  reig(IEIGN_W, 4) = css;
  reig(IEIGN_W, 5) = 0.0;
  reig(IEIGN_W, 6) = cff;

  reig(IEIGN_U, 0) = Qs*betx;
  reig(IEIGN_U, 1) = -bety;
  reig(IEIGN_U, 2) = -Qf*betx;
  reig(IEIGN_U, 3) = 0.0;
  reig(IEIGN_U, 4) = Qf*betx;
  reig(IEIGN_U, 5) = bety;
  reig(IEIGN_U, 6) = -Qs*betx;

  reig(IEIGN_V, 0) = Qs*bety;
  reig(IEIGN_V, 1) = betx;
  reig(IEIGN_V, 2) = -Qf*bety;
  reig(IEIGN_V, 3) = 0.0;
  reig(IEIGN_V, 4) = Qf*bety;
  reig(IEIGN_V, 5) = -betx;
  reig(IEIGN_V, 6) = -Qs*bety;

  reig(IEIGN_P, 0) = Q(QRHO)*as*alf;
  reig(IEIGN_P, 1) = 0.0;
  reig(IEIGN_P, 2) = Q(QRHO)*as*als;
  reig(IEIGN_P, 3) = 0.0;
  reig(IEIGN_P, 4) = Q(QRHO)*as*als;
  reig(IEIGN_P, 5) = 0.0;
  reig(IEIGN_P, 6) = Q(QRHO)*as*alf;

  reig(IEIGN_BT, 0) = AAs*betx;
  reig(IEIGN_BT, 1) = -bety*S*std::sqrt(Q(QRHO));
  reig(IEIGN_BT, 2) = -AAf*betx;
  reig(IEIGN_BT, 3) = 0.0;
  reig(IEIGN_BT, 4) = -AAf*betx;
  reig(IEIGN_BT, 5) = -bety*S*std::sqrt(Q(QRHO));
  reig(IEIGN_BT, 6) = AAs*betx;

  reig(IEIGN_BTT, 0) = AAs*bety;
  reig(IEIGN_BTT, 1) = betx*S*std::sqrt(Q(QRHO));
  reig(IEIGN_BTT, 2) = -AAf*bety;
  reig(IEIGN_BTT, 3) = 0.0;
  reig(IEIGN_BTT, 4) = -AAf*bety;
  reig(IEIGN_BTT, 5) = betx*S*std::sqrt(Q(QRHO));
  reig(IEIGN_BTT, 6) = AAs*bety;

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
check_evecs(Array2D<Real, 0, NEIGN-1, 0, NEIGN-1>& leig,
            Array2D<Real, 0, NEIGN-1, 0, NEIGN-1>& reig) {

  // this checks if the eigenvectors are properly orthornomal

  for (int n = 0; n < NEIGN; n++) {
    for (int m = 0; m < NEIGN; m++) {

      // dot product of l_m . r_n
      Real dot = 0.0;
      for (int p = 0; p < NEIGN; p++)
      {
        dot += leig(m, p) * reig(p, n);
      }

#ifndef AMREX_USE_GPU
      std::cout << "l(" << m << ") . r(" << n << ") = " << dot << std::endl;
#endif

    }
#ifndef AMREX_USE_GPU
    std::cout << std::endl;
#endif
  }
}
#endif
