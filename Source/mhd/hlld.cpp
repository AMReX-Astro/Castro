#include <Castro.H>

#include <mhd_util.H>

using namespace amrex;

void
Castro::hlld(const Box& bx,
             Array4<Real const> const& qleft,
             Array4<Real const> const& qright,
             Array4<Real> const& flx,
             const int dir) {

  // Riemann solve:

  // Main assumption, the normal velocity/Mag field is constant in the
  // Riemann fan, and is sM/Bn respectively.  Total Pressure is constant
  // throughout the Riemann fan, pst!


  // `n` here is the normal
  // `p` are the perpendicular

  int QMAGN, QMAGP1, QMAGP2;
  int QVELN, QVELP1, QVELP2;
  int UMN, UMP1, UMP2;
  int UMAGN, UMAGP1, UMAGP2;

  if (dir == 0) {
    QMAGN  = QMAGX;
    QMAGP1 = QMAGY;
    QMAGP2 = QMAGZ;
    QVELN  = QU;
    QVELP1 = QV;
    QVELP2 = QW;
    UMN    = UMX;
    UMP1   = UMY;
    UMP2   = UMZ;
    UMAGN  = UMAGX;
    UMAGP1 = UMAGY;
    UMAGP2 = UMAGZ;

  } else if (dir == 1) {
    QMAGN  = QMAGY;
    QMAGP1 = QMAGZ;
    QMAGP2 = QMAGX;
    QVELN  = QV;
    QVELP1 = QW;
    QVELP2 = QU;
    UMN    = UMY;
    UMP1   = UMZ;
    UMP2   = UMX;
    UMAGN  = UMAGY;
    UMAGP1 = UMAGZ;
    UMAGP2 = UMAGX;

  } else {  // dir == 2
    QMAGN  = QMAGZ;
    QMAGP1 = QMAGX;
    QMAGP2 = QMAGY;
    QVELN  = QW;
    QVELP1 = QU;
    QVELP2 = QV;
    UMN    = UMZ;
    UMP1   = UMX;
    UMP2   = UMY;
    UMAGN  = UMAGZ;
    UMAGP1 = UMAGX;
    UMAGP2 = UMAGY;
  }

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    // this is a loop over interfaces, so, e.g., for idir = 0 (x), we are seeing
    // q_{i-1/2,j,k,L} and q_{i-1/2,j,k,R}

    Array1D<Real, 0, NQ-1> qL;
    Array1D<Real, 0, NQ-1> qR;

    for (int n = 0; n < NQ; n++) {
      qL(n) = qleft(i,j,k,n);
      qR(n) = qright(i,j,k,n);
    }

    // check to not go lower than small_dens and small_press
    qL(QRHO) = amrex::max(small_dens, qL(QRHO));
    qR(QRHO) = amrex::max(small_dens, qR(QRHO));
    qL(QPRES) = amrex::max(small_pres, qL(QPRES));
    qR(QPRES) = amrex::max(small_pres, qR(QPRES));

    Array1D<Real, 0, NUM_STATE+2> uL;
    Real gam1_L;

    PToC(qL, uL, gam1_L);

    // Compute the fluxes.  Here use total p not just p_g eq.11 in Miniati
    // Also note that Miniati has an overall sign error in the B fluxes
    Real BL2 = qL(QMAGX) * qL(QMAGX) + qL(QMAGY) * qL(QMAGY) + qL(QMAGZ) * qL(QMAGZ);
    Real UBL = qL(QMAGX) * qL(QU) + qL(QMAGY) * qL(QV) + qL(QMAGZ) * qL(QW);
    Array1D<Real, 0, NUM_STATE+2> FL;

    FL(URHO) = qL(QRHO) * qL(QVELN);
    FL(UMN) = qL(QRHO) * qL(QVELN) * qL(QVELN) + (qL(QPRES) + 0.5_rt * BL2) - qL(QMAGN) * qL(QMAGN);
    FL(UMP1) = qL(QRHO) * qL(QVELN) * qL(QVELP1) - qL(QMAGN) * qL(QMAGP1);
    FL(UMP2) = qL(QRHO) * qL(QVELN) * qL(QVELP2) - qL(QMAGN) * qL(QMAGP2);
    FL(UEDEN) = qL(QVELN) * (uL(UEDEN) + (qL(QPRES) + 0.5_rt * BL2)) - qL(QMAGN) * UBL;
    FL(UMAGN) = 0.0;
    FL(UMAGP1) = qL(QVELN) * qL(QMAGP1) - qL(QVELP1) * qL(QMAGN);
    FL(UMAGP2) = qL(QVELN) * qL(QMAGP2) - qL(QVELP2) * qL(QMAGN);
    for (int n = 0; n < NumSpec; n++) {
      FL(UFS+n) = qL(QVELN) * uL(UFS+n);
    }
    FL(UEINT) = qL(QVELN) * uL(UEINT);
    FL(UTEMP) = 0.0_rt;

    Array1D<Real, 0, NUM_STATE+2> uR;
    Real gam1_R;

    PToC(qR, uR, gam1_R);

    Real BR2 = qR(QMAGX) * qR(QMAGX) + qR(QMAGY) * qR(QMAGY) + qR(QMAGZ) * qR(QMAGZ);
    Real UBR = qR(QMAGX) * qR(QU) + qR(QMAGY) * qR(QV) + qR(QMAGZ) * qR(QW);

    Array1D<Real, 0, NUM_STATE+2> FR;

    FR(URHO) = qR(QRHO) * qR(QVELN);
    FR(UMN) = qR(QRHO) * qR(QVELN) * qR(QVELN) + (qR(QPRES) + 0.5_rt * BR2) - qR(QMAGN) * qR(QMAGN);
    FR(UMP1) = qR(QRHO) * qR(QVELN) * qR(QVELP1) - qR(QMAGN) * qR(QMAGP1);
    FR(UMP2) = qR(QRHO) * qR(QVELN) * qR(QVELP2) - qR(QMAGN) * qR(QMAGP2);
    FR(UEDEN) = qR(QVELN) * (uR(UEDEN) + (qR(QPRES) + 0.5_rt * BR2)) - qR(QMAGN) * UBR;
    FR(UMAGN) = 0.0;
    FR(UMAGP1) = qR(QVELN) * qR(QMAGP1) - qR(QVELP1) * qR(QMAGN);
    FR(UMAGP2) = qR(QVELN) * qR(QMAGP2) - qR(QVELP2) * qR(QMAGN);
    for (int n = 0; n < NumSpec; n++) {
      FR(UFS+n) = qR(QVELN) * uR(UFS+n);
    }
    FR(UEINT) = qR(QVELN) * uR(UEINT);
    FR(UTEMP) = 0.0_rt;

    // From Miyoshi and Kusano paper eq.(3)

    Real asL = gam1_L * qL(QPRES) / qL(QRHO);
    Real asR = gam1_R * qR(QPRES) / qR(QRHO);

    Real caL  = BL2 / qL(QRHO); // Magnetic Speeds
    Real caR  = BR2 / qR(QRHO);

    Real canL = qL(QMAGN) * qL(QMAGN) / qL(QRHO);
    Real canR = qR(QMAGN) * qR(QMAGN) / qR(QRHO);

    // fast waves
    // cf as in equation (3) of HLLD paper (Miyoshi and Kusano)

    Real cfL2 = std::sqrt(0.5_rt * ((asL + caL) + std::sqrt((asL + caL)*(asL + caL) - 4.0_rt * asL * canL)));
    Real cfR = std::sqrt(0.5_rt * ((asR + caR) + std::sqrt((asR + caR)*(asR + caR) - 4.0_rt * asR * canR)));

    // Riemann Speeds
    // sL and sR, eq.(12) (Miyoshi and Kusano)

    Real sL = amrex::min(qL(QVELN) - cfL2, qR(QVELN) - cfR);
    Real sR = amrex::max(qL(QVELN) + cfL2, qR(QVELN) + cfR);

    // Pressures in the Riemann Fan

    Real ptL = qL(QPRES) + 0.5_rt * BL2;
    Real ptR = qR(QPRES) + 0.5_rt * BR2;

    // sM -- the entropy wave, eq.(38)

    Real sM  = ((sR - qR(QVELN)) * qR(QRHO) * qR(QVELN) - (sL - qL(QVELN)) * qL(QRHO) * qL(QVELN) -
           ptR + ptL) /
      ((sR - qR(QVELN)) * qR(QRHO) - (sL - qL(QVELN)) * qL(QRHO));


    // pst, eq.(41) (pstar_total)

    Real pst  = (sR - qR(QVELN)) * qR(QRHO) * ptL - (sL - qL(QVELN)) * qL(QRHO) * ptR +
      qL(QRHO) * qR(QRHO) * (sR - qR(QVELN)) * (sL - qL(QVELN)) * (qR(QVELN) - qL(QVELN));
    pst  = pst / ((sR - qR(QVELN)) * qR(QRHO) - (sL - qL(QVELN)) * qL(QRHO));

    // Density * states

    Array1D<Real, 0, NUM_STATE+2> UsL;
    Array1D<Real, 0, NUM_STATE+2> UsR;

    // Density eq.(43)

    UsL(URHO) = amrex::max(small_dens, qL(QRHO) * ((sL - qL(QVELN)) / (sL - sM)));
    UsR(URHO) = amrex::max(small_dens, qR(QRHO) * ((sR - qR(QVELN)) / (sR - sM)));

    // Species * states

    for (int n = 0; n < NumSpec; n++) {
      UsL(UFS+n) = qL(QFS+n) * UsL(URHO);
      UsR(UFS+n) = qR(QFS+n) * UsR(URHO);
    }

    // e

    UsL(UEINT) = uL(UEINT) / uL(URHO) * UsL(URHO);
    UsR(UEINT) = uR(UEINT) / uR(URHO) * UsR(URHO);

    // Vel * states

    // Normal dir (Eq. 39)

    UsL(UMN) = sM;
    UsR(UMN) = sM;

    // Perpendicular dir (Eq 44)
    // Second Perpendicular dir (Eq 46)

    Real denom_upL = qL(QRHO) * (sL - qL(QVELN)) * (sL - sM) - qL(QMAGN) * qL(QMAGN);
    Real denom_upR = qR(QRHO) * (sR - qR(QVELN)) * (sR - sM) - qR(QMAGN)*qR(QMAGN);

    if (std::abs(denom_upL) < 1.e-14_rt) {
      UsL(UMP1) = qL(QVELP1);
      UsL(UMP2) = qL(QVELP2);

    } else {
      UsL(UMP1) = qL(QVELP1) - qL(QMAGN)*qL(QMAGP1) *
        ((sM - qL(QVELN)) / denom_upL);
      UsL(UMP2) = qL(QVELP2) - qL(QMAGN) * qL(QMAGP2) *
        ((sM - qL(QVELN)) / denom_upL);

    }

    if (std::abs(denom_upR) < 1.e-14_rt) {
      UsR(UMP1) = qR(QVELP1);
      UsR(UMP2) = qR(QVELP2);

    } else {
      UsR(UMP1) = qR(QVELP1) - qR(QMAGN)*qR(QMAGP1) *
        ((sM - qR(QVELN)) / denom_upR);
      UsR(UMP2) = qR(QVELP2) - qR(QMAGN) * qR(QMAGP2) *
        ((sM - qR(QVELN)) / denom_upR);

    }


    UsL(UMX) = UsL(UMX) * UsL(URHO);
    UsL(UMY) = UsL(UMY) * UsL(URHO);
    UsL(UMZ) = UsL(UMZ) * UsL(URHO);

    UsR(UMX) = UsR(UMX) * UsR(URHO);
    UsR(UMY) = UsR(UMY) * UsR(URHO);
    UsR(UMZ) = UsR(UMZ) * UsR(URHO);

    // B * states

    // Normal dir

    UsL(UMAGN) = qL(QMAGN);
    UsR(UMAGN) = qR(QMAGN);

    // Perpendicular dir (Eq. 45)
    // Second Perpendicular dir (Eq. 47)

    Real denom_bpL = qL(QRHO) * (sL - qL(QVELN)) * (sL - sM) - qL(QMAGN) * qL(QMAGN);
    Real denom_bpR = qR(QRHO) * (sR - qR(QVELN)) * (sR - sM) - qR(QMAGN) * qR(QMAGN);

    if (std::abs(denom_bpL) < 1.e-14_rt) {
      UsL(UMAGP1) = 0.0_rt;
      UsL(UMAGP2) = 0.0_rt;

    } else {
      UsL(UMAGP1) = qL(QMAGP1) * (qL(QRHO) * (sL - qL(QVELN)) * (sL - qL(QVELN)) - qL(QMAGN) * qL(QMAGN)) /
                    denom_bpL;
      UsL(UMAGP2) = qL(QMAGP2) * (qL(QRHO) * (sL - qL(QVELN)) * (sL - qL(QVELN)) - qL(QMAGN) * qL(QMAGN)) /
                    denom_bpL;

    }

    if (std::abs(denom_bpR) < 1.e-14_rt) {
      UsR(UMAGP1) = 0.0_rt;
      UsR(UMAGP2) = 0.0_rt;

    } else {
      UsR(UMAGP1) = qR(QMAGP1) * (qR(QRHO) * (sR - qR(QVELN)) * (sR - qR(QVELN)) - qR(QMAGN) * qR(QMAGN)) /
                    denom_bpR;
      UsR(UMAGP2) = qR(QMAGP2) * (qR(QRHO) * (sR - qR(QVELN)) * (sR - qR(QVELN)) - qR(QMAGN) * qR(QMAGN)) /
                    denom_bpR;

    }


    // Energy, eq.(48)

    UsL(UEDEN) = (sL - qL(QVELN)) * uL(UEDEN) - ptL * qL(QVELN) + pst * sM +
      qL(QMAGN) * (UBL - (UsL(UMX) * UsL(UMAGX) + UsL(UMY) * UsL(UMAGY) + UsL(UMZ) * UsL(UMAGZ)) / UsL(URHO));
    UsL(UEDEN) = UsL(UEDEN) / (sL - sM);

    UsR(UEDEN) = (sR - qR(QVELN)) * uR(UEDEN) - ptR * qR(QVELN) + pst * sM +
      qR(QMAGN) * (UBR - (UsR(UMX) * UsR(UMAGX) + UsR(UMY) * UsR(UMAGY) + UsR(UMZ) * UsR(UMAGZ)) / UsR(URHO));
    UsR(UEDEN) = UsR(UEDEN) / (sR - sM);

    UsL(UTEMP) = 0.0_rt;
    UsR(UTEMP) = 0.0_rt;

    // speeds, eq.(51)

    Real ssL = sM - std::abs(qL(QMAGN)) / std::sqrt(UsL(URHO));
    Real ssR = sM + std::abs(qR(QMAGN)) / std::sqrt(UsR(URHO));

    // ** states

    Array1D<Real, 0, NUM_STATE+2> UssL;
    Array1D<Real, 0, NUM_STATE+2> UssR;

    // Dens (Eq. 49)

    UssL(URHO) = UsL(URHO);
    UssR(URHO) = UsR(URHO);

    // species

    for (int n = 0; n < NumSpec; n++) {
      UssL(UFS+n) = UsL(UFS+n);
      UssR(UFS+n) = UsR(UFS+n);
    }

    // e
    UssL(UEINT) = UsL(UEINT);
    UssR(UEINT) = UsR(UEINT);

    // Vel in normal direction (Eq. 39)

    UssL(UMN) = sM;
    UssR(UMN) = sM;

    // Vel in perpendicular direction, eq.(59)

    UssL(UMP1) = (std::sqrt(UsL(URHO)) * UsL(UMP1) / UsL(URHO) +
                  std::sqrt(UsR(URHO)) * UsR(UMP1) / UsR(URHO) +
                  (UsR(UMAGP1) - UsL(UMAGP1)) * std::copysign(1.0_rt, qL(QMAGN))) /
      (std::sqrt(UsL(URHO)) + std::sqrt(UsR(URHO)));
    UssR(UMP1) = UssL(UMP1);

    // Vel in second perpendicular direction, eq(60)

    UssL(UMP2) = (std::sqrt(UsL(URHO)) * UsL(UMP2) / UsL(URHO) +
                  std::sqrt(UsR(URHO)) * UsR(UMP2) / UsR(URHO) +
                  (UsR(UMAGP2) - UsL(UMAGP2)) * std::copysign(1.0_rt, qL(QMAGN))) /
      (std::sqrt(UsL(URHO)) + std::sqrt(UsR(URHO)));
    UssR(UMP2) = UssL(UMP2);

    UssL(UMX) = UssL(UMX) * UssL(URHO);
    UssL(UMY) = UssL(UMY) * UssL(URHO);
    UssL(UMZ) = UssL(UMZ) * UssL(URHO);

    UssR(UMX) = UssR(UMX) * UssR(URHO);
    UssR(UMY) = UssR(UMY) * UssR(URHO);
    UssR(UMZ) = UssR(UMZ) * UssR(URHO);

    // B in normal direction

    UssL(UMAGN) = UsL(UMAGN);
    UssR(UMAGN) = UsR(UMAGN);

    // B in perpendicular direction (Eq 61)

    UssL(UMAGP1) = (std::sqrt(UsL(URHO)) * UsR(UMAGP1) + std::sqrt(UsR(URHO)) * UsL(UMAGP1) +
                    std::sqrt(UsL(URHO) * UsR(URHO)) * (UsR(UMP1) / UsR(URHO) -
                                                        UsL(UMP1) / UsL(URHO)) * std::copysign(1.0_rt, qL(QMAGN))) /
      (std::sqrt(UsL(URHO)) + std::sqrt(UsR(URHO)));

    UssR(UMAGP1) = UssL(UMAGP1);

    // B in second perpendicular direction (Eq 62)

    UssL(UMAGP2) = (std::sqrt(UsL(URHO)) * UsR(UMAGP2) + std::sqrt(UsR(URHO)) * UsL(UMAGP2) +
                    std::sqrt(UsL(URHO) * UsR(URHO)) * (UsR(UMP2) / UsR(URHO) -
                                                        UsL(UMP2) / UsL(URHO)) * std::copysign(1.0_rt, qL(QMAGN))) /
      (std::sqrt(UsL(URHO)) + std::sqrt(UsR(URHO)));
    UssR(UMAGP2) = UssL(UMAGP2);

    // Energy , eq.(63)

    UssL(UEDEN) = UsL(UEDEN) - std::sqrt(UsL(URHO)) *
      ((UsL(UMX) * UsL(UMAGX) + UsL(UMY) * UsL(UMAGY) + UsL(UMZ) * UsL(UMAGZ)) / UsL(URHO) -
       (UssL(UMX) * UssL(UMAGX) + UssL(UMY) * UssL(UMAGY) + UssL(UMZ) * UssL(UMAGZ)) / UssL(URHO)) *
      std::copysign(1.0_rt, qL(QMAGN));
    UssR(UEDEN) = UsR(UEDEN) + std::sqrt(UsR(QRHO)) *
      ((UsR(UMX) * UsR(UMAGX) + UsR(UMY) * UsR(UMAGY) + UsR(UMZ) * UsR(UMAGZ)) / UsR(URHO) -
       (UssR(UMX) * UssR(UMAGX) + UssR(UMY) * UssR(UMAGY) + UssR(UMZ) * UssR(UMAGZ)) / UssR(URHO)) *
      std::copysign(1.0_rt, qR(QMAGN));

    UssL(UTEMP) = 0.0_rt;
    UssR(UTEMP) = 0.0_rt;

    // fluxes (eq 64 and 65)

    // Solve the RP

    if (sL > 0.0) {
      for (int n = 0; n < NUM_STATE+3; n++) {
        flx(i,j,k,n) = FL(n);
      }

    } else if (sL <= 0.0 && ssL > 0.0) {
      // FsL  = FL + sL*UsL - sL*uL
      for (int n = 0; n < NUM_STATE+3; n++) {
        flx(i,j,k,n) = FL(n) + sL * UsL(n) - sL * uL(n);
      }

    } else if (ssL <= 0.0 && sM > 0.0) {
      // FssL = FL + ssL*UssL - (ssL-sL)*UsL - sL*uL
      for (int n = 0; n < NUM_STATE+3; n++) {
        flx(i,j,k,n) = FL(n) + ssL * UssL(n) - (ssL - sL) * UsL(n) - sL * uL(n);
      }

    } else if (sM <= 0.0 && ssR > 0.0) {
      // FssR = FR + ssR*UssR - (ssR-sR)*UsR - sR*uR
      for (int n = 0; n < NUM_STATE+3; n++) {
        flx(i,j,k,n) = FR(n) + ssR * UssR(n) - (ssR - sR) * UsR(n) - sR * uR(n);
      }

    } else if (ssR <= 0.0 && sR > 0.0) {
      // FsR  = FR + sR*UsR - sR*uR
      for (int n = 0; n < NUM_STATE+3; n++) {
        flx(i,j,k,n) = FR(n) + sR * UsR(n) - sR * uR(n);
      }

    } else {
      for (int n = 0; n < NUM_STATE+3; n++) {
        flx(i,j,k,n) = FR(n);
      }
    }

    flx(i,j,k,UTEMP) = 0.0;
  });
}

