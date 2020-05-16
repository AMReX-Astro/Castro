#include "Castro.H"
#include "Castro_F.H"
#include "Castro_util.H"
#include "Castro_hydro_F.H"

#include "mhd_sound_speed.H"

using namespace amrex;

void
Castro::mhd_speeds(const Box& bx,
                   Array4<Real const> const& Bx,
                   Array4<Real const> const& By,
                   Array4<Real const> const& Bz,
                   Array4<Real const> const& q_arr,
                   Array4<Real const> const& qaux_arr,
                   Array4<Real> const& cx,
                   Array4<Real> const& cy,
                   Array4<Real> const& cz) {

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    //sound speed for ideal mhd
    Real as2 = qaux_arr(i,j,k,QC) * qaux_arr(i,j,k,QC);
    Real rho_inv = 1.0_rt / q_arr(i,j,k,QRHO);

    Real ca = (q_arr(i,j,k,QMAGX) * q_arr(i,j,k,QMAGX) +
               q_arr(i,j,k,QMAGY) * q_arr(i,j,k,QMAGY) +
               q_arr(i,j,k,QMAGZ) * q_arr(i,j,k,QMAGZ)) * rho_inv;

    Real ctmp = 0.0;

    Real cad = q_arr(i,j,k,QMAGX) * q_arr(i,j,k,QMAGX) * rho_inv;
    eos_soundspeed_mhd(ctmp, as2, ca, cad);
    cx(i,j,k) = ctmp;

    cad = q_arr(i,j,k,QMAGY) * q_arr(i,j,k,QMAGY) * rho_inv;
    eos_soundspeed_mhd(ctmp, as2, ca, cad);
    cy(i,j,k) = ctmp;

    cad = q_arr(i,j,k,QMAGZ) * q_arr(i,j,k,QMAGZ) * rho_inv;
    eos_soundspeed_mhd(ctmp, as2, ca, cad);
    cz(i,j,k) = ctmp;

  });

}


void
Castro::consup_mhd(const Box& bx,
                   Array4<Real> const& update,
                   Array4<Real const> const& flux0,
                   Array4<Real const> const& flux1,
                   Array4<Real const> const& flux2) {

  // do the conservative update and store - div{F} in update.  Note,
  // in contrast to the CTU hydro case, we don't add the pdivu term to
  // (rho e) here, because we included that in the hydro source terms.

  const auto dx = geom.CellSizeArray();

  Real dxinv = 1.0_rt/dx[0];
  Real dyinv = 1.0_rt/dx[1];
  Real dzinv = 1.0_rt/dx[2];

  amrex::ParallelFor(bx, NUM_STATE,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n) noexcept
  {

    if (n == UTEMP) {
      update(i,j,k,n) = 0.0_rt;
#ifdef SHOCK_VAR
    } else if (n == USHK) {
      update(i,j,k,n) = 0.0_rt;
#endif
    } else {
      update(i,j,k,n) =
        (flux0(i,j,k,n) - flux0(i+1,j,k,n)) * dxinv +
        (flux1(i,j,k,n) - flux1(i,j+1,k,n)) * dyinv +
        (flux2(i,j,k,n) - flux2(i,j,k+1,n)) * dzinv;
    }

  });

}
