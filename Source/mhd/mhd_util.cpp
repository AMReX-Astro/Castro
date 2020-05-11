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
