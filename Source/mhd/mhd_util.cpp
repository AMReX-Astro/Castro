#include <Castro.H>
#include <Castro_util.H>

#include <mhd_util.H>

using namespace amrex;

void
Castro::consup_mhd(const Box& bx, const Real dt,
                   Array4<Real> const& U_new,
                   Array4<Real const> const& flux0,
                   Array4<Real const> const& flux1,
                   Array4<Real const> const& flux2) {

  // do the conservative update.  Note, in contrast to the CTU hydro
  // case, we don't add the pdivu term to (rho e) here, because we
  // included that in the hydro source terms.

  const auto dx = geom.CellSizeArray();

  Real dxinv = 1.0_rt/dx[0];
#if AMREX_SPACEDIM >= 2
  Real dyinv = 1.0_rt/dx[1];
#endif
#if AMREX_SPACEDIM == 3
  Real dzinv = 1.0_rt/dx[2];
#endif

  amrex::ParallelFor(bx, NUM_STATE,
  [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {

    if (n == UTEMP) {
      U_new(i,j,k,n) = 0.0_rt;
#ifdef SHOCK_VAR
    } else if (n == USHK) {
      U_new(i,j,k,n) = 0.0_rt;
#endif
    } else {
      U_new(i,j,k,n) += dt * (flux0(i,j,k,n) - flux0(i+1,j,k,n)) * dxinv;
#if AMREX_SPACEDIM >= 2
      U_new(i,j,k,n) += dt * (flux1(i,j,k,n) - flux1(i,j+1,k,n)) * dyinv;
#endif
#if AMREX_SPACEDIM == 3
      U_new(i,j,k,n) += dt * (flux2(i,j,k,n) - flux2(i,j,k+1,n)) * dzinv;
#endif
    }

  });

}


void
Castro::PrimToCons(const Box& bx,
                   Array4<Real const> const& q_arr,
                   Array4<Real> const& u_arr) {

  // calculate the conserved variables from the primitive

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    u_arr(i,j,k,URHO) = q_arr(i,j,k,QRHO);
    u_arr(i,j,k,UMX) = q_arr(i,j,k,QRHO)*q_arr(i,j,k,QU);
    u_arr(i,j,k,UMY) = q_arr(i,j,k,QRHO)*q_arr(i,j,k,QV);
    u_arr(i,j,k,UMZ) = q_arr(i,j,k,QRHO)*q_arr(i,j,k,QW);

    eos_t eos_state;
    eos_state.rho = q_arr(i,j,k,QRHO);
    eos_state.p = q_arr(i,j,k,QPRES);
    eos_state.T = 100.0_rt; // dummy initial T.
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = q_arr(i,j,k,QFS+n);
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = q_arr(i,j,k,QFX+n);
    }
#endif

    eos(eos_input_rp, eos_state);

    u_arr(i,j,k,UEDEN) = eos_state.rho * eos_state.e +
      + 0.5_rt * q_arr(i,j,k,QRHO) * (q_arr(i,j,k,QU) * q_arr(i,j,k,QU) +
                                      q_arr(i,j,k,QV) * q_arr(i,j,k,QV) +
                                      q_arr(i,j,k,QW) * q_arr(i,j,k,QW)) +
      + 0.5_rt * (q_arr(i,j,k,QMAGX) * q_arr(i,j,k,QMAGX) +
                  q_arr(i,j,k,QMAGY) * q_arr(i,j,k,QMAGY) +
                  q_arr(i,j,k,QMAGZ) * q_arr(i,j,k,QMAGZ));

    u_arr(i,j,k,UEINT) = eos_state.rho * eos_state.e;
    u_arr(i,j,k,UTEMP) = eos_state.T;

    u_arr(i,j,k,UMAGX) = q_arr(i,j,k,QMAGX);
    u_arr(i,j,k,UMAGY) = q_arr(i,j,k,QMAGY);
    u_arr(i,j,k,UMAGZ) = q_arr(i,j,k,QMAGZ);

    // species
    for (int n = 0; n < NumSpec; n++) {
      u_arr(i,j,k,UFS+n) = q_arr(i,j,k,QRHO) * q_arr(i,j,k,QFS+n);
    }

  });
}


void
Castro::prim_half(const Box& bx,
                  Array4<Real> const& q2D,
                  Array4<Real const> const& q_arr,
                  Array4<Real const> const& flxx,
                  Array4<Real const> const& flxy,
                  Array4<Real const> const& flxz,
                  const Real dt) {

  // Find the 2D corrected primitive variables and n+1/2

  auto dx = geom.CellSizeArray();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    Real divF[NUM_STATE+3];
    Real divF_q[NQ];
    Real q_zone[NQ];

    for (int n = 0; n < NUM_STATE+3; n++) {
      divF[n] = (flxx(i+1,j,k,n) - flxx(i,j,k,n)) / dx[0];
#if AMREX_SPACEDIM >= 2
      divF[n] += (flxy(i,j+1,k,n) - flxy(i,j,k,n)) / dx[1];
#endif
#if AMREX_SPACEDIM == 3
      divF[n] += (flxz(i,j,k+1,n) - flxz(i,j,k,n)) / dx[2];
#endif
    }

    // that is a flux of conserved variables -- transform it to primitive
    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j,k,n);
    }

    qflux(divF_q, divF, q_zone);

    // Right below eq. 48
    for (int n = 0; n < NQ; n++) {
      q2D(i,j,k,n) = q_arr(i,j,k,n) - 0.5_rt * dt * divF_q[n];
    }
  });
}
