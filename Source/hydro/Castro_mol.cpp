#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#include <ppm.H>

using namespace amrex;


void
Castro::mol_ppm_reconstruct(const Box& bx,
                            const int idir,
                            Array4<Real const> const q_arr,
                            Array4<Real const> const flatn_arr,
                            Array4<Real> const qm,
                            Array4<Real> const qp) {

  AMREX_PARALLEL_FOR_4D(bx, NQ, i, j, k, n,
  {

    Real s[5];
    Real flat = flatn_arr(i,j,k);
    Real sm;
    Real sp;

    if (idir == 0) {
      s[im2] = q_arr(i-2,j,k,n);
      s[im1] = q_arr(i-1,j,k,n);
      s[i0]  = q_arr(i,j,k,n);
      s[ip1] = q_arr(i+1,j,k,n);
      s[ip2] = q_arr(i+2,j,k,n);

    } else if (idir == 1) {
      s[im2] = q_arr(i,j-2,k,n);
      s[im1] = q_arr(i,j-1,k,n);
      s[i0]  = q_arr(i,j,k,n);
      s[ip1] = q_arr(i,j+1,k,n);
      s[ip2] = q_arr(i,j+2,k,n);

    } else {
      s[im2] = q_arr(i,j,k-2,n);
      s[im1] = q_arr(i,j,k-1,n);
      s[i0]  = q_arr(i,j,k,n);
      s[ip1] = q_arr(i,j,k+1,n);
      s[ip2] = q_arr(i,j,k+2,n);

    }

    ppm_reconstruct(s, flat, sm, sp);

    if (idir == 1) {
      // right state at i-1/2
      qp(i,j,k,n) = sm;

      // left state at i+1/2
      qm(i+1,j,k,n) = sp;

    } else if (idir == 2) {
      // right state at j-1/2
      qp(i,j,k,n) = sm;

      // left state at j+1/2
      qm(i,j+1,k,n) = sp;

    } else {
      // right state at k-1/2
      qp(i,j,k,n) = sm;

      // left state at k+1/2
      qm(i,j,k+1,n) = sp;

    }

  });

}

