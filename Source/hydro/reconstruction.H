#ifndef CASTRO_RECONSTRUCTION_H
#define CASTRO_RECONSTRUCTION_H

namespace reconstruction {
   enum slope_indices {
       im2 = 0,
       im1,
       i0,
       ip1,
       ip2,
       nslp
   };
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
load_stencil(Array4<Real const> const& q_arr, const int idir,
             const int i, const int j, const int k, const int ncomp,
             Real* s) {

    using namespace reconstruction;

    if (idir == 0) {
        s[im2] = q_arr(i-2,j,k,ncomp);
        s[im1] = q_arr(i-1,j,k,ncomp);
        s[i0]  = q_arr(i,j,k,ncomp);
        s[ip1] = q_arr(i+1,j,k,ncomp);
        s[ip2] = q_arr(i+2,j,k,ncomp);

      } else if (idir == 1) {
        s[im2] = q_arr(i,j-2,k,ncomp);
        s[im1] = q_arr(i,j-1,k,ncomp);
        s[i0]  = q_arr(i,j,k,ncomp);
        s[ip1] = q_arr(i,j+1,k,ncomp);
        s[ip2] = q_arr(i,j+2,k,ncomp);

      } else {
        s[im2] = q_arr(i,j,k-2,ncomp);
        s[im1] = q_arr(i,j,k-1,ncomp);
        s[i0]  = q_arr(i,j,k,ncomp);
        s[ip1] = q_arr(i,j,k+1,ncomp);
        s[ip2] = q_arr(i,j,k+2,ncomp);

      }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
load_passive_stencil(Array4<Real const> const& U_arr, Array4<Real const> const& rho_inv_arr, const int idir,
                     const int i, const int j, const int k, const int ncomp,
                     Real* s) {

    using namespace reconstruction;

    if (idir == 0) {
        s[im2] = U_arr(i-2,j,k,ncomp) * rho_inv_arr(i-2,j,k);
        s[im1] = U_arr(i-1,j,k,ncomp) * rho_inv_arr(i-1,j,k);
        s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);
        s[ip1] = U_arr(i+1,j,k,ncomp) * rho_inv_arr(i+1,j,k);
        s[ip2] = U_arr(i+2,j,k,ncomp) * rho_inv_arr(i+2,j,k);

      } else if (idir == 1) {
        s[im2] = U_arr(i,j-2,k,ncomp) * rho_inv_arr(i,j-2,k);
        s[im1] = U_arr(i,j-1,k,ncomp) * rho_inv_arr(i,j-1,k);
        s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);
        s[ip1] = U_arr(i,j+1,k,ncomp) * rho_inv_arr(i,j+1,k);
        s[ip2] = U_arr(i,j+2,k,ncomp) * rho_inv_arr(i,j+2,k);

      } else {
        s[im2] = U_arr(i,j,k-2,ncomp) * rho_inv_arr(i,j,k-2);
        s[im1] = U_arr(i,j,k-1,ncomp) * rho_inv_arr(i,j,k-1);
        s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);
        s[ip1] = U_arr(i,j,k+1,ncomp) * rho_inv_arr(i,j,k+1);
        s[ip2] = U_arr(i,j,k+2,ncomp) * rho_inv_arr(i,j,k+2);

      }

}

#endif

