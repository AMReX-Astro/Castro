#ifndef CASTRO_RECONSTRUCTION_H
#define CASTRO_RECONSTRUCTION_H

constexpr int im2 = 0;
constexpr int im1 = 1;
constexpr int i0 = 2;
constexpr int ip1 = 3;
constexpr int ip2 = 4;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
load_stencil(Array4<Real const> const& q_arr, const int idir,
             const int i, const int j, const int k, const int ncomp,
             Real* s) {

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

#endif

