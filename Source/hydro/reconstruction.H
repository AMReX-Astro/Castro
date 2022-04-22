#ifndef CASTRO_RECONSTRUCTION_H
#define CASTRO_RECONSTRUCTION_H

constexpr int im3 = 0;
constexpr int im2 = 1;
constexpr int im1 = 2;
constexpr int i0 = 3;
constexpr int ip1 = 4;
constexpr int ip2 = 5;
constexpr int ip3 = 6;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
load_stencil(Array4<Real const> const& q_arr, const int idir,
             const int i, const int j, const int k, const int ncomp,
             Real* s) {

    if (idir == 0) {
        s[im3] = q_arr(i-3,j,k,ncomp);
        s[im2] = q_arr(i-2,j,k,ncomp);
        s[im1] = q_arr(i-1,j,k,ncomp);
        s[i0]  = q_arr(i,j,k,ncomp);
        s[ip1] = q_arr(i+1,j,k,ncomp);
        s[ip2] = q_arr(i+2,j,k,ncomp);
        s[ip3] = q_arr(i+3,j,k,ncomp);

      } else if (idir == 1) {
        s[im3] = q_arr(i,j-3,k,ncomp);
        s[im2] = q_arr(i,j-2,k,ncomp);
        s[im1] = q_arr(i,j-1,k,ncomp);
        s[i0]  = q_arr(i,j,k,ncomp);
        s[ip1] = q_arr(i,j+1,k,ncomp);
        s[ip2] = q_arr(i,j+2,k,ncomp);
        s[ip3] = q_arr(i,j+3,k,ncomp);

      } else {
        s[im3] = q_arr(i,j,k-3,ncomp);
        s[im2] = q_arr(i,j,k-2,ncomp);
        s[im1] = q_arr(i,j,k-1,ncomp);
        s[i0]  = q_arr(i,j,k,ncomp);
        s[ip1] = q_arr(i,j,k+1,ncomp);
        s[ip2] = q_arr(i,j,k+2,ncomp);
        s[ip3] = q_arr(i,j,k+3,ncomp);

      }

}

#endif

