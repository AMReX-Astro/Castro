#ifndef CASTRO_RECONSTRUCTION_H
#define CASTRO_RECONSTRUCTION_H

#include <state_indices.H>

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
load_stencil(amrex::Array4<amrex::Real const> const& q_arr, const int idir,
             const int i, const int j, const int k, const int ncomp,
             amrex::Real* s) {

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
load_passive_stencil(amrex::Array4<amrex::Real const> const& U_arr,
                     amrex::Array4<amrex::Real const> const& rho_inv_arr,
                     const int idir,
                     const int i, const int j, const int k, const int ncomp,
                     amrex::Real* s) {

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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
add_geometric_rho_source(amrex::Array4<amrex::Real const> const& q_arr,
                         amrex::Array4<amrex::Real const> const& dloga,
                         const int i, const int j, const int k,
                         const int ncomp, amrex::Real* s) {

    using namespace reconstruction;

    // For idir == 0, i.e. r-direction:
    // this takes the form: -alpha rho u / r
    // where alpha = 1 for cylindrical and 2 for spherical
    // note dloga(idir==0) = alpha/r

    // For idir == 1,
    // i.e. theta-direction for spherical and z-direction for cylindrical
    // this takes the form: -rho v cot(theta) / r for spherical and 0 for cylindrical
    // note: dloga(idir==1) = cot(theta)/r for spherical and 0 for cylindrical.

    // ncomp should be QU for idir == 0 and QV for idir == 1

    if (ncomp == QU) {

        s[im2] += -dloga(i-2,j,k) * q_arr(i-2,j,k,QRHO) * q_arr(i-2,j,k,ncomp);
        s[im1] += -dloga(i-1,j,k) * q_arr(i-1,j,k,QRHO) * q_arr(i-1,j,k,ncomp);
        s[i0] += -dloga(i,j,k) * q_arr(i,j,k,QRHO) * q_arr(i,j,k,ncomp);
        s[ip1] += -dloga(i+1,j,k) * q_arr(i+1,j,k,QRHO) * q_arr(i+1,j,k,ncomp);
        s[ip2] += -dloga(i+2,j,k) * q_arr(i+2,j,k,QRHO) * q_arr(i+2,j,k,ncomp);

    } else if (ncomp == QV) {

        s[im2] += -dloga(i,j-2,k) * q_arr(i,j-2,k,QRHO) * q_arr(i,j-2,k,ncomp);
        s[im1] += -dloga(i,j-1,k) * q_arr(i,j-1,k,QRHO) * q_arr(i,j-1,k,ncomp);
        s[i0] += -dloga(i,j,k) * q_arr(i,j,k,QRHO) * q_arr(i,j,k,ncomp);
        s[ip1] += -dloga(i,j+1,k) * q_arr(i,j+1,k,QRHO) * q_arr(i,j+1,k,ncomp);
        s[ip2] += -dloga(i,j+2,k) * q_arr(i,j+2,k,QRHO) * q_arr(i,j+2,k,ncomp);

    }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
add_geometric_rhoe_source(amrex::Array4<amrex::Real const> const& q_arr,
                          amrex::Array4<amrex::Real const> const& dloga,
                          const int i, const int j, const int k,
                          const int ncomp, amrex::Real* s) {

    using namespace reconstruction;

    // For idir == 0, i.e. r-direction:
    // this takes the form: -alpha (rho e + p ) u / r
    // where alpha = 1 for cylindrical and 2 for spherical
    // note dloga(idir==0) = alpha/r

    // For idir == 1,
    // i.e. theta-direction for spherical and z-direction for cylindrical
    // this takes the form: - (rho e + p) v cot(theta) / r for spherical and 0 for cylindrical
    // note: dloga(idir==1) = cot(theta)/r for spherical and 0 for cylindrical.

    // ncomp should be QU for idir == 0 and QV for idir == 1

    if (ncomp == QU) {

        s[im2] += -dloga(i-2,j,k) * (q_arr(i-2,j,k,QREINT) + q_arr(i-2,j,k,QPRES)) * q_arr(i-2,j,k,ncomp);
        s[im1] += -dloga(i-1,j,k) * (q_arr(i-1,j,k,QREINT) + q_arr(i-1,j,k,QPRES)) * q_arr(i-1,j,k,ncomp);
        s[i0] += -dloga(i,j,k) * (q_arr(i,j,k,QREINT) + q_arr(i,j,k,QPRES)) * q_arr(i,j,k,ncomp);
        s[ip1] += -dloga(i+1,j,k) * (q_arr(i+1,j,k,QREINT) + q_arr(i+1,j,k,QPRES)) * q_arr(i+1,j,k,ncomp);
        s[ip2] += -dloga(i+2,j,k) * (q_arr(i+2,j,k,QREINT) + q_arr(i+2,j,k,QPRES)) * q_arr(i+2,j,k,ncomp);

    } else if (ncomp == QV) {

        s[im2] += -dloga(i,j-2,k) * (q_arr(i,j-2,k,QREINT) + q_arr(i,j-2,k,QPRES)) * q_arr(i,j-2,k,ncomp);
        s[im1] += -dloga(i,j-1,k) * (q_arr(i,j-1,k,QREINT) + q_arr(i,j-1,k,QPRES)) * q_arr(i,j-1,k,ncomp);
        s[i0] += -dloga(i,j,k) * (q_arr(i,j,k,QREINT) + q_arr(i,j,k,QPRES)) * q_arr(i,j,k,ncomp);
        s[ip1] += -dloga(i,j+1,k) * (q_arr(i,j+1,k,QREINT) + q_arr(i,j+1,k,QPRES)) * q_arr(i,j+1,k,ncomp);
        s[ip2] += -dloga(i,j+2,k) * (q_arr(i,j+2,k,QREINT) + q_arr(i,j+2,k,QPRES)) * q_arr(i,j+2,k,ncomp);

    }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
add_geometric_p_source(amrex::Array4<amrex::Real const> const& q_arr,
                       amrex::Array4<amrex::Real const> const& qaux_arr,
                       amrex::Array4<amrex::Real const> const& dloga,
                       const int i, const int j, const int k,
                       const int ncomp, amrex::Real* s) {

    using namespace reconstruction;

    // For idir == 0, i.e. r-direction:
    // this takes the form: -alpha Gamma1 p u / r
    // where alpha = 1 for cylindrical and 2 for spherical
    // note: dloga(idir==0) = alpha/r

    // For idir == 1,
    // i.e. theta-direction for spherical and z-direction for cylindrical
    // this takes the form: - Gamma1 p v cot(theta) / r for spherical and 0 for cylindrical
    // note: dloga(idir==1) = cot(theta)/r for spherical and 0 for cylindrical.

    // ncomp should be QU for idir == 0 and QV for idir == 1


    if (ncomp == QU) {

        s[im2] += -dloga(i-2,j,k) * q_arr(i-2,j,k,QPRES) * qaux_arr(i-2,j,k,QGAMC) * q_arr(i-2,j,k,ncomp);
        s[im1] += -dloga(i-1,j,k) * q_arr(i-1,j,k,QPRES) * qaux_arr(i-1,j,k,QGAMC) * q_arr(i-1,j,k,ncomp);
        s[i0] += -dloga(i,j,k) * q_arr(i,j,k,QPRES) * qaux_arr(i,j,k,QGAMC) * q_arr(i,j,k,ncomp);
        s[ip1] += -dloga(i+1,j,k) * q_arr(i+1,j,k,QPRES) * qaux_arr(i+1,j,k,QGAMC) * q_arr(i+1,j,k,ncomp);
        s[ip2] += -dloga(i+2,j,k) * q_arr(i+2,j,k,QPRES) * qaux_arr(i+2,j,k,QGAMC) * q_arr(i+2,j,k,ncomp);

    } else if (ncomp == QV) {

        s[im2] += -dloga(i,j-2,k) * q_arr(i,j-2,k,QPRES) * qaux_arr(i,j-2,k,QGAMC) * q_arr(i,j-2,k,ncomp);
        s[im1] += -dloga(i,j-1,k) * q_arr(i,j-1,k,QPRES) * qaux_arr(i,j-1,k,QGAMC) * q_arr(i,j-1,k,ncomp);
        s[i0] += -dloga(i,j,k) * q_arr(i,j,k,QPRES) * qaux_arr(i,j,k,QGAMC) * q_arr(i,j,k,ncomp);
        s[ip1] += -dloga(i,j+1,k) * q_arr(i,j+1,k,QPRES) * qaux_arr(i,j+1,k,QGAMC) * q_arr(i,j+1,k,ncomp);
        s[ip2] += -dloga(i,j+2,k) * q_arr(i,j+2,k,QPRES) * qaux_arr(i,j+2,k,QGAMC) * q_arr(i,j+2,k,ncomp);
    }

}


#endif
