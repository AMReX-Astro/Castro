#ifndef CASTRO_RECONSTRUCTION_H
#define CASTRO_RECONSTRUCTION_H

#include <state_indices.H>

namespace reconstruction {
    enum slope_indices : std::uint8_t {
        im3 = 0,
        im2,
        im1,
        i0,
        ip1,
        ip2,
        ip3,
        nslp
    };

    enum class Centering : std::uint8_t {
        zone_on_left_bndry=0,
        zone_one_from_left_bndry,
        zone_centered,
        zone_one_from_right_bndry,
        zone_on_right_bndry
    };

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
reconstruction::Centering
get_centering(const int i, const int j, const int k, const int idir,
              const amrex::GpuArray<int, 3>& domlo, const amrex::GpuArray<int, 3>& domhi,
              const bool lo_bc_test, const bool hi_bc_test) {

    // look at the location of the zone we are reconstructing and
    // return a Centering that tells us if it is near the boundary
    // (for specific boundary conditions).  These are cases where we
    // might want to use one-sided stencils.

    reconstruction::Centering centering{reconstruction::Centering::zone_centered};

    if (lo_bc_test && ((idir == 0 && i == domlo[0]) ||
                       (idir == 1 && j == domlo[1]) ||
                       (idir == 2 && k == domlo[2]))) {
        centering = reconstruction::Centering::zone_on_left_bndry;

    } else if (lo_bc_test && ((idir == 0 && i == domlo[0]+1) ||
                              (idir == 1 && j == domlo[1]+1) ||
                              (idir == 2 && k == domlo[2]+1))) {
        centering = reconstruction::Centering::zone_one_from_left_bndry;

    } else if (hi_bc_test && ((idir == 0 && i == domhi[0]-1) ||
                              (idir == 1 && j == domhi[1]-1) ||
                              (idir == 2 && k == domhi[2]-1))) {
        centering = reconstruction::Centering::zone_one_from_right_bndry;

    } else if (hi_bc_test && ((idir == 0 && i == domhi[0]) ||
                              (idir == 1 && j == domhi[1]) ||
                              (idir == 2 && k == domhi[2]))) {
        centering = reconstruction::Centering::zone_on_right_bndry;
    }

    return centering;

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
load_stencil(amrex::Array4<amrex::Real const> const& q_arr,
             reconstruction::Centering centering, const int idir,
             const int i, const int j, const int k, const int ncomp,
             amrex::Real* s) {

    using namespace reconstruction;

    // care needs to be done when filling the stencil here, since we
    // might be passing in the source for q_arr, which has
    // NUM_GROW_SRC < NUM_GROW ghost cells.  But we know that we only
    // need the one-sided stencils when we are at a physical boundary,
    // and we also know that extra zones are all inside the domain, so
    // we don't access any cells we shouldn't if we fill according to
    // the centering.

    switch (centering) {
    case Centering::zone_on_left_bndry:

        // --- | --- | --- | --- | --- | --- | --- |
        // i-3   i-2   i-1 ^  i    i+1   i+2   i+3
        //                 ^
        //              boundary

        if (idir == 0) {
            s[i0]  = q_arr(i,j,k,ncomp);
            s[ip1] = q_arr(i+1,j,k,ncomp);
            s[ip2] = q_arr(i+2,j,k,ncomp);
            s[ip3] = q_arr(i+3,j,k,ncomp);

        } else if (idir == 1) {
            s[i0]  = q_arr(i,j,k,ncomp);
            s[ip1] = q_arr(i,j+1,k,ncomp);
            s[ip2] = q_arr(i,j+2,k,ncomp);
            s[ip3] = q_arr(i,j+3,k,ncomp);

        } else {
            s[i0]  = q_arr(i,j,k,ncomp);
            s[ip1] = q_arr(i,j,k+1,ncomp);
            s[ip2] = q_arr(i,j,k+2,ncomp);
            s[ip3] = q_arr(i,j,k+3,ncomp);

        }
        break;

    case Centering::zone_one_from_left_bndry:

        // --- | --- | --- | --- | --- | --- | --- |
        // i-3   i-2 ^ i-1    i    i+1   i+2   i+3
        //           ^
        //        boundary

        if (idir == 0) {
            s[im1] = q_arr(i-1,j,k,ncomp);
            s[i0]  = q_arr(i,j,k,ncomp);
            s[ip1] = q_arr(i+1,j,k,ncomp);
            s[ip2] = q_arr(i+2,j,k,ncomp);

        } else if (idir == 1) {
            s[im1] = q_arr(i,j-1,k,ncomp);
            s[i0]  = q_arr(i,j,k,ncomp);
            s[ip1] = q_arr(i,j+1,k,ncomp);
            s[ip2] = q_arr(i,j+2,k,ncomp);

        } else {
            s[im1] = q_arr(i,j,k-1,ncomp);
            s[i0]  = q_arr(i,j,k,ncomp);
            s[ip1] = q_arr(i,j,k+1,ncomp);
            s[ip2] = q_arr(i,j,k+2,ncomp);

        }
        break;

    case Centering::zone_centered:
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
        break;

    case Centering::zone_one_from_right_bndry:

        // --- | --- | --- | --- | --- | --- | --- |
        // i-3   i-2   i-1    i    i+1 ^ i+2   i+3
        //                             ^
        //                          boundary

        if (idir == 0) {
            s[im2] = q_arr(i-2,j,k,ncomp);
            s[im1] = q_arr(i-1,j,k,ncomp);
            s[i0]  = q_arr(i,j,k,ncomp);
            s[ip1] = q_arr(i+1,j,k,ncomp);

        } else if (idir == 1) {
            s[im2] = q_arr(i,j-2,k,ncomp);
            s[im1] = q_arr(i,j-1,k,ncomp);
            s[i0]  = q_arr(i,j,k,ncomp);
            s[ip1] = q_arr(i,j+1,k,ncomp);

        } else {
            s[im2] = q_arr(i,j,k-2,ncomp);
            s[im1] = q_arr(i,j,k-1,ncomp);
            s[i0]  = q_arr(i,j,k,ncomp);
            s[ip1] = q_arr(i,j,k+1,ncomp);

        }
        break;

    case Centering::zone_on_right_bndry:

        // --- | --- | --- | --- | --- | --- | --- |
        // i-3   i-2   i-1    i  ^ i+1   i+2   i+3
        //                       ^
        //                    boundary

        if (idir == 0) {
            s[im3] = q_arr(i-3,j,k,ncomp);
            s[im2] = q_arr(i-2,j,k,ncomp);
            s[im1] = q_arr(i-1,j,k,ncomp);
            s[i0]  = q_arr(i,j,k,ncomp);

        } else if (idir == 1) {
            s[im3] = q_arr(i,j-3,k,ncomp);
            s[im2] = q_arr(i,j-2,k,ncomp);
            s[im1] = q_arr(i,j-1,k,ncomp);
            s[i0]  = q_arr(i,j,k,ncomp);

        } else {
            s[im3] = q_arr(i,j,k-3,ncomp);
            s[im2] = q_arr(i,j,k-2,ncomp);
            s[im1] = q_arr(i,j,k-1,ncomp);
            s[i0]  = q_arr(i,j,k,ncomp);

        }
        break;
    }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
load_passive_stencil(amrex::Array4<amrex::Real const> const& U_arr,
                     amrex::Array4<amrex::Real const> const& rho_inv_arr,
                     reconstruction::Centering centering, const int idir,
                     const int i, const int j, const int k, const int ncomp,
                     amrex::Real* s) {

    using namespace reconstruction;

    switch (centering) {
    case Centering::zone_on_left_bndry:

        if (idir == 0) {
            s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);
            s[ip1] = U_arr(i+1,j,k,ncomp) * rho_inv_arr(i+1,j,k);
            s[ip2] = U_arr(i+2,j,k,ncomp) * rho_inv_arr(i+2,j,k);
            s[ip3] = U_arr(i+3,j,k,ncomp) * rho_inv_arr(i+3,j,k);

        } else if (idir == 1) {
            s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);
            s[ip1] = U_arr(i,j+1,k,ncomp) * rho_inv_arr(i,j+1,k);
            s[ip2] = U_arr(i,j+2,k,ncomp) * rho_inv_arr(i,j+2,k);
            s[ip3] = U_arr(i,j+3,k,ncomp) * rho_inv_arr(i,j+3,k);

        } else {
            s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);
            s[ip1] = U_arr(i,j,k+1,ncomp) * rho_inv_arr(i,j,k+1);
            s[ip2] = U_arr(i,j,k+2,ncomp) * rho_inv_arr(i,j,k+2);
            s[ip3] = U_arr(i,j,k+3,ncomp) * rho_inv_arr(i,j,k+3);

        }
        break;

    case Centering::zone_one_from_left_bndry:

        if (idir == 0) {
            s[im1] = U_arr(i-1,j,k,ncomp) * rho_inv_arr(i-1,j,k);
            s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);
            s[ip1] = U_arr(i+1,j,k,ncomp) * rho_inv_arr(i+1,j,k);
            s[ip2] = U_arr(i+2,j,k,ncomp) * rho_inv_arr(i+2,j,k);

        } else if (idir == 1) {
            s[im1] = U_arr(i,j-1,k,ncomp) * rho_inv_arr(i,j-1,k);
            s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);
            s[ip1] = U_arr(i,j+1,k,ncomp) * rho_inv_arr(i,j+1,k);
            s[ip2] = U_arr(i,j+2,k,ncomp) * rho_inv_arr(i,j+2,k);

        } else {
            s[im1] = U_arr(i,j,k-1,ncomp) * rho_inv_arr(i,j,k-1);
            s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);
            s[ip1] = U_arr(i,j,k+1,ncomp) * rho_inv_arr(i,j,k+1);
            s[ip2] = U_arr(i,j,k+2,ncomp) * rho_inv_arr(i,j,k+2);

        }
        break;

    case Centering::zone_centered:

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
        break;

    case Centering::zone_one_from_right_bndry:

        if (idir == 0) {
            s[im2] = U_arr(i-2,j,k,ncomp) * rho_inv_arr(i-2,j,k);
            s[im1] = U_arr(i-1,j,k,ncomp) * rho_inv_arr(i-1,j,k);
            s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);
            s[ip1] = U_arr(i+1,j,k,ncomp) * rho_inv_arr(i+1,j,k);

        } else if (idir == 1) {
            s[im2] = U_arr(i,j-2,k,ncomp) * rho_inv_arr(i,j-2,k);
            s[im1] = U_arr(i,j-1,k,ncomp) * rho_inv_arr(i,j-1,k);
            s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);
            s[ip1] = U_arr(i,j+1,k,ncomp) * rho_inv_arr(i,j+1,k);

        } else {
            s[im2] = U_arr(i,j,k-2,ncomp) * rho_inv_arr(i,j,k-2);
            s[im1] = U_arr(i,j,k-1,ncomp) * rho_inv_arr(i,j,k-1);
            s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);
            s[ip1] = U_arr(i,j,k+1,ncomp) * rho_inv_arr(i,j,k+1);

        }
        break;

    case Centering::zone_on_right_bndry:

        if (idir == 0) {
            s[im3] = U_arr(i-3,j,k,ncomp) * rho_inv_arr(i-3,j,k);
            s[im2] = U_arr(i-2,j,k,ncomp) * rho_inv_arr(i-2,j,k);
            s[im1] = U_arr(i-1,j,k,ncomp) * rho_inv_arr(i-1,j,k);
            s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);

        } else if (idir == 1) {
            s[im3] = U_arr(i,j-3,k,ncomp) * rho_inv_arr(i,j-3,k);
            s[im2] = U_arr(i,j-2,k,ncomp) * rho_inv_arr(i,j-2,k);
            s[im1] = U_arr(i,j-1,k,ncomp) * rho_inv_arr(i,j-1,k);
            s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);

        } else {
            s[im3] = U_arr(i,j,k-3,ncomp) * rho_inv_arr(i,j,k-3);
            s[im2] = U_arr(i,j,k-2,ncomp) * rho_inv_arr(i,j,k-2);
            s[im1] = U_arr(i,j,k-1,ncomp) * rho_inv_arr(i,j,k-1);
            s[i0]  = U_arr(i,j,k,ncomp) * rho_inv_arr(i,j,k);

        }
        break;
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

    // note, we know that dloga and q_arr all have NUM_GROW ghost
    // cells, so we can fill the entire stencil and not worry about
    // the centering.

    if (ncomp == QU) {

        s[im3] += -dloga(i-3,j,k) * q_arr(i-3,j,k,QRHO) * q_arr(i-3,j,k,ncomp);
        s[im2] += -dloga(i-2,j,k) * q_arr(i-2,j,k,QRHO) * q_arr(i-2,j,k,ncomp);
        s[im1] += -dloga(i-1,j,k) * q_arr(i-1,j,k,QRHO) * q_arr(i-1,j,k,ncomp);
        s[i0] += -dloga(i,j,k) * q_arr(i,j,k,QRHO) * q_arr(i,j,k,ncomp);
        s[ip1] += -dloga(i+1,j,k) * q_arr(i+1,j,k,QRHO) * q_arr(i+1,j,k,ncomp);
        s[ip2] += -dloga(i+2,j,k) * q_arr(i+2,j,k,QRHO) * q_arr(i+2,j,k,ncomp);
        s[ip3] += -dloga(i+3,j,k) * q_arr(i+3,j,k,QRHO) * q_arr(i+3,j,k,ncomp);

    } else if (ncomp == QV) {

        s[im3] += -dloga(i,j-3,k) * q_arr(i,j-3,k,QRHO) * q_arr(i,j-3,k,ncomp);
        s[im2] += -dloga(i,j-2,k) * q_arr(i,j-2,k,QRHO) * q_arr(i,j-2,k,ncomp);
        s[im1] += -dloga(i,j-1,k) * q_arr(i,j-1,k,QRHO) * q_arr(i,j-1,k,ncomp);
        s[i0] += -dloga(i,j,k) * q_arr(i,j,k,QRHO) * q_arr(i,j,k,ncomp);
        s[ip1] += -dloga(i,j+1,k) * q_arr(i,j+1,k,QRHO) * q_arr(i,j+1,k,ncomp);
        s[ip2] += -dloga(i,j+2,k) * q_arr(i,j+2,k,QRHO) * q_arr(i,j+2,k,ncomp);
        s[ip3] += -dloga(i,j+3,k) * q_arr(i,j+3,k,QRHO) * q_arr(i,j+3,k,ncomp);

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

    // note, we know that dloga and q_arr all have NUM_GROW ghost
    // cells, so we can fill the entire stencil and not worry about
    // the centering.

    if (ncomp == QU) {

        s[im3] += -dloga(i-3,j,k) * (q_arr(i-3,j,k,QREINT) + q_arr(i-3,j,k,QPRES)) * q_arr(i-3,j,k,ncomp);
        s[im2] += -dloga(i-2,j,k) * (q_arr(i-2,j,k,QREINT) + q_arr(i-2,j,k,QPRES)) * q_arr(i-2,j,k,ncomp);
        s[im1] += -dloga(i-1,j,k) * (q_arr(i-1,j,k,QREINT) + q_arr(i-1,j,k,QPRES)) * q_arr(i-1,j,k,ncomp);
        s[i0] += -dloga(i,j,k) * (q_arr(i,j,k,QREINT) + q_arr(i,j,k,QPRES)) * q_arr(i,j,k,ncomp);
        s[ip1] += -dloga(i+1,j,k) * (q_arr(i+1,j,k,QREINT) + q_arr(i+1,j,k,QPRES)) * q_arr(i+1,j,k,ncomp);
        s[ip2] += -dloga(i+2,j,k) * (q_arr(i+2,j,k,QREINT) + q_arr(i+2,j,k,QPRES)) * q_arr(i+2,j,k,ncomp);
        s[ip3] += -dloga(i+3,j,k) * (q_arr(i+3,j,k,QREINT) + q_arr(i+3,j,k,QPRES)) * q_arr(i+3,j,k,ncomp);

    } else if (ncomp == QV) {

        s[im3] += -dloga(i,j-3,k) * (q_arr(i,j-3,k,QREINT) + q_arr(i,j-3,k,QPRES)) * q_arr(i,j-3,k,ncomp);
        s[im2] += -dloga(i,j-2,k) * (q_arr(i,j-2,k,QREINT) + q_arr(i,j-2,k,QPRES)) * q_arr(i,j-2,k,ncomp);
        s[im1] += -dloga(i,j-1,k) * (q_arr(i,j-1,k,QREINT) + q_arr(i,j-1,k,QPRES)) * q_arr(i,j-1,k,ncomp);
        s[i0] += -dloga(i,j,k) * (q_arr(i,j,k,QREINT) + q_arr(i,j,k,QPRES)) * q_arr(i,j,k,ncomp);
        s[ip1] += -dloga(i,j+1,k) * (q_arr(i,j+1,k,QREINT) + q_arr(i,j+1,k,QPRES)) * q_arr(i,j+1,k,ncomp);
        s[ip2] += -dloga(i,j+2,k) * (q_arr(i,j+2,k,QREINT) + q_arr(i,j+2,k,QPRES)) * q_arr(i,j+2,k,ncomp);
        s[ip3] += -dloga(i,j+3,k) * (q_arr(i,j+3,k,QREINT) + q_arr(i,j+3,k,QPRES)) * q_arr(i,j+3,k,ncomp);

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

    // note, we know that dloga, q_arr, and qaux_arr all have NUM_GROW
    // ghost cells, so we can fill the entire stencil and not worry
    // about the centering.

    if (ncomp == QU) {

        s[im3] += -dloga(i-3,j,k) * q_arr(i-3,j,k,QPRES) * qaux_arr(i-3,j,k,QGAMC) * q_arr(i-3,j,k,ncomp);
        s[im2] += -dloga(i-2,j,k) * q_arr(i-2,j,k,QPRES) * qaux_arr(i-2,j,k,QGAMC) * q_arr(i-2,j,k,ncomp);
        s[im1] += -dloga(i-1,j,k) * q_arr(i-1,j,k,QPRES) * qaux_arr(i-1,j,k,QGAMC) * q_arr(i-1,j,k,ncomp);
        s[i0] += -dloga(i,j,k) * q_arr(i,j,k,QPRES) * qaux_arr(i,j,k,QGAMC) * q_arr(i,j,k,ncomp);
        s[ip1] += -dloga(i+1,j,k) * q_arr(i+1,j,k,QPRES) * qaux_arr(i+1,j,k,QGAMC) * q_arr(i+1,j,k,ncomp);
        s[ip2] += -dloga(i+2,j,k) * q_arr(i+2,j,k,QPRES) * qaux_arr(i+2,j,k,QGAMC) * q_arr(i+2,j,k,ncomp);
        s[ip3] += -dloga(i+3,j,k) * q_arr(i+3,j,k,QPRES) * qaux_arr(i+3,j,k,QGAMC) * q_arr(i+3,j,k,ncomp);

    } else if (ncomp == QV) {

        s[im2] += -dloga(i,j-3,k) * q_arr(i,j-3,k,QPRES) * qaux_arr(i,j-3,k,QGAMC) * q_arr(i,j-3,k,ncomp);
        s[im2] += -dloga(i,j-2,k) * q_arr(i,j-2,k,QPRES) * qaux_arr(i,j-2,k,QGAMC) * q_arr(i,j-2,k,ncomp);
        s[im1] += -dloga(i,j-1,k) * q_arr(i,j-1,k,QPRES) * qaux_arr(i,j-1,k,QGAMC) * q_arr(i,j-1,k,ncomp);
        s[i0] += -dloga(i,j,k) * q_arr(i,j,k,QPRES) * qaux_arr(i,j,k,QGAMC) * q_arr(i,j,k,ncomp);
        s[ip1] += -dloga(i,j+1,k) * q_arr(i,j+1,k,QPRES) * qaux_arr(i,j+1,k,QGAMC) * q_arr(i,j+1,k,ncomp);
        s[ip2] += -dloga(i,j+2,k) * q_arr(i,j+2,k,QPRES) * qaux_arr(i,j+2,k,QGAMC) * q_arr(i,j+2,k,ncomp);
        s[ip3] += -dloga(i,j+3,k) * q_arr(i,j+3,k,QPRES) * qaux_arr(i,j+3,k,QGAMC) * q_arr(i,j+3,k,ncomp);
    }

}


#endif
