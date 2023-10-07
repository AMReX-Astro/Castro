#ifndef filter_H
#define filter_H

// reference: R. J. Purser (J. of Clim. and Apld. Meteorology, 1987)

namespace filter
{
    // 3-point filter: T=1, R=S=0
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff1 (int i)
    {
        if (i == 0) {
            return 0.50_rt;
        }
        else if (i == 1) {
            return 0.25_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    // For boundary cell
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff1b (int i)
    {
        if (i == 0) {
            return 0.75_rt;
        }
        else if (i == 1) {
            return 0.25_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    // 5-point filter: T=2, R+S+1=T
    // S = 0, 1
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff2 (int i, int j)
    {
        if (i == 0 && j == 0) {
            return 0.625_rt;
        }
        else if (i == 1 && j == 0) {
            return 0.25_rt;
        }
        else if (i == 2 && j == 0) {
            return -0.0625_rt;
        }
        else if (i == 0 && j == 1) {
            return 0.375_rt;
        }
        else if (i == 1 && j == 1) {
            return 0.25_rt;
        }
        else if (i == 2 && j == 1) {
            return 0.0625_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    // For first boundary cell
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff2b0 (int i)
    {
        if (i == 0) {
            return 17.0_rt / 16.0_rt;
        }
        else if (i == 1) {
            return -2.0_rt / 16.0_rt;
        }
        else if (i == 2) {
            return 1.0_rt / 16.0_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    // For second boundary cell
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff2b1 (int i)
    {
        if (i == -1) {
            return -2.0_rt / 16.0_rt;
        }
        else if (i == 0) {
            return 21.0_rt / 16.0_rt;
        }
        else if (i == 1) {
            return -4.0_rt / 16.0_rt;
        }
        else if (i == 2) {
            return 1.0_rt / 16.0_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    // 7-point filter: T=3, R+S+1=T
    // S = 0, 1, 2
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff3 (int i, int j)
    {
        if (i == 0 && j == 0) {
            return 44.0_rt / 64.0_rt;
        }
        else if (i == 1 && j == 0) {
            return 15.0_rt / 64.0_rt;
        }
        else if (i == 2 && j == 0) {
            return -6.0_rt / 64.0_rt;
        }
        else if (i == 3 && j == 0) {
            return 1.0_rt / 64.0_rt;
        }
        else if (i == 0 && j == 1) {
            return 32.0_rt / 64.0_rt;
        }
        else if (i == 1 && j == 1) {
            return 18.0_rt / 64.0_rt;
        }
        else if (i == 2 && j == 1) {
            return 0.e0_rt;
        }
        else if (i == 3 && j == 1) {
            return -2.0_rt / 64.0_rt;
        }
        else if (i == 0 && j == 2) {
            return 20.0_rt / 64.0_rt;
        }
        else if (i == 1 && j == 2) {
            return 15.0_rt / 64.0_rt;
        }
        else if (i == 2 && j == 2) {
            return 6.0_rt / 64.0_rt;
        }
        else if (i == 3 && j == 2) {
            return 1.0_rt / 64.0_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    // For boundary cells
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff3b0 (int i)
    {
        if (i == 0) {
            return 63.0_rt / 64.0_rt;
        }
        else if (i == 1) {
            return 3.0_rt / 64.0_rt;
        }
        else if (i == 2) {
            return -3.0_rt / 64.0_rt;
        }
        else if (i == 3) {
            return 1.0_rt / 64.0_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff3b1 (int i)
    {
        if (i == -1) {
            return 3.0_rt / 64.0_rt;
        }
        else if (i == 0) {
            return 54.0_rt / 64.0_rt;
        }
        else if (i == 1) {
            return 12.0_rt / 64.0_rt;
        }
        else if (i == 2) {
            return -6.e0_rt/64.e0_rt;
        }
        else if (i == 3) {
            return 1.0_rt / 64.0_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff3b2 (int i)
    {
        if (i == -2) {
            return -3.0_rt / 64.0_rt;
        }
        else if (i == -1) {
            return 12.0_rt / 64.0_rt;
        }
        else if (i == 0) {
            return 45.0_rt / 64.0_rt;
        }
        else if (i == 1) {
            return 15.0_rt / 64.0_rt;
        }
        else if (i == 2) {
            return -6.0_rt / 64.0_rt;
        }
        else if (i == 3) {
            return 1.0_rt / 64.0_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    // 9-point filter: T=4, R+S+1=T
    // S = 0, 1, 2, 3
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff4 (int i, int j)
    {
        if (i == 0 && j == 0) {
            return 186.0_rt / 256.0_rt;
        }
        else if (i == 1 && j == 0) {
            return 56.0_rt / 256.0_rt;
        }
        else if (i == 2 && j == 0) {
            return -28.0_rt / 256.0_rt;
        }
        else if (i == 3 && j == 0) {
            return 8.0_rt / 256.0_rt;
        }
        else if (i == 4 && j == 0) {
            return -1.0_rt / 256.0_rt;
        }
        else if (i == 0 && j == 1) {
            return 146.0_rt / 256.0_rt;
        }
        else if (i == 1 && j == 1) {
            return 72.0_rt / 256.0_rt;
        }
        else if (i == 2 && j == 1) {
            return -12.0_rt / 256.0_rt;
        }
        else if (i == 3 && j == 1) {
            return -8.0_rt / 256.0_rt;
        }
        else if (i == 4 && j == 1) {
            return 3.0_rt / 256.0_rt;
        }
        else if (i == 0 && j == 2) {
            return 110.0_rt / 256.0_rt;
        }
        else if (i == 1 && j == 2) {
            return 72.0_rt / 256.0_rt;
        }
        else if (i == 2 && j == 2) {
            return 12.0_rt / 256.0_rt;
        }
        else if (i == 3 && j == 2) {
            return -8.0_rt / 256.0_rt;
        }
        else if (i == 4 && j == 2) {
            return -3.0_rt / 256.0_rt;
        }
        else if (i == 0 && j == 3) {
            return 70.0_rt / 256.0_rt;
        }
        else if (i == 1 && j == 3) {
            return 56.0_rt / 256.0_rt;
        }
        else if (i == 2 && j == 3) {
            return 28.0_rt / 256.0_rt;
        }
        else if (i == 3 && j == 3) {
            return 8.0_rt / 256.0_rt;
        }
        else if (i == 4 && j == 3) {
            return 1.0_rt / 256.0_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    // For boundary cells
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff4b0 (int i)
    {
        if (i == 0) {
            return 257.0_rt / 256.0_rt;
        }
        else if (i == 1) {
            return -4.0_rt / 256.0_rt;
        }
        else if (i == 2) {
            return 6.0_rt / 256.0_rt;
        }
        else if (i == 3) {
            return -4.0_rt / 256.0_rt;
        }
        else if (i == 4) {
            return 1.0_rt / 256.0_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff4b1 (int i)
    {
        if (i == -1) {
            return -4.0_rt / 256.0_rt;
        }
        else if (i == 0) {
            return 273.0_rt / 256.0_rt;
        }
        else if (i == 1) {
            return -28.0_rt / 256.0_rt;
        }
        else if (i == 2) {
            return 22.0_rt / 256.0_rt;
        }
        else if (i == 3) {
            return -8.0_rt / 256.0_rt;
        }
        else if (i == 4) {
            return 1.0_rt / 256.0_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff4b2 (int i)
    {
        if (i == -2) {
            return 6.0_rt / 256.0_rt;
        }
        else if (i == -1) {
            return -28.0_rt / 256.0_rt;
        }
        else if (i == 0) {
            return 309.0_rt / 256.0_rt;
        }
        else if (i == 1) {
            return -52.0_rt / 256.0_rt;
        }
        else if (i == 2) {
            return 28.0_rt / 256.0_rt;
        }
        else if (i == 3) {
            return -8.0_rt / 256.0_rt;
        }
        else if (i == 4) {
            return 1.0_rt / 256.0_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real ff4b3 (int i)
    {
        if (i == -3) {
            return -4.0_rt / 256.0_rt;
        }
        else if (i == -2) {
            return 22.0_rt / 256.0_rt;
        }
        else if (i == -1) {
            return -52.0_rt / 256.0_rt;
        }
        else if (i == 0) {
            return 325.0_rt / 256.0_rt;
        }
        else if (i == 1) {
            return -56.0_rt / 256.0_rt;
        }
        else if (i == 2) {
            return 28.0_rt / 256.0_rt;
        }
        else if (i == 3) {
            return -8.0_rt / 256.0_rt;
        }
        else if (i == 4) {
            return 1.0_rt / 256.0_rt;
        }
        else {
            return 0.0_rt;
        }
    }

    template<int width>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void load_filter_stencil(Array4<Real const> const input,
                             Array1D<Real, -width, width>& stencil,
                             int dir, int i, int j, int k, int g)
    {
        stencil(0) = input(i,j,k,g);

        if constexpr (width >= 1) {
            if (dir == 0) {
                stencil(-1) = input(i-1,j,k,g);
                stencil( 1) = input(i+1,j,k,g);
            }
            else if (dir == 1) {
                stencil(-1) = input(i,j-1,k,g);
                stencil( 1) = input(i,j+1,k,g);
            }
            else {
                stencil(-1) = input(i,j,k-1,g);
                stencil( 1) = input(i,j,k+1,g);
            }
        }

        if constexpr (width >= 2) {
            if (dir == 0) {
                stencil(-2) = input(i-2,j,k,g);
                stencil( 2) = input(i+2,j,k,g);
            }
            else if (dir == 1) {
                stencil(-2) = input(i,j-2,k,g);
                stencil( 2) = input(i,j+2,k,g);
            }
            else {
                stencil(-2) = input(i,j,k-2,g);
                stencil( 2) = input(i,j,k+2,g);
            }
        }

        if constexpr (width >= 3) {
            if (dir == 0) {
                stencil(-3) = input(i-3,j,k,g);
                stencil( 3) = input(i+3,j,k,g);
            }
            else if (dir == 1) {
                stencil(-3) = input(i,j-3,k,g);
                stencil( 3) = input(i,j+3,k,g);
            }
            else {
                stencil(-3) = input(i,j,k-3,g);
                stencil( 3) = input(i,j,k+3,g);
            }
        }

        if constexpr (width >= 4) {
            if (dir == 0) {
                stencil(-4) = input(i-4,j,k,g);
                stencil( 4) = input(i+4,j,k,g);
            }
            else if (dir == 1) {
                stencil(-4) = input(i,j-4,k,g);
                stencil( 4) = input(i,j+4,k,g);
            }
            else {
                stencil(-4) = input(i,j,k-4,g);
                stencil( 4) = input(i,j,k+4,g);
            }
        }
    }

    template<int width>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real apply_filter (Array4<Real const> const Er,
                       Array4<Real const> const lam,
                       int dir, int S, int i, int j, int k, int g)
    {
        Real lamfil;

        // In whatever direction we're filtering, retrieve a stencil for both
        // lambda and the radiation energy, Er. The energy does not directly
        // factor into the calculation, but helps us understand when we have
        // gone over the edge of the grid. In the code below, we'll look at
        // the stencil only for those zones where Er != -1, which constitutes
        // our valid data. The approach we will use is to look progressively
        // one zone further out to the left or right; if we detect that any of
        // the zones we would normally use for our centered stencil are an
        // invalid zone, we use a one-sided or skewed stencil with only zones
        // that are available. Note that this relies on the assumption that if
        // Er(-1) is -1, then all zones to the left of that are also -1, so we
        // don't bother actually examining the contents of those zones (and
        // the analogous assumption holds on the right edge of the grid).

        Array1D<Real, -width, width> Er_stencil, lam_stencil;

        load_filter_stencil<width>(Er, Er_stencil, dir, i, j, k, g);
        load_filter_stencil<width>(lam, lam_stencil, dir, i, j, k, g);

        if constexpr (width == 1) {
            if (Er_stencil(-1) == -1.0_rt) {
                lamfil = ff1b(0) * lam_stencil(0) +
                         ff1b(1) * lam_stencil(1);
            }
            else if (Er_stencil(1) == -1.0_rt) {
                lamfil = ff1b(1) * lam_stencil(-1) +
                         ff1b(0) * lam_stencil( 0);
            }
            else {
                lamfil = ff1(0) * lam_stencil(0) +
                         ff1(1) * (lam_stencil(-1) + lam_stencil(1));
            }
        }
        else if constexpr (width == 2) {
            if (Er_stencil(-1) == -1.0_rt) {
                lamfil = ff2b0(0) * lam_stencil(0) +
                         ff2b0(1) * lam_stencil(1) +
                         ff2b0(2) * lam_stencil(2);
            }
            else if (Er_stencil(-2) == -1.0_rt) {
                lamfil = ff2b1(-1) * lam_stencil(-1) +
                         ff2b1( 0) * lam_stencil( 0) +
                         ff2b1( 1) * lam_stencil( 1) +
                         ff2b1( 2) * lam_stencil( 2);
            }
            else if (Er_stencil(1) == -1.0_rt) {
                lamfil = ff2b0(2) * lam_stencil(-2) +
                         ff2b0(1) * lam_stencil(-1) +
                         ff2b0(0) * lam_stencil( 0);
            }
            else if (Er_stencil(2) == -1.0_rt) {
                lamfil = ff2b1( 2) * lam_stencil(-2) +
                         ff2b1( 1) * lam_stencil(-1) +
                         ff2b1( 0) * lam_stencil( 0) +
                         ff2b1(-1) * lam_stencil( 1);
            }
            else {
                lamfil = ff2(0,S) * lam_stencil(0) +
                         ff2(1,S) * (lam_stencil(-1) + lam_stencil(1)) +
                         ff2(2,S) * (lam_stencil(-2) + lam_stencil(2));
            }
        }
        else if constexpr (width == 3) {
            if (Er_stencil(-1) == -1.0_rt) {
                lamfil = ff3b0(0) * lam_stencil(0) +
                         ff3b0(1) * lam_stencil(1) +
                         ff3b0(2) * lam_stencil(2) +
                         ff3b0(3) * lam_stencil(3);
            }
            else if (Er_stencil(-2) == -1.0_rt) {
                lamfil = ff3b1(-1) * lam_stencil(-1) +
                         ff3b1( 0) * lam_stencil( 0) +
                         ff3b1( 1) * lam_stencil( 1) +
                         ff3b1( 2) * lam_stencil( 2) +
                         ff3b1( 3) * lam_stencil( 3);
            }
            else if (Er_stencil(-3) == -1.0_rt) {
                lamfil = ff3b2(-2) * lam_stencil(-2) +
                         ff3b2(-1) * lam_stencil(-1) +
                         ff3b2( 0) * lam_stencil( 0) +
                         ff3b2( 1) * lam_stencil( 1) +
                         ff3b2( 2) * lam_stencil( 2) +
                         ff3b2( 3) * lam_stencil( 3);
            }
            else if (Er_stencil(1) == -1.0_rt) {
                lamfil = ff3b0(3) * lam_stencil(-3) +
                         ff3b0(2) * lam_stencil(-2) +
                         ff3b0(1) * lam_stencil(-1) +
                         ff3b0(0) * lam_stencil( 0);
            }
            else if (Er_stencil(2) == -1.0_rt) {
                lamfil = ff3b1( 3) * lam_stencil(-3) +
                         ff3b1( 2) * lam_stencil(-2) +
                         ff3b1( 1) * lam_stencil(-1) +
                         ff3b1( 0) * lam_stencil( 0) +
                         ff3b1(-1) * lam_stencil( 1);
            }
            else if (Er_stencil(3) == -1.0_rt) {
                lamfil = ff3b2( 3) * lam_stencil(-3) +
                         ff3b2( 2) * lam_stencil(-2) +
                         ff3b2( 1) * lam_stencil(-1) +
                         ff3b2( 0) * lam_stencil( 0) +
                         ff3b2(-1) * lam_stencil( 1) +
                         ff3b2(-2) * lam_stencil( 2);
            }
            else {
                lamfil = ff3(0,S) * lam_stencil(0) +
                         ff3(1,S) * (lam_stencil(-1) + lam_stencil(1)) +
                         ff3(2,S) * (lam_stencil(-2) + lam_stencil(2)) +
                         ff3(3,S) * (lam_stencil(-3) + lam_stencil(3));
            }
        }
        else if constexpr (width == 4) {
            if (Er_stencil(-1) == -1.0_rt) {
                lamfil = ff4b0(0) * lam_stencil(0) +
                         ff4b0(1) * lam_stencil(1) +
                         ff4b0(2) * lam_stencil(2) +
                         ff4b0(3) * lam_stencil(3) +
                         ff4b0(4) * lam_stencil(4);
            }
            else if (Er_stencil(-2) == -1.0_rt) {
                lamfil = ff4b1(-1) * lam_stencil(-1) +
                         ff4b1( 0) * lam_stencil( 0) +
                         ff4b1( 1) * lam_stencil( 1) +
                         ff4b1( 2) * lam_stencil( 2) +
                         ff4b1( 3) * lam_stencil( 3) +
                         ff4b1( 4) * lam_stencil( 4);
            }
            else if (Er_stencil(-3) == -1.0_rt) {
                lamfil = ff4b2(-2) * lam_stencil(-2) +
                         ff4b2(-1) * lam_stencil(-1) +
                         ff4b2( 0) * lam_stencil( 0) +
                         ff4b2( 1) * lam_stencil( 1) +
                         ff4b2( 2) * lam_stencil( 2) +
                         ff4b2( 3) * lam_stencil( 3) +
                         ff4b2( 4) * lam_stencil( 4);
            }
            else if (Er_stencil(-4) == -1.0_rt) {
                lamfil = ff4b3(-3) * lam_stencil(-3) +
                         ff4b3(-2) * lam_stencil(-2) +
                         ff4b3(-1) * lam_stencil(-1) +
                         ff4b3( 0) * lam_stencil( 0) +
                         ff4b3( 1) * lam_stencil( 1) +
                         ff4b3( 2) * lam_stencil( 2) +
                         ff4b3( 3) * lam_stencil( 3) +
                         ff4b3( 4) * lam_stencil( 4);
            }
            else if (Er_stencil(1) == -1.0_rt) {
                lamfil = ff4b0(4) * lam_stencil(-4) +
                         ff4b0(3) * lam_stencil(-3) +
                         ff4b0(2) * lam_stencil(-2) +
                         ff4b0(1) * lam_stencil(-1) +
                         ff4b0(0) * lam_stencil( 0);
            }
            else if (Er_stencil(2) == -1.0_rt) {
                lamfil = ff4b1( 4) * lam_stencil(-4) +
                         ff4b1( 3) * lam_stencil(-3) +
                         ff4b1( 2) * lam_stencil(-2) +
                         ff4b1( 1) * lam_stencil(-1) +
                         ff4b1( 0) * lam_stencil( 0) +
                         ff4b1(-1) * lam_stencil( 1);
            }
            else if (Er_stencil(3) == -1.0_rt) {
                lamfil = ff4b2( 4) * lam_stencil(-4) +
                         ff4b2( 3) * lam_stencil(-3) +
                         ff4b2( 2) * lam_stencil(-2) +
                         ff4b2( 1) * lam_stencil(-1) +
                         ff4b2( 0) * lam_stencil( 0) +
                         ff4b2(-1) * lam_stencil( 1) +
                         ff4b2(-2) * lam_stencil( 2);
            }
            else if (Er_stencil(4) == -1.0_rt) {
                lamfil = ff4b3( 4) * lam_stencil(-4) +
                         ff4b3( 3) * lam_stencil(-3) +
                         ff4b3( 2) * lam_stencil(-2) +
                         ff4b3( 1) * lam_stencil(-1) +
                         ff4b3( 0) * lam_stencil( 0) +
                         ff4b3(-1) * lam_stencil( 1) +
                         ff4b3(-2) * lam_stencil( 2) +
                         ff4b3(-3) * lam_stencil( 3);
            }
            else {
                lamfil = ff4(0,S) * lam_stencil(0) +
                         ff4(1,S) * (lam_stencil(-1) + lam_stencil(1)) +
                         ff4(2,S) * (lam_stencil(-2) + lam_stencil(2)) +
                         ff4(3,S) * (lam_stencil(-3) + lam_stencil(3)) +
                         ff4(4,S) * (lam_stencil(-4) + lam_stencil(4));
            }
        }

        lamfil = std::min(1.0_rt / 3.0_rt, std::max(1.0e-25_rt, lamfil));

        return lamfil;
    }
}

#endif
