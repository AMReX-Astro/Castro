#ifndef FOURTH_CENTER_AVERAGE_H
#define FOURTH_CENTER_AVERAGE_H

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
compute_laplacian(int i, int j, int k, int n,
                  Array4<Real const> const& a,
                  GpuArray<bool, AMREX_SPACEDIM> const& lo_periodic, GpuArray<bool, AMREX_SPACEDIM> const& hi_periodic,
                  GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi) {

  Real lapx = 0.0;
  Real lapy = 0.0;
  Real lapz = 0.0;

  // we use 2nd-order accurate one-sided stencils at the physical
  // boundaries note: this differs from the suggestion in MC2011 --
  // they just suggest using the Laplacian from +1 off the interior.
  // I like the one-sided better.

  if (i == domlo[0] && !lo_periodic[0]) {
    lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i+1,j,k,n) + 4.0_rt*a(i+2,j,k,n) - a(i+3,j,k,n);

  } else if (i == domhi[0] && !hi_periodic[0]) {
    lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i-1,j,k,n) + 4.0_rt*a(i-2,j,k,n) - a(i-3,j,k,n);

  } else {
    lapx = a(i+1,j,k,n) - 2.0_rt*a(i,j,k,n) + a(i-1,j,k,n);
  }

#if AMREX_SPACEDIM >= 2
  if (j == domlo[1] && !lo_periodic[1]) {
    lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j+1,k,n) + 4.0_rt*a(i,j+2,k,n) - a(i,j+3,k,n);

  } else if (j == domhi[1] && !hi_periodic[1]) {
    lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j-1,k,n) + 4.0_rt*a(i,j-2,k,n) - a(i,j-3,k,n);

  } else {
    lapy = a(i,j+1,k,n) - 2.0_rt*a(i,j,k,n) + a(i,j-1,k,n);
  }
#endif

#if AMREX_SPACEDIM == 3
  if (k == domlo[2] && !lo_periodic[2]) {
    lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k+1,n) + 4.0_rt*a(i,j,k+2,n) - a(i,j,k+3,n);

  } else if (k == domhi[2] && !hi_periodic[2]) {
    lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k-1,n) + 4.0_rt*a(i,j,k-2,n) - a(i,j,k-3,n);

  } else {
    lapz = a(i,j,k+1,n) - 2.0_rt*a(i,j,k,n) + a(i,j,k-1,n);
  }
#endif

  return lapx + lapy + lapz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
trans_laplacian(int i, int j, int k, int n,
                int idir,
                Array4<Real const> const& a,
                GpuArray<bool, AMREX_SPACEDIM> const& lo_periodic, GpuArray<bool, AMREX_SPACEDIM> const& hi_periodic,
                GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi) {

  Real lapx = 0.0;
  Real lapy = 0.0;
  Real lapz = 0.0;

  // we use 2nd-order accurate one-sided stencils at the physical boundaries
  if (idir != 0) {

    if (i == domlo[0] && !lo_periodic[0]) {
      lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i+1,j,k,n) + 4.0_rt*a(i+2,j,k,n) - a(i+3,j,k,n);

    } else if (i == domhi[0] && !hi_periodic[0]) {
      lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i-1,j,k,n) + 4.0_rt*a(i-2,j,k,n) - a(i-3,j,k,n);

    } else {
      lapx = a(i+1,j,k,n) - 2.0_rt*a(i,j,k,n) + a(i-1,j,k,n);
    }
  }

  if (idir != 1) {

    if (j == domlo[1] && !lo_periodic[1]) {
      lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j+1,k,n) + 4.0_rt*a(i,j+2,k,n) - a(i,j+3,k,n);

    } else if (j == domhi[1] && !hi_periodic[1]) {
      lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j-1,k,n) + 4.0_rt*a(i,j-2,k,n) - a(i,j-3,k,n);

    } else {
      lapy = a(i,j+1,k,n) - 2.0_rt*a(i,j,k,n) + a(i,j-1,k,n);
    }
  }

#if AMREX_SPACEDIM == 3
  if (idir != 2) {

    if (k == domlo[2] && !lo_periodic[2]) {
      lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k+1,n) + 4.0_rt*a(i,j,k+2,n) - a(i,j,k+3,n);

    } else if (k == domhi[2] && !hi_periodic[2]) {
      lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k-1,n) + 4.0_rt*a(i,j,k-2,n) - a(i,j,k-3,n);

    } else {
      lapz = a(i,j,k+1,n) - 2.0_rt*a(i,j,k,n) + a(i,j,k-1,n);
    }
  }
#endif

  return lapx + lapy + lapz;
}


#endif
