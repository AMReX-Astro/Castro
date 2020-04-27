

// Note: pretty much all of these routines below assume that dx(1) = dx(2) = dx(3)

Real compute_laplacian(int i, int j, int k, int n,
                       Array4<Real const> const a,
                       GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi) {

  Real lapx = 0.0;
  Real lapy = 0.0;
  Real lapz = 0.0;

  // we use 2nd-order accurate one-sided stencils at the physical
  // boundaries note: this differs from the suggestion in MC2011 --
  // they just suggest using the Laplacian from +1 off the interior.
  // I like the one-sided better.

  if (i == domlo[0] && physbc_lo[0] != Interior) {
    lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i+1,j,k,n) + 4.0_rt*a(i+2,j,k,n) - a(i+3,j,k,n);

  } else if (i == domhi[0] && physbc_hi[0] != Interior) {
    lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i-1,j,k,n) + 4.0_rt*a(i-2,j,k,n) - a(i-3,j,k,n);

  } else {
    lapx = a(i+1,j,k,n) - TWO*a(i,j,k,n) + a(i-1,j,k,n);
  }

#if AMREX_SPACEDIM >= 2
  if (j == domlo[1] && physbc_lo[1] != Interior) {
    lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j+1,k,n) + 4.0_rt*a(i,j+2,k,n) - a(i,j+3,k,n);

  } else if (j == domhi[1] && physbc_hi[1] != Interior) {
    lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j-1,k,n) + 4.0_rt*a(i,j-2,k,n) - a(i,j-3,k,n);

  } else {
    lapy = a(i,j+1,k,n) - TWO*a(i,j,k,n) + a(i,j-1,k,n);
  }
#endif

#if AMREX_SPACEDIM == 3
  if (k == domlo[2] && physbc_lo[2] != Interior) {
    lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k+1,n) + 4.0_rt*a(i,j,k+2,n) - a(i,j,k+3,n);

  } else if (k == domhi[2] && physbc_hi[2] != Interior) {
    lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k-1,n) + 4.0_rt*a(i,j,k-2,n) - a(i,j,k-3,n);

  } else {
    lapz = a(i,j,k+1,n) - TWO*a(i,j,k,n) + a(i,j,k-1,n);
  }
#endif

  return lapx + lapy + lapz;
}


Real trans_laplacian(int i, int j, int k, int n,
                     int idir,
                     Array4<Real const> const a,
                     GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi) {

  Real lapx = 0.0;
  Real lapy = 0.0;
  Real lapz = 0.0;

  // we use 2nd-order accurate one-sided stencils at the physical boundaries
  if (idir != 0) {

    if (i == domlo[0] && physbc_lo[0] != Interior) {
      lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i+1,j,k,n) + 4.0_rt*a(i+2,j,k,n) - a(i+3,j,k,n);

    } else if (i == domhi[0] && physbc_hi[0] != Interior) {
      lapx = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i-1,j,k,n) + 4.0_rt*a(i-2,j,k,n) - a(i-3,j,k,n);

    } else {
      lapx = a(i+1,j,k,n) - TWO*a(i,j,k,n) + a(i-1,j,k,n);
    }
  }

  if (idir != 1) {

    if (j == domlo[1] && physbc_lo[1] != Interior) {
      lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j+1,k,n) + 4.0_rt*a(i,j+2,k,n) - a(i,j+3,k,n);

    } else if (j == domhi[1] && physbc_hi[1] != Interior) {
      lapy = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j-1,k,n) + 4.0_rt*a(i,j-2,k,n) - a(i,j-3,k,n);

    } else {
      lapy = a(i,j+1,k,n) - TWO*a(i,j,k,n) + a(i,j-1,k,n);
    }
  }

#if AMREX_SPACEDIM == 3
  if (idir != 2) {

    if (k == domlo[2] && physbc_lo[2] != Interior) {
      lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k+1,n) + 4.0_rt*a(i,j,k+2,n) - a(i,j,k+3,n);

    } else if (k == domhi[2] && physbc_hi[2] != Interior) {
      lapz = 2.0_rt*a(i,j,k,n) - 5.0_rt*a(i,j,k-1,n) + 4.0_rt*a(i,j,k-2,n) - a(i,j,k-3,n);

    } else {
      lapz = a(i,j,k+1,n) - TWO*a(i,j,k,n) + a(i,j,k-1,n);
    }
  }
#endif

  return lapx + lapy + lapz;
}


void make_cell_center(const Box& bx,
                      Array4<Real const> const& U,
                      Array4<Real> const& U_cc,
                      GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi) {

  // Take a cell-average state U and a convert it to a cell-center
  // state U_cc via U_cc = U - 1/24 L U

  U_lo = U.lbound();
  U_hi = U.ubound();

  AMREX_ASSERT(U_lo[0] <= lo[0]-1 && U_hi[0] >= hi[0]+1 &&
               (AMREX_SPACEDIM >= 2 && (U_lo[1] <= lo[1]-1 && U_hi[1] >= hi[1]+1)) &&
               (AMREX_SPACEDIM == 3 && (U_lo[2] <= lo[2]-1 && U_hi[2] >= hi[2]+1)));

  amrex::ParallelFor(bx, U.nComp(),
  [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real lap = compute_laplacian(i, j, k, n, U, domlo, domhi);

    U_cc(i,j,k,n) = U(i,j,k,n) - (1.0_rt/24.0_rt) * lap;

  });
}

void make_cell_center_in_place(const Box& bx,
                               Array4<Real> const& U,
                               Array4<Real> const& tmp,
                               GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi) {

  // Take a cell-average state U and make it cell-centered in place
  // via U <- U - 1/24 L U.  Note that this operation is not tile
  // safe.

  // here, tmp is a temporary memory space we use to store one-component's Laplacian

  for (int n = 0; n < U.nComp(); n++) {

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      tmp(i,j,k) = compute_laplacian(i, j, k, n, U, domlo, domhi);
    });

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      U(i,j,k,n) = U(i,j,k,n) - (1.0_rt/24.0_rt) * tmp(i,j,k);
    });
  }
}


void compute_lap_term(const Box& bx,
                      Array4<Real const> const& Y,
                      Array4<Real> const& lap, const int ncomp,
                      GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi) {

  // Computes the h**2/24 L U term that is used in correcting
  // cell-center to averages (and back)

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    lap(i,j,k) = (1.0_rt/24.0_rt) *
      compute_laplacian(i, j, k, ncomp, U, domlo, domhi);
  });

}


void make_fourth_average(const Box& bx,
                         Array4<Real> const& q,
                         Array4<Real const> const& q_bar,
                         GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi) {

  // Take the cell-center state q and another state q_bar (e.g.,
  // constructed from the cell-average U) and replace the cell-center
  // q with a 4th-order accurate cell-average, q <- q + 1/24 L q_bar

  amrex::ParallelFor(bx, q.nComp(),
  [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real lap = compute_laplacian(i, j, k, n, q_bar, domlo, domhi);

    q(i,j,k,n) += (1.0_rt/24.0_rt) * lap;
  });
}


void make_fourth_in_place(const Box& bx,
                          Array4<Real> const& q,
                          Array4<Real> const& tmp,
                          GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi) {

  // Take the cell-center q and makes it a cell-average q, in place
  // (e.g. q is overwritten by its average), q <- q + 1/24 L q.
  // Note: this routine is not tile safe.

  for (int n = 0; n < q.nComp(); n++) {
    make_fourth_in_place_n(bx, q, n, tmp, domlo, domhi);
  }
}


void make_fourth_in_place_n(const Box& bx,
                            Array4<Real> const& q, const int ncomp,
                            Array4<Real> const& tmp,
                            GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi) {

  // Take the cell-center q and makes it a cell-average q, in place
  // (e.g. q is overwritten by its average), q <- q + 1/24 L q.
  // Note: this routine is not tile safe.
  //
  // This version operates on a single component.  Here ncomp is the
  // component to update

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    tmp(i,j,k) = compute_laplacian(i, j, k, ncomp, q, domlo, domhi);
  });

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    q(i,j,k,ncomp) += (1.0_rt/24.0_rt) * lap(i,j,k);
  });

}
