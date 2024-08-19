#include <Castro.H>
#include <fourth_center_average.H>

using namespace amrex;

// Note: pretty much all of these routines below assume that dx(1) = dx(2) = dx(3)


void
Castro::make_cell_center(const Box& bx,
                         Array4<Real const> const& U,
                         Array4<Real> const& U_cc,
                         GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi) {

  // Take a cell-average state U and a convert it to a cell-center
  // state U_cc via U_cc = U - 1/24 L U

  auto U_lo = lbound(U);
  auto U_hi = ubound(U);

  const auto *lo = bx.loVect();
  const auto *hi = bx.hiVect();

  AMREX_ASSERT(U_lo.x <= lo[0]-1 && U_hi.x >= hi[0]+1);
#if AMREX_SPACEDIM >= 2
  AMREX_ASSERT(U_lo.y <= lo[1]-1 && U_hi.y >= hi[1]+1);
#endif
#if AMREX_SPACEDIM == 3
  AMREX_ASSERT(U_lo.z <= lo[2]-1 && U_hi.z >= hi[2]+1);
#endif

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  GpuArray<bool, AMREX_SPACEDIM> lo_periodic;
  GpuArray<bool, AMREX_SPACEDIM> hi_periodic;
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    lo_periodic[idir] = lo_bc[idir] == amrex::PhysBCType::interior;
    hi_periodic[idir] = hi_bc[idir] == amrex::PhysBCType::interior;
  }

  amrex::ParallelFor(bx, U.nComp(),
  [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real lap = compute_laplacian(i, j, k, n, U,
                                 lo_periodic, hi_periodic, domlo, domhi);

    U_cc(i,j,k,n) = U(i,j,k,n) - (1.0_rt/24.0_rt) * lap;

  });
}

void
Castro::make_cell_center_in_place(const Box& bx,
                                  Array4<Real> const& U,
                                  Array4<Real> const& tmp,
                                  GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi) {

  // Take a cell-average state U and make it cell-centered in place
  // via U <- U - 1/24 L U.  Note that this operation is not tile
  // safe.

  // here, tmp is a temporary memory space we use to store one-component's Laplacian

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  GpuArray<bool, AMREX_SPACEDIM> lo_periodic;
  GpuArray<bool, AMREX_SPACEDIM> hi_periodic;
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    lo_periodic[idir] = lo_bc[idir] == amrex::PhysBCType::interior;
    hi_periodic[idir] = hi_bc[idir] == amrex::PhysBCType::interior;
  }

  for (int n = 0; n < U.nComp(); n++) {

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      tmp(i,j,k) = compute_laplacian(i, j, k, n, U,
                                     lo_periodic, hi_periodic, domlo, domhi);
    });

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      U(i,j,k,n) = U(i,j,k,n) - (1.0_rt/24.0_rt) * tmp(i,j,k);
    });
  }
}


void
Castro::compute_lap_term(const Box& bx,
                         Array4<Real const> const& U,
                         Array4<Real> const& lap, const int ncomp,
                         GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi) {

  // Computes the h**2/24 L U term that is used in correcting
  // cell-center to averages (and back)

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  GpuArray<bool, AMREX_SPACEDIM> lo_periodic;
  GpuArray<bool, AMREX_SPACEDIM> hi_periodic;
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    lo_periodic[idir] = lo_bc[idir] == amrex::PhysBCType::interior;
    hi_periodic[idir] = hi_bc[idir] == amrex::PhysBCType::interior;
  }

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    lap(i,j,k) = (1.0_rt/24.0_rt) *
      compute_laplacian(i, j, k, ncomp, U,
                        lo_periodic, hi_periodic, domlo, domhi);
  });

}


void
Castro::make_fourth_average(const Box& bx,
                            Array4<Real> const& q,
                            Array4<Real const> const& q_bar,
                            GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi) {

  // Take the cell-center state q and another state q_bar (e.g.,
  // constructed from the cell-average U) and replace the cell-center
  // q with a 4th-order accurate cell-average, q <- q + 1/24 L q_bar

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  GpuArray<bool, AMREX_SPACEDIM> lo_periodic;
  GpuArray<bool, AMREX_SPACEDIM> hi_periodic;
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    lo_periodic[idir] = lo_bc[idir] == amrex::PhysBCType::interior;
    hi_periodic[idir] = hi_bc[idir] == amrex::PhysBCType::interior;
  }

  amrex::ParallelFor(bx, q.nComp(),
  [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real lap = compute_laplacian(i, j, k, n, q_bar,
                                 lo_periodic, hi_periodic, domlo, domhi);

    q(i,j,k,n) += (1.0_rt/24.0_rt) * lap;
  });
}


void
Castro::make_fourth_in_place(const Box& bx,
                             Array4<Real> const& q,
                             Array4<Real> const& tmp,
                             GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi) {

  // Take the cell-center q and makes it a cell-average q, in place
  // (e.g. q is overwritten by its average), q <- q + 1/24 L q.
  // Note: this routine is not tile safe.

  for (int n = 0; n < q.nComp(); n++) {
    make_fourth_in_place_n(bx, q, n, tmp, domlo, domhi);
  }
}


void
Castro::make_fourth_in_place_n(const Box& bx,
                               Array4<Real> const& q, const int ncomp,
                               Array4<Real> const& tmp,
                               GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi) {

  // Take the cell-center q and makes it a cell-average q, in place
  // (e.g. q is overwritten by its average), q <- q + 1/24 L q.
  // Note: this routine is not tile safe.
  //
  // This version operates on a single component.  Here ncomp is the
  // component to update

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  GpuArray<bool, AMREX_SPACEDIM> lo_periodic;
  GpuArray<bool, AMREX_SPACEDIM> hi_periodic;
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    lo_periodic[idir] = lo_bc[idir] == amrex::PhysBCType::interior;
    hi_periodic[idir] = hi_bc[idir] == amrex::PhysBCType::interior;
  }

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    tmp(i,j,k) = compute_laplacian(i, j, k, ncomp, q,
                                   lo_periodic, hi_periodic, domlo, domhi);
  });

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    q(i,j,k,ncomp) += (1.0_rt/24.0_rt) * tmp(i,j,k);
  });

}
