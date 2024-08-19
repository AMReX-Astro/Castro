#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

using namespace amrex;

void
apply_metric(const Box& bx,
             Array4<Real> const rhs, const Box& rbx,
             Array4<Real> const ecx, const Box& xbx,
#if AMREX_SPACEDIM >= 2
             Array4<Real> const ecy, const Box& ybx,
#endif
             GpuArray<Real, AMREX_SPACEDIM> dx,
             const int coord_type)
{
    // r-z
    if (coord_type == 1) {

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            IntVect idx(AMREX_D_DECL(i, j, k));

            // at centers
            if (rbx.contains(idx)) {
                Real r = (static_cast<Real>(i) + 0.5_rt) * dx[0];
                rhs(i,j,k) *= r;
            }

            // On x-edges
            if (xbx.contains(idx)) {
                Real r = static_cast<Real>(i) * dx[0];
                ecx(i,j,k) *= r;
            }

#if AMREX_SPACEDIM >= 2
            // On y-edges
            if (ybx.contains(idx)) {
                Real r = (static_cast<Real>(i) + 0.5_rt) * dx[0];
                ecy(i,j,k) *= r;
            }
#endif

        });

    } else {
        amrex::Print() << "Bogus coord_type in apply_metric " << coord_type << std::endl;
        amrex::Error("Error:: MGutils.cpp :: ca_apply_metric");
    }
}


void
do_weight_cc(const Box& bx,
             Array4<Real> const cc,
             GpuArray<Real, AMREX_SPACEDIM> dx,
             const int coord_type)
{
    // r-z
    if (coord_type == 1) {

        // At centers
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real r = (static_cast<Real>(i) + 0.5_rt) * dx[0];
            cc(i,j,k) *= r;
        });

    } else {
      amrex::Print() << "Bogus coord_type in weight_cc " << coord_type << std::endl;
      amrex::Error("Error:: MGutils.cpp :: ca_weight_cc");
    }
}


void
do_unweight_cc(const Box& bx,
               Array4<Real> const cc,
               GpuArray<Real, AMREX_SPACEDIM> dx,
               const int coord_type)
{
    // r-z
    if (coord_type == 1) {

        // At centers
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real r = (static_cast<Real>(i) + 0.5_rt) * dx[0];
            cc(i,j,k) /= r;
        });

    } else {
        amrex::Print() << "Bogus coord_type in unweight_cc " << coord_type << std::endl;
        amrex::Error("Error:: MGutils.cpp :: ca_unweight_cc");
    }
}


void
do_unweight_edges(const Box& bx,
                  Array4<Real> const ec,
                  const int idir,
                  GpuArray<Real, AMREX_SPACEDIM> dx,
                  const int coord_type)
{
    // r-z
    if (coord_type == 1) {

        if (idir == 0) {

            // On x-edges
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i != 0) {
                    Real r = static_cast<Real>(i) * dx[0];
                    ec(i,j,k) /= r;
                }
            });

        } else {

            // On y-edges
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real r = (static_cast<Real>(i) + 0.5_rt) * dx[0];
                ec(i,j,k) /= r;
            });
        }

    } else {
        amrex::Print() << "Bogus coord_type in unweight_edges " << coord_type << std::endl;
        amrex::Error("Error:: MGutils.cpp :: ca_unweight_edges");
    }
}

