#include "Castro.H"
#include <Castro_F.H>
#include "prob_parameters.H"

using namespace amrex;

void
Castro::fill_ext_source (const Real time, const Real dt, 
                         const MultiFab& state_old, 
                         const MultiFab& state_new, MultiFab& ext_src)
{
    // Compute the external sources for all the conservative equations.
    //
    // This is called twice in the evolution:
    //
    // First, for the predictor, it is called with (old, old) states.
    //
    // This is also used in the first pass of the conservative update
    // (adding dt * S there).
    //
    // Next we correct the source terms in the conservative update to
    // time-center them.  Here we call ext_src(old, new), and then
    // in time_center_source_terms we subtract off 1/2 of the first S
    // and add 1/2 of the new S.
    //
    // Therefore, to get a properly time-centered source, generally
    // speaking, you always want to use the "new" state here.  That
    // will be the time n state in the first call and the n+1 in the
    // second call.

    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    auto domain_lo = geom.Domain().loVect3d();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    const auto H_0 = heating_peak;
    const auto W_0 = heating_sigma;
    const auto r_0 = heating_rad;
    const auto t_stop = heating_time;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(ext_src, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        
        const auto lo = bx.loVect3d();

        if (lo[0] == 0 && lo[1] == 0 && lo[2] == 0) {
            AllPrint() << "TIME vs TSTOP " << time << " " << t_stop;
        }

        Array4<Real const> const sold = state_old.array(mfi);
        Array4<Real> const src = ext_src.array(mfi);

        if (prob_type == 1) {

            // For heating at the center

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                auto x = domain_lo[0] + (Real(i) + 0.5_rt)*dx[0] - center[0];
                auto y = domain_lo[1] + (Real(j) + 0.5_rt)*dx[1] - center[1];
                auto z = domain_lo[2] + (Real(k) + 0.5_rt)*dx[2] - center[2];

                auto dist = std::sqrt(x*x + y*y + z*z);

                auto Hext = H_0 * std::exp(-((dist - r_0)*(dist - r_0))/(W_0*W_0));

                src(i,j,k,UEINT) = sold(i,j,k,URHO) * Hext;
                src(i,j,k,UEDEN) = sold(i,j,k,URHO) * Hext;
            });

        } else if (prob_type == 3) {

            // sub-chandra heating -- modulate by He

            const auto ihe4 = network_spec_index("helium-4");

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                auto x = domain_lo[0] + (Real(i) + 0.5_rt)*dx[0] - center[0];
                auto y = domain_lo[1] + (Real(j) + 0.5_rt)*dx[1] - center[1];
                auto z = domain_lo[2] + (Real(k) + 0.5_rt)*dx[2] - center[2];

                auto dist = std::sqrt(x*x + y*y + z*z);

                auto Hext = H_0 * std::exp(-((dist - r_0)*(dist - r_0))/(W_0*W_0)) * 
                    sold(i,j,k,UFS-1+ihe4) / sold(i,j,k,URHO);

                src(i,j,k,UEINT) = sold(i,j,k,URHO) * Hext;
                src(i,j,k,UEDEN) = sold(i,j,k,URHO) * Hext;
            });
        }
    }
}
