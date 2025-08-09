#include <Castro.H>
#include <global.H>
#include <timestep.H>

#ifdef DIFFUSION
#include <conductivity.H>
#endif

#ifdef MHD
#include <mhd_util.H>
#endif

#ifdef ROTATION
#include <Rotation.H>
#endif

#ifdef REACTIONS
#include <actual_network.H>
#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>
#else
#include <actual_rhs.H>
#endif
#endif

#ifdef RADIATION
#include <Radiation.H>
#endif

using namespace amrex;

#ifdef RADIATION
Real
timestep::estdt_rad (const MultiFab& stateMF, const MultiFab& radMF, const GeometryData& geomdata)
{
    // Compute radiation + hydro limited timestep.

    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(stateMF, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbox = mfi.tilebox();
        const Box& vbox = mfi.validbox();

        FArrayBox gPr;
        gPr.resize(tbox);
        global::the_radiation_ptr->estimate_gamrPr(stateMF[mfi], radMF[mfi], gPr, geomdata.CellSize(), vbox);

        auto u = stateMF[mfi].array();
        auto gPr_arr = gPr.array();

        // Note that we synchronize in the call to estimate_gamrPr. Ideally
        // we would merge these into one loop later (probably by making that
        // call occur on a zone-by-zone basis) so that we can simplify.

        reduce_op.eval(tbox, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real rhoInv = 1.0_rt / u(i,j,k,URHO);

            const auto* dx = geomdata.CellSize();
#if AMREX_SPACEDIM >= 2
            const auto* problo = geomdata.ProbLo();
            const auto coord = geomdata.Coord();
#endif

            eos_t eos_state;
            eos_state.rho = u(i,j,k,URHO);
            eos_state.T   = u(i,j,k,UTEMP);
            eos_state.e   = u(i,j,k,UEINT) * rhoInv;
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = u(i,j,k,UFS+n) * rhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; ++n) {
                eos_state.aux[n] = u(i,j,k,UFX+n) * rhoInv;
            }
#endif

            eos(eos_input_re, eos_state);

            Real c = eos_state.cs;
            c = std::sqrt(c * c + gPr_arr(i,j,k) * rhoInv);

            Real ux = u(i,j,k,UMX) * rhoInv;
            Real uy = u(i,j,k,UMY) * rhoInv;
            Real uz = u(i,j,k,UMZ) * rhoInv;

            Real dt1 = dx[0] / (c + std::abs(ux));
#if AMREX_SPACEDIM >= 2
            Real dt2 = dx[1] / (c + std::abs(uy));
            if (coord == 2) {
                dt2 *= problo[0] + 0.5_rt * dx[0];
            }
#else
            Real dt2 = std::numeric_limits<Real>::max();
#endif
#if AMREX_SPACEDIM == 3
            Real dt3 = dx[2] / (c + std::abs(uz));
#else
            Real dt3 = std::numeric_limits<Real>::max();
#endif

            Real dt_min = amrex::min(dt1, dt2, dt3);

            return {dt_min};
        });

        Gpu::synchronize();
    }

    ReduceTuple hv = reduce_data.value();
    Real estdt = amrex::get<0>(hv);
    estdt = std::min(estdt, max_dt / castro::cfl);

    return estdt;
}
#endif
