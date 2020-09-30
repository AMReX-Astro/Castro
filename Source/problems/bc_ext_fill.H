#ifndef BC_EXT_FILL_H
#define BC_EXT_FILL_H

#include <model_parser.H>

AMREX_INLINE AMREX_GPU_HOST_DEVICE
void
ext_denfill_c(const Box& bx, Array4<Real> const& adv, const int dcomp,
              Geometry const& geom, const BCRec& bcr,
              const Real time) {

    // Note: this function should not be needed, technically, but is
    // provided to filpatch because there are many times in the algorithm
    // when just the density is needed.  We try to rig up the filling so
    // that the same function is called here and in hypfill where all the
    // states are filled.

    auto domlo = geom.Domain().loVect3d();
    auto domhi = geom.Domain().hiVect3d();

    auto problo = geom.ProbLoArray();

    auto lo = bx.loVect();
    auto hi = bx.hiVect();

    auto adv_bx = Box(adv);
    auto adv_lo = adv_bx.loVect();
    auto adv_hi = adv_bx.hiVect();

    auto dx = geom.CellSizeArray();

    // XLO

    if (bcr.lo(0) == EXT_DIR && lo[0] < domlo[0]) {
        Box gbx(IntVect(D_DECL(adv_lo[0], lo[1], lo[2])),
                IntVect(D_DECL(domlo[0]-1, hi[1], hi[2])));

        amrex::ParallelFor(gbx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            Real x = problo[0] + dx[0]*(static_cast<Real>(i) + 0.5_rt);
            adv(i,j,k,dcomp) = interpolate(x, model::idens);
        });
    }


    // XHI

    if (bcr.hi(0) == EXT_DIR && hi[0] > domhi[0]) {
        Box gbx(IntVect(D_DECL(domhi[0]+1, lo[1], lo[2])),
                IntVect(D_DECL(adv_hi[0], hi[1], hi[2])));

        amrex::ParallelFor(gbx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            Real x = problo[0] + dx[0]*(static_cast<Real>(i) + 0.5_rt);
            adv(i,j,k,dcomp) = interpolate(x, model::idens);
        });
    }


#if AMREX_SPACEDIM >= 2

    // YLO

    if (bcr.lo(1) == EXT_DIR && lo[1] < domlo[1]) {
        Box gbx(IntVect(D_DECL(lo[0], adv_lo[1], lo[2])),
                IntVect(D_DECL(hi[0], domlo[1]-1, hi[2])));

        amrex::ParallelFor(gbx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            Real y = problo[1] + dx[1]*(static_cast<Real>(j) + 0.5_rt);
            adv(i,j,k,dcomp) = interpolate(y, model::idens);
        });
    }


    // YHI

    if (bcr.hi(1) == EXT_DIR && hi[1] > domhi[1]) {
        Box gbx(IntVect(D_DECL(lo[0], domhi[1]+1, lo[2])),
                IntVect(D_DECL(hi[0], adv_hi[1], hi[2])));

        amrex::ParallelFor(gbx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            Real y = problo[1] + dx[1]*(static_cast<Real>(j) + 0.5_rt);
            adv(i,j,k,dcomp) = interpolate(y, model::idens);
        });
    }
#endif

#if AMREX_SPACEDIM == 3

    // ZLO

    if (bcr.lo(2) == EXT_DIR && lo[2] < domlo[2]) {
        Box gbx(IntVect(D_DECL(lo[0], lo[1], adv_lo[2])),
                IntVect(D_DECL(hi[0], hi[1], domlo[2]-1)));

        amrex::ParallelFor(gbx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            Real z = problo[2] + dx[2]*(static_cast<Real>(k) + 0.5_rt);
            adv(i,j,k,dcomp) = interpolate(z, model::idens);
        });
    }

    // ZHI

    if (bcr.hi(2) == EXT_DIR && hi[2] > domhi[2]) {
        Box gbx(IntVect(D_DECL(lo[0], lo[1], domhi[2]+1)),
                IntVect(D_DECL(hi[0], hi[1], adv_hi[2])));

        amrex::ParallelFor(gbx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            Real z = problo[2] + dx[2]*(static_cast<Real>(k) + 0.5_rt);
            adv(i,j,k,dcomp) = interpolate(z, model::idens);
        });
    }
#endif
#endif

