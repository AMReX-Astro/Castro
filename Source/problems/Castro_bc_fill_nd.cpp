
#include <AMReX_BLFort.H>
#include <Castro.H>
#include <Castro_bc_fill_nd.H>
#include <Castro_bc_fill_nd_F.H>
#include <Castro_bc_ext_fill_nd_F.H>
#include <Castro_generic_fill.H>
#include <bc_ext_fill.H>

using namespace amrex;

void ca_statefill(Box const& bx, FArrayBox& data,
                  const int dcomp, const int numcomp,
                  Geometry const& geom, const Real time,
                  const Vector<BCRec>& bcr, const int bcomp,
                  const int scomp)
{

    // Here dcomp is the component in the destination array that we
    // are filling and bcr is a vector of length ncomp which are the
    // BC values corresponding to components dcomp to dcomp + ncomp -
    // 1

    // Make a copy of the raw BCRec data in the format
    // our BC routines can handle (a contiguous array
    // of integers).

    Vector<int> bcrs(2 * AMREX_SPACEDIM * numcomp);

    for (int n = 0; n < numcomp; ++n) {
        for (int k = 0; k < 2 * AMREX_SPACEDIM; ++k) {
            bcrs[2 * AMREX_SPACEDIM * n + k] = bcr[n].vect()[k];
        }
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bcrs.data(), numcomp);
#else
    const int* bc_f = bcrs.data();
#endif

    if (Gpu::inLaunchRegion()) {
        GpuBndryFuncFab<CastroGenericFill> gpu_bndry_func(castro_generic_fill_func);
        gpu_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
    }
    else {
        CpuBndryFuncFab cpu_bndry_func(nullptr);
        cpu_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
    }

    // This routine either comes in with one component or all NUM_STATE.

    if (numcomp == 1) {

#pragma gpu box(bx)
        denfill(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                BL_TO_FORTRAN_N_ANYD(data, dcomp),
                AMREX_INT_ANYD(geom.Domain().loVect()), AMREX_INT_ANYD(geom.Domain().hiVect()),
                AMREX_REAL_ANYD(geom.CellSize()), AMREX_REAL_ANYD(geom.ProbLo()), time, bc_f);

    }
    else {

        AMREX_ALWAYS_ASSERT(numcomp == NUM_STATE);

#pragma gpu box(bx)
        hypfill(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                BL_TO_FORTRAN_ANYD(data),
                AMREX_INT_ANYD(geom.Domain().loVect()), AMREX_INT_ANYD(geom.Domain().hiVect()),
                AMREX_REAL_ANYD(geom.CellSize()), AMREX_REAL_ANYD(geom.ProbLo()), time, bc_f);

    }

    // we just did the standard BC fills (reflect, outflow, ...)  now
    // we consider the external ones (HSE).  Note, if we are at a
    // corner where two (or three) faces want to do HSE, we may run
    // into a situation that the data is not valid in the corner where
    // we start the integration.  We'll abort, for now, if we run into
    // this case.
    //
    // The future fix is to first call ext_fill on the ghost cells
    // that are not corners and then call it a second time on just the
    // corners.

#if AMREX_SPACEDIM == 2
    if ((bcr[URHO].lo(0) == EXT_DIR && bcr[URHO].lo(1) == EXT_DIR) ||
        (bcr[URHO].lo(0) == EXT_DIR && bcr[URHO].hi(1) == EXT_DIR) ||
        (bcr[URHO].hi(0) == EXT_DIR && bcr[URHO].lo(1) == EXT_DIR) ||
        (bcr[URHO].hi(0) == EXT_DIR && bcr[URHO].hi(1) == EXT_DIR)) {
        amrex::Error("Error: external boundaries meeting at a corner not supported");
    }
#endif

#if AMREX_SPACEDIM == 3
    if ((bcr[URHO].lo(0) == EXT_DIR &&           // xl, yl, zl corner
         (bcr[URHO].lo(1) == EXT_DIR || bcr[URHO].lo(2) == EXT_DIR)) ||
        (bcr[URHO].lo(1) == EXT_DIR && bcr[URHO].lo(2) == EXT_DIR) ||
        (bcr[URHO].lo(0) == EXT_DIR &&           // xl, yr, zl corner
         (bcr[URHO].hi(1) == EXT_DIR || bcr[URHO].lo(2) == EXT_DIR)) ||
        (bcr[URHO].hi(1) == EXT_DIR && bcr[URHO].lo(2) == EXT_DIR) ||
        (bcr[URHO].lo(0) == EXT_DIR &&           // xl, yl, zr corner
         (bcr[URHO].lo(1) == EXT_DIR || bcr[URHO].hi(2) == EXT_DIR)) ||
        (bcr[URHO].lo(1) == EXT_DIR && bcr[URHO].hi(2) == EXT_DIR) ||
        (bcr[URHO].lo(0) == EXT_DIR &&           // xl, yr, zr corner
         (bcr[URHO].hi(1) == EXT_DIR || bcr[URHO].hi(2) == EXT_DIR)) ||
        (bcr[URHO].hi(1) == EXT_DIR && bcr[URHO].hi(2) == EXT_DIR) ||
        (bcr[URHO].hi(0) == EXT_DIR &&           // xr, yl, zl corner
         (bcr[URHO].lo(1) == EXT_DIR || bcr[URHO].lo(2) == EXT_DIR)) ||
        (bcr[URHO].lo(1) == EXT_DIR && bcr[URHO].lo(2) == EXT_DIR) ||
        (bcr[URHO].hi(0) == EXT_DIR &&           // xr, yr, zl corner
         (bcr[URHO].hi(1) == EXT_DIR || bcr[URHO].lo(2) == EXT_DIR)) ||
        (bcr[URHO].hi(1) == EXT_DIR && bcr[URHO].lo(2) == EXT_DIR) ||
        (bcr[URHO].hi(0) == EXT_DIR &&           // xr, yl, zr corner
         (bcr[URHO].lo(1) == EXT_DIR || bcr[URHO].hi(2) == EXT_DIR)) ||
        (bcr[URHO].lo(1) == EXT_DIR && bcr[URHO].hi(2) == EXT_DIR) ||
        (bcr[URHO].hi(0) == EXT_DIR &&           // xr, yr, zr corner
         (bcr[URHO].hi(1) == EXT_DIR || bcr[URHO].hi(2) == EXT_DIR)) ||
        (bcr[URHO].hi(1) == EXT_DIR && bcr[URHO].hi(2) == EXT_DIR)) {
        amrex::Error("Error: external boundaries meeting at a corner not supported");
    }
#endif

    if (numcomp == 1) {

#ifndef CXX_MODEL_PARSER
#pragma gpu box(bx)
        ext_denfill(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                    BL_TO_FORTRAN_N_ANYD(data, dcomp),
                    AMREX_INT_ANYD(geom.Domain().loVect()), AMREX_INT_ANYD(geom.Domain().hiVect()),
                    AMREX_REAL_ANYD(geom.CellSize()), AMREX_REAL_ANYD(geom.ProbLo()), time, bc_f);
#else
        ext_denfill_c(bx, data.array(dcomp), geom, bcr[0], time);
#endif

    }
    else {

        AMREX_ALWAYS_ASSERT(numcomp == NUM_STATE);

#pragma gpu box(bx)
        ext_fill(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                 BL_TO_FORTRAN_ANYD(data),
                 AMREX_INT_ANYD(geom.Domain().loVect()), AMREX_INT_ANYD(geom.Domain().hiVect()),
                 AMREX_REAL_ANYD(geom.CellSize()), AMREX_REAL_ANYD(geom.ProbLo()), time, bc_f);

    }

#ifdef AMREX_USE_CUDA
    clean_bc(bc_f);
#endif

  }



#ifdef MHD
  void ca_face_fillx(Real* var, const int* var_lo, const int* var_hi,
                     const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                     const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = var_lo[i];
      hi[i] = var_hi[i];
    }

    const int* bc_f = bc;


    face_fillx(AMREX_ARLIM_ANYD(lo), AMREX_ARLIM_ANYD(hi),
               var, AMREX_ARLIM_ANYD(var_lo), AMREX_ARLIM_ANYD(var_hi),
               AMREX_ARLIM_ANYD(domlo), AMREX_ARLIM_ANYD(domhi),
               AMREX_ZFILL(dx), AMREX_ZFILL(xlo), *time, bc_f);

  }

  void ca_face_filly(Real* var, const int* var_lo, const int* var_hi,
                     const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                     const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = var_lo[i];
      hi[i] = var_hi[i];
    }

    const int* bc_f = bc;


    face_filly(AMREX_ARLIM_ANYD(lo), AMREX_ARLIM_ANYD(hi),
               var, AMREX_ARLIM_ANYD(var_lo), AMREX_ARLIM_ANYD(var_hi),
               AMREX_ARLIM_ANYD(domlo), AMREX_ARLIM_ANYD(domhi),
               AMREX_ZFILL(dx), AMREX_ZFILL(xlo), *time, bc_f);

  }

  void ca_face_fillz(Real* var, const int* var_lo, const int* var_hi,
                     const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                     const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = var_lo[i];
      hi[i] = var_hi[i];
    }

    const int* bc_f = bc;


    face_fillz(AMREX_ARLIM_ANYD(lo), AMREX_ARLIM_ANYD(hi),
               var, AMREX_ARLIM_ANYD(var_lo), AMREX_ARLIM_ANYD(var_hi),
               AMREX_ARLIM_ANYD(domlo), AMREX_ARLIM_ANYD(domhi),
               AMREX_ZFILL(dx), AMREX_ZFILL(xlo), *time, bc_f);

  }
#endif  

