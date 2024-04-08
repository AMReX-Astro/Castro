
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>

#include <HypreABec.H>
#include <HABEC.H>
#include <rad_util.H>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <_hypre_struct_mv.h>

using namespace amrex;

static int ispow2(int i)
{
  return (i == 1) ? 1 : (((i <= 0) || (i & 1)) ? 0 : ispow2(i / 2));
}

Real HypreABec::flux_factor = 1.0;

#if (AMREX_SPACEDIM == 1)
static int vl[2] = { 0, 0 };
static int vh[2] = { 0, 0 };
#endif

static int* loV(const Box& b) {
#if (AMREX_SPACEDIM == 1)
  vl[0] = b.smallEnd(0);
  return vl;
#else
  return (int*) b.loVect();
#endif
}

static int* hiV(const Box& b) {
#if (AMREX_SPACEDIM == 1)
  vh[0] = b.bigEnd(0);
  return vh;
#else
  return (int*) b.hiVect();
#endif
}

HypreABec::HypreABec(const BoxArray& grids,
                     const DistributionMapping& dmap,
                     const Geometry& _geom,
                     int _solver_flag)
  : geom(_geom), solver_flag(_solver_flag)
{
  ParmParse pp("habec");

  pfmg_relax_type = 1; pp.query("pfmg_relax_type", pfmg_relax_type);
  verbose = 0; pp.query("v", verbose); pp.query("verbose", verbose);
  verbose_threshold = 0; pp.query("verbose_threshold", verbose_threshold);

  static int first = 1;
  if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
    first = 0;
    if (solver_flag == 1 || solver_flag == 3 || solver_flag == 5) {
      std::cout << "habec.pfmg_relax_type           = " << pfmg_relax_type << std::endl;
    }
    std::cout << "habec.verbose                   = " << verbose << std::endl;
    std::cout << "habec.verbose_threshold         = " << verbose_threshold << std::endl;
  }
  bho = 0; // higher order boundaries don't work with symmetric matrices

  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    dx[i] = geom.CellSize(i);
  }

#if (AMREX_SPACEDIM == 1)

  // Hypre doesn't support 1D directly, so we use 2D Hypre with
  // the second dimension collapsed.
  // (SMG reduces to cyclic reduction in this case, so it's an exact solve.)
  // (PFMG will not work.)

  HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &hgrid);

  if (geom.isAnyPeriodic()) {
    BL_ASSERT(geom.isPeriodic(0));
    BL_ASSERT(geom.Domain().smallEnd(0) == 0);

    int is_periodic[2];
    is_periodic[0] = geom.period(0);
    is_periodic[1] = 0;
    BL_ASSERT(ispow2(is_periodic[0]));

    HYPRE_StructGridSetPeriodic(hgrid, is_periodic);
  }

#else

  HYPRE_StructGridCreate(MPI_COMM_WORLD, AMREX_SPACEDIM, &hgrid);

  if (geom.isAnyPeriodic()) {
    int is_periodic[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      is_periodic[i] = 0;
      if (geom.isPeriodic(i)) {
        is_periodic[i] = geom.period(i);
        BL_ASSERT(ispow2(is_periodic[i]));
        BL_ASSERT(geom.Domain().smallEnd(i) == 0);
      }
    }
    HYPRE_StructGridSetPeriodic(hgrid, is_periodic);
  }
#endif

  if (ParallelDescriptor::NProcs() != 1) {
    // parallel section:
    for (int i = 0; i < grids.size(); i++) {
      if (dmap[i] == ParallelDescriptor::MyProc()) {
        HYPRE_StructGridSetExtents(hgrid, loV(grids[i]), hiV(grids[i]));
      }
    }
  }
  else {
    for (int i = 0; i < grids.size(); i++) {
      HYPRE_StructGridSetExtents(hgrid, loV(grids[i]), hiV(grids[i]));
    }
  }

  HYPRE_StructGridAssemble(hgrid);

#if (AMREX_SPACEDIM == 1)
  // if we were really 1D:
/*
  int offsets[2][1] = {{-1},
                       { 0}};
*/
  // fake 1D as a 2D problem:
  int offsets[2][2] = {{-1,  0},
                       { 0,  0}};
#elif (AMREX_SPACEDIM == 2)
  int offsets[3][2] = {{-1,  0},
                       { 0, -1},
                       { 0,  0}};
#elif (AMREX_SPACEDIM == 3)
  int offsets[4][3] = {{-1,  0,  0},
                       { 0, -1,  0},
                       { 0,  0, -1},
                       { 0,  0,  0}};
#endif

#if   (AMREX_SPACEDIM == 1)
  int A_num_ghost[6] = { 1, 1, 0, 0, 0, 0 };
#elif (AMREX_SPACEDIM == 2)
  //int A_num_ghost[4] = { 1, 1, 1, 1 };
  int A_num_ghost[6] = { 1, 1, 1, 1, 0, 0 };
#elif (AMREX_SPACEDIM == 3)
  int A_num_ghost[6] = { 1, 1, 1, 1, 1, 1 };
#endif

  HYPRE_StructStencil stencil;

#if (AMREX_SPACEDIM == 1)
  HYPRE_StructStencilCreate(2, 2, &stencil);
#else
  HYPRE_StructStencilCreate(AMREX_SPACEDIM, AMREX_SPACEDIM + 1, &stencil);
#endif

  for (int i = 0; i < AMREX_SPACEDIM + 1; i++) {
    HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
  }

  HYPRE_StructMatrixCreate(MPI_COMM_WORLD, hgrid, stencil, &A);
  HYPRE_StructMatrixSetSymmetric(A, 1);
  HYPRE_StructMatrixSetNumGhost(A, A_num_ghost);
  HYPRE_StructMatrixInitialize(A);

  HYPRE_StructMatrixCreate(MPI_COMM_WORLD, hgrid, stencil, &A0);
  HYPRE_StructMatrixSetSymmetric(A0, 1);
  HYPRE_StructMatrixSetNumGhost(A0, A_num_ghost);
  HYPRE_StructMatrixInitialize(A0);

  //HYPRE_StructVectorCreate(MPI_COMM_WORLD, hgrid, stencil, &b);
  //HYPRE_StructVectorCreate(MPI_COMM_WORLD, hgrid, stencil, &x);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, hgrid, &b);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, hgrid, &x);

  HYPRE_StructStencilDestroy(stencil); // no longer needed

  HYPRE_StructVectorInitialize(b);
  HYPRE_StructVectorInitialize(x);

  Gpu::synchronize();

  int ncomp=1;
  int ngrow=0;
  acoefs.reset(new MultiFab(grids, dmap, ncomp, ngrow));
  acoefs->setVal(0.0);

  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    BoxArray edge_boxes(grids);
    edge_boxes.surroundingNodes(i);
    bcoefs[i].reset(new MultiFab(edge_boxes, dmap, ncomp, ngrow));
  }
}

HypreABec::~HypreABec()
{
  HYPRE_StructVectorDestroy(b);
  HYPRE_StructVectorDestroy(x);

  HYPRE_StructMatrixDestroy(A);
  HYPRE_StructMatrixDestroy(A0);

  HYPRE_StructGridDestroy(hgrid);
}

void HypreABec::setScalars(Real Alpha, Real Beta)
{
  alpha = Alpha;
  beta  = Beta;
}

void HypreABec::aCoefficients(const MultiFab &a)
{
  BL_ASSERT( a.ok() );
  BL_ASSERT( a.boxArray() == acoefs->boxArray() );
  MultiFab::Copy(*acoefs, a, 0, 0, 1, 0);
}

void HypreABec::bCoefficients(const MultiFab &b, int dir)
{
  BL_ASSERT( b.ok() );
  BL_ASSERT( b.boxArray() == bcoefs[dir]->boxArray() );
  MultiFab::Copy(*bcoefs[dir], b, 0, 0, 1, 0);
}

void HypreABec::SPalpha(const MultiFab& a)
{
  BL_ASSERT( a.ok() );
  if (SPa == 0) {
    const BoxArray& grids = a.boxArray();
    const DistributionMapping& dmap = a.DistributionMap();
    SPa.reset(new MultiFab(grids,dmap,1,0));
  }
  MultiFab::Copy(*SPa, a, 0, 0, 1, 0);
}

void HypreABec::boundaryFlux(MultiFab* Flux, MultiFab& Soln, int icomp,
                             BC_Mode inhom)
{
    BL_PROFILE("HypreABec::boundaryFlux");

    const BoxArray &grids = Soln.boxArray();

    const NGBndry& bd = getBndry();
    const Box& domain = bd.getDomain();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Real foo=1.e200;

        for (MFIter si(Soln); si.isValid(); ++si) {
            int i = si.index();
            const Box &reg = grids[i];
            for (OrientationIter oitr; oitr; oitr++) {
                int cdir(oitr());
                int idim = oitr().coordDir();
                const RadBoundCond &bct = bd.bndryConds(oitr())[i];
                const Real      &bcl = bd.bndryLocs(oitr())[i];
                const FArrayBox       &fs  = bd.bndryValues(oitr())[si];
                const Mask      &msk = bd.bndryMasks(oitr(),i);

                if (reg[oitr()] == domain[oitr()]) {
                    int bctype = bct;
                    Array4<int const> tf_arr;
                    if (bd.mixedBndry(oitr())) {
                        const BaseFab<int> &tf = *(bd.bndryTypes(oitr())[i]);
                        tf_arr = tf.array();
                        bctype = -1;
                    }
                    // In normal code operation only the fluxes at internal
                    // Dirichlet boundaries are used.  Some diagnostics use the
                    // fluxes computed at domain boundaries but these do not
                    // influence the evolution of the interior solution.
                    Array4<Real const> sp_arr;
                    if (SPa != 0) {
                        sp_arr = (*SPa)[si].array();
                    }
                    HABEC::hbflx3(Flux[idim][si].array(),
                                  Soln[si].array(icomp),
                                  reg,
                                  cdir, bctype,
                                  tf_arr,
                                  bho, bcl,
                                  fs.array(bdcomp),
                                  msk.array(),
                                  (*bcoefs[idim])[si].array(),
                                  beta, dx, flux_factor,
                                  oitr(), geom.data(), inhom,
                                  sp_arr);
                }
                else {
                    HABEC::hbflx(Flux[idim][si].array(),
                                 Soln[si].array(icomp),
                                 reg,
                                 cdir, bct, bho, bcl,
                                 fs.array(bdcomp),
                                 msk.array(),
                                 (*bcoefs[idim])[si].array(),
                                 beta, dx, inhom);
                }
            }
        }
    }
}

void HypreABec::getFaceMetric(Vector<Real>& r,
                              const Box& reg,
                              const Orientation& ori,
                              const Geometry& geom)
{
  if (ori.coordDir() == 0) {
    if (geom.IsCartesian()) {
      r.resize(1, 1.0);
    }
    else { // RZ or Spherical
      r.resize(1);
      if (ori.isLow()) {
        r[0] = geom.LoEdge(reg.smallEnd(0), 0);
      }
      else {
        r[0] = geom.HiEdge(reg.bigEnd(0), 0);
      }
      if (geom.IsSPHERICAL()) {
        r[0] *= r[0];
      }
    }
  }
  else {
    if (geom.IsCartesian()) {
      r.resize(reg.length(0), 1.0);
    }
    else { // RZ
      // We only support spherical coordinates in 1D
      BL_ASSERT(geom.IsRZ());
      geom.GetCellLoc(r, reg, 0);
    }
  }
}

void HypreABec::hacoef (const Box& bx,
                        Array4<GpuArray<Real, AMREX_SPACEDIM+1>> const& mat,
                        Array4<Real const> const& a,
                        Real alpha)
{
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        if (alpha == 0.e0_rt) {
            mat(i,j,k)[AMREX_SPACEDIM] = 0.e0_rt;
        }
        else {
            mat(i,j,k)[AMREX_SPACEDIM] = alpha * a(i,j,k);
        }
    });

    Gpu::synchronize();
}

void HypreABec::hbcoef (const Box& bx,
                        Array4<GpuArray<Real, AMREX_SPACEDIM+1>> const& mat,
                        Array4<Real const> const& b,
                        Real beta, const Real* dx,
                        int idir)
{

    if (idir == 0) {

        const Real fac = beta / (dx[0] * dx[0]);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            mat(i,j,k)[0] = -fac * b(i,j,k);
            mat(i,j,k)[AMREX_SPACEDIM] += fac * (b(i,j,k) + b(i+1,j,k));
        });

    }
    else if (idir == 1) {

        const Real fac = beta / (dx[1] * dx[1]);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            mat(i,j,k)[1] = -fac * b(i,j,k);
            mat(i,j,k)[AMREX_SPACEDIM] += fac * (b(i,j,k) + b(i,j+1,k));
        });

    }
    else {

        const Real fac = beta / (dx[2] * dx[2]);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            mat(i,j,k)[2] = -fac * b(i,j,k);
            mat(i,j,k)[AMREX_SPACEDIM] += fac * (b(i,j,k) + b(i,j,k+1));
        });

    }

    Gpu::synchronize();
}

void HypreABec::hbmat (const Box& bx,
                       Array4<GpuArray<Real, AMREX_SPACEDIM+1>> const& mat,
                       int cdir, int bct, Real bcl,
                       Array4<int const> const& mask,
                       Array4<Real const> const& b,
                       Real beta, const Real* dx)
{
    bool xlo = false;
    bool ylo = false;
    bool zlo = false;

    bool xhi = false;
    bool yhi = false;
    bool zhi = false;

    Real h;

    if (AMREX_SPACEDIM == 1) {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            xhi = true;
            h = dx[0];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }
    else if (AMREX_SPACEDIM == 2) {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 2) {
            xhi = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            ylo = true;
            h = dx[1];
        }
        else if (cdir == 3) {
            yhi = true;
            h = dx[1];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }
    else {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 3) {
            xhi = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            ylo = true;
            h = dx[1];
        }
        else if (cdir == 4) {
            yhi = true;
            h = dx[1];
        }
        else if (cdir == 2) {
            zlo = true;
            h = dx[2];
        }
        else if (cdir == 5) {
            zhi = true;
            h = dx[2];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }

    const Real fac = beta / (h * h);;

    Real bfv, bfm;

    if (bct == AMREX_LO_DIRICHLET) {
        bfv = fac * h / (0.5e0_rt * h + bcl);
        bfm = bfv - fac;
    }
    else if (bct == AMREX_LO_NEUMANN) {
        bfv = beta / h;
        bfm = -fac;
    }
    else {
        amrex::Error("hbmat: unsupported boundary type");
    }

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        if (mask.contains(i-1,j,k)) {

            if (xlo && mask(i-1,j,k) > 0) {
                mat(i,j,k)[AMREX_SPACEDIM] += bfm * b(i,j,k);
                mat(i,j,k)[0] = 0.e0_rt;
            }

        }
        else if (mask.contains(i+1,j,k)) {

            if (xhi && mask(i+1,j,k) > 0) {
                mat(i,j,k)[AMREX_SPACEDIM] += bfm * b(i+1,j,k);
            }

        }
        else if (mask.contains(i,j-1,k)) {

            if (ylo && mask(i,j-1,k) > 0) {
                mat(i,j,k)[AMREX_SPACEDIM] += bfm * b(i,j,k);
                mat(i,j,k)[1] = 0.e0_rt;
            }

        }
        else if (mask.contains(i,j+1,k)) {

            if (yhi && mask(i,j+1,k) > 0) {
                mat(i,j,k)[AMREX_SPACEDIM] += bfm * b(i,j+1,k);
            }

        }
        else if (mask.contains(i,j,k-1)) {

            if (zlo && mask(i,j,k-1) > 0) {
                mat(i,j,k)[AMREX_SPACEDIM] += bfm * b(i,j,k);
                mat(i,j,k)[2] = 0.e0_rt;
            }

        }
        else if (mask.contains(i,j,k+1)) {

            if (zhi && mask(i,j,k+1) > 0) {
                mat(i,j,k)[3] += bfm * b(i,j,k+1);
            }

        }
    });

    Gpu::synchronize();
}

void HypreABec::hbmat3 (const Box& bx,
                        int ori_lo, int idir,
                        Array4<GpuArray<Real, AMREX_SPACEDIM+1>> const& mat,
                        int cdir, int bctype,
                        Array4<int const> const& tf,
                        Real bcl,
                        Array4<int const> const& mask,
                        Array4<Real const> const& b,
                        Real beta, const GeometryData& geomdata,
                        Real c,
                        Array4<Real const> const& spa)
{
    bool xlo = false;
    bool ylo = false;
    bool zlo = false;

    bool xhi = false;
    bool yhi = false;
    bool zhi = false;

    const auto dx = geomdata.CellSize();

    Real h;

    if (AMREX_SPACEDIM == 1) {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            xhi = true;
            h = dx[0];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }
    else if (AMREX_SPACEDIM == 2) {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 2) {
            xhi = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            ylo = true;
            h = dx[1];
        }
        else if (cdir == 3) {
            yhi = true;
            h = dx[1];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }
    else {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 3) {
            xhi = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            ylo = true;
            h = dx[1];
        }
        else if (cdir == 4) {
            yhi = true;
            h = dx[1];
        }
        else if (cdir == 2) {
            zlo = true;
            h = dx[2];
        }
        else if (cdir == 5) {
            zhi = true;
            h = dx[2];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }

    const Real fac = beta / (h * h);

    // The -fac * b(i,j,k) term applied to the matrix diagonal is the contribution
    // from the interior stencil which must be removed at the boundary.

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        Real r;
        face_metric(i, j, k, bx.loVect()[0], bx.hiVect()[0], geomdata, idir, ori_lo, r);

        int bct;
        Real bfm, bfv;

        if (mask.contains(i-1,j,k)) {

            if (xlo && mask(i-1,j,k) > 0) {

                if (bctype == -1) {
                    bct = tf(i-1,j,k);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {
                    bfv = fac * h / (0.5e0_rt * h + bcl);
                    bfm = bfv * b(i,j,k);
                }
                else if (bct == AMREX_LO_NEUMANN) {
                    bfm = 0.e0_rt;
                }
                else if (bct == AMREX_LO_MARSHAK) {
                    bfv = 2.e0_rt * beta * r / h;
                    bfm = 0.25e0_rt * c * bfv;
                }
                else if (bct == AMREX_LO_SANCHEZ_POMRANING) {
                    bfv = 2.e0_rt * beta * r / h;
                    bfm = spa(i,j,k) * c * bfv;
                }
#ifndef AMREX_USE_GPU
                else {
                    amrex::Error("hbmat3: unsupported boundary type");
                }
#endif

                mat(i,j,k)[AMREX_SPACEDIM] += bfm - fac * b(i,j,k);
                mat(i,j,k)[0] = 0.e0_rt;

            }

        }
        else if (mask.contains(i+1,j,k)) {

            if (xhi && mask(i+1,j,k) > 0) {

                if (bctype == -1) {
                    bct = tf(i+1,j,k);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {
                    bfv = fac * h / (0.5e0_rt * h + bcl);
                    bfm = bfv * b(i+1,j,k);
                }
                else if (bct == AMREX_LO_NEUMANN) {
                    bfm = 0.e0_rt;
                }
                else if (bct == AMREX_LO_MARSHAK) {
                    bfv = 2.e0_rt * beta * r / h;
                    bfm = 0.25e0_rt * c * bfv;
                }
                else if (bct == AMREX_LO_SANCHEZ_POMRANING) {
                    bfv = 2.e0_rt * beta * r / h;
                    bfm = spa(i,j,k) * c * bfv;
                }
#ifndef AMREX_USE_GPU
                else {
                    amrex::Error("hbmat3: unsupported boundary type");
                }
#endif

                mat(i,j,k)[AMREX_SPACEDIM] += bfm - fac * b(i+1,j,k);

            }

        }
        else if (mask.contains(i,j-1,k)) {

            if (ylo && mask(i,j-1,k) > 0) {

                if (bctype == -1) {
                    bct = tf(i,j-1,k);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {
                    bfv = fac * h / (0.5e0_rt * h + bcl);
                    bfm = bfv * b(i,j,k);
                }
                else if (bct == AMREX_LO_NEUMANN) {
                    bfm = 0.e0_rt;
                }
                else if (bct == AMREX_LO_MARSHAK) {
                    bfv = 2.e0_rt * beta * r / h;
                    bfm = 0.25e0_rt * c * bfv;
                }
                else if (bct == AMREX_LO_SANCHEZ_POMRANING) {
                    bfv = 2.e0_rt * beta * r / h;
                    bfm = spa(i,j,k) * c * bfv;
                }
#ifndef AMREX_USE_GPU
                else {
                    amrex::Error("hbmat3: unsupported boundary type");
                }
#endif

                mat(i,j,k)[AMREX_SPACEDIM] += bfm - fac * b(i,j,k);
                mat(i,j,k)[1] = 0.e0_rt;

            }

        }
        else if (mask.contains(i,j+1,k)) {

            if (yhi && mask(i,j+1,k) > 0) {

                if (bctype == -1) {
                    bct = tf(i,j+1,k);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {
                    bfv = fac * h / (0.5e0_rt * h + bcl);
                    bfm = bfv * b(i,j+1,k);
                }
                else if (bct == AMREX_LO_NEUMANN) {
                    bfm = 0.e0_rt;
                }
                else if (bct == AMREX_LO_MARSHAK) {
                    bfv = 2.e0_rt * beta * r / h;
                    bfm = 0.25e0_rt * c * bfv;
                }
                else if (bct == AMREX_LO_SANCHEZ_POMRANING) {
                    bfv = 2.e0_rt * beta * r / h;
                    bfm = spa(i,j,k) * c * bfv;
                }
#ifndef AMREX_USE_GPU
                else {
                    amrex::Error("hbmat3: unsupported boundary type");
                }
#endif

                mat(i,j,k)[AMREX_SPACEDIM] += bfm - fac * b(i,j+1,k);

            }


        }
        else if (mask.contains(i,j,k-1)) {

            if (zlo && mask(i,j,k-1) > 0) {

                if (bctype == -1) {
                    bct = tf(i,j,k-1);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {
                    bfv = fac * h / (0.5e0_rt * h + bcl);
                    bfm = bfv * b(i,j,k);
                }
                else if (bct == AMREX_LO_NEUMANN) {
                    bfm = 0.e0_rt;
                }
                else if (bct == AMREX_LO_MARSHAK) {
                    bfv = 2.e0_rt * beta * r / h;
                    bfm = 0.25e0_rt * c * bfv;
                }
                else if (bct == AMREX_LO_SANCHEZ_POMRANING) {
                    bfv = 2.e0_rt * beta * r / h;
                    bfm = spa(i,j,k) * c * bfv;
                }
#ifndef AMREX_USE_GPU
                else {
                    amrex::Error("hbmat3: unsupported boundary type");
                }
#endif

                mat(i,j,k)[AMREX_SPACEDIM] += bfm - fac * b(i,j,k);
                mat(i,j,k)[2] = 0.e0_rt;

            }

        }
        else if (mask.contains(i,j,k+1)) {

            if (zhi && mask(i,j,k+1) > 0) {

                if (bctype == -1) {
                    bct = tf(i,j,k+1);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {
                    bfv = fac * h / (0.5e0_rt * h + bcl);
                    bfm = bfv * b(i,j,k+1);
                }
                else if (bct == AMREX_LO_NEUMANN) {
                    bfm = 0.e0_rt;
                }
                else if (bct == AMREX_LO_MARSHAK) {
                    bfv = 2.e0_rt * beta * r / h;
                    bfm = 0.25e0_rt * c * bfv;
                }
                else if (bct == AMREX_LO_SANCHEZ_POMRANING) {
                    bfv = 2.e0_rt * beta * r / h;
                    bfm = spa(i,j,k) * c * bfv;
                }
#ifndef AMREX_USE_GPU
                else {
                    amrex::Error("hbmat3: unsupported boundary type");
                }
#endif

                mat(i,j,k)[AMREX_SPACEDIM] += bfm - fac * b(i,j,k+1);

            }

        }
    });

    Gpu::synchronize();
}

void HypreABec::setupSolver(Real _reltol, Real _abstol, int maxiter)
{
  BL_PROFILE("HypreABec::setupSolver");

  const BoxArray& grids = acoefs->boxArray();

  const int size = AMREX_SPACEDIM + 1;
  int i, idim;

  int stencil_indices[size];

  for (i = 0; i < size; i++) {
    stencil_indices[i] = i;
  }

  BaseFab<GpuArray<Real, size>> matfab; // AoS indexing
  for (MFIter ai(*acoefs); ai.isValid(); ++ai) {
    i = ai.index();
    const Box &reg = grids[i];

    matfab.resize(reg);
    Elixir matfab_elix = matfab.elixir();
    Real* mat = (Real*) matfab.dataPtr();

    // build matrix interior

    hacoef(reg, matfab.array(), (*acoefs)[ai].array(), alpha);

    for (idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        hbcoef(reg, matfab.array(), (*bcoefs[idim])[ai].array(), beta, dx, idim);
    }

    // add b.c.'s to matrix diagonal, and
    // zero out offdiag values at low domain boundaries (high done by symmetry)

    const NGBndry& bd = getBndry();
    const Box& domain = bd.getDomain();
    for (OrientationIter oitr; oitr; oitr++) {
      int cdir(oitr());
      idim = oitr().coordDir();
      const RadBoundCond &bct = bd.bndryConds(oitr())[i];
      const Real      &bcl = bd.bndryLocs(oitr())[i];
      const Mask      &msk = bd.bndryMasks(oitr(),i);
      const Box &bbox = (*bcoefs[idim])[ai].box();
      const Box &msb  = msk.box();
      if (reg[oitr()] == domain[oitr()]) {
        Array4<int const> tfp{};
        int bctype = bct;
        if (bd.mixedBndry(oitr())) {
          const BaseFab<int> &tf = *(bd.bndryTypes(oitr())[i]);
          tfp = tf.array();
          bctype = -1;
        }
        const Box &fsb = bd.bndryValues(oitr())[ai].box();
        Array4<Real const> pSPa{};
        if (SPa != 0) {
          pSPa = (*SPa)[ai].array();
        }

        hbmat3(reg, oitr().isLow(), idim,
               matfab.array(), cdir, bctype,
               tfp, bcl, msk.array(),
               (*bcoefs[idim])[ai].array(),
               beta, geom.data(), flux_factor,
               pSPa);
      }
      else {
          hbmat(reg, matfab.array(), cdir, bct, bcl, msk.array(), (*bcoefs[idim])[ai].array(), beta, dx);
      }
    }

    // initialize matrix

    HYPRE_StructMatrixSetBoxValues(A, loV(reg), hiV(reg),
                                   size, stencil_indices, mat);
    Gpu::synchronize();
  }

  HYPRE_StructMatrixAssemble(A);

  HYPRE_StructVectorAssemble(b); // currently a no-op
  HYPRE_StructVectorAssemble(x); // currently a no-op

  reltol = _reltol;
  abstol = _abstol; // may be used to change tolerance for solve

  if (solver_flag == 0) {
    HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_StructSMGSetMemoryUse(solver, 0);
    HYPRE_StructSMGSetMaxIter(solver, maxiter);
    HYPRE_StructSMGSetRelChange(solver, 0);
    HYPRE_StructSMGSetTol(solver, reltol);
    HYPRE_StructSMGSetNumPreRelax(solver, 1);
    HYPRE_StructSMGSetNumPostRelax(solver, 1);
    HYPRE_StructSMGSetLogging(solver, 1);
    HYPRE_StructSMGSetup(solver, A, b, x);
  }
  else if (solver_flag == 1) {
    HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &solver);
    //HYPRE_StructPFMGSetMemoryUse(solver, 0);
    HYPRE_StructPFMGSetSkipRelax(solver, 0);
    HYPRE_StructPFMGSetMaxIter(solver, maxiter);
    HYPRE_StructPFMGSetRelChange(solver, 0);
    HYPRE_StructPFMGSetTol(solver, reltol);
// in following line, Falgout says use 1 as relax type, not 2 (rbp, 9/27/05)
// weighted Jacobi = 1; red-black GS = 2
    HYPRE_StructPFMGSetRelaxType(solver, pfmg_relax_type);
    HYPRE_StructPFMGSetNumPreRelax(solver, 1);
    HYPRE_StructPFMGSetNumPostRelax(solver, 1);
    HYPRE_StructPFMGSetLogging(solver, 1);
    HYPRE_StructPFMGSetup(solver, A, b, x);
  }
  else if (solver_flag == 2) {
    HYPRE_StructJacobiCreate(MPI_COMM_WORLD, &solver);
    //HYPRE_StructPFMGSetMemoryUse(solver, 0);
    //HYPRE_StructPFMGSetSkipRelax(solver, 0);
    HYPRE_StructJacobiSetMaxIter(solver, maxiter);
    //HYPRE_StructPFMGSetRelChange(solver, 0);
    //HYPRE_StructPFMGSetTol(solver, reltol);
    //HYPRE_StructPFMGSetNumPreRelax(solver, 1);
    //HYPRE_StructPFMGSetNumPostRelax(solver, 1);
    //HYPRE_StructPFMGSetLogging(solver, 1);
    HYPRE_StructJacobiSetup(solver, A, b, x);
  }
  else if (solver_flag == 3 || solver_flag == 4) {
    HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_StructPCGSetMaxIter(solver, maxiter);
    HYPRE_StructPCGSetRelChange(solver, 0);
    HYPRE_StructPCGSetTol(solver, reltol);

    if (solver_flag == 3) {
// pfmg pre-conditioned cg
      HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &precond);
      HYPRE_StructPFMGSetMaxIter(precond, 1);
      HYPRE_StructPFMGSetTol(precond, 0.0);
      HYPRE_StructPFMGSetZeroGuess(precond);
// weighted Jacobi = 1; red-black GS = 2
      HYPRE_StructPFMGSetRelaxType(precond, pfmg_relax_type);
      HYPRE_StructPFMGSetNumPreRelax(precond, 1);
      HYPRE_StructPFMGSetNumPostRelax(precond, 1);
      HYPRE_StructPFMGSetSkipRelax(precond, 0);
      HYPRE_StructPFMGSetLogging(precond, 0);
      HYPRE_StructPCGSetPrecond(solver,
                                HYPRE_StructPFMGSolve,
                                HYPRE_StructPFMGSetup,
                                precond);
    }
    else if (solver_flag == 4) {
      HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
      HYPRE_StructSMGSetMemoryUse(precond, 0);
      HYPRE_StructSMGSetMaxIter(precond, 1);
      HYPRE_StructSMGSetRelChange(precond, 0);
      HYPRE_StructSMGSetTol(precond, 0.0);
      HYPRE_StructSMGSetNumPreRelax(precond, 1);
      HYPRE_StructSMGSetNumPostRelax(precond, 1);
      HYPRE_StructSMGSetLogging(precond, 0);
      HYPRE_StructPCGSetPrecond(solver,
                                HYPRE_StructSMGSolve,
                                HYPRE_StructSMGSetup,
                                precond);
    }

    HYPRE_StructPCGSetLogging(solver, 1);
    HYPRE_StructPCGSetup(solver, A, b, x);
  }
  else if (solver_flag == 5 || solver_flag == 6) {
    HYPRE_StructHybridCreate(MPI_COMM_WORLD, &solver);
    HYPRE_StructHybridSetDSCGMaxIter(solver, maxiter);
    HYPRE_StructHybridSetPCGMaxIter(solver, maxiter);
    HYPRE_StructHybridSetTol(solver, reltol);

    HYPRE_StructHybridSetConvergenceTol(solver, 0.90);
    HYPRE_StructHybridSetRelChange(solver, 0);
    HYPRE_StructHybridSetLogging(solver, 1);
    HYPRE_StructHybridSetSolverType(solver, 1); /* pcg */

    /* pfmg preconditioning */
    if (solver_flag == 5) {
      HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &precond);
      HYPRE_StructPFMGSetMaxIter(precond, 1);
      HYPRE_StructPFMGSetTol(precond, 0.0);
      HYPRE_StructPFMGSetZeroGuess(precond);
// weighted Jacobi = 1; red-black GS = 2
      HYPRE_StructPFMGSetRelaxType(precond, pfmg_relax_type);
      HYPRE_StructPFMGSetNumPreRelax(precond, 1);
      HYPRE_StructPFMGSetNumPostRelax(precond, 1);
      HYPRE_StructPFMGSetSkipRelax(precond, 0);
      HYPRE_StructPFMGSetLogging(precond, 0);
      HYPRE_StructHybridSetPrecond(solver,
                                   HYPRE_StructPFMGSolve,
                                   HYPRE_StructPFMGSetup,
                                   precond);
    }
    else if (solver_flag == 6) {
      HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
      HYPRE_StructSMGSetMemoryUse(precond, 0);
      HYPRE_StructSMGSetMaxIter(precond, 1);
      HYPRE_StructSMGSetRelChange(precond, 0);
      HYPRE_StructSMGSetTol(precond, 0.0);
      HYPRE_StructSMGSetNumPreRelax(precond, 1);
      HYPRE_StructSMGSetNumPostRelax(precond, 1);
      HYPRE_StructSMGSetLogging(precond, 0);
      HYPRE_StructHybridSetPrecond(solver,
                                   HYPRE_StructSMGSolve,
                                   HYPRE_StructSMGSetup,
                                   precond);
    }

    HYPRE_StructHybridSetup(solver, A, b, x);
  }
  else {
      amrex::Error("HypreABec: no such solver");
  }
  Gpu::synchronize();
}

void HypreABec::clearSolver()
{
  BL_PROFILE("HypreABec::clearSolver");

  if (solver_flag == 0) {
    HYPRE_StructSMGDestroy(solver);
  }
  else if (solver_flag == 1) {
    HYPRE_StructPFMGDestroy(solver);
  }
  else if(solver_flag == 2) {
    HYPRE_StructJacobiDestroy(solver);
  }
  else if(solver_flag == 3 || solver_flag == 4) {
    HYPRE_StructPCGDestroy(solver);
    if (solver_flag == 3)
    {
       HYPRE_StructPFMGDestroy(precond);
    }
    else if (solver_flag == 4)
    {
       HYPRE_StructSMGDestroy(precond);
    }
  }
  else if(solver_flag == 5 || solver_flag == 6) {
    HYPRE_StructHybridDestroy(solver);
    if(solver_flag == 5) {
       HYPRE_StructPFMGDestroy(precond);
    }
    if(solver_flag == 6) {
       HYPRE_StructSMGDestroy(precond);
    }
  }
}

void HypreABec::hbvec3 (const Box& bx,
                        int ori_lo, int idir,
                        Array4<Real> const& vec,
                        int cdir, int bctype,
                        Array4<int const> const& tf,
                        int bho, Real bcl,
                        Array4<Real const> const& bcval,
                        Array4<int const> const& mask,
                        Array4<Real const> const& b,
                        Real beta, const GeometryData& geomdata)
{
    bool xlo = false;
    bool ylo = false;
    bool zlo = false;

    bool xhi = false;
    bool yhi = false;
    bool zhi = false;

    Real h;

    const auto dx = geomdata.CellSize();

    if (AMREX_SPACEDIM == 1) {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            xhi = true;
            h = dx[0];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }
    else if (AMREX_SPACEDIM == 2) {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 2) {
            xhi = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            ylo = true;
            h = dx[1];
        }
        else if (cdir == 3) {
            yhi = true;
            h = dx[1];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }
    else {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 3) {
            xhi = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            ylo = true;
            h = dx[1];
        }
        else if (cdir == 4) {
            yhi = true;
            h = dx[1];
        }
        else if (cdir == 2) {
            zlo = true;
            h = dx[2];
        }
        else if (cdir == 5) {
            zhi = true;
            h = dx[2];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        Real r;
        face_metric(i, j, k, bx.loVect()[0], bx.hiVect()[0], geomdata, idir, ori_lo, r);

        int bct;
        Real bfv, h2, th2;

        if (mask.contains(i-1,j,k)) {

            if (xlo && mask(i-1,j,k) > 0) {

                if (bctype == -1) {
                    bct = tf(i-1,j,k);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {
                    if (bho >= 1) {
                        h2 = 0.5e0_rt * h;
                        th2 = 3.e0_rt * h2;
                        bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2));
                    }
                    else {
                        bfv = (beta / h) / (0.5e0_rt * h + bcl);
                    }

                    bfv = bfv * b(i,j,k);
                }
                else if (bct == AMREX_LO_NEUMANN) {
                    bfv = beta * r / h;
                }
                else if (bct == AMREX_LO_MARSHAK || bct == AMREX_LO_SANCHEZ_POMRANING) {
                    bfv = 2.e0_rt * beta * r / h;
                }
#ifndef AMREX_USE_GPU
                else {
                    amrex::Error("hbvec3: unsupported boundary type");
                }
#endif

                vec(i,j,k) = vec(i,j,k) + bfv * bcval(i-1,j,k);

            }

        }
        else if (mask.contains(i+1,j,k)) {

            if (xhi && mask(i+1,j,k) > 0) {

                if (bctype == -1) {
                    bct = tf(i+1,j,k);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {
                    if (bho >= 1) {
                        h2 = 0.5e0_rt * h;
                        th2 = 3.e0_rt * h2;
                        bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2));
                    }
                    else {
                        bfv = (beta / h) / (0.5e0_rt * h + bcl);
                    }

                    bfv = bfv * b(i+1,j,k);
                }
                else if (bct == AMREX_LO_NEUMANN) {
                    bfv = beta * r / h;
                }
                else if (bct == AMREX_LO_MARSHAK || bct == AMREX_LO_SANCHEZ_POMRANING) {
                    bfv = 2.e0_rt * beta * r / h;
                }
#ifndef AMREX_USE_GPU
                else {
                    amrex::Error("hbvec3: unsupported boundary type");
                }
#endif

                vec(i,j,k) = vec(i,j,k) + bfv * bcval(i+1,j,k);

            }

        }
        else if (mask.contains(i,j-1,k)) {

            if (ylo && mask(i,j-1,k) > 0) {

                if (bctype == -1) {
                    bct = tf(i,j-1,k);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {
                    if (bho >= 1) {
                        h2 = 0.5e0_rt * h;
                        th2 = 3.e0_rt * h2;
                        bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2));
                    }
                    else {
                        bfv = (beta / h) / (0.5e0_rt * h + bcl);
                    }

                    bfv = bfv * b(i,j,k);
                }
                else if (bct == AMREX_LO_NEUMANN) {
                    bfv = beta * r / h;
                }
                else if (bct == AMREX_LO_MARSHAK || bct == AMREX_LO_SANCHEZ_POMRANING) {
                    bfv = 2.e0_rt * beta * r / h;
                }
#ifndef AMREX_USE_GPU
                else {
                    amrex::Error("hbvec3: unsupported boundary type");
                }
#endif

                vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j-1,k);

            }

        }
        else if (mask.contains(i,j+1,k)) {

            if (yhi && mask(i,j+1,k) > 0) {

                if (bctype == -1) {
                    bct = tf(i,j+1,k);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {
                    if (bho >= 1) {
                        h2 = 0.5e0_rt * h;
                        th2 = 3.e0_rt * h2;
                        bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2));
                    }
                    else {
                        bfv = (beta / h) / (0.5e0_rt * h + bcl);
                    }

                    bfv = bfv * b(i,j+1,k);
                }
                else if (bct == AMREX_LO_NEUMANN) {
                    bfv = beta * r / h;
                }
                else if (bct == AMREX_LO_MARSHAK || bct == AMREX_LO_SANCHEZ_POMRANING) {
                    bfv = 2.e0_rt * beta * r / h;
                }
#ifndef AMREX_USE_GPU
                else {
                    amrex::Error("hbvec3: unsupported boundary type");
                }
#endif

                vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j+1,k);

            }

        }
        else if (mask.contains(i,j,k-1)) {

            if (zlo && mask(i,j,k-1) > 0) {

                if (bctype == -1) {
                    bct = tf(i,j,k-1);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {
                    if (bho >= 1) {
                        h2 = 0.5e0_rt * h;
                        th2 = 3.e0_rt * h2;
                        bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2));
                    }
                    else {
                        bfv = (beta / h) / (0.5e0_rt * h + bcl);
                    }

                    bfv = bfv * b(i,j,k);
                }
                else if (bct == AMREX_LO_NEUMANN) {
                    bfv = beta * r / h;
                }
                else if (bct == AMREX_LO_MARSHAK || bct == AMREX_LO_SANCHEZ_POMRANING) {
                    bfv = 2.e0_rt * beta * r / h;
                }
#ifndef AMREX_USE_GPU
                else {
                    amrex::Error("hbvec3: unsupported boundary type");
                }
#endif

                vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j,k-1);

            }

        }
        else if (mask.contains(i,j,k+1)) {

            if (zhi && mask(i,j,k+1) > 0) {

                if (bctype == -1) {
                    bct = tf(i,j,k+1);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {
                    if (bho >= 1) {
                        h2 = 0.5e0_rt * h;
                        th2 = 3.e0_rt * h2;
                        bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2));
                    }
                    else {
                        bfv = (beta / h) / (0.5e0_rt * h + bcl);
                    }

                    bfv = bfv * b(i,j,k+1);
                }
                else if (bct == AMREX_LO_NEUMANN) {
                    bfv = beta * r / h;
                }
                else if (bct == AMREX_LO_MARSHAK || bct == AMREX_LO_SANCHEZ_POMRANING) {
                    bfv = 2.e0_rt * beta * r / h;
                }
#ifndef AMREX_USE_GPU
                else {
                    amrex::Error("hbvec3: unsupported boundary type");
                }
#endif

                vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j,k+1);

            }

        }
    });

    Gpu::synchronize();
}

void HypreABec::hbvec (const Box& bx,
                       Array4<Real> const& vec,
                       int cdir, int bct, int bho, Real bcl,
                       Array4<Real const> const& bcval,
                       Array4<int const> const& mask,
                       Array4<Real const> const& b,
                       Real beta, const Real* dx)
{
    bool xlo = false;
    bool ylo = false;
    bool zlo = false;

    bool xhi = false;
    bool yhi = false;
    bool zhi = false;

    Real h, bfv, h2, th2;

    if (AMREX_SPACEDIM == 1) {
        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            xhi = true;
            h = dx[0];
        }
        else {
            amrex::Error("Unknown cdir");
        }
    }
    else if (AMREX_SPACEDIM == 2) {
        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 2) {
            xhi = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            ylo = true;
            h = dx[1];
        }
        else if (cdir == 3) {
            yhi = true;
            h = dx[1];
        }
        else {
            amrex::Error("Unknown cdir");
        }
    }
    else {
        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 3) {
            xhi = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            ylo = true;
            h = dx[1];
        }
        else if (cdir == 4) {
            yhi = true;
            h = dx[1];
        }
        else if (cdir == 2) {
            zlo = true;
            h = dx[2];
        }
        else if (cdir == 5) {
            zhi = true;
            h = dx[2];
        }
        else {
            amrex::Error("Unknown cdir");
        }
    }

    if (bct == AMREX_LO_DIRICHLET) {
        if (bho >= 1) {
            h2 = 0.5e0_rt * h;
            th2 = 3.e0_rt * h2;
            bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2));
        }
        else {
            bfv = (beta / h) / (0.5e0_rt * h + bcl);
        }
    }
    else if (bct == AMREX_LO_NEUMANN) {
        bfv = beta / h;
    }
    else {
        amrex::Error("hbvec: unsupported boundary type");
    }

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        if (mask.contains(i-1,j,k)) {

            if (xlo && mask(i-1,j,k) > 0) {

                vec(i,j,k) = vec(i,j,k) + bfv * b(i,j,k) * bcval(i-1,j,k);

            }

        }
        else if (mask.contains(i+1,j,k)) {

            if (xhi && mask(i+1,j,k) > 0) {

                vec(i,j,k) = vec(i,j,k) + bfv * b(i+1,j,k) * bcval(i+1,j,k);

            }

        }
        else if (mask.contains(i,j-1,k)) {

            if (ylo && mask(i,j-1,k) > 0) {

                vec(i,j,k) = vec(i,j,k) + bfv * b(i,j,k) * bcval(i,j-1,k);

            }

        }
        else if (mask.contains(i,j+1,k)) {

            if (yhi && mask(i,j+1,k) > 0) {

                vec(i,j,k) = vec(i,j,k) + bfv * b(i,j+1,k) * bcval(i,j+1,k);

            }

        }
        else if (mask.contains(i,j,k-1)) {

            if (zlo && mask(i,j,k-1) > 0) {

                vec(i,j,k) = vec(i,j,k) + bfv * b(i,j,k) * bcval(i,j,k-1);

            }

        }
        else if (mask.contains(i,j,k+1)) {

            if (zhi && mask(i,j,k+1) > 0) {

                vec(i,j,k) = vec(i,j,k) + bfv * b(i,j,k+1) * bcval(i,j,k+1);

            }

        }

    });

    Gpu::synchronize();
}

void HypreABec::solve(MultiFab& dest, int icomp, MultiFab& rhs, BC_Mode inhom)
{
  BL_PROFILE("HypreABec::solve");

  const BoxArray& grids = dest.boxArray();

  int i, idim;

  Real *vec;
  FArrayBox fnew;
  for (MFIter di(dest); di.isValid(); ++di) {
    i = di.index();
    const Box &reg = grids[i];

    // initialize dest, since we will reuse the space to set up rhs below:

    FArrayBox *f;
    int fcomp;

    Array4<Real> const d_arr = dest.array(di);
    Array4<Real> const r_arr = rhs.array(di);

    if (dest.nGrow() == 0) { // need a temporary if dest is the wrong size
      f = &dest[di];
      fcomp = icomp;
    }
    else {
      f = &fnew;
      f->resize(reg);

      Array4<Real> const f_arr = f->array();

      fcomp = 0;

      AMREX_PARALLEL_FOR_3D(reg, i, j, k, { f_arr(i,j,k,fcomp) = d_arr(i,j,k,icomp); });
      Gpu::streamSynchronize();
    }
    Elixir f_elix = fnew.elixir();

    vec = f->dataPtr(fcomp); // sharing space, dest will be overwritten below

    HYPRE_StructVectorSetBoxValues(x, loV(reg), hiV(reg), vec);

    Gpu::streamSynchronize();

    Array4<Real> const f_arr = f->array();

    AMREX_PARALLEL_FOR_3D(reg, i, j, k, { f_arr(i,j,k,fcomp) = r_arr(i,j,k,0); });
    Gpu::streamSynchronize();

    // add b.c.'s to rhs

    if (inhom) {
      const NGBndry& bd = getBndry();
      const Box& domain = bd.getDomain();
      for (OrientationIter oitr; oitr; oitr++) {
        int cdir(oitr());
        idim = oitr().coordDir();
        const RadBoundCond &bct = bd.bndryConds(oitr())[i];
        const Real      &bcl = bd.bndryLocs(oitr())[i];
        const FArrayBox       &fs  = bd.bndryValues(oitr())[di];
        const Mask      &msk = bd.bndryMasks(oitr(),i);
        const Box &bbox = (*bcoefs[idim])[di].box();

        if (reg[oitr()] == domain[oitr()]) {
          Array4<const int> tfp{};
          int bctype = bct;
          if (bd.mixedBndry(oitr())) {
            const BaseFab<int> &tf = *(bd.bndryTypes(oitr())[i]);
            tfp = tf.array();
            bctype = -1;
          }
          hbvec3(reg,
                 oitr().isLow(), idim,
                 f->array(fcomp),
                 cdir, bct,
                 tfp,
                 bho, bcl,
                 fs.array(bdcomp),
                 msk.array(),
                 (*bcoefs[idim])[di].array(),
                 beta, geom.data());
        }
        else {
            hbvec(reg, f->array(fcomp),
                  cdir, bct, bho, bcl,
                  fs.array(bdcomp), msk.array(),
                  (*bcoefs[idim])[di].array(),
                  beta, dx);
        }
      }
    }

    Gpu::streamSynchronize();

    // initialize rhs

    HYPRE_StructVectorSetBoxValues(b, loV(reg), hiV(reg), vec);
  }

  HYPRE_StructVectorAssemble(b); // currently a no-op
  HYPRE_StructVectorAssemble(x); // currently a no-op
  Gpu::synchronize();

  if (abstol > 0.0) {
    Real bnorm;
    bnorm = hypre_StructInnerProd((hypre_StructVector *) b,
                                  (hypre_StructVector *) b);
    bnorm = std::sqrt(bnorm);

    const BoxArray& grids = acoefs->boxArray();
    Real volume = 0.0;
    for (int i = 0; i < grids.size(); i++) {
      volume += grids[i].numPts();
    }

    Real reltol_new = (bnorm > 0.0
                       ? abstol / bnorm * std::sqrt(volume)
                       : reltol);

    if (reltol_new > reltol) {
      if (solver_flag == 0) {
        HYPRE_StructSMGSetTol(solver, reltol_new);
      }
      else if(solver_flag == 1) {
        HYPRE_StructPFMGSetTol(solver, reltol_new);
      }
      else if(solver_flag == 2) {
        // nothing for this option
      }
      else if(solver_flag == 3 || solver_flag == 4) {
        HYPRE_StructPCGSetTol(solver, reltol_new);
      }
    }
  }

  if (solver_flag == 0) {
    HYPRE_StructSMGSolve(solver, A, b, x);
    //HYPRE_StructVectorPrint("Xsmg", x, 0);
    //HYPRE_StructVectorPrint("Bsmg", b, 0);
    //cin.get();
  }
  else if (solver_flag == 1) {
    HYPRE_StructPFMGSolve(solver, A, b, x);
  }
  else if (solver_flag == 2) {
    HYPRE_StructJacobiSolve(solver, A, b, x);
  }
  else if (solver_flag == 3 || solver_flag == 4) {
    HYPRE_StructPCGSolve(solver, A, b, x);
  }
  else if (solver_flag == 5 || solver_flag == 6) {
    HYPRE_StructHybridSolve(solver, A, b, x);
  }

  Gpu::synchronize();

  for (MFIter di(dest); di.isValid(); ++di) {
    i = di.index();
    const Box &reg = grids[i];

    FArrayBox *f;
    int fcomp;
    if (dest.nGrow() == 0) { // need a temporary if dest is the wrong size
      f = &dest[di];
      fcomp = icomp;
    }
    else {
      f = &fnew;
      f->resize(reg);
      fcomp = 0;
    }
    Elixir f_elix = fnew.elixir();

    vec = f->dataPtr(fcomp);
    HYPRE_StructVectorGetBoxValues(x, loV(reg), hiV(reg),
                                   vec);
    Gpu::synchronize();

    if (dest.nGrow() != 0) {
        Array4<Real> const f_arr = f->array();
        Array4<Real> const d_arr = dest.array(di);
        AMREX_PARALLEL_FOR_3D(reg, i, j, k, { d_arr(i,j,k,icomp) = f_arr(i,j,k); });
    }
  }

  if (verbose >= 2 && ParallelDescriptor::IOProcessor()) {
    int num_iterations;
    Real res;
    if (solver_flag == 0) {
      HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
      HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &res);
    }
    else if(solver_flag == 1) {
      HYPRE_StructPFMGGetNumIterations(solver, &num_iterations);
      HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &res);
    }
    else if(solver_flag == 2) {
      HYPRE_StructJacobiGetNumIterations(solver, &num_iterations);
      HYPRE_StructJacobiGetFinalRelativeResidualNorm(solver, &res);
    }
    else if(solver_flag == 3 || solver_flag == 4) {
      HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
      HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &res);
    }
    else if(solver_flag == 5 || solver_flag == 6) {
      HYPRE_StructHybridGetNumIterations(solver, &num_iterations);
      HYPRE_StructHybridGetFinalRelativeResidualNorm(solver, &res);
    }
    if (num_iterations >= verbose_threshold) {
      int oldprec = std::cout.precision(20);
      std::cout << num_iterations
           << " Hypre Multigrid Iterations, Relative Residual "
           << res << std::endl;
      std::cout.precision(oldprec);
    }
  }
  Gpu::synchronize();
}

Real HypreABec::getAbsoluteResidual()
{
  BL_PROFILE("HypreABec::getAbsoluteResidual");

  Real bnorm;
  bnorm = hypre_StructInnerProd((hypre_StructVector *) b,
                                (hypre_StructVector *) b);
  bnorm = std::sqrt(bnorm);

  Real res;
  if (solver_flag == 0) {
    HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &res);
  }
  else if(solver_flag == 1) {
    HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &res);
  }
  else if(solver_flag == 2) {
    HYPRE_StructJacobiGetFinalRelativeResidualNorm(solver, &res);
  }
  else if(solver_flag == 3 || solver_flag == 4) {
    HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &res);
  }
  else if(solver_flag == 5 || solver_flag == 6) {
    HYPRE_StructHybridGetFinalRelativeResidualNorm(solver, &res);
  }

  Gpu::synchronize();

  const BoxArray& grids = acoefs->boxArray();
  Real volume = 0.0;
  for (int i = 0; i < grids.size(); i++) {
    volume += grids[i].numPts();
  }

  return bnorm * res / std::sqrt(volume);
}
