
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>

#include "HypreABec.H"
#include "HABEC_F.H"

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "_hypre_struct_mv.h"

using namespace amrex;

static int ispow2(int i)
{
  return (i == 1) ? 1 : (((i <= 0) || (i & 1)) ? 0 : ispow2(i / 2));
}

Real HypreABec::flux_factor = 1.0;

#if (BL_SPACEDIM == 1)
static int vl[2] = { 0, 0 };
static int vh[2] = { 0, 0 };
#endif

static int* loV(const Box& b) {
#if (BL_SPACEDIM == 1)
  vl[0] = b.smallEnd(0);
  return vl;
#else
  return (int*) b.loVect();
#endif
}

static int* hiV(const Box& b) {
#if (BL_SPACEDIM == 1)
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

  int i;
#if defined(BL_USE_MPI) || !(defined(BL_BGL) || defined(chaos_3_x86_64_ib) || defined(chaos_3_x86_64))
  MPI_Initialized(&i);
#else
  i=1;
#endif
  if (!i) {
    int   argc   = 1;
    const char *argv[] = { "mf" };
    // arguments must be set, though not used for anything important
    MPI_Init(&argc, (char***)&argv);
  }

  int num_procs, myid;

  MPI_Comm_size(MPI_COMM_WORLD, &num_procs );
  MPI_Comm_rank(MPI_COMM_WORLD, &myid );

  for (i = 0; i < BL_SPACEDIM; i++) {
    dx[i] = geom.CellSize(i);
  }

#if (BL_SPACEDIM == 1)

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

  HYPRE_StructGridCreate(MPI_COMM_WORLD, BL_SPACEDIM, &hgrid);

  if (geom.isAnyPeriodic()) {
    int is_periodic[BL_SPACEDIM];
    for (i = 0; i < BL_SPACEDIM; i++) {
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

  if (num_procs != 1) {
    // parallel section:
    BL_ASSERT(ParallelDescriptor::NProcs() == num_procs);
    BL_ASSERT(ParallelDescriptor::MyProc() == myid);

    for (i = 0; i < grids.size(); i++) {
      if (dmap[i] == myid) {
	HYPRE_StructGridSetExtents(hgrid, loV(grids[i]), hiV(grids[i]));
      }
    }
  }
  else {
    for (i = 0; i < grids.size(); i++) {
      HYPRE_StructGridSetExtents(hgrid, loV(grids[i]), hiV(grids[i]));
    }
  }

  HYPRE_StructGridAssemble(hgrid);

#if (BL_SPACEDIM == 1)
  // if we were really 1D:
/*
  int offsets[2][1] = {{-1},
		       { 0}};
*/
  // fake 1D as a 2D problem:
  int offsets[2][2] = {{-1,  0},
		       { 0,  0}};
#elif (BL_SPACEDIM == 2)
  int offsets[3][2] = {{-1,  0},
		       { 0, -1},
		       { 0,  0}};
#elif (BL_SPACEDIM == 3)
  int offsets[4][3] = {{-1,  0,  0},
		       { 0, -1,  0},
		       { 0,  0, -1},
		       { 0,  0,  0}};
#endif

#if   (BL_SPACEDIM == 1)
  int A_num_ghost[6] = { 1, 1, 0, 0, 0, 0 };
#elif (BL_SPACEDIM == 2)
  //int A_num_ghost[4] = { 1, 1, 1, 1 };
  int A_num_ghost[6] = { 1, 1, 1, 1, 0, 0 };
#elif (BL_SPACEDIM == 3)
  int A_num_ghost[6] = { 1, 1, 1, 1, 1, 1 };
#endif

  HYPRE_StructStencil stencil;

#if (BL_SPACEDIM == 1)
  HYPRE_StructStencilCreate(2, 2, &stencil);
#else
  HYPRE_StructStencilCreate(BL_SPACEDIM, BL_SPACEDIM + 1, &stencil);
#endif

  for (i = 0; i < BL_SPACEDIM + 1; i++) {
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
 
  for (i = 0; i < BL_SPACEDIM; i++) {
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
	Vector<Real> r;
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
		    const int *tfp = NULL;
		    int bctype = bct;
		    if (bd.mixedBndry(oitr())) {
			const BaseFab<int> &tf = *(bd.bndryTypes(oitr())[i]);
			tfp = tf.dataPtr();
			bctype = -1;
		    }
		    // In normal code operation only the fluxes at internal
		    // Dirichlet boundaries are used.  Some diagnostics use the
		    // fluxes computed at domain boundaries but these do not
		    // influence the evolution of the interior solution.
		    Real* pSPa;
		    Box SPabox; 
		    if (SPa != 0) {
			pSPa = (*SPa)[si].dataPtr();
			SPabox = (*SPa)[si].box();
		    }
		    else {
			pSPa = &foo;
			SPabox = Box(IntVect::TheZeroVector(),IntVect::TheZeroVector());
		    }
		    getFaceMetric(r, reg, oitr(), geom);
		    hbflx3(BL_TO_FORTRAN(Flux[idim][si]),
			   BL_TO_FORTRAN_N(Soln[si], icomp),
			   ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
			   cdir, bctype, tfp, bho, bcl,
			   BL_TO_FORTRAN_N(fs, bdcomp),
			   BL_TO_FORTRAN(msk),
			   BL_TO_FORTRAN((*bcoefs[idim])[si]),
			   beta, dx, flux_factor, r.dataPtr(), inhom,
			   pSPa, ARLIM(SPabox.loVect()), ARLIM(SPabox.hiVect()));
		}
		else {
		    hbflx(BL_TO_FORTRAN(Flux[idim][si]),
			  BL_TO_FORTRAN_N(Soln[si], icomp),
			  ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
			  cdir, bct, bho, bcl,
			  BL_TO_FORTRAN_N(fs, bdcomp),
			  BL_TO_FORTRAN(msk),
			  BL_TO_FORTRAN((*bcoefs[idim])[si]),
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

void HypreABec::setupSolver(Real _reltol, Real _abstol, int maxiter)
{
  BL_PROFILE("HypreABec::setupSolver");

  const BoxArray& grids = acoefs->boxArray();

  const int size = BL_SPACEDIM + 1;
  int i, idim;

  int stencil_indices[size];

  for (i = 0; i < size; i++) {
    stencil_indices[i] = i;
  }

  Real foo=1.e200;

  FArrayBox matfab;
  for (MFIter ai(*acoefs); ai.isValid(); ++ai) {
    i = ai.index();
    const Box &reg = grids[i];

    matfab.resize(reg,size);
    Elixir matfab_elix = matfab.elixir();
    Real* mat = matfab.dataPtr();

    // build matrix interior

    // Note that we are using AoS indexing of matfab inside these functions.

#pragma gpu box(reg) sync
    hacoef(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
           BL_TO_FORTRAN_ANYD(matfab), 
           BL_TO_FORTRAN_ANYD((*acoefs)[ai]),
           alpha);

    for (idim = 0; idim < BL_SPACEDIM; ++idim) {
#pragma gpu box(reg) sync
        hbcoef(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
               BL_TO_FORTRAN_ANYD(matfab),
               BL_TO_FORTRAN_ANYD((*bcoefs[idim])[ai]),
               beta, AMREX_REAL_ANYD(dx), idim);
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
        const int *tfp = NULL;
        int bctype = bct;
        if (bd.mixedBndry(oitr())) {
          const BaseFab<int> &tf = *(bd.bndryTypes(oitr())[i]);
          tfp = tf.dataPtr();
          bctype = -1;
        }
        const Box &fsb = bd.bndryValues(oitr())[ai].box();
	Real* pSPa;
	Box SPabox; 
	if (SPa != 0) {
	  pSPa = (*SPa)[ai].dataPtr();
	  SPabox = (*SPa)[ai].box();
	}
	else {
	  pSPa = &foo;
	  SPabox = Box(IntVect::TheZeroVector(),IntVect::TheZeroVector());
	}

#pragma gpu box(reg) sync
        hbmat3(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
               reg.loVect()[0], reg.hiVect()[0],
               oitr().isLow(), idim+1,
               BL_TO_FORTRAN_ANYD(matfab),
               cdir, bctype,
               tfp, AMREX_INT_ANYD(fsb.loVect()), AMREX_INT_ANYD(fsb.hiVect()),
               bcl,
               msk.dataPtr(), AMREX_INT_ANYD(msk.loVect()), AMREX_INT_ANYD(msk.hiVect()),
               BL_TO_FORTRAN_ANYD((*bcoefs[idim])[ai]),
               beta, AMREX_REAL_ANYD(dx), flux_factor,
               pSPa, AMREX_INT_ANYD(SPabox.loVect()), AMREX_INT_ANYD(SPabox.hiVect()));
      }
      else {
#pragma gpu box(reg) sync
        hbmat(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
              BL_TO_FORTRAN_ANYD(matfab),
              cdir, bct, bcl,
              msk.dataPtr(), AMREX_INT_ANYD(msk.loVect()), AMREX_INT_ANYD(msk.hiVect()),
              BL_TO_FORTRAN_ANYD((*bcoefs[idim])[ai]),
              beta, AMREX_REAL_ANYD(dx));
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

void HypreABec::solve(MultiFab& dest, int icomp, MultiFab& rhs, BC_Mode inhom)
{
  BL_PROFILE("HypreABec::solve");

  const BoxArray& grids = dest.boxArray();

  int i, idim;

  //dest.setVal(0.0);

  Vector<Real> r;

  Real *vec;
  FArrayBox fnew;
  for (MFIter di(dest); di.isValid(); ++di) {
    i = di.index();
    const Box &reg = grids[i];

    // initialize dest, since we will reuse the space to set up rhs below:

    FArrayBox *f;
    int fcomp;
    if (dest.nGrow() == 0) { // need a temporary if dest is the wrong size
      f = &dest[di];
      fcomp = icomp;
    }
    else {
      f = &fnew;
      f->resize(reg);
      f->copy(dest[di], icomp, 0, 1);
      fcomp = 0;
    }
    Elixir f_elix = fnew.elixir();

    vec = f->dataPtr(fcomp); // sharing space, dest will be overwritten below

    HYPRE_StructVectorSetBoxValues(x, loV(reg), hiV(reg), vec);
    Gpu::synchronize();

    f->copy(rhs[di], 0, fcomp, 1);

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
          const int *tfp = NULL;
          int bctype = bct;
          if (bd.mixedBndry(oitr())) {
            const BaseFab<int> &tf = *(bd.bndryTypes(oitr())[i]);
            tfp = tf.dataPtr();
            bctype = -1;
          }
          getFaceMetric(r, reg, oitr(), geom);
	  hbvec3(vec, ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
		 cdir, bctype, tfp, bho, bcl,
		 BL_TO_FORTRAN_N(fs, bdcomp),
		 BL_TO_FORTRAN(msk),
		 BL_TO_FORTRAN((*bcoefs[idim])[di]),
		 beta, dx, r.dataPtr());
	}
	else {
	  hbvec(vec, ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
		cdir, bct, bho, bcl,
		BL_TO_FORTRAN_N(fs, bdcomp),
		BL_TO_FORTRAN(msk),
		BL_TO_FORTRAN((*bcoefs[idim])[di]),
		beta, dx);
	}
      }
    }

    // initialize rhs

    HYPRE_StructVectorSetBoxValues(b, loV(reg), hiV(reg), vec);
    Gpu::synchronize();
  }

  HYPRE_StructVectorAssemble(b); // currently a no-op
  HYPRE_StructVectorAssemble(x); // currently a no-op
  Gpu::synchronize();

  if (abstol > 0.0) {
    Real bnorm;
    bnorm = hypre_StructInnerProd((hypre_StructVector *) b,
				  (hypre_StructVector *) b);
    bnorm = sqrt(bnorm);

    const BoxArray& grids = acoefs->boxArray();
    Real volume = 0.0;
    for (int i = 0; i < grids.size(); i++) {
      volume += grids[i].numPts();
    }

    Real reltol_new = (bnorm > 0.0
		       ? abstol / bnorm * sqrt(volume)
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
      dest[di].copy(*f, 0, icomp, 1);
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
  bnorm = sqrt(bnorm);

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

  return bnorm * res / sqrt(volume);
}
