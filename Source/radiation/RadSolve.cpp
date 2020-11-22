
#include <AMReX_ParmParse.H>
#include <AMReX_AmrLevel.H>

#include <AMReX_LO_BCTYPES.H>
//#include <CompSolver.H>

#include <RadSolve.H>
#include <Radiation.H>  // for access to static physical constants only

#include <rad_util.H>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <RAD_F.H>

#include <HABEC_F.H>    // only for nonsymmetric flux; may be changed?

using namespace amrex;

RadSolve::RadSolve (Amr* Parent, int level, const BoxArray& grids, const DistributionMapping& dmap)
    : parent(Parent)
{
    read_params();

    if (radsolve::level_solver_flag < 100) {
        hd.reset(new HypreABec(grids, dmap, parent->Geom(level), radsolve::level_solver_flag));
    }
    else {
        if (radsolve::use_hypre_nonsymmetric_terms == 0) {
            hm.reset(new HypreMultiABec(level, level, radsolve::level_solver_flag));
            hm->addLevel(level, parent->Geom(level), grids, dmap,
                         IntVect::TheUnitVector());
            hm->buildMatrixStructure();
        }
        else {
            hem.reset(new HypreExtMultiABec(level, level, radsolve::level_solver_flag));
            cMulti  = hem->cMultiplier();
            d1Multi = hem->d1Multiplier();
            d2Multi = hem->d2Multiplier();
            hem->addLevel(level, parent->Geom(level), grids, dmap,
                          IntVect::TheUnitVector());
            hem->buildMatrixStructure();
        }
    }
}

void
RadSolve::read_params ()
{
    ParmParse pp("radsolve");

    // Override some defaults manually.

    if (BL_SPACEDIM == 1) {
        // pfmg will not work in 1D
        radsolve::level_solver_flag = 0;
    }

    if (Radiation::SolverType == Radiation::SGFLDSolver
        && Radiation::Er_Lorentz_term) { 
        radsolve::use_hypre_nonsymmetric_terms = 1;
    }

    if (Radiation::SolverType == Radiation::MGFLDSolver && 
        Radiation::accelerate == 2 && Radiation::nGroups > 1) {
        radsolve::use_hypre_nonsymmetric_terms = 1;
    }

    if (Radiation::SolverType == Radiation::SGFLDSolver ||
        Radiation::SolverType == Radiation::MGFLDSolver) {
        radsolve::abstol = 0.0;
    }

#include <radsolve_queries.H>

    // Check for unsupported options.

    if (BL_SPACEDIM == 1) {
        if (radsolve::level_solver_flag == 1) {
            amrex::Error("radsolve.level_solver_flag = 1 is not supported in 1D");
        }
    }

    if (Radiation::SolverType == Radiation::SGFLDSolver
        && Radiation::Er_Lorentz_term) { 

        if (radsolve::level_solver_flag < 100) {
            amrex::Error("To do Lorentz term implicitly level_solver_flag must be >= 100.");
        }
    }

    if (Radiation::SolverType == Radiation::MGFLDSolver && 
        Radiation::accelerate == 2 && Radiation::nGroups > 1) {

        if (radsolve::level_solver_flag < 100) {
            amrex::Error("When accelerate is 2, level_solver_flag must be >= 100.");
        }
    }

}

void RadSolve::levelInit(int level)
{
  BL_PROFILE("RadSolve::levelInit");
}

void RadSolve::levelBndry(RadBndry& bd)
{
  BL_PROFILE("RadSolve::levelBndry");

  if (hd) {
    hd->setBndry(bd);
  }
  else if (hm) {
    hm->setBndry(hm->crseLevel(), bd);
  }
  else if (hem) {
    hem->setBndry(hem->crseLevel(), bd);
  }
}

// update multigroup version
void RadSolve::levelBndry(MGRadBndry& mgbd, const int comp)
{
  BL_PROFILE("RadSolve::levelBndryMG (updated)");

  if (hd) {
    hd->setBndry(mgbd, comp);
  }
  else if (hm) {
    hm->setBndry(hm->crseLevel(), mgbd, comp);
  }
  else if (hem) {
    hem->setBndry(hem->crseLevel(), mgbd, comp);
  }
}

void RadSolve::cellCenteredApplyMetrics(int level, MultiFab& cc)
{
    BL_PROFILE("RadSolve::cellCenteredApplyMetrics");
    BL_ASSERT(cc.nGrow() == 0);

    auto geomdata = parent->Geom(level).data();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(cc, TilingIfNotGPU()); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.tilebox();

        auto cc_arr = cc[mfi].array();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            Real r, s;
            cell_center_metric(i, j, k, geomdata, r, s);

            cc_arr(i,j,k) *= r * s;
        });
    }
}

void RadSolve::setLevelACoeffs(int level, const MultiFab& acoefs)
{
    if (hd) {
        hd->aCoefficients(acoefs);
    }
    else if (hm) {
        hm->aCoefficients(level, acoefs);
    }
    else if (hem) {
        hem->aCoefficients(level, acoefs);
    }
}

void RadSolve::setLevelBCoeffs(int level, const MultiFab& bcoefs, int dir)
{
    if (hd) {
        hd->bCoefficients(bcoefs, dir);
    }
    else if (hm) {
        hm->bCoefficients(level, bcoefs, dir);
    }
    else if (hem) {
        hem->bCoefficients(level, bcoefs, dir);
    }
}

void RadSolve::setLevelCCoeffs(int level, const MultiFab& ccoefs, int dir)
{
    if (hem) {
      hem->cCoefficients(level, ccoefs, dir);
    }
}

void RadSolve::levelACoeffs(int level,
                            MultiFab& fkp, MultiFab& eta, MultiFab& etainv,
                            Real c, Real delta_t, Real theta)
{
  BL_PROFILE("RadSolve::levelACoeffs");
  const BoxArray& grids = parent->boxArray(level);
  const DistributionMapping& dmap = parent->DistributionMap(level);

  auto geomdata = parent->Geom(level).data();

  // Allocate space for ABecLapacian acoeffs, fill with values

  int Ncomp  = 1;
  int Nghost = 0;

  MultiFab acoefs(grids, dmap, Ncomp, Nghost);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(fkp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box &bx = mfi.tilebox();

      auto a = acoefs[mfi].array();
      auto fkp_arr = fkp[mfi].array();
      auto eta_arr = eta[mfi].array();
      auto etainv_arr = etainv[mfi].array();

      const Real dtm = 1.e0_rt / delta_t;

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {
          Real r, s;
          cell_center_metric(i, j, k, geomdata, r, s);

          if (AMREX_SPACEDIM == 1) {
              s = 1.e0_rt;
          }

          a(i,j,k) = r * s *
                     (fkp_arr(i,j,k) * etainv_arr(i,j,k) * c + dtm) /
                     (1.e0_rt - (1.e0_rt - theta) * eta_arr(i,j,k));
      });
  }

  if (hd) {
    hd->aCoefficients(acoefs);
  }
  else if (hm) {
    hm->aCoefficients(level, acoefs);
  }
  else if (hem) {
    hem->aCoefficients(level, acoefs);
  }
}

void RadSolve::levelSPas(int level, Array<MultiFab, BL_SPACEDIM>& lambda, int igroup, 
                         int lo_bc[3], int hi_bc[3])
{
  const BoxArray& grids = parent->boxArray(level);
  const DistributionMapping& dmap = parent->DistributionMap(level);
  const Geometry& geom = parent->Geom(level);
  const Box& domainBox = geom.Domain();

  MultiFab spa(grids, dmap, 1, 0);
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(spa,true); mfi.isValid(); ++mfi) {
      const Box& reg  = mfi.tilebox();
    
      spa[mfi].setVal<RunOn::Host>(1.e210,reg,0);
    
      bool nexttoboundary=false;
      for (int idim=0; idim<BL_SPACEDIM; idim++) {
          if (lo_bc[idim] == LO_SANCHEZ_POMRANING &&
              reg.smallEnd(idim) == domainBox.smallEnd(idim)) {
              nexttoboundary=true;
              break;
          }
          if (hi_bc[idim] == LO_SANCHEZ_POMRANING &&
              reg.bigEnd(idim) == domainBox.bigEnd(idim)) {
              nexttoboundary=true;
              break;
          }
      }
    
      if (nexttoboundary) {
          ca_spalpha(reg.loVect(), reg.hiVect(),
                     BL_TO_FORTRAN(spa[mfi]),
                     D_DECL(BL_TO_FORTRAN(lambda[0][mfi]),
                            BL_TO_FORTRAN(lambda[1][mfi]),
                            BL_TO_FORTRAN(lambda[2][mfi])),
                     &igroup);
      }
  }

  if (hm) {
    hm->SPalpha(level, spa);
  }
  else if (hem) {
    hem->SPalpha(level, spa);
  }
  else if (hd) {
    hd->SPalpha(spa);
  }
  else {
    amrex::Abort("Should not be in RadSolve::levelSPas");    
  }
}

void RadSolve::levelBCoeffs(int level,
                            Array<MultiFab, BL_SPACEDIM>& lambda,
                            MultiFab& kappa_r, int kcomp,
                            Real c, int lamcomp)
{
  BL_PROFILE("RadSolve::levelBCoeffs");
  BL_ASSERT(kappa_r.nGrow() == 1);

  const Geometry& geom = parent->Geom(level);
  const Real* dx       = geom.CellSize();

  for (int idim = 0; idim < BL_SPACEDIM; ++idim) {

    MultiFab bcoefs(lambda[idim].boxArray(), lambda[idim].DistributionMap(), 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(lambda[idim], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
        bclim(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
              BL_TO_FORTRAN_ANYD(bcoefs[mfi]), 
              BL_TO_FORTRAN_N_ANYD(lambda[idim][mfi], lamcomp),
              idim, 
              BL_TO_FORTRAN_N_ANYD(kappa_r[mfi], kcomp), 
              c, AMREX_REAL_ANYD(dx));
    }

    if (hd) {
        hd->bCoefficients(bcoefs, idim);
    }
    else if (hm) {
      hm->bCoefficients(level, bcoefs, idim);
    }
    else if (hem) {
      hem->bCoefficients(level, bcoefs, idim);
    }
  } // -->> over dimension
}

void RadSolve::levelDCoeffs(int level, Array<MultiFab, BL_SPACEDIM>& lambda,
                            MultiFab& vel, MultiFab& dcf)
{
    BL_PROFILE("RadSolve::levelDCoeffs");
    const Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
    const DistributionMapping& dm = castro->DistributionMap();
    const Geometry& geom = parent->Geom(level);
    const Real* dx       = geom.CellSize();

    for (int idim=0; idim<BL_SPACEDIM; idim++) {

        MultiFab dcoefs(castro->getEdgeBoxArray(idim), dm, 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(dcoefs, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
            ca_compute_dcoefs(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                              BL_TO_FORTRAN_ANYD(dcoefs[mfi]),
                              BL_TO_FORTRAN_ANYD(lambda[idim][mfi]),
                              BL_TO_FORTRAN_ANYD(vel[mfi]),
                              BL_TO_FORTRAN_ANYD(dcf[mfi]), 
                              AMREX_REAL_ANYD(dx), idim);

        }

        hem->d2Coefficients(level, dcoefs, idim);
        hem->d2Multiplier() = 1.0;
    }
}

void RadSolve::levelRhs(int level, MultiFab& rhs,
                        MultiFab& temp,
                        MultiFab& fkp, MultiFab& eta, MultiFab& etainv,
                        MultiFab& rhoem, MultiFab& rhoes,
                        MultiFab& dflux_old, MultiFab& Er_old, MultiFab& Edot,
                        Real delta_t, Real sigma, Real c, Real theta,
                        FluxRegister* fine_corr, Real scale,
                        int igroup, Real nu, Real dnu)
{
  BL_PROFILE("RadSolve::levelRhs");
  BL_ASSERT(rhs.nGrow() == 0);

  const Geometry& geom = parent->Geom(level);
  const Real* dx       = geom.CellSize();

  rhs.setVal(0.0);
  if (fine_corr) {
    // This works trivially for a multilevel solve since the finer level is
    // present in the rhs to overwrite the junk produced under it.
    // In a single-level version we have to be sure that fine_corr
    // has been cleaned up using ClearInternalBorders.

    // Hack:  For the single group case igroup defaults to -1, which is
    // significant later in this routine.  So we have to construct the
    // correct component number here:

    int igrouptmp = (igroup < 0) ? 0 : igroup;

    fine_corr->Reflux(rhs, scale, igrouptmp, 0, 1, parent->Geom(level));
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
      lrhs(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           BL_TO_FORTRAN_ANYD(rhs[mfi]), 
           BL_TO_FORTRAN_ANYD(temp[mfi]),
           BL_TO_FORTRAN_ANYD(fkp[mfi]),
           BL_TO_FORTRAN_ANYD(eta[mfi]),
           BL_TO_FORTRAN_ANYD(etainv[mfi]),
           BL_TO_FORTRAN_ANYD(rhoem[mfi]),
           BL_TO_FORTRAN_ANYD(rhoes[mfi]),
           BL_TO_FORTRAN_ANYD(dflux_old[mfi]),
           BL_TO_FORTRAN_N_ANYD(Er_old[mfi], 0), 
           BL_TO_FORTRAN_ANYD(Edot[mfi]),
           delta_t, AMREX_REAL_ANYD(dx), sigma, c, theta);
  }
}


void RadSolve::levelSolve(int level,
                          MultiFab& Er, int igroup, MultiFab& rhs,
                          Real sync_absres_factor)
{
  BL_PROFILE("RadSolve::levelSolve");

  // Set coeffs, build solver, solve
  if (hd) {
    hd->setScalars(radsolve::alpha, radsolve::beta);
  }
  else if (hm) {
    hm->setScalars(radsolve::alpha, radsolve::beta);
  }
  else if (hem) {
    hem->setScalars(radsolve::alpha, radsolve::beta);
  }

  if (hd) {
    hd->setupSolver(radsolve::reltol, radsolve::abstol, radsolve::maxiter);
    hd->solve(Er, igroup, rhs, Inhomogeneous_BC);
    Real res = hd->getAbsoluteResidual();
    if (verbose >= 2 && ParallelDescriptor::IOProcessor()) {
      int oldprec = std::cout.precision(20);
      std::cout << "Absolute residual = " << res << std::endl;
      std::cout.precision(oldprec);
    }
    res *= sync_absres_factor;
    hd->clearSolver();
  }
  else if (hm) {
    hm->loadMatrix();
    hm->finalizeMatrix();
    hm->loadLevelVectors(level, Er, igroup, rhs, Inhomogeneous_BC);
    hm->finalizeVectors();
    hm->setupSolver(radsolve::reltol, radsolve::abstol, radsolve::maxiter);
    hm->solve();
    hm->getSolution(level, Er, igroup);
    Real res = hm->getAbsoluteResidual();
    if (verbose >= 2 && ParallelDescriptor::IOProcessor()) {
      int oldprec = std::cout.precision(20);
      std::cout << "Absolute residual = " << res << std::endl;
      std::cout.precision(oldprec);
    }
    res *= sync_absres_factor;
    hm->clearSolver();
  }
  else if (hem) {
    hem->loadMatrix();
    hem->finalizeMatrix();
    hem->loadLevelVectors(level, Er, igroup, rhs, Inhomogeneous_BC);
    hem->finalizeVectors();
    hem->setupSolver(radsolve::reltol, radsolve::abstol, radsolve::maxiter);
    hem->solve();
    hem->getSolution(level, Er, igroup);
    Real res = hem->getAbsoluteResidual();
    if (verbose >= 2 && ParallelDescriptor::IOProcessor()) {
      int oldprec = std::cout.precision(20);
      std::cout << "Absolute residual = " << res << std::endl;
      std::cout.precision(oldprec);
    }
    res *= sync_absres_factor;
    hem->clearSolver();
  }
}

void RadSolve::levelFluxFaceToCenter(int level, const Array<MultiFab, BL_SPACEDIM>& Flux,
                                     MultiFab& flx, int iflx)
{
    int nflx = flx.nComp();
    
    const Geometry& geom = parent->Geom(level);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<Real> r, s;
    
        for (int idim = 0; idim < BL_SPACEDIM; idim++) {
            for (MFIter mfi(flx,true); mfi.isValid(); ++mfi) 
            {
                const Box &ccbx  = mfi.tilebox();
                const Box &ndbx = amrex::surroundingNodes(ccbx, idim);

                getEdgeMetric(idim, geom, ndbx, r, s);

                int rlo = ndbx.smallEnd(0);
                int rhi = rlo + r.size() - 1;

                ca_flux_face2center(ccbx.loVect(), ccbx.hiVect(),
                                    BL_TO_FORTRAN(flx[mfi]),
                                    BL_TO_FORTRAN(Flux[idim][mfi]),
                                    r.dataPtr(), &rlo, &rhi, 
                                    &nflx, &idim, &iflx);
            }
        }
    }
}

void RadSolve::levelFlux(int level,
                         Array<MultiFab, BL_SPACEDIM>& Flux,
                         MultiFab& Er, int igroup)
{
  BL_PROFILE("RadSolve::levelFlux");
  const BoxArray& grids = parent->boxArray(level);
  const DistributionMapping& dmap = parent->DistributionMap(level);

  // grow a larger MultiFab to hold Er so we can difference across faces
  MultiFab Erborder(grids, dmap, 1, 1);
  Erborder.setVal(0.0);
  MultiFab::Copy(Erborder, Er, igroup, 0, 1, 0);

  Erborder.FillBoundary(parent->Geom(level).periodicity()); // zeroes left in off-level boundaries

  const Real* dx = parent->Geom(level).CellSize();

  for (int n = 0; n < BL_SPACEDIM; n++) {

      const MultiFab *bp;

      if (hd) {
          bp = &hd->bCoefficients(n);
      }
      else if (hm) {
          bp = &hm->bCoefficients(level, n);
      }
      else if (hem) {
          bp = &hem->bCoefficients(level, n);
      }

      MultiFab &bcoef = *(MultiFab*)bp;

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(Flux[n], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
          set_abec_flux(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                        n,
                        BL_TO_FORTRAN_ANYD(Erborder[mfi]), 
                        BL_TO_FORTRAN_ANYD(bcoef[mfi]), 
                        radsolve::beta, AMREX_REAL_ANYD(dx),
                        BL_TO_FORTRAN_ANYD(Flux[n][mfi]));
      }

  }

  // Correct fluxes at physical and coarse-fine boundaries.

  // Note: It would be good to move much of this section into the
  // linear solver drivers themselves (one layer down), to enable
  // an implementation that does not compute fluxes everywhere in
  // an EdgeVar.  We can't just move the nonsymmetric pieces easily
  // by themselves, though, because the current implementation
  // trashes the boundary fluxes before fixing them.

  if (hd) {
    hd->boundaryFlux(&Flux[0], Er, igroup, Inhomogeneous_BC);
  }
  else if (hm) {
    hm->boundaryFlux(level, &Flux[0], Er, igroup, Inhomogeneous_BC);
  }
}

void RadSolve::levelFluxReg(int level,
                            FluxRegister* flux_in, FluxRegister* flux_out,
                            const Array<MultiFab, BL_SPACEDIM>& Flux,
                            int igroup)
{
  BL_PROFILE("RadSolve::levelFluxReg");

  const Real* dx = parent->Geom(level).CellSize();

  const Real volume = D_TERM(dx[0], * dx[1], * dx[2]);

  if (flux_in) {
    for (int n = 0; n < BL_SPACEDIM; n++) {
      const Real scale = volume / dx[n];
      flux_in->CrseInit(Flux[n], n, 0, igroup, 1, scale);
    }
  }
  if (flux_out) {
    for (OrientationIter face; face; ++face) {
      Orientation ori = face();
      (*flux_out)[ori].setVal(0.0, igroup, 1);
    }
    for (int n = 0; n < BL_SPACEDIM; n++) {
      const Real scale = volume / dx[n];
      flux_out->FineAdd(Flux[n], n, 0, igroup, 1, scale);
    }
  }
}

void RadSolve::levelDterm(int level, MultiFab& Dterm, MultiFab& Er, int igroup)
{
  BL_PROFILE("RadSolve::levelDterm");
  const BoxArray& grids = parent->boxArray(level);
  const DistributionMapping& dmap = parent->DistributionMap(level);
  const Geometry& geom = parent->Geom(level);
  const Real* dx = parent->Geom(level).CellSize();
  const Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));

  Array<MultiFab, BL_SPACEDIM> Dterm_face;
  for (int idim=0; idim<BL_SPACEDIM; idim++) {
      Dterm_face[idim].define(castro->getEdgeBoxArray(idim), dmap, 1, 0);
  }

  // grow a larger MultiFab to hold Er so we can difference across faces
  MultiFab Erborder(grids, dmap, 1, 1);
  Erborder.setVal(0.0);
  MultiFab::Copy(Erborder, Er, igroup, 0, 1, 0);

  Erborder.FillBoundary(parent->Geom(level).periodicity()); // zeroes left in off-level boundaries

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (int n = 0; n < BL_SPACEDIM; n++) {
      const MultiFab *dp;

      dp = &hem->d2Coefficients(level, n);
      MultiFab &dcoef = *(MultiFab*)dp;
      
      for (MFIter fi(dcoef,true); fi.isValid(); ++fi) {
          const Box& bx = fi.tilebox();
          ca_set_dterm_face(bx.loVect(), bx.hiVect(),
                            BL_TO_FORTRAN(Erborder[fi]),
                            BL_TO_FORTRAN(dcoef[fi]), 
                            BL_TO_FORTRAN(Dterm_face[n][fi]), 
                            dx, &n);
      }
  }

  // Correct D terms at physical and coarse-fine boundaries.
  hem->boundaryDterm(level, &Dterm_face[0], Er, igroup);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      Vector<Real> rc, re, s;
      
      if (geom.IsSPHERICAL()) {
          for (MFIter fi(Dterm_face[0]); fi.isValid(); ++fi) {  // omp over boxes
              int i = fi.index();
              const Box &reg = grids[i];
              parent->Geom(level).GetEdgeLoc(re, reg, 0);
              parent->Geom(level).GetCellLoc(rc, reg, 0);
              parent->Geom(level).GetCellLoc(s, reg, 0);
              const Box &dbox = Dterm_face[0][fi].box();
              sphe(re.dataPtr(), s.dataPtr(), 0,
                   ARLIM(dbox.loVect()), ARLIM(dbox.hiVect()), dx);
              
              ca_correct_dterm(D_DECL(BL_TO_FORTRAN(Dterm_face[0][fi]),
                                      BL_TO_FORTRAN(Dterm_face[1][fi]),
                                      BL_TO_FORTRAN(Dterm_face[2][fi])),
                               re.dataPtr(), rc.dataPtr());
          }
#ifdef _OPENMP
#pragma omp barrier
#endif
      }
      else if (geom.IsRZ()) {
          for (MFIter fi(Dterm_face[0]); fi.isValid(); ++fi) {  // omp over boxes
              int i = fi.index();
              const Box &reg = grids[i];
              parent->Geom(level).GetEdgeLoc(re, reg, 0);
              parent->Geom(level).GetCellLoc(rc, reg, 0);
              
              ca_correct_dterm(D_DECL(BL_TO_FORTRAN(Dterm_face[0][fi]),
                                      BL_TO_FORTRAN(Dterm_face[1][fi]),
                                      BL_TO_FORTRAN(Dterm_face[2][fi])),
                               re.dataPtr(), rc.dataPtr());
          }
#ifdef _OPENMP
#pragma omp barrier
#endif
      }

      for (MFIter fi(Dterm,true); fi.isValid(); ++fi) {
          const Box& bx = fi.tilebox();
          int scomp = 0;
          int dcomp = 0;
          int ncomp = 1;
          int nf = 1;
          int nc = 1;
          ca_face2center(bx.loVect(), bx.hiVect(), scomp, dcomp, ncomp, nf, nc,
                         D_DECL(BL_TO_FORTRAN(Dterm_face[0][fi]),
                                BL_TO_FORTRAN(Dterm_face[1][fi]),
                                BL_TO_FORTRAN(Dterm_face[2][fi])),
                         BL_TO_FORTRAN(Dterm[fi]));
      }
  }
}

// <MGFLD routines>
void RadSolve::computeBCoeffs(MultiFab& bcoefs, int idim,
                              MultiFab& kappa_r, int kcomp,
                              MultiFab& lambda, int lamcomp,
                              Real c, const Geometry& geom)
{
  BL_PROFILE("RadSolve::computeBCoeffs (MGFLD)");
  BL_ASSERT(kappa_r.nGrow() == 1);

  const Real* dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(lambda, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
      bclim(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
            BL_TO_FORTRAN_ANYD(bcoefs[mfi]),
            BL_TO_FORTRAN_N_ANYD(lambda[mfi], lamcomp),
            idim, 
            BL_TO_FORTRAN_N_ANYD(kappa_r[mfi], kcomp), 
            c, AMREX_REAL_ANYD(dx));
  }
}

void RadSolve::levelACoeffs(int level, MultiFab& kpp, 
                            Real delta_t, Real c, int igroup, Real ptc_tau)
{
  BL_PROFILE("RadSolve::levelACoeffs (MGFLD)");
  const BoxArray& grids = parent->boxArray(level);
  const DistributionMapping& dmap = parent->DistributionMap(level);
  const Real* dx = parent->Geom(level).CellSize();

  // allocate space for ABecLaplacian acoeffs, fill with values

  int Ncomp = 1;
  int Nghost = 0;
  MultiFab acoefs(grids, dmap, Ncomp, Nghost);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(kpp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();
          
      Real dt_ptc = delta_t / (1.0 + ptc_tau);

#pragma gpu box(bx)
      lacoefmgfld(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                  BL_TO_FORTRAN_ANYD(acoefs[mfi]), 
                  BL_TO_FORTRAN_N_ANYD(kpp[mfi], igroup), 
                  AMREX_REAL_ANYD(dx), dt_ptc, c);
  }

  // set a coefficients
  if (hd) {
    hd->aCoefficients(acoefs);
  }
  else if (hm) {
    hm->aCoefficients(level,acoefs);
  }
  else if (hem) {
    hem->aCoefficients(level,acoefs);
  }
}


void RadSolve::levelRhs(int level, MultiFab& rhs, const MultiFab& jg, 
                        const MultiFab& mugT,
                        const MultiFab& coupT,
                        const MultiFab& etaT,
                        const MultiFab& Er_step, const MultiFab& rhoe_step,
                        const MultiFab& Er_star, const MultiFab& rhoe_star,
                        Real delta_t, int igroup, int it, Real ptc_tau)
{
  BL_PROFILE("RadSolve::levelRhs (MGFLD version)");
  Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
  Real time = castro->get_state_data(Rad_Type).curTime();
  const Real* dx = parent->Geom(level).CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter ri(rhs, TilingIfNotGPU()); ri.isValid(); ++ri) {

      const Box& bx = ri.tilebox();

#pragma gpu box(bx)
      ca_compute_rhs(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                     BL_TO_FORTRAN_ANYD(rhs[ri]),
                     BL_TO_FORTRAN_ANYD(jg[ri]),
                     BL_TO_FORTRAN_ANYD(mugT[ri]),
                     BL_TO_FORTRAN_ANYD(coupT[ri]),
                     BL_TO_FORTRAN_ANYD(etaT[ri]),
                     BL_TO_FORTRAN_ANYD(Er_step[ri]),
                     BL_TO_FORTRAN_ANYD(rhoe_step[ri]),
                     BL_TO_FORTRAN_ANYD(Er_star[ri]),
                     BL_TO_FORTRAN_ANYD(rhoe_star[ri]),
                     AMREX_REAL_ANYD(dx), delta_t, igroup, ptc_tau);

#pragma gpu box(bx)
      ca_rad_source(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                    BL_TO_FORTRAN_ANYD(rhs[ri]),
                    AMREX_REAL_ANYD(dx), delta_t, time, igroup);

  }
}

// </ MGFLD routines>

void RadSolve::setHypreMulti(Real cMul, Real d1Mul, Real d2Mul)
{
  if (hem) {
    hem-> cMultiplier() =  cMul;
    hem->d1Multiplier() = d1Mul;
    hem->d2Multiplier() = d2Mul;
  }
}

void RadSolve::restoreHypreMulti()
{
  if (hem) {
    hem-> cMultiplier() =  cMulti;
    hem->d1Multiplier() = d1Multi;
    hem->d2Multiplier() = d2Multi;  
  }
}

void RadSolve::getCellCenterMetric(const Geometry& geom, const Box& reg, Vector<Real>& r, Vector<Real>& s)
{
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (geom.IsCartesian()) {
        r.resize(reg.length(0), 1);
        s.resize(reg.length(I), 1);
    }
    else if (geom.IsRZ()) {
        geom.GetCellLoc(r, reg, 0);
        s.resize(reg.length(I), 1);
    }
    else {
        geom.GetCellLoc(r, reg, 0);
        geom.GetCellLoc(s, reg, I);
        const Real *dx = geom.CellSize();
        sphc(r.dataPtr(), s.dataPtr(),
             ARLIM(reg.loVect()), ARLIM(reg.hiVect()), dx);
    }
}
        
void RadSolve::getEdgeMetric(int idim, const Geometry& geom, const Box& edgebox, 
                             Vector<Real>& r, Vector<Real>& s)
{
    const Box& reg = amrex::enclosedCells(edgebox);
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (geom.IsCartesian()) {
        r.resize(reg.length(0)+1, 1);
        s.resize(reg.length(I)+1, 1);
    }
    else if (geom.IsRZ()) {
        if (idim == 0) {
            geom.GetEdgeLoc(r, reg, 0);
        }
        else {
            geom.GetCellLoc(r, reg, 0);
        }
        s.resize(reg.length(I)+1, 1);
    }
    else {
      if (idim == 0) {
        geom.GetEdgeLoc(r, reg, 0);
        geom.GetCellLoc(s, reg, I);
      }
      else {
        geom.GetCellLoc(r, reg, 0);
        geom.GetEdgeLoc(s, reg, I);
      }
      const Real *dx = geom.CellSize();
      sphe(r.dataPtr(), s.dataPtr(), idim,
           ARLIM(edgebox.loVect()), ARLIM(edgebox.hiVect()), dx);
    }
}

