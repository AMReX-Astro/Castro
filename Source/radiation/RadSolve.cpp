
#include <AMReX_ParmParse.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_LO_BCTYPES.H>

#include <RadSolve.H>
#include <Radiation.H>  // for access to static physical constants only
#include <rad_util.H>
#include <problem_rad_source.H>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

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

#include <radsolve_queries.H>

    if (Radiation::SolverType == Radiation::SGFLDSolver ||
        Radiation::SolverType == Radiation::MGFLDSolver) {
        radsolve::abstol = 0.0;
    }

    // Check for unsupported options.

    if (AMREX_SPACEDIM == 1) {
        if (radsolve::level_solver_flag == 1) {
            amrex::Error("radsolve.level_solver_flag = 1 is not supported in 1D");
        }
    }

    if (Radiation::SolverType == Radiation::SGFLDSolver
        && Radiation::Er_Lorentz_term) {

        if (radsolve::level_solver_flag < 100) {
            amrex::Error("To do Lorentz term implicitly level_solver_flag must be >= 100.");
        }

        if (radsolve::use_hypre_nonsymmetric_terms == 0) {
            amrex::Error("To do Lorentz term implicitly use_hypre_nonsymmetric_terms must be 1.");
        }
    }

    if (Radiation::SolverType == Radiation::MGFLDSolver &&
        Radiation::accelerate == 2 && Radiation::nGroups > 1) {

        if (radsolve::level_solver_flag < 100) {
            amrex::Error("When accelerate is 2, level_solver_flag must be >= 100.");
        }

        if (radsolve::use_hypre_nonsymmetric_terms == 0) {
            amrex::Error("When accelerate is 2, use_hypre_nonsymmetric_terms must be 1.");
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
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
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
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
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

void RadSolve::levelSPas(int level, Array<MultiFab, AMREX_SPACEDIM>& lambda, int igroup,
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
      for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
          if (lo_bc[idim] == AMREX_LO_SANCHEZ_POMRANING &&
              reg.smallEnd(idim) == domainBox.smallEnd(idim)) {
              nexttoboundary=true;
              break;
          }
          if (hi_bc[idim] == AMREX_LO_SANCHEZ_POMRANING &&
              reg.bigEnd(idim) == domainBox.bigEnd(idim)) {
              nexttoboundary=true;
              break;
          }
      }

      if (nexttoboundary) {
          auto spa_arr = spa[mfi].array();

          auto lmx = lambda[0][mfi].array();
#if AMREX_SPACEDIM >= 2
          auto lmy = lambda[1][mfi].array();
#endif
#if AMREX_SPACEDIM == 3
          auto lmz = lambda[2][mfi].array();
#endif

          amrex::ParallelFor(reg,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
          {
              Real lam;

              if (i == reg.loVect()[0] || i == reg.hiVect()[0] ||
                  j == reg.loVect()[1] || j == reg.hiVect()[1] ||
                  k == reg.hiVect()[2] || k == reg.hiVect()[2]) {
#if AMREX_SPACEDIM == 1
                  if (i == reg.loVect()[0]) {
                      lam = lmx(i,j,k,igroup);
                  }
                  else {
                      lam = lmx(i+1,j,k,igroup);
                  }
#elif AMREX_SPACEDIM == 2
                  lam = 0.25_rt * (lmx(i,j,k,igroup) + lmx(i+1,j  ,k,igroup) +
                                   lmy(i,j,k,igroup) + lmy(i  ,j+1,k,igroup));
#else
                  lam = (lmx(i,j,k,igroup) + lmx(i+1,j  ,k  ,igroup) +
                         lmy(i,j,k,igroup) + lmy(i  ,j+1,k  ,igroup) +
                         lmz(i,j,k,igroup) + lmz(i  ,j  ,k+1,igroup)) / 6.e0_rt;
#endif

                  spa_arr(i,j,k) = FLDalpha(lam);
              }
          });
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
                            Array<MultiFab, AMREX_SPACEDIM>& lambda,
                            MultiFab& kappa_r, int kcomp,
                            Real c, int lamcomp)
{
  BL_PROFILE("RadSolve::levelBCoeffs");
  BL_ASSERT(kappa_r.nGrow() == 1);

  auto geomdata = parent->Geom(level).data();
  auto dx = parent->Geom(level).CellSizeArray();

  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {

    MultiFab bcoefs(lambda[idim].boxArray(), lambda[idim].DistributionMap(), 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(lambda[idim], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto bcoefs_arr = bcoefs[mfi].array();
        auto lambda_arr = lambda[idim][mfi].array(lamcomp);
        auto kappa_r_arr = kappa_r[mfi].array(kcomp);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            if (idim == 0) {

                Real r, s;
                edge_center_metric(i, j, k, idim, geomdata, r, s);

                if (AMREX_SPACEDIM == 1) {
                    s = 1.e0_rt;
                }

                Real kap = kavg(kappa_r_arr(i-1,j,k), kappa_r_arr(i,j,k), dx[0], -1);
                bcoefs_arr(i,j,k) = r * s * c * lambda_arr(i,j,k) / kap;

            }
            else if (idim == 1) {

                Real r, s;
                edge_center_metric(i, j, k, idim, geomdata, r, s);

                if (AMREX_SPACEDIM == 1) {
                    s = 1.e0_rt;
                }

                Real kap = kavg(kappa_r_arr(i,j-1,k), kappa_r_arr(i,j,k), dx[1], -1);
                bcoefs_arr(i,j,k) = r * s * c * lambda_arr(i,j,k) / kap;

            }
            else {

                Real r, s;
                edge_center_metric(i, j, k, idim, geomdata, r, s);

                if (AMREX_SPACEDIM == 1) {
                    s = 1.e0_rt;
                }

                Real kap = kavg(kappa_r_arr(i,j,k-1), kappa_r_arr(i,j,k), dx[2], -1);
                bcoefs_arr(i,j,k) = r * s * c * lambda_arr(i,j,k) / kap;

            }
        });
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

void RadSolve::levelDCoeffs(int level, Array<MultiFab, AMREX_SPACEDIM>& lambda,
                            MultiFab& vel, MultiFab& dcf)
{
    BL_PROFILE("RadSolve::levelDCoeffs");
    const Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
    const DistributionMapping& dm = castro->DistributionMap();
    const Geometry& geom = parent->Geom(level);
    const auto dx = geom.CellSizeArray();
    const auto geomdata = geom.data();

    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {

        MultiFab dcoefs(castro->getEdgeBoxArray(idim), dm, 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(dcoefs, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            auto dcoefs_arr = dcoefs[mfi].array();
            auto lambda_arr = lambda[idim][mfi].array();
            auto vel_arr = vel[mfi].array();
            auto dcf_arr = dcf[mfi].array();

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {
                Real r, s;

                if (idim == 0) {

                    edge_center_metric(i, j, k, idim, geomdata, r, s);

                    if (vel_arr(i-1,j,k,0) + vel_arr(i,j,k,0) > 0.e0_rt) {
                        dcoefs_arr(i,j,k) = dcf_arr(i-1,j,k) * vel_arr(i-1,j,k,0) * lambda_arr(i,j,k);
                    }
                    else if (vel_arr(i-1,j,k,0) + vel_arr(i,j,k,0) < 0.e0_rt) {
                        dcoefs_arr(i,j,k) = dcf_arr(i,j,k) * vel_arr(i,j,k,0) * lambda_arr(i,j,k);
                    }
                    else {
                        dcoefs_arr(i,j,k) = 0.e0_rt;
                    }

                    dcoefs_arr(i,j,k) = dcoefs_arr(i,j,k) * r;

                }
                else if (idim == 1) {

                    edge_center_metric(i, j, k, idim, geomdata, r, s);

                    if (vel_arr(i,j-1,k,1) + vel_arr(i,j,k,1) > 0.e0_rt) {
                        dcoefs_arr(i,j,k) = dcf_arr(i,j-1,k) * vel_arr(i,j-1,k,1) * lambda_arr(i,j,k);
                    }
                    else if (vel_arr(i,j-1,k,1) + vel_arr(i,j,k,1) < 0.e0_rt) {
                        dcoefs_arr(i,j,k) = dcf_arr(i,j,k) * vel_arr(i,j,k,1) * lambda_arr(i,j,k);
                    }
                    else {
                        dcoefs_arr(i,j,k) = 0.e0_rt;
                    }

                    dcoefs_arr(i,j,k) = dcoefs_arr(i,j,k) * r;

                }
                else {

                    edge_center_metric(i, j, k, idim, geomdata, r, s);

                    if (vel_arr(i,j,k-1,2) + vel_arr(i,j,k,2) > 0.e0_rt) {
                        dcoefs_arr(i,j,k) = dcf_arr(i,j,k-1) * vel_arr(i,j,k-1,2) * lambda_arr(i,j,k);
                    }
                    else if (vel_arr(i,j,k-1,2) + vel_arr(i,j,k,2) < 0.e0_rt) {
                        dcoefs_arr(i,j,k) = dcf_arr(i,j,k) * vel_arr(i,j,k,2) * lambda_arr(i,j,k);
                    }
                    else {
                        dcoefs_arr(i,j,k) = 0.e0_rt;
                    }

                    dcoefs_arr(i,j,k) = dcoefs_arr(i,j,k) * r;

                }

            });
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
  auto geomdata = geom.data();

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

      auto rhs_arr = rhs[mfi].array();
      auto temp_arr = temp[mfi].array();
      auto fkp_arr = fkp[mfi].array();
      auto eta_arr = eta[mfi].array();
      auto etainv_arr = etainv[mfi].array();
      auto rhoem_arr = rhoem[mfi].array();
      auto rhoes_arr = rhoes[mfi].array();
      auto dflux_old_arr = dflux_old[mfi].array();
      auto Er_old_arr = Er_old[mfi].array(0);
      auto Edot_arr = Edot[mfi].array();

      const Real dtm = 1.0_rt / delta_t;

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
      {
          Real ek = fkp_arr(i,j,k) * eta_arr(i,j,k);
          Real bs = etainv_arr(i,j,k) * 4.e0_rt * sigma * fkp_arr(i,j,k) * std::pow(temp_arr(i,j,k), 4);
          Real es = eta_arr(i,j,k) * (rhoem_arr(i,j,k) - rhoes_arr(i,j,k));
          Real ekt = (1.0_rt - theta) * eta_arr(i,j,k);

          Real r, s;
          cell_center_metric(i, j, k, geomdata, r, s);

          if (AMREX_SPACEDIM == 1) {
              s = 1.0_rt;
          }

          rhs_arr(i,j,k) = (rhs_arr(i,j,k) + r * s *
                            (bs + dtm * (Er_old_arr(i,j,k) + es) +
                             ek * c * Edot_arr(i,j,k) -
                             ekt * dflux_old_arr(i,j,k))) / (1.0_rt - ekt);
      });
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

void RadSolve::levelFluxFaceToCenter(int level, const Array<MultiFab, AMREX_SPACEDIM>& Flux,
                                     MultiFab& flx, int iflx)
{
    int nflx = flx.nComp();

    const Geometry& geom = parent->Geom(level);
    auto geomdata = geom.data();

    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(flx,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            auto t = flx[mfi].array();
            auto f = Flux[idim][mfi].array();

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {
                int it = idim * NGROUPS + iflx;

                Real r_left, s_left;
                edge_center_metric(i, j, k, idim, geomdata, r_left, s_left);

                Real r_right, s_right;
                edge_center_metric(i+1, j, k, idim, geomdata, r_right, s_right);

                if (idim == 0) {
                    t(i,j,k,it) = (f(i,j,k) / (r_left + 1.e-50_rt) + f(i+1,j,k) / r_right) * 0.5_rt;
                }
                else if (idim == 1) {
                    t(i,j,k,it) = (f(i,j,k) / r_left + f(i,j+1,k) / r_left) * 0.5_rt;
                }
                else {
                    t(i,j,k,it) = (f(i,j,k) + f(i,j,k+1)) * 0.5_rt;
                }
            });
        }
    }
}

void RadSolve::levelFlux(int level,
                         Array<MultiFab, AMREX_SPACEDIM>& Flux,
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

  auto dx = parent->Geom(level).CellSizeArray();

  for (int n = 0; n < AMREX_SPACEDIM; n++) {

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

          auto Erborder_arr = Erborder[mfi].array();
          auto bcoef_arr = bcoef[mfi].array();
          auto Flux_arr = Flux[n][mfi].array();

          Real beta = radsolve::beta;

          amrex::ParallelFor(bx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
          {
              if (n == 0) {

                  const Real fac = -beta / dx[0];

                  Flux_arr(i,j,k) = bcoef_arr(i,j,k) * (Erborder_arr(i,j,k) - Erborder_arr(i-1,j,k)) * fac;

              }
              else if (n == 1) {

                  const Real fac = -beta / dx[1];

                  Flux_arr(i,j,k) = bcoef_arr(i,j,k) * (Erborder_arr(i,j,k) - Erborder_arr(i,j-1,k)) * fac;

              }
              else {

                  const Real fac = -beta / dx[2];

                  Flux_arr(i,j,k) = bcoef_arr(i,j,k) * (Erborder_arr(i,j,k) - Erborder_arr(i,j,k-1)) * fac;

              }
          });
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
                            const Array<MultiFab, AMREX_SPACEDIM>& Flux,
                            int igroup)
{
  BL_PROFILE("RadSolve::levelFluxReg");

  const Real* dx = parent->Geom(level).CellSize();

  const Real volume = AMREX_D_TERM(dx[0], * dx[1], * dx[2]);

  if (flux_in) {
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
      const Real scale = volume / dx[n];
      flux_in->CrseInit(Flux[n], n, 0, igroup, 1, scale);
    }
  }
  if (flux_out) {
    for (OrientationIter face; face; ++face) {
      Orientation ori = face();
      (*flux_out)[ori].setVal(0.0, igroup, 1);
    }
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
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
  const GeometryData& geomdata = geom.data();
  auto dx = parent->Geom(level).CellSizeArray();
  const Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));

  Array<MultiFab, AMREX_SPACEDIM> Dterm_face;
  for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
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
  for (int n = 0; n < AMREX_SPACEDIM; n++) {
      const MultiFab *dp;

      dp = &hem->d2Coefficients(level, n);
      MultiFab &dcoef = *(MultiFab*)dp;

      for (MFIter fi(dcoef,true); fi.isValid(); ++fi) {
          const Box& bx = fi.tilebox();

          auto Er = Erborder[fi].array();
          auto dc = dcoef[fi].array();
          auto dtf = Dterm_face[n][fi].array();

          amrex::ParallelFor(bx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
          {
              if (n == 0) {
                  dtf(i,j,k) = (Er(i,j,k) - Er(i-1,j,k)) / dx[0] * dc(i,j,k);
              }
              else if (n == 1) {
                  dtf(i,j,k) = (Er(i,j,k) - Er(i,j-1,k)) / dx[1] * dc(i,j,k);
              }
              else {
                  dtf(i,j,k) = (Er(i,j,k) - Er(i,j,k-1)) / dx[2] * dc(i,j,k);
              }
          });
      }
  }

  // Correct D terms at physical and coarse-fine boundaries.
  hem->boundaryDterm(level, &Dterm_face[0], Er, igroup);

  // Correct for metric terms (only has an effect in non-Cartesian geometries).
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(Dterm_face[dir], true); mfi.isValid(); ++mfi) {
          const Box& box = mfi.tilebox();
          Array4<Real> const d = Dterm_face[dir][mfi].array();

          amrex::ParallelFor(box,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              if (dir == 0 && (geomdata.Coord() == CoordSys::SPHERICAL || geomdata.Coord() == CoordSys::RZ)) {
                  Real r, s;
                  edge_center_metric(i, j, k, dir, geomdata, r, s);

                  d(i,j,k) = d(i,j,k) / (r + 1.0e-50_rt);
              }
              else if (dir == 1 && geomdata.Coord() == CoordSys::RZ) {
                  Real r, s;
                  cell_center_metric(i, j, k, geomdata, r, s);

                  d(i,j,k) = d(i,j,k) / r;
              }
          });
      }
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter fi(Dterm,true); fi.isValid(); ++fi) {
      const Box& bx = fi.tilebox();

      auto Dx = Dterm_face[0][fi].array();
#if AMREX_SPACEDIM >= 2
      auto Dy = Dterm_face[1][fi].array();
#endif
#if AMREX_SPACEDIM == 3
      auto Dz = Dterm_face[2][fi].array();
#endif

      auto D = Dterm[fi].array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
      {
#if AMREX_SPACEDIM == 1
          D(i,j,k) = (Dx(i,j,k) + Dx(i+1,j,k)) * 0.5_rt;
#elif AMREX_SPACEDIM == 2
          D(i,j,k) = (Dx(i,j,k) + Dx(i+1,j,k) +
                      Dy(i,j,k) + Dy(i,j+1,k)) * 0.25_rt;
#else
          D(i,j,k) = (Dx(i,j,k) + Dx(i+1,j,k) +
                      Dy(i,j,k) + Dy(i,j+1,k) +
                      Dz(i,j,k) + Dz(i,j,k+1)) * (1.0_rt / 6.0_rt);
#endif
      });
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

  auto geomdata = geom.data();
  auto dx = geom.CellSizeArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(lambda, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      auto bcoefs_arr = bcoefs[mfi].array();
      auto lambda_arr = lambda[mfi].array(lamcomp);
      auto kappa_r_arr = kappa_r[mfi].array(kcomp);

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
      {
          if (idim == 0) {

              Real r, s;
              edge_center_metric(i, j, k, idim, geomdata, r, s);

              if (AMREX_SPACEDIM == 1) {
                  s = 1.e0_rt;
              }

              Real kap = kavg(kappa_r_arr(i-1,j,k), kappa_r_arr(i,j,k), dx[0], -1);
              bcoefs_arr(i,j,k) = r * s * c * lambda_arr(i,j,k) / kap;

          }
          else if (idim == 1) {

              Real r, s;
              edge_center_metric(i, j, k, idim, geomdata, r, s);

              if (AMREX_SPACEDIM == 1) {
                  s = 1.e0_rt;
              }

              Real kap = kavg(kappa_r_arr(i,j-1,k), kappa_r_arr(i,j,k), dx[1], -1);
              bcoefs_arr(i,j,k) = r * s * c * lambda_arr(i,j,k) / kap;

          }
          else {

              Real r, s;
              edge_center_metric(i, j, k, idim, geomdata, r, s);

              if (AMREX_SPACEDIM == 1) {
                  s = 1.e0_rt;
              }

              Real kap = kavg(kappa_r_arr(i,j,k-1), kappa_r_arr(i,j,k), dx[2], -1);
              bcoefs_arr(i,j,k) = r * s * c * lambda_arr(i,j,k) / kap;

          }
      });
  }
}

void RadSolve::levelACoeffs(int level, MultiFab& kpp,
                            Real delta_t, Real c, int igroup, Real ptc_tau)
{
  BL_PROFILE("RadSolve::levelACoeffs (MGFLD)");
  const BoxArray& grids = parent->boxArray(level);
  const DistributionMapping& dmap = parent->DistributionMap(level);
  const auto geomdata = parent->Geom(level).data();

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

      auto acoefs_arr = acoefs[mfi].array();
      auto kpp_arr = kpp[mfi].array(igroup);

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
      {
          Real r, s;
          cell_center_metric(i, j, k, geomdata, r, s);

          acoefs_arr(i,j,k) = c * kpp_arr(i,j,k) + 1.e0_rt / dt_ptc;
          acoefs_arr(i,j,k) = r * s * acoefs_arr(i,j,k);
      });
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
  auto geomdata = parent->Geom(level).data();

  const Real dt1 = 1.0_rt / delta_t;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter ri(rhs, TilingIfNotGPU()); ri.isValid(); ++ri) {

      const Box& bx = ri.tilebox();

      auto rhs_arr = rhs[ri].array();
      auto jg_arr = jg[ri].array();
      auto mugT_arr = mugT[ri].array();
      auto coupT_arr = coupT[ri].array();
      auto etaT_arr = etaT[ri].array();
      auto Er_step_arr = Er_step[ri].array();
      auto rhoe_step_arr = rhoe_step[ri].array();
      auto Er_star_arr = Er_star[ri].array();
      auto rhoe_star_arr = rhoe_star[ri].array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
      {
          Real Hg = mugT_arr(i,j,k,igroup) * etaT_arr(i,j,k);

          rhs_arr(i,j,k) = C::c_light * (jg_arr(i,j,k,igroup) + Hg * coupT_arr(i,j,k))
                           + dt1 * (Er_step_arr(i,j,k,igroup) - Hg * (rhoe_star_arr(i,j,k) - rhoe_step_arr(i,j,k))
                                    + ptc_tau * Er_star_arr(i,j,k,igroup));

          Real r, s;
          cell_center_metric(i, j, k, geomdata, r, s);

          rhs_arr(i,j,k) *= r;

          problem_rad_source(i, j, k, rhs_arr, geomdata, time, delta_t, igroup);
      });
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

void RadSolve::getEdgeMetric(int idim, const Geometry& geom, const Box& edgebox,
                             Vector<Real>& r, Vector<Real>& s)
{
    const Box& reg = amrex::enclosedCells(edgebox);
    const int I = (AMREX_SPACEDIM >= 2) ? 1 : 0;
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

      if (idim == 0) {
          for (int i = edgebox.loVect()[0]; i <= edgebox.hiVect()[0]; ++i) {
              r[i] *= r[i];
          }
#if AMREX_SPACEDIM >= 2
          Real h2 = 0.5e0_rt * dx[1];
          Real d2 = 1.e0_rt / dx[1];
          for (int j = edgebox.loVect()[1]; j <= edgebox.hiVect()[1]; ++j) {
              s[j] = d2 * (std::cos(s[j] - h2) - std::cos(s[j] + h2));
          }
#endif
      }
      else {
          Real h1 = 0.5e0_rt * dx[0];
          Real d1 = 1.e0_rt / (3.e0_rt * dx[0]);
          for (int i = edgebox.loVect()[0]; i <= edgebox.hiVect()[0]; ++i) {
              r[i] = d1 * (std::pow(r[i] + h1, 3) - std::pow(r[i] - h1, 3));
          }
#if AMREX_SPACEDIM >= 2
          for (int j = edgebox.loVect()[1]; j <= edgebox.hiVect()[1]; ++j) {
              s[j] = std::sin(s[j]);
          }
#endif
      }
    }
}

