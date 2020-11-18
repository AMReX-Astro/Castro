#include <Radiation.H>
#include <Castro_F.H>
#include <RAD_F.H>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

void Radiation::check_convergence_er(Real& relative, Real& absolute, Real& err_er,
                                     const MultiFab& Er_new, const MultiFab& Er_pi, 
                                     const MultiFab& kappa_p,
                                     const MultiFab& etaTz,
                                     const MultiFab& temp_new,
                                     Real delta_t)
{
  BL_PROFILE("Radiation::check_convergence_er (MGFLD)");

  ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpMax> reduce_op;
  ReduceData<Real, Real, Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(Er_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      const auto Ern = Er_new[mfi].array();
      const auto Erl = Er_pi[mfi].array();
      const auto kap = kappa_p[mfi].array();
      const auto etTz = etaTz[mfi].array();
      const auto temp = temp_new[mfi].array();

      reduce_op.eval(bx, reduce_data,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
      {
          Real chg = 0.e0_rt;
          Real tot = 0.e0_rt;
          Real kde = 0.e0_rt;

          for (int g = 0; g < NGROUPS; ++g) {
              Real der = Ern(i,j,k,g) - Erl(i,j,k,g);
              chg = chg + std::abs(der);
              tot = tot + std::abs(Ern(i,j,k,g));
              kde = kde + kap(i,j,k,g) * der;
          }

          Real abso = chg;
          Real rela = chg / (tot + 1.e-50_rt);

          Real err_T = etTz(i,j,k) * kde;
          Real errr = std::abs(err_T / (temp(i,j,k) + 1.e-50_rt));

          return {rela, abso, errr};
      });
  }

  ReduceTuple hv = reduce_data.value();

  relative = amrex::get<0>(hv);
  absolute = amrex::get<1>(hv);
  err_er   = amrex::get<2>(hv);

  Real data[3] = {relative, absolute, err_er};

  ParallelDescriptor::ReduceRealMax(data, 3);

  relative = data[0];
  absolute = data[1];
  err_er   = data[2];
}


void Radiation::check_convergence_matt(const MultiFab& rhoe_new, const MultiFab& rhoe_star, 
                                       const MultiFab& rhoe_step, const MultiFab& Er_new, 
                                       const MultiFab& temp_new, const MultiFab& temp_star, 
                                       const MultiFab& rho,
                                       const MultiFab& kappa_p, const MultiFab& jg,
                                       const MultiFab& dedT,
                                       Real& rel_rhoe, Real& abs_rhoe, 
                                       Real& rel_FT,   Real& abs_FT, 
                                       Real& rel_T,    Real& abs_T, 
                                       Real delta_t)
{
  ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpMax, ReduceOpMax, ReduceOpMax, ReduceOpMax> reduce_op;
  ReduceData<Real, Real, Real, Real, Real, Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(rhoe_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      auto ren = rhoe_new[mfi].array();
      auto res = rhoe_star[mfi].array();
      auto re2 = rhoe_step[mfi].array();
      auto Ern = Er_new[mfi].array();
      auto Tmn = temp_new[mfi].array();
      auto Tms = temp_star[mfi].array();
      auto rho_arr = rho[mfi].array();
      auto kap = kappa_p[mfi].array();
      auto jg_arr = jg[mfi].array();
      auto deT = dedT[mfi].array();

      int ng = Radiation::nGroups;
      Real cdt = C::c_light * delta_t;

      reduce_op.eval(bx, reduce_data,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
      {
          Real abschg_re = std::abs(ren(i,j,k) - res(i,j,k));
          Real relchg_re = std::abs(abschg_re / (ren(i,j,k) + 1.e-50_rt));

          Real abschg_T = std::abs(Tmn(i,j,k) - Tms(i,j,k));
          Real relchg_T = std::abs(abschg_T / (Tmn(i,j,k) + 1.e-50_rt));

          Real FT = ren(i,j,k) - re2(i,j,k);
          for (int g = 0; g < ng; ++g) {
              FT -= cdt * (kap(i,j,k,g) * Ern(i,j,k,g) - jg_arr(i,j,k,g));
          }
          FT = std::abs(FT);

          Real dTe = Tmn(i,j,k);
          Real FTdenom = rho_arr(i,j,k) * std::abs(deT(i,j,k) * dTe);
          //Real FTdenom = amrex::max(std::abs(ren(i,j,k) - re2(i,j,k)), std::abs(ren(i,j,k) * 1.e-15_rt))

          Real abschg_FT = FT;
          Real relchg_FT = FT / (FTdenom + 1.e-50_rt);

          return {relchg_re, abschg_re, relchg_FT, abschg_FT, relchg_T, abschg_T};
      });
  }

  ReduceTuple hv = reduce_data.value();
  rel_rhoe = amrex::get<0>(hv);
  abs_rhoe = amrex::get<1>(hv);
  rel_FT   = amrex::get<2>(hv);
  abs_FT   = amrex::get<3>(hv);
  rel_T    = amrex::get<4>(hv);
  abs_T    = amrex::get<5>(hv);

  int ndata = 6;
  Real data[6] = {rel_rhoe, abs_rhoe, rel_FT, abs_FT, rel_T, abs_T};

  ParallelDescriptor::ReduceRealMax(data, ndata);

  rel_rhoe = data[0]; 
  abs_rhoe = data[1];
  rel_FT   = data[2];
  abs_FT   = data[3];
  rel_T    = data[4];
  abs_T    = data[5];
}


void Radiation::compute_coupling(MultiFab& coupT,
                                 const MultiFab& kpp, const MultiFab& Eg,
                                 const MultiFab& jg)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(kpp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto coupT_arr = coupT[mfi].array();
        auto kpp_arr = kpp[mfi].array();
        auto Eg_arr = Eg[mfi].array();
        auto jg_arr = jg[mfi].array();

        int ng = Radiation::nGroups;

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            coupT_arr(i,j,k) = 0.0_rt;

            for (int g = 0; g < ng; ++g) {
                coupT_arr(i,j,k) = coupT_arr(i,j,k) + (kpp_arr(i,j,k,g) * Eg_arr(i,j,k,g) - jg_arr(i,j,k,g));
            }
        });
    }
}


void Radiation::compute_etat(MultiFab& etaT, MultiFab& etaTz, 
                             MultiFab& eta1, MultiFab& djdT,
                             const MultiFab& dkdT, const MultiFab& dedT,
                             const MultiFab& Er_star, const MultiFab& rho,
                             Real delta_t, Real ptc_tau)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rho, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto etaT_arr = etaT[mfi].array();
        auto etaTz_arr = etaTz[mfi].array();
        auto eta1_arr = eta1[mfi].array();
        auto djdT_arr = djdT[mfi].array();
        auto dkdT_arr = dkdT[mfi].array();
        auto dedT_arr = dedT[mfi].array();
        auto Er_star_arr = Er_star[mfi].array();
        auto rho_arr = rho[mfi].array();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            Real sigma = 1.e0_rt + ptc_tau;
            Real cdt = C::c_light * delta_t;

            Real dZdT[NGROUPS];
            Real sumdZdT = 0.0_rt;
            for (int g = 0; g < NGROUPS; ++g) {
                dZdT[g] = djdT_arr(i,j,k,g) - dkdT_arr(i,j,k,g) * Er_star_arr(i,j,k,g);
                sumdZdT += dZdT[g];
            }

            if (sumdZdT == 0.0_rt) {
                sumdZdT = 1.e-50_rt;
            }

            Real foo = cdt * sumdZdT;
            Real bar = sigma * rho_arr(i,j,k) * dedT_arr(i,j,k);
            etaT_arr(i,j,k) = foo / (foo + bar);
            etaTz_arr(i,j,k) = etaT_arr(i,j,k) / sumdZdT;
            eta1_arr(i,j,k) = bar / (foo + bar);
            for (int g = 0; g < NGROUPS; ++g) {
                djdT_arr(i,j,k,g) = dZdT[g] / sumdZdT;
            }
        });
    }
}


void Radiation::eos_opacity_emissivity(const MultiFab& S_new, 
                                       const MultiFab& temp_new,
                                       const MultiFab& temp_star,
                                       MultiFab& kappa_p, MultiFab& kappa_r, MultiFab& jg, 
                                       MultiFab& djdT, MultiFab& dkdT, MultiFab& dedT,
                                       int level, int it, int ngrow)
{
  int star_is_valid = 1 - ngrow;

  int lag_opac;
  if (it == 1) {
    lag_opac = 0;
  }
  else if (it <= update_opacity) {
    lag_opac = 0;
  }
  else {
    lag_opac = 1;
  }

  const Geometry& geom = parent->Geom(level);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& box = mfi.tilebox();

#pragma gpu box(box)
      ca_compute_c_v
          (AMREX_INT_ANYD(box.loVect()), AMREX_INT_ANYD(box.hiVect()),
           BL_TO_FORTRAN_ANYD(dedT[mfi]),
           BL_TO_FORTRAN_ANYD(temp_new[mfi]),
           BL_TO_FORTRAN_ANYD(S_new[mfi]));
  }

  if (dedT_fac > 1.0) {
    dedT.mult(dedT_fac);
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.growntilebox(ngrow);

#pragma gpu box(bx)
      ca_opacs
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           BL_TO_FORTRAN_ANYD(S_new[mfi]),
           BL_TO_FORTRAN_ANYD(temp_new[mfi]),
           BL_TO_FORTRAN_ANYD(temp_star[mfi]),
           BL_TO_FORTRAN_ANYD(kappa_p[mfi]),
           BL_TO_FORTRAN_ANYD(kappa_r[mfi]),
           BL_TO_FORTRAN_ANYD(dkdT[mfi]),
           use_dkdT, star_is_valid, lag_opac);

      const Box& reg = mfi.tilebox();

#pragma gpu box(reg)
      ca_compute_emissivity
          (AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
           BL_TO_FORTRAN_ANYD(jg[mfi]),  
           BL_TO_FORTRAN_ANYD(djdT[mfi]),  
           BL_TO_FORTRAN_ANYD(temp_new[mfi]),
           BL_TO_FORTRAN_ANYD(kappa_p[mfi]),
           BL_TO_FORTRAN_ANYD(dkdT[mfi]));
  }    

  if (ngrow == 0 && !lag_opac) {
      kappa_r.FillBoundary(geom.periodicity());
  }
}


void Radiation::gray_accel(MultiFab& Er_new, MultiFab& Er_pi, 
                           MultiFab& kappa_p, MultiFab& kappa_r,
                           MultiFab& etaT, MultiFab& eta1,
                           MultiFab& mugT,
                           Array<MultiFab, BL_SPACEDIM>& lambda,
                           RadSolve* solver, MGRadBndry& mgbd, 
                           const BoxArray& grids, int level, Real time, 
                           Real delta_t, Real ptc_tau)
{
  const Geometry& geom = parent->Geom(level);
  const Real* dx = parent->Geom(level).CellSize();
  const Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
  const DistributionMapping& dmap = castro->DistributionMap();

  if (nGroups > 1) {
    solver->setHypreMulti(1.0);
  }
  else {
    solver->setHypreMulti(0.0);
  }
  mgbd.setCorrection();

  MultiFab Er_zero(grids, dmap, 1, 0);
  Er_zero.setVal(0.0);
  getBndryDataMG_ga(mgbd, Er_zero, level);

  MultiFab spec(grids, dmap, nGroups, 1);

#ifdef _OPENMP
#pragma omp parallel
#endif 
  for (MFIter mfi(spec, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();

    auto kappa_p_arr = kappa_p[mfi].array();
    auto mugT_arr = mugT[mfi].array();
    auto spec_arr = spec[mfi].array();

    Real cdt1 = 1.e0_rt / (C::c_light * delta_t);

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        Real epsilon[NGROUPS];
        Real sumeps = 0.0;
        for (int g = 0; g < NGROUPS; ++g) {
            Real kapt = kappa_p_arr(i,j,k,g) + (1.e0_rt + ptc_tau) * cdt1;
            epsilon[g] = mugT_arr(i,j,k,g) / kapt;
            sumeps += epsilon[g];
        }

        if (sumeps == 0.e0_rt) {
            for (int g = 0; g < NGROUPS; ++g) {
                spec_arr(i,j,k,g) = 0.e0_rt;
            }
        } else {
            for (int g = 0; g < NGROUPS; ++g) {
                spec_arr(i,j,k,g) = epsilon[g] / sumeps;
            }
        }
    });
  }

  // Extrapolate spectrum out one cell
  for (int indx = 0; indx < nGroups; indx++) {
    extrapolateBorders(spec, indx);
  }
  // Overwrite all extrapolated components with values from
  // neighboring fine grids where they exist:
  spec.FillBoundary(parent->Geom(level).periodicity());

  // set boundary condition
  solver->levelBndry(mgbd,0);

  // A coefficients
  MultiFab acoefs(grids, dmap, 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(acoefs, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
      ca_accel_acoe
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           BL_TO_FORTRAN_ANYD(eta1[mfi]),
           BL_TO_FORTRAN_ANYD(spec[mfi]),
           BL_TO_FORTRAN_ANYD(kappa_p[mfi]),
           BL_TO_FORTRAN_ANYD(acoefs[mfi]),
           delta_t, ptc_tau);    
  }  

  solver->cellCenteredApplyMetrics(level, acoefs);
  solver->setLevelACoeffs(level, acoefs);

  const DistributionMapping& dm = castro->DistributionMap();

  // B & C coefficients
  Array<MultiFab, BL_SPACEDIM> bcoefs, ccoefs, bcgrp;
  for (int idim = 0; idim < BL_SPACEDIM; idim++) {
    const BoxArray& edge_boxes = castro->getEdgeBoxArray(idim);

    bcoefs[idim].define(edge_boxes, dm, 1, 0);
    bcoefs[idim].setVal(0.0);

    bcgrp [idim].define(edge_boxes, dm, 1, 0);

    if (nGroups > 1) {
        ccoefs[idim].define(edge_boxes, dm, 2, 0);
        ccoefs[idim].setVal(0.0);
    }
  }

  for (int igroup = 0; igroup < nGroups; igroup++) {
    for (int idim=0; idim<BL_SPACEDIM; idim++) {
      solver->computeBCoeffs(bcgrp[idim], idim, kappa_r, igroup,
                            lambda[idim], igroup, c, geom);
      // metrics is already in bcgrp

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(spec, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box&  bx  = mfi.nodaltilebox(idim);

#pragma gpu box(bx)
          lbcoefna(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                   BL_TO_FORTRAN_ANYD(bcoefs[idim][mfi]),
                   BL_TO_FORTRAN_ANYD(bcgrp[idim][mfi]),
                   BL_TO_FORTRAN_N_ANYD(spec[mfi], igroup), 
                   idim);
          
          if (nGroups > 1) {
#pragma gpu box(bx)
              ca_accel_ccoe
                  (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                   BL_TO_FORTRAN_ANYD(bcgrp[idim][mfi]),
                   BL_TO_FORTRAN_ANYD(spec[mfi]),
                   BL_TO_FORTRAN_ANYD(ccoefs[idim][mfi]),
                   AMREX_REAL_ANYD(dx), idim, igroup);
          }
      }
    }
  }

  for (int idim = 0; idim < BL_SPACEDIM; idim++) {
    solver->setLevelBCoeffs(level, bcoefs[idim], idim);

    if (nGroups > 1) {
      solver->setLevelCCoeffs(level, ccoefs[idim], idim);
    }
  }

  // rhs
  MultiFab rhs(grids,dmap,1,0);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      auto Ern = Er_new[mfi].array();
      auto Erl = Er_pi[mfi].array();
      auto kap = kappa_p[mfi].array();
      auto etaT_arr = etaT[mfi].array();
      auto rhs_arr = rhs[mfi].array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {
          Real rt_term = 0.0;
          for (int g = 0; g < NGROUPS; ++g) {
              rt_term += kap(i,j,k,g) * (Ern(i,j,k,g) - Erl(i,j,k,g));
          }

          Real H = etaT_arr(i,j,k);

          rhs_arr(i,j,k) = C::c_light * H * rt_term;
      });
  }

  // must apply metrics to rhs here
  solver->cellCenteredApplyMetrics(level, rhs);

  // solve
  MultiFab accel(grids,dmap,1,0);
  accel.setVal(0.0);
  solver->levelSolve(level, accel, 0, rhs, 0.01);

  // update Er_new
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(spec, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& reg  = mfi.tilebox();

      auto Er_new_arr = Er_new[mfi].array();
      auto spec_arr = spec[mfi].array();
      auto accel_arr = accel[mfi].array();

      amrex::ParallelFor(reg,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {
          for (int n = 0; n < NGROUPS; ++n) {
              Er_new_arr(i,j,k,n) = Er_new_arr(i,j,k,n) + spec_arr(i,j,k,n) * accel_arr(i,j,k);
          }
      });
  }

  mgbd.unsetCorrection();
  getBndryDataMG(mgbd, Er_new, time, level);

  solver->restoreHypreMulti();
}


void Radiation::local_accel(MultiFab& Er_new, const MultiFab& Er_pi,
                            const MultiFab& kappa_p, 
                            const MultiFab& etaT,
                            const MultiFab& mugT,
                            Real delta_t, Real ptc_tau)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Er_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto Ern = Er_new[mfi].array();
        auto Erl = Er_pi[mfi].array();
        auto kap = kappa_p[mfi].array();
        auto etaT_arr = etaT[mfi].array();
        auto mugT_arr = mugT[mfi].array();

        Real cdt1 = 1.0_rt / (C::c_light * delta_t);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            Real Hg[NGROUPS], kapt[NGROUPS];
            Real rt_term = 0.0_rt;
            Real p = 1.0_rt;

            for (int g = 0; g < NGROUPS; ++g) {
                Hg[g] = mugT_arr(i,j,k,g) * etaT_arr(i,j,k);
                kapt[g] = kap(i,j,k,g) + (1.e0_rt + ptc_tau) * cdt1;
                Real kk = kap(i,j,k,g) / kapt[g];

                p -= Hg[g] * kk;
                rt_term += kap(i,j,k,g) * (Ern(i,j,k,g) - Erl(i,j,k,g));
            }

            for (int g = 0; g < NGROUPS; ++g) {
                Ern(i,j,k,g) = Ern(i,j,k,g) + (Hg[g] * rt_term) / (kapt[g] * p + 1.e-50_rt);
            }
        });
    }
}


void Radiation::state_energy_update(MultiFab& state, const MultiFab& rhoe, 
                                    const MultiFab& temp, 
                                    const BoxArray& grids,
                                    Real& derat, Real& dT, int level)
{
  BL_PROFILE("Radiation::state_energy_update (MGFLD)");

  const DistributionMapping& dmap = state.DistributionMap();

  MultiFab msk(grids,dmap,1,0);

  {
      BoxArray baf;

      if (level < parent->finestLevel()) {
          baf = parent->boxArray(level+1);
          baf.coarsen(parent->refRatio(level));
      }

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(msk); mfi.isValid(); ++mfi) 
      {
          msk[mfi].setVal<RunOn::Device>(1.0);

          // Now mask off finer level also:
          if (level < parent->finestLevel()) {
              std::vector< std::pair<int,Box> > isects = baf.intersections(mfi.validbox());
              for (int ii = 0; ii < isects.size(); ii++) {
                  msk[mfi].setVal<RunOn::Device>(0.0, isects[ii].second, 0);
              }
          } 
      }
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(state,true); mfi.isValid(); ++mfi) 
  {
      const Box& reg  = mfi.tilebox();

      BL_FORT_PROC_CALL(CA_STATE_UPDATE, ca_state_update)
          (reg.loVect(), reg.hiVect(),
           BL_TO_FORTRAN(state[mfi]),
           BL_TO_FORTRAN(rhoe[mfi]),
           BL_TO_FORTRAN(temp[mfi]),
           BL_TO_FORTRAN(msk[mfi]),
           &derat, &dT);
  }

  ParallelDescriptor::ReduceRealMax(derat);
  ParallelDescriptor::ReduceRealMax(dT);
}


void Radiation::update_matter(MultiFab& rhoe_new, MultiFab& temp_new, 
                              const MultiFab& Er_new, const MultiFab& Er_pi,
                              const MultiFab& rhoe_star,
                              const MultiFab& rhoe_step,
                              const MultiFab& etaT, const MultiFab& etaTz,
                              const MultiFab& eta1,
                              const MultiFab& coupT,
                              const MultiFab& kappa_p, const MultiFab& jg, 
                              const MultiFab& mugT,
                              const MultiFab& S_new, 
                              int level, Real delta_t, Real ptc_tau,
                              int it, bool conservative_update)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rhoe_new,true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto temp_new_arr = temp_new[mfi].array();
        auto S_new_arr = S_new[mfi].array();

        if (conservative_update) {
#pragma gpu box(bx) sync
            ca_update_matter
                (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                 BL_TO_FORTRAN_ANYD(rhoe_new[mfi]),
                 BL_TO_FORTRAN_ANYD(Er_new[mfi]),
                 BL_TO_FORTRAN_ANYD(Er_pi[mfi]),
                 BL_TO_FORTRAN_ANYD(rhoe_star[mfi]),
                 BL_TO_FORTRAN_ANYD(rhoe_step[mfi]),
                 BL_TO_FORTRAN_ANYD(eta1[mfi]),
                 BL_TO_FORTRAN_ANYD(coupT[mfi]),
                 BL_TO_FORTRAN_ANYD(kappa_p[mfi]),
                 delta_t, ptc_tau);

            temp_new[mfi].copy<RunOn::Device>(rhoe_new[mfi],bx);
            Gpu::synchronize();

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                // Get T from rhoe (temp_new comes in with rhoe)

                if (temp_new_arr(i,j,k) <= 0.e0_rt)
                {
                    temp_new_arr(i,j,k) = small_temp;
                }
                else
                {
                    Real rhoInv = 1.e0_rt / S_new_arr(i,j,k,URHO);

                    eos_t eos_state;
                    eos_state.rho = S_new_arr(i,j,k,URHO);
                    eos_state.T   = S_new_arr(i,j,k,UTEMP);
                    eos_state.e   = temp_new_arr(i,j,k) * rhoInv;
                    for (int n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] = S_new_arr(i,j,k,UFS+n) * rhoInv;
                    }
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; ++n) {
                        eos_state.aux[n] = S_new_arr(i,j,k,UFX+n) * rhoInv;
                    }
#endif

                    eos(eos_input_re, eos_state);

                    temp_new_arr(i,j,k) = eos_state.T;
                }
            });
            Gpu::synchronize();
        }
        else {
            BL_FORT_PROC_CALL(CA_NCUPDATE_MATTER, ca_ncupdate_matter)
                (bx.loVect(), bx.hiVect(),
                 BL_TO_FORTRAN(temp_new[mfi]),
                 BL_TO_FORTRAN(Er_new[mfi]),
                 BL_TO_FORTRAN(rhoe_star[mfi]),
                 BL_TO_FORTRAN(rhoe_step[mfi]),
                 BL_TO_FORTRAN(etaTz[mfi]),
                 BL_TO_FORTRAN(kappa_p[mfi]),
                 BL_TO_FORTRAN(jg[mfi]),
                 &delta_t);

#pragma gpu box(bx) sync
            ca_get_rhoe
                (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                 BL_TO_FORTRAN_ANYD(rhoe_new[mfi]),
                 BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
                 BL_TO_FORTRAN_ANYD(S_new[mfi]));
        }
    }  
}

// ========================================================================
// for the hyperbolic solver

void Radiation::compute_limiter(int level, const BoxArray& grids,
                                const MultiFab &Sborder, 
                                const MultiFab &Erborder,
                                MultiFab &lamborder)
{ // it works for both single- and multi-group
  int ngrow = lamborder.nGrow();
  const DistributionMapping& dmap = lamborder.DistributionMap();

  if (filter_lambda_T) {
    BL_ASSERT(ngrow == 4);
  }

  if (limiter == 0) {

    lamborder.setVal(1./3., ngrow);

  }
  else {

    MultiFab kpr(grids,dmap,Radiation::nGroups,ngrow);  

    if (do_multigroup) {
        MGFLD_compute_rosseland(kpr, Sborder); 
    }
    else {
      SGFLD_compute_rosseland(kpr, Sborder); 
    }

    MultiFab Er_wide(grids, dmap, nGroups, ngrow+1);
    Er_wide.setVal(-1.0);
    MultiFab::Copy(Er_wide, Erborder, 0, 0, nGroups, 0);
    
    Er_wide.FillBoundary(parent->Geom(level).periodicity());
    
    const Real* dx = parent->Geom(level).CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(Er_wide,false); mfi.isValid(); ++mfi) {
      BL_FORT_PROC_CALL(CA_COMPUTE_LAMBORDER, ca_compute_lamborder)
        (BL_TO_FORTRAN(Er_wide[mfi]), 
         BL_TO_FORTRAN(kpr[mfi]),
         BL_TO_FORTRAN(lamborder[mfi]), 
         dx, &ngrow, &limiter, &filter_lambda_T, &filter_lambda_S);
    }

    if (filter_lambda_T) {
        lamborder.FillBoundary(parent->Geom(level).periodicity());
    }
  }
}


void Radiation::estimate_gamrPr(const FArrayBox& state, const FArrayBox& Er, 
                                FArrayBox& gPr, const Real*dx, const Box& box)
{
  if (limiter == 0) {
    BL_FORT_PROC_CALL(CA_EST_GPR0, ca_est_gpr0)
      (BL_TO_FORTRAN(Er),
       BL_TO_FORTRAN(gPr));
  }
  else {

    FArrayBox kappa_r(gPr.box(), nGroups);

    if (do_multigroup) {
      MGFLD_compute_rosseland(kappa_r, state);
    }
    else {
      SGFLD_compute_rosseland(kappa_r, state);
    }

    BL_FORT_PROC_CALL(CA_EST_GPR2, ca_est_gpr2)
      (BL_TO_FORTRAN(kappa_r),
       BL_TO_FORTRAN(Er),
       BL_TO_FORTRAN(gPr), 
       box.loVect(),box.hiVect(),
       dx, &limiter, &comoving);
  }
}


void Radiation::MGFLD_compute_rosseland(FArrayBox& kappa_r, const FArrayBox& state)
{
  BL_PROFILE("Radiation::MGFLD_compute_rosseland (FArrayBox)");

  const Box& kbox = kappa_r.box();

#pragma gpu box(kbox) sync
  ca_compute_rosseland(AMREX_INT_ANYD(kbox.loVect()), AMREX_INT_ANYD(kbox.hiVect()),
                       BL_TO_FORTRAN_ANYD(kappa_r),
                       BL_TO_FORTRAN_ANYD(state),
                       0, nGroups-1, nGroups);

}


void Radiation::MGFLD_compute_rosseland(MultiFab& kappa_r, const MultiFab& state)
{
    BL_PROFILE("Radiation::MGFLD_compute_rosseland (MultiFab)");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(kappa_r, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.growntilebox();

#pragma gpu box(bx)
        ca_compute_rosseland(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                             BL_TO_FORTRAN_ANYD(kappa_r[mfi]),
                             BL_TO_FORTRAN_ANYD(state[mfi]),
                             0, nGroups-1, nGroups);
    }
}

void Radiation::MGFLD_compute_scattering(FArrayBox& kappa_s, const FArrayBox& state)
{
    BL_PROFILE("Radiation::MGFLD_compute_scattering");

    const Box& kbox = kappa_s.box();

#pragma gpu box(kbox) sync
    ca_compute_scattering
        (AMREX_INT_ANYD(kbox.loVect()), AMREX_INT_ANYD(kbox.hiVect()),
         BL_TO_FORTRAN_ANYD(kappa_s), 
         BL_TO_FORTRAN_ANYD(state));
}

void Radiation::bisect_matter(MultiFab& rhoe_new, MultiFab& temp_new, 
                              const MultiFab& rhoe_star, const MultiFab& temp_star, 
                              const MultiFab& S_new, const BoxArray& grids, int level)
{
  temp_new.plus(temp_star, 0, 1, 0);
  temp_new.mult(0.5, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(rhoe_new,true); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
      ca_get_rhoe
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           BL_TO_FORTRAN_ANYD(rhoe_new[mfi]),
           BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
           BL_TO_FORTRAN_ANYD(S_new[mfi]));
  }
}


void Radiation::rhstoEr(MultiFab& rhs, Real dt, int level)
{
    const Real* dx = parent->Geom(level).CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter ri(rhs, TilingIfNotGPU()); ri.isValid(); ++ri) 
    {
        const Box& bx = ri.tilebox();       

#pragma gpu box(bx)
        ca_rhstoer(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                   BL_TO_FORTRAN_ANYD(rhs[ri]),
                   AMREX_REAL_ANYD(dx), dt);
    }
}

void Radiation::inelastic_scattering(int level)
{
    if (do_inelastic_scattering) {
        Real dt = parent->dtLevel(level);
        Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
        MultiFab& S_new = castro->get_new_data(State_Type);
        MultiFab& Er_new = castro->get_new_data(Rad_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            FArrayBox kps;
            for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();

                kps.resize(bx,1); // we assume scattering is independent of nu
                MGFLD_compute_scattering(kps, S_new[mfi]);

                ca_inelastic_sct(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                                 BL_TO_FORTRAN_ANYD(S_new[mfi]),
                                 BL_TO_FORTRAN_ANYD(Er_new[mfi]),
                                 BL_TO_FORTRAN_ANYD(kps),
                                 dt);           
            }
        }

        castro->computeTemp(S_new, castro->state[State_Type].curTime(), S_new.nGrow());
    }
}
