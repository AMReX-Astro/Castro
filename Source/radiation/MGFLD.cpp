#include <Radiation.H>
#include <Castro_F.H>
#include <RAD_F.H>
#include <rad_util.H>
#include <blackbody.H>
#include <opacity.H>
#include <problem_emissivity.H>

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

      auto dedT_arr = dedT[mfi].array();
      auto temp_arr = temp_new[mfi].array();
      auto S_new_arr = S_new[mfi].array();

      amrex::ParallelFor(box,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
      {
          Real rhoInv = 1.e0_rt / S_new_arr(i,j,k,URHO);

          eos_t eos_state;
          eos_state.rho = S_new_arr(i,j,k,URHO);
          eos_state.T   = temp_arr(i,j,k);
          for (int n = 0; n < NumSpec; ++n) {
              eos_state.xn[n] = S_new_arr(i,j,k,UFS+n) * rhoInv;
          }
#if NAUX_NET > 0
          for (int n = 0; n < NumAux; ++n) {
              eos_state.aux[n] = S_new_arr(i,j,k,UFX+n) * rhoInv;
          }
#endif

          eos(eos_input_rt, eos_state);

          dedT_arr(i,j,k) = eos_state.cv;
      });
  }

  if (dedT_fac > 1.0) {
    dedT.mult(dedT_fac);
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.growntilebox(ngrow);

      auto S_new_arr = S_new[mfi].array();
      auto temp_new_arr = temp_new[mfi].array();
      auto temp_star_arr = temp_star[mfi].array();
      auto kappa_p_arr = kappa_p[mfi].array();
      auto kappa_r_arr = kappa_r[mfi].array();
      auto dkdT_arr = dkdT[mfi].array();
      auto jg_arr = jg[mfi].array();
      auto djdT_arr = djdT[mfi].array();

      bool use_dkdT_loc = use_dkdT;

      GpuArray<Real, NGROUPS> nugroup_loc;
      for (int g = 0; g < NGROUPS; ++g) {
          nugroup_loc[g] = nugroup[g];
      }
      GpuArray<Real, NGROUPS+1> xnu_loc;
      for (int g = 0; g < NGROUPS+1; ++g) {
          xnu_loc[g] = xnu[g];
      }

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {
          const Real fac = 0.5e0_rt;
          const Real minfrac = 1.e-8_rt;

          if (lag_opac) {
              dkdT_arr(i,j,k) = 0.0_rt;
              return;
          }

          Real rho = S_new_arr(i,j,k,URHO);
          Real temp = temp_new_arr(i,j,k);

          Real Ye;
          if (NumAux > 0) {
              Real Ye = S_new_arr(i,j,k,UFX);
          } else {
              Ye = 0.e0_rt;
          }

          Real dT;
          if (star_is_valid > 0) {
              dT = fac * std::abs(temp_star_arr(i,j,k) - temp_new_arr(i,j,k));
              dT = amrex::max(dT, minfrac * temp_new_arr(i,j,k));
          } else {
              dT = temp_new_arr(i,j,k) * 1.e-3_rt + 1.e-50_rt;
          }

          for (int g = 0; g < NGROUPS; ++g) {
              Real nu = nugroup_loc[g];

              bool comp_kp = true;
              bool comp_kr = true;

              Real kp, kr;

              opacity(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr);

              kappa_p_arr(i,j,k,g) = kp;
              kappa_r_arr(i,j,k,g) = kr;

              if (use_dkdT_loc == 0) {

                  dkdT_arr(i,j,k,g) = 0.e0_rt;

              } else {

                  bool comp_kp = true;
                  bool comp_kr = false;

                  Real kp1, kr1, kp2, kr2;

                  opacity(kp1, kr1, rho, temp-dT, Ye, nu, comp_kp, comp_kr);
                  opacity(kp2, kr2, rho, temp+dT, Ye, nu, comp_kp, comp_kr);

                  dkdT_arr(i,j,k,g) = (kp2 - kp1) / (2.e0_rt * dT);
              }
          }
      });

      const Box& reg = mfi.tilebox();

      amrex::ParallelFor(reg,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {
          // Integrate the Planck distribution upward from zero frequency.
          // This handles both the single-group and multi-group cases.

          Real Teff = amrex::max(temp_new_arr(i,j,k), 1.e-50_rt);

          Real B1, dBdT1;
          BdBdTIndefInteg(Teff, 0.0_rt, B1, dBdT1);

          for (int g = 0; g < NGROUPS; ++g) {

              Real xnup = xnu_loc[g+1];

              // For the last group, make sure that we complete
              // the integral up to "infinity".

              if (g == NGROUPS - 1) {
                  xnup = amrex::max(xnup, 1.e25_rt);
              }

              Real B0 = B1;
              Real dBdT0 = dBdT1;
              BdBdTIndefInteg(Teff, xnup, B1, dBdT1);
              Real Bg = B1 - B0;
              Real dBdT = dBdT1 - dBdT0;

              jg_arr(i,j,k,g) = Bg * kappa_p_arr(i,j,k,g);
              djdT_arr(i,j,k,g) = dkdT_arr(i,j,k,g) * Bg + dBdT * kappa_p_arr(i,j,k,g);

              // Allow a problem to override this emissivity.

              problem_emissivity(i, j, k, g,
                                 nugroup_loc, xnu_loc,
                                 temp_new_arr(i,j,k), kappa_p_arr(i,j,k,g),
                                 dkdT_arr(i,j,k,g), jg_arr(i,j,k,g), djdT_arr(i,j,k,g));
          }
      });
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
  auto dx = parent->Geom(level).CellSizeArray();
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

      auto eta1_arr = eta1[mfi].array();
      auto spec_arr = spec[mfi].array();
      auto kappa_p_arr = kappa_p[mfi].array();
      auto acoefs_arr = acoefs[mfi].array();

      const Real dt1 = (1.e0_rt + ptc_tau) / delta_t;

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {
          Real kbar = 0.0;
          for (int g = 0; g < NGROUPS; ++g) {
              kbar += spec_arr(i,j,k,g) * kappa_p_arr(i,j,k,g);
          }

          Real H1 = eta1_arr(i,j,k);

          acoefs_arr(i,j,k) = H1 * kbar * C::c_light + dt1;
      });
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

          auto bcoefs_arr = bcoefs[idim][mfi].array();
          auto bcgrp_arr = bcgrp[idim][mfi].array();
          auto spec_arr = spec[mfi].array(igroup);
          auto ccoefs_arr = ccoefs[idim][mfi].array();

          amrex::ParallelFor(bx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
          {
              if (idim == 0) {
                  bcoefs_arr(i,j,k) += 0.5e0_rt * (spec_arr(i-1,j,k) + spec_arr(i,j,k)) * bcgrp_arr(i,j,k);
              }
              else if (idim == 1) {
                  bcoefs_arr(i,j,k) += 0.5e0_rt * (spec_arr(i,j-1,k) + spec_arr(i,j,k)) * bcgrp_arr(i,j,k);
              }
              else {
                  bcoefs_arr(i,j,k) += 0.5e0_rt * (spec_arr(i,j,k-1) + spec_arr(i,j,k)) * bcgrp_arr(i,j,k);
              }
          });
          
          if (nGroups > 1) {
              amrex::ParallelFor(bx,
              [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
              {
                  int ioff, joff, koff;
                  Real h1;

                  if (idim == 0) {
                      ioff = 1;
                      joff = 0;
                      koff = 0;
                      h1 = 1.e0_rt / dx[0];
                  }
                  else if (idim == 1) {
                      ioff = 0;
                      joff = 1;
                      koff = 0;
                      h1 = 1.e0_rt / dx[1];
                  }
                  else {
                      ioff = 0;
                      joff = 0;
                      koff = 1;
                      h1 = 1.e0_rt / dx[2];
                  }

                  Real grad_spec = (spec_arr(i,j,k) - spec_arr(i-ioff,j-joff,k-koff)) * h1;
                  Real foo = -0.5e0_rt * bcgrp_arr(i,j,k) * grad_spec;
                  ccoefs_arr(i,j,k,0) += foo;
                  ccoefs_arr(i,j,k,1) += foo;
              });
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
    const Real cdt = C::c_light * delta_t;
    const Real cdt1 = 1.e0_rt / (C::c_light * delta_t);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rhoe_new,true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto re_n = rhoe_new[mfi].array();
        auto S_new_arr = S_new[mfi].array();
        auto Tp_n = temp_new[mfi].array();
        auto Er_n = Er_new[mfi].array();
        auto Er_l = Er_pi[mfi].array();
        auto re_s = rhoe_star[mfi].array();
        auto re_2 = rhoe_step[mfi].array();
        auto eta1_arr = eta1[mfi].array();
        auto cpT = coupT[mfi].array();
        auto etTz = etaTz[mfi].array();
        auto kpp  = kappa_p[mfi].array();
        auto jg_arr = jg[mfi].array();

        if (conservative_update) {

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                Real H1 = eta1_arr(i,j,k);

                Real dkEE = 0.0;
                for (int g = 0; g < NGROUPS; ++g) {
                    dkEE += kpp(i,j,k,g) * (Er_n(i,j,k,g) - Er_l(i,j,k,g));
                }

                Real chg = cdt * dkEE + H1 * ((re_2(i,j,k) - re_s(i,j,k)) + cdt * cpT(i,j,k));

                re_n(i,j,k) = re_s(i,j,k) + chg;

                re_n(i,j,k) = (re_n(i,j,k) + ptc_tau * re_s(i,j,k)) / (1.e0_rt + ptc_tau);
            });

            temp_new[mfi].copy<RunOn::Device>(rhoe_new[mfi],bx);
            Gpu::synchronize();

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                // Get T from rhoe (temp_new comes in with rhoe)

                if (Tp_n(i,j,k) <= 0.e0_rt)
                {
                    Tp_n(i,j,k) = small_temp;
                }
                else
                {
                    Real rhoInv = 1.e0_rt / S_new_arr(i,j,k,URHO);

                    eos_t eos_state;
                    eos_state.rho = S_new_arr(i,j,k,URHO);
                    eos_state.T   = S_new_arr(i,j,k,UTEMP);
                    eos_state.e   = Tp_n(i,j,k) * rhoInv;
                    for (int n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] = S_new_arr(i,j,k,UFS+n) * rhoInv;
                    }
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; ++n) {
                        eos_state.aux[n] = S_new_arr(i,j,k,UFX+n) * rhoInv;
                    }
#endif

                    eos(eos_input_re, eos_state);

                    Tp_n(i,j,k) = eos_state.T;
                }
            });
            Gpu::synchronize();
        }
        else {

            const Real fac = 0.01_rt;

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                Real cpT = 0.e0_rt;

                for (int g = 0; g < NGROUPS; ++g) {
                    cpT = cpT + kpp(i,j,k,g) * Er_n(i,j,k,g) - jg_arr(i,j,k,g);
                }

                Real scrch_re = cpT - (re_s(i,j,k) - re_2(i,j,k)) * cdt1;

                Real dTemp = etTz(i,j,k) * scrch_re;

                if (std::abs(dTemp / (Tp_n(i,j,k) + 1.e-50_rt)) > fac) {
                    dTemp = std::copysign(fac * Tp_n(i,j,k), dTemp);
                }

                Tp_n(i,j,k) = Tp_n(i,j,k) + dTemp;

                Real rhoInv = 1.e0_rt / S_new_arr(i,j,k,URHO);

                eos_t eos_state;
                eos_state.rho = S_new_arr(i,j,k,URHO);
                eos_state.T   = Tp_n(i,j,k);
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = S_new_arr(i,j,k,UFS+n) * rhoInv;
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; ++n) {
                    eos_state.aux[n] = S_new_arr(i,j,k,UFX+n) * rhoInv;
                }
#endif

                eos(eos_input_rt, eos_state);

                re_n(i,j,k) = eos_state.rho * eos_state.e;
            });
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

  auto state_arr = state.array();
  auto kpr = kappa_r.array();

  GpuArray<Real, NGROUPS> nugroup_loc;
  for (int g = 0; g < NGROUPS; ++g) {
      nugroup_loc[g] = nugroup[g];
  }

  amrex::ParallelFor(kbox,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {
      Real rho = state_arr(i,j,k,URHO);
      Real temp = state_arr(i,j,k,UTEMP);
      Real Ye;
      if (NumAux > 0) {
          Ye = state_arr(i,j,k,UFX);
      } else {
          Ye = 0.e0_rt;
      }

      Real kp, kr;
      bool comp_kp = false;
      bool comp_kr = true;

      for (int g = 0; g < NGROUPS; ++g) {
          opacity(kp, kr, rho, temp, Ye, nugroup_loc[g], comp_kp, comp_kr);
          kpr(i,j,k,g) = kr;
      }
  });
  Gpu::synchronize();
}


void Radiation::MGFLD_compute_rosseland(MultiFab& kappa_r, const MultiFab& state)
{
    BL_PROFILE("Radiation::MGFLD_compute_rosseland (MultiFab)");

    GpuArray<Real, NGROUPS> nugroup_loc;
    for (int g = 0; g < NGROUPS; ++g) {
        nugroup_loc[g] = nugroup[g];
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(kappa_r, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.growntilebox();

        auto state_arr = state[mfi].array();
        auto kpr = kappa_r[mfi].array();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            Real rho = state_arr(i,j,k,URHO);
            Real temp = state_arr(i,j,k,UTEMP);
            Real Ye;
            if (NumAux > 0) {
                Ye = state_arr(i,j,k,UFX);
            } else {
                Ye = 0.e0_rt;
            }

            Real kp, kr;
            bool comp_kp = false;
            bool comp_kr = true;

            for (int g = 0; g < NGROUPS; ++g) {
                opacity(kp, kr, rho, temp, Ye, nugroup_loc[g], comp_kp, comp_kr);
                kpr(i,j,k,g) = kr;
            }
        });
    }
}

void Radiation::MGFLD_compute_scattering(FArrayBox& kappa_s, const FArrayBox& state)
{
    BL_PROFILE("Radiation::MGFLD_compute_scattering");

    const Box& kbox = kappa_s.box();

    auto kps = kappa_s.array();
    auto sta = state.array();

    // scattering is assumed to be independent of nu.
    const Real nu = nugroup[0];

    amrex::ParallelFor(kbox,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        Real rho = sta(i,j,k,URHO);
        Real temp = sta(i,j,k,UTEMP);
        Real Ye;
        if (NumAux > 0) {
            Ye = sta(i,j,k,UFX);
        }
        else {
            Ye = 0.e0_rt;
        }

        Real kp, kr;
        bool comp_kp = true;
        bool comp_kr = true;
        opacity(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr);

        kps(i,j,k) = amrex::max(kr - kp, 0.e0_rt);
    });
    Gpu::synchronize();
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

      auto rhoe = rhoe_new[mfi].array();
      auto temp = temp_new[mfi].array();
      auto state = S_new[mfi].array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {
          Real rhoInv = 1.e0_rt / state(i,j,k,URHO);

          eos_t eos_state;
          eos_state.rho = state(i,j,k,URHO);
          eos_state.T   =  temp(i,j,k);
          for (int n = 0; n < NumSpec; ++n) {
              eos_state.xn[n] = state(i,j,k,UFS+n) * rhoInv;
          }
#if NAUX_NET > 0
          for (int n = 0; n < NumAux; ++n) {
              eos_state.aux[n] = state(i,j,k,UFX+n) * rhoInv;
          }
#endif

          eos(eos_input_rt, eos_state);

          rhoe(i,j,k) = eos_state.rho * eos_state.e;
      });
  }
}


void Radiation::rhstoEr(MultiFab& rhs, Real dt, int level)
{
    auto geomdata = parent->Geom(level).data();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter ri(rhs, TilingIfNotGPU()); ri.isValid(); ++ri) 
    {
        const Box& bx = ri.tilebox();

        auto rhs_arr = rhs[ri].array();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            Real r, s;
            cell_center_metric(i, j, k, geomdata, r, s);

            rhs_arr(i,j,k) *= dt / r;
        });
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
