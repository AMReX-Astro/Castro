#include <Castro.H>

#ifdef ROTATION
#include <Rotation.H>
#endif

#include <Castro_util.H>
#include <advection_util.H>

#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif

#ifdef RADIATION
#include <Radiation.H>
#include <fluxlimiter.H>
#include <rad_util.H>
#endif

#include <eos.H>

using namespace amrex;


void
Castro::ctoprim(const Box& bx,
                const Real time,
                Array4<Real const> const& uin,
#ifdef MHD
                Array4<Real const> const& Bx,
                Array4<Real const> const& By,
                Array4<Real const> const& Bz,
#endif
#ifdef RADIATION
                Array4<Real const> const& Erin,
                Array4<Real const> const& lam,
#endif
                Array4<Real> const& q_arr,
                Array4<Real> const& qaux_arr) {

  amrex::ignore_unused(time);

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
      auto q = [&] (int n) -> Real& { return q_arr(i,j,k,n); };
      auto qaux = [&] (int n) -> Real& { return qaux_arr(i,j,k,n); };

      hydro::conservative_to_primitive(i, j, k, uin,
#ifdef MHD
                                       Bx, By, Bz,
#endif
#ifdef RADIATION
                                       Erin, lam,
#endif
                                       q, qaux, q_arr.nComp() == NQ);
  });
}


void
Castro::shock(const Box& bx,
              Array4<Real const> const& q_arr,
              Array4<Real const> const& U_src_arr,
              Array4<Real> const& shk) {

  // This is a basic multi-dimensional shock detection algorithm.
  // we look for |grad P . dx| / P > 2/3 and div u < 0
  // This is basically the method in Gronow et al. 2020

  const auto dx = geom.CellSizeArray();
  const int coord_type = geom.Coord();

  Real dxinv = 1.0_rt / dx[0];
#if AMREX_SPACEDIM >= 2
  Real dyinv = 1.0_rt / dx[1];
#endif
#if AMREX_SPACEDIM == 3
  Real dzinv = 1.0_rt / dx[2];
#endif

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real div_u = 0.0_rt;

    // construct div{U}
    if (coord_type == 0) {

      // Cartesian
      div_u += 0.5_rt * (q_arr(i+1,j,k,QU) - q_arr(i-1,j,k,QU)) * dxinv;
#if (AMREX_SPACEDIM >= 2)
      div_u += 0.5_rt * (q_arr(i,j+1,k,QV) - q_arr(i,j-1,k,QV)) * dyinv;
#endif
#if (AMREX_SPACEDIM == 3)
      div_u += 0.5_rt * (q_arr(i,j,k+1,QW) - q_arr(i,j,k-1,QW)) * dzinv;
#endif

#if AMREX_SPACEDIM <= 2
   } else if (coord_type == 1) {

     // r-z
     Real rc = (i + 0.5_rt) * dx[0];
     Real rm = (i - 1 + 0.5_rt) * dx[0];
     Real rp = (i + 1 + 0.5_rt) * dx[0];

#if (AMREX_SPACEDIM == 1)
     div_u += 0.5_rt * (rp * q_arr(i+1,j,k,QU) - rm * q_arr(i-1,j,k,QU)) / (rc * dx[0]);
#endif
#if (AMREX_SPACEDIM == 2)
     div_u += 0.5_rt * (rp * q_arr(i+1,j,k,QU) - rm * q_arr(i-1,j,k,QU)) / (rc * dx[0]) +
              0.5_rt * (q_arr(i,j+1,k,QV) - q_arr(i,j-1,k,QV)) * dyinv;
#endif
#endif

#if AMREX_SPACEDIM == 1
    } else if (coord_type == 2) {

      // 1-d spherical
      Real rc = (i + 0.5_rt) * dx[0];
      Real rm = (i - 1 + 0.5_rt) * dx[0];
      Real rp = (i + 1 + 0.5_rt) * dx[0];

      div_u += 0.5_rt * (rp * rp * q_arr(i+1,j,k,QU) - rm * rm * q_arr(i-1,j,k,QU)) / (rc * rc * dx[0]);
#endif

#ifndef AMREX_USE_GPU

    } else {
      amrex::Error("ERROR: invalid coord_type in shock");
#endif
    }


    // now compute (grad P - rho g) . dx
    // We subtract off the hydrostatic force, since the pressure that
    // balances that is not available to make a shock.
    // We'll use a centered diff for the pressure gradient.
    Real dP_x = 0.5_rt * (q_arr(i+1,j,k,QPRES) - q_arr(i-1,j,k,QPRES));
    if (shock_detection_include_sources == 1) {
        dP_x += -0.25_rt * dx[0] * (U_src_arr(i+1,j,k,UMX) + 2.0_rt * U_src_arr(i,j,k,UMX) + U_src_arr(i-1,j,k,UMX));
    }
    Real dP_y = 0.0_rt;
    Real dP_z = 0.0_rt;
#if AMREX_SPACEDIM >= 2
    dP_y = 0.5_rt * (q_arr(i,j+1,k,QPRES) - q_arr(i,j-1,k,QPRES));
    if (shock_detection_include_sources == 1) {
        dP_y += -0.25_rt * dx[1] * (U_src_arr(i,j+1,k,UMY) + 2.0_rt * U_src_arr(i,j,k,UMY) + U_src_arr(i,j-1,k,UMY));
    }
#endif
#if AMREX_SPACEDIM == 3
    dP_z = 0.5_rt * (q_arr(i,j,k+1,QPRES) - q_arr(i,j,k-1,QPRES));
    if (shock_detection_include_sources == 1) {
        dP_z += -0.25_rt * dx[2] * (U_src_arr(i,j,k+1,UMZ) + 2.0_rt * U_src_arr(i,j,k,UMZ) + U_src_arr(i,j,k-1,UMZ));
    }
#endif

    //Real gradPdx_over_P = std::sqrt(dP_x * dP_x + dP_y * dP_y + dP_z * dP_z) / q_arr(i,j,k,QPRES);
    Real vel = std::sqrt(q_arr(i,j,k,QU) * q_arr(i,j,k,QU) +
                         q_arr(i,j,k,QV) * q_arr(i,j,k,QV) +
                         q_arr(i,j,k,QW) * q_arr(i,j,k,QW));

    Real gradPdx_over_P{0.0_rt};
    if (vel != 0.0) {
        gradPdx_over_P = std::abs(dP_x * q_arr(i,j,k,QU) +
                                  dP_y * q_arr(i,j,k,QV) +
                                  dP_z * q_arr(i,j,k,QW)) / vel;
    }
    gradPdx_over_P /= q_arr(i,j,k,QPRES);

    if (gradPdx_over_P > castro::shock_detection_threshold && div_u < 0.0_rt) {
      shk(i,j,k) = 1.0_rt;
    } else {
      shk(i,j,k) = 0.0_rt;
    }
  });

}


void
Castro::divu(const Box& bx,
             Array4<Real const> const& q_arr,
             Array4<Real> const& div) {
  // this computes the *node-centered* divergence

  const auto dx = geom.CellSizeArray();

#if AMREX_SPACEDIM <= 2
  const int coord_type = geom.Coord();
  const auto problo = geom.ProbLoArray();
#endif

  Real dxinv = 1.0_rt / dx[0];
#if AMREX_SPACEDIM >= 2
  Real dyinv = 1.0_rt / dx[1];
#endif
#if AMREX_SPACEDIM == 3
  Real dzinv = 1.0_rt / dx[2];
#endif

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

#if AMREX_SPACEDIM == 1
    if (coord_type == 0) {
      div(i,j,k) = (q_arr(i,j,k,QU) - q_arr(i-1,j,k,QU)) * dxinv;

    } else if (coord_type == 1) {
      // axisymmetric
      if (i == 0) {
        div(i,j,k) = 0.0_rt;
      } else {
        Real rl = (i - 0.5_rt) * dx[0] + problo[0];
        Real rr = (i + 0.5_rt) * dx[0] + problo[0];
        Real rc = (i) * dx[0] + problo[0];

        div(i,j,k) = (rr * q_arr(i,j,k,QU) - rl * q_arr(i-1,j,k,QU)) * dxinv / rc;
      }
    } else {
      // spherical
      if (i == 0) {
        div(i,j,k) = 0.0_rt;
      } else {
        Real rl = (i - 0.5_rt) * dx[0] + problo[0];
        Real rr = (i + 0.5_rt) * dx[0] + problo[0];
        Real rc = (i) * dx[0] + problo[0];

        div(i,j,k) = (rr * rr * q_arr(i,j,k,QU) - rl * rl * q_arr(i-1,j,k,QU)) * dxinv / (rc * rc);
      }
    }
#endif

#if AMREX_SPACEDIM == 2
    Real ux = 0.0_rt;
    Real vy = 0.0_rt;

    if (coord_type == 0) {
      ux = 0.5_rt * (q_arr(i,j,k,QU) - q_arr(i-1,j,k,QU) + q_arr(i,j-1,k,QU) - q_arr(i-1,j-1,k,QU)) * dxinv;
      vy = 0.5_rt * (q_arr(i,j,k,QV) - q_arr(i,j-1,k,QV) + q_arr(i-1,j,k,QV) - q_arr(i-1,j-1,k,QV)) * dyinv;

    } else {
      if (i == 0) {
        ux = 0.0_rt;
        vy = 0.0_rt;  // is this part correct?
      } else {
        Real rl = (i - 0.5_rt) * dx[0] + problo[0];
        Real rr = (i + 0.5_rt) * dx[0] + problo[0];
        Real rc = (i) * dx[0] + problo[0];

        // These are transverse averages in the y-direction
        Real ul = 0.5_rt * (q_arr(i-1,j,k,QU) + q_arr(i-1,j-1,k,QU));
        Real ur = 0.5_rt * (q_arr(i,j,k,QU) + q_arr(i,j-1,k,QU));

        // Take 1/r d/dr(r*u)
        ux = (rr * ur - rl * ul) * dxinv / rc;

        // These are transverse averages in the x-direction
        Real vb = 0.5_rt * (q_arr(i,j-1,k,QV) + q_arr(i-1,j-1,k,QV));
        Real vt = 0.5_rt * (q_arr(i,j,k,QV) + q_arr(i-1,j,k,QV));

        vy = (vt - vb) * dyinv;
      }
    }

    div(i,j,k) = ux + vy;
#endif

#if AMREX_SPACEDIM == 3
    Real ux = 0.25_rt * (q_arr(i,j,k,QU) - q_arr(i-1,j,k,QU) +
                         q_arr(i,j,k-1,QU) - q_arr(i-1,j,k-1,QU) +
                         q_arr(i,j-1,k,QU) - q_arr(i-1,j-1,k,QU) +
                         q_arr(i,j-1,k-1,QU) - q_arr(i-1,j-1,k-1,QU)) * dxinv;

    Real vy = 0.25_rt * (q_arr(i,j,k,QV) - q_arr(i,j-1,k,QV) +
                         q_arr(i,j,k-1,QV) - q_arr(i,j-1,k-1,QV) +
                         q_arr(i-1,j,k,QV) - q_arr(i-1,j-1,k,QV) +
                         q_arr(i-1,j,k-1,QV) - q_arr(i-1,j-1,k-1,QV)) * dyinv;

    Real wz = 0.25_rt * (q_arr(i,j,k,QW) - q_arr(i,j,k-1,QW) +
                         q_arr(i,j-1,k,QW) - q_arr(i,j-1,k-1,QW) +
                         q_arr(i-1,j,k,QW) - q_arr(i-1,j,k-1,QW) +
                         q_arr(i-1,j-1,k,QW) - q_arr(i-1,j-1,k-1,QW)) * dzinv;

    div(i,j,k) = ux + vy + wz;
#endif

  });

}


void
Castro::apply_av(const Box& bx,
                 const int idir,
                 Array4<Real const> const& div,
                 Array4<Real const> const& uin,
                 Array4<Real> const& flux) {

  const auto dx = geom.CellSizeArray();

  Real diff_coeff = difmag;

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

      Real div1;
      if (idir == 0) {
          div1 = 0.25_rt * (div(i,j,k) + div(i,j+dg1,k) +
                            div(i,j,k+dg2) + div(i,j+dg1,k+dg2));
      } else if (idir == 1) {
          div1 = 0.25_rt * (div(i,j,k) + div(i+1,j,k) +
                            div(i,j,k+dg2) + div(i+1,j,k+dg2));
      } else {
          div1 = 0.25_rt * (div(i,j,k) + div(i+1,j,k) +
                            div(i,j+dg1,k) + div(i+1,j+dg1,k));
      }

      div1 = diff_coeff * std::min(0.0_rt, div1);

      if (div1 == 0.0_rt) {
          return;
      }

      // we want to prevent the artificial viscosity from breaking
      // species conservation (sum X_k = 1).  This can happen if the
      // artificial viscosity flips the sign of some of the fluxes.
      // Here we scale div1 to ensure that this is not an issue.

      if (limit_avisc_on_species) {

          amrex::Real limit_div1{div1};

          for (int n = 0; n < NumSpec; ++n) {

              Real du{};

              if (idir == 0) {
                  du = uin(i,j,k,UFS+n) - uin(i-1,j,k,UFS+n);
              } else if (idir == 1) {
                  du = uin(i,j,k,UFS+n) - uin(i,j-dg1,k,UFS+n);
              } else {
                  du = uin(i,j,k,UFS+n) - uin(i,j,k-dg2,UFS+n);
              }

              Real fnew = flux(i,j,k,UFS+n) + dx[idir] * div1 * du;
              if (fnew * flux(i,j,k,UFS+n) < 0) {
                  limit_div1 = std::min(limit_div1,
                                        std::abs(flux(i,j,k,UFS+n) / (dx[idir] * du)));
              }
          }

          div1 = limit_div1;

      }


      for (int n = 0; n < NUM_STATE; ++n) {

          if (n == UTEMP) {
              continue;
          }
#ifdef SHOCK_VAR
          if (n == USHK) {
              continue;
          }
#endif
#ifdef NSE_NET
          if (n == UMUP || n == UMUN) {
              continue;
          }
#endif

          Real div_var{};

          if (idir == 0) {
              div_var = div1 * (uin(i,j,k,n) - uin(i-1,j,k,n));
          } else if (idir == 1) {
              div_var = div1 * (uin(i,j,k,n) - uin(i,j-dg1,k,n));
          } else {
              div_var = div1 * (uin(i,j,k,n) - uin(i,j,k-dg2,n));
          }

          flux(i,j,k,n) += dx[idir] * div_var;
      }
  });
}


#ifdef RADIATION
void
Castro::apply_av_rad(const Box& bx,
                     const int idir,
                     Array4<Real const> const& div,
                     Array4<Real const> const& Erin,
                     Array4<Real> const& radflux) {

  const auto dx = geom.CellSizeArray();

  Real diff_coeff = difmag;

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

      Real div1;
      if (idir == 0) {
          div1 = 0.25_rt * (div(i,j,k) + div(i,j+dg1,k) +
                            div(i,j,k+dg2) + div(i,j+dg1,k+dg2));
      } else if (idir == 1) {
          div1 = 0.25_rt * (div(i,j,k) + div(i+1,j,k) +
                            div(i,j,k+dg2) + div(i+1,j,k+dg2));
      } else {
          div1 = 0.25_rt * (div(i,j,k) + div(i+1,j,k) +
                            div(i,j+dg1,k) + div(i+1,j+dg1,k));
      }

      div1 = diff_coeff * std::min(0.0_rt, div1);

      if (div1 == 0.0_rt) {
          return;
      }

      for (int n = 0; n < Radiation::nGroups; ++n) {

          Real div_var{};

          if (idir == 0) {
              div_var = div1 * (Erin(i,j,k,n) - Erin(i-1,j,k,n));
          } else if (idir == 1) {
              div_var = div1 * (Erin(i,j,k,n) - Erin(i,j-dg1,k,n));
          } else {
              div_var = div1 * (Erin(i,j,k,n) - Erin(i,j,k-dg2,n));
          }

          radflux(i,j,k,n) += dx[idir] * div_var;
      }
  });
}
#endif


void
Castro::normalize_species_fluxes(const Box& bx,
                                 Array4<Real> const& flux) {

  // Normalize the fluxes of the mass fractions so that
  // they sum to 0.  This is essentially the CMA procedure that is
  // defined in Plewa & Muller, 1999, A&A, 342, 179.

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    Real sum = 0.0_rt;

    for (int n = UFS; n < UFS+NumSpec; n++) {
      sum += flux(i,j,k,n);
    }

    Real fac = 1.0_rt;

    // We skip the normalization if the sum is zero or within epsilon.
    // There can be numerical problems here if the density flux is
    // approximately zero at the interface but not exactly, resulting in
    // division by a small number and/or resulting in one of the species
    // fluxes being negative because of roundoff error. There are also other
    // terms like artificial viscosity which can cause these problems.
    // So checking that sum is sufficiently large helps avoid this.

    if (std::abs(sum) > std::numeric_limits<Real>::epsilon() * std::abs(flux(i,j,k,URHO))) {
      fac = flux(i,j,k,URHO) / sum;
    }

    for (int n = UFS; n < UFS+NumSpec; n++) {
      flux(i,j,k,n) = flux(i,j,k,n) * fac;
    }
  });
}


void  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
Castro::scale_flux(const Box& bx,
#if AMREX_SPACEDIM == 1
                   Array4<Real const> const& qint,
#endif
                   Array4<Real> const& flux,
                   Array4<Real const> const& area_arr,
                   const Real dt) {

#if AMREX_SPACEDIM == 1
  const int coord_type = geom.Coord();
#endif

  amrex::ParallelFor(bx, NUM_STATE,
  [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {

    flux(i,j,k,n) = dt * flux(i,j,k,n) * area_arr(i,j,k);
#if AMREX_SPACEDIM == 1
    // Correct the momentum flux with the grad p part.
    if (coord_type == 0 && n == UMX) {
      flux(i,j,k,n) += dt * area_arr(i,j,k) * qint(i,j,k,GDPRES);
    }
#endif
  });
}


#ifdef RADIATION
void
Castro::scale_rad_flux(const Box& bx,
                       Array4<Real> const& rflux,
                       Array4<Real const> const& area_arr,
                       const Real dt) {

  amrex::ParallelFor(bx, Radiation::nGroups,
  [=] AMREX_GPU_DEVICE (int i, int j, int k, int g) noexcept
  {
    rflux(i,j,k,g) = dt * rflux(i,j,k,g) * area_arr(i,j,k);
  });
}
#endif


#ifndef MHD
void
Castro::limit_hydro_fluxes_on_small_dens(const Box& bx,
                                         int idir,
                                         Array4<Real const> const& u,
                                         Array4<Real const> const& vol,
                                         Array4<Real> const& flux,
                                         Array4<Real const> const& area_arr,
                                         Real dt,
                                         bool scale_by_dAdt)
{
    // Hu, Adams, and Shu (2013), JCP, 242, 169, "Positivity-preserving method for
    // high-order conservative schemes solving compressible Euler equations," proposes
    // a positivity-preserving advection scheme. That algorithm blends the actual
    // (second-order) hydro flux with the first-order Lax-Friedrichs flux to ensure
    // positivity. However, we demand a stronger requirement, that rho > small_dens.
    // Additionally, the blending approach can cause problems with multiple advecting
    // species, since adding the Lax-Friedrichs flux does not guarantee physical consistency
    // between the species advection and the density advection. So instead of trying the
    // blending approach, we simply apply a linear scaling to each flux such that it does
    // not violate the density floor.

    const Real density_floor_tolerance = 1.1_rt;

    // The density floor is the small density, modified by a small factor.
    // In practice numerical error can cause the density that is created
    // by this flux limiter to be slightly lower than the target density,
    // so we set the target to be slightly larger than the real density floor
    // to avoid density resets.

    Real density_floor = small_dens * density_floor_tolerance;

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Grab the states on either side of the interface we are working with,
        // depending on which dimension we're currently calling this with.

        Real rhoR = u(i,j,k,URHO);
        Real volR = vol(i,j,k);

        Real rhoL, volL;

        if (idir == 0) {
            rhoL = u(i-1,j,k,URHO);
            volL = vol(i-1,j,k);
        }
        else if (idir == 1) {
            rhoL = u(i,j-1,k,URHO);
            volL = vol(i,j-1,k);
        }
        else {
            rhoL = u(i,j,k-1,URHO);
            volL = vol(i,j,k-1);
        }

        // Coefficients of fluxes on either side of the interface.

        Real flux_coefR = 1.0_rt / volR;
        Real flux_coefL = 1.0_rt / volL;

        if (scale_by_dAdt) {
            flux_coefR *= dt * area_arr(i,j,k);
            flux_coefL *= dt * area_arr(i,j,k);
        }

        // Updates to the zones on either side of the interface.

        Real drhoR = flux_coefR * flux(i,j,k,URHO);
        Real drhoL = flux_coefL * flux(i,j,k,URHO);

        // Limit all fluxes such that the zone does not go negative in density.
        // We apply this flux limiter on a per-edge basis. So we can guarantee
        // that any individual flux cannot cause a small density in one step,
        // but with the above floor we cannot guarantee that the sum of the
        // fluxes will enforce this constraint. A simple way to guarantee that
        // is to note that in the worst case, all six faces adjacent to a zone
        // could be evacuating that zone (more generally, 2 * AMREX_SPACEDIM).
        // So instead of checking whether rho + drho is less than the floor,
        // we can check whether rho + 6 * drho is less than the floor, and
        // strengthen the limiter by that same factor of 6, so that even in the
        // worst case, the zone will remain positive. This is simple to implement
        // at the risk of being overly aggressive on flux limiting. A more
        // sophisticated version of this might check the ratio of the fluxes
        // (or, similarly, the ratio of the edge velocities coming from the
        // Riemann solver) and scale the larger fluxes more strongly than the
        // weaker fluxes, and perhaps only consider fluxes that subtract from
        // the zone, but this would also be more complicated to implement.

        if (rhoR >= density_floor && std::abs(drhoR) > 0.0_rt && rhoR + 2 * AMREX_SPACEDIM * drhoR < density_floor) {
            Real limiting_factor = std::abs((density_floor - rhoR) / (2 * AMREX_SPACEDIM * drhoR));
            for (int n = 0; n < NUM_STATE; ++n) {
                flux(i,j,k,n) = flux(i,j,k,n) * limiting_factor;
            }
        }
        else if (rhoL >= density_floor && std::abs(drhoL) > 0.0_rt && rhoL - 2 * AMREX_SPACEDIM * drhoL < density_floor) {
            Real limiting_factor = std::abs((density_floor - rhoL) / (2 * AMREX_SPACEDIM * drhoL));
            for (int n = 0; n < NUM_STATE; ++n) {
                flux(i,j,k,n) = flux(i,j,k,n) * limiting_factor;
            }
        }

    });
}
#endif


void  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
Castro::do_enforce_minimum_density(const Box& bx,
                                   Array4<Real> const& state_arr,
                                   const int verbose_warnings) {

#ifdef HYBRID_MOMENTUM
  GeometryData geomdata = geom.data();
#endif

  amrex::ignore_unused(verbose_warnings);

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    if (state_arr(i,j,k,URHO) < small_dens) {

#ifndef AMREX_USE_GPU
      if (verbose_warnings > 1 ||
          (verbose_warnings > 0 && state_arr(i,j,k,URHO) > castro::retry_small_density_cutoff)) {
        std::cout << " " << std::endl;
        if (state_arr(i,j,k,URHO) < 0.0_rt) {
          std::cout << ">>> RESETTING NEG.  DENSITY AT " << i << ", " << j << ", " << k << std::endl;
        }
        else if (state_arr(i,j,k,URHO) == 0.0_rt) {
          // If the density is *exactly* zero, that almost certainly means something has gone wrong,
          // like we failed to properly fill the state data on grid creation.
          amrex::Error("Density exactly zero at " + std::to_string(i) + ", " +
                                                    std::to_string(j) + ", " +
                                                    std::to_string(k));
        }
        else {
          std::cout << ">>> RESETTING SMALL DENSITY AT " << i << ", " << j << ", " << k << std::endl;
        }
        std::cout << ">>> FROM " << state_arr(i,j,k,URHO) << " TO " << small_dens << std::endl;
        std::cout << ">>> IN GRID " << bx << std::endl;
        std::cout << " " << std::endl;
      }
#endif

      for (int ipassive = 0; ipassive < npassive; ipassive++) {
          const int n = upassmap(ipassive);
          state_arr(i,j,k,n) *= (small_dens / state_arr(i,j,k,URHO));
      }

      eos_re_t eos_state;
      eos_state.rho = small_dens;
      eos_state.T = small_temp;
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state_arr(i,j,k,UFS+n) / small_dens;
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = state_arr(i,j,k,UFX+n) / small_dens;
      }
#endif

      eos(eos_input_rt, eos_state);

      state_arr(i,j,k,URHO ) = eos_state.rho;
      state_arr(i,j,k,UTEMP) = eos_state.T;

      state_arr(i,j,k,UMX) = 0.0_rt;
      state_arr(i,j,k,UMY) = 0.0_rt;
      state_arr(i,j,k,UMZ) = 0.0_rt;

      state_arr(i,j,k,UEINT) = eos_state.rho * eos_state.e;
      state_arr(i,j,k,UEDEN) = state_arr(i,j,k,UEINT);

#ifdef HYBRID_MOMENTUM
      GpuArray<Real, 3> loc;

      position(i, j, k, geomdata, loc);

      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        loc[dir] -= problem::center[dir];
      }

      GpuArray<Real, 3> linear_mom;

      for (int dir = 0; dir < 3; ++dir) {
        linear_mom[dir] = state_arr(i,j,k,UMX+dir);
      }

      GpuArray<Real, 3> hybrid_mom;

      linear_to_hybrid(loc, linear_mom, hybrid_mom);

      for (int dir = 0; dir < 3; ++dir) {
        state_arr(i,j,k,UMR+dir) = hybrid_mom[dir];
      }
#endif
    }
  });
}


void
Castro::enforce_reflect_states(const Box& bx, const int idir,
                               Array4<Real> const& qm,
                               Array4<Real> const& qp) {

    // special care for reflecting BCs
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();

    const auto domlo = geom.Domain().loVect3d();
    const auto domhi = geom.Domain().hiVect3d();

    bool lo_bc_test = lo_bc[idir] == amrex::PhysBCType::symmetry;
    bool hi_bc_test = hi_bc[idir] == amrex::PhysBCType::symmetry;

    // normal velocity
    const int QUN = QU + idir;

    // this is a loop over interfaces

    if (lo_bc_test) {

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            // reset the left state at domlo if needed -- it is outside the domain
            if ((idir == 0 && i == domlo[0]) ||
                (idir == 1 && j == domlo[1]) ||
                (idir == 2 && k == domlo[2])) {
                for (int n = 0; n < NQ; n++) {
                    if (n == QUN) {
                        qm(i,j,k,QUN) = -qp(i,j,k,QUN);
                    } else {
                        qm(i,j,k,n) = qp(i,j,k,n);
                    }
                }
            }
        });
    }

    if (hi_bc_test) {

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            // reset the right state at domhi+1 if needed -- it is outside the domain
            if ((idir == 0 && i == domhi[0]+1) ||
                (idir == 1 && j == domhi[1]+1) ||
                (idir == 2 && k == domhi[2]+1)) {
                for (int n = 0; n < NQ; n++) {
                    if (n == QUN) {
                        qp(i,j,k,QUN) = -qm(i,j,k,QUN);
                    } else {
                        qp(i,j,k,n) = qm(i,j,k,n);
                    }
                }
            }
        });

    }
}
