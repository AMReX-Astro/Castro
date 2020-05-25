#include "Castro.H"

void
Castro::fourth_interfaces(const Box& bx,
                          const int idir, const int ncomp,
                          Array4<Real const> const& a,
                          Array4<Real> const& a_int) {

  // this just computes the unlimited single-value interface state
  // for the 4th order method.
  //
  // Note: this needs to be run on lo-1 to [hi(1)+2, hi(2)+1, hi(3)+1] for x,
  // and analogously for y and z

  // our convention here is that:
  //     al(i,j,k)   will be al_{i-1/2,j,k),
  //     al(i+1,j,k) will be al_{i+1/2,j,k)

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  auto domlo = geom.Domain().loVect3d();
  auto domhi = geom.Domain().hiVect3d();

  if (idir == 0) {

    // this loop is over interfaces

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {

      // interpolate to the edges -- this is a_{i-1/2}
      // note for non-periodic physical boundaries, we use a special stencil

      if (i == domlo[0]+1 && lo_bc[0] != Interior) {
        // use a stencil for the interface that is one zone
        // from the left physical boundary, MC Eq. 22
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i-1,j,k,ncomp) + 13.0_rt*a(i,j,k,ncomp) -
                                         5.0_rt*a(i+1,j,k,ncomp) + a(i+2,j,k,ncomp));

      } else if (i == domlo[0] && lo_bc[0] != Interior) {
        // use a stencil for when the interface is on the
        // left physical boundary MC Eq. 21
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j,k,ncomp) - 23.0_rt*a(i+1,j,k,ncomp) +
                                         13.0_rt*a(i+2,j,k,ncomp) - 3.0_rt*a(i+3,j,k,ncomp));

      } else if (i == domhi[0] && hi_bc[0] != Interior) {
        // use a stencil for the interface that is one zone
        // from the right physical boundary, MC Eq. 22
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j,k,ncomp) + 13.0_rt*a(i-1,j,k,ncomp) -
                                         5.0_rt*a(i-2,j,k,ncomp) + a(i-3,j,k,ncomp));

      } else if (i == domhi[0]+1 && hi_bc[0] != Interior) {
        // use a stencil for when the interface is on the
        // right physical boundary MC Eq. 21
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i-1,j,k,ncomp) - 23.0_rt*a(i-2,j,k,ncomp) +
                                         13.0_rt*a(i-3,j,k,ncomp) - 3.0_rt*a(i-4,j,k,ncomp));

      } else if (i == domlo[0]-1 && lo_bc[0] == Outflow) {
        // extrapolate to the domlo[0]-1 cell using a
        // conservative cubic polynomial averaged over
        // the cell
        Real ac = 4.0_rt*a(domlo[0],j,k,ncomp) - 6.0_rt*a(domlo[0]+1,j,k,ncomp) + 4.0_rt*a(domlo[0]+2,j,k,ncomp) - a(domlo[0]+3,j,k,ncomp);

        // now use the 1-sided stencil from above with
        // this extrapolated value
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*ac - 23.0_rt*a(i+1,j,k,ncomp) +
                                         13.0_rt*a(i+2,j,k,ncomp) - 3.0_rt*a(i+3,j,k,ncomp));

      } else if (i == domhi[0]+2 && hi_bc[0] == Outflow) {
        // extrapolate to the domhi[0]+1 cell using a
        // conservative cubic polynomial averaged over
        // the cell
        Real ac = 4.0_rt*a(domhi[0],j,k,ncomp) - 6.0_rt*a(domhi[0]-1,j,k,ncomp) + 4.0_rt*a(domhi[0]-2,j,k,ncomp) - a(domhi[0]-3,j,k,ncomp);

        // now use the 1-sided stencil from above with
        // this extrapolated value
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*ac - 23.0_rt*a(i-2,j,k,ncomp) +
                                         13.0_rt*a(i-3,j,k,ncomp) - 3.0_rt*a(i-4,j,k,ncomp));

      } else {
        // regular stencil
        a_int(i,j,k) = (7.0_rt/12.0_rt)*(a(i-1,j,k,ncomp) + a(i,j,k,ncomp)) -
                       (1.0_rt/12.0_rt)*(a(i-2,j,k,ncomp) + a(i+1,j,k,ncomp));
      }
    });

  } else if (idir == 1) {

    // this loop is over interfaces

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {

      // interpolate to the edges

      if (j == domlo[1]+1 && lo_bc[1] != Interior) {
        // use a stencil for the interface that is one zone
        // from the left physical boundary, MC Eq. 22
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j-1,k,ncomp) + 13.0_rt*a(i,j,k,ncomp) -
                                         5.0_rt*a(i,j+1,k,ncomp) + a(i,j+2,k,ncomp));

      } else if (j == domlo[1] && lo_bc[1] != Interior) {
        // use a stencil for when the interface is on the
        // left physical boundary MC Eq. 21
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j,k,ncomp) - 23.0_rt*a(i,j+1,k,ncomp) +
                                         13.0_rt*a(i,j+2,k,ncomp) - 3.0_rt*a(i,j+3,k,ncomp));

      } else if (j == domhi[1] && hi_bc[1] != Interior) {
        // use a stencil for the interface that is one zone
        // from the right physical boundary, MC Eq. 22
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j,k,ncomp) + 13.0_rt*a(i,j-1,k,ncomp) -
                                         5.0_rt*a(i,j-2,k,ncomp) + a(i,j-3,k,ncomp));

      } else if (j == domhi[1]+1 && hi_bc[1] != Interior) {
        // use a stencil for when the interface is on the
        // right physical boundary MC Eq. 21
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j-1,k,ncomp) - 23.0_rt*a(i,j-2,k,ncomp) +
                                         13.0_rt*a(i,j-3,k,ncomp) - 3.0_rt*a(i,j-4,k,ncomp));

      } else if (j == domlo[1]-1 && lo_bc[1] == Outflow) {
        // extrapolate to the domlo[1]-1 cell using a
        // conservative cubic polynomial averaged over
        // the cell
        Real ac = 4.0_rt*a(i,domlo[1],k,ncomp) - 6.0_rt*a(i,domlo[1]+1,k,ncomp) + 4.0_rt*a(i,domlo[1]+2,k,ncomp) - a(i,domlo[1]+3,k,ncomp);

        // now use the 1-sided stencil from above with
        // this extrapolated value
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*ac - 23.0_rt*a(i,j+1,k,ncomp) +
                                         13.0_rt*a(i,j+2,k,ncomp) - 3.0_rt*a(i,j+3,k,ncomp));

      } else if (j == domhi[1]+2 && hi_bc[1] == Outflow) {
        // extrapolate to the domhi[1]+1 cell using a
        // conservative cubic polynomial averaged over
        // the cell
        Real ac = 4.0_rt*a(i,domhi[1],k,ncomp) - 6.0_rt*a(i,domhi[1]-1,k,ncomp) + 4.0_rt*a(i,domhi[1]-2,k,ncomp) - a(i,domhi[1]-3,k,ncomp);

        // now use the 1-sided stencil from above with
        // this extrapolated value
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*ac - 23.0_rt*a(i,j-2,k,ncomp) +
                                         13.0_rt*a(i,j-3,k,ncomp) - 3.0_rt*a(i,j-4,k,ncomp));

      } else {
        // regular stencil
        a_int(i,j,k) = (7.0_rt/12.0_rt)*(a(i,j-1,k,ncomp) + a(i,j,k,ncomp)) -
                       (1.0_rt/12.0_rt)*(a(i,j-2,k,ncomp) + a(i,j+1,k,ncomp));
      }

    });

  } else if (idir == 2) {

    // this loop is over interfaces

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {

      // interpolate to the edges

      if (k == domlo[2]+1 && lo_bc[2] != Interior) {
        // use a stencil for the interface that is one zone
        // from the left physical boundary, MC Eq. 22
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j,k-1,ncomp) + 13.0_rt*a(i,j,k,ncomp) -
                                         5.0_rt*a(i,j,k+1,ncomp) + a(i,j,k+2,ncomp));

      } else if (k == domlo[2] && lo_bc[2] != Interior) {
        // use a stencil for when the interface is on the
        // left physical boundary MC Eq. 21
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j,k,ncomp) - 23.0_rt*a(i,j,k+1,ncomp) +
                                         13.0_rt*a(i,j,k+2,ncomp) - 3.0_rt*a(i,j,k+3,ncomp));

      } else if (k == domhi[2] && hi_bc[2] != Interior) {
        // use a stencil for the interface that is one zone
        // from the right physical boundary, MC Eq. 22
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(3.0_rt*a(i,j,k,ncomp) + 13.0_rt*a(i,j,k-1,ncomp) -
                                         5.0_rt*a(i,j,k-2,ncomp) + a(i,j,k-3,ncomp));

      } else if (k == domhi[2]+1 && hi_bc[2] != Interior) {
        // use a stencil for when the interface is on the
        // right physical boundary MC Eq. 21
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*a(i,j,k-1,ncomp) - 23.0_rt*a(i,j,k-2,ncomp) +
                                         13.0_rt*a(i,j,k-3,ncomp) - 3.0_rt*a(i,j,k-4,ncomp));

      } else if (k == domlo[2]-1 && lo_bc[2] == Outflow) {
        // extrapolate to the domlo[2]-1 cell using a
        // conservative cubic polynomial averaged over
        // the cell
        Real ac = 4.0_rt*a(i,j,domlo[2],ncomp) - 6.0_rt*a(i,j,domlo[2]+1,ncomp) + 4.0_rt*a(i,j,domlo[2]+2,ncomp) - a(i,j,domlo[2]+3,ncomp);

        // now use the 1-sided stencil from above with
        // this extrapolated value
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*ac - 23.0_rt*a(i,j,k+1,ncomp) +
                                         13.0_rt*a(i,j,k+2,ncomp) - 3.0_rt*a(i,j,k+3,ncomp));

      } else if (k == domhi[2]+2 && hi_bc[2] == Outflow) {
        // extrapolate to the domhi[2]+1 cell using a
        // conservative cubic polynomial averaged over
        // the cell
        Real ac = 4.0_rt*a(i,j,domhi[2],ncomp) - 6.0_rt*a(i,j,domhi[2]-1,ncomp) + 4.0_rt*a(i,j,domhi[2]-2,ncomp) - a(i,j,domhi[2]-3,ncomp);

        // now use the 1-sided stencil from above with
        // this extrapolated value
        a_int(i,j,k) = (1.0_rt/12.0_rt)*(25.0_rt*ac - 23.0_rt*a(i,j,k-2,ncomp) +
                                         13.0_rt*a(i,j,k-3,ncomp) - 3.0_rt*a(i,j,k-4,ncomp));

      } else {
        // regular stencil
        a_int(i,j,k) = (7.0_rt/12.0_rt)*(a(i,j,k-1,ncomp) + a(i,j,k,ncomp)) -
                       (1.0_rt/12.0_rt)*(a(i,j,k-2,ncomp) + a(i,j,k+1,ncomp));
      }

    });
  }
}


void
Castro::states(const Box& bx,
               const int idir, const int ncomp,
               Array4<Real const> const& a,
               Array4<Real const> const& a_int,
               Array4<Real const> const& flatn,
               Array4<Real> const& al,
               Array4<Real> const& ar) {

  // our convention here is that:
  //     al(i,j,k)   will be al_{i-1/2,j,k),
  //     al(i+1,j,k) will be al_{i+1/2,j,k)
  //
  // Note, this needs to be run on lo-1, hi+1 in all directions

  // we need interface values on all faces of the domain

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  auto domlo = geom.Domain().loVect3d();
  auto domhi = geom.Domain().hiVect3d();

  constexpr Real C2 = 1.25_rt;
  constexpr Real C3 = 0.1_rt;

  if (idir == 0) {

    if (limit_fourth_order == 0) {

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {
        al(i+1,j,k,ncomp) = a_int(i+1,j,k);
        ar(i,j,k,ncomp) = a_int(i,j,k);
      });

    } else {

      // the limiting loops are now over zones

      // this is a loop over cell centers, affecting
      // i-1/2,R and i+1/2,L

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {

        al(i+1,j,k,ncomp) = a_int(i+1,j,k);
        ar(i,j,k,ncomp) = a_int(i,j,k);

        // these live on cell-centers
        Real dafm = a(i,j,k,ncomp) - a_int(i,j,k);
        Real dafp = a_int(i+1,j,k) - a(i,j,k,ncomp);

        // these live on cell-centers
        Real d2af = 6.0_rt*(a_int(i,j,k) - 2.0_rt*a(i,j,k,ncomp) + a_int(i+1,j,k));

        Real d2acm2 = a(i-3,j,k,ncomp) - 2.0_rt*a(i-2,j,k,ncomp) + a(i-1,j,k,ncomp);
        Real d2acm1 = a(i-2,j,k,ncomp) - 2.0_rt*a(i-1,j,k,ncomp) + a(i,j,k,ncomp);
        Real d2ac0  = a(i-1,j,k,ncomp) - 2.0_rt*a(i,j,k,ncomp) + a(i+1,j,k,ncomp);
        Real d2acp1 = a(i,j,k,ncomp) - 2.0_rt*a(i+1,j,k,ncomp) + a(i+2,j,k,ncomp);
        Real d2acp2 = a(i+1,j,k,ncomp) - 2.0_rt*a(i+2,j,k,ncomp) + a(i+3,j,k,ncomp);

        // limit? MC Eq. 24 and 25
        if (dafm * dafp <= 0.0_rt ||
            (a(i,j,k,ncomp) - a(i-2,j,k,ncomp))*(a(i+2,j,k,ncomp) - a(i,j,k,ncomp)) <= 0.0_rt) {

          // we are at an extrema

          Real s = std::copysign(1.0_rt, d2ac0);

          Real d2a_lim;

          if (s == std::copysign(1.0_rt, d2acm1) &&
              s == std::copysign(1.0_rt, d2acp1) &&
              s == std::copysign(1.0_rt, d2af)) {
            // MC Eq. 26
            d2a_lim = s * amrex::min(std::abs(d2af), C2*std::abs(d2acm1),
                                     C2*std::abs(d2ac0), C2*std::abs(d2acp1));
          } else {
            d2a_lim = 0.0_rt;
          }

          Real rho;
          if (std::abs(d2af) <= 1.e-12_rt *
              amrex::max(std::abs(a(i-2,j,k,ncomp)), std::abs(a(i-1,j,k,ncomp)),
                         std::abs(a(i,j,k,ncomp)), std::abs(a(i+1,j,k,ncomp)),
                         std::abs(a(i+2,j,k,ncomp)))) {
            rho = 0.0_rt;
          } else {
            // MC Eq. 27
            rho = d2a_lim/d2af;
          }

          if (rho < 1.0_rt - 1.e-12_rt) {
            // we may need to limit -- these quantities are at cell-centers
            Real d3am1 = d2acm1 - d2acm2;
            Real d3a0  = d2ac0 - d2acm1;
            Real d3ap1 = d2acp1 - d2ac0;
            Real d3ap2 = d2acp2 - d2acp1;

            Real d3a_min = amrex::min(d3am1, d3a0, d3ap1, d3ap2);
            Real d3a_max = amrex::max(d3am1, d3a0, d3ap1, d3ap2);

            if (C3 * amrex::max(std::abs(d3a_min), std::abs(d3a_max)) <=
                (d3a_max - d3a_min)) {
              // limit
              if (dafm*dafp < 0.0_rt) {
                // Eqs. 29, 30
                ar(i,j,k,ncomp) = a(i,j,k,ncomp) - rho*dafm;  // note: typo in Eq 29
                al(i+1,j,k,ncomp) = a(i,j,k,ncomp) + rho*dafp;
              } else if (std::abs(dafm) >= 2.0_rt*std::abs(dafp)) {
                // Eq. 31
                ar(i,j,k,ncomp) = a(i,j,k,ncomp) - 2.0_rt*(1.0_rt - rho)*dafp - rho*dafm;
              } else if (std::abs(dafp) >= 2.0_rt*std::abs(dafm)) {
                // Eq. 32
                al(i+1,j,k,ncomp) = a(i,j,k,ncomp) + 2.0_rt*(1.0_rt - rho)*dafm + rho*dafp;
              }

            }
          }

        } else {
          // if Eqs. 24 or 25 didn't hold we still may need to limit
          if (std::abs(dafm) >= 2.0_rt*std::abs(dafp)) {
            ar(i,j,k,ncomp) = a(i,j,k,ncomp) - 2.0_rt*dafp;
          }
          if (std::abs(dafp) >= 2.0_rt*std::abs(dafm)) {
            al(i+1,j,k,ncomp) = a(i,j,k,ncomp) + 2.0_rt*dafm;
          }
        }

        // apply flattening

        al(i+1,j,k,ncomp) = flatn(i,j,k)*al(i+1,j,k,ncomp) + (1.0_rt - flatn(i,j,k))*a(i,j,k,ncomp);
        ar(i,j,k,ncomp) = flatn(i,j,k)*ar(i,j,k,ncomp) + (1.0_rt - flatn(i,j,k))*a(i,j,k,ncomp);


        // now handle any physical boundaries here by modifying the interface values

        if (i == domlo[0]) {

          // reset the left state at domlo[0] if needed -- it is outside the domain

          if (lo_bc[0] == Outflow) {
            //al(domlo[0],j,k,:) = ar(domlo[0],j,k,:)

          } else if (lo_bc[0] == Symmetry) {
            if (ncomp == QU) {
              al(domlo[0],j,k,QU) = -ar(domlo[0],j,k,QU);
            } else {
              al(domlo[0],j,k,ncomp) = ar(domlo[0],j,k,ncomp);
            }

          } else if (lo_bc[0] == Interior) {
            // we don't need to do anything here

          } else {
            // not supported
            amrex::Error("ERROR: boundary conditions not supported for -X in states");
          }

        }

        if (i == domhi[0]+1) {

          // reset the right state at domhi[0]+1 if needed -- it is outside the domain

          if (hi_bc[0] == Outflow) {
            //ar(domhi[0]+1,j,k,:) = al(domhi[0]+1,j,k,:)

          } else if (hi_bc[0] == Symmetry) {
            if (ncomp == QU) {
              ar(domhi[0]+1,j,k,QU) = -al(domhi[0]+1,j,k,QU);
            } else {
              ar(domhi[0]+1,j,k,ncomp) = al(domhi[0]+1,j,k,ncomp);
            }

          } else if (hi_bc[0] == Interior) {
            // we don't need to do anything here

          } else {
            // not supported
            amrex::Error("ERROR: boundary conditions not supported for +X in states");
          }

        }

      });

    }

  } else if (idir == 1) {

    if (limit_fourth_order == 0) {

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {
        al(i,j+1,k,ncomp) = a_int(i,j+1,k);
        ar(i,j,k,ncomp) = a_int(i,j,k);
      });

    } else {

      // the limiting loops are now over zones

      // this is a loop over cell centers, affecting
      // j-1/2,R and j+1/2,L

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {
        al(i,j+1,k,ncomp) = a_int(i,j+1,k);
        ar(i,j,k,ncomp) = a_int(i,j,k);

        // these live on cell-centers
        Real dafm = a(i,j,k,ncomp) - a_int(i,j,k);
        Real dafp = a_int(i,j+1,k) - a(i,j,k,ncomp);

        // these live on cell-centers
        Real d2af = 6.0_rt*(a_int(i,j,k) - 2.0_rt*a(i,j,k,ncomp) + a_int(i,j+1,k));

        Real d2acm2 = a(i,j-3,k,ncomp) - 2.0_rt*a(i,j-2,k,ncomp) + a(i,j-1,k,ncomp);
        Real d2acm1 = a(i,j-2,k,ncomp) - 2.0_rt*a(i,j-1,k,ncomp) + a(i,j,k,ncomp);
        Real d2ac0  = a(i,j-1,k,ncomp) - 2.0_rt*a(i,j,k,ncomp) + a(i,j+1,k,ncomp);
        Real d2acp1 = a(i,j,k,ncomp) - 2.0_rt*a(i,j+1,k,ncomp) + a(i,j+2,k,ncomp);
        Real d2acp2 = a(i,j+1,k,ncomp) - 2.0_rt*a(i,j+2,k,ncomp) + a(i,j+3,k,ncomp);

        // limit? MC Eq. 24 and 25
        if (dafm * dafp <= 0.0_rt ||
            (a(i,j,k,ncomp) - a(i,j-2,k,ncomp))*(a(i,j+2,k,ncomp) - a(i,j,k,ncomp)) <= 0.0_rt) {

          // we are at an extrema

          Real s = std::copysign(1.0_rt, d2ac0);

          Real d2a_lim;

          if (s == std::copysign(1.0_rt, d2acm1) &&
              s == std::copysign(1.0_rt, d2acp1) &&
              s == std::copysign(1.0_rt, d2af)) {
            // MC Eq. 26
            d2a_lim = s * amrex::min(std::abs(d2af), C2*std::abs(d2acm1),
                                     C2*std::abs(d2ac0), C2*std::abs(d2acp1));
          } else {
            d2a_lim = 0.0_rt;
          }

          Real rho;
          if (std::abs(d2af) <= 1.e-12_rt *
              amrex::max(std::abs(a(i,j-2,k,ncomp)), std::abs(a(i,j-1,k,ncomp)),
                         std::abs(a(i,j,k,ncomp)), std::abs(a(i,j+1,k,ncomp)),
                         std::abs(a(i,j+2,k,ncomp)))) {
            rho = 0.0_rt;
          } else {
            // MC Eq. 27
            rho = d2a_lim/d2af;
          }

          if (rho < 1.0_rt - 1.e-12_rt) {
            // we may need to limit -- these quantities are at cell-centers
            Real d3am1 = d2acm1 - d2acm2;
            Real d3a0  = d2ac0 - d2acm1;
            Real d3ap1 = d2acp1 - d2ac0;
            Real d3ap2 = d2acp2 - d2acp1;

            Real d3a_min = amrex::min(d3am1, d3a0, d3ap1, d3ap2);
            Real d3a_max = amrex::max(d3am1, d3a0, d3ap1, d3ap2);

            if (C3 * amrex::max(std::abs(d3a_min), std::abs(d3a_max)) <=
                (d3a_max - d3a_min)) {
              // limit
              if (dafm*dafp < 0.0_rt) {
                // Eqs. 29, 30
                ar(i,j,k,ncomp) = a(i,j,k,ncomp) - rho*dafm;  // note: typo in Eq 29
                al(i,j+1,k,ncomp) = a(i,j,k,ncomp) + rho*dafp;
              } else if (std::abs(dafm) >= 2.0_rt*std::abs(dafp)) {
                // Eq. 31
                ar(i,j,k,ncomp) = a(i,j,k,ncomp) - 2.0_rt*(1.0_rt - rho)*dafp - rho*dafm;
              } else if (std::abs(dafp) >= 2.0_rt*std::abs(dafm)) {
                // Eq. 32
                al(i,j+1,k,ncomp) = a(i,j,k,ncomp) + 2.0_rt*(1.0_rt - rho)*dafm + rho*dafp;
              }

            }
          }

        } else {
          // if Eqs. 24 or 25 didn't hold we still may need to limit
          if (std::abs(dafm) >= 2.0_rt*std::abs(dafp)) {
            ar(i,j,k,ncomp) = a(i,j,k,ncomp) - 2.0_rt*dafp;
          }
          if (std::abs(dafp) >= 2.0_rt*std::abs(dafm)) {
            al(i,j+1,k,ncomp) = a(i,j,k,ncomp) + 2.0_rt*dafm;
          }
        }

        // apply flattening

        al(i,j+1,k,ncomp) = flatn(i,j,k)*al(i,j+1,k,ncomp) + (1.0_rt - flatn(i,j,k))*a(i,j,k,ncomp);
        ar(i,j,k,ncomp) = flatn(i,j,k)*ar(i,j,k,ncomp) + (1.0_rt - flatn(i,j,k))*a(i,j,k,ncomp);

        // now handle any physical boundaries here by modifying the interface values

        if (j == domlo[1]) {

          // reset the left state at domlo[1] if needed -- it is outside the domain

          if (lo_bc[1] == Outflow) {
            //al(i,domlo[1],k,:) = ar(i,domlo[1],k,:)

          } else if (lo_bc[1] == Symmetry) {
            if (ncomp == QV) {
              al(i,domlo[1],k,QV) = -ar(i,domlo[1],k,QV);
            } else {
              al(i,domlo[1],k,ncomp) = ar(i,domlo[1],k,ncomp);
            }

          } else if (lo_bc[1] == Interior) {
            // we don't need to do anything here

          } else {
            // not supported
            amrex::Error("ERROR: boundary conditions not supported for -Y in states");
          }

        }

        if (j == domhi[1]+1) {

          // reset the right state at domhi[1]+1 if needed -- it is outside the domain

          if (hi_bc[1] == Outflow) {
            //ar(i,domhi[1]+1,k,:) = al(i,domhi[1]+1,k,:)

          } else if (hi_bc[1] == Symmetry) {
            if (ncomp == QV) {
              ar(i,domhi[1]+1,k,QV) = -al(i,domhi[1]+1,k,QV);
            } else {
              ar(i,domhi[1]+1,k,ncomp) = al(i,domhi[1]+1,k,ncomp);
            }

          } else if (hi_bc[1] == Interior) {
            // we don't need to do anything here

          } else {
            // not supported
            amrex::Error("ERROR: boundary conditions not supported for +Y in states");
          }

        }

      });

    }

  } else if (idir == 2) {

    if (limit_fourth_order == 0) {

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {
        al(i,j,k+1,ncomp) = a_int(i,j,k+1);
        ar(i,j,k,ncomp) = a_int(i,j,k);
      });

    } else {

      // this is a loop over cell centers, affecting
      // k-1/2,R and k+1/2,L

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
      {

        al(i,j,k+1,ncomp) = a_int(i,j,k+1);
        ar(i,j,k,ncomp) = a_int(i,j,k);

        // these live on cell-centers
        Real dafm = a(i,j,k,ncomp) - a_int(i,j,k);
        Real dafp = a_int(i,j,k+1) - a(i,j,k,ncomp);

        // these live on cell-centers
        Real d2af = 6.0_rt*(a_int(i,j,k) - 2.0_rt*a(i,j,k,ncomp) + a_int(i,j,k+1));

        Real d2acm2 = a(i,j,k-3,ncomp) - 2.0_rt*a(i,j,k-2,ncomp) + a(i,j,k-1,ncomp);
        Real d2acm1 = a(i,j,k-2,ncomp) - 2.0_rt*a(i,j,k-1,ncomp) + a(i,j,k,ncomp);
        Real d2ac0  = a(i,j,k-1,ncomp) - 2.0_rt*a(i,j,k,ncomp) + a(i,j,k+1,ncomp);
        Real d2acp1 = a(i,j,k,ncomp) - 2.0_rt*a(i,j,k+1,ncomp) + a(i,j,k+2,ncomp);
        Real d2acp2 = a(i,j,k+1,ncomp) - 2.0_rt*a(i,j,k+2,ncomp) + a(i,j,k+3,ncomp);

        // limit? MC Eq. 24 and 25
        if (dafm * dafp <= 0.0_rt ||
            (a(i,j,k,ncomp) - a(i,j,k-2,ncomp))*(a(i,j,k+2,ncomp) - a(i,j,k,ncomp)) <= 0.0_rt) {

          // we are at an extrema

          Real s = std::copysign(1.0_rt, d2ac0);

          Real d2a_lim;

          if (s == std::copysign(1.0_rt, d2acm1) &&
              s == std::copysign(1.0_rt, d2acp1) &&
              s == std::copysign(1.0_rt, d2af)) {
            // MC Eq. 26
            d2a_lim = s * amrex::min(std::abs(d2af), C2*std::abs(d2acm1),
                                     C2*std::abs(d2ac0), C2*std::abs(d2acp1));
          } else {
            d2a_lim = 0.0_rt;
          }

          Real rho;
          if (std::abs(d2af) <= 1.e-12_rt *
              amrex::max(std::abs(a(i,j,k-2,ncomp)), std::abs(a(i,j,k-1,ncomp)),
                         std::abs(a(i,j,k,ncomp)), std::abs(a(i,j,k+1,ncomp)),
                         std::abs(a(i,j,k+2,ncomp)))) {
            rho = 0.0_rt;
          } else {
            // MC Eq. 27
            rho = d2a_lim/d2af;
          }

          if (rho < 1.0_rt - 1.e-12_rt) {
            // we may need to limit -- these quantities are at cell-centers
            Real d3am1 = d2acm1 - d2acm2;
            Real d3a0  = d2ac0 - d2acm1;
            Real d3ap1 = d2acp1 - d2ac0;
            Real d3ap2 = d2acp2 - d2acp1;

            Real d3a_min = amrex::min(d3am1, d3a0, d3ap1, d3ap2);
            Real d3a_max = amrex::max(d3am1, d3a0, d3ap1, d3ap2);

            if (C3 * amrex::max(std::abs(d3a_min), std::abs(d3a_max)) <=
                (d3a_max - d3a_min)) {
              // limit
              if (dafm*dafp < 0.0_rt) {
                // Eqs. 29, 30
                ar(i,j,k,ncomp) = a(i,j,k,ncomp) - rho*dafm;  // note: typo in Eq 29
                al(i,j,k+1,ncomp) = a(i,j,k,ncomp) + rho*dafp;
              } else if (std::abs(dafm) >= 2.0_rt*std::abs(dafp)) {
                // Eq. 31
                ar(i,j,k,ncomp) = a(i,j,k,ncomp) - 2.0_rt*(1.0_rt - rho)*dafp - rho*dafm;
              } else if (std::abs(dafp) >= 2.0_rt*std::abs(dafm)) {
                // Eq. 32
                al(i,j,k+1,ncomp) = a(i,j,k,ncomp) + 2.0_rt*(1.0_rt - rho)*dafm + rho*dafp;
              }

            }
          }

        } else {
          // if Eqs. 24 or 25 didn't hold we still may need to limit
          if (std::abs(dafm) >= 2.0_rt*std::abs(dafp)) {
            ar(i,j,k,ncomp) = a(i,j,k,ncomp) - 2.0_rt*dafp;
          }
          if (std::abs(dafp) >= 2.0_rt*std::abs(dafm)) {
            al(i,j,k+1,ncomp) = a(i,j,k,ncomp) + 2.0_rt*dafm;
          }
        }

        // apply flattening
        al(i,j,k+1,ncomp) = flatn(i,j,k)*al(i,j,k+1,ncomp) + (1.0_rt - flatn(i,j,k))*a(i,j,k,ncomp);
        ar(i,j,k,ncomp) = flatn(i,j,k)*ar(i,j,k,ncomp) + (1.0_rt - flatn(i,j,k))*a(i,j,k,ncomp);


        // now handle any physical boundaries here by modifying the interface values

        if (k == domlo[2]) {
          // reset the left state at domlo[2] if needed -- it is outside the domain

          if (lo_bc[2] == Outflow) {
            //al(i,j,domlo[2],:) = ar(i,j,domlo[2],:)

          } else if (lo_bc[2] == Symmetry) {
            if (ncomp == QW) {
              al(i,j,domlo[2],QW) = -ar(i,j,domlo[2],QW);
            } else {
              al(i,j,domlo[2],ncomp) = ar(i,j,domlo[2],ncomp);
            }

          } else if (lo_bc[2] == Interior) {
            // we don't need to do anything here

          } else {
            // not supported
            amrex::Error("ERROR: boundary conditions not supported for -Z in states");
          }

        }

        if (k == domhi[2]+1) {
          // reset the right state at domhi[2]+1 if needed -- it is outside the domain

          if (hi_bc[2] == Outflow) {
            //ar(i,j,domhi[2]+1,:) = al(i,j,domhi[2]+1,:)

          } else if (hi_bc[2] == Symmetry) {
            if (ncomp == QW) {
              ar(i,j,domhi[2]+1,QW) = -al(i,j,domhi[2]+1,QW);
            } else {
              ar(i,j,domhi[2]+1,ncomp) = al(i,j,domhi[2]+1,ncomp);
            }

          } else if (lo_bc[2] == Interior) {
            // we don't need to do anything here

          } else {
            // not supported
            amrex::Error("ERROR: boundary conditions not supported for +Z in states");
          }

        }

      });

    }

  }

}


void
Castro::fourth_avisc(const Box& bx,
                     Array4<Real const> const& q,
                     Array4<Real const> const& qaux,
                     Array4<Real> const& avis,
                     const int idir) {

  // this computes the *face-centered* artifical viscosity using the
  // 4th order expression from McCorquodale & Colella (Eq. 35)

  constexpr Real beta = 0.3_rt;

  const auto dx = geom.CellSizeArray();

  Real dxinv = 1.0_rt / dx[0];
  Real dyinv = 1.0_rt / dx[1];
  Real dzinv = 1.0_rt / dx[2];

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    Real cmin;

    if (idir == 0) {

      // normal direction

      avis(i,j,k) = (q(i,j,k,QU) - q(i-1,j,k,QU)) * dxinv;
#if AMREX_SPACEDIM >= 2
      avis(i,j,k) += 0.25_rt*(q(i,j+1,k,QV) - q(i,j-1,k,QV) +
                              q(i-1,j+1,k,QV) - q(i-1,j-1,k,QV)) * dyinv;
#endif
#if AMREX_SPACEDIM >= 3
      avis(i,j,k) += 0.25_rt*(q(i,j,k+1,QW) - q(i,j,k-1,QW) +
                              q(i-1,j,k+1,QW) - q(i-1,j,k-1,QW)) * dzinv;
#endif

      cmin = amrex::min(qaux(i,j,k,QC), qaux(i-1,j,k,QC));

    } else if (idir == 1) {

      // normal direction

      avis(i,j,k) = (q(i,j,k,QV) - q(i,j-1,k,QV)) * dyinv;

      avis(i,j,k) += 0.25_rt*(q(i+1,j,k,QU) - q(i-1,j,k,QU) +
                              q(i+1,j-1,k,QU) - q(i-1,j-1,k,QU)) * dxinv;

#if AMREX_SPACEDIM >= 3
      avis(i,j,k) += 0.25_rt*(q(i,j,k+1,QW) - q(i,j,k-1,QW) +
                              q(i,j-1,k+1,QW) - q(i,j-1,k-1,QW)) * dzinv;
#endif

      cmin = amrex::min(qaux(i,j,k,QC), qaux(i,j-1,k,QC));

    } else {

      // normal direction

      avis(i,j,k) = (q(i,j,k,QW) - q(i,j,k-1,QW)) * dzinv;

      avis(i,j,k) += 0.25_rt*(q(i,j+1,k,QV) - q(i,j-1,k,QV) +
                              q(i,j+1,k-1,QV) - q(i,j-1,k-1,QV)) * dyinv;

      avis(i,j,k) += 0.25_rt*(q(i+1,j,k,QU) - q(i-1,j,k,QU) +
                              q(i+1,j,k-1,QU) - q(i-1,j,k-1,QU)) * dxinv;

      cmin = amrex::min(qaux(i,j,k,QC), qaux(i,j,k-1,QC));

    }

    // MC Eq. 36

    Real coeff = amrex::min(1.0_rt, (dx[idir] * avis(i,j,k)) * (dx[idir] * avis(i,j,k)) /
                                    (beta * cmin * cmin));

    if (avis(i,j,k) < 0.0_rt) {
      avis(i,j,k) = dx[idir] * avis(i,j,k) * coeff;
    } else {
      avis(i,j,k) = 0.0_rt;
    }

  });

}


#ifdef DIFFUSION
  subroutine add_diffusive_flux(lo, hi, &
                                q, q_lo, q_hi, ncomp, temp_comp, &
                                qint, qi_lo, qi_hi, &
                                F, F_lo, F_hi, &
                                dx, idir, is_avg) bind(C, name="add_diffusive_flux")

    ! add the diffusive flux to the energy fluxes
    !

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NQ, NVAR, &
                                   URHO, &
                                   UEDEN, UEINT, UTEMP, &
                                   QRHO, QREINT, QFS
    use eos_type_module, only : eos_t, eos_input_re
    use conductivity_module, only : conducteos
    use network, only : nspec
    use castro_error_module, only : castro_error

    integer, intent(in), value :: idir
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qi_lo(3), qi_hi(3)
    integer, intent(in) :: F_lo(3), F_hi(3)
    integer, intent(in), value :: ncomp, temp_comp

    real(rt), intent(in) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), ncomp)
    real(rt), intent(in) :: qint(qi_lo(1):qi_hi(1), qi_lo(2):qi_hi(2), qi_lo(3):qi_hi(3), NQ)
    real(rt), intent(out) :: F(F_lo(1):F_hi(1), F_lo(2):F_hi(2), F_lo(3):F_hi(3), NVAR)
    integer, intent(in) :: lo(3), hi(3)
    real(rt), intent(in) :: dx(3)
    integer, intent(in), value :: is_avg

    integer :: i, j, k

    type(eos_t) :: eos_state
    real(rt) :: dTdx

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state % rho = qint(i,j,k,QRHO)
             eos_state % T = q(i,j,k,temp_comp)   ! initial guess
             eos_state % e = qint(i,j,k,QREINT) / qint(i,j,k,QRHO)
             eos_state % xn(:) = qint(i,j,k,QFS:QFS-1+nspec)

             call conducteos(eos_input_re, eos_state)

             if (idir == 1) then

                if (is_avg == 0) then
                   ! we are working with the cell-center state
                   dTdx = (-q(i+1,j,k,temp_comp) + 27*q(i,j,k,temp_comp) - &
                        27*q(i-1,j,k,temp_comp) + q(i-2,j,k,temp_comp))/(24.0_rt * dx(1))

                else
                   ! we are working with the cell-average state
                   dTdx = (-q(i+1,j,k,temp_comp) + 15*q(i,j,k,temp_comp) - &
                        15*q(i-1,j,k,temp_comp) + q(i-2,j,k,temp_comp))/(12.0_rt * dx(1))
                end if

             else if (idir == 2) then

                if (is_avg == 0) then
                   ! we are working with the cell-center state
                   dTdx = (-q(i,j+1,k,temp_comp) + 27*q(i,j,k,temp_comp) - &
                        27*q(i,j-1,k,temp_comp) + q(i,j-2,k,temp_comp))/(24.0_rt * dx(2))

                else
                   ! we are working with the cell-average state
                   dTdx = (-q(i,j+1,k,temp_comp) + 15*q(i,j,k,temp_comp) - &
                        15*q(i,j-1,k,temp_comp) + q(i,j-2,k,temp_comp))/(12.0_rt * dx(2))
                end if

             else

                if (is_avg == 0) then
                   ! we are working with the cell-center state
                   dTdx = (-q(i,j,k+1,temp_comp) + 27*q(i,j,k,temp_comp) - &
                        27*q(i,j,k-1,temp_comp) + q(i,j,k-2,temp_comp))/(24.0_rt * dx(3))

                else
                   ! we are working with the cell-average state
                   dTdx = (-q(i,j,k+1,temp_comp) + 15*q(i,j,k,temp_comp) - &
                        15*q(i,j,k-1,temp_comp) + q(i,j,k-2,temp_comp))/(12.0_rt * dx(3))
                end if

             endif

             F(i,j,k,UEINT) = F(i,j,k,UEINT) - eos_state % conductivity * dTdx
             F(i,j,k,UEDEN) = F(i,j,k,UEDEN) - eos_state % conductivity * dTdx

          end do
       end do
    end do

  end subroutine add_diffusive_flux
