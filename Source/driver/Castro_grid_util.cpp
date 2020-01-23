#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

amrex::Vector<Real>
Castro::position_c(const int i, const int j, const int k,
                   const bool ccx, const bool ccy, const bool ccz) {

  amrex::Vector<Real> loc;
  loc.resize(3);

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

  const Real *dx = geom.CellSize();

  const Real* prob_lo = geom.ProbLo();
  const Real* prob_hi = geom.ProbHi();

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  amrex::Vector<int> idx = {i, j, k};

  for (int idir = 0; idir < 3; idir++) {
    int correction = 0;
    Real offset = prob_lo[idir];

    if ((idir == 0 && ccx) || (idir == 1 && ccy) || (idir == 2 && ccz)) {

      // If we're cell-centered, we want to be in the middle of the zone.
      offset += 0.5 * dx[idir];

    } else {
      // Take care of the fact that for edge-centered indexing,
      // we actually range from (domlo, domhi+1).
      correction = 1;
    }

    // Be careful when using periodic boundary conditions. In that case,
    // we need to loop around to the other side of the domain.

    if (lo_bc[idir] == Interior && idx[idir] < domain_lo[idir]) {
      offset += (prob_hi[idir] - prob_lo[idir]);

    } else if (hi_bc[idir] == Interior && idx[idir] > domain_hi[idir] + correction) {
      offset += (prob_lo[idir] - prob_hi[idir]);
    }

    loc[idir] = offset + static_cast<double> (idx[idir]) * dx[idir];
  }

  return loc;

}
