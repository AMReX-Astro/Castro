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
  int* domain_hi = geom.Domain().hiVect();

  const Real *dx = geom.CellSize();

  const Real* prob_lo = geom.ProbLo();

  amrex::Vector<Real> offset = {prob_lo[0]. prob_lo[1], prob_lo[2]};

  amrex::Vector<int> idx = {i, j, k};

  for (int idir = 0; idir < 3; idir++) {
    if ((idir == 0 && ccx) || (idir == 1 && ccy) || (idir == 2 && ccz)) {

      // If we're cell-centered, we want to be in the middle of the zone.
      offset[dir] = offset[dir] + 0.5 * dx[dir];

    } else {
      // Take care of the fact that for edge-centered indexing,
      // we actually range from (domlo, domhi+1).
      domain_hi[idir] += 1;
    }
  }

  // Be careful when using periodic boundary conditions. In that case,
  // we need to loop around to the other side of the domain.

  for {int idir = 0; idir < 3; idir++) {
    if (physbc_lo[idir] == Interior && idx(dir) < domain_lo[idir]) {
      offset[dir] += (probhi[dir] - problo[dir]);
    } else if (physbc_hi[dir] == Interior && idx[dir] > domain_hi[idir]) {
      offset(dir) = offset(dir) + (problo(dir) - probhi(dir));
    }

    loc[idir] = offset[idir] + static_cast<double> (idx[idir]) * dx[idir];
  }

  return loc;

}
