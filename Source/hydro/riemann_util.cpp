#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"
#include "Castro_util.H"

#ifdef RADIATION
#include "Radiation.H"
#include "fluxlimiter.H"
#endif

#ifdef HYBRID_MOMENTUM
#include "hybrid.H"
#endif

#include <eos.H>

#include <cmath>

using namespace amrex;

#include <riemann.H>

void
Castro::store_godunov_state(const Box& bx,
                            Array4<Real const> const qint,
#ifdef RADIATION
                            Array4<Real const> const lambda,
#endif
                            Array4<Real> const qgdnv) {

  // this copies the full interface state (NQ -- one for each primitive
  // variable) over to a smaller subset of size NGDNV for use later in the
  // hydro advancement.

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    qgdnv(i,j,k,GDRHO) = qint(i,j,k,QRHO);
    qgdnv(i,j,k,GDU) = qint(i,j,k,QU);
    qgdnv(i,j,k,GDV) = qint(i,j,k,QV);
    qgdnv(i,j,k,GDW) = qint(i,j,k,QW);
    qgdnv(i,j,k,GDPRES) = qint(i,j,k,QPRES);
#ifdef RADIATION
    for (int g = 0; g < NGROUPS; g++) {
      qgdnv(i,j,k,GDLAMS+g) = lambda(i,j,k,g);
      qgdnv(i,j,k,GDERADS+g) = qint(i,j,k,QRAD+g);
    }
#endif
  });
}
