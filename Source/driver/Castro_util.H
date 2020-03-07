#ifndef _Castro_util_H_
#define _Castro_util_H_

#include "AMReX_Geometry.H"

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void position(int i, int j, int k,
              GeometryData const& geomdata, Real* loc,
              bool ccx = true, bool ccy = true, bool ccz = true)
{
    // Given 3D indices (i,j,k), return the cell-centered spatial position.
    // Optionally we can also be edge-centered in any of the directions.

    int idx[3] = {i, j, k};

    bool cc[3] = {ccx, ccy, ccz};

    Real offset[AMREX_SPACEDIM];

    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbLo();
    const Real* dx     = geomdata.CellSize();
    const int*  domlo  = geomdata.domain.loVect();
    const int*  domhi  = geomdata.domain.hiVect();

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (cc[dir]) {
            // If we're cell-centered, we want to be in the middle of the zone.
            offset[dir] = problo[dir] + 0.5_rt * dx[dir];
        }
        else {
            // We are edge-centered.
            offset[dir] = problo[dir];
        }
    }

    // Be careful when using periodic boundary conditions. In that case,
    // we need to loop around to the other side of the domain.

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {

        if (geomdata.isPeriodic(dir)) {
            if (idx[dir] < domlo[dir]) {
                offset[dir] += probhi[dir] - problo[dir];
            }

            // For the hi end we need to deal correctly with
            // edge-centered indexing.
            if ((cc[dir] && idx[dir] > domhi[dir]) ||
                (!cc[dir] && idx[dir] > domhi[dir]+1)) {
                offset[dir] += problo[dir] - probhi[dir];
            }
        }

    }

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        loc[dir] = offset[dir] + idx[dir] * dx[dir];
    }

    for (int dir = AMREX_SPACEDIM; dir < 3; ++dir) {
        loc[dir] = 0.0_rt;
    }

}

#endif
