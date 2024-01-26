#include <limits.h>
#include <math.h>
#include <AMReX_LO_BCTYPES.H>
#include <RadInterpBndryData.H>

using namespace amrex;

#if (AMREX_SPACEDIM == 2)
#define NUMDERIV 2
#endif

#if (AMREX_SPACEDIM == 3)
#define NUMDERIV 5
#endif

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
REAL* fabdat = (fab).dataPtr();
#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const REAL* fabdat = (fab).dataPtr();

RadInterpBndryData::RadInterpBndryData(const BoxArray& _grids, const DistributionMapping& _dmap,
                                       int _ncomp, const Geometry& geom)
    : RadBndryData(_grids,_dmap,_ncomp,geom)
{
}

void
RadInterpBndryData::setBndryValues(Real bv)
{
  for (OrientationIter fi; fi; ++fi) {
      bndry[fi()].setVal(bv);
  }
}

// At the coarsest level the bndry values are taken from adjacent grids.
void
RadInterpBndryData::setBndryValues(const MultiFab& mf, int mf_start,
                                int bnd_start, int num_comp,
                                const BCRec& bc )
{
      // check that boxarrays are identical
    BL_ASSERT( grids.size() );
    BL_ASSERT( grids == mf.boxArray() );

      // set bndry flags and locations
    IntVect ref_ratio = IntVect::TheUnitVector();
    setBndryConds(bc, geom, ref_ratio);

    //for (int grd = 0; grd < ngrd; grd++) {
    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
        BL_ASSERT(grids[mfi.index()] == mfi.validbox());
        int grd = mfi.index();
        const Box& bx = grids[grd];
        for (OrientationIter fi; fi; ++fi) {
            Orientation face(fi());
            if (bx[face] == geom.Domain()[face]) {
                  // physical bndry, copy from grid
                FArrayBox& bnd_fab = bndry[face][mfi];
                bnd_fab.copy<RunOn::Host>(mf[mfi],mf_start,bnd_start,num_comp);
            }
        }
    }

    // now copy boundary values stored in ghost cells of fine
    // into bndry.  This only does something for physical boundaries,
    // we don't need to make it periodic aware
    for (OrientationIter fi; fi; ++fi) {
        bndry[fi()].copyFrom(mf,0,mf_start,bnd_start,num_comp);
    }
}

// (1) set bndry type and location of bndry value on each face of
//     each grid
// (2) set actual bndry value by:
//     (A) Interpolate from crse bndryRegister at crse/fine interface
//     (B) Copy from ghost region of MultiFab at physical bndry
//     (C) Copy from valid region of MultiFab at fine/fine interface
void
RadInterpBndryData::setBndryValues(BndryRegister& crse, int c_start,
                                const MultiFab& fine, int f_start,
                                int bnd_start, int num_comp, IntVect& ratio,
                                const BCRec& bc)
{
      // check that boxarrays are identical
    BL_ASSERT( grids.size() );
    BL_ASSERT( grids == fine.boxArray() );

      // set bndry types and bclocs
    setBndryConds(bc, geom, ratio);

      // first interpolate from coarse to fine on bndry
    const Box& fine_domain = geom.Domain();
      // mask turned off if covered by fine grid
    Real *derives = 0;
    int  tmplen = 0;
    //for (int grd = 0; grd < ngrd; grd++) {
    for(MFIter mfi(fine); mfi.isValid(); ++mfi) {
        BL_ASSERT(grids[mfi.index()] == mfi.validbox());
        int grd = mfi.index();
        const Box& fine_bx = grids[grd];
        Box crse_bx = amrex::coarsen(fine_bx,ratio);
        auto cblo = crse_bx.loVect3d();
        auto cbhi = crse_bx.hiVect3d();
        int mxlen = crse_bx.longside() + 2;
        if (std::pow(mxlen,(float)AMREX_SPACEDIM-1) > tmplen) {
            delete [] derives;
#if (AMREX_SPACEDIM == 1)
            derives = new Real[1];
#else
            tmplen = mxlen;
#  if (AMREX_SPACEDIM > 2)
            tmplen *= mxlen;
#  endif
            derives = new Real[tmplen*NUMDERIV];
#endif
        }
        auto lo = fine_bx.loVect3d();
        auto hi = fine_bx.hiVect3d();
        const FArrayBox& fine_grd = fine[mfi];

        for (OrientationIter fi; fi; ++fi) {
            Orientation face(fi());
            int dir = face.coordDir();
            if (fine_bx[face] != fine_domain[face] ||
                geom.isPeriodic(dir)) {
                  // internal or periodic edge, interpolate from crse data
                const Mask& mask = *masks[face][grd];
                Array4<int const> const mask_arr = mask.array();

                const FArrayBox& crse_fab = crse[face][mfi];
                Array4<Real const> const crse_arr = crse_fab.array(c_start);

                int clo[3], chi[3], flo[3], fhi[3];

                for (int d = 0; d < 3; ++d) {
                    clo[d] = cblo[d];
                    chi[d] = cbhi[d];
                    flo[d] = lo[d];
                    fhi[d] = hi[d];

                    // For the face we're operating on, we want to make
                    // sure that we're only operating on the ghost zone
                    // immediately adjacent to the boundary.

                    if (d == dir) {
                        if (face.isLow()) {
                            clo[d] -= 1;
                            flo[d] -= 1;
                            chi[d] = clo[d];
                            fhi[d] = flo[d];
                        }
                        else {
                            chi[d] += 1;
                            fhi[d] += 1;
                            clo[d] = chi[d];
                            flo[d] = fhi[d];
                        }
                    }
                }

                int ratiox = 1;
                int ratioy = 1;
                int ratioz = 1;

                if (dir != 0) {
                    ratiox = ratio[0];
                }

#if AMREX_SPACEDIM >= 2
                if (dir != 1) {
                    ratioy = ratio[1];
                }
#endif

#if AMREX_SPACEDIM == 3
                if (dir != 2) {
                    ratioz = ratio[2];
                }
#endif

                FArrayBox& bnd_fab = bndry[face][mfi];
                Array4<Real> const bdry = bnd_fab.array(bnd_start);

                int is_not_covered = RadBndryData::not_covered;

                for (int n = 0; n < num_comp; ++n) {

                    // Note that only two of these three loops will do something
                    // nontrivial, depending on which face we are working on.

                    for (int kc = clo[2]; kc <= chi[2]; ++kc) {
                        int k = (dir == 2) ? flo[2] : ratioz * kc;

                        for (int jc = clo[1]; jc <= chi[1]; ++jc) {
                            int j = (dir == 1) ? flo[1] : ratioy * jc;

                            for (int ic = clo[0]; ic <= chi[0]; ++ic) {
                                int i = (dir == 0) ? flo[0] : ratiox * ic;

                                Real dcdx = 0.0_rt;
                                Real dcdy = 0.0_rt;
                                Real dcdz = 0.0_rt;

                                Real dcdx2 = 0.0_rt;
                                Real dcdy2 = 0.0_rt;
                                Real dcdz2 = 0.0_rt;

                                Real dcdxy = 0.0_rt;
                                Real dcdxz = 0.0_rt;
                                Real dcdyz = 0.0_rt;

                                if (dir != 0) {
                                    dcdx = 0.5_rt * (crse_arr(ic+1,jc,kc,n) - crse_arr(ic-1,jc,kc,n));
                                    dcdx2 = 0.5_rt * (crse_arr(ic+1,jc,kc,n) - 2.0_rt * crse_arr(ic,jc,kc,n) + crse_arr(ic-1,jc,kc,n));
                                }

#if AMREX_SPACEDIM >= 2
                                if (dir != 1) {
                                    dcdy = 0.5_rt * (crse_arr(ic,jc+1,kc,n) - crse_arr(ic,jc-1,kc,n));
                                    dcdy2 = 0.5_rt * (crse_arr(ic,jc+1,kc,n) - 2.0_rt * crse_arr(ic,jc,kc,n) + crse_arr(ic,jc-1,kc,n));
                                }
#endif

#if AMREX_SPACEDIM == 3
                                if (dir != 2) {
                                    dcdz = 0.5_rt * (crse_arr(ic,jc,kc+1,n) - crse_arr(ic,jc,kc-1,n));
                                    dcdz2 = 0.5_rt * (crse_arr(ic,jc,kc+1,n) - 2.0_rt * crse_arr(ic,jc,kc,n) + crse_arr(ic,jc,kc-1,n));
                                }
#endif

#if AMREX_SPACEDIM == 3
                                if (dir == 0) {
                                    dcdyz = 0.25_rt * (crse_arr(ic,jc+1,kc+1,n) - crse_arr(ic,jc-1,kc+1,n) +
                                                       crse_arr(ic,jc-1,kc-1,n) - crse_arr(ic,jc+1,kc-1,n));
                                }

                                if (dir == 1) {
                                    dcdxz = 0.25_rt * (crse_arr(ic+1,jc,kc+1,n) - crse_arr(ic-1,jc,kc+1,n) +
                                                       crse_arr(ic-1,jc,kc-1,n) - crse_arr(ic+1,jc,kc-1,n));
                                }

                                if (dir == 2) {
                                    dcdxy = 0.25_rt * (crse_arr(ic+1,jc+1,kc,n) - crse_arr(ic-1,jc+1,kc,n) +
                                                       crse_arr(ic-1,jc-1,kc,n) - crse_arr(ic+1,jc-1,kc,n));
                                }
#endif

                                if (dir != 0) {
                                    if (mask_arr(i-1,j,k) != is_not_covered) {
                                        dcdx = crse_arr(ic+1,jc,kc,n) - crse_arr(ic,jc,kc,n);
                                        dcdx2 = 0.0_rt;
                                    }
                                    if (mask_arr(i+ratiox,j,k) != is_not_covered) {
                                        dcdx = crse_arr(ic,jc,kc,n) - crse_arr(ic-1,jc,kc,n);
                                        dcdx2 = 0.0_rt;
                                    }
                                    if (mask_arr(i-1,j,k) != is_not_covered && mask_arr(i+ratiox,j,k) != is_not_covered) {
                                        dcdx = 0.0_rt;
                                    }
                                }

#if AMREX_SPACEDIM >= 2
                                if (dir != 1) {
                                    if (mask_arr(i,j-1,k) != is_not_covered) {
                                        dcdy = crse_arr(ic,jc+1,kc,n) - crse_arr(ic,jc,kc,n);
                                        dcdy2 = 0.0_rt;
                                    }
                                    if (mask_arr(i,j+ratioy,k) != is_not_covered) {
                                        dcdy = crse_arr(ic,jc,kc,n) - crse_arr(ic,jc-1,kc,n);
                                        dcdy2 = 0.0_rt;
                                    }
                                    if (mask_arr(i,j-1,k) != is_not_covered && mask_arr(i,j+ratioy,k) != is_not_covered) {
                                        dcdy = 0.0_rt;
                                    }
                                }
#endif

#if AMREX_SPACEDIM == 3
                                if (dir != 2) {
                                    if (mask_arr(i,j,k-1) != is_not_covered) {
                                        dcdz = crse_arr(ic,jc,kc+1,n) - crse_arr(ic,jc,kc,n);
                                        dcdz2 = 0.0_rt;
                                    }
                                    if (mask_arr(i,j,k+ratioz) != is_not_covered) {
                                        dcdz = crse_arr(ic,jc,kc,n) - crse_arr(ic,jc,kc-1,n);
                                        dcdz2 = 0.0_rt;
                                    }
                                    if (mask_arr(i,j,k-1) != is_not_covered && mask_arr(i,j,k+ratioz) != is_not_covered) {
                                        dcdz = 0.0_rt;
                                    }
                                }
#endif

#if AMREX_SPACEDIM == 3
                                if (dir == 0) {
                                    if ((mask_arr(i,j+ratioy,k+ratioz) != is_not_covered) ||
                                        (mask_arr(i,j-1     ,k+ratioz) != is_not_covered) ||
                                        (mask_arr(i,j+ratioy,k-1     ) != is_not_covered) ||
                                        (mask_arr(i,j-1     ,k-1     ) != is_not_covered)) {
                                        dcdyz = 0.0_rt;
                                    }
                                }

                                if (dir == 1) {
                                    if ((mask_arr(i+ratiox,j,k+ratioz) != is_not_covered) ||
                                        (mask_arr(i-1     ,j,k+ratioz) != is_not_covered) ||
                                        (mask_arr(i+ratiox,j,k-1     ) != is_not_covered) ||
                                        (mask_arr(i-1     ,j,k-1     ) != is_not_covered)) {
                                        dcdxz = 0.0_rt;
                                    }
                                }

                                if (dir == 2) {
                                    if ((mask_arr(i+ratiox,j+ratioy,k) != is_not_covered) ||
                                        (mask_arr(i-1     ,j+ratioy,k) != is_not_covered) ||
                                        (mask_arr(i+ratiox,j-1     ,k) != is_not_covered) ||
                                        (mask_arr(i-1     ,j-1     ,k) != is_not_covered)) {
                                        dcdxy = 0.0_rt;
                                    }
                                }
#endif

                                for (int koff = 0; koff < ratioz; ++koff) {
                                    Real zz = (koff - 0.5_rt * ratioz + 0.5_rt) / ratioz;
                                    int kk = (dir == 2) ? flo[2] : ratioz * kc + koff;

                                    for (int joff = 0; joff < ratioy; ++joff) {
                                        Real yy = (joff - 0.5_rt * ratioy + 0.5_rt) / ratioy;
                                        int jj = (dir == 1) ? flo[1] : ratioy * jc + joff;

                                        for (int ioff = 0; ioff < ratiox; ++ioff) {
                                            Real xx = (ioff - 0.5_rt * ratiox + 0.5_rt) / ratiox;
                                            int ii = (dir == 0) ? flo[0] : ratiox * ic + ioff;

                                            bdry(ii,jj,kk,n) = crse_arr(ic,jc,kc,n) +
                                                               xx * dcdx + xx * xx * dcdx2 +
                                                               yy * dcdy + yy * yy * dcdy2 +
                                                               zz * dcdz + zz * zz * dcdz2 +
                                                               xx * yy * dcdxy + xx * zz * dcdxz + yy * zz * dcdyz;
                                        }
                                    }
                                }
                            }
                        }
                    }

                } // num_comp
            } else {
                  // physical bndry, copy from ghost region of
                  // corresponding grid
                FArrayBox& bnd_fab = bndry[face][mfi];
                bnd_fab.copy<RunOn::Host>(fine_grd,f_start,bnd_start,num_comp);
            }
        }
    }
    delete[] derives;

    // now copy boundary values stored in ghost cells of fine
    // into bndry.  This only does something for physical boundaries,
    // we don't need to make it periodic aware
    for (OrientationIter face; face; ++face) {
        bndry[face()].copyFrom(fine,0,f_start,bnd_start,num_comp);
    }
}

void RadInterpBndryData::setBndryConds(const BCRec& phys_bc,
                                    const Geometry& geom, int ratio)
{

    IntVect ratio_vect = ratio * IntVect::TheUnitVector();
    setBndryConds(phys_bc, geom, ratio_vect);

}

