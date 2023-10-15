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
        const int* cblo = crse_bx.loVect();
        const int* cbhi = crse_bx.hiVect();
        int mxlen = crse_bx.longside() + 2;
        if (pow(mxlen,(float)AMREX_SPACEDIM-1) > tmplen) {
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
        const int* lo = fine_bx.loVect();
        const int* hi = fine_bx.hiVect();
        const FArrayBox& fine_grd = fine[mfi];

        for (OrientationIter fi; fi; ++fi) {
            Orientation face(fi());
            int dir = face.coordDir();
            if (fine_bx[face] != fine_domain[face] ||
                geom.isPeriodic(dir)) {
                  // internal or periodic edge, interpolate from crse data
                const Mask& mask = *masks[face][grd];
                const int* mlo = mask.loVect();
                const int* mhi = mask.hiVect();
                Array4<int const> const mask_arr = mask.array();

                const FArrayBox& crse_fab = crse[face][mfi];
                const int* clo = crse_fab.loVect();
                const int* chi = crse_fab.hiVect();
                Array4<Real const> const crse = crse_fab.array(c_start);

                int iclo = 0;
                int ichi = 0;
                int jclo = 0;
                int jchi = 0;
                int kclo = 0;
                int kchi = 0;

                if (face.coordDir() != 0) {
                    iclo = cblo[0];
                    ichi = cbhi[1];
                }

#if AMREX_SPACEDIM >= 2
                if (face.coordDir() != 1) {
                    jclo = cblo[1];
                    jchi = cbhi[1];
                }
#endif

#if AMREX_SPACEDIM == 3
                if (face.coordDir() != 2) {
                    kclo = cblo[2];
                    kchi = cbhi[2];
                }
#endif

                int ratiox = 1;
                int ratioy = 1;
                int ratioz = 1;

                if (face.coordDir() != 0) {
                    ratiox = ratio[0];
                }

#if AMREX_SPACEDIM >= 2
                if (face.coordDir() != 1) {
                    ratioy = ratio[1];
                }
#endif

#if AMREX_SPACEDIM == 3
                if (face.coordDir() != 2) {
                    ratioz = ratio[2];
                }
#endif

                FArrayBox& bnd_fab = bndry[face][mfi];
                const int* blo = bnd_fab.loVect();
                const int* bhi = bnd_fab.hiVect();
                Array4<Real> const bdry = bnd_fab.array(bnd_start);

                int is_not_covered = RadBndryData::not_covered;

                for (int n = 0; n < num_comp; ++n) {

                    // Note that only two of these three loops will do something
                    // nontrivial, depending on which face we are working on.

                    for (int koff = 0; koff < ratioz; ++koff) {
                        Real zz = (koff - 0.5_rt * ratioz + 0.5_rt) / ratioz;

                        for (int kc = kclo; kc <= kchi; ++kc) {
                            int k = ratioz * kc + koff;

                            for (int joff = 0; joff < ratioy; ++joff) {
                                Real yy = (joff - 0.5_rt * ratioy + 0.5_rt) / ratioy;

                                for (int jc = jclo; jc <= jchi; ++jc) {
                                    int j = ratioy * jc + joff;

                                    for (int ioff = 0; ioff < ratiox; ++ioff) {
                                        Real xx = (ioff - 0.5_rt * ratiox + 0.5_rt) / ratiox;

                                        for (int ic = iclo; ic <= ichi; ++ic) {
                                            int i = ratiox * ic + ioff;

                                            Real xderiv = 0.0_rt;
                                            Real yderiv = 0.0_rt;
                                            Real zderiv = 0.0_rt;

                                            Real xxderiv = 0.0_rt;
                                            Real yyderiv = 0.0_rt;
                                            Real zzderiv = 0.0_rt;

                                            Real xyderiv = 0.0_rt;
                                            Real xzderiv = 0.0_rt;
                                            Real yzderiv = 0.0_rt;

                                            if (face.coordDir() != 0) {
                                                xderiv = 0.5_rt * (crse(ic+1,jc,kc,n) - crse(ic-1,jc,kc,n));
                                                xxderiv = 0.5_rt * (crse(ic+1,jc,kc,n) - 2.0_rt * crse(ic,jc,kc,n) + crse(ic-1,jc,kc,n));
                                            }

#if AMREX_SPACEDIM >= 2
                                            if (face.coordDir() != 1) {
                                                yderiv = 0.5_rt * (crse(ic,jc+1,kc,n) - crse(ic,jc-1,kc,n));
                                                yyderiv = 0.5_rt * (crse(ic,jc+1,kc,n) - 2.0_rt * crse(ic,jc,kc,n) + crse(ic,jc-1,kc,n));
                                            }
#endif

#if AMREX_SPACEDIM == 3
                                            if (face.coordDir() != 2) {
                                                zderiv = 0.5_rt * (crse(ic,jc,kc+1,n) - crse(ic,jc,kc-1,n));
                                                zzderiv = 0.5_rt * (crse(ic,jc,kc+1,n) - 2.0_rt * crse(ic,jc,kc,n) + crse(ic,jc,kc-1,n));
                                            }
#endif

#if AMREX_SPACEDIM == 3
                                            if (face.coordDir() == 0) {
                                                yzderiv = 0.25_rt * (crse(ic,jc+1,kc+1,n) - crse(ic,jc-1,kc+1,n) +
                                                                     crse(ic,jc-1,kc-1,n) - crse(ic,jc+1,kc-1,n));
                                            }

                                            if (face.coordDir() == 1) {
                                                xzderiv = 0.25_rt * (crse(ic+1,jc,kc+1,n) - crse(ic-1,jc,kc+1,n) +
                                                                     crse(ic-1,jc,kc-1,n) - crse(ic+1,jc,kc-1,n));
                                            }

                                            if (face.coordDir() == 2) {
                                                xyderiv = 0.25_rt * (crse(ic+1,jc+1,kc,n) - crse(ic-1,jc+1,kc,n) +
                                                                     crse(ic-1,jc-1,kc,n) - crse(ic+1,jc-1,kc,n));
                                            }
#endif

                                            if (mask_arr(i-1,j,k) != not_covered) {
                                                xderiv = crse(ic+1,jc,kc,n) - crse(ic,jc,kc,n);
                                                xxderiv = 0.0_rt;
                                            }
                                            if (mask_arr(i+ratiox,j,k) != not_covered) {
                                                xderiv = crse(ic,jc,kc,n) - crse(ic-1,jc,kc,n);
                                                xxderiv = 0.0_rt;
                                            }
                                            if (mask_arr(i-1,j,k) != not_covered && mask_arr(i+ratiox,j,k) != not_covered) {
                                                xderiv = 0.0_rt;
                                            }

#if AMREX_SPACEDIM >= 2
                                            if (mask_arr(i,j-1,k) != not_covered) {
                                                yderiv = crse(ic,jc+1,kc,n) - crse(ic,jc,kc,n);
                                                yyderiv = 0.0_rt;
                                            }
                                            if (mask_arr(i,j+ratioy,k) != not_covered) {
                                                yderiv = crse(ic,jc,kc,n) - crse(ic,jc-1,kc,n);
                                                yyderiv = 0.0_rt;
                                            }
                                            if (mask_arr(i,j-1,k) != not_covered && mask_arr(i,j+ratioy,k) != not_covered) {
                                                yderiv = 0.0_rt;
                                            }
#endif

#if AMREX_SPACEDIM == 3
                                            if (mask_arr(i,j,k-1) != not_covered) {
                                                yderiv = crse(ic,jc,kc+1,n) - crse(ic,jc,kc,n);
                                                yyderiv = 0.0_rt;
                                            }
                                            if (mask_arr(i,j,k+ratioz) != not_covered) {
                                                yderiv = crse(ic,jc,kc,n) - crse(ic,jc,kc-1,n);
                                                yyderiv = 0.0_rt;
                                            }
                                            if (mask_arr(i,j,k-1) != not_covered && mask_arr(i,j,k+ratioz) != not_covered) {
                                                yderiv = 0.0_rt;
                                            }
#endif

#if AMREX_SPACEDIM == 3
                                            if ((mask_arr(i,j+ratioy,k+ratioz) != not_covered) ||
                                                (mask_arr(i,j-1     ,k+ratioz) != not_covered) ||
                                                (mask_arr(i,j+ratioy,k-1     ) != not_covered) ||
                                                (mask_arr(i,j-1     ,k-1     ) != not_covered)) {
                                                yzderiv = 0.0_rt;
                                            }

                                            if ((mask_arr(i+ratiox,j,k+ratioz) != not_covered) ||
                                                (mask_arr(i-1     ,j,k+ratioz) != not_covered) ||
                                                (mask_arr(i+ratiox,j,k-1     ) != not_covered) ||
                                                (mask_arr(i-1     ,j,k-1     ) != not_covered)) {
                                                xzderiv = 0.0_rt;
                                            }

                                            if ((mask_arr(i+ratiox,j+ratioy,k) != not_covered) ||
                                                (mask_arr(i-1     ,j+ratioy,k) != not_covered) ||
                                                (mask_arr(i+ratiox,j-1     ,k) != not_covered) ||
                                                (mask_arr(i-1     ,j-1     ,k) != not_covered)) {
                                                xyderiv = 0.0_rt;
                                            }
#endif

                                            bdry(i,j,k,n) = crse(ic,jc,kc,n) +
                                                            xx * xderiv + xx * xx * xxderiv +
                                                            yy * yderiv + yy * yy * yyderiv +
                                                            zz * zderiv + zz * zz * zzderiv +
                                                            xx * yy * xyderiv + xx * zz * xzderiv + yy * zz * yzderiv;
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

