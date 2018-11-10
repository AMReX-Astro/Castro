
#include <climits>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_INTERP_F.H>
#include <AMReX_Interpolater.H>

#include "DG_Interpolater.H"

namespace amrex {

DG_Interpolater dg_interp;

DG_Interpolater::DG_Interpolater ()
{}

DG_Interpolater::~DG_Interpolater ()
{}

Box
DG_Interpolater::CoarseBox (const Box&     fine,
                                   const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

Box
DG_Interpolater::CoarseBox (const Box& fine,
                                   int        ratio)
{
    Box crse(amrex::coarsen(fine,ratio));
    crse.grow(1);
    return crse;
}

void
DG_Interpolater::interp (const FArrayBox& crse,
                                int              crse_comp,
                                FArrayBox&       fine,
                                int              fine_comp,
                                int              ncomp,
                                const Box&       fine_region,
                                const IntVect&   ratio,
                                const Geometry&  crse_geom,
                                const Geometry&  fine_geom,
                                Vector<BCRec>&    bcr,
                                int              actual_comp,
                                int              actual_state)
{
    BL_PROFILE("DG_Interpolater::interp()");
    BL_ASSERT(bcr.size() >= ncomp);

    amrex::Print() << " ************** NOTE NOTE ************** " << std::endl;
    amrex::Print() << " We are calling DG_interpolater but it is really just CellConsInterpolater for now" << std::endl;

    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();
    //
    // crse_bx is coarsening of target_fine_region, grown by 1.
    //
    Box crse_bx = CoarseBox(target_fine_region,ratio);
    //
    // Slopes are needed only on coarsening of target_fine_region.
    //
    Box cslope_bx(crse_bx);
    cslope_bx.grow(-1);
    //
    // Make a refinement of cslope_bx
    //
    Box fine_version_of_cslope_bx = amrex::refine(cslope_bx,ratio);
    //
    // Get coarse and fine edge-centered volume coordinates.
    //
    Vector<Real> fvc[AMREX_SPACEDIM];
    Vector<Real> cvc[AMREX_SPACEDIM];
    int dir;
    for (dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        fine_geom.GetEdgeVolCoord(fvc[dir],fine_version_of_cslope_bx,dir);
        crse_geom.GetEdgeVolCoord(cvc[dir],crse_bx,dir);
    }
    //
    // alloc tmp space for slope calc.
    //
    // In ucc_slopes and lcc_slopes , there is a slight abuse of 
    // the number of compenents argument
    // --> there is a slope for each component in each coordinate 
    //     direction
    //
    FArrayBox ucc_slopes(cslope_bx,ncomp*AMREX_SPACEDIM);
    FArrayBox lcc_slopes(cslope_bx,ncomp*AMREX_SPACEDIM);
    FArrayBox slope_factors(cslope_bx,AMREX_SPACEDIM);

    FArrayBox  cmax(cslope_bx,ncomp);
    FArrayBox  cmin(cslope_bx,ncomp);
    FArrayBox alpha(cslope_bx,ncomp);

    Real* fdat       = fine.dataPtr(fine_comp);
    const Real* cdat = crse.dataPtr(crse_comp);
    Real* ucc_xsldat = ucc_slopes.dataPtr(0);
    Real* lcc_xsldat = lcc_slopes.dataPtr(0);
    Real* xslfac_dat = slope_factors.dataPtr(0);
#if (AMREX_SPACEDIM>=2)
    Real* ucc_ysldat = ucc_slopes.dataPtr(ncomp);
    Real* lcc_ysldat = lcc_slopes.dataPtr(ncomp);
    Real* yslfac_dat = slope_factors.dataPtr(1);
#endif
#if (AMREX_SPACEDIM==3)
    Real* ucc_zsldat = ucc_slopes.dataPtr(2*ncomp);
    Real* lcc_zsldat = lcc_slopes.dataPtr(2*ncomp);
    Real* zslfac_dat = slope_factors.dataPtr(2);
#endif
    
    const int* flo    = fine.loVect();
    const int* fhi    = fine.hiVect();
    const int* clo    = crse.loVect();
    const int* chi    = crse.hiVect();
    const int* fblo   = target_fine_region.loVect();
    const int* fbhi   = target_fine_region.hiVect();
    const int* csbhi  = cslope_bx.hiVect();
    const int* csblo  = cslope_bx.loVect();
    int lin_limit     = 0;
    const int* cvcblo = crse_bx.loVect();
    const int* fvcblo = fine_version_of_cslope_bx.loVect();
    int slope_flag    = 1;

    int cvcbhi[AMREX_SPACEDIM];
    int fvcbhi[AMREX_SPACEDIM];

    for (dir=0; dir<AMREX_SPACEDIM; dir++)
    {
        cvcbhi[dir] = cvcblo[dir] + cvc[dir].size() - 1;
        fvcbhi[dir] = fvcblo[dir] + fvc[dir].size() - 1;
    }

    AMREX_D_TERM(Real* voffx = new Real[fvc[0].size()];,
           Real* voffy = new Real[fvc[1].size()];,
           Real* voffz = new Real[fvc[2].size()];);

    Vector<int> bc     = GetBCArray(bcr);
    const int* ratioV = ratio.getVect();

    amrex_linccinterp (fdat,AMREX_ARLIM(flo),AMREX_ARLIM(fhi),
                      fblo, fbhi,
                      AMREX_ARLIM(fvcblo), AMREX_ARLIM(fvcbhi),
                      cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),
                      AMREX_ARLIM(cvcblo), AMREX_ARLIM(cvcbhi),
                      ucc_xsldat, lcc_xsldat, xslfac_dat,
#if (AMREX_SPACEDIM>=2)
                      ucc_ysldat, lcc_ysldat, yslfac_dat,
#endif
#if (AMREX_SPACEDIM==3)
                      ucc_zsldat, lcc_zsldat, zslfac_dat,
#endif
                      AMREX_ARLIM(csblo), AMREX_ARLIM(csbhi),
                      csblo, csbhi,
                      &ncomp,AMREX_D_DECL(&ratioV[0],&ratioV[1],&ratioV[2]),
                      bc.dataPtr(), &slope_flag, &lin_limit,
                      AMREX_D_DECL(fvc[0].dataPtr(),fvc[1].dataPtr(),fvc[2].dataPtr()),
                      AMREX_D_DECL(cvc[0].dataPtr(),cvc[1].dataPtr(),cvc[2].dataPtr()),
                      AMREX_D_DECL(voffx,voffy,voffz),
                      alpha.dataPtr(),cmax.dataPtr(),cmin.dataPtr(),
                      &actual_comp,&actual_state);

    AMREX_D_TERM(delete [] voffx;, delete [] voffy;, delete [] voffz;);

}

}
