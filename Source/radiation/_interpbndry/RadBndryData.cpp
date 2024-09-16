#include <limits.h>
#include <RadBndryData.H>
#include <AMReX_Utility.H>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
REAL* fabdat = (fab).dataPtr();
#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const REAL* fabdat = (fab).dataPtr();

RadBndryData::RadBndryData(const BoxArray& _grids, const DistributionMapping& _dmap,
                           int _ncomp, const ProxyGeometry& _geom)
    : BndryRegister(), geom(_geom)
{
    define(_grids,_dmap,_ncomp,_geom);
}

// // copy constructor
// RadBndryData::RadBndryData(const RadBndryData& src)
// {
//     (*this) = src;
// }

std::ostream& operator << (std::ostream& os, const RadBndryData &mgb)
{
    const BoxArray& grds = mgb.boxes();
    int ngrds = grds.size();
    os << "[RadBndryData with " << ngrds << " grids:\n";
    for (int grd = 0; grd < ngrds; grd++){
        for (OrientationIter face; face; ++face) {
            Orientation f = face();
            os << "::: face " << f << " of grid " << grds[grd] << "\n";
            os << "BC = " << mgb.bcond[f][grd]
               << " LOC = " << mgb.bcloc[f][grd] << "\n";
            os << *mgb.masks[f][grd];
            os << mgb.bndry[f][grd];
        }
        os << "--------------------------------------------------" << std::endl;
    }
    return os;
}

void
RadBndryData::define(const BoxArray& _grids, const DistributionMapping& _dmap,
                     int _ncomp, const ProxyGeometry& _geom)
{
    geom = _geom;
    BndryRegister::setBoxes(_grids);
    int len = grids.size();
    BL_ASSERT( len > 0 );

    for (OrientationIter fi; fi; ++fi) {
        Orientation face = fi();
        int coord_dir = face.coordDir();
        masks[face].resize(len);
        bcloc[face].resize(len);
        bcond[face].resize(len);

        BndryRegister::define(face,IndexType::TheCellType(),0,1,0,_ncomp,_dmap);

        // alloc mask and set to quad_interp value
        //for (int k = 0; k < len; k++) {
        for (FabSetIter bi(bndry[face]); bi.isValid(); ++bi) {
            int k = bi.index();
            Box face_box = amrex::adjCell(grids[k],face,1);

            // extend box in directions orthogonal to face normal
            for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
                if (dir == coord_dir) continue;
                face_box.grow(dir,1);
            }
            masks[face][k].reset(new Mask(face_box));
            Mask *m = masks[face][k].get();
            m->setVal<RunOn::Host>(outside_domain,0);
            Box dbox(geom.Domain());
            dbox &= face_box;
            m->setVal<RunOn::Host>(not_covered,dbox,0);
            // now have to set as not_covered the periodic translates as well
            if( geom.isAnyPeriodic() ){
              Box dombox(geom.Domain());
              Vector<IntVect> pshifts(27);
              geom.periodicShift( dombox, face_box, pshifts );
              for( int iiv=0; iiv<pshifts.size(); iiv++){
                IntVect iv = pshifts[iiv];
                m->shift(iv);
                Box target(dombox);
                target &= m->box();
                if (target.ok()) m->setVal<RunOn::Host>(not_covered,target,0);
                m->shift(-iv);
              }
            }
              // turn mask off on intersection with grids at this level
            for (int g = 0; g < len; g++) {
                Box ovlp(grids[g]);
                ovlp &= face_box;
                if (ovlp.ok()) m->setVal<RunOn::Host>(covered,ovlp,0);
            }
            // handle special cases if is periodic
            if( geom.isAnyPeriodic() &&
                !geom.Domain().contains(face_box) ){
              Vector<IntVect> pshifts(27);
              geom.periodicShift( geom.Domain(), face_box, pshifts);
              for( int iiv=0; iiv<pshifts.size(); iiv++ ){
                IntVect iv = pshifts[iiv];
                m->shift(iv);
                for( int g=0; g<len; g++){
                  Box ovlp(grids[g]);
                  ovlp &= m->box();
                  if( ovlp.ok() ) m->setVal<RunOn::Host>(covered,ovlp,0);
                }
                m->shift(-iv);
              }
            }
        }
    }
}
