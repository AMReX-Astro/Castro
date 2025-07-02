#ifndef CASTRO_NGBNDRY_H
#define CASTRO_NGBNDRY_H

#include <RadInterpBndryData.H>
#include <AMReX_BC_TYPES.H>

enum BC_Mode { Homogeneous_BC = 0, Inhomogeneous_BC };

///
/// @class NGBndry
/// @brief
///
class NGBndry : public RadInterpBndryData
{

///
/// @param _grids
/// @param _dmap
/// @param _ncomp
/// @param _geom
///
public:
    NGBndry(const amrex::BoxArray& _grids, const amrex::DistributionMapping& _dmap,
          int _ncomp, const amrex::Geometry& _geom) :
        RadInterpBndryData(_grids,_dmap,_ncomp,_geom) { }



///
/// @param bc
/// @param phys_bc_mode
///
  virtual void setBndryFluxConds(const amrex::BCRec& bc,
                       const BC_Mode phys_bc_mode = Inhomogeneous_BC) = 0;


///
/// @param _face
///
  const amrex::Vector< std::unique_ptr<amrex::BaseFab<int> > >& bndryTypes(const amrex::Orientation& _face) const {
    return bctypearray[_face];
  }


///
/// @param _face
///
  virtual int mixedBndry(const amrex::Orientation& /* _face */) const {
    return 0;
  }

protected:

///
/// If used, these arrays must be built/deleted by a derived class:
///
  amrex::Vector< std::unique_ptr<amrex::BaseFab<int> > > bctypearray[2*AMREX_SPACEDIM];

};

///
/// @class NGBndryBld
/// @brief
///
class NGBndryBld
{
public:
  virtual NGBndry* operator()(const amrex::BoxArray& _grids,
                              const amrex::DistributionMapping& _dmap,
                              int _ncomp,
                              const amrex::Geometry& _geom) const = 0;
};

#endif
