#ifndef CASTRO_RADBNDRY_H
#define CASTRO_RADBNDRY_H

#include <NGBndry.H>

///
/// @class RadBndry
/// @brief
///
class RadBndry : public NGBndry
{

///
/// @param _grids
/// @param _dmap
/// @param _geom
///
public:
  RadBndry(const amrex::BoxArray& _grids, const amrex::DistributionMapping& _dmap,
           const amrex::Geometry& _geom);

///
/// @param _grids
/// @param _dmap
/// @param _geom
/// @param bv
///
  RadBndry(const amrex::BoxArray& _grids, const amrex::DistributionMapping& _dmap,
           const amrex::Geometry& _geom, amrex::Real bv);


///
/// @param phys_bc
/// @param geom
/// @param ratio
///
  void setBndryConds(const amrex::BCRec& phys_bc,
                     const amrex::Geometry& geom, amrex::IntVect& ratio) override;


///
/// @param bc
/// @param phys_bc_mode
///
  void setBndryFluxConds(const amrex::BCRec& bc,
                         const BC_Mode phys_bc_mode = Inhomogeneous_BC) override;


///
/// @param _face
///
  int mixedBndry(const amrex::Orientation& _face) const override {
    return (bcflag[_face] > 1) ? 1 : 0;
  }


///
/// @param bc
/// @param ratio
///
  void setHomogValues(const amrex::BCRec& bc, amrex::IntVect& ratio);


///
/// @param Time
///
  static void setTime(amrex::Real Time) {
    time = Time;
  }
///
/// @warning hidden state change, use carefully!
///
  static void setCorrection() {
    correction = 1;
  }
  static void unsetCorrection() {
    correction = 0;
  }

protected:
  static void init();

///
/// @param bv
///
  static void init(amrex::Real bv);

  static int first;            ///< only set up bcval once
  static amrex::Vector<int> bcflag;    ///< 0 = use constant bcval; 1 = use function
  static amrex::Vector<amrex::Real> bcval;    ///< constant bcval across each face
  static amrex::Real time;            ///< time value passed to bc function
  static int correction;       ///< ignore bcval when this is set
};

///
/// @class RadBndryBld
/// @brief
///
class RadBndryBld : public NGBndryBld
{
public:

///
/// @param _grids
/// @param _dmap
/// @param _ncomp
/// @param _geom
///
  NGBndry* operator()(const amrex::BoxArray& _grids,
                      const amrex::DistributionMapping& _dmap,
                      int /* _ncomp */,
                      const amrex::Geometry& _geom) const override {
      return new RadBndry(_grids, _dmap, _geom);
  }
};

#endif
