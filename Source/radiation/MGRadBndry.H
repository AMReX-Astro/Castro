#ifndef CASTRO_MGRADBNDRY_H
#define CASTRO_MGRADBNDRY_H

#include <NGBndry.H>

///
/// @class MGRadBndry
/// @brief
///
class MGRadBndry : public NGBndry
{

///
/// @param _grids
/// @param _dmap
/// @param _ngroups
/// @param _geom
///
public:
  MGRadBndry(const amrex::BoxArray& _grids,
             const amrex::DistributionMapping& _dmap,
             const int _ngroups,
             const amrex::Geometry& _geom);

///
/// @param phys_bc
/// @param geom
/// @param ratio
///
  virtual void setBndryConds(const amrex::BCRec& phys_bc,
                             const amrex::Geometry& geom, amrex::IntVect& ratio);


///
/// @param bc
/// @param phys_bc_mode
///
  virtual void setBndryFluxConds(const amrex::BCRec& bc,
                       const BC_Mode phys_bc_mode = Inhomogeneous_BC);


///
/// @param _face
///
  virtual int mixedBndry(const amrex::Orientation& _face) const {
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


///
/// @param _ngroups
///
protected:
  static void init(const int _ngroups);
  static int ngroups;

  static int first;            ///< only set up bcval once
  static amrex::Vector<int> bcflag;    ///< 0 = use constant bcval; 1 = use function
  static amrex::Vector< amrex::Vector<amrex::Real> > bcval;
  static amrex::Real time;            ///< time value passed to bc function
  static int correction;       ///< ignore bcval when this is set
};

///
/// @class MGRadBndryBld
/// @brief
///
class MGRadBndryBld : public NGBndryBld
{
public:
  NGBndry* operator()(const amrex::BoxArray& _grids,
                      const amrex::DistributionMapping& _dmap,
                      int _ngroups,
                      const amrex::Geometry& _geom) const {

///
/// @param _grids
/// @param _dmap
/// @param _ngroups
/// @param _geom
///
      return new MGRadBndry(_grids, _dmap, _ngroups, _geom);
  }
};

#endif
