#ifndef _NGBndry_H_
#define _NGBndry_H_

#include <RadInterpBndryData.H>
#include <BC_TYPES.H>
#include "LHH.H"

enum BC_Mode { Homogeneous_BC = 0, Inhomogeneous_BC };

class NGBndry : public RadInterpBndryData
{
public:
  NGBndry(const BoxArray& _grids, int _ncomp, const Geometry& _geom) :
    RadInterpBndryData(_grids,_ncomp,_geom) { }

  virtual void setBndryFluxConds(const BCRec& bc,
		       const BC_Mode phys_bc_mode = Inhomogeneous_BC) = 0;

  const Array< std::unique_ptr<BaseFab<int> > >& bndryTypes(const Orientation& _face) const {
    return bctypearray[_face];
  }

  virtual int mixedBndry(const Orientation& _face) const {
    return 0;
  }

protected:
  // If used, these arrays must be built/deleted by a derived class:
  Array< std::unique_ptr<BaseFab<int> > > bctypearray[2*BL_SPACEDIM];

private:
  NGBndry(const NGBndry& src);
  NGBndry& operator=(const NGBndry& src);
};

class NGBndryBld
{
public:
  virtual NGBndry* operator()(const BoxArray& _grids, int _ncomp,
			      const Geometry& _geom) const = 0;
};

#endif
