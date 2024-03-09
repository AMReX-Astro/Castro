#ifndef CASTRO_RAD_BNDRYDATA_H
#define CASTRO_RAD_BNDRYDATA_H

#include <AMReX_BLassert.H>
#include <AMReX_BoxArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Orientation.H>
#include <AMReX_BndryRegister.H>
#include <AMReX_Mask.H>
#include <RadBoundCond.H>
#include <AMReX_Geometry.H>

#include <iostream>

//@Man:
/*@Memo:
        A ProxyGeometry is a class constructed solely to prevent automatic
        type conversion with the Geometry::Geometry(const Box&) constructor.
        This class can easily make a Geometry object, but does not have a
        ProxyGeometry::ProxyGeometry(const Box&) constructor.
*/
/*@Doc:
        A ProxyGeometry object is, for most purposes, merely a wrapper to the
        Geometry class.  The Geometry class contains a single-argument
        constructor taking a const Box &, and therefore allows implicit
        type conversion.  However, that constructor requires parm-parsed
        input for the coordinate system and prob lo and hi.  This has led to
        some some mysterious error messages, by correctly compiled code
        passing in a Box where a Geometry was called for.  This little class
        is an attempt to circumvnet this little funny business.
*/
class ProxyGeometry
{
protected:
//@ManMemo: internal data
//@ManDoc: The underlying Geometry
    const amrex::Geometry& m_g;

public:
//@ManMemo: administrative functions
//@ManDoc: default constructor
    ProxyGeometry(const amrex::Geometry& g) : m_g( g ) {}
//@ManDoc: type conversion operator
    operator amrex::Geometry() const
    {  return m_g; }
};


//@Man:
/*@Memo:
        A BndryData stores and manipulates boundary data
        information on each side of each box in a BoxArray.
*/
/*@Doc:
        A BndryData contains a BndryRegister about each side of each grid in
        a Boxarray.  These data are used to store information along the
        outer layer of each grid (at the same level of coarseness), such
        as the value of boundary conditions there.  Also, for each
        boundary, this class maintains a BoundCond identifier and a
        location.  This "location" often specifies where, in physical
        space, a value, which is stored in the boundary registers, will
        be applied by a solver (although interpretation of this information
        is application-specific).

        In addition to boundary condition values, types and locations, and
        BndryDate object maintains a mask for each boundary value.  Given
        a Geometry describing the physical "domain" of application, a BndryData
        object fills the mask arrays with a mask value to indicate if the
        node is outside the domain (outside_domain), or, if not, whether it
        lays within the valid region of an adjacent grid (covered), or not
        (not_covered).  This mask data is created upon non-default class
        instantiation.
*/
class RadBndryData : public amrex::BndryRegister
{
public:

//@ManMemo: class enumeration
//@ManDoc: mask values enumeration
    enum MaskVal { covered = 0, not_covered = 1, outside_domain = 2 };

protected:
//@ManMemo: internal data
//@ManDoc: Array (on orientation) of boundary condition type specifiers
    amrex::Vector<RadBoundCond>   bcond[2*AMREX_SPACEDIM];
//@ManDoc: Array (on orientation) of boundary condition locations
    amrex::Vector<amrex::Real>    bcloc[2*AMREX_SPACEDIM];
//@ManDoc: Array (on orientation) of boundary condition mask arrays (FAB-like)
    amrex::Vector<std::unique_ptr<amrex::Mask> > masks[2*AMREX_SPACEDIM];
//@ManDoc: Domain used for mask definitions
    amrex::Geometry geom;

public:
//@ManMemo: administrative functions
//@ManDoc: default constructor
    RadBndryData() : amrex::BndryRegister() {};
//@ManDoc: constructor specifying number of components and box of physical domain (cell-centered)
    RadBndryData(const amrex::BoxArray& _grids, const amrex::DistributionMapping& _dmap,
                 int _ncomp, const ProxyGeometry& geom);
//@ManDoc: destructor
    virtual ~RadBndryData();
//@ManDoc: allocate bndry fabs along given face
    void define(const amrex::BoxArray& _grids, const amrex::DistributionMapping& _dmap,
                int _ncomp, const ProxyGeometry& geom);

//@ManMemo: I/O functions
//@ManDoc: write to output stream
    friend std::ostream& operator << (std::ostream&, const RadBndryData&);

//@ManMemo: access functions
//@ManDoc: return FabSet on given face
    const amrex::FabSet& bndryValues(const amrex::Orientation& _face) const
    {
        return bndry[_face];
    }
//@ManDoc: return boundary location on given face
    const amrex::Vector<amrex::Real>& bndryLocs(const amrex::Orientation& _face) const
    {
        return bcloc[_face];
    }
//@ManDoc: return boundary type specifyer on given face
    const amrex::Vector<RadBoundCond>& bndryConds(const amrex::Orientation& _face) const
    {
        return bcond[_face];
    }
//@ManDoc: return boundary mask on given face
    const amrex::Mask& bndryMasks (const amrex::Orientation& _face, int i) const
    {
        return *masks[_face][i];
    }
//@ManDoc: return domain used to define masks
    const amrex::Box& getDomain() const
    {
        return geom.Domain();
    }
//@ManDoc: return geometry used to define masks
    const amrex::Geometry& getGeom() const
    {
        return geom;
    }
//@ManDoc: set values of boundary Fab for given orientation on nth grid
    void setValue(const amrex::Orientation& _face, int _n, amrex::Real _val)
    {
        bndry[_face][_n].setVal<amrex::RunOn::Host>(_val);
    }
//@ManDoc: set mask values for given orientation on nth grid
    void setMaskValue(const amrex::Orientation& _face, int _n, int _val)
    {
        masks[_face][_n]->setVal<amrex::RunOn::Host>(_val);
    }
//@ManDoc: set boundary type specifier for given orientation on nth grid
    void setBoundCond(const amrex::Orientation& _face, int _n, const RadBoundCond& _bcn)
    {
        bcond[_face][_n] = _bcn;
    }
//@ManDoc: set boundary location for given orientation on nth grid
    void setBoundLoc(const amrex::Orientation& _face, int _n, amrex::Real _val)
    {
        bcloc[_face][_n] = _val;
    }
//@ManDoc: implement public access to const BndryRegister::operator[]
    const amrex::FabSet &operator[](const amrex::Orientation &_face) const {
      return amrex::BndryRegister::bndry[_face];
    }
//@ManDoc: implement public access to BndryRegister::operator[]
    amrex::FabSet &operator[](const amrex::Orientation &_face) {
      return amrex::BndryRegister::bndry[_face];
    }


private:
    //
    // Disabled!
    //
//@ManDoc: copy constructor
    RadBndryData(const RadBndryData& src);
//@ManDoc: copy operator
    RadBndryData& operator = (const RadBndryData& src);
};

#endif

