#ifndef CASTRO_RAD_INTERPBNDRYDATA_H
#define CASTRO_RAD_INTERPBNDRYDATA_H

#include <iostream>

#include <AMReX_BLassert.H>
#include <AMReX_BoxArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Orientation.H>
#include <AMReX_Mask.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <RadBoundCond.H>
#include <RadBndryData.H>


//@Man:
/*@Memo:
        An InterpBndryData object adds to a BndryData object the ability to
        manipulate and set the data stored in the boundary cells.
*/
/*@Doc:
        The "Interpbndrydata" class is a virtual base class derived from
        BndryData.  It is intended to provide a more physical method for
        filling boundary-related data.  Boundary values in a BndryData object
        are stored in FabSets around each grid in the domain, and the
        InterpBndryData class provides a mechanism for filling these FabSets,
        consistent with AMR-like numerical discretizations.  When asked to
        set it's boundary values, an InterpBndryData object:
        \begin{enumerate}
        \item Fills with physical boundary values if the FAB is on the
        domain boundary ({\it the corresponding values are presumed to be
        stored in the ghost cells of a MultiFab given to the boundary filling
        routine})
        \item Fills on intersection with data from the VALID region of the
        input MultiFab, ensuring that adjacent FABs use consistent data at
        their intersection, and otherwise,
        \item Fills with values interpolated from a coarser FAB that
        bounds the cells that do not meet the above two criteria
        \end{enumerate}

        This class does NOT provide a copy constructor or assignment operator.

*/

class RadInterpBndryData : public RadBndryData
{
public:
//@ManMemo: Constructor
//@ManDoc: default constructor
    RadInterpBndryData() : RadBndryData() {};
//@ManDoc: constructor for given BoxArray, etc
    RadInterpBndryData(const amrex::BoxArray& _grids, const amrex::DistributionMapping& _dmap,
                       int _ncomp, const amrex::Geometry& geom);
//@ManDoc: interpret user BCs to those expected by LinOp: takes integer ratio and
// is a wrapper to the IntVect version
    void setBndryConds(const amrex::BCRec& phys_bc,
                       const amrex::Geometry& geom, int ratio);
//@ManDoc: interpret user BCs to those expected by LinOp: takes IntVect ratio
    virtual void setBndryConds(const amrex::BCRec& phys_bc,
                               const amrex::Geometry& geom, amrex::IntVect& ratio) = 0;
//@ManDoc: set bndry values at coarse level (non interpolation performed)
    void setBndryValues(const amrex::MultiFab& mf, int mf_start,
                        int bnd_start, int num_comp,
                        const amrex::BCRec& phys_bc);
//@ManDoc: set bndry values at fine level, performing necessary interpolations
    void setBndryValues(amrex::BndryRegister& crse, int c_start,
                        const amrex::MultiFab& fine, int f_start,
                        int bnd_start, int num_comp, amrex::IntVect& ratio,
                        const amrex::BCRec& phys_bc);
//@ManDoc: set bndry values to provided value
    void setBndryValues(amrex::Real bv);
};


#endif

