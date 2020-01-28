#ifndef _HypreMultiABec_H_
#define _HypreMultiABec_H_

#include "HypreABec.H"

#include <list>
#include "HYPRE_sstruct_ls.h"

///
/// An AuxVar is an intermediate variable used in building up a
/// linear system, defined as a linear combination of other
/// variables which are either true independent degrees of freedom
/// or other AuxVars.
///
/// @class AuxVar
/// @brief
///
class AuxVar {

 protected:

  class Connex {
   public:
    Connex() {
      other = NULL;
    }

///
/// @param p
/// @param r
///
    Connex(AuxVar* p, amrex::Real r) {
      val   = r;
      other = p;
    }

///
/// @param lev
/// @param iv
/// @param r
///
    Connex(int lev, const amrex::IntVect& iv, amrex::Real r) {
      val   = r;
      index = iv;
      level = lev;
      other = NULL;
    }

///
/// @param c
///
    bool same_target(const Connex& c) {
      return ((other != NULL)
              ? (other == c.other)
              : (c.other == NULL && level == c.level && index == c.index));
    }

    amrex::Real val;
    amrex::IntVect index;
    int level; //, grid; // may not need grid
    AuxVar *other;
  };

 public:

  AuxVar() : slave_flag(0) {
  }


///
/// @param p
/// @param r
///
  void push(AuxVar* p, amrex::Real r) {
    a.push_back(Connex(p,r));
  }


///
/// @param lev
/// @param iv
/// @param r
///
  void push(int lev, const amrex::IntVect& iv, amrex::Real r) {
    a.push_back(Connex(lev,iv,r));
  }


///
/// @param p
///
  void push_slave(AuxVar* p) {

///
/// @param p->slave_flag
///
    if (p->slave_flag == 0) { // don't count the same slave twice
      a.push_back(Connex(p,1.0));
      p->slave_flag = 1;
    }
  }

  bool empty() {
    return a.empty();
  }

  bool slave() {
    return slave_flag;
  }

  void collapse();


///
/// @param levels
/// @param cells
///
  int get_locations(amrex::Vector<int>& levels, amrex::Vector<amrex::IntVect>& cells);


///
/// @param values
///
  int get_coeffs(amrex::Vector<amrex::Real>& values);

  void clear();

 protected:

  std::list<Connex> a;
  int slave_flag;
};

///
/// AuxVarBox is like a FAB, except that an AuxVar is stored at every
/// grid point rather than a double.  It is a minimal container class
/// with much less functionality than even a BaseFab.
/// Could a BaseFab<AuxVar> be used instead? (AuxVar does not have arithmetic)
///
/// @class AuxVarBox
/// @brief
///
class AuxVarBox {

 public:


///
/// @param bx
///
  AuxVarBox(const amrex::Box& bx) : domain(bx) {
    int numpts = domain.numPts();
    dptr = new AuxVar[numpts];
  }

  ~AuxVarBox() {
    delete[] dptr;
  }

  AuxVar& operator()(const amrex::IntVect& p) {
    BL_ASSERT(!(dptr == 0));
    BL_ASSERT(domain.contains(p));
    return dptr[domain.index(p)];
  }

  const amrex::Box& box() const {
    return domain;
  }

 protected:
  amrex::Box domain;
  AuxVar *dptr;
};

///
/// @class BndryAuxVarBase
/// @brief
///
class BndryAuxVarBase {

 public:


///
/// @param dm
///
  BndryAuxVarBase (const amrex::DistributionMapping& dm)
      : distributionMap(dm) {}

///
/// Location type explanation ---
///  INTERIOR: cells along insides of grid faces, overlaps marked "slave".
///  EXTERIOR: cells along outsides of grid faces, overlaps (from other
///            grids this time) marked "slave".
///  GHOST: cells along outsides of grid faces, overlaps not marked.
///  (Note that GHOST type is located at cell centers, but the overlap
///  behavior more logically resembles that of face-centered data.  It
///  is also analogous to ghost cells, in that ghost cells from
///  different grids are distinct memory locations even if they have
///  the same cell index.)
///
  enum Location { INTERIOR, EXTERIOR, GHOST };

///
/// This base class provides parallel iterator functionality,
/// since I'm too lazy to write true iterator classes:
///
  int firstLocal();

///
/// @param i
///
  int nextLocal(int i);

///
/// @param i
///
  bool isValid(int i) {
    return (i < distributionMap.ProcessorMap().size());
  }

 protected:

  amrex::DistributionMapping        distributionMap;
};


///
/// @class BndryAuxVar

/// @brief BndryAuxVar is like a BndryRegister with AuxVars in it.  This too
/// is a rather minimal container class.  It provides its own
/// pseudo-iterator instead of building on something like FabSet.
///
class BndryAuxVar : public BndryAuxVarBase {

 public:


///
/// @param _grids
/// @param _dmap
/// @param loc
///
  BndryAuxVar(const amrex::BoxArray& _grids,
	      const amrex::DistributionMapping& _dmap,
	      Location loc);

  //~BndryAuxVar(); destructor not needed

  AuxVarBox & operator() (const amrex::Orientation ori, int i) {
    return *aux[ori][i];
  }

 protected:

  amrex::BoxArray            grids;
  amrex::Vector<std::unique_ptr<AuxVarBox> >   aux[2*BL_SPACEDIM];
};

///
/// @class CrseBndryAuxVar
///
/// @brief CrseBndryAuxVar is like a BndryAuxVar organized by coarse grid
/// rather than by fine grid.  Entries only exist where the coarse grid
/// is owned by the local processor.  Each entry has subarrays
/// containing information for all fine grid edges falling within
/// the local coarse grid.  Note that while cgrids defines the
/// coarse grid locations and fgrids defines the fine grid locations,
/// both of these BoxArrays must be at the same level of refinement!
/// This class does not know the refinement ratio or problem domain.
///
/// Mask arrays are also built, analogous to those in BndryData.
///
class CrseBndryAuxVar : public BndryAuxVarBase {

 public:


///
/// @param _cgrids
/// @param _cdmap
/// @param _fgrids
/// @param loc
///
  CrseBndryAuxVar(const amrex::BoxArray& _cgrids,
		  const amrex::DistributionMapping& _cdmap,
                  const amrex::BoxArray& _fgrids, Location loc);


///
/// @param other
/// @param loc
///
  CrseBndryAuxVar(const CrseBndryAuxVar& other, Location loc);


///
/// @param _cgrids
/// @param _fgrids
/// @param other
/// @param loc
///
  CrseBndryAuxVar(const amrex::BoxArray& _cgrids,
                  const amrex::BoxArray& _fgrids,
                  const CrseBndryAuxVar& other, Location loc);


///
/// @param loc
///
  void reinitialize_connections(Location loc); ///< use for new discretization


///
/// @param rat
/// @param ncomp
///
  void buildFaceData(amrex::IntVect& rat, int ncomp = 1);


///
/// @param rat
/// @param ncomp
///
  void rebuildFaceData(amrex::IntVect& rat, int ncomp = 1); ///< use when ncomp changes


///
/// @param ori
/// @param src
/// @param srccomp
/// @param destcomp
/// @param numcomp
///
  void loadFaceData(const amrex::Orientation ori,
                    amrex::MultiFab& src,
                    int srccomp,
                    int destcomp,
                    int numcomp = 1);

  //~CrseBndryAuxVar(); destructor not needed


///
/// @param ori
/// @param i
///
  int size (const amrex::Orientation ori, int i) const { return aux[ori][i].size(); }

  AuxVarBox& operator() (const amrex::Orientation ori, int i, int j) {
    return *aux[ori][i][j];
  }


///
/// @param ori
/// @param i
/// @param j
///
  amrex::Mask& mask (const amrex::Orientation ori, int i, int j) {
    return *msk[ori][i][j];
  }


///
/// @param ori
/// @param i
/// @param j
///
  amrex::FArrayBox& faceData (const amrex::Orientation ori, int i, int j) {
    return *face_data[ori][i][j];
  }

 protected:


///
/// @param loc
///
  void initialize_slaves(Location loc);

  amrex::BoxArray                    cgrids, fgrids;
  amrex::Vector< amrex::Vector<std::unique_ptr<AuxVarBox> > >  aux[2*BL_SPACEDIM];
  amrex::Vector< amrex::Vector<std::unique_ptr<amrex::Mask> > >       msk[2*BL_SPACEDIM];
  amrex::Vector< amrex::Vector<int> >                          fine_index[2*BL_SPACEDIM];
  amrex::Vector< amrex::Vector<std::unique_ptr<amrex::FArrayBox> > >  face_data[2*BL_SPACEDIM];
};

///
/// @class HypreMultiABec
/// @brief
///
class HypreMultiABec {

 public:

///
/// solver_flag should be set directly from the inputs file as
/// hmabec.solver_flag.  The value passed in to the constructor will
/// be used as a default if that variable does not occur in inputs.
///
/// @param _crse_level
/// @param _fine_level
/// @param _solver_flag
///
  HypreMultiABec(int _crse_level, int _fine_level, int _solver_flag);

  virtual ~HypreMultiABec();


///
/// @param level
/// @param _geom
/// @param _grids
/// @param _dmap
/// @param _crse_ratio
///
  void addLevel(int level,
		const amrex::Geometry& _geom,
		const amrex::BoxArray& _grids,
		const amrex::DistributionMapping& _dmap,
		amrex::IntVect         _crse_ratio);

  int crseLevel() {
    return crse_level;
  }
  int fineLevel() {
    return fine_level;
  }


///
/// @param level
/// @param _bd
/// @param _comp
///
  void setBndry(int level, NGBndry& _bd, int _comp = 0) {
      bd[level] = &_bd;
      bdcomp = _comp;
  }

///
/// @param level
///
  NGBndry& bndryData(int level) {
    return *bd[level];
  }
  static amrex::Real& fluxFactor() {
    return flux_factor;
  }


///
/// @param r
/// @param reg
/// @param ori
/// @param geom
///
  static void getFaceMetric(amrex::Vector<amrex::Real>& r,
                            const amrex::Box& reg,
                            const amrex::Orientation& ori,
                            const amrex::Geometry& geom) {

///
/// @param r
/// @param reg
/// @param ori
/// @param geom
///
    HypreABec::getFaceMetric(r, reg, ori, geom);
  }


///
/// @param v
///
  void setVerbose(int v) {
    verbose = v;
  }

///
/// buildMatrixStructure builds the graph structure for the matrix.
/// It can be called once addLevel and setBndry have been called for
/// every level.  No numerical data is used---the only info used from
/// the boundary object are the masks at the coarse-fine interfaces.
/// If there is only one level then it is not necessary for the call
/// to setBndry to preceed this one.
///
/// buildMatrixStructure also allocates the coefficient arrays and
/// initializes the vector objects as well as the matrix object, so
/// it must be called before any numerical data can be loaded into
/// those structures.
///
  void buildMatrixStructure();


///
/// @param alpha
/// @param beta
///
  void setScalars(amrex::Real alpha, amrex::Real beta);

  amrex::Real getAlpha() const {
    return alpha;
  }
  amrex::Real getBeta() const {
    return beta;
  }


///
/// @param level
/// @param &a
///
  void aCoefficients(int level, const amrex::MultiFab &a);

///
/// @param level
/// @param &b
/// @param dir
///
  void bCoefficients(int level, const amrex::MultiFab &b, int dir);


///
/// @param level
/// @param &Spa
///
  void SPalpha(int level, const amrex::MultiFab &Spa);


///
/// @param level
///
  const amrex::MultiFab& aCoefficients(int level) {
    return *acoefs[level];
  }

///
/// @param level
/// @param dir
///
  const amrex::MultiFab& bCoefficients(int level, int dir) {
      return (*bcoefs[level])[dir];
  }

  virtual void loadMatrix(); ///< once all level coeffs and scalars have been set
  void finalizeMatrix(); ///< once matrix load completed in derived classes also

///
/// The argument inhom in the following methods formerly defaulted
/// to 1.  For greater type safety (to avoid confusion with icomp) it
/// is now an enum with no default.  The argument icomp is always a
/// component number for the independent variable, whether it is
/// called Er, vector, or dest.
///
/// Call loadLevelVectors if dest and rhs are available at the
/// same time, otherwise call the other two.  In the latter case
/// the rhs MultiFab will be altered with boundary condition info
/// if inhom is true (but the main reason to call the two separate
/// functions is if rhs is a temporary, so this shouldn't be a problem).
///
/// The B vector depends on not only rhs, but also on boundary condition
/// information and bCoefficients, so these must be set before the B vector
/// can be initialized.  finalizeMatrix need not have been called yet.
///
  virtual void loadLevelVectors(int level,
                                amrex::MultiFab& dest,
                                int icomp,
                                amrex::MultiFab& rhs, // will not be altered
                                BC_Mode inhom);

///
/// @param level
/// @param dest
/// @param icomp
///
          void loadLevelVectorX(int level,
                                amrex::MultiFab& dest,
                                int icomp);
  virtual void loadLevelVectorB(int level,
                                amrex::MultiFab& rhs, // will be altered
                                BC_Mode inhom);

  void finalizeVectors(); ///< once all level values have been loaded

///
/// setup, solve and clear phases separated so that
/// multiple calls to solve can be made:
///
/// finalizeMatrix should be called before setupSolver.  The calls
/// to set and finalize the vectors can be put off until after
/// setupSolver if desired.  It should be ok to call setupSolver
/// once and then call solve multiple times with different
/// right hand sides.
///
/// @param _reltol
/// @param _abstol
/// @param maxiter
///
  void setupSolver(amrex::Real _reltol, amrex::Real _abstol, int maxiter);
  void solve();

///
/// @param level
/// @param dest
/// @param icomp
///
  void getSolution(int level, amrex::MultiFab& dest, int icomp);

///
/// This is the 2-norm of the complete rhs, including b.c. contributions
///
  amrex::Real getAbsoluteResidual();

  void clearSolver();


///
/// @param level
/// @param Flux
/// @param Er
/// @param icomp
/// @param inhom
///
  void boundaryFlux(int level,
		    amrex::MultiFab* Flux,
		    amrex::MultiFab& Er,
                    int icomp,
		    BC_Mode inhom);


///
/// @param level
/// @param prod
///
  void getProduct(int level, amrex::MultiFab& prod);

 protected:

  int crse_level, fine_level, solver_flag;

  amrex::Vector<amrex::Geometry> geom;
  amrex::Vector<amrex::BoxArray> grids;
  amrex::Vector<amrex::DistributionMapping> dmap;
  amrex::Vector<amrex::IntVect> fine_ratio;
  amrex::Vector<NGBndry*> bd;
  int bdcomp; ///< component number used for bd

  amrex::Vector< amrex::Vector<amrex::BoxArray> > subgrids;

  amrex::Vector<std::unique_ptr<amrex::MultiFab> > acoefs;
  amrex::Vector<std::unique_ptr<amrex::Array<amrex::MultiFab, BL_SPACEDIM> > > bcoefs;
  amrex::Real alpha, beta;
  amrex::Real reltol, abstol;

  amrex::Vector<std::unique_ptr<amrex::MultiFab> > SPa; ///< LO_SANCHEZ_POMRANING alpha

  amrex::Vector<std::unique_ptr<BndryAuxVar> >    cintrp, ederiv;
  amrex::Vector<std::unique_ptr<CrseBndryAuxVar> > c_cintrp, c_ederiv, c_entry;

  int verbose, verbose_threshold, bho, use_subgrids;

  HYPRE_SStructGrid     hgrid;
  HYPRE_SStructStencil  stencil;
  HYPRE_SStructGraph    graph;
  HYPRE_SStructMatrix   A, A0;
  HYPRE_SStructVector   b;
  HYPRE_SStructVector   x;
  HYPRE_SStructSolver   sstruct_solver;
  HYPRE_SStructSolver   sstruct_precond;
  HYPRE_Solver          solver;
  HYPRE_Solver          precond;
  int                   ObjectType;

  static amrex::Real flux_factor;

  // static utility functions follow:


///
/// @param x
/// @param part
/// @param reg
/// @param sgr
/// @param *vec
///
  static void vectorSetBoxValues(HYPRE_SStructVector x,
                                 int part,
                                 const amrex::Box& reg,
                                 const amrex::BoxArray& sgr,
                                 amrex::Real *vec);


///
/// @param x
/// @param part
/// @param reg
/// @param sgr
/// @param f
/// @param fcomp
///
  static void vectorGetBoxValues(HYPRE_SStructVector x,
                                 int part,
                                 const amrex::Box& reg,
                                 const amrex::BoxArray& sgr,
                                 amrex::FArrayBox& f, int fcomp);

#if (BL_SPACEDIM == 1)
  static int vl[2];
  static int vh[2];
#endif

///
/// The implicit assumption in the following is that vl[1], vh[1] are
/// never changed from 0:
///
/// @param i
///
  static int* getV1(const amrex::IntVect& i) {
#if (BL_SPACEDIM == 1)
    vl[0] = i[0];
    return vl;
#else

///
/// @param int*
///
    return (int*) i.getVect();
#endif
  }


///
/// @param i
///
  static int* getV2(const amrex::IntVect& i) {
#if (BL_SPACEDIM == 1)
    vh[0] = i[0];
    return vh;
#else

///
/// @param int*
///
    return (int*) i.getVect();
#endif
  }


///
/// @param b
///
  static int* loV(const amrex::Box& b) {
#if (BL_SPACEDIM == 1)
    vl[0] = b.smallEnd(0);
    return vl;
#else

///
/// @param int*
///
    return (int*) b.loVect();
#endif
  }


///
/// @param b
///
  static int* hiV(const amrex::Box& b) {
#if (BL_SPACEDIM == 1)
    vh[0] = b.bigEnd(0);
    return vh;
#else

///
/// @param int*
///
    return (int*) b.hiVect();
#endif
  }

};

#endif
