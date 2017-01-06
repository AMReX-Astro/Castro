#ifndef _HypreMultiABec_H_
#define _HypreMultiABec_H_

#include "HypreABec.H"

#include "HYPRE_sstruct_ls.h"

// An AuxVar is an intermediate variable used in building up a
// linear system, defined as a linear combination of other
// variables which are either true independent degrees of freedom
// or other AuxVars.

class AuxVar {

 protected:

  class Connex {
   public:
    Connex() {
      other = NULL;
    }
    Connex(AuxVar* p, Real r) {
      val   = r;
      other = p;
    }
    Connex(int lev, const IntVect& iv, Real r) {
      val   = r;
      index = iv;
      level = lev;
      other = NULL;
    }
    bool same_target(const Connex& c) {
      return ((other != NULL)
              ? (other == c.other)
              : (c.other == NULL && level == c.level && index == c.index));
    }

    Real val;
    IntVect index;
    int level; //, grid; // may not need grid
    AuxVar *other;
  };

 public:

  AuxVar() : slave_flag(0) {
  }

  void push(AuxVar* p, Real r) {
    a.push_back(Connex(p,r));
  }

  void push(int lev, const IntVect& iv, Real r) {
    a.push_back(Connex(lev,iv,r));
  }

  void push_slave(AuxVar* p) {
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

  int get_locations(Array<int>& levels, Array<IntVect>& cells);

  int get_coeffs(Array<Real>& values);

  void clear();

 protected:

  std::list<Connex> a;
  int slave_flag;
};

// AuxVarBox is like a FAB, except that an AuxVar is stored at every
// grid point rather than a double.  It is a minimal container class
// with much less functionality than even a BaseFab.

// Could a BaseFab<AuxVar> be used instead? (AuxVar does not have arithmetic)

class AuxVarBox {

 public:

  AuxVarBox(const Box& bx) : domain(bx) {
    int numpts = domain.numPts();
    dptr = new AuxVar[numpts];
  }

  ~AuxVarBox() {
    delete[] dptr;
  }

  AuxVar& operator()(const IntVect& p) {
    BL_ASSERT(!(dptr == 0));
    BL_ASSERT(domain.contains(p));
    return dptr[domain.index(p)];
  }

  const Box& box() const {
    return domain;
  }

 protected:
  Box domain;
  AuxVar *dptr;
};

class BndryAuxVarBase {

 public:

  // Location type explanation ---
  //  INTERIOR: cells along insides of grid faces, overlaps marked "slave".
  //  EXTERIOR: cells along outsides of grid faces, overlaps (from other
  //            grids this time) marked "slave".
  //  GHOST: cells along outsides of grid faces, overlaps not marked.
  //  (Note that GHOST type is located at cell centers, but the overlap
  //  behavior more logically resembles that of face-centered data.  It
  //  is also analogous to ghost cells, in that ghost cells from
  //  different grids are distinct memory locations even if they have
  //  the same cell index.)

  enum Location { INTERIOR, EXTERIOR, GHOST };

  // This base class provides parallel iterator functionality,
  // since I'm too lazy to write true iterator classes:

  int firstLocal();
  int nextLocal(int i);
  bool isValid(int i) {
    return (i < distributionMap.ProcessorMap().size() - 1);
  }

 protected:

  DistributionMapping        distributionMap;
};

// BndryAuxVar is like a BndryRegister with AuxVars in it.  This too
// is a rather minimal container class.  It provides its own
// pseudo-iterator instead of building on something like FabSet.

class BndryAuxVar : public BndryAuxVarBase {

 public:

  BndryAuxVar(const BoxArray& _grids, Location loc);

  //~BndryAuxVar(); destructor not needed, PArrayManage is used

  PArray<AuxVarBox>& operator()(const Orientation ori) {
    return aux[ori];
  }

 protected:

  BoxArray            grids;
  PArray<AuxVarBox>   aux[2*BL_SPACEDIM];
};

// CrseBndryAuxVar is like a BndryAuxVar organized by coarse grid
// rather than by fine grid.  Entries only exist where the coarse grid
// is owned by the local processor.  Each entry has subarrays
// containing information for all fine grid edges falling within
// the local coarse grid.  Note that while cgrids defines the
// coarse grid locations and fgrids defines the fine grid locations,
// both of these BoxArrays must be at the same level of refinement!
// This class does not know the refinement ratio or problem domain.

// Mask arrays are also built, analogous to those in BndryData.

class CrseBndryAuxVar : public BndryAuxVarBase {

 public:

  CrseBndryAuxVar(const BoxArray& _cgrids,
                  const BoxArray& _fgrids, Location loc);

  CrseBndryAuxVar(const CrseBndryAuxVar& other, Location loc);

  CrseBndryAuxVar(const BoxArray& _cgrids,
                  const BoxArray& _fgrids,
                  const CrseBndryAuxVar& other, Location loc);

  void reinitialize_connections(Location loc); // use for new discretization

  void buildFaceData(IntVect& rat, int ncomp = 1);

  void rebuildFaceData(IntVect& rat, int ncomp = 1); // use when ncomp changes

  void loadFaceData(const Orientation ori,
                    MultiFab& src,
                    int srccomp,
                    int destcomp,
                    int numcomp = 1);

  //~CrseBndryAuxVar(); destructor not needed, PArrayManage is used

  Array< PArray<AuxVarBox> >& operator()(const Orientation ori) {
    return aux[ori];
  }

  Array< PArray<Mask> >& mask(const Orientation ori) {
    return msk[ori];
  }

  Array< PArray<FArrayBox> >& faceData(const Orientation ori) {
    return face_data[ori];
  }

 protected:

  void initialize_slaves(Location loc);

  BoxArray                    cgrids, fgrids;
  Array< PArray<AuxVarBox> >  aux[2*BL_SPACEDIM];
  Array< PArray<Mask> >       msk[2*BL_SPACEDIM];
  Array< Array<int> >         fine_index[2*BL_SPACEDIM];
  Array< PArray<FArrayBox> >  face_data[2*BL_SPACEDIM];
};

class HypreMultiABec {

 public:

  // solver_flag should be set directly from the inputs file as
  // hmabec.solver_flag.  The value passed in to the constructor will
  // be used as a default if that variable does not occur in inputs.

  HypreMultiABec(int _crse_level, int _fine_level, int _solver_flag);

  virtual ~HypreMultiABec();

  void addLevel(int level,
		const Geometry& _geom,
		const BoxArray& _grids,
		IntVect         _crse_ratio);

  int crseLevel() {
    return crse_level;
  }
  int fineLevel() {
    return fine_level;
  }

  void setBndry(int level, NGBndry& _bd, int _comp = 0) {
    // clear must be done before set if already has value;
    // will not deallocate previous value because policy is NoManage.
    bd.clear(level);
    bd.set(level, &_bd);
    bdcomp = _comp;
  }
  NGBndry& bndryData(int level) {
    return bd[level];
  }
  void bndryClear(int level) {
    delete bd.remove(level);
  }
  static Real& fluxFactor() {
    return flux_factor;
  }

  static void getFaceMetric(Array<Real>& r,
                            const Box& reg,
                            const Orientation& ori,
                            const Geometry& geom) {
    HypreABec::getFaceMetric(r, reg, ori, geom);
  }

  void setVerbose(int v) {
    verbose = v;
  }

  // buildMatrixStructure builds the graph structure for the matrix.
  // It can be called once addLevel and setBndry have been called for
  // every level.  No numerical data is used---the only info used from
  // the boundary object are the masks at the coarse-fine interfaces.
  // If there is only one level then it is not necessary for the call
  // to setBndry to preceed this one.

  // buildMatrixStructure also allocates the coefficient arrays and
  // initializes the vector objects as well as the matrix object, so
  // it must be called before any numerical data can be loaded into
  // those structures.

  void buildMatrixStructure();

  void setScalars(Real alpha, Real beta);

  Real getAlpha() const {
    return alpha;
  }
  Real getBeta() const {
    return beta;
  }

  void aCoefficients(int level, const MultiFab &a);
  void bCoefficients(int level, const MultiFab &b, int dir);

  void SPalpha(int level, const MultiFab &Spa);

  const MultiFab& aCoefficients(int level) {
    return acoefs[level];
  }
  const MultiFab& bCoefficients(int level, int dir) {
    return bcoefs[level][dir];
  }

  virtual void loadMatrix(); // once all level coeffs and scalars have been set
  void finalizeMatrix(); // once matrix load completed in derived classes also

  // The argument inhom in the following methods formerly defaulted
  // to 1.  For greater type safety (to avoid confusion with icomp) it
  // is now an enum with no default.  The argument icomp is always a
  // component number for the independent variable, whether it is
  // called Er, vector, or dest.

  // Call loadLevelVectors if dest and rhs are available at the
  // same time, otherwise call the other two.  In the latter case
  // the rhs MultiFab will be altered with boundary condition info
  // if inhom is true (but the main reason to call the two separate
  // functions is if rhs is a temporary, so this shouldn't be a problem).

  // The B vector depends on not only rhs, but also on boundary condition
  // information and bCoefficients, so these must be set before the B vector
  // can be initialized.  finalizeMatrix need not have been called yet.

  virtual void loadLevelVectors(int level,
                                MultiFab& dest,
                                int icomp,
                                MultiFab& rhs, // will not be altered
                                BC_Mode inhom);
          void loadLevelVectorX(int level,
                                MultiFab& dest,
                                int icomp);
  virtual void loadLevelVectorB(int level,
                                MultiFab& rhs, // will be altered
                                BC_Mode inhom);

  void finalizeVectors(); // once all level values have been loaded

  // setup, solve and clear phases separated so that
  // multiple calls to solve can be made:

  // finalizeMatrix should be called before setupSolver.  The calls
  // to set and finalize the vectors can be put off until after
  // setupSolver if desired.  It should be ok to call setupSolver
  // once and then call solve multiple times with different
  // right hand sides.

  void setupSolver(Real _reltol, Real _abstol, int maxiter);
  void solve();
  void getSolution(int level, MultiFab& dest, int icomp);

  // This is the 2-norm of the complete rhs, including b.c. contributions
  Real getAbsoluteResidual();

  void clearSolver();

  void boundaryFlux(int level,
		    MultiFab* Flux,
		    MultiFab& Er,
                    int icomp,
		    BC_Mode inhom);

  // Matrix apply (matrix-vector multiply) only implemented for single level

  void initializeApplyLevel(int level,
			    MultiFab& product,
			    MultiFab& vector,
                            int icomp,
			    BC_Mode inhom);
  void apply();
  void getProduct(int level, MultiFab& prod);

 protected:

  int crse_level, fine_level, solver_flag;

  Array<Geometry> geom;
  Array<BoxArray> grids;
  Array<IntVect> fine_ratio;
  PArray<NGBndry> bd;
  int bdcomp; // component number used for bd

  Array< Array<BoxArray> > subgrids;

  PArray<MultiFab> acoefs;
  PArray< Tuple<MultiFab, BL_SPACEDIM> > bcoefs;
  Real alpha, beta;
  Real reltol, abstol;

  PArray<MultiFab> SPa; // LO_SANCHEZ_POMRANING alpha

  PArray<BndryAuxVar>     cintrp, ederiv;
  PArray<CrseBndryAuxVar> c_cintrp, c_ederiv, c_entry;

  int verbose, verbose_threshold, bho, use_subgrids;

  HYPRE_SStructGrid     grid;
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

  static Real flux_factor;

  // static utility functions follow:

  static void vectorSetBoxValues(HYPRE_SStructVector x,
                                 int part,
                                 const Box& reg,
                                 const BoxArray& sgr,
                                 Real *vec);

  static void vectorGetBoxValues(HYPRE_SStructVector x,
                                 int part,
                                 const Box& reg,
                                 const BoxArray& sgr,
                                 FArrayBox& f, int fcomp);

#if (BL_SPACEDIM == 1)
  static int vl[2];
  static int vh[2];
#endif

// The implicit assumption in the following is that vl[1], vh[1] are
// never changed from 0:

  static int* getV1(const IntVect& i) {
#if (BL_SPACEDIM == 1)
    vl[0] = i[0];
    return vl;
#else
    return (int*) i.getVect();
#endif
  }

  static int* getV2(const IntVect& i) {
#if (BL_SPACEDIM == 1)
    vh[0] = i[0];
    return vh;
#else
    return (int*) i.getVect();
#endif
  }

  static int* loV(const Box& b) {
#if (BL_SPACEDIM == 1)
    vl[0] = b.smallEnd(0);
    return vl;
#else
    return (int*) b.loVect();
#endif
  }

  static int* hiV(const Box& b) {
#if (BL_SPACEDIM == 1)
    vh[0] = b.bigEnd(0);
    return vh;
#else
    return (int*) b.hiVect();
#endif
  }

};

#endif
