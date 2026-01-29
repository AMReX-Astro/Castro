#include <Radiation.H>

#include <AMReX_ParmParse.H>

#include <HypreMultiABec.H>
#include <HABEC.H>
#include <rad_util.H>
#include <AMReX_LO_BCTYPES.H>

#include <_hypre_sstruct_mv.h>
#include <HYPRE_krylov.h>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

static int ispow2(int i)
{
  return (i == 1) ? 1 : (((i <= 0) || (i & 1)) ? 0 : ispow2(i / 2));
}

Real HypreMultiABec::flux_factor = 1.0;

#if (AMREX_SPACEDIM == 1)
int HypreMultiABec::vl[2] = { 0, 0 };
int HypreMultiABec::vh[2] = { 0, 0 };
#endif

void AuxVar::collapse()
{
  // "Flattens" the dependency list.  Any entry that points to another
  // AuxVar is replaced by a copy of the dependency list from that AuxVar.
  // The values in the inserted elements are multiplied by the value from
  // the entry being replaced.  AuxVar pointers that were in the dependency
  // list are removed from it, but the AuxVars they point to are not
  // themselves modified.  Entries pointing to the same cell are then
  // combined.

  // The loops in this function are tricky, since we're modifying
  // a list at the same time we're looping through it.  for statements
  // are used in somewhat unusual ways.  Erasing an entry invalidates
  // its iterator, so before we can do that we need to have another
  // iterator pointing to an adjacent entry.  While "past-the-end" is
  // a valid target for an iterator, we assume there is no corresponding
  // "before-the-beginning" position.

  // There is one property that can be useful in dealing with the hypre
  // semi-structured interface:  Call two AuxVar objects or two Connex
  // lists "similar" if they point to the same locations in the same order.
  // If two AuxVar objects are built whose Connex lists are similar up to
  // a point, and all the latter entries point only to locations already
  // referenced in the similar initial sections, then after being
  // collapsed the two AuxVars will be similar.  If one AuxVar is then
  // used for establishing the graph entries, and the other is used for
  // initializing the matrix coefficients, then it is ok for the second
  // AuxVar to include additional entries so long as they don't point
  // to any new locations not already established in the graph.

  for (std::list<Connex>::iterator it = a.begin(); it != a.end(); ) {
    AuxVar *o = it->other;
    if (o) {
      std::list<Connex>::iterator kt = it;
      a.insert(++kt, o->a.begin(), o->a.end());

      // multiply inserted entries by it->var here:
      for (std::list<Connex>::iterator jt = it; ++jt != kt; ) {
        jt->val *= it->val;
      }

      kt = it; // entry to be erased
      ++it;
      a.erase(kt);
    }
    else {
      ++it;
    }
  }

  // Combine entries that point to same cell:

  for (std::list<Connex>::iterator it = a.begin(); it != a.end(); ++it) {
    for (std::list<Connex>::iterator jt = it; ++jt != a.end(); ) {
      if (it->same_target(*jt)) {
        it->val += jt->val;
        std::list<Connex>::iterator kt = jt;
        --jt;
        a.erase(kt);
      }
    }
  }
}

void AuxVar::clear()
{
  a.clear();
  secondary_flag = 0;
}

int AuxVar::get_locations(Vector<int>& levels, Vector<IntVect>& cells)
{
  if (secondary()) {
    return 1; // failure
  }

  int n = a.size();
  levels.resize(n);
  cells.resize(n);

  int i = 0;
  for (std::list<Connex>::iterator it = a.begin(); it != a.end(); ++it) {
    if (it->other != NULL) {
      return 2; // failure
    }
    levels[i] = it->level;
    cells[i]  = it->index;
    i++;
  }

  return 0; // success
}

int AuxVar::get_coeffs(Vector<Real>& values)
{
  if (secondary()) {
    return 1; // failure
  }

  int n = a.size();
  values.resize(n);

  int i = 0;
  for (std::list<Connex>::iterator it = a.begin(); it != a.end(); ++it) {
    if (it->other != NULL) {
      return 2; // failure
    }
    values[i] = it->val;
    i++;
  }

  return 0; // success
}

int BndryAuxVarBase::firstLocal()
{
  const int MyProc = ParallelDescriptor::MyProc();
  int i = 0;
  while (distributionMap[i] != MyProc) {
    i++;
    if (i >= distributionMap.size()) {
        i = -1;
        break;
    }
  }
  return i;
}

int BndryAuxVarBase::nextLocal(int i)
{
  const int MyProc = ParallelDescriptor::MyProc();
  do {
    i++;
    if (i >= distributionMap.size()) {
        i = -1;
        break;
    }
  } while (distributionMap[i] != MyProc);
  return i;
}

BndryAuxVar::BndryAuxVar(const BoxArray& _grids,
                         const DistributionMapping& _dmap,
                         Location loc)
    : BndryAuxVarBase(_dmap), grids(_grids)
{
  // For Location type EXTERIOR use CrseBndryAuxVar instead:
  BL_ASSERT(loc == INTERIOR || loc == GHOST);

  int inormal = (loc == INTERIOR) ? 0 : 1;

  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();
    aux[ori].resize(grids.size());
    int ishift = ori.isLow() ? 1-inormal : inormal-1;
    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      Box b = amrex::adjCell(grids[i], ori);
      aux[ori][i].reset(new AuxVarBox(b.shift(ori.coordDir(), ishift)));
    }
  }

  // Make primary-secondary connections:

  if (loc == INTERIOR) {
    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
#if 1
      // This is the current default implementation:
      for (OrientationIter omitr; omitr; ++omitr) {
        Orientation om = omitr();
        const Box& bm = aux[om][i]->box();
        for (OrientationIter ositr; ositr; ++ositr) {
          Orientation os = ositr();
          if (os.coordDir() > om.coordDir()) {
            const Box& bs = aux[os][i]->box();
            if (bm.intersects(bs)) {
              Box reg = (bm & bs);
              for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
                  (*aux[om][i])(v).push_secondary(&(*aux[os][i])(v));
              }
            }
          }
        }
      }
#elif 0
      // This was the original implementation, 2D only:
      Orientation oxlo(0, Orientation::low);
      Orientation oxhi(0, Orientation::high);
      Orientation oylo(1, Orientation::low);
      Orientation oyhi(1, Orientation::high);
      IntVect p = aux[oxlo][i]->box().smallEnd();
      (*aux[oxlo][i])(p).push_secondary(&(*aux[oylo][i])(p));
      p = aux[oxlo][i]->box().bigEnd();
      (*aux[oxlo][i])(p).push_secondary(&(*aux[oyhi][i])(p));
      p = aux[oxhi][i]->box().smallEnd();
      (*aux[oxhi][i])(p).push_secondary(&(*aux[oylo][i])(p));
      p = aux[oxhi][i]->box().bigEnd();
      (*aux[oxhi][i])(p).push_secondary(&(*aux[oyhi][i])(p));
#elif 0
      // This version is like the new default, except that
      // it loops through orientations in a different order.
      // Some primary/secondary pairs are therefore flipped, and in
      // the end the solvers return slightly different numbers.
      // (Results should be the same within the solver tolerance.)
      for (OrientationIter omitr; omitr; ++omitr) {
        Orientation om = omitr();
        const Box& bm = aux[om][i]->box();
        for (OrientationIter ositr(om); ++ositr; ) {
          Orientation os = ositr();
          const Box& bs = aux[os][i]->box();
          if (bm.intersects(bs)) {
            Box reg = (bm & bs);
            for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
                (*aux[om][i])(v).push_secondary(&(*aux[os][i])(v));
            }
          }
        }
      }
#endif
    }
  }
}

CrseBndryAuxVar::CrseBndryAuxVar(const BoxArray& _cgrids,
                                 const DistributionMapping& _cdmap,
                                 const BoxArray& _fgrids, Location loc)
    : BndryAuxVarBase(_cdmap), cgrids(_cgrids), fgrids(_fgrids)
{
  // For Location type INTERIOR use BndryAuxVar instead:
  BL_ASSERT(loc == EXTERIOR || loc == GHOST);

  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();

    aux[ori].resize(cgrids.size());
    msk[ori].resize(cgrids.size());

    fine_index[ori].resize(cgrids.size());

    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      std::list<Box> bl;
      std::list<int> fi;
      for (int j = 0; j < fgrids.size(); j++) {
        Box face = amrex::adjCell(fgrids[j], ori);
        if (cgrids[i].intersects(face)) {
          bl.push_back(face & cgrids[i]);
          fi.push_back(j);
        }
      }

      // bl now has every fine grid face in it, even those entirely
      // covered by other fine grids.  We might optimize by removing
      // the covered ones here.

      int n = bl.size();
      aux[ori][i].resize(n);
      msk[ori][i].resize(n);

      fine_index[ori][i].resize(n);

      int j = 0;
      for (std::list<int>::iterator it = fi.begin(); it != fi.end(); ++it) {
        fine_index[ori][i][j++] = *it;
      }

      j = 0;
      for (std::list<Box>::iterator it = bl.begin(); it != bl.end(); ++it) {
        aux[ori][i][j].reset(new AuxVarBox(*it));
        Box mask_box = *it;
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          if (dir == ori.coordDir())
            continue;
          mask_box.grow(dir,1);
        }
        msk[ori][i][j].reset(new Mask(mask_box));

        msk[ori][i][j]->setVal<RunOn::Host>(RadBndryData::outside_domain);

        // If we had the Geometry, we would have the problem domain
        // and would not have to loop over coarse grids below.  Also,
        // if we ever care to do this right for periodic domains we
        // will need the Geometry.  (See RadBndryData.cpp)

        for (int k = 0; k < cgrids.size(); k++) {
          if (cgrids[k].intersects(mask_box)) {
              msk[ori][i][j]->setVal<RunOn::Host>(RadBndryData::not_covered,
                                                  (cgrids[k] & mask_box), 0);
          }
        }

        for (int k = 0; k < fgrids.size(); k++) {
          if (fgrids[k].intersects(mask_box)) {
              msk[ori][i][j]->setVal<RunOn::Host>(RadBndryData::covered,
                                                  (fgrids[k] & mask_box), 0);
          }
        }

        j++;
      }
    }
  }

  initialize_secondaries(loc);
}

CrseBndryAuxVar::CrseBndryAuxVar(const CrseBndryAuxVar& other, Location loc)
    : BndryAuxVarBase(other.distributionMap), cgrids(other.cgrids), fgrids(other.fgrids)
{
  // For Location type INTERIOR use BndryAuxVar instead:
  BL_ASSERT(loc == EXTERIOR || loc == GHOST);

  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();

    aux[ori].resize(cgrids.size());
    msk[ori].resize(cgrids.size());

    fine_index[ori].resize(cgrids.size());

    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      int n = other.aux[ori][i].size();

      aux[ori][i].resize(n);
      msk[ori][i].resize(n);

      fine_index[ori][i].resize(n);

      for (int j = 0; j < n; j++) {
        fine_index[ori][i][j] = other.fine_index[ori][i][j];

        aux[ori][i][j].reset(new AuxVarBox(other.aux[ori][i][j]->box()));

        msk[ori][i][j].reset(new Mask(other.msk[ori][i][j]->box()));
        msk[ori][i][j]->copy<RunOn::Host>(*other.msk[ori][i][j]);
      }
    }
  }

  initialize_secondaries(loc);
}

CrseBndryAuxVar::CrseBndryAuxVar(const BoxArray& _cgrids,
                                 const BoxArray& _fgrids,
                                 const CrseBndryAuxVar& other, Location loc)
    : BndryAuxVarBase(other.distributionMap), cgrids(_cgrids), fgrids(_fgrids)
{
  // For Location type INTERIOR use BndryAuxVar instead:
  BL_ASSERT(loc == EXTERIOR || loc == GHOST);

  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();

    aux[ori].resize(cgrids.size());
    msk[ori].resize(cgrids.size());

    fine_index[ori].resize(cgrids.size());

    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      int n = other.aux[ori][i].size();

      aux[ori][i].resize(n);
      msk[ori][i].resize(n);

      fine_index[ori][i].resize(n);

      for (int j = 0; j < n; j++) {
        fine_index[ori][i][j] = other.fine_index[ori][i][j];

        Box face = amrex::adjCell(fgrids[fine_index[ori][i][j]], ori);
        aux[ori][i][j].reset(new AuxVarBox(face & cgrids[i]));

        Box mask_box = aux[ori][i][j]->box();
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          if (dir == ori.coordDir())
            continue;
          mask_box.grow(dir,1);
        }
        msk[ori][i][j].reset(new Mask(mask_box));

        msk[ori][i][j]->setVal<RunOn::Host>(RadBndryData::outside_domain);

        // If we had the Geometry, we would have the problem domain
        // and would not have to loop over coarse grids below.  Also,
        // if we ever care to do this right for periodic domains we
        // will need the Geometry.  (See RadBndryData.cpp)

        for (int k = 0; k < cgrids.size(); k++) {
          if (cgrids[k].intersects(mask_box)) {
              msk[ori][i][j]->setVal<RunOn::Host>(RadBndryData::not_covered,
                                                  (cgrids[k] & mask_box), 0);
          }
        }

        for (int k = 0; k < fgrids.size(); k++) {
          if (fgrids[k].intersects(mask_box)) {
              msk[ori][i][j]->setVal<RunOn::Host>(RadBndryData::covered,
                                                  (fgrids[k] & mask_box), 0);
          }
        }

      }
    }
  }

  initialize_secondaries(loc);
}

void CrseBndryAuxVar::reinitialize_connections(Location loc)
{
  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();
    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      int n = aux[ori][i].size();
      for (int j = 0; j < n; j++) {
        const Box& reg = aux[ori][i][j]->box();
        for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
            (*aux[ori][i][j])(v).clear();
        }
      }
    }
  }

  initialize_secondaries(loc);
}

void CrseBndryAuxVar::initialize_secondaries(Location loc)
{
  // Make primary-secondary connections:

  if (loc == EXTERIOR) {
    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      for (OrientationIter omitr; omitr; ++omitr) {
        Orientation om = omitr();
        for (OrientationIter ositr(om); ++ositr; ) {
          Orientation os = ositr();

          for (int jm = 0; jm < aux[om][i].size(); jm++) {
            const Box& bm = aux[om][i][jm]->box();
            for (int js = 0; js < aux[os][i].size(); js++) {
                const Box& bs = aux[os][i][js]->box();

              if (bm.intersects(bs)) {
                Box reg = (bm & bs);
                for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
                    (*aux[om][i][jm])(v).push_secondary(&(*aux[os][i][js])(v));
                }
              }

            }
          }

        }
      }
    }
  }
}

void CrseBndryAuxVar::buildFaceData(IntVect& rat, int ncomp)
{
  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();

    face_data[ori].resize(cgrids.size());

    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      int n = aux[ori][i].size();
      face_data[ori][i].resize(n);
      for (int j = 0; j < n; j++) {
          Box face_box = aux[ori][i][j]->box();
        // face_box is known to be "adjacent cell", convert to face:
        face_box.shiftHalf(ori.coordDir(), (ori.isLow() ? 1 : -1));
        face_box.refine(rat);
        face_data[ori][i][j].reset(new FArrayBox(face_box, ncomp));
      }
    }
  }
}

void CrseBndryAuxVar::rebuildFaceData(IntVect& rat, int ncomp)
{
  for (OrientationIter oitr; oitr; ++oitr) {
    Orientation ori = oitr();

    for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
      int n = face_data[ori][i].size();
      for (int j = 0; j < n; j++) {
        const Box& face_box = face_data[ori][i][j]->box();
        face_data[ori][i][j]->resize(face_box, ncomp);
      }
    }
  }
}

void CrseBndryAuxVar::loadFaceData(const Orientation ori,
                                   MultiFab& src,
                                   int srccomp,
                                   int destcomp,
                                   int numcomp)
{
  MultiFabCopyDescriptor mfcd;

  MultiFabId mfid = mfcd.RegisterMultiFab(&src);

  Vector< Vector<FillBoxId> > fbid(face_data[ori].size());
  for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
    int n = face_data[ori][i].size();
    fbid[i].resize(n);
    for (int j = 0; j < n; j++) {
      fbid[i][j] = mfcd.AddBox(mfid, face_data[ori][i][j]->box(), NULL,
                               fine_index[ori][i][j],
                               srccomp, destcomp, numcomp);
    }
  }

  mfcd.CollectData();

  for (int i = firstLocal(); isValid(i); i = nextLocal(i)) {
    int n = face_data[ori][i].size();
    for (int j = 0; j < n; j++) {
      mfcd.FillFab(mfid, fbid[i][j], *face_data[ori][i][j]);
    }
  }
}

void HypreMultiABec::vectorSetBoxValues(HYPRE_SStructVector x,
                                        int part,
                                        const Box& reg,
                                        const BoxArray& sgr,
                                        Real *vec)
{
  BL_PROFILE("HypreMultiABec::vectorSetBoxValues");

  if (sgr.size() > 0) {
    FArrayBox svecfab;
    for (int j = 0; j < sgr.size(); j++) {
      const Box& sreg = sgr[j];
      svecfab.resize(sreg);
      Real *svec = svecfab.dataPtr();
      for (IntVect v = sreg.smallEnd(); v <= sreg.bigEnd(); sreg.next(v)) {
        int is = sreg.index(v);
        int ir =  reg.index(v);
        svec[is] = vec[ir];
      }
      HYPRE_SStructVectorSetBoxValues(x, part, loV(sreg), hiV(sreg), 0, svec);
    }
  }
  else {
    HYPRE_SStructVectorSetBoxValues(x, part, loV(reg), hiV(reg), 0, vec);
  }
}

void HypreMultiABec::vectorGetBoxValues(HYPRE_SStructVector x,
                                        int part,
                                        const Box& reg,
                                        const BoxArray& sgr,
                                        FArrayBox& f, int fcomp)
{
  BL_PROFILE("HypreMultiABec::vectorGetBoxValues");

  BL_ASSERT(f.box() == reg);
  Real* vec = f.dataPtr(fcomp);
  if (sgr.size() > 0) {
    f.setVal<RunOn::Host>(0.0, fcomp);
    FArrayBox svecfab;
    for (int j = 0; j < sgr.size(); j++) {
      const Box& sreg = sgr[j];
      svecfab.resize(sreg);
      Real *svec = svecfab.dataPtr();
      HYPRE_SStructVectorGetBoxValues(x, part, loV(sreg), hiV(sreg), 0, svec);
      for (IntVect v = sreg.smallEnd(); v <= sreg.bigEnd(); sreg.next(v)) {
        int is = sreg.index(v);
        int ir =  reg.index(v);
        vec[ir] = svec[is];
      }
    }
  }
  else {
    HYPRE_SStructVectorGetBoxValues(x, part, loV(reg), hiV(reg), 0, vec);
  }
}

HypreMultiABec::HypreMultiABec(int _crse_level, int _fine_level,
                               int _solver_flag)
  : crse_level(_crse_level), fine_level(_fine_level),
    solver_flag(_solver_flag),
    geom(fine_level+1),
    grids(fine_level+1),
    dmap(fine_level+1),
    fine_ratio(fine_level+1),
    bd(fine_level+1),
    subgrids(fine_level+1),
    acoefs(fine_level+1),
    bcoefs(fine_level+1),
    SPa(fine_level+1),
    cintrp(fine_level+1),
    ederiv(fine_level+1),
    c_cintrp(fine_level+1),
    c_ederiv(fine_level+1),
    c_entry(fine_level+1),
    hgrid(NULL), stencil(NULL), graph(NULL),
    A(NULL), A0(NULL), b(NULL), x(NULL),
    sstruct_solver(NULL), solver(NULL), precond(NULL)
{
  ParmParse pp("hmabec");

  verbose = 1;       pp.query("v", verbose); pp.query("verbose", verbose);
  verbose_threshold = 0; pp.query("verbose_threshold", verbose_threshold);
  bho = 0;           pp.query("bho", bho);
  use_subgrids = 0;  pp.query("use_subgrids", use_subgrids);

  static int first = 1;
  if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
    first = 0;
    std::cout << "hmabec.bho               = " << bho << std::endl;
    std::cout << "hmabec.use_subgrids      = " << use_subgrids << std::endl;
    std::cout << "hmabec.verbose           = " << verbose << std::endl;
    std::cout << "hmabec.verbose_threshold = " << verbose_threshold << std::endl;
  }

  if (solver_flag == 100 || solver_flag == 102 || solver_flag == 104 ||
      solver_flag == 105 ||
      solver_flag == 150 || solver_flag == 151 || solver_flag == 153 ||
      solver_flag == 1002) {
    ObjectType = HYPRE_PARCSR;
  }
  else if (solver_flag == 101 || solver_flag == 103 ||
           solver_flag == 106 || solver_flag == 107 ||
           solver_flag == 152 ||
           solver_flag == 1003) {
    ObjectType = HYPRE_SSTRUCT;
  }
  else if (solver_flag == 108 || solver_flag == 109) {
    ObjectType = HYPRE_STRUCT;
  }
  else {
    std::cout << "HypreMultiABec: no such solver" << std::endl;
    exit(1);
  }

  int nparts = fine_level - crse_level + 1;

#if (AMREX_SPACEDIM == 1)

  // Hypre doesn't support 1D directly, so we use 2D Hypre with
  // the second dimension collapsed.
  // (SMG reduces to cyclic reduction in this case, so it's an exact solve.)
  // (PFMG will not work.)

  HYPRE_SStructGridCreate(MPI_COMM_WORLD, 2, nparts, &hgrid);

#else

  HYPRE_SStructGridCreate(MPI_COMM_WORLD, AMREX_SPACEDIM, nparts, &hgrid);

#endif
}

HypreMultiABec::~HypreMultiABec()
{
  HYPRE_SStructVectorDestroy(b);
  HYPRE_SStructVectorDestroy(x);

  HYPRE_SStructMatrixDestroy(A);
  HYPRE_SStructMatrixDestroy(A0);

  HYPRE_SStructGraphDestroy(graph);
  HYPRE_SStructStencilDestroy(stencil);
  HYPRE_SStructGridDestroy(hgrid);
}

void HypreMultiABec::addLevel(int             level,
                              const Geometry& _geom,
                              const BoxArray& _grids,
                              const DistributionMapping& _dmap,
                              IntVect         _fine_ratio)
{
  int part = level - crse_level;

  geom[level]  = _geom;
  grids[level] = _grids;
  dmap[level]  = _dmap;
  fine_ratio[level] = _fine_ratio;

#if (AMREX_SPACEDIM == 1)

  if (geom[level].isAnyPeriodic()) {
    BL_ASSERT(geom[level].isPeriodic(0));
    BL_ASSERT(geom[level].Domain().smallEnd(0) == 0);

    int is_periodic[2];
    is_periodic[0] = geom[level].period(0);
    is_periodic[1] = 0;
    BL_ASSERT(ispow2(is_periodic[0]));

    HYPRE_SStructGridSetPeriodic(hgrid, part, is_periodic);
  }

#else

  if (geom[level].isAnyPeriodic()) {
    int is_periodic[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      is_periodic[i] = 0;
      if (geom[level].isPeriodic(i)) {
        is_periodic[i] = geom[level].period(i);
        BL_ASSERT(ispow2(is_periodic[i]));
        BL_ASSERT(geom[level].Domain().smallEnd(i) == 0);
      }
    }
    HYPRE_SStructGridSetPeriodic(hgrid, part, is_periodic);
  }

#endif

  int myid = ParallelDescriptor::MyProc();

  subgrids[level].resize(grids[level].size());

  if (!use_subgrids || level == fine_level || grids[level+1].size() == 0) {
    for (int i = 0; i < grids[level].size(); i++) {
      if (dmap[level][i] == myid) {
        HYPRE_SStructGridSetExtents(hgrid, part,
                                    loV(grids[level][i]),
                                    hiV(grids[level][i]));
      }
    }
  }
  else {
    BoxArray mask = grids[level+1];
    mask.coarsen(fine_ratio[level]);
    for (int i = 0; i < grids[level].size(); i++) {
      if (dmap[level][i] == myid) {
        subgrids[level][i] = amrex::complementIn(grids[level][i], mask);
#if 0
        std::cout << "Coarse level, grid " << i << " subgrids are "
             << subgrids[level][i] << std::endl;
        std::cout << "Box numPts = " << grids[level][i].numPts()
             << "Subgrids numPts = " << subgrids[level][i].numPts() << std::endl;
#endif
        for (int j = 0; j < subgrids[level][i].size(); j++) {
          HYPRE_SStructGridSetExtents(hgrid, part,
                                      loV(subgrids[level][i][j]),
                                      hiV(subgrids[level][i][j]));
        }
      }
    }
  }

  // All variables are cell-centered
  HYPRE_SStructVariable vars[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};
  HYPRE_SStructGridSetVariables(hgrid, part, 1, vars);
}

static void
TransverseInterpolant(AuxVarBox& cintrp, const Mask& msk,
                      const Box& reg, const Box& creg,
                      AMREX_D_DECL(const IntVect& rat, const IntVect& vj1, const IntVect& vk1),
                      AMREX_D_DECL(const IntVect& ve,  const IntVect& vjr, const IntVect& vkr),
                      AMREX_D_DECL(int idir, int jdir, int kdir),
                      int clevel)
{
  for (IntVect vc = creg.smallEnd(); vc <= creg.bigEnd(); creg.next(vc)) {
    IntVect vf = rat * vc;
    vf[idir] = reg.smallEnd(idir); // same as bigEnd(idir)
    Box face(vf, vf + ve);
    if (msk(vf) == RadBndryData::not_covered) {
#if (0)
      // force piecewise constant interpolation for debugging:
      for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
        cintrp(v).push(clevel, vc,     1.0);
      }
#elif (AMREX_SPACEDIM == 1)
      cintrp(vf).push(clevel, vc, 1.0);
#elif (AMREX_SPACEDIM == 2)
      if (msk(vf-vj1) != RadBndryData::not_covered &&
          msk(vf+vjr) == RadBndryData::not_covered) {
        // low direction not available, use linear interp upwards:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          cintrp(v).push(clevel, vc+vj1, xx);
          cintrp(v).push(clevel, vc,     1.0 - xx);
        }
      }
      else if (msk(vf-vj1) == RadBndryData::not_covered &&
               msk(vf+vjr) == RadBndryData::not_covered) {
        // use piecewise quadratic interpolation whenever possible:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          cintrp(v).push(clevel, vc+vj1, 0.5*xx*(xx+1));
          cintrp(v).push(clevel, vc,     1.0-xx*xx);
          cintrp(v).push(clevel, vc-vj1, 0.5*xx*(xx-1));
        }
      }
      else if (msk(vf-vj1) == RadBndryData::not_covered &&
               msk(vf+vjr) != RadBndryData::not_covered) {
        // high direction not available, use linear interp downwards:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          cintrp(v).push(clevel, vc,     1.0 + xx);
          cintrp(v).push(clevel, vc-vj1, -xx);
        }
      }
      else {
        // neither direction available, drop back to piecewise const:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          cintrp(v).push(clevel, vc,     1.0);
        }
        //amrex::Error("Case not implemented");
      }
#elif (AMREX_SPACEDIM == 3)

      // First do the jdir direction, including piecewise-constant term:

      if (msk(vf-vj1) != RadBndryData::not_covered &&
          msk(vf+vjr) == RadBndryData::not_covered) {
        // low direction not available, use linear interp upwards:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          cintrp(v).push(clevel, vc+vj1, xx);
          cintrp(v).push(clevel, vc,     1.0 - xx);
        }
      }
      else if (msk(vf-vj1) == RadBndryData::not_covered &&
               msk(vf+vjr) == RadBndryData::not_covered) {
        // use piecewise quadratic interpolation whenever possible:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          cintrp(v).push(clevel, vc+vj1, 0.5*xx*(xx+1));
          cintrp(v).push(clevel, vc,     1.0-xx*xx);
          cintrp(v).push(clevel, vc-vj1, 0.5*xx*(xx-1));
        }
      }
      else if (msk(vf-vj1) == RadBndryData::not_covered &&
               msk(vf+vjr) != RadBndryData::not_covered) {
        // high direction not available, use linear interp downwards:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          cintrp(v).push(clevel, vc,     1.0 + xx);
          cintrp(v).push(clevel, vc-vj1, -xx);
        }
      }
      else {
        // neither direction available, drop back to piecewise const:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          cintrp(v).push(clevel, vc,     1.0);
        }
      }

      // Then add on contributions from the kdir direction:

      if (msk(vf-vk1) != RadBndryData::not_covered &&
          msk(vf+vkr) == RadBndryData::not_covered) {
        // low direction not available, use linear interp upwards:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real yy = (v[kdir] - vf[kdir] - 0.5 * ve[kdir]) / rat[kdir];
          cintrp(v).push(clevel, vc+vk1, yy);
          cintrp(v).push(clevel, vc,    -yy);
        }
      }
      else if (msk(vf-vk1) == RadBndryData::not_covered &&
               msk(vf+vkr) == RadBndryData::not_covered) {
        // use piecewise quadratic interpolation whenever possible:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real yy = (v[kdir] - vf[kdir] - 0.5 * ve[kdir]) / rat[kdir];
          cintrp(v).push(clevel, vc+vk1, 0.5*yy*(yy+1));
          cintrp(v).push(clevel, vc,     -yy*yy);
          cintrp(v).push(clevel, vc-vk1, 0.5*yy*(yy-1));
        }
      }
      else if (msk(vf-vk1) == RadBndryData::not_covered &&
               msk(vf+vkr) != RadBndryData::not_covered) {
        // high direction not available, use linear interp downwards:
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real yy = (v[kdir] - vf[kdir] - 0.5 * ve[kdir]) / rat[kdir];
          cintrp(v).push(clevel, vc,      yy);
          cintrp(v).push(clevel, vc-vk1, -yy);
        }
      }
      else {
        // neither direction available, drop back to piecewise const:
        // (do nothing, no need to explicitly add a zero contribution here)
      }

      // Finally add in cross derivative terms:

      if (msk(vf-vj1-vk1) == RadBndryData::not_covered &&
          msk(vf-vj1+vkr) == RadBndryData::not_covered &&
          msk(vf+vjr-vk1) == RadBndryData::not_covered &&
          msk(vf+vjr+vkr) == RadBndryData::not_covered) {
        for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
          Real xx = (v[jdir] - vf[jdir] - 0.5 * ve[jdir]) / rat[jdir];
          Real yy = (v[kdir] - vf[kdir] - 0.5 * ve[kdir]) / rat[kdir];
          cintrp(v).push(clevel, vc-vj1-vk1,  0.25*xx*yy);
          cintrp(v).push(clevel, vc-vj1+vk1, -0.25*xx*yy);
          cintrp(v).push(clevel, vc+vj1-vk1, -0.25*xx*yy);
          cintrp(v).push(clevel, vc+vj1+vk1,  0.25*xx*yy);
        }
      }
#endif
    }
  }
}

static void
NormalDerivative(AuxVarBox& ederiv, AuxVarBox& cintrp,
                 const Mask& msk, const Box& reg,
                 const IntVect& vin, Real h, int r, int bho, int flevel)
{
  if (bho == 1) {
    Real efacb = 8.0 / (h * (1 + r) * (3 + r));
    Real efac1 = (r - 3) / (h * (1 + r));
    Real efac2 = (1 - r) / (h * (3 + r));
    IntVect vi2 = 2 * vin;
    for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
      if (msk(v) == RadBndryData::not_covered) {
        ederiv(v).push(&cintrp(v),    efacb);
        ederiv(v).push(flevel, v-vin, efac1);
        ederiv(v).push(flevel, v-vi2, efac2);
      }
    }
  }
  else {
    Real efac = 2.0 / (h * (1 + r)); // normal derivative factor
    for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
      if (msk(v) == RadBndryData::not_covered) {
        ederiv(v).push(&cintrp(v),     efac);
        ederiv(v).push(flevel, v-vin, -efac);
      }
    }
  }
}

void HypreMultiABec::buildMatrixStructure()
{
  // Build MultiFabs first, so that distribution maps are available
  // when building the graph below:

  for (int level = crse_level; level <= fine_level; level++) {
    int ncomp=1;
    int ngrow=0;
    acoefs[level].reset(new MultiFab(grids[level], dmap[level], ncomp, ngrow));
    acoefs[level]->setVal(0.0);

    bcoefs[level].reset(new Array<MultiFab, AMREX_SPACEDIM>);

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      BoxArray edge_boxes(grids[level]);
      edge_boxes.surroundingNodes(i);
      (*bcoefs[level])[i].define(edge_boxes, dmap[level], ncomp, ngrow);
      (*bcoefs[level])[i].setVal(0.0);
    }
  }

  // This can be done now that addLevel has been called for each level:

  HYPRE_SStructGridAssemble(hgrid);

  // Setup stencils:

#if (AMREX_SPACEDIM == 1)
  // if we were really 1D:
/*
  int offsets[3][1] = {{ 0}
                       {-1},
                       { 1}};
*/
  // fake 1D as a 2D problem:
  int offsets[3][2] = {{ 0,  0},
                       {-1,  0},
                       { 1,  0}};
#elif (AMREX_SPACEDIM == 2)
  int offsets[5][2] = {{ 0,  0},
                       {-1,  0},
                       { 1,  0},
                       { 0, -1},
                       { 0,  1}};
#elif (AMREX_SPACEDIM == 3)
  int offsets[7][3] = {{ 0,  0,  0},
                       {-1,  0,  0},
                       { 1,  0,  0},
                       { 0, -1,  0},
                       { 0,  1,  0},
                       { 0,  0, -1},
                       { 0,  0,  1}};
#endif

  BL_ASSERT(stencil == NULL);
#if (AMREX_SPACEDIM == 1)
  HYPRE_SStructStencilCreate(2, 3, &stencil);
#else
  HYPRE_SStructStencilCreate(AMREX_SPACEDIM, 2 * AMREX_SPACEDIM + 1, &stencil);
#endif

  for (int i = 0; i < 2 * AMREX_SPACEDIM + 1; i++) {
    HYPRE_SStructStencilSetEntry(stencil, i, offsets[i], 0);
  }

  BL_ASSERT(graph == NULL);
  HYPRE_SStructGraphCreate(MPI_COMM_WORLD, hgrid, &graph);
  HYPRE_SStructGraphSetObjectType(graph, ObjectType);

  for (int level = crse_level; level <= fine_level; level++) {
    int part = level - crse_level;
    HYPRE_SStructGraphSetStencil(graph, part, 0, stencil);
  }

  // Add non-stencil entries to the graph here:

  for (int level = crse_level + 1; level <= fine_level; level++) {
    int part = level - crse_level;
    cintrp[level].reset(new BndryAuxVar(grids[level], dmap[level], BndryAuxVar::GHOST));
    ederiv[level].reset(new BndryAuxVar(grids[level], dmap[level],  BndryAuxVar::GHOST));
    BndryAuxVar entry(grids[level], dmap[level], BndryAuxVar::INTERIOR);
    IntVect rat = fine_ratio[level-1];

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      IntVect vin = amrex::BASISV(idir);
      vin = (ori.isLow() ? -vin : vin); // outward normal unit vector
      Real h = geom[level].CellSize(idir); // normal fine grid spacing
      IntVect ve; // default constructor initializes to zero
#if (AMREX_SPACEDIM >= 2)
      int jdir = (idir + 1) % AMREX_SPACEDIM;
      IntVect vj1 = amrex::BASISV(jdir); // tangential unit vector
      IntVect vjr = rat * vj1;
      ve += (vjr - vj1);
#endif
#if (AMREX_SPACEDIM == 3)
      int kdir = (idir + 2) % 3;
      IntVect vk1 = amrex::BASISV(kdir);
      IntVect vkr = rat * vk1;
      ve += (vkr - vk1);
#endif
      for (int i = cintrp[level]->firstLocal(); cintrp[level]->isValid(i);
           i = cintrp[level]->nextLocal(i)) {
        Box reg = amrex::adjCell(grids[level][i], ori);
        Box creg = amrex::coarsen(reg, rat); // coarse adjacent cells
        const Mask &msk = bd[level]->bndryMasks(ori,i);

        TransverseInterpolant((*cintrp[level])(ori,i), msk, reg, creg,
                              AMREX_D_DECL(rat,  vj1,  vk1),
                              AMREX_D_DECL(ve,   vjr,  vkr),
                              AMREX_D_DECL(idir, jdir, kdir),
                              level-1);

        NormalDerivative((*ederiv[level])(ori,i),
                         (*cintrp[level])(ori,i),
                         msk, reg, vin, h, rat[idir], bho, level);

        // cintrp and ederiv are now complete.  entry will be done in
        // draft form to establish the graph connections, but must be
        // done again later when the edge coefficients are available.

        reg.shift(-vin); // fine interior cells
        for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
          if (msk(v+vin) == RadBndryData::not_covered) {
            // value not important, since coefficient not known.
              entry(ori,i)(v).push(&(*ederiv[level])(ori,i)(v+vin), 1.0);
          }
        }
      }
    }

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      IntVect vin = amrex::BASISV(idir);
      vin = (ori.isLow() ? -vin : vin); // outward normal unit vector
      for (int i = entry.firstLocal(); entry.isValid(i);
           i = entry.nextLocal(i)) {
        Box reg = amrex::adjCell(grids[level][i], ori);
        reg.shift(-vin); // fine interior cells
//        const Mask &msk = bd[level]->bndryMasks(ori,i);
        for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
#if (0 && !defined(NDEBUG))
          if (msk(v+vin) == RadBndryData::not_covered &&
              entry(ori,i)(v).secondary()) {
            std::cout << v << " is secondary in orientation " << ori
                 << " on processor " << ParallelDescriptor::MyProc()
                 << std::endl;
          }
#endif
          // Even if this entry is covered, it could have a
          // not_covered secondary:
          if (!entry(ori,i)(v).empty() &&
              !entry(ori,i)(v).secondary()) {
            entry(ori,i)(v).collapse();
            Vector<int> levels;
            Vector<IntVect> cells;
            int retval = entry(ori,i)(v).get_locations(levels, cells);
            BL_ASSERT(retval == 0);
            for (int j = 0; j < levels.size(); j++) {
              // eliminate stencil-like connections:
              int not_stencil = 1;
              if (levels[j] == level) {
                IntVect d = cells[j] - v;
                for (int k = 0; k < 2 * AMREX_SPACEDIM + 1; k++) {
                  if (d == IntVect(offsets[k])) {
                    not_stencil = 0;
                  }
                }
              }
              if (not_stencil) {
                int to_part = levels[j] - crse_level;
                HYPRE_SStructGraphAddEntries(graph,
                                             part,    getV1(v),        0,
                                             to_part, getV2(cells[j]), 0);
              }
            }
          }
        }
      }
    }

    // Now add the graph entries seen by the coarse cells adjacent
    // to the coarse-fine interface.  These are averages of ederiv
    // over the fine faces making up each coarse face.  Since we
    // have to do this from the processor owning the coarse grid,
    // we recompute information using CrseBndryAuxVar

    const BoxArray& f_fgrids(grids[level]);
    const BoxArray& c_cgrids(grids[level-1]);
    BoxArray f_cgrids(c_cgrids);
    f_cgrids.refine(rat);
    BoxArray c_fgrids(f_fgrids);
    c_fgrids.coarsen(rat);

    c_cintrp[level].reset(new CrseBndryAuxVar(f_cgrids, dmap[level-1],
                                              f_fgrids, BndryAuxVar::GHOST));
    c_ederiv[level].reset(new CrseBndryAuxVar(*c_cintrp[level],
                                            BndryAuxVar::GHOST));
    c_entry[level].reset (new CrseBndryAuxVar(c_cgrids, c_fgrids,
                                            *c_cintrp[level],
                                            BndryAuxVar::EXTERIOR));

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      IntVect vin = amrex::BASISV(idir);
      vin = (ori.isLow() ? -vin : vin); // outward normal unit vector
      Real h = geom[level].CellSize(idir); // normal fine grid spacing
      IntVect ve; // default constructor initializes to zero
#if (AMREX_SPACEDIM >= 2)
      int jdir = (idir + 1) % AMREX_SPACEDIM;
      IntVect vj1 = amrex::BASISV(jdir); // tangential unit vector
      IntVect vjr = rat * vj1;
      ve += (vjr - vj1);
#endif
#if (AMREX_SPACEDIM == 3)
      int kdir = (idir + 2) % 3;
      IntVect vk1 = amrex::BASISV(kdir);
      IntVect vkr = rat * vk1;
      ve += (vkr - vk1);
#endif
      for (int i = c_cintrp[level]->firstLocal(); c_cintrp[level]->isValid(i);
         i = c_cintrp[level]->nextLocal(i)) {
         for (int j = 0; j < (*c_cintrp[level]).size(ori,i); j++) {
           const Box& reg = (*c_cintrp[level])(ori,i,j).box(); // adjacent cells
           const Box& creg = (*c_entry[level])(ori,i,j).box(); // adjacent cells
           const Mask &msk = c_cintrp[level]->mask(ori,i,j); // fine mask

         TransverseInterpolant((*c_cintrp[level])(ori,i,j), msk, reg, creg,
                                AMREX_D_DECL(rat,  vj1,  vk1),
                                AMREX_D_DECL(ve,   vjr,  vkr),
                                AMREX_D_DECL(idir, jdir, kdir),
                                level-1);

         NormalDerivative((*c_ederiv[level])(ori,i,j),
                          (*c_cintrp[level])(ori,i,j),
                          msk, reg, vin, h, rat[idir], bho, level);

          // c_cintrp and c_ederiv are now complete.  c_entry will be done
          // in draft form to establish the graph connections, but must be
          // done again later when the edge coefficients are available.

          for (IntVect vc = creg.smallEnd(); vc <= creg.bigEnd(); creg.next(vc)) {
            IntVect vf = rat * vc;
            vf[idir] = reg.smallEnd(idir); // same as bigEnd(idir)
            Box face(vf, vf + ve);
            if (msk(vf) == RadBndryData::not_covered) {
              for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
                // value not important, since coefficient not known.
                  (*c_entry[level])(ori,i,j)(vc)
                      .push(&(*c_ederiv[level])(ori,i,j)(v), 1.0);
              }
            }
          }
        }
      }
    }

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      for (int i = c_cintrp[level]->firstLocal(); c_cintrp[level]->isValid(i);
           i = c_cintrp[level]->nextLocal(i)) {
        for (int j = 0; j < (*c_cintrp[level]).size(ori,i); j++) {
          const Box& reg = (*c_cintrp[level])(ori,i,j).box(); // adjacent cells
          const Box& creg = (*c_entry[level])(ori,i,j).box(); // adjacent cells
          const Mask &msk = c_cintrp[level]->mask(ori,i,j); // fine mask
          for (IntVect vc = creg.smallEnd(); vc <= creg.bigEnd(); creg.next(vc)) {
            IntVect vf = rat * vc;
            vf[idir] = reg.smallEnd(idir); // same as bigEnd(idir)
            // Unlike fine entry, it should not be possible for this
            // entry to be covered but have a not_covered secondary:
            if (msk(vf) == RadBndryData::not_covered &&
                !(*c_entry[level])(ori,i,j)(vc).secondary()) {
              (*c_entry[level])(ori,i,j)(vc).collapse();
              Vector<int> levels;
              Vector<IntVect> cells;
              int retval = (*c_entry[level])(ori,i,j)(vc)
                .get_locations(levels, cells);
              BL_ASSERT(retval == 0);
              for (int jj = 0; jj < levels.size(); jj++) {
                // eliminate stencil-like connections:
                int not_stencil = 1;
                if (levels[jj] == level-1) {
                  IntVect d = cells[jj] - vc;
                  for (int k = 0; k < 2 * AMREX_SPACEDIM + 1; k++) {
                    if (d == IntVect(offsets[k])) {
                      not_stencil = 0;
                    }
                  }
                }
                if (not_stencil) {
                  int to_part = levels[jj] - crse_level;
                  HYPRE_SStructGraphAddEntries(graph,
                                               part-1,  getV1(vc),        0,
                                               to_part, getV2(cells[jj]), 0);
                }
              }
            }
          }
        }
      }
    }
  }

  HYPRE_SStructGraphAssemble(graph);

  BL_ASSERT(A == NULL);
  HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A);
  HYPRE_SStructMatrixSetObjectType(A, ObjectType);
  //HYPRE_StructMatrixSetSymmetric(A, 1);
  //HYPRE_StructMatrixSetNumGhost(A, A_num_ghost);
  HYPRE_SStructMatrixInitialize(A);

  BL_ASSERT(A0 == NULL);
  HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A0);
  HYPRE_SStructMatrixSetObjectType(A0, ObjectType);
  //HYPRE_StructMatrixSetSymmetric(A0, 1);
  //HYPRE_StructMatrixSetNumGhost(A0, A_num_ghost);
  HYPRE_SStructMatrixInitialize(A0);

  BL_ASSERT(b == NULL);
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, hgrid, &b);
  HYPRE_SStructVectorSetObjectType(b, ObjectType);

  BL_ASSERT(x == NULL);
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, hgrid, &x);
  HYPRE_SStructVectorSetObjectType(x, ObjectType);

  HYPRE_SStructVectorInitialize(b);
  HYPRE_SStructVectorInitialize(x);

  // According to Rob, the following is necessary in some cases before
  // we call the solver setup routine (and we may do that before loading
  // the vectors with data):

  HYPRE_SStructVectorAssemble(b);
  HYPRE_SStructVectorAssemble(x);
}

void HypreMultiABec::setScalars(Real Alpha, Real Beta)
{
  alpha = Alpha;
  beta  = Beta;
}

void HypreMultiABec::aCoefficients(int level, const MultiFab &a)
{
  BL_PROFILE("HypreMultiABec::aCoefficients");

  BL_ASSERT( a.ok() );
  BL_ASSERT( a.boxArray() == acoefs[level]->boxArray() );
  MultiFab::Copy(*acoefs[level], a, 0, 0, 1, 0);
}

void HypreMultiABec::bCoefficients(int level, const MultiFab &b, int dir)
{
  BL_PROFILE("HypreMultiABec::bCoefficients");

  BL_ASSERT( b.ok() );
  BL_ASSERT( b.boxArray() == (*bcoefs[level])[dir].boxArray() );
  MultiFab::Copy((*bcoefs[level])[dir], b, 0, 0, 1, 0);
}

void HypreMultiABec::SPalpha(int level, const MultiFab& a)
{
  BL_ASSERT( a.ok() );
  if (! SPa[level]) {
    SPa[level].reset(new MultiFab(grids[level],dmap[level],1,0));
    BL_ASSERT( a.boxArray() == SPa[level]->boxArray() );
    BL_ASSERT( a.DistributionMap() == SPa[level]->DistributionMap() );
  }
  MultiFab::Copy(*SPa[level], a, 0, 0, 1, 0);
}

void HypreMultiABec::hmac (const Box& bx,
                           Array4<GpuArray<Real, 2 * AMREX_SPACEDIM + 1>> const& mat,
                           Array4<Real const> const& a,
                           Real alpha)
{
    BL_PROFILE("HypreMultiABec::hmac");

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        if (alpha == 0.e0_rt) {
            mat(i,j,k)[0] = 0.e0_rt;
        }
        else {
            mat(i,j,k)[0] = alpha * a(i,j,k);
        }
    });

    Gpu::synchronize();
}

void HypreMultiABec::hmbc (const Box& bx,
                           Array4<GpuArray<Real, 2 * AMREX_SPACEDIM + 1>> const& mat,
                           Array4<Real const> const& b,
                           Real beta, const Real* dx, int n)
{
    BL_PROFILE("HypreMultiABec::hmbc");

    if (n == 0) {

        const Real fac = beta / (dx[0] * dx[0]);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            mat(i,j,k)[0] += fac * (b(i,j,k) + b(i+1,j,k));
            mat(i,j,k)[1] = -fac * b(i,j,k);
            mat(i,j,k)[2] = -fac * b(i+1,j,k);
        });

    }
    else if (n == 1) {

        const Real fac = beta / (dx[1] * dx[1]);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            mat(i,j,k)[0] += fac * (b(i,j,k) + b(i,j+1,k));
            mat(i,j,k)[3] = -fac * b(i,j,k);
            mat(i,j,k)[4] = -fac * b(i,j+1,k);
        });

    }
    else {

        const Real fac = beta / (dx[2] * dx[2]);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            mat(i,j,k)[0] += fac * (b(i,j,k) + b(i,j,k+1));
            mat(i,j,k)[5] = -fac * b(i,j,k);
            mat(i,j,k)[6] = -fac * b(i,j,k+1);
        });

    }

    Gpu::synchronize();
}

void
HypreMultiABec::hmmat (const Box& bx,
                       Array4<GpuArray<Real, 2 * AMREX_SPACEDIM + 1>> const& mat,
                       int cdir, int bct, int bho, Real bcl,
                       Array4<int const> const& mask,
                       Array4<Real const> const& b,
                       Real beta, const Real* dx)
{
    bool xlo = false;
    bool ylo = false;
    bool zlo = false;

    bool xhi = false;
    bool yhi = false;
    bool zhi = false;

    Real h;

    if (AMREX_SPACEDIM == 1) {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            xhi = true;
            h = dx[0];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }
    else if (AMREX_SPACEDIM == 2) {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 2) {
            xhi = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            ylo = true;
            h = dx[1];
        }
        else if (cdir == 3) {
            yhi = true;
            h = dx[1];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }
    else {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 3) {
            xhi = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            ylo = true;
            h = dx[1];
        }
        else if (cdir == 4) {
            yhi = true;
            h = dx[1];
        }
        else if (cdir == 2) {
            zlo = true;
            h = dx[2];
        }
        else if (cdir == 5) {
            zhi = true;
            h = dx[2];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }

    const Real fac = beta / (h * h);

    Real bfm, bfv;
    Real bfm2, h2, th2;

    if (bct == AMREX_LO_DIRICHLET) {

        if (bho >= 1) {

            h2 = 0.5e0_rt * h;
            th2 = 3.e0_rt * h2;
            bfm = fac * (th2 - bcl) / (bcl + h2) - fac;
            bfm2 = fac * (bcl - h2) / (bcl + th2);

        }
        else {

            bfv = (beta / h) / (0.5e0_rt * h + bcl);
            bfm = bfv - fac;

        }

    }
    else if (bct == AMREX_LO_NEUMANN) {

        bfm = -fac;
        bfm2 = 0.e0_rt;

    }
    else {

        amrex::Error("hmmat: unsupported boundary type");

    }

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        if (mask.contains(i-1,j,k)) {

            if (xlo && mask(i-1,j,k) > 0) {

                mat(i,j,k)[0] += bfm * b(i,j,k);
                mat(i,j,k)[1] = 0.e0_rt;
                if (bho >= 1) {
                    mat(i,j,k)[2] += bfm2 * b(i,j,k);
                }

            }

        }
        else if (mask.contains(i+1,j,k)) {

            if (xhi && mask(i+1,j,k) > 0) {

                mat(i,j,k)[0] += bfm * b(i+1,j,k);
                mat(i,j,k)[2] = 0.e0_rt;
                if (bho >= 1) {
                    mat(i,j,k)[1] += bfm2 * b(i+1,j,k);
                }

            }

        }
        else if (mask.contains(i,j-1,k)) {

            if (ylo && mask(i,j-1,k) > 0) {

                mat(i,j,k)[0] += bfm * b(i,j,k);
                mat(i,j,k)[3] = 0.e0_rt;
                if (bho >= 1) {
                    mat(i,j,k)[4] += bfm2 * b(i,j,k);
                }

            }

        }
        else if (mask.contains(i,j+1,k)) {

            if (yhi && mask(i,j+1,k) > 0) {

                mat(i,j,k)[0] += bfm * b(i,j+1,k);
                mat(i,j,k)[4] = 0.e0_rt;
                if (bho >= 1) {
                    mat(i,j,k)[3] += bfm2 * b(i,j+1,k);
                }

            }

        }
        else if (mask.contains(i,j,k-1)) {

            if (zlo && mask(i,j,k-1) > 0) {

                mat(i,j,k)[0] += bfm * b(i,j,k);
                mat(i,j,k)[5] = 0.e0_rt;
                if (bho >= 1) {
                    mat(i,j,k)[6] += bfm2 * b(i,j,k);
                }

            }

        }
        else if (mask.contains(i,j,k+1)) {

            if (zhi && mask(i,j,k+1) > 0) {

                mat(i,j,k)[0] += bfm * b(i,j,k+1);
                mat(i,j,k)[6] = 0.e0_rt;
                if (bho >= 1) {
                    mat(i,j,k)[5] += bfm2 * b(i,j,k+1);
                }

            }

        }
    });

    Gpu::synchronize();
}

void
HypreMultiABec::hmmat3 (const Box& bx,
                        int ori_lo, int idir,
                        Array4<GpuArray<Real, 2 * AMREX_SPACEDIM + 1>> const& mat,
                        int cdir, int bctype,
                        Array4<int const> const& tf,
                        int bho, Real bcl,
                        Array4<int const> const& mask,
                        Array4<Real const> const& b,
                        Real beta, const Real* dx,
                        const GeometryData& geomdata, Real c,
                        Array4<Real const> const& spa)
{
    BL_PROFILE("HypreMultiABec::hmmat3");

    bool xlo = false;
    bool ylo = false;
    bool zlo = false;

    bool xhi = false;
    bool yhi = false;
    bool zhi = false;

    Real h;

    if (AMREX_SPACEDIM == 1) {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            xhi = true;
            h = dx[0];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }
    else if (AMREX_SPACEDIM == 2) {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 2) {
            xhi = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            ylo = true;
            h = dx[1];
        }
        else if (cdir == 3) {
            yhi = true;
            h = dx[1];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }
    else {

        if (cdir == 0) {
            xlo = true;
            h = dx[0];
        }
        else if (cdir == 3) {
            xhi = true;
            h = dx[0];
        }
        else if (cdir == 1) {
            ylo = true;
            h = dx[1];
        }
        else if (cdir == 4) {
            yhi = true;
            h = dx[1];
        }
        else if (cdir == 2) {
            zlo = true;
            h = dx[2];
        }
        else if (cdir == 5) {
            zhi = true;
            h = dx[2];
        }
        else {
            amrex::Error("Unknown cdir");
        }

    }

    const Real fac = beta / (h * h);

    // The -fac * b(i,j,k) term applied to the matrix diagonal is the contribution
    // from the interior stencil which must be removed at the boundary.

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        Real r;
        face_metric(i, j, k, bx.loVect()[0], bx.hiVect()[0], geomdata, idir, ori_lo, r);

        Real bfm, bfv;
        Real bfm2, h2, th2;
        int bct;

        if (mask.contains(i-1,j,k)) {

            if (xlo && mask(i-1,j,k) > 0) {

                if (bctype == -1) {
                    bct = tf(i-1,j,k);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {

                    if (bho >= 1) {
                        h2 = 0.5e0_rt * h;
                        th2 = 3.e0_rt * h2;
                        bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k);
                        bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k);
                    }
                    else  {
                        bfv = (beta / h) / (0.5e0_rt * h + bcl);
                        bfm = bfv * b(i,j,k);
                    }

                }
                else if (bct == AMREX_LO_NEUMANN) {

                    bfm  = 0.e0_rt;
                    bfm2 = 0.e0_rt;

                }
                else if (bct == AMREX_LO_MARSHAK) {

                    bfv = 2.e0_rt * c * beta * r / h;

                    if (bho >= 1) {
                        bfm  =  0.375e0_rt * bfv;
                        bfm2 = -0.125e0_rt * bfv;
                    }
                    else {
                        bfm = 0.25e0_rt * bfv;
                    }

                }
                else if (bct == AMREX_LO_SANCHEZ_POMRANING) {

                    bfv = 2.e0_rt * c * beta * r / h;

                    if (bho >= 1) {
                        bfm  =  1.5e0_rt * spa(i,j,k) * bfv;
                        bfm2 = -0.5e0_rt * spa(i,j,k) * bfv;
                    }
                    else {
                        bfm = spa(i,j,k) * bfv;
                    }

                }
#ifndef AMREX_USE_GPU
                else {

                    amrex::Error("hmmat3: unsupported boundary type");

                }
#endif

                mat(i,j,k)[0] += bfm - fac * b(i,j,k);
                mat(i,j,k)[1] = 0.e0_rt;
                if (bho >= 1) {
                    mat(i,j,k)[2] += bfm2;
                }

            }

        }
        else if (mask.contains(i+1,j,k)) {

            if (xhi && mask(i+1,j,k) > 0) {

                if (bctype == -1) {
                    bct = tf(i+1,j,k);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {

                    if (bho >= 1) {
                        h2 = 0.5e0_rt * h;
                        th2 = 3.e0_rt * h2;
                        bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i+1,j,k);
                        bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i+1,j,k);
                    }
                    else {
                        bfv = (beta / h) / (0.5e0_rt * h + bcl);
                        bfm = bfv * b(i+1,j,k);
                    }

                }
                else if (bct == AMREX_LO_NEUMANN) {

                    bfm  = 0.e0_rt;
                    bfm2 = 0.e0_rt;

                }
                else if (bct == AMREX_LO_MARSHAK) {

                    bfv = 2.e0_rt * c * beta * r / h;

                    if (bho >= 1) {
                        bfm  =  0.375e0_rt * bfv;
                        bfm2 = -0.125e0_rt * bfv;
                    }
                    else {
                        bfm = 0.25e0_rt * bfv;
                    }

                }
                else if (bct == AMREX_LO_SANCHEZ_POMRANING) {

                    bfv = 2.e0_rt * c * beta * r / h;

                    if (bho >= 1) {
                        bfm  =  1.5e0_rt * spa(i,j,k) * bfv;
                        bfm2 = -0.5e0_rt * spa(i,j,k) * bfv;
                    }
                    else {
                        bfm = spa(i,j,k) * bfv;
                    }

                }
#ifndef AMREX_USE_GPU
                else {

                    amrex::Error("hmmat3: unsupported boundary type");

                }
#endif

                mat(i,j,k)[0] += bfm - fac * b(i+1,j,k);
                mat(i,j,k)[2] = 0.e0_rt;
                if (bho >= 1) {
                    mat(i,j,k)[1] += bfm2;
                }

            }

        }
        else if (mask.contains(i,j-1,k)) {

            if (ylo && mask(i,j-1,k) > 0) {

                if (bctype == -1) {
                    bct = tf(i,j-1,k);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {

                    if (bho >= 1) {
                        h2 = 0.5e0_rt * h;
                        th2 = 3.e0_rt * h2;
                        bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k);
                        bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k);
                    }
                    else {
                        bfv = (beta / h) / (0.5e0_rt * h + bcl);
                        bfm = bfv * b(i,j,k);
                    }

                }
                else if (bct == AMREX_LO_NEUMANN) {

                    bfm  = 0.e0_rt;
                    bfm2 = 0.e0_rt;

                }
                else if (bct == AMREX_LO_MARSHAK) {

                    bfv = 2.e0_rt * c * beta * r / h;

                    if (bho >= 1) {
                        bfm  =  0.375e0_rt * bfv;
                        bfm2 = -0.125e0_rt * bfv;
                    }
                    else {
                        bfm = 0.25e0_rt * bfv;
                    }

                }
                else if (bct == AMREX_LO_SANCHEZ_POMRANING) {

                    bfv = 2.e0_rt * c * beta * r / h;

                    if (bho >= 1) {
                        bfm  =  1.5e0_rt * spa(i,j,k) * bfv;
                        bfm2 = -0.5e0_rt * spa(i,j,k) * bfv;
                    }
                    else {
                        bfm = spa(i,j,k) * bfv;
                    }

                }
#ifndef AMREX_USE_GPU
                else {

                    amrex::Error("hmmat3: unsupported boundary type");

                }
#endif

                mat(i,j,k)[0] += bfm - fac * b(i,j,k);
                mat(i,j,k)[3] = 0.e0_rt;
                if (bho >= 1) {
                    mat(i,j,k)[4] += bfm2;
                }

            }

        }
        else if (mask.contains(i,j+1,k)) {

            if (yhi && mask(i,j+1,k) > 0) {

                if (bctype == -1) {
                    bct = tf(i,j+1,k);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {

                    if (bho >= 1) {
                        h2 = 0.5e0_rt * h;
                        th2 = 3.e0_rt * h2;
                        bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j+1,k);
                        bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j+1,k);
                    }
                    else {
                        bfv = (beta / h) / (0.5e0_rt * h + bcl);
                        bfm = bfv * b(i,j+1,k);
                    }

                }
                else if (bct == AMREX_LO_NEUMANN) {

                    bfm  = 0.e0_rt;
                    bfm2 = 0.e0_rt;

                }
                else if (bct == AMREX_LO_MARSHAK) {

                    bfv = 2.e0_rt * c * beta * r / h;

                    if (bho >= 1) {
                        bfm  =  0.375e0_rt * bfv;
                        bfm2 = -0.125e0_rt * bfv;
                    }
                    else {
                        bfm = 0.25e0_rt * bfv;
                    }

                }
                else if (bct == AMREX_LO_SANCHEZ_POMRANING) {

                    bfv = 2.e0_rt * c * beta * r / h;

                    if (bho >= 1) {
                        bfm  =  1.5e0_rt * spa(i,j,k) * bfv;
                        bfm2 = -0.5e0_rt * spa(i,j,k) * bfv;
                    }
                    else {
                        bfm = spa(i,j,k) * bfv;
                    }

                }
#ifndef AMREX_USE_GPU
                else {

                    amrex::Error("hmmat3: unsupported boundary type");

                }
#endif

                mat(i,j,k)[0] += bfm - fac * b(i,j+1,k);
                mat(i,j,k)[4] = 0.e0_rt;
                if (bho >= 1) {
                    mat(i,j,k)[3] += bfm2;
                }

            }

        }
        else if (mask.contains(i,j,k-1)) {

            if (zlo && mask(i,j,k-1) > 0) {

                if (bctype == -1) {
                    bct = tf(i,j,k-1);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {

                    if (bho >= 1) {
                        h2 = 0.5e0_rt * h;
                        th2 = 3.e0_rt * h2;
                        bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k);
                        bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k);
                    }
                    else {
                        bfv = (beta / h) / (0.5e0_rt * h + bcl);
                        bfm = bfv * b(i,j,k);
                    }

                }
                else if (bct == AMREX_LO_NEUMANN) {

                    bfm  = 0.e0_rt;
                    bfm2 = 0.e0_rt;

                }
                else if (bct == AMREX_LO_MARSHAK) {

                    bfv = 2.e0_rt * c * beta * r / h;

                    if (bho >= 1) {
                        bfm  =  0.375e0_rt * bfv;
                        bfm2 = -0.125e0_rt * bfv;
                    }
                    else {
                        bfm = 0.25e0_rt * bfv;
                    }

                }
                else if (bct == AMREX_LO_SANCHEZ_POMRANING) {

                    bfv = 2.e0_rt * c * beta * r / h;

                    if (bho >= 1) {
                        bfm  =  1.5e0_rt * spa(i,j,k) * bfv;
                        bfm2 = -0.5e0_rt * spa(i,j,k) * bfv;
                    }
                    else {
                        bfm = spa(i,j,k) * bfv;
                    }

                }
#ifndef AMREX_USE_GPU
                else {

                    amrex::Error("hmmat3: unsupported boundary type");

                }
#endif

                mat(i,j,k)[0] += bfm - fac * b(i,j,k);
                mat(i,j,k)[5] = 0.e0_rt;
                if (bho >= 1) {
                    mat(i,j,k)[6] += bfm2;
                }

            }

        }
        else if (mask.contains(i,j,k+1)) {

            if (zhi && mask(i,j,k+1) > 0) {

                if (bctype == -1) {
                    bct = tf(i,j,k+1);
                }
                else {
                    bct = bctype;
                }

                if (bct == AMREX_LO_DIRICHLET) {

                    if (bho >= 1) {
                        h2 = 0.5e0_rt * h;
                        th2 = 3.e0_rt * h2;
                        bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k+1);
                        bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k+1);
                    }
                    else {
                        bfv = (beta / h) / (0.5e0_rt * h + bcl);
                        bfm = bfv * b(i,j,k+1);
                    }

                }
                else if (bct == AMREX_LO_NEUMANN) {

                    bfm  = 0.e0_rt;
                    bfm2 = 0.e0_rt;

                }
                else if (bct == AMREX_LO_MARSHAK) {

                    bfv = 2.e0_rt * c * beta * r / h;

                    if (bho >= 1) {
                        bfm  =  0.375e0_rt * bfv;
                        bfm2 = -0.125e0_rt * bfv;
                    }
                    else {
                        bfm = 0.25e0_rt * bfv;
                    }

                }
                else if (bct == AMREX_LO_SANCHEZ_POMRANING) {

                    bfv = 2.e0_rt * c * beta * r / h;

                    if (bho >= 1) {
                        bfm  =  1.5e0_rt * spa(i,j,k) * bfv;
                        bfm2 = -0.5e0_rt * spa(i,j,k) * bfv;
                    }
                    else {
                        bfm = spa(i,j,k) * bfv;
                    }

                }
#ifndef AMREX_USE_GPU
                else {

                    amrex::Error("hmmat3: unsupported boundary type");

                }
#endif

                mat(i,j,k)[0] += bfm - fac * b(i,j,k+1);
                mat(i,j,k)[6] = 0.e0_rt;
                if (bho >= 1) {
                    mat(i,j,k)[5] += bfm2;
                }

            }

        }
    });

    Gpu::synchronize();
}

void HypreMultiABec::loadMatrix()
{
  BL_PROFILE("HypreMultiABec::loadMatrix");

#if (AMREX_SPACEDIM == 1)
  // if we were really 1D:
/*
  int offsets[3][1] = {{ 0}
                       {-1},
                       { 1}};
*/
  // fake 1D as a 2D problem:
  int offsets[3][2] = {{ 0,  0},
                       {-1,  0},
                       { 1,  0}};
#elif (AMREX_SPACEDIM == 2)
  int offsets[5][2] = {{ 0,  0},
                       {-1,  0},
                       { 1,  0},
                       { 0, -1},
                       { 0,  1}};
#elif (AMREX_SPACEDIM == 3)
  int offsets[7][3] = {{ 0,  0,  0},
                       {-1,  0,  0},
                       { 1,  0,  0},
                       { 0, -1,  0},
                       { 0,  1,  0},
                       { 0,  0, -1},
                       { 0,  0,  1}};
#endif

  const int size = 2 * AMREX_SPACEDIM + 1;
  int i, idim;

  int stencil_indices[size];

  for (i = 0; i < size; i++) {
    stencil_indices[i] = i;
  }

  BaseFab<GpuArray<Real, size>> matfab; // AoS indexing
  FArrayBox smatfab;
  for (int level = crse_level; level <= fine_level; level++) {
    int part = level - crse_level;

    for (MFIter mfi(*acoefs[level]); mfi.isValid(); ++mfi) {
      i = mfi.index();
      const Box &reg = grids[level][i];

      matfab.resize(reg);
      Real* mat = (Real*) matfab.dataPtr();
      Elixir mat_elix = matfab.elixir();

      // build matrix interior

      hmac(reg, matfab.array(), (*acoefs[level])[mfi].array(), alpha);

      for (idim = 0; idim < AMREX_SPACEDIM; idim++) {
          hmbc(reg, matfab.array(), (*bcoefs[level])[idim][mfi].array(),
               beta, geom[level].CellSize(), idim);
      }

      // add b.c.'s to matrix diagonal, and
      // zero out offdiag values at domain boundaries

      const Box& domain = bd[level]->getDomain();
      for (OrientationIter oitr; oitr; oitr++) {
        int cdir(oitr());
        idim = oitr().coordDir();
        const RadBoundCond &bct = bd[level]->bndryConds(oitr())[i];
        const Real      &bcl = bd[level]->bndryLocs(oitr())[i];
        const Mask      &msk = bd[level]->bndryMasks(oitr(),i);
        const Box &bbox = (*bcoefs[level])[idim][mfi].box();
        const Box &msb  = msk.box();
        if (reg[oitr()] == domain[oitr()] || level == crse_level) {

          // Treat an exposed grid edge here as a boundary condition
          // for the linear solver:

          if (reg[oitr()] == domain[oitr()]) {
            Array4<const int> tfp{};
            int bctype = bct;
            if (bd[level]->mixedBndry(oitr())) {
              const BaseFab<int> &tf = *(bd[level]->bndryTypes(oitr())[i]);
              tfp = tf.array();
              bctype = -1;
            }
            const FArrayBox &fs = bd[level]->bndryValues(oitr())[mfi];
            Array4<const Real> pSPa{};
            if (SPa[level]) {
                pSPa = (*SPa[level])[mfi].array();
            }
            hmmat3(reg, oitr().isLow(), idim, matfab.array(), cdir, bctype,
                   tfp, bho, bcl, msk.array(), (*bcoefs[level])[idim][mfi].array(),
                   beta, geom[level].CellSize(), geom[level].data(),
                   flux_factor, pSPa);
          }
          else {
              hmmat(reg, matfab.array(), cdir, bct, bho, bcl,
                    msk.array(), (*bcoefs[level])[idim][mfi].array(),
                    beta, geom[level].CellSize());
          }
        }
        else {

          // An exposed grid edge here actually borders the next coarser
          // level in the current linear system.  Zero out the interior
          // stencil using Neumann BC:

          const RadBoundCond bct_coarse = AMREX_LO_NEUMANN;
          hmmat(reg, matfab.array(), cdir, bct_coarse, bho, bcl,
                msk.array(), (*bcoefs[level])[idim][mfi].array(),
                beta, geom[level].CellSize());
        }
      }

      // initialize matrix

      if (subgrids[level][i].size() > 0) {
        for (int j = 0; j < subgrids[level][i].size(); j++) {
          const Box& sreg = subgrids[level][i][j];
          smatfab.resize(sreg,size);
          Real* smat = smatfab.dataPtr();
          for (IntVect v = sreg.smallEnd(); v <= sreg.bigEnd(); sreg.next(v)) {
            int is = sreg.index(v);
            int ir =  reg.index(v);
            for (int s = 0; s < size; s++) {
              smat[is * size + s] = mat[ir * size + s];
            }
          }
          HYPRE_SStructMatrixSetBoxValues(A, part, loV(sreg), hiV(sreg), 0,
                                          size, stencil_indices, smat);
        }
      }
      else {
        HYPRE_SStructMatrixSetBoxValues(A, part, loV(reg), hiV(reg), 0,
                                        size, stencil_indices, mat);
      }
    }

    // Add coarse-fine interface entries to the matrix here:
    if (level == crse_level) {
      continue;
    }

    // First we do the entries as seen by the fine cells adjacent
    // to the interface, working on the fine processor:

    BndryAuxVar entry(grids[level], dmap[level], BndryAuxVar::INTERIOR);
    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      IntVect vin = amrex::BASISV(idir), ves;
      vin = (ori.isLow() ? -vin : vin);    // outward normal unit vector
      ves = (ori.isLow() ?  ves : vin);    // edge shift vector
      Real h = geom[level].CellSize(idir); // normal fine grid spacing
      Real ffac = (-beta / h);             // divergence factor
      for (int i = cintrp[level]->firstLocal(); cintrp[level]->isValid(i);
           i = cintrp[level]->nextLocal(i)) {
        Box reg = amrex::adjCell(grids[level][i], ori);
        reg.shift(-vin); // fine interior cells
        const Mask &msk = bd[level]->bndryMasks(ori,i);
        for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
          if (msk(v+vin) == RadBndryData::not_covered) {
              entry(ori,i)(v).push(&(*ederiv[level])(ori,i)(v+vin),
                                    ffac * (*bcoefs[level])[idir][i](v+ves));
          }
        }
      }
    }

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      IntVect vin = amrex::BASISV(idir);
      vin = (ori.isLow() ? -vin : vin);    // outward normal unit vector
      for (int i = cintrp[level]->firstLocal(); cintrp[level]->isValid(i);
           i = cintrp[level]->nextLocal(i)) {
        Box reg = amrex::adjCell(grids[level][i], ori);
        reg.shift(-vin); // fine interior cells
//        const Mask &msk = bd[level]->bndryMasks(ori,i);
        for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
          if (!entry(ori,i)(v).empty() &&
              !entry(ori,i)(v).secondary()) {
            entry(ori,i)(v).collapse();
            Vector<int> levels;
            Vector<IntVect> cells;
            int retval = entry(ori,i)(v).get_locations(levels, cells);
            BL_ASSERT(retval == 0);
            Vector<Real> values;
            retval = entry(ori,i)(v).get_coeffs(values);
            BL_ASSERT(retval == 0);
            int ientry = 2 * AMREX_SPACEDIM + 1;
            for (int j = 0; j < levels.size(); j++) {
              // identify stencil-like connections for separate treatment:
              int not_stencil = 1;
              if (levels[j] == level) {
                IntVect d = cells[j] - v;
                for (int k = 0; k < 2 * AMREX_SPACEDIM + 1; k++) {
                  if (d == IntVect(offsets[k])) {
                    not_stencil = 0;
                    HYPRE_SStructMatrixAddToValues(A, part, getV1(v), 0,
                                                   1, &k, &values[j]);
                  }
                }
              }
              if (not_stencil) {
                HYPRE_SStructMatrixSetValues(A, part, getV1(v), 0,
                                             1, &ientry, &values[j]);
                ientry++;
              }
            }
          }
        }
      }
    }

    // Now add the matrix values seen by the coarse cells adjacent
    // to the coarse-fine interface.  These are averages of ederiv
    // over the fine faces making up each coarse face.  We use
    // CrseBndryAuxVar now because we need to work on the coarse
    // processor.

    IntVect rat = fine_ratio[level-1];
    const BoxArray& f_fgrids(grids[level]);
    const BoxArray& c_cgrids(grids[level-1]);
    BoxArray f_cgrids(c_cgrids);
    f_cgrids.refine(rat);
    BoxArray c_fgrids(f_fgrids);
    c_fgrids.coarsen(rat);
    //CrseBndryAuxVar c_entry(c_cgrids, c_fgrids, BndryAuxVar::EXTERIOR);
    c_entry[level]->reinitialize_connections(BndryAuxVar::EXTERIOR);
    c_entry[level]->buildFaceData(rat);

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      c_entry[level]->loadFaceData(ori, (*bcoefs[level])[idir], 0, 0, 1);
      IntVect vin = amrex::BASISV(idir), ves;
      vin = (ori.isLow() ? -vin : vin); // outward normal unit vector
      ves = (ori.isLow() ? -vin : ves); // edge shift vector (diff from above)
      Real hc = geom[level-1].CellSize(idir); // normal coarse grid spacing
      Real cfac = (beta / hc);                // divergence factor
      //cfac = 0.0;
      Real zfac = (-cfac / hc);               // factor for covered cell
      //cfac = 0.0;
      IntVect ve; // default constructor initializes to zero
#if (AMREX_SPACEDIM >= 2)
      int jdir = (idir + 1) % AMREX_SPACEDIM;
      ve += (rat[jdir] - 1) * amrex::BASISV(jdir);
      cfac /= rat[jdir]; // will average over fine cells in tangential dir
#endif
#if (AMREX_SPACEDIM == 3)
      int kdir = (idir + 2) % 3;
      ve += (rat[kdir] - 1) * amrex::BASISV(kdir);
      cfac /= rat[kdir]; // will average over fine cells in tangential dir
#endif
      for (int i = c_entry[level]->firstLocal(); c_entry[level]->isValid(i);
           i = c_entry[level]->nextLocal(i)) {
        // parallel loop is tied to coarse grids
          for (int j = 0; j < (*c_entry[level]).size(ori,i); j++) {
          const Box& reg = (*c_cintrp[level])(ori,i,j).box(); // adjacent cells
          const Box& creg = (*c_entry[level])(ori,i,j).box(); // adjacent cells
          const Mask& msk = c_cintrp[level]->mask(ori,i,j); // fine mask
          const FArrayBox& fbcoefs = c_entry[level]->faceData(ori,i,j);
          for (IntVect vc = creg.smallEnd(); vc <= creg.bigEnd(); creg.next(vc)) {
            IntVect vf = rat * vc;
            vf[idir] = reg.smallEnd(idir); // same as bigEnd(idir)
            Box face(vf, vf + ve);
            if (msk(vf) == RadBndryData::not_covered) {
              // Zero out connection to covered coarse cell:
              (*c_entry[level])(ori,i,j)(vc).push(-1, vc-vin, 0.0);
              (*c_entry[level])(ori,i,j)(vc).push(level-1, vc,
                                                       zfac * (*bcoefs[level-1])[idir][i](vc+ves));
              // Add fine fluxes over face of coarse cell:
              for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
                  (*c_entry[level])(ori,i,j)(vc)
                      .push(&(*c_ederiv[level])(ori,i,j)(v), cfac * fbcoefs(v+ves));
              }
            }
          }
        }
      }
    }

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      for (int i = c_entry[level]->firstLocal(); c_entry[level]->isValid(i);
           i = c_entry[level]->nextLocal(i)) {
        // parallel loop is tied to coarse grids
          for (int j = 0; j < (*c_entry[level]).size(ori,i); j++) {
          const Box& reg = (*c_cintrp[level])(ori,i,j).box(); // adjacent cells
          const Box& creg = (*c_entry[level])(ori,i,j).box(); // adjacent cells
          const Mask &msk = c_cintrp[level]->mask(ori,i,j); // fine mask
          for (IntVect vc = creg.smallEnd(); vc <= creg.bigEnd(); creg.next(vc)) {
            IntVect vf = rat * vc;
            vf[idir] = reg.smallEnd(idir); // same as bigEnd(idir)
            if (msk(vf) == RadBndryData::not_covered &&
                !(*c_entry[level])(ori,i,j)(vc).secondary()) {
              (*c_entry[level])(ori,i,j)(vc).collapse();
              Vector<int> levels;
              Vector<IntVect> cells;
              int retval = (*c_entry[level])(ori,i,j)(vc)
                .get_locations(levels, cells);
              BL_ASSERT(retval == 0);
              Vector<Real> values;
              retval = (*c_entry[level])(ori,i,j)(vc).get_coeffs(values);
              BL_ASSERT(retval == 0);
              int ientry = 2 * AMREX_SPACEDIM + 1;
              for (int jj = 0; jj < levels.size(); jj++) {
                // identify stencil-like connections for separate treatment:
                int not_stencil = 1;
                if (levels[jj] == -1) {
                  // connection to covered coarse cell to be zeroed out:
                  IntVect d = cells[jj] - vc;
                  for (int k = 0; k < 2 * AMREX_SPACEDIM + 1; k++) {
                    if (d == IntVect(offsets[k])) {
                      not_stencil = 0;
                      HYPRE_SStructMatrixSetValues(A, part-1, getV1(vc), 0,
                                                   1, &k, &values[jj]);
                      BL_ASSERT(values[jj] == 0.0);
                    }
                  }
                  BL_ASSERT(not_stencil == 0);
                }
                if (levels[jj] == level-1) {
                  // other coarse-level entry, may or may not be in stencil:
                  IntVect d = cells[jj] - vc;
                  for (int k = 0; k < 2 * AMREX_SPACEDIM + 1; k++) {
                    if (d == IntVect(offsets[k])) {
                      not_stencil = 0;
                      HYPRE_SStructMatrixAddToValues(A, part-1, getV1(vc), 0,
                                                     1, &k, &values[jj]);
                    }
                  }
                }
                if (not_stencil) {
                  HYPRE_SStructMatrixSetValues(A, part-1, getV1(vc), 0,
                                               1, &ientry, &values[jj]);
                  ientry++;
                }
              }
            }
          }
        }
      }
    }
  }
}

void HypreMultiABec::finalizeMatrix()
{
  BL_PROFILE("HypreMultiABec::finalizeMatrix");

  HYPRE_SStructMatrixAssemble(A);
}

void HypreMultiABec::loadLevelVectors(int level,
                                      MultiFab& dest,
                                      int icomp,
                                      MultiFab& rhs,
                                      BC_Mode inhom)
{
  BL_PROFILE("HypreMultiABec::loadLevelVectors");

  int part = level - crse_level;

  Real *vec;
  FArrayBox fnew;
  for (MFIter mfi(dest); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = grids[level][i];

    // initialize dest, since we will reuse the space to set up rhs below:

    FArrayBox *f;
    int fcomp;
    if (dest.nGrow() == 0) { // need a temporary if dest is the wrong size
      f = &dest[mfi];
      fcomp = icomp;
    }
    else {
      f = &fnew;
      f->resize(reg);
      f->copy<RunOn::Device>(dest[mfi], icomp, 0, 1);
      fcomp = 0;
    }
    Elixir f_elix = fnew.elixir();

    vec = f->dataPtr(fcomp); // sharing space, dest will be overwritten below

    vectorSetBoxValues(x, part, reg, subgrids[level][i], vec);

    f->copy<RunOn::Device>(rhs[mfi], 0, fcomp, 1);

    // add b.c.'s to rhs

    if (inhom) {
      const Box& domain = bd[level]->getDomain();
      for (OrientationIter oitr; oitr; oitr++) {
        int cdir(oitr());
        int idim = oitr().coordDir();
        const RadBoundCond &bct = bd[level]->bndryConds(oitr())[i];
        const Real      &bcl = bd[level]->bndryLocs(oitr())[i];
        const FArrayBox       &fs  = bd[level]->bndryValues(oitr())[mfi];
        const Mask      &msk = bd[level]->bndryMasks(oitr(), i);
        const Box &bbox = (*bcoefs[level])[idim][mfi].box();

        if (reg[oitr()] == domain[oitr()] || level == crse_level) {

          // Treat an exposed grid edge here as a boundary condition
          // for the linear solver:

          if (reg[oitr()] == domain[oitr()]) {
            Array4<const int> tfp{};
            int bctype = bct;
            if (bd[level]->mixedBndry(oitr())) {
              const BaseFab<int> &tf = *(bd[level]->bndryTypes(oitr())[i]);
              tfp = tf.array();
              bctype = -1;
            }
            HypreABec::hbvec3(reg,
                              oitr().isLow(), idim,
                              f->array(fcomp),
                              cdir, bctype,
                              tfp,
                              bho, bcl,
                              fs.array(bdcomp),
                              msk.array(),
                              (*bcoefs[level])[idim][mfi].array(),
                              beta, geom[level].data());
          }
          else {
              HypreABec::hbvec(reg, f->array(fcomp),
                               cdir, bct, bho, bcl,
                               fs.array(bdcomp), msk.array(),
                               (*bcoefs[level])[idim][mfi].array(),
                               beta, geom[level].CellSize());
          }
        }
        // There is no else here, since we would then be at an
        // interior edge and would not need to add anything to vec.
      }
    }

    // initialize rhs

    vectorSetBoxValues(b, part, reg, subgrids[level][i], vec);
  }
}

void HypreMultiABec::loadLevelVectorX(int level,
                                      MultiFab& dest,
                                      int icomp)
{
  BL_PROFILE("HypreMultiABec::loadLevelVectorX");

  int part = level - crse_level;

  FArrayBox fnew;
  for (MFIter mfi(dest); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = grids[level][i];

    FArrayBox *f;
    int fcomp;
    if (dest.nGrow() == 0) { // need a temporary if dest is the wrong size
      f = &dest[mfi];
      fcomp = icomp;
    }
    else {
      f = &fnew;
      f->resize(reg);
      f->copy<RunOn::Device>(dest[mfi], icomp, 0, 1);
      fcomp = 0;
    }
    Elixir f_elix = fnew.elixir();

    Real* vec = f->dataPtr(fcomp);

    vectorSetBoxValues(x, part, reg, subgrids[level][i], vec);
  }
}

void HypreMultiABec::loadLevelVectorB(int level,
                                      MultiFab& rhs, // will be altered
                                      BC_Mode inhom)
{
  BL_PROFILE("HypreMultiABec::loadLevelVectorB");

  int part = level - crse_level;

  FArrayBox fnew;
  for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = grids[level][i];

    FArrayBox *f;
    if (rhs.nGrow() == 0) { // need a temporary if rhs is the wrong size
      f = &rhs[mfi];
    }
    else {
      f = &fnew;
      f->resize(reg);
      f->copy<RunOn::Device>(rhs[mfi]);
    }
    Elixir f_elix = fnew.elixir();

    Real* vec = f->dataPtr();

    // add b.c.'s to rhs

    if (inhom) {
      const Box& domain = bd[level]->getDomain();
      for (OrientationIter oitr; oitr; oitr++) {
        int cdir(oitr());
        int idim = oitr().coordDir();
        const RadBoundCond &bct = bd[level]->bndryConds(oitr())[i];
        const Real      &bcl = bd[level]->bndryLocs(oitr())[i];
        const FArrayBox       &fs  = bd[level]->bndryValues(oitr())[mfi];
        const Mask      &msk = bd[level]->bndryMasks(oitr(), i);
        const Box &bbox = (*bcoefs[level])[idim][mfi].box();
        const Box &msb  = msk.box();
        if (reg[oitr()] == domain[oitr()] || level == crse_level) {

          // Treat an exposed grid edge here as a boundary condition
          // for the linear solver:

          if (reg[oitr()] == domain[oitr()]) {
            Array4<const int> tfp{};
            int bctype = bct;
            if (bd[level]->mixedBndry(oitr())) {
              const BaseFab<int> &tf = *(bd[level]->bndryTypes(oitr())[i]);
              tfp = tf.array();
              bctype = -1;
            }
            HypreABec::hbvec3(reg,
                              oitr().isLow(), idim,
                              f->array(),
                              cdir, bctype,
                              tfp,
                              bho, bcl,
                              fs.array(bdcomp),
                              msk.array(),
                              (*bcoefs[level])[idim][mfi].array(),
                              beta, geom[level].data());
          }
          else {
              HypreABec::hbvec(reg, f->array(),
                               cdir, bct, bho, bcl,
                               fs.array(bdcomp), msk.array(),
                               (*bcoefs[level])[idim][mfi].array(),
                               beta, geom[level].CellSize());
          }
        }
        // There is no else here, since we would then be at an
        // interior edge and would not need to add anything to vec.
      }
    }

    // initialize rhs

    vectorSetBoxValues(b, part, reg, subgrids[level][i], vec);
  }
}

void HypreMultiABec::finalizeVectors()
{
  BL_PROFILE("HypreMultiABec::finalizeVectors");

  HYPRE_SStructVectorAssemble(b);
  HYPRE_SStructVectorAssemble(x);
}

void HypreMultiABec::setupSolver(Real _reltol, Real _abstol, int maxiter)
{
  BL_PROFILE("HypreMultiABec::setupSolver");

  reltol = _reltol;
  abstol = _abstol; // may be used to change tolerance for solve

  BL_ASSERT(sstruct_solver == NULL);
  BL_ASSERT(solver         == NULL);
  BL_ASSERT(precond        == NULL);

  if (solver_flag == 100) {
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);

    HYPRE_BoomerAMGCreate(&solver);
    HYPRE_BoomerAMGSetMinIter(solver, 1);
    HYPRE_BoomerAMGSetMaxIter(solver, maxiter);
    //HYPRE_BoomerAMGSetMaxIter(solver, 1);
    HYPRE_BoomerAMGSetTol(solver, reltol);
#if 0
    // Barry used these four settings, at least for level solves:
    HYPRE_BoomerAMGSetInterpType(solver, 6);
    HYPRE_BoomerAMGSetCoarsenType(solver, 8);
    HYPRE_BoomerAMGSetStrongThreshold(solver, 0.25);
    HYPRE_BoomerAMGSetTruncFactor(solver, 0.);
#endif
    //HYPRE_BoomerAMGSetStrongThreshold(solver, 0.6);
    //HYPRE_BoomerAMGSetPrintLevel(solver, 2);
    //HYPRE_BoomerAMGSetTruncFactor(solver, 0.5);
    //HYPRE_BoomerAMGSetCoarsenType(solver, 6);
    //HYPRE_BoomerAMGSetLogging(solver, 2);
    HYPRE_BoomerAMGSetup(solver, par_A, par_b, par_x);
#if 0
    // The default num_sweeps seems to work best.  Note that changing
    // the default can't be done with a temp array on the stack, as the
    // array will not be copied.  The following code will leak since
    // num_sweeps is never deleted, but at least it will work:
    int num_sweeps[4] = {2, 2, 2, 1};
    HYPRE_BoomerAMGSetNumGridSweeps(solver, num_sweeps);
#endif
    //HYPRE_BoomerAMGSetIOutDat(solver, 2);
  }
  else if (solver_flag == 101) {
    int nparts = fine_level - crse_level + 1;
    int *plevels = new int[nparts];
    int (*prefinements)[3] = new int[nparts][3];
    for (int i = 0; i < nparts; i++) {
      plevels[i] = i;
      if (i == 0) {
        prefinements[i][0] = 1;
        prefinements[i][1] = 1;
        prefinements[i][2] = 1;
      }
      else {
        prefinements[i][0] = fine_ratio[crse_level+i-1][0];
#if (AMREX_SPACEDIM == 1)
        prefinements[i][1] = 1;
        prefinements[i][2] = 1;
#elif (AMREX_SPACEDIM == 2)
        prefinements[i][1] = fine_ratio[crse_level+i-1][1];
        prefinements[i][2] = 1;
#elif (AMREX_SPACEDIM == 3)
        prefinements[i][1] = fine_ratio[crse_level+i-1][1];
        prefinements[i][2] = fine_ratio[crse_level+i-1][2];
#endif
      }
    }
    int n_pre  = prefinements[nparts-1][0] - 1;
    int n_post = prefinements[nparts-1][0] - 1;

    n_pre  = 3;
    n_post = 3;

    HYPRE_SStructFACCreate(MPI_COMM_WORLD, &sstruct_solver);
    HYPRE_SStructFACSetMaxLevels(sstruct_solver, nparts);
    HYPRE_SStructFACSetMaxIter(sstruct_solver, maxiter);
    HYPRE_SStructFACSetTol(sstruct_solver, reltol);
    HYPRE_SStructFACSetPLevels(sstruct_solver, nparts, plevels);
    HYPRE_SStructFACSetPRefinements(sstruct_solver, nparts, prefinements);
    HYPRE_SStructFACSetRelChange(sstruct_solver, 0);
    HYPRE_SStructFACSetRelaxType(sstruct_solver, 1);
    //HYPRE_SStructFACSetRelaxType(sstruct_solver, 2);
    HYPRE_SStructFACSetRelaxType(sstruct_solver, 2);
    HYPRE_SStructFACSetNumPreRelax(sstruct_solver, n_pre);
    HYPRE_SStructFACSetNumPostRelax(sstruct_solver, n_post);
    HYPRE_SStructFACSetCoarseSolverType(sstruct_solver, 2);
    HYPRE_SStructFACSetLogging(sstruct_solver, 1);

    for (int i = 1; i < nparts; i++) {
      HYPRE_SStructFACZeroCFSten(A, hgrid, i, prefinements[i]);
      HYPRE_SStructFACZeroFCSten(A, hgrid, i);
    }
    for (int i = 0; i < nparts-1; i++) {
      HYPRE_SStructFACZeroAMRMatrixData(A, i, prefinements[i+1]);
    }
    HYPRE_SStructFACZeroAMRVectorData(b, plevels, prefinements);
    finalizeMatrix();
    finalizeVectors();

    HYPRE_SStructFACSetup2(sstruct_solver, A, b, x);
    std::cin.get();

    delete[] plevels;
    delete[] prefinements;
  }
  else if (solver_flag == 102) {
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);

    ParmParse pp("hmabec");
    int kdim = 5; pp.query("kdim", kdim);
    static int first = 1;
    if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
      first = 0;
      std::cout << "hmabec.kdim            = " << kdim << std::endl;
    }

    HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_ParCSRGMRESSetKDim(solver, kdim);
    HYPRE_ParCSRGMRESSetMaxIter(solver, maxiter);
    HYPRE_ParCSRGMRESSetTol(solver, reltol);
    //HYPRE_ParCSRGMRESSetLogging(solver, 1);

    HYPRE_ParCSRGMRESSetup(solver, par_A, par_b, par_x);
  }
  else if (solver_flag == 103) {
    ParmParse pp("hmabec");
    int kdim = 5; pp.query("kdim", kdim);
    static int first = 1;
    if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
      first = 0;
      std::cout << "hmabec.kdim            = " << kdim << std::endl;
    }

    HYPRE_SStructGMRESCreate(MPI_COMM_WORLD, &sstruct_solver);
    HYPRE_SStructGMRESSetKDim(sstruct_solver, kdim);
    HYPRE_SStructGMRESSetMaxIter(sstruct_solver, maxiter);
    HYPRE_SStructGMRESSetTol(sstruct_solver, reltol);
    //HYPRE_SStructGMRESSetLogging(sstruct_solver, 1);

    HYPRE_SStructGMRESSetup(sstruct_solver, A, b, x);
  }
  else if (solver_flag == 1002) {
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);

    HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_ParCSRPCGSetMaxIter(solver, maxiter);
    HYPRE_ParCSRPCGSetTol(solver, reltol);
    //HYPRE_ParCSRPCGSetLogging(solver, 1);

    HYPRE_ParCSRPCGSetup(solver, par_A, par_b, par_x);
  }
  else if (solver_flag == 1003) {
    HYPRE_SStructPCGCreate(MPI_COMM_WORLD, &sstruct_solver);
    HYPRE_SStructPCGSetMaxIter(sstruct_solver, maxiter);
    HYPRE_SStructPCGSetTol(sstruct_solver, reltol);
    //HYPRE_SStructPCGSetLogging(sstruct_solver, 1);

    HYPRE_SStructPCGSetup(sstruct_solver, A, b, x);
  }
  else if (solver_flag == 104) {
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);

    ParmParse pp("hmabec");
    int kdim = 5; pp.query("kdim", kdim);
    static int first = 1;
    if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
      first = 0;
      std::cout << "hmabec.kdim            = " << kdim << std::endl;
    }

    HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_ParCSRGMRESSetKDim(solver, kdim);
    HYPRE_ParCSRGMRESSetMaxIter(solver, maxiter);
    HYPRE_ParCSRGMRESSetTol(solver, reltol);
    //HYPRE_ParCSRGMRESSetLogging(solver, 1);

    HYPRE_BoomerAMGCreate(&precond);
    HYPRE_BoomerAMGSetMaxIter(precond, 1);
    HYPRE_BoomerAMGSetTol(precond, reltol);
    //HYPRE_BoomerAMGSetStrongThreshold(precond, 0.6);
    //HYPRE_BoomerAMGSetPrintLevel(precond, 2);
    //HYPRE_BoomerAMGSetTruncFactor(precond, 0.5);
    //HYPRE_BoomerAMGSetCoarsenType(precond, 6);
    //HYPRE_BoomerAMGSetLogging(precond, 2);
    HYPRE_ParCSRGMRESSetPrecond(solver,
                                (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSolve,
                                (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSetup,
                                precond);

    HYPRE_ParCSRGMRESSetup(solver, par_A, par_b, par_x);
  }
  else if (solver_flag == 105) {
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);

    ParmParse pp("hmabec");
    int kdim    = 5;    pp.query("kdim", kdim);
    int interp  = 6;    pp.query("interp", interp);
    int coarsen = 8;    pp.query("coarsen", coarsen);
    int relax   = 6;    pp.query("relax", relax);
    Real strong = 0.25; pp.query("strong", strong);
    Real trunc  = 0.0;  pp.query("trunc", trunc);
    int p_max_elmts    = -1; pp.query("p_max_elmts", p_max_elmts);
    int agg_num_levels = -1; pp.query("agg_num_levels", agg_num_levels);
    int num_paths      = -1; pp.query("num_paths", num_paths);
    int print_level    = 0;  pp.query("print_level", print_level);
    static int first = 1;
    if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
      first = 0;
      std::cout << "hmabec.kdim            = " << kdim << std::endl;
      std::cout << "hmabec.interp          = " << interp << std::endl;
      std::cout << "hmabec.coarsen         = " << coarsen << std::endl;
      std::cout << "hmabec.relax           = " << relax << std::endl;
      std::cout << "hmabec.strong          = " << strong << std::endl;
      std::cout << "hmabec.trunc           = " << trunc << std::endl;
      if (p_max_elmts >= 0) {
        std::cout << "hmabec.p_max_elmts     = " << p_max_elmts << std::endl;
      }
      if (agg_num_levels >= 0) {
        std::cout << "hmabec.agg_num_levels  = " << agg_num_levels << std::endl;
      }
      if (num_paths >= 0) {
        std::cout << "hmabec.num_paths       = " << num_paths << std::endl;
      }
    }

    HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_ParCSRGMRESSetKDim(solver, kdim);
    HYPRE_ParCSRGMRESSetMaxIter(solver, maxiter);
    HYPRE_ParCSRGMRESSetTol(solver, reltol);
    //HYPRE_ParCSRGMRESSetLogging(solver, 1);

    HYPRE_BoomerAMGCreate(&precond);
    HYPRE_BoomerAMGSetMaxIter(precond, 1);
    HYPRE_BoomerAMGSetTol(precond, reltol);

    HYPRE_BoomerAMGSetInterpType(precond, interp);
    HYPRE_BoomerAMGSetCoarsenType(precond, coarsen);
    HYPRE_BoomerAMGSetRelaxType(precond, relax);
    HYPRE_BoomerAMGSetStrongThreshold(precond, strong);
    HYPRE_BoomerAMGSetTruncFactor(precond, trunc);
    if (p_max_elmts >= 0) {
      HYPRE_BoomerAMGSetPMaxElmts(precond, p_max_elmts);
    }
    if (agg_num_levels >= 0) {
      HYPRE_BoomerAMGSetAggNumLevels(precond, agg_num_levels);
    }
    if (num_paths >= 0) {
      HYPRE_BoomerAMGSetNumPaths(precond, num_paths);
    }
    HYPRE_BoomerAMGSetPrintLevel(precond, print_level);
    //HYPRE_BoomerAMGSetLogging(precond, 2);
    HYPRE_ParCSRGMRESSetPrecond(solver,
                                (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSolve,
                                (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSetup,
                                precond);

    HYPRE_ParCSRGMRESSetup(solver, par_A, par_b, par_x);
  }
  else if (solver_flag == 106) {
    // split solver
    ParmParse pp("hmabec");
    int struct_iter = 1; pp.query("struct_iter", struct_iter);
#if (AMREX_SPACEDIM == 1)
    int struct_flag = 0;
#else
    int struct_flag = 1;
#endif
    pp.query("struct_flag", struct_flag);
    static int first = 1;
    if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
      first = 0;
      std::cout << "hmabec.struct_iter     = " << struct_iter << std::endl;
      if (struct_flag == 0) {
        std::cout << "hmabec.struct_flag     = HYPRE_SMG" << std::endl;
      }
      else if (struct_flag == 1) {
        std::cout << "hmabec.struct_flag     = HYPRE_PFMG" << std::endl;
      }
      else {
        std::cout << "hmabec.struct_flag     = HYPRE_Jacobi" << std::endl;
      }
    }

    HYPRE_SStructSplitCreate(MPI_COMM_WORLD, &sstruct_solver);
    HYPRE_SStructSplitSetMaxIter(sstruct_solver, maxiter);
#if 0
    HYPRE_SStructSplitSetStructSolverNumIterations(sstruct_solver, struct_iter);
#endif
    HYPRE_SStructSplitSetTol(sstruct_solver, reltol);
    //HYPRE_SStructSplitSetZeroGuess(sstruct_solver);
    //HYPRE_SStructSplitSetLogging(sstruct_solver, 1);
    if (struct_flag == 0) {
      HYPRE_SStructSplitSetStructSolver(sstruct_solver, HYPRE_SMG);
    }
    else if (struct_flag == 1) {
      HYPRE_SStructSplitSetStructSolver(sstruct_solver, HYPRE_PFMG);
    }
    else {
#if 0
      HYPRE_SStructSplitSetStructSolver(sstruct_solver, HYPRE_Jacobi);
#endif
    }

    HYPRE_SStructSplitSetup(sstruct_solver, A, b, x);
  }
  else if (solver_flag == 107) {
    // split solver as preconditioner to GMRES
    ParmParse pp("hmabec");
    int kdim = 5; pp.query("kdim", kdim);
    int split_iter  = 1; pp.query("split_iter",  split_iter);
    int struct_iter = 1; pp.query("struct_iter", struct_iter);
#if (AMREX_SPACEDIM == 1)
    int struct_flag = 0;
#else
    int struct_flag = 1;
#endif
    pp.query("struct_flag", struct_flag);
    static int first = 1;
    if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
      first = 0;
      std::cout << "hmabec.kdim            = " << kdim << std::endl;
      std::cout << "hmabec.split_iter      = " << split_iter << std::endl;
      std::cout << "hmabec.struct_iter     = " << struct_iter << std::endl;
      if (struct_flag == 0) {
        std::cout << "hmabec.struct_flag     = HYPRE_SMG" << std::endl;
      }
      else if (struct_flag == 1) {
        std::cout << "hmabec.struct_flag     = HYPRE_PFMG" << std::endl;
      }
      else {
        std::cout << "hmabec.struct_flag     = HYPRE_Jacobi" << std::endl;
      }
    }

    HYPRE_SStructGMRESCreate(MPI_COMM_WORLD, &sstruct_solver);
    HYPRE_SStructGMRESSetKDim(sstruct_solver, kdim);
    HYPRE_SStructGMRESSetMaxIter(sstruct_solver, maxiter);
    HYPRE_SStructGMRESSetTol(sstruct_solver, reltol);

    HYPRE_SStructSplitCreate(MPI_COMM_WORLD, &sstruct_precond);
    HYPRE_SStructSplitSetMaxIter(sstruct_precond, split_iter);
#if 0
    HYPRE_SStructSplitSetStructSolverNumIterations(sstruct_precond, struct_iter);
#endif
    HYPRE_SStructSplitSetTol(sstruct_precond, 0.0);
    //HYPRE_SStructSplitSetMaxIter(sstruct_precond, 1);
    //HYPRE_SStructSplitSetTol(sstruct_precond, reltol);
    //HYPRE_SStructSplitSetLogging(sstruct_precond, 1);
    //HYPRE_SStructSplitSetZeroGuess(sstruct_precond);
    if (struct_flag == 0) {
      HYPRE_SStructSplitSetStructSolver(sstruct_precond, HYPRE_SMG);
    }
    else if (struct_flag == 1) {
      HYPRE_SStructSplitSetStructSolver(sstruct_precond, HYPRE_PFMG);
    }
    else {
#if 0
      HYPRE_SStructSplitSetStructSolver(sstruct_precond, HYPRE_Jacobi);
#endif
    }

    HYPRE_SStructGMRESSetPrecond(sstruct_solver,
                        (HYPRE_PtrToSStructSolverFcn) HYPRE_SStructSplitSolve,
                        (HYPRE_PtrToSStructSolverFcn) HYPRE_SStructSplitSetup,
                        sstruct_precond);
    HYPRE_SStructGMRESSetup(sstruct_solver, A, b, x);
  }
  else if (solver_flag == 108) {
    // SMG or PFMG (one level only)
    HYPRE_StructMatrix s_A;
    HYPRE_StructVector s_b;
    HYPRE_StructVector s_x;
    HYPRE_SStructMatrixGetObject(A, (void**) &s_A);
    HYPRE_SStructVectorGetObject(b, (void**) &s_b);
    HYPRE_SStructVectorGetObject(x, (void**) &s_x);
    HYPRE_StructSolver& struct_solver = *(HYPRE_StructSolver*)&solver;

    ParmParse pp("hmabec");
#if (AMREX_SPACEDIM == 1)
    int struct_flag = 0;
#else
    int struct_flag = 1;
#endif
    pp.query("struct_flag", struct_flag);
    static int first = 1;
    if (struct_flag == 0) {
      if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
        first = 0;
        std::cout << "hmabec.struct_flag     = HYPRE_SMG" << std::endl;
      }
      HYPRE_StructSMGCreate(MPI_COMM_WORLD, &struct_solver);
      HYPRE_StructSMGSetMemoryUse(struct_solver, 0);
      HYPRE_StructSMGSetMaxIter(struct_solver, maxiter);
      HYPRE_StructSMGSetRelChange(struct_solver, 0);
      HYPRE_StructSMGSetTol(struct_solver, reltol);
      HYPRE_StructSMGSetNumPreRelax(struct_solver, 1);
      HYPRE_StructSMGSetNumPostRelax(struct_solver, 1);
      HYPRE_StructSMGSetLogging(struct_solver, 1);
      HYPRE_StructSMGSetup(struct_solver, s_A, s_b, s_x);
    }
    else if (struct_flag == 1) {
      int pfmg_relax_type = 1; pp.query("pfmg_relax_type", pfmg_relax_type);
      if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
        first = 0;
        std::cout << "hmabec.struct_flag     = HYPRE_PFMG" << std::endl;
        std::cout << "hmabec.pfmg_relax_type = " << pfmg_relax_type << std::endl;
      }
      HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &struct_solver);
      //HYPRE_StructPFMGSetMemoryUse(struct_solver, 0);
      HYPRE_StructPFMGSetSkipRelax(struct_solver, 0);
      HYPRE_StructPFMGSetMaxIter(struct_solver, maxiter);
      HYPRE_StructPFMGSetRelChange(struct_solver, 0);
      HYPRE_StructPFMGSetTol(struct_solver, reltol);
      HYPRE_StructPFMGSetRelaxType(struct_solver, pfmg_relax_type);
      HYPRE_StructPFMGSetNumPreRelax(struct_solver, 1);
      HYPRE_StructPFMGSetNumPostRelax(struct_solver, 1);
      HYPRE_StructPFMGSetLogging(struct_solver, 1);
      HYPRE_StructPFMGSetup(struct_solver, s_A, s_b, s_x);
    }
  }
  else if (solver_flag == 109) {
    // SMG or PFMG as preconditioner to GMRES (one level only)
    HYPRE_StructMatrix s_A;
    HYPRE_StructVector s_b;
    HYPRE_StructVector s_x;
    HYPRE_SStructMatrixGetObject(A, (void**) &s_A);
    HYPRE_SStructVectorGetObject(b, (void**) &s_b);
    HYPRE_SStructVectorGetObject(x, (void**) &s_x);
    HYPRE_StructSolver& struct_precond = *(HYPRE_StructSolver*)&precond;
    HYPRE_StructSolver& struct_solver  = *(HYPRE_StructSolver*)&solver;

    ParmParse pp("hmabec");
    int kdim = 5; pp.query("kdim", kdim);
#if (AMREX_SPACEDIM == 1)
    int struct_flag = 0;
#else
    int struct_flag = 1;
#endif
    pp.query("struct_flag", struct_flag);

    HYPRE_StructGMRESCreate(MPI_COMM_WORLD, &struct_solver);
    //HYPRE_StructGMRESSetKDim(struct_solver, kdim);
    HYPRE_GMRESSetKDim((HYPRE_Solver) struct_solver, kdim);
    HYPRE_StructGMRESSetMaxIter(struct_solver, maxiter);
    HYPRE_StructGMRESSetTol(struct_solver, reltol);

    static int first = 1;
    if (struct_flag == 0) {
      if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
        first = 0;
        std::cout << "hmabec.kdim            = " << kdim << std::endl;
        std::cout << "hmabec.struct_flag     = HYPRE_SMG" << std::endl;
      }
      HYPRE_StructSMGCreate(MPI_COMM_WORLD, &struct_precond);
      HYPRE_StructSMGSetMemoryUse(struct_precond, 0);
      HYPRE_StructSMGSetMaxIter(struct_precond, 1);
      HYPRE_StructSMGSetRelChange(struct_precond, 0);
      HYPRE_StructSMGSetTol(struct_precond, reltol);
      HYPRE_StructSMGSetNumPreRelax(struct_precond, 1);
      HYPRE_StructSMGSetNumPostRelax(struct_precond, 1);
      HYPRE_StructSMGSetLogging(struct_precond, 1);
      HYPRE_StructGMRESSetPrecond(struct_solver,
                       (HYPRE_PtrToStructSolverFcn) HYPRE_StructSMGSolve,
                       (HYPRE_PtrToStructSolverFcn) HYPRE_StructSMGSetup,
                       struct_precond);
    }
    else if (struct_flag == 1) {
      int pfmg_relax_type = 1; pp.query("pfmg_relax_type", pfmg_relax_type);
      int pfmg_pre_relax  = 1; pp.query("pfmg_pre_relax",  pfmg_pre_relax);
      int pfmg_post_relax = 1; pp.query("pfmg_post_relax", pfmg_post_relax);
      if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
        first = 0;
        std::cout << "hmabec.kdim            = " << kdim << std::endl;
        std::cout << "hmabec.struct_flag     = HYPRE_PFMG" << std::endl;
        std::cout << "hmabec.pfmg_relax_type = " << pfmg_relax_type << std::endl;
        std::cout << "hmabec.pfmg_pre_relax  = " << pfmg_pre_relax << std::endl;
        std::cout << "hmabec.pfmg_post_relax = " << pfmg_post_relax << std::endl;
      }
      HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &struct_precond);
      //HYPRE_StructPFMGSetMemoryUse(struct_precond, 0);
      HYPRE_StructPFMGSetSkipRelax(struct_precond, 0);
      HYPRE_StructPFMGSetMaxIter(struct_precond, 1);
      HYPRE_StructPFMGSetRelChange(struct_precond, 0);
      HYPRE_StructPFMGSetTol(struct_precond, reltol);
      HYPRE_StructPFMGSetRelaxType(struct_precond, pfmg_relax_type);
      HYPRE_StructPFMGSetNumPreRelax(struct_precond, pfmg_pre_relax);
      HYPRE_StructPFMGSetNumPostRelax(struct_precond, pfmg_post_relax);
      HYPRE_StructPFMGSetLogging(struct_precond, 1);
      HYPRE_StructGMRESSetPrecond(struct_solver,
                       (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,
                       (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup,
                       struct_precond);
    }

    HYPRE_StructGMRESSetup(struct_solver, s_A, s_b, s_x);
  }

  // The next one is AMG again with Barry's different settings.
  // The three after that added by Barry, appear to be for symmetric matrices.

  else if (solver_flag == 150) {
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);
#if 0
    HYPRE_SStructMatrixPrint("A", A, 0);
    HYPRE_SStructVectorPrint("B", b, 0);
    cin.get();
    std::cout << "HypreMultiABec: creating solver" << std::endl;
#endif
    HYPRE_BoomerAMGCreate(&solver);
    HYPRE_BoomerAMGSetMinIter(solver, 1);
    HYPRE_BoomerAMGSetMaxIter(solver, maxiter);
    //HYPRE_BoomerAMGSetMaxIter(solver, 1);
    HYPRE_BoomerAMGSetTol(solver, reltol);
#if 1
    // Barry used these four settings, at least for level solves:
    HYPRE_BoomerAMGSetInterpType(solver, 6);
    HYPRE_BoomerAMGSetCoarsenType(solver, 8);
    HYPRE_BoomerAMGSetStrongThreshold(solver, 0.25);
    HYPRE_BoomerAMGSetTruncFactor(solver, 0.);
#endif
    //HYPRE_BoomerAMGSetPrintLevel(solver, 2);
    //HYPRE_BoomerAMGSetTruncFactor(solver, 0.5);
    //HYPRE_BoomerAMGSetCoarsenType(solver, 6);
    //HYPRE_BoomerAMGSetLogging(solver, 2);
    HYPRE_BoomerAMGSetup(solver, par_A, par_b, par_x);
#if 0
    // The default num_sweeps seems to work best.  Note that changing
    // the default can't be done with a temp array on the stack, as the
    // array will not be copied.  The following code will leak since
    // num_sweeps is never deleted, but at least it will work:
    int num_sweeps[4] = {2,2,2,1};
    HYPRE_BoomerAMGSetNumGridSweeps(solver, num_sweeps);
#endif
    //HYPRE_BoomerAMGSetIOutDat(solver, 2);
  }
  else if (solver_flag == 151 || solver_flag == 153) {
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);

    HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_PCGSetMaxIter(solver, maxiter);
    HYPRE_PCGSetTol(solver, reltol);
    HYPRE_PCGSetTwoNorm(solver, 1 );
    HYPRE_PCGSetRelChange(solver, 0 );

    HYPRE_BoomerAMGCreate(&precond);
    HYPRE_BoomerAMGSetTol(precond, 0.0);
    HYPRE_BoomerAMGSetMaxIter(precond, 1);
    HYPRE_BoomerAMGSetCoarsenType(precond, 8);
    HYPRE_BoomerAMGSetInterpType(precond, 6);
    HYPRE_BoomerAMGSetStrongThreshold(precond, 0.25);
    HYPRE_BoomerAMGSetTruncFactor(precond, 0.);
    if (solver_flag == 151) {
       HYPRE_BoomerAMGSetAggNumLevels(precond, 1);
       HYPRE_BoomerAMGSetNumPaths(precond, 2);
    }
    HYPRE_PCGSetPrecond(solver,
                        (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                        (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
                        precond );

    HYPRE_PCGSetup(solver, (HYPRE_Matrix)par_A,
                   (HYPRE_Vector)par_b, (HYPRE_Vector)par_x);
  }

  else if (solver_flag == 152) {
    HYPRE_SStructPCGCreate(MPI_COMM_WORLD, (HYPRE_SStructSolver *) &solver);
    HYPRE_PCGSetMaxIter(solver, maxiter);
    HYPRE_PCGSetTol(solver, reltol);
    HYPRE_PCGSetTwoNorm(solver, 1 );
    HYPRE_PCGSetRelChange(solver, 0 );

    HYPRE_SStructSplitCreate(MPI_COMM_WORLD, (HYPRE_SStructSolver *) &precond);
    HYPRE_SStructSplitSetMaxIter((HYPRE_SStructSolver) precond, 1);
    HYPRE_SStructSplitSetTol((HYPRE_SStructSolver) precond, 0.0);
    HYPRE_SStructSplitSetZeroGuess((HYPRE_SStructSolver) precond);
    HYPRE_SStructSplitSetStructSolver((HYPRE_SStructSolver) precond, HYPRE_PFMG);
    HYPRE_PCGSetPrecond(solver,
                        (HYPRE_PtrToSolverFcn) HYPRE_SStructSplitSolve,
                        (HYPRE_PtrToSolverFcn) HYPRE_SStructSplitSetup,
                        precond);

    HYPRE_PCGSetup(solver, (HYPRE_Matrix) A,
                   (HYPRE_Vector) b, (HYPRE_Vector) x);

  }
  else {
    std::cout << "HypreMultiABec: no such solver" << std::endl;
    exit(1);
  }
}

void HypreMultiABec::clearSolver()
{
  BL_PROFILE("HypreMultiABec::clearSolver");

  if (solver_flag == 100) {
    HYPRE_BoomerAMGDestroy(solver);
  }
  else if (solver_flag == 101) {
    HYPRE_SStructFACDestroy2(sstruct_solver);
  }
  else if (solver_flag == 102) {
    HYPRE_ParCSRGMRESDestroy(solver);
  }
  else if (solver_flag == 103) {
    HYPRE_SStructGMRESDestroy(sstruct_solver);
  }
  else if (solver_flag == 1002) {
    HYPRE_ParCSRPCGDestroy(solver);
  }
  else if (solver_flag == 1003) {
    HYPRE_SStructPCGDestroy(sstruct_solver);
  }
  else if (solver_flag == 104 || solver_flag == 105) {
    HYPRE_ParCSRGMRESDestroy(solver);
    HYPRE_BoomerAMGDestroy(precond);
  }
  else if (solver_flag == 106) {
    HYPRE_SStructSplitDestroy(sstruct_solver);
  }
  else if (solver_flag == 107) {
    HYPRE_SStructGMRESDestroy(sstruct_solver);
    HYPRE_SStructSplitDestroy(sstruct_precond);
  }
  else if (solver_flag == 108) {
    ParmParse pp("hmabec");
#if (AMREX_SPACEDIM == 1)
    int struct_flag = 0;
#else
    int struct_flag = 1;
#endif
    pp.query("struct_flag", struct_flag);
    if (struct_flag == 0) {
      HYPRE_StructSMGDestroy((HYPRE_StructSolver) solver);
    }
    else {
      HYPRE_StructPFMGDestroy((HYPRE_StructSolver) solver);
    }
  }
  else if (solver_flag == 109) {
    ParmParse pp("hmabec");
#if (AMREX_SPACEDIM == 1)
    int struct_flag = 0;
#else
    int struct_flag = 1;
#endif
    pp.query("struct_flag", struct_flag);
    HYPRE_StructGMRESDestroy((HYPRE_StructSolver) solver);
    if (struct_flag == 0) {
      HYPRE_StructSMGDestroy((HYPRE_StructSolver) precond);
    }
    else {
      HYPRE_StructPFMGDestroy((HYPRE_StructSolver) precond);
    }
  }

  if (solver_flag == 150) {
    HYPRE_BoomerAMGDestroy(solver);
  }
  else if (solver_flag == 151 || solver_flag == 153) {
    HYPRE_ParCSRPCGDestroy(solver);
    HYPRE_BoomerAMGDestroy(precond);
  }
  else if (solver_flag == 152) {
    HYPRE_SStructPCGDestroy((HYPRE_SStructSolver) solver);
    HYPRE_SStructSplitDestroy((HYPRE_SStructSolver) precond);
  }

  sstruct_solver = NULL;
  solver         = NULL;
  precond        = NULL;
}

void HypreMultiABec::solve()
{
  BL_PROFILE("HypreMultiABec::solve");

  if (abstol > 0.0) {
    Real bnorm;
    hypre_SStructInnerProd((hypre_SStructVector *) b,
                           (hypre_SStructVector *) b,
                           &bnorm);
    bnorm = std::sqrt(bnorm);

    Real volume = 0.0;
    for (int level = crse_level; level <= fine_level; level++) {
      for (int i = 0; i < grids[level].size(); i++) {
        volume += grids[level][i].numPts();
      }
    }

    Real reltol_new = (bnorm > 0.0
                       ? abstol / bnorm * std::sqrt(volume)
                       : reltol);

    if (reltol_new > reltol) {
      if (solver_flag == 100) {
        HYPRE_BoomerAMGSetTol(solver, reltol_new);
      }
      else if(solver_flag == 101) {
        HYPRE_SStructFACSetTol(sstruct_solver, reltol_new);
      }
      else if(solver_flag == 102) {
        HYPRE_ParCSRGMRESSetTol(solver, reltol_new);
      }
      else if (solver_flag == 103 || solver_flag == 107) {
        HYPRE_SStructGMRESSetTol(sstruct_solver, reltol_new);
      }
      else if(solver_flag == 1002) {
        HYPRE_ParCSRPCGSetTol(solver, reltol_new);
      }
      else if (solver_flag == 1003) {
        HYPRE_SStructPCGSetTol(sstruct_solver, reltol_new);
      }
      else if (solver_flag == 104 || solver_flag == 105) {
        HYPRE_ParCSRGMRESSetTol(solver, reltol_new);
        HYPRE_BoomerAMGSetTol(precond, reltol_new);
      }
      else if (solver_flag == 106) {
        HYPRE_SStructSplitSetTol(sstruct_solver, reltol_new);
      }
      else if (solver_flag == 108) {
        ParmParse pp("hmabec");
#if (AMREX_SPACEDIM == 1)
        int struct_flag = 0;
#else
        int struct_flag = 1;
#endif
        pp.query("struct_flag", struct_flag);
        if (struct_flag == 0) {
          HYPRE_StructSMGSetTol((HYPRE_StructSolver) solver, reltol_new);
        }
        else {
          HYPRE_StructPFMGSetTol((HYPRE_StructSolver) solver, reltol_new);
        }
      }
      else if (solver_flag == 109) {
        HYPRE_StructGMRESSetTol((HYPRE_StructSolver) solver, reltol_new);
      }
      else if (solver_flag == 150) {
        HYPRE_BoomerAMGSetTol(solver, reltol_new);
      }
      else if(solver_flag == 151 || solver_flag == 153) {
        HYPRE_PCGSetTol(solver, reltol_new);
      }
      else if(solver_flag == 152) {
        HYPRE_PCGSetTol(solver, reltol_new);
      }
    }
  }

  if (solver_flag == 100) {
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);
    //std::cout << "HypreMultiABec: starting solver..." << std::endl;
    HYPRE_BoomerAMGSolve(solver, par_A, par_b, par_x);
    //std::cout << "                                  done." << std::endl;
    //HYPRE_SStructVectorPrint("Xamg", x, 0);
    //HYPRE_SStructVectorPrint("Bamg", b, 0);
    //cin.get();
  }
  else if (solver_flag == 101) {
    HYPRE_SStructFACSolve3(sstruct_solver, A, b, x);
  }
  else if (solver_flag == 102) {
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);

    HYPRE_ParCSRGMRESSolve(solver, par_A, par_b, par_x);
  }
  else if (solver_flag == 103 || solver_flag == 107) {
    HYPRE_SStructGMRESSolve(sstruct_solver, A, b, x);
  }
  else if (solver_flag == 1002) {
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);

    HYPRE_ParCSRPCGSolve(solver, par_A, par_b, par_x);
  }
  else if (solver_flag == 1003) {
    HYPRE_SStructPCGSolve(sstruct_solver, A, b, x);
  }
  else if (solver_flag == 104 || solver_flag == 105) {
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);

    HYPRE_ParCSRGMRESSolve(solver, par_A, par_b, par_x);
  }
  else if (solver_flag == 106) {
    HYPRE_SStructSplitSolve(sstruct_solver, A, b, x);
  }
  else if (solver_flag == 108) {
    HYPRE_StructMatrix s_A;
    HYPRE_StructVector s_b;
    HYPRE_StructVector s_x;
    HYPRE_SStructMatrixGetObject(A, (void**) &s_A);
    HYPRE_SStructVectorGetObject(b, (void**) &s_b);
    HYPRE_SStructVectorGetObject(x, (void**) &s_x);
    ParmParse pp("hmabec");
#if (AMREX_SPACEDIM == 1)
    int struct_flag = 0;
#else
    int struct_flag = 1;
#endif
    pp.query("struct_flag", struct_flag);
    if (struct_flag == 0) {
      HYPRE_StructSMGSolve((HYPRE_StructSolver) solver, s_A, s_b, s_x);
    }
    else {
      HYPRE_StructPFMGSolve((HYPRE_StructSolver) solver, s_A, s_b, s_x);
    }
  }
  else if (solver_flag == 109) {
    HYPRE_StructMatrix s_A;
    HYPRE_StructVector s_b;
    HYPRE_StructVector s_x;
    HYPRE_SStructMatrixGetObject(A, (void**) &s_A);
    HYPRE_SStructVectorGetObject(b, (void**) &s_b);
    HYPRE_SStructVectorGetObject(x, (void**) &s_x);
    HYPRE_StructGMRESSolve((HYPRE_StructSolver) solver, s_A, s_b, s_x);
  }
  else if (solver_flag == 150) {
    int num_iterations, myid;
    Real res;

    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);
    //std::cout << "HypreMultiABec: starting solver..." << std::endl;
    HYPRE_BoomerAMGSolve(solver, par_A, par_b, par_x);
    HYPRE_SStructVectorGather(x);

    HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &res);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid );
    if (!myid)
    {
       std::cout << num_iterations << " Hypre Multigrid Iterations_inside, Relative Residual "
            << res << std::endl;
    }
    //std::cout << "                                  done." << std::endl;
    //HYPRE_SStructVectorPrint("Xamg", x, 0);
    //HYPRE_SStructVectorPrint("Bamg", b, 0);
    //cin.get();
  }
  else if (solver_flag == 151 || solver_flag == 153) {
    int num_iterations, myid;
    Real res;

    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);
    HYPRE_PCGSolve(solver, (HYPRE_Matrix)par_A,
                   (HYPRE_Vector)par_b, (HYPRE_Vector)par_x);
    HYPRE_SStructVectorGather(x);

    HYPRE_PCGGetNumIterations(solver, &num_iterations);
    HYPRE_PCGGetFinalRelativeResidualNorm(solver, &res);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid );
    if (!myid)
    {
       std::cout << num_iterations << " Hypre Multigrid Iterations_inside, Relative Residual "
            << res << std::endl;
    }
  }
  else if (solver_flag == 152) {
    int num_iterations, myid;
    Real res;

    HYPRE_PCGSolve(solver, (HYPRE_Matrix) A,
                   (HYPRE_Vector) b, (HYPRE_Vector) x);
    HYPRE_SStructVectorGather(x);

    HYPRE_PCGGetNumIterations(solver, &num_iterations);
    HYPRE_PCGGetFinalRelativeResidualNorm(solver, &res);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid );
    if (!myid)
    {
       std::cout << num_iterations << " Hypre Multigrid Iterations_inside, Relative Residual "
            << res << std::endl;
    }
  }

  HYPRE_SStructVectorGather(x);

  if (verbose >= 2 && ParallelDescriptor::IOProcessor()) {
    int num_iterations;
    Real res;
    if (solver_flag == 100) {
      HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &res);
    }
    else if (solver_flag == 101) {
      HYPRE_SStructFACGetNumIterations(sstruct_solver, &num_iterations);
      HYPRE_SStructFACGetFinalRelativeResidualNorm(sstruct_solver, &res);
    }
    else if (solver_flag == 102) {
      HYPRE_ParCSRGMRESGetNumIterations(solver, &num_iterations);
      HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver, &res);
    }
    else if (solver_flag == 103 || solver_flag == 107) {
      HYPRE_SStructGMRESGetNumIterations(sstruct_solver, &num_iterations);
      HYPRE_SStructGMRESGetFinalRelativeResidualNorm(sstruct_solver, &res);
    }
    else if (solver_flag == 1002) {
      HYPRE_ParCSRPCGGetNumIterations(solver, &num_iterations);
      HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(solver, &res);
    }
    else if (solver_flag == 1003) {
      HYPRE_SStructPCGGetNumIterations(sstruct_solver, &num_iterations);
      HYPRE_SStructPCGGetFinalRelativeResidualNorm(sstruct_solver, &res);
    }
    else if (solver_flag == 104 || solver_flag == 105) {
      HYPRE_ParCSRGMRESGetNumIterations(solver, &num_iterations);
      HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver, &res);
    }
    else if (solver_flag == 106) {
      HYPRE_SStructSplitGetNumIterations(sstruct_solver, &num_iterations);
      HYPRE_SStructSplitGetFinalRelativeResidualNorm(sstruct_solver, &res);
    }
    else if (solver_flag == 108) {
      ParmParse pp("hmabec");
#if (AMREX_SPACEDIM == 1)
      int struct_flag = 0;
#else
      int struct_flag = 1;
#endif
      pp.query("struct_flag", struct_flag);
      HYPRE_StructSolver& struct_solver = *(HYPRE_StructSolver*)&solver;
      if (struct_flag == 0) {
        HYPRE_StructSMGGetNumIterations(struct_solver, &num_iterations);
        HYPRE_StructSMGGetFinalRelativeResidualNorm(struct_solver, &res);
      }
      else {
        HYPRE_StructPFMGGetNumIterations(struct_solver, &num_iterations);
        HYPRE_StructPFMGGetFinalRelativeResidualNorm(struct_solver, &res);
      }
    }
    else if (solver_flag == 109) {
      HYPRE_StructSolver& struct_solver = *(HYPRE_StructSolver*)&solver;
      HYPRE_StructGMRESGetNumIterations(struct_solver, &num_iterations);
      HYPRE_StructGMRESGetFinalRelativeResidualNorm(struct_solver, &res);
    }
    else if (solver_flag == 150) {
      HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &res);
    }
    else if (solver_flag == 151 || solver_flag == 152 || solver_flag == 153) {
      HYPRE_PCGGetNumIterations(solver, &num_iterations);
      HYPRE_PCGGetFinalRelativeResidualNorm(solver, &res);
    }

    if (num_iterations >= verbose_threshold) {
      int oldprec = std::cout.precision(20);
      if (Radiation::current_group_number >= 0) {
        std::cout << Radiation::current_group_name << " Group "
             << Radiation::current_group_number << ": ";
      }
      std::cout << num_iterations
           << " Hypre Multigrid Iterations, Relative Residual "
           << res << std::endl;
      std::cout.precision(oldprec);
    }

  }
}

void HypreMultiABec::getSolution(int level, MultiFab& dest, int icomp)
{
  BL_PROFILE("HypreMultiABec::getSolution");

  int part = level - crse_level;

  FArrayBox fnew;
  for (MFIter mfi(dest); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = grids[level][i];

    FArrayBox *f;
    int fcomp;
    if (dest.nGrow() == 0) { // need a temporary if dest is the wrong size
      f = &dest[mfi];
      fcomp = icomp;
    }
    else {
      f = &fnew;
      f->resize(reg);

      fcomp = 0;
    }
    Elixir f_elix = fnew.elixir();

    vectorGetBoxValues(x, part, reg, subgrids[level][i], *f, fcomp);

    if (dest.nGrow() != 0) {
        dest[mfi].copy<RunOn::Device>(*f, 0, icomp, 1);
    }
  }
}

Real HypreMultiABec::getAbsoluteResidual()
{
  BL_PROFILE("HypreMultiABec::getAbsoluteResidual");

  Real bnorm;
  hypre_SStructInnerProd((hypre_SStructVector *) b,
                         (hypre_SStructVector *) b,
                         &bnorm);
  bnorm = std::sqrt(bnorm);

  Real res;
  if (solver_flag == 100) {
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &res);
  }
  else if (solver_flag == 101) {
    HYPRE_SStructFACGetFinalRelativeResidualNorm(sstruct_solver, &res);
  }
  else if (solver_flag == 102) {
    HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver, &res);
  }
  else if (solver_flag == 103 || solver_flag == 107) {
    HYPRE_SStructGMRESGetFinalRelativeResidualNorm(sstruct_solver, &res);
  }
  else if (solver_flag == 1002) {
    HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(solver, &res);
  }
  else if (solver_flag == 1003) {
    HYPRE_SStructPCGGetFinalRelativeResidualNorm(sstruct_solver, &res);
  }
  else if (solver_flag == 104 || solver_flag == 105) {
    HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver, &res);
  }
  else if (solver_flag == 106) {
    HYPRE_SStructSplitGetFinalRelativeResidualNorm(sstruct_solver, &res);
  }
  else if (solver_flag == 108) {
    ParmParse pp("hmabec");
#if (AMREX_SPACEDIM == 1)
    int struct_flag = 0;
#else
    int struct_flag = 1;
#endif
    pp.query("struct_flag", struct_flag);
    HYPRE_StructSolver& struct_solver = *(HYPRE_StructSolver*)&solver;
    if (struct_flag == 0) {
      HYPRE_StructSMGGetFinalRelativeResidualNorm(struct_solver, &res);
    }
    else {
      HYPRE_StructPFMGGetFinalRelativeResidualNorm(struct_solver, &res);
    }
  }
  else if (solver_flag == 109) {
    HYPRE_StructSolver& struct_solver = *(HYPRE_StructSolver*)&solver;
    HYPRE_StructGMRESGetFinalRelativeResidualNorm(struct_solver, &res);
  }
  else if (solver_flag == 150) {
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &res);
  }
  else if (solver_flag == 151 || solver_flag == 152 || solver_flag == 153) {
    HYPRE_PCGGetFinalRelativeResidualNorm(solver, &res);
  }

  Real volume = 0.0;
  for (int level = crse_level; level <= fine_level; level++) {
    for (int i = 0; i < grids[level].size(); i++) {
      volume += grids[level][i].numPts();
    }
  }

  return bnorm * res / std::sqrt(volume);
}

void HypreMultiABec::boundaryFlux(int level,
                                  MultiFab* Flux,
                                  MultiFab& Soln,
                                  int icomp,
                                  BC_Mode inhom)
{
    BL_PROFILE("HypreMultiABec::boundaryFlux");

    const Box& domain = bd[level]->getDomain();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(Soln); mfi.isValid(); ++mfi) {
            int i = mfi.index();
            const Box &reg = grids[level][i];
            for (OrientationIter oitr; oitr; oitr++) {
                int cdir(oitr());
                int idim = oitr().coordDir();
                const RadBoundCond &bct = bd[level]->bndryConds(oitr())[i];
                const Real      &bcl = bd[level]->bndryLocs(oitr())[i];
                const FArrayBox       &fs  = bd[level]->bndryValues(oitr())[mfi];
                const Mask      &msk = bd[level]->bndryMasks(oitr(), i);
                const Box &fbox = Flux[idim][mfi].box();
                const Box &sbox = Soln[mfi].box();
                const Box &msb  = msk.box();
                const Box &bbox = (*bcoefs[level])[idim][mfi].box();
                if (reg[oitr()] == domain[oitr()]) {
                    int bctype = bct;
                    Array4<int const> tf_arr;
                    if (bd[level]->mixedBndry(oitr())) {
                        const BaseFab<int> &tf = *(bd[level]->bndryTypes(oitr())[i]);
                        tf_arr = tf.array();
                        bctype = -1;
                    }
                    // In normal code operation only the fluxes at internal
                    // Dirichlet boundaries are used.  Some diagnostics use the
                    // fluxes computed at domain boundaries but these do not
                    // influence the evolution of the interior solution.
                    Array4<Real const> sp_arr;
                    if (SPa[level]) {
                        sp_arr = (*SPa[level])[mfi].array();
                    }

                    HABEC::hbflx3(Flux[idim][mfi].array(),
                                  Soln[mfi].array(icomp),
                                  reg,
                                  cdir, bctype,
                                  tf_arr,
                                  bho, bcl,
                                  fs.array(bdcomp),
                                  msk.array(),
                                  (*bcoefs[level])[idim][mfi].array(),
                                  beta, geom[level].CellSize(),
                                  flux_factor, oitr(),
                                  geom[level].data(), inhom,
                                  sp_arr);
                }
                else {
                    HABEC::hbflx(Flux[idim][mfi].array(),
                                 Soln[mfi].array(icomp),
                                 reg,
                                 cdir, bct, bho, bcl,
                                 fs.array(bdcomp),
                                 msk.array(),
                                 (*bcoefs[level])[idim][mfi].array(),
                                 beta, geom[level].CellSize(), inhom);
                }
            }
        }
    }
}

void HypreMultiABec::getProduct(int level, MultiFab& product)
{
  BL_PROFILE("HypreMultiABec::getProduct");

  int part = level - crse_level;

  for (MFIter mfi(product); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = grids[level][i];

    vectorGetBoxValues(b, part, reg, subgrids[level][i], product[i], 0);
  }
}
