
#include <HypreExtMultiABec.H>
#include <HABEC.H>
#include <AMReX_LO_BCTYPES.H>

#include <_hypre_sstruct_mv.h>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

HypreExtMultiABec::~HypreExtMultiABec()
{
}

void HypreExtMultiABec::a2Coefficients(int level, const MultiFab &a2, int dir)
{
  BL_PROFILE("HypreExtMultiABec::a2Coefficients");

  BL_ASSERT( a2.ok() );

  int ncomp=1;
  int ngrow=0;

  if (!a2coefs[level]) {
    a2coefs[level].reset(new Array<MultiFab, AMREX_SPACEDIM>);

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      BoxArray edge_boxes(grids[level]);
      edge_boxes.surroundingNodes(i);
      (*a2coefs[level])[i].define(edge_boxes, dmap[level], ncomp, ngrow);
      (*a2coefs[level])[i].setVal(0.0);
    }
  }

  BL_ASSERT( a2.boxArray() == (*a2coefs[level])[dir].boxArray() );

  MultiFab::Copy((*a2coefs[level])[dir], a2, 0, 0, ncomp, ngrow);
}

void HypreExtMultiABec::cCoefficients(int level, const MultiFab &c, int dir)
{
  BL_PROFILE("HypreExtMultiABec::cCoefficients");

  BL_ASSERT( c.ok() );

  int ncomp=2; // coeffs are 2-sided for upwinding
  int ngrow=0;

  if (!ccoefs[level]) {
    ccoefs[level].reset(new Array<MultiFab, AMREX_SPACEDIM>);

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      BoxArray edge_boxes(grids[level]);
      edge_boxes.surroundingNodes(i);
      (*ccoefs[level])[i].define(edge_boxes, dmap[level], ncomp, ngrow);
      (*ccoefs[level])[i].setVal(0.0);
    }
  }

  BL_ASSERT( c.boxArray() == (*ccoefs[level])[dir].boxArray() );

  MultiFab::Copy((*ccoefs[level])[dir], c, 0, 0, ncomp, ngrow);
}

void HypreExtMultiABec::d1Coefficients(int level, const MultiFab &d1, int dir)
{
  BL_PROFILE("HypreExtMultiABec::d1Coefficients");

  BL_ASSERT( d1.ok() );

  int ncomp=1;
  int ngrow=0;

  if (!d1coefs[level]) {
    d1coefs[level].reset(new Array<MultiFab, AMREX_SPACEDIM>);

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      (*d1coefs[level])[i].define(grids[level], dmap[level], ncomp, ngrow);
      (*d1coefs[level])[i].setVal(0.0);
    }
  }

  BL_ASSERT( d1.boxArray() == (*d1coefs[level])[dir].boxArray() );

  MultiFab::Copy((*d1coefs[level])[dir], d1, 0, 0, ncomp, ngrow);
}

void HypreExtMultiABec::d2Coefficients(int level, const MultiFab &d2, int dir)
{
  BL_PROFILE("HypreExtMultiABec::d2Coefficients");

  BL_ASSERT( d2.ok() );

  int ncomp=1;
  int ngrow=0;

  if (!d2coefs[level]) {
    d2coefs[level].reset(new Array<MultiFab, AMREX_SPACEDIM>);

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
      BoxArray edge_boxes(grids[level]);
      edge_boxes.surroundingNodes(i);
      (*d2coefs[level])[i].define(edge_boxes, dmap[level], ncomp, ngrow);
      (*d2coefs[level])[i].setVal(0.0);
    }
  }

  BL_ASSERT( d2.boxArray() == (*d2coefs[level])[dir].boxArray() );

  MultiFab::Copy((*d2coefs[level])[dir], d2, 0, 0, ncomp, ngrow);
}

static void
FaceValue(AuxVarBox& evalue, AuxVarBox& cintrp,
          const Mask& msk, const Box& reg,
          const IntVect& vin, int r, int bho, int flevel)
{
  if (bho == 1) {
    Real efacb = 3.0 / ((1 + r) * (3 + r));
    Real efac1 = ( 1.5 * r) / (1 + r);
    Real efac2 = (-0.5 * r) / (3 + r);
    IntVect vi2 = 2 * vin;
    for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
      if (msk(v) == RadBndryData::not_covered) {
        evalue(v).push(&cintrp(v),    efacb);
        evalue(v).push(flevel, v-vin, efac1);
        evalue(v).push(flevel, v-vi2, efac2);
      }
    }
  }
  else {
    Real efacb = 1.0 / (1 + r);
    Real efac1 = r / (1.0 + r);
    for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
      if (msk(v) == RadBndryData::not_covered) {
        evalue(v).push(&cintrp(v),    efacb);
        evalue(v).push(flevel, v-vin, efac1);
      }
    }
  }
}

void HypreExtMultiABec::loadMatrix()
{
  BL_PROFILE("HypreExtMultiABec::loadMatrix");

  HypreMultiABec::loadMatrix();

  if (0 && verbose >= 1 && ParallelDescriptor::IOProcessor()) {
    std::cout << "In HypreExtMultiABec::loadMatrix(), the multipliers are:" << std::endl;
    std::cout << "  HypreExtMultiABec::alpha2 = " << alpha2 << std::endl;
    std::cout << "  HypreExtMultiABec::gamma  = " << gamma  << std::endl;
    std::cout << "  HypreExtMultiABec::delta1 = " << delta1 << std::endl;
    std::cout << "  HypreExtMultiABec::delta2 = " << delta2 << std::endl;
  }

  // These really ought to be members and already defined:
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

  int stencil_indices[size];

  for (int i = 0; i < size; i++) {
    stencil_indices[i] = i;
  }

  FArrayBox matfab;
  FArrayBox mat_tmpfab;
  for (int level = crse_level; level <= fine_level; level++) {
    int part = level - crse_level;

    for (MFIter mfi(*acoefs[level]); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      const Box &reg = grids[level][i];

      matfab.resize(reg,size);
      auto mat = matfab.array();
      Elixir mat_elix = matfab.elixir();

      matfab.setVal<RunOn::Device>(0.0);

      // build matrix interior

      if (a2coefs[level]) {
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            const Real fac = 0.25_rt * alpha2;

            auto a2 = (*a2coefs[level])[idim][mfi].array();

            amrex::ParallelFor(reg,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {
                if (idim == 0) {
                    mat(i,j,k,0) += fac * (a2(i,j,k) + a2(i+1,j,k));
                    mat(i,j,k,1) += fac * a2(i,j,k);
                    mat(i,j,k,2) += fac * a2(i+1,j,k);
                }
                else if (idim == 1) {
                    mat(i,j,k,0) += fac * (a2(i,j,k) + a2(i,j+1,k));
                    mat(i,j,k,3) += fac * a2(i,j,k);
                    mat(i,j,k,4) += fac * a2(i,j+1,k);
                }
                else {
                    mat(i,j,k,0) += fac * (a2(i,j,k) + a2(i,j,k+1));
                    mat(i,j,k,5) += fac * a2(i,j,k);
                    mat(i,j,k,6) += fac * a2(i,j,k+1);
                }
            });
        }
      }

      if (ccoefs[level]) {
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            const Real fac = 0.5_rt * gamma / geom[level].CellSize(idim);

            auto c = (*ccoefs[level])[idim][mfi].array();

            amrex::ParallelFor(reg,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {
                if (idim == 0) {
                    mat(i,j,k,0) += -fac * (c(i,j,k) - c(i+1,j,k));
                    mat(i,j,k,1) += -fac * c(i,j,k);
                    mat(i,j,k,2) += fac * c(i+1,j,k);
                }
                else if (idim == 1) {
                    mat(i,j,k,0) += -fac * (c(i,j,k) - c(i,j+1,k));
                    mat(i,j,k,3) += -fac * c(i,j,k);
                    mat(i,j,k,4) += fac * c(i,j+1,k);
                }
                else {
                    mat(i,j,k,0) += -fac * (c(i,j,k) - c(i,j,k+1));
                    mat(i,j,k,5) += -fac * c(i,j,k);
                    mat(i,j,k,6) += fac * c(i,j,k+1);
                }
            });
        }
      }

      if (d1coefs[level]) {
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            const Real fac = 0.5_rt * delta1 / geom[level].CellSize(idim);

            auto d1 = (*d1coefs[level])[idim][mfi].array();

            amrex::ParallelFor(reg,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {
                if (idim == 0) {
                    mat(i,j,k,1) += -fac * d1(i,j,k);
                    mat(i,j,k,2) += fac * d1(i,j,k);
                }
                else if (idim == 1) {
                    mat(i,j,k,3) += -fac * d1(i,j,k);
                    mat(i,j,k,4) += fac * d1(i,j,k);
                }
                else {
                    mat(i,j,k,5) += -fac * d1(i,j,k);
                    mat(i,j,k,6) += fac * d1(i,j,k);
                }
            });
        }
      }

      if (d2coefs[level]) {
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            const Real fac = 0.5_rt * delta2 / geom[level].CellSize(idim);

            auto d2 = (*d2coefs[level])[idim][mfi].array();

            amrex::ParallelFor(reg,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {
                if (idim == 0) {
                    mat(i,j,k,0) += fac * (d2(i,j,k) - d2(i+1,j,k));
                    mat(i,j,k,1) += -fac * d2(i,j,k);
                    mat(i,j,k,2) += fac * d2(i+1,j,k);
                }
                else if (idim == 1) {
                    mat(i,j,k,0) += fac * (d2(i,j,k) - d2(i,j+1,k));
                    mat(i,j,k,3) += -fac * d2(i,j,k);
                    mat(i,j,k,4) += fac * d2(i,j+1,k);
                }
                else {
                    mat(i,j,k,0) += fac * (d2(i,j,k) - d2(i,j,k+1));
                    mat(i,j,k,5) += -fac * d2(i,j,k);
                    mat(i,j,k,6) += fac * d2(i,j,k+1);
                }
            });
        }
      }

      Gpu::streamSynchronize();

      // Boundary conditions will be corrected below.

      // update matrix
      mat_tmpfab.resize(reg);
      Real* mat_tmp = mat_tmpfab.dataPtr();
      int volume = reg.numPts();
      for (int s = 0; s < size; s++) {
        for (int k = 0; k < volume; k++) {
          mat_tmp[k] = matfab.dataPtr()[s * volume + k];
        }
        HYPRE_SStructMatrixAddToBoxValues(A, part, loV(reg), hiV(reg), 0,
                                          1, &stencil_indices[s], mat_tmp);
      }
    }

    Gpu::streamSynchronize();

    // At this point we begin adding matrix entries for points around
    // the edges of the current level.  This now includes boundary
    // conditions as well as the coarse-fine interface discretization
    // as seen from the fine side of each interface.

    // (Later, we'll add entries corresponding to the coarse side of
    // each coarse-fine interface.  That will be a separate section
    // because the action will be happening on a different set of
    // processors.)

    BndryAuxVar evalue(grids[level], dmap[level], BndryAuxVar::GHOST);
    BndryAuxVar entry( grids[level], dmap[level], BndryAuxVar::INTERIOR);
    if (level == crse_level) {
      // HypreMultiABec creates this for finer levels but not this one.
      // We will delete it manually below.
      ederiv[level].reset(new BndryAuxVar(grids[level], dmap[level], BndryAuxVar::GHOST));
    }

    // Add matrix entries corresponding to physical boundary conditions
    // and all exposed outer edges of crse_level:

    const Box& domain = bd[level]->getDomain();
    for (OrientationIter oitr; oitr; oitr++) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      IntVect vin = amrex::BASISV(idir), ves;
      vin = (ori.isLow() ? -vin : vin);    // outward normal unit vector
      ves = (ori.isLow() ?  ves : vin);    // edge shift vector
      Real h    = geom[level].CellSize(idir); // normal fine grid spacing
      Real h2   = 0.5 * h;
      Real th2  = 1.5 * h;
      Real ofac = (ori.isLow() ? -1.0 : 1.0);
      Real ofh  = ofac / h;
      Real of2h = ofac * 0.5 / h;
      for (int i = entry.firstLocal(); entry.isValid(i);
           i = entry.nextLocal(i)) {
        Box reg = amrex::adjCell(grids[level][i], ori);
        if (grids[level][i][ori] == domain[ori] || level == crse_level) {
          RadBoundCond     bct = bd[level]->bndryConds(ori)[i];
          // bct may be changed below if this is a mixed boundary
          const Real      &bcl = bd[level]->bndryLocs(ori)[i];
          const Mask      &msk = bd[level]->bndryMasks(ori,i);
//          const Box &bbox = (*bcoefs[level])[idir][i].box();
//          const Box &msb  = msk.box();

          // Treat an exposed grid edge here as a boundary condition
          // for the linear solver:

          reg.shift(-vin); // fine interior cells
          for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
            if (msk(v+vin) != RadBndryData::covered) {

              // zero out interior stencil since this is a boundary cell
              if (a2coefs[level]) {
                Real fac = alpha2 * (*a2coefs[level])[idir][i](v+ves);
                entry(ori,i)(v).push(level, v,     -0.25 * fac);
                entry(ori,i)(v).push(level, v+vin, -0.25 * fac);
              }
              if (ccoefs[level]) {
                // upwinding is determined by the components of ccoefs
                Real fac = gamma * ofh;
                int  i0  = ori.isLow();
                int  i1  = ori.isHigh();
                entry(ori,i)(v).push(level, v,
                                      -fac * (*ccoefs[level])[idir][i](v+ves,i0));
                entry(ori,i)(v).push(level, v+vin,
                                      -fac * (*ccoefs[level])[idir][i](v+ves,i1));
              }
              if (d1coefs[level]) {
                Real fac = delta1 * of2h * (*d1coefs[level])[idir][i](v);
                entry(ori,i)(v).push(level, v,      fac);
                entry(ori,i)(v).push(level, v+vin, -fac);
              }
              if (d2coefs[level]) {
                Real fac = delta2 * of2h * (*d2coefs[level])[idir][i](v+ves);
                entry(ori,i)(v).push(level, v,      fac);
                entry(ori,i)(v).push(level, v+vin, -fac);
              }

              // determine what type of boundary this is and act accordingly
              if (reg[ori] == domain[ori] && bd[level]->mixedBndry(ori)) {
                bct = (*(bd[level]->bndryTypes(ori)[i]))(v+vin);
              }
              if (bct == AMREX_LO_DIRICHLET) {
                if (bho == 1) {
                  evalue(ori,i)(v+vin).push(level, v,
                                            (bcl * th2) / (h * (bcl + h2)));
                  evalue(ori,i)(v+vin).push(level, v-vin,
                                            (-bcl * h2) / (h * (bcl + th2)));
                  (*ederiv[level])(ori,i)(v+vin).push(level, v,
                                            (bcl - th2) / (h*(bcl + h2)));
                  (*ederiv[level])(ori,i)(v+vin).push(level, v-vin,
                                            (h2 - bcl) / (h*(bcl + th2)));
                }
                else {
                  evalue(ori,i)(v+vin).push(level, v,
                                             bcl / (bcl + h2));
                  (*ederiv[level])(ori,i)(v+vin).push(level, v,
                                            -1.0 / (bcl + h2));
                }
                if (a2coefs[level]) {
                  Real fac = 0.5 * alpha2 * (*a2coefs[level])[idir][i](v+ves);
                  entry(ori,i)(v).push(&evalue(ori,i)(v+vin), fac);
                }
                if (ccoefs[level]) {

                  // This was the centered difference version:
                  //Real fac = gamma * ofh * (*ccoefs[level])[idir][i](v+ves);
                  //entry(ori,i)(v).push(&evalue(ori,i)(v+vin), fac);

                  // For upwinding, we can't use evalue since it's centered.
                  // The good side is we don't need to worry about being
                  // high-order, since it's a low-order stencil anyway.
                  // If we're upwinding from the interior, just use that value.
                  // If we're upwinding from the exterior, interpolate
                  // (linearly) to the ghost cell center position.

                  Real fac = gamma * ofh;
                  int  i0  = ori.isLow();
                  int  i1  = ori.isHigh();
                  // upwinding from interior:
                  entry(ori,i)(v).push(level, v,
                                        fac
                                        * (*ccoefs[level])[idir][i](v+ves,i0));
                  // upwinding from exterior:
                  entry(ori,i)(v).push(level, v,
                                        fac * (bcl - h2) / (bcl + h2)
                                        * (*ccoefs[level])[idir][i](v+ves,i1));
                }
                if (d1coefs[level]) {
                  Real fac = 0.5*delta1*ofac * (*d1coefs[level])[idir][i](v);
                  entry(ori,i)(v).push(&(*ederiv[level])(ori,i)(v+vin), fac);
                }
                if (d2coefs[level]) {
                  Real fac = 0.5*delta2*ofac * (*d2coefs[level])[idir][i](v+ves);
                  entry(ori,i)(v).push(&(*ederiv[level])(ori,i)(v+vin), fac);
                }
              }
              else if (bct == AMREX_LO_NEUMANN) {
                // no more action required here
              }
              else if (bct == AMREX_LO_MARSHAK || bct == AMREX_LO_SANCHEZ_POMRANING) {
                if (bho == 1) {
                  evalue(ori,i)(v+vin).push(level, v,      1.5);
                  evalue(ori,i)(v+vin).push(level, v-vin, -0.5);
                }
                else {
                  evalue(ori,i)(v+vin).push(level, v,      1.0);
                }
                // cmult should be c for photons, 1 for neutrinos
                Real cmult = flux_factor;
                Real xi = 0.0; // xi should be passed in through RadBndry?
                //Real tmp = beta / h; // already done in HypreMultiABec
                Real tmp = 0.0;
                Real fac = 0.5 * (tmp - 0.5 * (delta1 + delta2) * ofac * xi);
                entry(ori,i)(v).push(&evalue(ori,i)(v+vin), cmult * fac);
              }
              else {
                amrex::Error("HypreExtMultiABec: unsupported boundary type");
              }
            }
          }
        }
      }
    }

    // If there is a coarser active level, add matrix entries here
    // that correspond to the fine side of the interface to the next
    // coarser level:

    if (level > crse_level) {

      // First we do the entries as seen by the fine cells adjacent
      // to the interface, working on the fine processor:

      IntVect rat = fine_ratio[level-1];
      for (OrientationIter oitr; oitr; ++oitr) {
        Orientation ori = oitr();
        int idir = ori.coordDir();
        IntVect vin = amrex::BASISV(idir), ves;
        vin = (ori.isLow() ? -vin : vin);    // outward normal unit vector
        ves = (ori.isLow() ?  ves : vin);    // edge shift vector
        Real h = geom[level].CellSize(idir); // normal fine grid spacing
        Real ofac = (ori.isLow() ? -1.0 : 1.0);
        Real of2h = ofac * 0.5 / h;
        for (int i = cintrp[level]->firstLocal(); cintrp[level]->isValid(i);
             i = cintrp[level]->nextLocal(i)) {
          Box reg = amrex::adjCell(grids[level][i], ori);
          const Mask &msk = bd[level]->bndryMasks(ori,i);
          FaceValue(evalue(ori,i),
                    (*cintrp[level])(ori,i),
                    msk, reg, vin, rat[idir], bho, level);
          reg.shift(-vin); // fine interior cells
          for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
            if (msk(v+vin) == RadBndryData::not_covered) {
              if (a2coefs[level]) {
                Real fac = alpha2 * (*a2coefs[level])[idir][i](v+ves);
                entry(ori,i)(v).push(level, v,     -0.25 * fac);
                entry(ori,i)(v).push(level, v+vin, -0.25 * fac);
                entry(ori,i)(v).push(&evalue(ori,i)(v+vin), 0.5 * fac);
              }
              if (ccoefs[level]) {
                // Remove this assertion when upwinding of ccoefs term
                // is implemented for multilevel linear systems:
                BL_ASSERT(level == crse_level);

                // not done yet
                Real fac = gamma * of2h * (*ccoefs[level])[idir][i](v+ves);
                entry(ori,i)(v).push(level, v,     -fac);
                entry(ori,i)(v).push(level, v+vin, -fac);
                entry(ori,i)(v).push(&evalue(ori,i)(v+vin), 2.0 * fac);
              }
              if (d1coefs[level]) {
                Real fac = delta1 * of2h * (*d1coefs[level])[idir][i](v);
                entry(ori,i)(v).push(level, v,      fac);
                entry(ori,i)(v).push(level, v+vin, -fac);
                entry(ori,i)(v).push(&(*ederiv[level])(ori,i)(v+vin), h * fac);
              }
              if (d2coefs[level]) {
                Real fac = delta2 * of2h * (*d2coefs[level])[idir][i](v+ves);
                entry(ori,i)(v).push(level, v,      fac);
                entry(ori,i)(v).push(level, v+vin, -fac);
                entry(ori,i)(v).push(&(*ederiv[level])(ori,i)(v+vin), h * fac);
              }
            }
          }
        }
      }
    }

    // Now we finish adding the entries from the above two sections.

    // The following section is identical to version in HypreMultiABec
    // except for use of AddTo in place of Set:
    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      IntVect vin = amrex::BASISV(idir);
      vin = (ori.isLow() ? -vin : vin);    // outward normal unit vector
      for (int i = entry.firstLocal(); entry.isValid(i);
           i = entry.nextLocal(i)) {
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
                HYPRE_SStructMatrixAddToValues(A, part, getV1(v), 0,
                                               1, &ientry, &values[j]);
                ientry++;
              }
            }
          }
        }
      }
    }

    if (level == crse_level) {
        ederiv[level].reset();
        continue;
    }

    // Now add the matrix values seen by the coarse cells adjacent
    // to the coarse-fine interface.  These are averages of ederiv
    // over the fine faces making up each coarse face.  We use
    // CrseExtBndryAuxVar now because we need to work on the coarse
    // processor.

    IntVect rat = fine_ratio[level-1];
    const BoxArray& f_fgrids(grids[level]);
    const BoxArray& c_cgrids(grids[level-1]);
    BoxArray f_cgrids(c_cgrids);
    f_cgrids.refine(rat);
    BoxArray c_fgrids(f_fgrids);
    c_fgrids.coarsen(rat);
    //CrseBndryAuxVar c_evalue(f_cgrids, f_fgrids, BndryAuxVar::GHOST);
    //CrseBndryAuxVar c_entry( c_cgrids, c_fgrids, BndryAuxVar::EXTERIOR);
    CrseBndryAuxVar c_evalue(*c_ederiv[level], BndryAuxVar::GHOST);

    int nc = 0, pc[3]; // pc is index translation array
    if (a2coefs[level]) {
      pc[0] = nc++;
    }
    if (ccoefs[level]) {
      // Remove this assertion when upwinding of ccoefs term
      // is implemented for multilevel linear systems:
      BL_ASSERT(level == crse_level);

      // not done yet
      pc[1] = nc++;
    }
    if (d2coefs[level]) {
      pc[2] = nc++;
    }
    c_entry[level]->reinitialize_connections(BndryAuxVar::EXTERIOR);
    c_entry[level]->rebuildFaceData(rat, nc);

    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      if (a2coefs[level]) {
        c_entry[level]->loadFaceData(ori, (*a2coefs[level])[idir], 0, pc[0], 1);
      }
      if (ccoefs[level]) {
        // Remove this assertion when upwinding of ccoefs term
        // is implemented for multilevel linear systems:
        BL_ASSERT(level == crse_level);

        // not done yet
        c_entry[level]->loadFaceData(ori,  (*ccoefs[level])[idir], 0, pc[1], 1);
      }
      if (d2coefs[level]) {
        c_entry[level]->loadFaceData(ori, (*d2coefs[level])[idir], 0, pc[2], 1);
      }
      IntVect vin = amrex::BASISV(idir), ves;
      vin = (ori.isLow() ? -vin : vin); // outward normal unit vector
      ves = (ori.isLow() ? -vin : ves); // edge shift vector (diff from above)
      Real hc = geom[level-1].CellSize(idir); // normal coarse grid spacing
      Real ofac = (ori.isLow() ? -1.0 : 1.0);
      Real ofhc = ofac / hc;
      Real rfac = 1.0;
      IntVect ve; // default constructor initializes to zero
#if (AMREX_SPACEDIM >= 2)
      int jdir = (idir + 1) % AMREX_SPACEDIM;
      ve += (rat[jdir] - 1) * amrex::BASISV(jdir);
      rfac /= rat[jdir]; // will average over fine cells in tangential dir
#endif
#if (AMREX_SPACEDIM == 3)
      int kdir = (idir + 2) % 3;
      ve += (rat[kdir] - 1) * amrex::BASISV(kdir);
      rfac /= rat[kdir]; // will average over fine cells in tangential dir
#endif
      for (int i = c_entry[level]->firstLocal(); c_entry[level]->isValid(i);
           i = c_entry[level]->nextLocal(i)) {
        // parallel loop is tied to coarse grids
        for (int j = 0; j < (*c_entry[level]).size(ori,i); j++) {
          const Box& reg = (*c_cintrp[level])(ori,i,j).box(); // adjacent cells
          const Box& creg = (*c_entry[level])(ori,i,j).box(); // adjacent cells
          const Mask& msk = (*c_cintrp[level]).mask(ori,i,j); // fine mask
          FaceValue(c_evalue(ori,i,j),
                    (*c_cintrp[level])(ori,i,j),
                    msk, reg, vin, rat[idir], bho, level);
          // fcoefs contains a2coefs, then ccoefs, then d2coefs, in order, but
          // only the ones that exist.  pc array does component translation.
          const FArrayBox& fcoefs = c_entry[level]->faceData(ori,i,j);
          for (IntVect vc = creg.smallEnd(); vc <= creg.bigEnd(); creg.next(vc)) {
            IntVect vf = rat * vc;
            vf[idir] = reg.smallEnd(idir); // same as bigEnd(idir)
            Box face(vf, vf + ve);
            if (msk(vf) == RadBndryData::not_covered) {
              // Zero out connection to covered coarse cell:
                (*c_entry[level])(ori,i,j)(vc).push(-1, vc-vin, 0.0);
              if (a2coefs[level]) {
                  (*c_entry[level])(ori,i,j)(vc)
                  .push(level-1, vc,
                        -0.25 * alpha2 * (*a2coefs[level-1])[idir][i](vc+ves));
                // Add fine fluxes over face of coarse cell:
                for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
                    (*c_entry[level])(ori,i,j)(vc)
                        .push(&((c_evalue(ori,i,j))(v)),
                              0.5 * alpha2 * rfac * fcoefs(v+ves, pc[0]));
                }
              }
              if (ccoefs[level]) {
                // not done yet
                  (*c_entry[level])(ori,i,j)(vc)
                      .push(level-1, vc,
                            0.5 * gamma * ofhc * (*ccoefs[level-1])[idir][i](vc+ves));
                // Add fine fluxes over face of coarse cell:
                for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
                    (*c_entry[level])(ori,i,j)(vc)
                        .push(&((c_evalue(ori,i,j))(v)),
                              -gamma * ofhc * rfac * fcoefs(v+ves, pc[1]));
                }
              }
              if (d1coefs[level]) {
                Real d1coef = (*d1coefs[level-1])[idir][i](vc);
                (*c_entry[level])(ori,i,j)(vc)
                  .push(level-1, vc,
                        -0.5 * delta1 * ofhc * d1coef);
                // Add fine fluxes over face of coarse cell:
                for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
                  // note that there are no fd1coefs:
                  (*c_entry[level])(ori,i,j)(vc)
                      .push(&(*c_ederiv[level])(ori,i,j)(v),
                              0.5 * delta1 * ofac * rfac * d1coef);
                }
              }
              if (d2coefs[level]) {
                (*c_entry[level])(ori,i,j)(vc)
                  .push(level-1, vc,
                        -0.5 * delta2 * ofhc * (*d2coefs[level-1])[idir][i](vc+ves));
                // Add fine fluxes over face of coarse cell:
                for (IntVect v = vf; v <= face.bigEnd(); face.next(v)) {
                  (*c_entry[level])(ori,i,j)(vc)
                      .push(&(*c_ederiv[level])(ori,i,j)(v),
                      0.5 * delta2 * ofac * rfac * fcoefs(v+ves, pc[2]));
                }
              }
            }
          }
        }
      }
    }

    // The following section is identical to version in HypreMultiABec
    // except for use of AddTo in place of Set:
    for (OrientationIter oitr; oitr; ++oitr) {
      Orientation ori = oitr();
      int idir = ori.coordDir();
      for (int i = c_entry[level]->firstLocal(); c_entry[level]->isValid(i);
           i = c_entry[level]->nextLocal(i)) {
        // parallel loop is tied to coarse grids
        for (int j = 0; j < (*c_entry[level]).size(ori,i); j++) {
          const Box& reg = (*c_cintrp[level])(ori,i,j).box(); // adjacent cells
          const Box& creg = (*c_entry[level])(ori,i,j).box(); // adjacent cells
          const Mask &msk = (*c_cintrp[level]).mask(ori,i,j); // fine mask
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
                  HYPRE_SStructMatrixAddToValues(A, part-1, getV1(vc), 0,
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

void HypreExtMultiABec::loadLevelVectors(int level,
                                         MultiFab& dest,
                                         int icomp,
                                         MultiFab& rhs, // will not be altered
                                         BC_Mode inhom)
{
  BL_PROFILE("HypreExtMultiABec::loadLevelVectors");

  loadLevelVectorX(level, dest, icomp);
  loadLevelVectorB(level, rhs, inhom);
}

void HypreExtMultiABec::loadLevelVectorB(int level,
                                         MultiFab& rhs, // will not be altered
                                         BC_Mode inhom)
{
  {
    static int first = 1;
    if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
      first = 0;
      std::cout << "In HypreExtMultiABec::loadLevelVectorB(), the multipliers are:"
           << std::endl;
      std::cout << "  HypreExtMultiABec::alpha  = " << alpha  << std::endl;
      std::cout << "  HypreExtMultiABec::alpha2 = " << alpha2 << std::endl;
      std::cout << "  HypreExtMultiABec::beta   = " << beta   << std::endl;
      std::cout << "  HypreExtMultiABec::gamma  = " << gamma  << std::endl;
      std::cout << "  HypreExtMultiABec::delta1 = " << delta1 << std::endl;
      std::cout << "  HypreExtMultiABec::delta2 = " << delta2 << std::endl;
    }

    // This block is here for debugging purposes only.  I think the code
    // below should be used in place of the HypreMultiABec version, not
    // in addition to it.
    //HypreMultiABec::loadLevelVectorB(level, rhs, inhom);
    //return;
  }

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
      f->copy<RunOn::Host>(rhs[mfi]);
    }
    Real* vec = f->dataPtr();

    // initialize rhs

    HYPRE_SStructVectorSetBoxValues(b, part, loV(reg), hiV(reg), 0, vec);
  }

  if (!inhom) {
    return;
  }

  // This version includes a new implementation of bcoefs separate from
  // the one in the base class.

  const Box& domain = bd[level]->getDomain();
  for (OrientationIter oitr; oitr; oitr++) {
    Orientation ori = oitr();
    int idir = ori.coordDir();
    IntVect vin = amrex::BASISV(idir), ves;
    vin = (ori.isLow() ? -vin : vin);    // outward normal unit vector
    ves = (ori.isLow() ?  ves : vin);    // edge shift vector
    Real h    = geom[level].CellSize(idir); // normal fine grid spacing
    Real h2   = 0.5 * h;
    Real th2  = 1.5 * h;
    Real ofac = (ori.isLow() ? -1.0 : 1.0);
    Real ofh  = ofac / h;
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      Box reg = amrex::adjCell(grids[level][i], ori);
      if (grids[level][i][ori] == domain[ori] || level == crse_level) {
        RadBoundCond     bct = bd[level]->bndryConds(ori)[i];
        // bct may be changed below if this is a mixed boundary
        const Real      &bcl = bd[level]->bndryLocs(ori)[i];
        const FArrayBox &fs  = bd[level]->bndryValues(ori)[mfi];
        const Mask      &msk = bd[level]->bndryMasks(ori,i);
//        const Box &bbox = (*bcoefs[level])[idir][i].box();
//        const Box &msb  = msk.box();

        // Treat an exposed grid edge here as a boundary condition
        // for the linear solver:

        reg.shift(-vin); // fine interior cells
        for (IntVect v = reg.smallEnd(); v <= reg.bigEnd(); reg.next(v)) {
          if (msk(v+vin) != RadBndryData::covered) {
            // determine what type of boundary this is and act accordingly
            if (reg[ori] == domain[ori] && bd[level]->mixedBndry(ori)) {
              bct = (*(bd[level]->bndryTypes(ori)[i]))(v+vin);
            }
            if (bct == AMREX_LO_DIRICHLET) {
              Real dfac, vfac;
              if (bho == 1) {
                dfac = 1.0 / ((bcl + h2) * (bcl + th2));
                vfac = h2 * th2 * dfac;
                dfac *= 2.0 * h;
              }
              else {
                dfac = 1.0 / (bcl + h2);
                vfac = h2 * dfac;
              }
              if (bcoefs[level]) {
                Real tmp = ((dfac / h) * fs(v+vin,bdcomp) *
                            beta * (*bcoefs[level])[idir][mfi](v+ves));
                HYPRE_SStructVectorAddToValues(b, part, getV1(v), 0, &tmp);
              }
              if (a2coefs[level]) {
                Real tmp = ((-0.5 * vfac) * fs(v+vin,bdcomp) *
                            alpha2 * (*a2coefs[level])[idir][mfi](v+ves));
                HYPRE_SStructVectorAddToValues(b, part, getV1(v), 0, &tmp);
              }
              if (ccoefs[level]) {
                // This was the centered difference version:
                //Real tmp = ((-vfac / h) * fs(v+vin,bdcomp) *
                //            gamma * ofac * (*ccoefs[level])[idir][mfi](v+ves));

                // For upwinding, vfac is not the correct interpolation factor.
                // If we're upwinding from the interior, ignore the bndry val.
                // If we're upwinding from the exterior, interpolate
                // (linearly) to the ghost cell center position.

                // The extra minus sign in both centered and upwind versions
                // is from moving this term to the rhs.

                Real fac = gamma * ofh;
                int  i0  = ori.isLow();
//                int  i1  = ori.isHigh();
                // upwinding from exterior is i1 direction:
                Real tmp = (-fac * (2.0 * h2) / (bcl + h2) * fs(v+vin,bdcomp)
                            * (*ccoefs[level])[idir][mfi](v+ves,i0));

                HYPRE_SStructVectorAddToValues(b, part, getV1(v), 0, &tmp);
              }
              if (d1coefs[level]) {
                Real tmp = ((-0.5 * dfac) * fs(v+vin,bdcomp) *
                            delta1 * ofac * (*d1coefs[level])[idir][mfi](v));
                HYPRE_SStructVectorAddToValues(b, part, getV1(v), 0, &tmp);
              }
              if (d2coefs[level]) {
                Real tmp = ((-0.5 * dfac) * fs(v+vin,bdcomp) *
                            delta2 * ofac * (*d2coefs[level])[idir][mfi](v+ves));
                HYPRE_SStructVectorAddToValues(b, part, getV1(v), 0, &tmp);
              }
            }
            else if (bct == AMREX_LO_NEUMANN) {
              // cmult should be c for photons, 1 for neutrinos
              Real cmult = 1.0;
              Real xi = 0.0; // xi should be passed in through RadBndry?
              Real tmp = beta / h;
              tmp = ((tmp - 0.5 * (delta1 + delta2) * ofac * xi) *
                     cmult * fs(v+vin,bdcomp));
              HYPRE_SStructVectorAddToValues(b, part, getV1(v), 0, &tmp);
            }
            else if (bct == AMREX_LO_MARSHAK || bct == AMREX_LO_SANCHEZ_POMRANING) {
              // cmult should be c for photons, 1 for neutrinos
              Real cmult = 1.0;
              Real xi = 0.0; // xi should be passed in through RadBndry?
              Real tmp = beta / h;
              tmp = 2.0 * ((tmp - 0.5 * (delta1 + delta2) * ofac * xi) *
                           cmult * fs(v+vin,bdcomp));
              HYPRE_SStructVectorAddToValues(b, part, getV1(v), 0, &tmp);
            }
          }
        }
      }
    }
  }
}


void HypreExtMultiABec::boundaryDterm(int level,
                                      MultiFab* Dterm,
                                      MultiFab& Soln,
                                      int icomp)
{
  BL_PROFILE("HypreExtMultiABec::boundaryDterm");

  const Box& domain = bd[level]->getDomain();
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(Soln); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = grids[level][i];
    for (OrientationIter oitr; oitr; oitr++) {
      int cdir(oitr());
      int idim = oitr().coordDir();
      const RadBoundCond &bct = bd[level]->bndryConds(oitr())[i];
      const Real      &bcl = bd[level]->bndryLocs(oitr())[i];
      const FArrayBox       &bcv  = bd[level]->bndryValues(oitr())[mfi];
      const Mask      &msk = bd[level]->bndryMasks(oitr(), i);

      if (reg[oitr()] == domain[oitr()]) {
        Array4<int const> tf_arr;
        int bctype = bct;
        if (bd[level]->mixedBndry(oitr())) {
          const BaseFab<int> &tf = *(bd[level]->bndryTypes(oitr())[i]);
          tf_arr = tf.array();
          bctype = -1;
        }
        HABEC::hdterm3(Dterm[idim][mfi].array(),
                       Soln[mfi].array(icomp),
                       reg,
                       cdir, bctype, tf_arr, bcl,
                       bcv.array(bdcomp),
                       msk.array(),
                       (*d2coefs[level])[idim][mfi].array(),
                       geom[level].CellSize());
      }
      else {
          HABEC::hdterm(Dterm[idim][mfi].array(),
                        Soln[mfi].array(icomp),
                        reg,
                        cdir, bct, bcl,
                        bcv.array(bdcomp),
                        msk.array(),
                        (*d2coefs[level])[idim][mfi].array(),
                        geom[level].CellSize());
      }
    }
  }
}
