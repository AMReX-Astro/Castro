#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;


void
Castro::src_to_prim(const Box& bx,
                    Array4<Real> const q,
                    Array4<Real> const qaux,
                    Array4<Real> const src,
                    Array4<Real> const srcQ)
{

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {


      for (int n = 0; n < NQSRC; ++n) {
        srcQ(i,j,k,n) = 0.0;
      }

      Real rhoinv = 1.0 / q(i,j,k,QRHO);

      srcQ(i,j,k,QRHO) = src(i,j,k,Density);
      srcQ(i,j,k,QU) = (src(i,j,k,Xmom) - q(i,j,k,QU) * srcQ(i,j,k,QRHO)) * rhoinv;
      srcQ(i,j,k,QV) = (src(i,j,k,Ymom) - q(i,j,k,QV) * srcQ(i,j,k,QRHO)) * rhoinv;
      srcQ(i,j,k,QW) = (src(i,j,k,Zmom) - q(i,j,k,QW) * srcQ(i,j,k,QRHO)) * rhoinv;
      srcQ(i,j,k,QREINT) = src(i,j,k,Eint);
      srcQ(i,j,k,QPRES ) = qaux(i,j,k,QDPDE) *
        (srcQ(i,j,k,QREINT) - q(i,j,k,QREINT)*srcQ(i,j,k,QRHO)*rhoinv) *
        rhoinv + qaux(i,j,k,QDPDR)*srcQ(i,j,k,QRHO);

#ifdef PRIM_SPECIES_HAVE_SOURCES
      for (int ipassive = 0; ipassive < npassive; ++ipassive) {
        int n = upass_map(ipassive);
        int iq = qpass_map(ipassive);

       // we may not be including the ability to have species sources,
       //  so check to make sure that we are < NQSRC
        srcQ(i,j,k,iq) = (src(i,j,k,n) - q(i,j,k,iq) * srcQ(i,j,k,QRHO) ) /
          q(i,j,k,QRHO);
      }
#endif
  });

}
