#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#ifdef GRAVITY
#include "Gravity.H"
#endif

#include <ppm.H>

using namespace amrex;

void
Castro::mol_plm_reconstruct(const Box& bx,
                            const int idir,
                            Array4<Real const> const q_arr,
                            Array4<Real const> const flatn_arr,
                            Array4<Real> const dq,
                            Array4<Real> const qm,
                            Array4<Real> const qp) {

  const auto dx = geom.CellSizeArray();

#ifdef GRAVITY
  Real lconst_grav =  gravity->get_const_grav();
#else
  Real lconst_grav = 0.0_rt;
#endif

  int lplm_well_balanced = plm_well_balanced;

  for (int n = 0; n < NQ; n++) {
    // piecewise linear slopes
    uslope(bx, idir,
           q_arr, n,
           flatn_arr, dq);
  }


  AMREX_PARALLEL_FOR_4D(bx, NQ, i, j, k, n,
  {


   // this is a loop over zones.  For each slope in the zone, fill the
   // two adjacent edge states (e.g., the right state at i-1/2 and the
   // left state at i+1/2

   if (idir == 0) {

     if (lplm_well_balanced == 1 && n == QPRES && idir == AMREX_SPACEDIM-1) {

       // left state at i+1/2 interface
       qm(i+1,j,k,n) = q_arr(i,j,k,n) + 0.5_rt*dq(i,j,k,n) +
         0.5_rt*dx[0]*q_arr(i,j,k,QRHO)*lconst_grav;

       // right state at i-1/2 interface
       qp(i,j,k,n) = q_arr(i,j,k,n) - 0.5_rt*dq(i,j,k,n) -
         0.5_rt*dx[0]*q_arr(i,j,k,QRHO)*lconst_grav;

     } else {
       // left state at i+1/2 interface
       qm(i+1,j,k,n) = q_arr(i,j,k,n) + 0.5_rt*dq(i,j,k,n);

       // right state at i-1/2 interface
       qp(i,j,k,n) = q_arr(i,j,k,n) - 0.5_rt*dq(i,j,k,n);

     }

#if AMREX_SPACEDIM >= 2
   } else if (idir == 1) {

     if (lplm_well_balanced == 1 && n == QPRES && idir == AMREX_SPACEDIM-1) {

       // left state at i+1/2 interface
       qm(i,j+1,k,n) = q_arr(i,j,k,n) + 0.5_rt*dq(i,j,k,n) +
         0.5_rt*dx[1]*q_arr(i,j,k,QRHO)*lconst_grav;

       // right state at i-1/2 interface
       qp(i,j,k,n) = q_arr(i,j,k,n) - 0.5_rt*dq(i,j,k,n) -
         0.5_rt*dx[1]*q_arr(i,j,k,QRHO)*lconst_grav;

     } else {

       // left state at j+1/2 interface
       qm(i,j+1,k,n) = q_arr(i,j,k,n) + 0.5_rt*dq(i,j,k,n);

       // right state at j-1/2 interface
       qp(i,j,k,n) = q_arr(i,j,k,n) - 0.5_rt*dq(i,j,k,n);

     }
#endif

#if AMREX_SPACEDIM == 3
   } else {

     if (lplm_well_balanced == 1 && n == QPRES && idir == AMREX_SPACEDIM-1) {

       // left state at i+1/2 interface
       qm(i,j,k+1,n) = q_arr(i,j,k,n) + 0.5_rt*dq(i,j,k,n) +
         0.5_rt*dx[2]*q_arr(i,j,k,QRHO)*lconst_grav;

       // right state at i-1/2 interface
       qp(i,j,k,n) = q_arr(i,j,k,n) - 0.5_rt*dq(i,j,k,n) -
         0.5_rt*dx[2]*q_arr(i,j,k,QRHO)*lconst_grav;

     } else {

       // left state at k+1/2 interface
       qm(i,j,k+1,n) = q_arr(i,j,k,n) + 0.5_rt*dq(i,j,k,n);

       // right state at k-1/2 interface
       qp(i,j,k,n) = q_arr(i,j,k,n) - 0.5_rt*dq(i,j,k,n);

     }
#endif
   }

  });


  // special care for reflecting BCs
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

  bool lo_bc_test = lo_bc[idir] == Symmetry;
  bool hi_bc_test = hi_bc[idir] == Symmetry;

  // we have to do this after the loops above, since here we will
  // consider interfaces, not zones

  if (idir == 0) {
    if (lo_bc_test) {

      AMREX_PARALLEL_FOR_3D(bx, i, j, k,
      {

       // reset the left state at domlo(0) if needed -- it is outside the domain
       if (i == domlo[0]) {
         for (int n = 0; n < NQ; n++) {
           if (n == QU) {
             qm(i,j,k,QU) = -qp(i,j,k,QU);
           } else {
             qm(i,j,k,n) = qp(i,j,k,n);
           }
         }
       }
      });

    }


    if (hi_bc_test) {

      AMREX_PARALLEL_FOR_3D(bx, i, j, k,
      {

       // reset the right state at domhi(0)+1 if needed -- it is outside the domain
       if (i == domhi[0]+1) {
         for (int n = 0; n < NQ; n++) {
           if (n == QU) {
             qp(i,j,k,QU) = -qm(i,j,k,QU);
           } else {
             qp(i,j,k,n) = qm(i,j,k,n);
           }
         }
       }
      });

    }

#if AMREX_SPACEDIM >= 2
  } else if (idir == 1) {

    if (lo_bc_test) {

      AMREX_PARALLEL_FOR_3D(bx, i, j, k,
      {

       // reset the left state at domlo(0) if needed -- it is outside the domain
       if (j == domlo[1]) {
         for (int n = 0; n < NQ; n++) {
           if (n == QV) {
             qm(i,j,k,QV) = -qp(i,j,k,QV);
           } else {
             qm(i,j,k,n) = -qp(i,j,k,n);
           }
         }
       }
      });

    }


    if (hi_bc_test) {

      AMREX_PARALLEL_FOR_3D(bx, i, j, k,
      {

       // reset the right state at domhi(0)+1 if needed -- it is outside the domain
       if (j == domhi[1]+1) {
         for (int n = 0; n < NQ; n++) {
           if (n == QV) {
             qp(i,j,k,QV) = -qm(i,j,k,QV);
           } else {
             qp(i,j,k,n) = qm(i,j,k,n);
           }
         }
       }
      });

    }

#endif
#if AMREX_SPACEDIM == 3
  } else {

    if (lo_bc_test) {

      AMREX_PARALLEL_FOR_3D(bx, i, j, k,
      {

       // reset the left state at domlo(0) if needed -- it is outside the domain
       if (k == domlo[2]) {
         for (int n = 0; n < NQ; n++) {
           if (n == QW) {
             qm(i,j,k,QW) = -qp(i,j,k,QW);
           } else {
             qm(i,j,k,n) = qp(i,j,k,n);
           }
         }
       }
      });

    }


    if (hi_bc_test) {

      AMREX_PARALLEL_FOR_3D(bx, i, j, k,
      {

       // reset the right state at domhi(0)+1 if needed -- it is outside the domain
       if (k == domhi[2]+1) {
         for (int n = 0; n < NQ; n++) {
           if (n == QW) {
             qp(i,j,k,QW) = -qm(i,j,k,QW);
           } else {
             qp(i,j,k,n) = qm(i,j,k,n);
           }
         }
       }
      });

    }

#endif

  }
}


void
Castro::mol_ppm_reconstruct(const Box& bx,
                            const int idir,
                            Array4<Real const> const q_arr,
                            Array4<Real const> const flatn_arr,
                            Array4<Real> const qm,
                            Array4<Real> const qp) {

  // special care for reflecting BCs
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

  bool lo_bc_test = lo_bc[idir] == Symmetry;
  bool hi_bc_test = hi_bc[idir] == Symmetry;


  AMREX_PARALLEL_FOR_4D(bx, NQ, i, j, k, n,
  {

    Real s[7];
    Real flat = flatn_arr(i,j,k);
    Real sm;
    Real sp;

    if (idir == 0) {
      s[im3] = q_arr(i-3,j,k,n);
      s[im2] = q_arr(i-2,j,k,n);
      s[im1] = q_arr(i-1,j,k,n);
      s[i0]  = q_arr(i,j,k,n);
      s[ip1] = q_arr(i+1,j,k,n);
      s[ip2] = q_arr(i+2,j,k,n);
      s[ip3] = q_arr(i+3,j,k,n);

    } else if (idir == 1) {
      s[im3] = q_arr(i,j-3,k,n);
      s[im2] = q_arr(i,j-2,k,n);
      s[im1] = q_arr(i,j-1,k,n);
      s[i0]  = q_arr(i,j,k,n);
      s[ip1] = q_arr(i,j+1,k,n);
      s[ip2] = q_arr(i,j+2,k,n);
      s[ip3] = q_arr(i,j+3,k,n);

    } else {
      s[im3] = q_arr(i,j,k-3,n);
      s[im2] = q_arr(i,j,k-2,n);
      s[im1] = q_arr(i,j,k-1,n);
      s[i0]  = q_arr(i,j,k,n);
      s[ip1] = q_arr(i,j,k+1,n);
      s[ip2] = q_arr(i,j,k+2,n);
      s[ip3] = q_arr(i,j,k+3,n);

    }

    ppm_reconstruct(s, i, j, k, idir,
                    lo_bc_test, hi_bc_test, domlo, domhi,
                    flat, sm, sp);

    if (idir == 0) {
      // right state at i-1/2
      qp(i,j,k,n) = sm;

      // left state at i+1/2
      qm(i+1,j,k,n) = sp;

    } else if (idir == 1) {
      // right state at j-1/2
      qp(i,j,k,n) = sm;

      // left state at j+1/2
      qm(i,j+1,k,n) = sp;

    } else {
      // right state at k-1/2
      qp(i,j,k,n) = sm;

      // left state at k+1/2
      qm(i,j,k+1,n) = sp;

    }

  });

}

