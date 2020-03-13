#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#include <cmath>

#include <ppm.H>

using namespace amrex;

void
Castro::trace_plm(const Box& bx, const int idir,
                  Array4<Real const> const q_arr,
                  Array4<Real const> const qaux_arr,
                  Array4<Real const> const dq,
                  Array4<Real> const qm,
                  Array4<Real> const qp,
#if AMREX_SPACEDIM < 3
                  Array4<Real const> const dloga,
#endif
                  Array4<Real const> const srcQ,
                  const Box& vbx,
                  const Real dt) {

  // here, bx is the box we loop over -- this can include ghost cells
  // vbx is the valid box (no ghost cells)

  const auto dx = geom.CellSizeArray();

  const Real dtdx = dt/dx[idir];

  auto lo = bx.loVect3d();
  auto hi = bx.hiVect3d();

  auto vlo = vbx.loVect3d();
  auto vhi = vbx.hiVect3d();

#ifndef AMREX_USE_CUDA
  if (ppm_type != 0) {
    std::cout << "Oops -- shouldnt be in tracexy with ppm_type != 0" << std::endl;
    amrex::Error("Error:: trace_3d.f90 :: tracexy");
  }
#endif

  int QUN;
  int QUT;
  int QUTT;

  if (idir == 0) {
    QUN = QU;
    QUT = QV;
    QUTT = QW;

  } else if (idir == 1) {
    QUN = QV;
    QUT = QW;
    QUTT = QU;

  } else {
    QUN = QW;
    QUT = QU;
    QUTT = QV;
  }

  Real lsmall_dens = small_dens;
  Real lsmall_pres = small_pres;

  GpuArray<int, npassive> qpass_map_p;
  for (int n = 0; n < npassive; n++){
    qpass_map_p[n] = qpass_map[n];
  }

  // Compute left and right traced states

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    Real cc = qaux_arr(i,j,k,QC);
    Real csq = cc*cc;
    Real rho = q_arr(i,j,k,QRHO);

    Real un = q_arr(i,j,k,QUN);
    Real ut = q_arr(i,j,k,QUT);
    Real utt = q_arr(i,j,k,QUTT);

    Real p = q_arr(i,j,k,QPRES);
    Real rhoe = q_arr(i,j,k,QREINT);
    Real enth = (rhoe+p)/(rho*csq);

    Real drho = dq(i,j,k,QRHO);
    Real dun = dq(i,j,k,QUN);
    Real dut = dq(i,j,k,QUT);
    Real dutt = dq(i,j,k,QUTT);
    Real dp = dq(i,j,k,QPRES);
    Real drhoe = dq(i,j,k,QREINT);

    Real alpham = 0.5_rt*(dp/(rho*cc) - dun)*(rho/cc);
    Real alphap = 0.5_rt*(dp/(rho*cc) + dun)*(rho/cc);
    Real alpha0r = drho - dp/csq;
    Real alpha0e = drhoe - dp*enth;
    Real alpha0ut = dut;
    Real alpha0utt = dutt;

    Real e[3];
    e[0] = un - cc;
    e[1] = un;
    e[2] = un + cc;

    // construct the right state on the i interface

    Real ref_fac = 0.5_rt*(1.0_rt + dtdx*amrex::min(e[0], 0.0_rt));
    Real rho_ref = rho - ref_fac*drho;
    Real un_ref = un - ref_fac*dun;
    Real ut_ref = ut - ref_fac*dut;
    Real utt_ref = utt - ref_fac*dutt;
    Real p_ref = p - ref_fac*dp;
    Real rhoe_ref = rhoe - ref_fac*drhoe;

    // this is -(1/2) ( 1 + dt/dx lambda) (l . dq) r
    Real trace_fac0 = 0.0_rt; //  FOURTH*dtdx*(e(1) - e(1))*(1.0_rt - sign(1.0_rt, e(1)))
    Real trace_fac1 = 0.25_rt*dtdx*(e[0] - e[1])*(1.0_rt - std::copysign(1.0_rt, e[1]));
    Real trace_fac2 = 0.25_rt*dtdx*(e[0] - e[2])*(1.0_rt - std::copysign(1.0_rt, e[2]));

    Real apright = trace_fac2*alphap;
    Real amright = trace_fac0*alpham;

    Real azrright = trace_fac1*alpha0r;
    Real azeright = trace_fac1*alpha0e;
    Real azut1rght = trace_fac1*alpha0ut;
    Real azutt1rght = trace_fac1*alpha0utt;

    if ((idir == 0 && i >= vlo[0]) ||
        (idir == 1 && j >= vlo[1]) ||
        (idir == 2 && k >= vlo[2])) {

      qp(i,j,k,QRHO) = amrex::max(lsmall_dens, rho_ref + apright + amright + azrright);
      qp(i,j,k,QUN) = un_ref + (apright - amright)*cc/rho;
      qp(i,j,k,QUT) = ut_ref + azut1rght;
      qp(i,j,k,QUTT) = utt_ref + azutt1rght;
      qp(i,j,k,QPRES) = amrex::max(lsmall_pres, p_ref + (apright + amright)*csq);
      qp(i,j,k,QREINT) = rhoe_ref + (apright + amright)*enth*csq + azeright;

      // add the source terms 
      qp(i,j,k,QRHO  ) += 0.5_rt*dt*srcQ(i,j,k,QRHO);
      qp(i,j,k,QRHO  ) = std::max(lsmall_dens, qp(i,j,k,QRHO));
      qp(i,j,k,QUN   ) += 0.5_rt*dt*srcQ(i,j,k,QUN);
      qp(i,j,k,QUT   ) += 0.5_rt*dt*srcQ(i,j,k,QUT);
      qp(i,j,k,QUTT  ) += 0.5_rt*dt*srcQ(i,j,k,QUTT);
      qp(i,j,k,QREINT) += 0.5_rt*dt*srcQ(i,j,k,QREINT);
      qp(i,j,k,QPRES ) += 0.5_rt*dt*srcQ(i,j,k,QPRES);

    }

    // now construct the left state on the i+1 interface

    ref_fac = 0.5_rt*(1.0_rt - dtdx*amrex::max(e[2], 0.0_rt));
    rho_ref = rho + ref_fac*drho;
    un_ref = un + ref_fac*dun;
    ut_ref = ut + ref_fac*dut;
    utt_ref = utt + ref_fac*dutt;
    p_ref = p + ref_fac*dp;
    rhoe_ref = rhoe + ref_fac*drhoe;

    trace_fac0 = 0.25_rt*dtdx*(e[2] - e[0])*(1.0_rt + std::copysign(1.0_rt, e[0]));
    trace_fac1 = 0.25_rt*dtdx*(e[2] - e[1])*(1.0_rt + std::copysign(1.0_rt, e[1]));
    trace_fac2 = 0.0_rt; //  FOURTH*dtdx*(e(3) - e(3))*(1.0_rt + sign(1.0_rt, e(3)))

    Real apleft = trace_fac2*alphap;
    Real amleft = trace_fac0*alpham;

    Real azrleft = trace_fac1*alpha0r;
    Real azeleft = trace_fac1*alpha0e;
    Real azut1left = trace_fac1*alpha0ut;
    Real azutt1left = trace_fac1*alpha0utt;

    if (idir == 0 && i <= vhi[0]) {
      qm(i+1,j,k,QRHO) = amrex::max(lsmall_dens, rho_ref + apleft + amleft + azrleft);
      qm(i+1,j,k,QUN) = un_ref + (apleft - amleft)*cc/rho;
      qm(i+1,j,k,QUT) = ut_ref + azut1left;
      qm(i+1,j,k,QUTT) = utt_ref + azutt1left;
      qm(i+1,j,k,QPRES) = amrex::max(lsmall_pres, p_ref + (apleft + amleft)*csq);
      qm(i+1,j,k,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft;

      // add the source terms
      qm(i+1,j,k,QRHO) = amrex::max(lsmall_dens, qm(i+1,j,k,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO));
      qm(i+1,j,k,QUN) += 0.5_rt*dt*srcQ(i,j,k,QUN);
      qm(i+1,j,k,QUT) += 0.5_rt*dt*srcQ(i,j,k,QUT);
      qm(i+1,j,k,QUTT) += 0.5_rt*dt*srcQ(i,j,k,QUTT);
      qm(i+1,j,k,QREINT) += 0.5_rt*dt*srcQ(i,j,k,QREINT);
      qm(i+1,j,k,QPRES) += 0.5_rt*dt*srcQ(i,j,k,QPRES);

    } else if (idir == 1 && j <= vhi[1]) {
      qm(i,j+1,k,QRHO) = amrex::max(lsmall_dens, rho_ref + apleft + amleft + azrleft);
      qm(i,j+1,k,QUN) = un_ref + (apleft - amleft)*cc/rho;
      qm(i,j+1,k,QUT) = ut_ref + azut1left;
      qm(i,j+1,k,QUTT) = utt_ref + azutt1left;
      qm(i,j+1,k,QPRES) = max(lsmall_pres, p_ref + (apleft + amleft)*csq);
      qm(i,j+1,k,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft;

      // add the source terms
      qm(i,j+1,k,QRHO) = amrex::max(lsmall_dens, qm(i,j+1,k,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO));
      qm(i,j+1,k,QUN) += 0.5_rt*dt*srcQ(i,j,k,QUN);
      qm(i,j+1,k,QUT) += 0.5_rt*dt*srcQ(i,j,k,QUT);
      qm(i,j+1,k,QUTT) += 0.5_rt*dt*srcQ(i,j,k,QUTT);
      qm(i,j+1,k,QREINT) += 0.5_rt*dt*srcQ(i,j,k,QREINT);
      qm(i,j+1,k,QPRES) += 0.5_rt*dt*srcQ(i,j,k,QPRES);

    } else if (idir == 2 && k <= vhi[2]) {
      qm(i,j,k+1,QRHO) = amrex::max(lsmall_dens, rho_ref + apleft + amleft + azrleft);
      qm(i,j,k+1,QUN) = un_ref + (apleft - amleft)*cc/rho;
      qm(i,j,k+1,QUT) = ut_ref + azut1left;
      qm(i,j,k+1,QUTT) = utt_ref + azutt1left;
      qm(i,j,k+1,QPRES) = amrex::max(lsmall_pres, p_ref + (apleft + amleft)*csq);
      qm(i,j,k+1,QREINT) = rhoe_ref + (apleft + amleft)*enth*csq + azeleft;

      // add the source terms
      qm(i,j,k+1,QRHO) = amrex::max(lsmall_dens, qm(i,j,k+1,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO));
      qm(i,j,k+1,QUN) += 0.5_rt*dt*srcQ(i,j,k,QUN);
      qm(i,j,k+1,QUT) += 0.5_rt*dt*srcQ(i,j,k,QUT);
      qm(i,j,k+1,QUTT) += 0.5_rt*dt*srcQ(i,j,k,QUTT);
      qm(i,j,k+1,QREINT) += 0.5_rt*dt*srcQ(i,j,k,QREINT);
      qm(i,j,k+1,QPRES) += 0.5_rt*dt*srcQ(i,j,k,QPRES);

    }

#if (AMREX_SPACEDIM < 3)
    // geometry source terms -- these only apply to the x-states
    if (idir == 0 && dloga(i,j,k) != 0.0_rt) {
      Real courn = dtdx*(cc + abs(un));
      Real eta = (1.0_rt-courn)/(cc*dt*abs(dloga(i,j,k)));
      Real dlogatmp = amrex::min(eta, 1.0_rt)*dloga(i,j,k);
      Real sourcr = -0.5_rt*dt*rho*dlogatmp*un;
      Real sourcp = sourcr*csq;
      Real source = sourcp*enth;

      if (i <= vhi[0]) {
        qm(i+1,j,k,QRHO) += sourcr;
        qm(i+1,j,k,QRHO) = amrex::max(qm(i+1,j,k,QRHO), lsmall_dens);
        qm(i+1,j,k,QPRES) += sourcp;
        qm(i+1,j,k,QREINT) += source;
      }

      if (i >= vlo[0]) {
        qp(i,j,k,QRHO) += sourcr;
        qp(i,j,k,QRHO) = amrex::max(qp(i,j,k,QRHO), lsmall_dens);
        qp(i,j,k,QPRES) += sourcp;
        qp(i,j,k,QREINT) += source;
      }
    }
#endif

    for (int ipassive = 0; ipassive < npassive; ipassive++) {
      int n = qpass_map_p[ipassive];

      // Right state
      if ((idir == 0 && i >= vlo[0]) ||
          (idir == 1 && j >= vlo[1]) ||
          (idir == 2 && k >= vlo[2])) {

        Real un = q_arr(i,j,k,QUN);
        Real spzero = un >= 0.0_rt ? -1.0_rt : un*dtdx;
        qp(i,j,k,n) = q_arr(i,j,k,n) + 0.5_rt*(-1.0_rt - spzero)*dq(i,j,k,n);
#if  PRIM_SPECIES_HAVE_SOURCES
        qp(i,j,k,n) += 0.5_rt*dt*srcQ(i,j,k,n);
#endif

        // Left state
        un = q_arr(i,j,k,QUN);
        spzero = un >= 0.0_rt ? un*dtdx : 1.0_rt;
        Real acmpleft = 0.5_rt*(1.0_rt - spzero )*dq(i,j,k,n);

        if (idir == 0 && i <= vhi[0]) {
          qm(i+1,j,k,n) = q_arr(i,j,k,n) + acmpleft;
#if  PRIM_SPECIES_HAVE_SOURCES
          qm(i+1,j,k,n) += 0.5_rt*dt*srcQ(i,j,k,n);
#endif

        } else if (idir == 1 && j <= vhi[1]) {
          qm(i,j+1,k,n) = q_arr(i,j,k,n) + acmpleft;
#if  PRIM_SPECIES_HAVE_SOURCES
          qm(i,j+1,k,n) += 0.5_rt*dt*srcQ(i,j,k,n);
#endif

        } else if (idir == 2 && k <= vhi[2]) {
          qm(i,j,k+1,n) = q_arr(i,j,k,n) + acmpleft;
#if  PRIM_SPECIES_HAVE_SOURCES
          qm(i,j,k+1,n) += 0.5_rt*dt*srcQ(i,j,k,n);
#endif
        }
      }
    }


  });
}
