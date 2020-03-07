#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#include "eos.H"

using namespace amrex;


void
Castro::ctoprim(const Box& bx,
                Array4<Real const> const uin,
#ifdef RADIATION
                Array4<Real const> const Erin,
                Array4<Real const> const lam,
#endif
                Array4<Real> const q_arr,
                Array4<Real> const qaux_arr) {


  GpuArray<int, npassive> upass_map_p;
  GpuArray<int, npassive> qpass_map_p;
  for (int n = 0; n < npassive; ++n) {
    upass_map_p[n] = upass_map[n];
    qpass_map_p[n] = qpass_map[n];
  }

  Real lsmall_dens = small_dens;
  Real ldual_energy_eta1 = dual_energy_eta1;

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {


#ifndef AMREX_USE_CUDA
    if (uin(i,j,k,URHO) <= 0.0_rt) {
      std::cout << std::endl;
      std::cout << ">>> Error: advection_util_nd.F90::ctoprim " << i << " " << j << " " << k << std::endl;
      std::cout << ">>> ... negative density " << uin(i,j,k,URHO) << std::endl;
      amrex::Error("Error:: advection_util_nd.f90 :: ctoprim");
    } else if (uin(i,j,k,URHO) < lsmall_dens) {
      std::cout << std::endl;
      std::cout << ">>> Error: advection_util_nd.F90::ctoprim " << i << " " << j << " " << k << std::endl;
      std::cout << ">>> ... small density " << uin(i,j,k,URHO) << std::endl;
      amrex::Error("Error:: advection_util_nd.f90 :: ctoprim");
    }
#endif

    q_arr(i,j,k,QRHO) = uin(i,j,k,URHO);
    Real rhoinv = 1.0_rt/q_arr(i,j,k,QRHO);

    q_arr(i,j,k,QU) = uin(i,j,k,UMX) * rhoinv;
    q_arr(i,j,k,QV) = uin(i,j,k,UMY) * rhoinv;
    q_arr(i,j,k,QW) = uin(i,j,k,UMZ) * rhoinv;

    // Get the internal energy, which we'll use for
    // determining the pressure.  We use a dual energy
    // formalism. If (E - K) < eta1 and eta1 is suitably
    // small, then we risk serious numerical truncation error
    // in the internal energy.  Therefore we'll use the result
    // of the separately updated internal energy equation.
    // Otherwise, we'll set e = E - K.

    Real kineng = 0.5_rt * q_arr(i,j,k,QRHO) * (q_arr(i,j,k,QU)*q_arr(i,j,k,QU) +
                                                q_arr(i,j,k,QV)*q_arr(i,j,k,QV) +
                                                q_arr(i,j,k,QW)*q_arr(i,j,k,QW));

    if ((uin(i,j,k,UEDEN) - kineng) > ldual_energy_eta1*uin(i,j,k,UEDEN)) {
      q_arr(i,j,k,QREINT) = (uin(i,j,k,UEDEN) - kineng) * rhoinv;
    } else {
      q_arr(i,j,k,QREINT) = uin(i,j,k,UEINT) * rhoinv;
    }

    // If we're advecting in the rotating reference frame,
    // then subtract off the rotation component here.

#ifdef ROTATION
    if (do_rotation == 1 && state_in_rotating_frame != 1) {
      Real vel[3];
      for (int n = 0; n < 3; n++) {
        vel[n] = uin(i,j,k,UMX+n) * rhoinv;
      }

      call inertial_to_rotational_velocity([i, j, k], amr_time, vel);
      q_arr(i,j,k,QU) = vel[0];
      q_arr(i,j,k,QV) = vel[1];
      q_arr(i,j,k,QW) = vel[2];
    }
#endif

    q_arr(i,j,k,QTEMP) = uin(i,j,k,UTEMP);
#ifdef RADIATION
    for (int g = 0; g < NGROUPS; g++) {
      q_arr(i,j,k,QRAD+g) = Erin(i,j,k,g);
    }
#endif

    // Load passively advected quatities into q
    for (int ipassive = 0; ipassive < npassive; ipassive++) {
      int n  = upass_map_p[ipassive];
      int iq = qpass_map_p[ipassive];
      q_arr(i,j,k,iq) = uin(i,j,k,n) * rhoinv;
    }

    // get gamc, p, T, c, csml using q state
    eos_t eos_state;
    eos_state.T = q_arr(i,j,k,QTEMP);
    eos_state.rho = q_arr(i,j,k,QRHO);
    eos_state.e = q_arr(i,j,k,QREINT);
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n]  = q_arr(i,j,k,QFS+n);
    }
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = q_arr(i,j,k,QFX+n);
    }

    eos(eos_input_re, eos_state);

    q_arr(i,j,k,QTEMP) = eos_state.T;
    q_arr(i,j,k,QREINT) = eos_state.e * q_arr(i,j,k,QRHO);
    q_arr(i,j,k,QPRES) = eos_state.p;
    q_arr(i,j,k,QGC) = eos_state.gam1;

    qaux_arr(i,j,k,QDPDR) = eos_state.dpdr_e;
    qaux_arr(i,j,k,QDPDE) = eos_state.dpde;

#ifdef RADIATION
    qaux_arr(i,j,k,QGAMCG) = eos_state.gam1;
    qaux_arr(i,j,k,QCG) = eos_state.cs;

    Real lams[NGROUPS];
    for (int g = 0; g < NGROUPS; g++) {
      lams[g] = lam(i,j,k,g);
    }
    Real qs[NQ];
    for (int n = 0; g < NQ; g++) {
      qs[n] = q_arr(i,j,k,n);
    }
    Real ptot;
    Real ctot;
    Real gamc_tot;
    compute_ptot_ctot(lams, qs, qaux_arr(i,j,k,QCG), ptot, ctot, gamc_tot);

    q_arr(i,j,k,QPTOT) = ptot;

    qaux_arr(i,j,k,QC) = ctot;
    qaux_arr(i,j,k,QGAMC) = gamc_tot;

    q_arr(i,j,k,qreitot) = q_arr(i,j,k,QREINT);
    for (int g = 0; g < NGROUPS; g++) {
      qaux_arr(i,j,k,QLAMS+g) = lam(i,j,k,g);
      q_arr(i,j,k,qreitot) += q_arr(i,j,k,QRAD+g);
    }

#else
    qaux_arr(i,j,k,QGAMC) = eos_state.gam1;
    qaux_arr(i,j,k,QC) = eos_state.cs;
#endif

  });
}


void
Castro::src_to_prim(const Box& bx,
                    Array4<Real const> const q_arr,
                    Array4<Real const> const qaux_arr,
                    Array4<Real const> const src,
                    Array4<Real> const srcQ)
{

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {


      for (int n = 0; n < NQSRC; ++n) {
        srcQ(i,j,k,n) = 0.0_rt;
      }

      Real rhoinv = 1.0_rt / q_arr(i,j,k,QRHO);

      srcQ(i,j,k,QRHO) = src(i,j,k,URHO);
      srcQ(i,j,k,QU) = (src(i,j,k,UMX) - q_arr(i,j,k,QU) * srcQ(i,j,k,QRHO)) * rhoinv;
      srcQ(i,j,k,QV) = (src(i,j,k,UMY) - q_arr(i,j,k,QV) * srcQ(i,j,k,QRHO)) * rhoinv;
      srcQ(i,j,k,QW) = (src(i,j,k,UMZ) - q_arr(i,j,k,QW) * srcQ(i,j,k,QRHO)) * rhoinv;
      srcQ(i,j,k,QREINT) = src(i,j,k,UEINT);
      srcQ(i,j,k,QPRES ) = qaux_arr(i,j,k,QDPDE) *
        (srcQ(i,j,k,QREINT) - q_arr(i,j,k,QREINT)*srcQ(i,j,k,QRHO)*rhoinv) *
        rhoinv + qaux_arr(i,j,k,QDPDR)*srcQ(i,j,k,QRHO);

#ifdef PRIM_SPECIES_HAVE_SOURCES
      for (int ipassive = 0; ipassive < npassive; ++ipassive) {
        int n = upass_map[ipassive];
        int iq = qpass_map[ipassive];

       // we may not be including the ability to have species sources,
       //  so check to make sure that we are < NQSRC
        srcQ(i,j,k,iq) = (src(i,j,k,n) - q_arr(i,j,k,iq) * srcQ(i,j,k,QRHO) ) /
          q_arr(i,j,k,QRHO);
      }
#endif
  });

}


void
Castro::shock(const Box& bx,
              Array4<Real const> const q_arr,
              Array4<Real> const shk) {

  // This is a basic multi-dimensional shock detection algorithm.
  // This implementation follows Flash, which in turn follows
  // AMRA and a Woodward (1995) (supposedly -- couldn't locate that).
  //
  // The spirit of this follows the shock detection in Colella &
  // Woodward (1984)
  //

  constexpr Real small = 1.e-10_rt;
  constexpr Real eps = 0.33e0_rt;

  const auto dx = geom.CellSizeArray();
  const int coord_type = geom.Coord();

  Real dxinv = 1.0_rt / dx[0];
#if AMREX_SPACEDIM >= 2
  Real dyinv = 1.0_rt / dx[1];
#endif
#if AMREX_SPACEDIM == 3
  Real dzinv = 1.0_rt / dx[2];
#endif

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {
    Real div_u = 0.0_rt;

    // construct div{U}
    if (coord_type == 0) {

      // Cartesian
      div_u += 0.5_rt * (q_arr(i+1,j,k,QU) - q_arr(i-1,j,k,QU)) * dxinv;
#if (AMREX_SPACEDIM >= 2)
      div_u += 0.5_rt * (q_arr(i,j+1,k,QV) - q_arr(i,j-1,k,QV)) * dyinv;
#endif
#if (AMREX_SPACEDIM == 3)
      div_u += 0.5_rt * (q_arr(i,j,k+1,QW) - q_arr(i,j,k-1,QW)) * dzinv;
#endif

   } else if (coord_type == 1) {

     // r-z
     Real rc = (i + 0.5_rt) * dx[0];
     Real rm = (i - 1 + 0.5_rt) * dx[0];
     Real rp = (i + 1 + 0.5_rt) * dx[0];

#if (AMREX_SPACEDIM == 1)
     div_u += 0.5_rt * (rp * q_arr(i+1,j,k,QU) - rm * q_arr(i-1,j,k,QU)) / (rc * dx[0]);
#endif
#if (AMREX_SPACEDIM == 2)
     div_u += 0.5_rt * (rp * q_arr(i+1,j,k,QU) - rm * q_arr(i-1,j,k,QU)) / (rc * dx[0]) +
              0.5_rt * (q_arr(i,j+1,k,QV) - q_arr(i,j-1,k,QV)) * dyinv;
#endif

    } else if (coord_type == 2) {

      // 1-d spherical
      Real rc = (i + 0.5_rt) * dx[0];
      Real rm = (i - 1 + 0.5_rt) * dx[0];
      Real rp = (i + 1 + 0.5_rt) * dx[0];

      div_u += 0.5_rt * (rp * rp * q_arr(i+1,j,k,QU) - rm * rm * q_arr(i-1,j,k,QU)) / (rc * rc * dx[0]);

#ifndef AMREX_USE_CUDA

    } else {
      amrex::Error("ERROR: invalid coord_type in shock");
#endif
    }

    // find the pre- and post-shock pressures in each direction
    Real px_pre;
    Real px_post;
    Real e_x;

    if (q_arr(i+1,j,k,QPRES) - q_arr(i-1,j,k,QPRES) < 0.0_rt) {
      px_pre = q_arr(i+1,j,k,QPRES);
      px_post = q_arr(i-1,j,k,QPRES);
    } else {
      px_pre = q_arr(i-1,j,k,QPRES);
      px_post = q_arr(i+1,j,k,QPRES);
    }

    // use compression to create unit vectors for the shock direction
    e_x = std::pow(q_arr(i+1,j,k,QU) - q_arr(i-1,j,k,QU), 2);

    Real py_pre;
    Real py_post;
    Real e_y;

#if (AMREX_SPACEDIM >= 2)
    if (q_arr(i,j+1,k,QPRES) - q_arr(i,j-1,k,QPRES) < 0.0_rt) {
      py_pre = q_arr(i,j+1,k,QPRES);
      py_post = q_arr(i,j-1,k,QPRES);
    } else {
      py_pre = q_arr(i,j-1,k,QPRES);
      py_post = q_arr(i,j+1,k,QPRES);
    }

    e_y = std::pow(q_arr(i,j+1,k,QV) - q_arr(i,j-1,k,QV), 2);

#else
    py_pre = 0.0_rt;
    py_post = 0.0_rt;

    e_y = 0.0_rt;
#endif

    Real pz_pre;
    Real pz_post;
    Real e_z;

#if (AMREX_SPACEDIM == 3)
    if (q_arr(i,j,k+1,QPRES) - q_arr(i,j,k-1,QPRES) < 0.0_rt) {
      pz_pre  = q_arr(i,j,k+1,QPRES);
      pz_post = q_arr(i,j,k-1,QPRES);
    } else {
      pz_pre  = q_arr(i,j,k-1,QPRES);
      pz_post = q_arr(i,j,k+1,QPRES);
    }

    e_z = std::pow(q_arr(i,j,k+1,QW) - q_arr(i,j,k-1,QW), 2);

#else
    pz_pre = 0.0_rt;
    pz_post = 0.0_rt;

    e_z = 0.0_rt;
#endif

    Real denom = 1.0_rt / (e_x + e_y + e_z + small);

    e_x = e_x * denom;
    e_y = e_y * denom;
    e_z = e_z * denom;

    // project the pressures onto the shock direction
    Real p_pre  = e_x * px_pre + e_y * py_pre + e_z * pz_pre;
    Real p_post = e_x * px_post + e_y * py_post + e_z * pz_post;

    // test for compression + pressure jump to flag a shock
    // this avoid U = 0, so e_x, ... = 0
    Real pjump = p_pre == 0 ? 0.0_rt : eps - (p_post - p_pre) / p_pre;

    if (pjump < 0.0 && div_u < 0.0_rt) {
      shk(i,j,k) = 1.0_rt;
    } else {
      shk(i,j,k) = 0.0_rt;
    }
  });

}


void
Castro::divu(const Box& bx,
             Array4<Real const> const q_arr,
             Array4<Real> const div) {
  // this computes the *node-centered* divergence

  const auto dx = geom.CellSizeArray();
  const int coord_type = geom.Coord();

  const auto problo = geom.ProbLoArray();

  Real dxinv = 1.0_rt / dx[0];
#if AMREX_SPACEDIM >= 2
  Real dyinv = 1.0_rt / dx[1];
#else
  Real dyinv = 0.0_rt;
#endif
#if AMREX_SPACEDIM == 3
  Real dzinv = 1.0_rt / dx[2];
#else
  Real dzinv = 0.0_rt;
#endif

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

#if AMREX_SPACEDIM == 1
    if (coord_type == 0) {
      div(i,j,k) = (q_arr(i,j,k,QU) - q_arr(i-1,j,k,QU)) * dxinv;

    } else if (coord_type == 1) {
      // axisymmetric
      if (i == 0) {
        div(i,j,k) = 0.0_rt;
      } else {
        Real rl = (i - 0.5_rt) * dx[0] + problo[0];
        Real rr = (i + 0.5_rt) * dx[0] + problo[0];
        Real rc = (i) * dx[0] + problo[0];

        div(i,j,k) = (rr * q_arr(i,j,k,QU) - rl * q_arr(i-1,j,k,QU)) * dxinv / rc;
      }
    } else {
      // spherical
      if (i == 0) {
        div(i,j,k) = 0.0_rt;
      } else {
        Real rl = (i - 0.5_rt) * dx[0] + problo[0];
        Real rr = (i + 0.5_rt) * dx[0] + problo[0];
        Real rc = (i) * dx[0] + problo[0];

        div(i,j,k) = (rr * rr * q_arr(i,j,k,QU) - rl * rl * q_arr(i-1,j,k,QU)) * dxinv / (rc * rc);
      }
    }
#endif

#if AMREX_SPACEDIM == 2
    Real ux = 0.0_rt;
    Real vy = 0.0_rt;

    if (coord_type == 0) {
      ux = 0.5_rt * (q_arr(i,j,k,QU) - q_arr(i-1,j,k,QU) + q_arr(i,j-1,k,QU) - q_arr(i-1,j-1,k,QU)) * dxinv;
      vy = 0.5_rt * (q_arr(i,j,k,QV) - q_arr(i,j-1,k,QV) + q_arr(i-1,j,k,QV) - q_arr(i-1,j-1,k,QV)) * dyinv;

    } else {
      if (i == 0) {
        ux = 0.0_rt;
        vy = 0.0_rt;  // is this part correct?
      } else {
        Real rl = (i - 0.5_rt) * dx[0] + problo[0];
        Real rr = (i + 0.5_rt) * dx[0] + problo[0];
        Real rc = (i) * dx[0] + problo[0];

        // These are transverse averages in the y-direction
        Real ul = 0.5_rt * (q_arr(i-1,j,k,QU) + q_arr(i-1,j-1,k,QU));
        Real ur = 0.5_rt * (q_arr(i,j,k,QU) + q_arr(i,j-1,k,QU));

        // Take 1/r d/dr(r*u)
        ux = (rr * ur - rl * ul) * dxinv / rc;

        // These are transverse averages in the x-direction
        Real vb = 0.5_rt * (q_arr(i,j-1,k,QV) + q_arr(i-1,j-1,k,QV));
        Real vt = 0.5_rt * (q_arr(i,j,k,QV) + q_arr(i-1,j,k,QV));

        vy = (vt - vb) * dyinv;
      }
    }

    div(i,j,k) = ux + vy;
#endif

#if AMREX_SPACEDIM == 3
    Real ux = 0.25_rt * (q_arr(i,j,k,QU) - q_arr(i-1,j,k,QU) +
                         q_arr(i,j,k-1,QU) - q_arr(i-1,j,k-1,QU) +
                         q_arr(i,j-1,k,QU) - q_arr(i-1,j-1,k,QU) +
                         q_arr(i,j-1,k-1,QU) - q_arr(i-1,j-1,k-1,QU)) * dxinv;

    Real vy = 0.25_rt * (q_arr(i,j,k,QV) - q_arr(i,j-1,k,QV) +
                         q_arr(i,j,k-1,QV) - q_arr(i,j-1,k-1,QV) +
                         q_arr(i-1,j,k,QV) - q_arr(i-1,j-1,k,QV) +
                         q_arr(i-1,j,k-1,QV) - q_arr(i-1,j-1,k-1,QV)) * dyinv;

    Real wz = 0.25_rt * (q_arr(i,j,k,QW) - q_arr(i,j,k-1,QW) +
                         q_arr(i,j-1,k,QW) - q_arr(i,j-1,k-1,QW) +
                         q_arr(i-1,j,k,QW) - q_arr(i-1,j,k-1,QW) +
                         q_arr(i-1,j-1,k,QW) - q_arr(i-1,j-1,k-1,QW)) * dzinv;

    div(i,j,k) = ux + vy + wz;
#endif

  });

}


void
Castro::apply_av(const Box& bx,
                 const int idir,
                 Array4<Real const> const div,
                 Array4<Real const> const uin,
                 Array4<Real> const flux) {

  const auto dx = geom.CellSizeArray();

  Real diff_coeff = difmag;

  AMREX_PARALLEL_FOR_4D(bx, NUM_STATE, i, j, k, n,
  {

    if (n == UTEMP) continue;
#ifdef SHOCK_VAR
    if (n == USHK) continue;
#endif

    Real div1;
    if (idir == 0) {

      div1 = 0.25_rt * (div(i,j,k) + div(i,j+dg1,k) +
                        div(i,j,k+dg2) + div(i,j+dg1,k+dg2));
      div1 = diff_coeff * amrex::min(0.0_rt, div1);
      div1 = div1 * (uin(i,j,k,n) - uin(i-1,j,k,n));

    } else if (idir == 1) {

      div1 = 0.25_rt * (div(i,j,k) + div(i+1,j,k) +
                        div(i,j,k+dg2) + div(i+1,j,k+dg2));
      div1 = diff_coeff * amrex::min(0.0_rt, div1);
      div1 = div1 * (uin(i,j,k,n) - uin(i,j-dg1,k,n));

    } else {

      div1 = 0.25_rt * (div(i,j,k) + div(i+1,j,k) +
                        div(i,j+dg1,k) + div(i+1,j+dg1,k));
      div1 = diff_coeff * amrex::min(0.0_rt, div1);
      div1 = div1 * (uin(i,j,k,n) - uin(i,j,k-dg2,n));

    }

    flux(i,j,k,n) += dx[idir] * div1;
  });
}


#ifdef RADIATION
void
Castro::apply_av_rad(const Box& bx,
                     const int idir,
                     Array4<Real const> const div,
                     Array4<Real const> const Erin,
                     Array4<Real> const radflux) {

  const auto dx = geom.CellSizeArray();

  Real diff_coeff = difmag;

  AMREX_PARALLEL_FOR_4D(bx, Radiation::nGroups, i, j, k, n,
  {

    Real div1;
    if (idir == 0) {

      div1 = 0.25_rt * (div(i,j,k) + div(i,j+dg1,k) +
                        div(i,j,k+dg2) + div(i,j+dg1,k+dg2));
      div1 = diff_coeff * amrex::min(0.0_rt, div1);
      div1 = div1 * (Erin(i,j,k,n) - Erin(i-1,j,k,n));

    } else if (idir == 1) {

      div1 = 0.25_rt * (div(i,j,k) + div(i+1,j,k) +
                        div(i,j,k+dg2) + div(i+1,j,k+dg2));
      div1 = diff_coeff * amrex::min(0.0_rt, div1);
      div1 = div1 * (Erin(i,j,k,n) - Erin(i,j-dg1,k,n));

    } else {

      div1 = 0.25_rt * (div(i,j,k) + div(i+1,j,k) +
                        div(i,j+dg1,k) + div(i+1,j+dg1,k));
      div1 = diff_coeff * amrex::min(0.0_rt, div1);
      div1 = div1 * (Erin(i,j,k,n) - Erin(i,j,k-dg2,n));

    }

    radflux(i,j,k,n) += dx[idir] * div1;
  });
}
#endif


void
Castro::normalize_species_fluxes(const Box& bx,
                                 Array4<Real> const flux) {

  // Normalize the fluxes of the mass fractions so that
  // they sum to 0.  This is essentially the CMA procedure that is
  // defined in Plewa & Muller, 1999, A&A, 342, 179.

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    Real sum = 0.0_rt;

    for (int n = UFS; n < UFS+NumSpec; n++) {
      sum += flux(i,j,k,n);
    }

    Real fac = 1.0_rt;
    if (sum != 0.0_rt) {
      fac = flux(i,j,k,URHO) / sum;
    }

    for (int n = UFS; n < UFS+NumSpec; n++) {
      flux(i,j,k,n) = flux(i,j,k,n) * fac;
    }
  });
}


void
Castro::scale_flux(const Box& bx,
#if AMREX_SPACEDIM == 1
                   Array4<Real const> const qint,
#endif
                   Array4<Real> const flux,
                   Array4<Real const> const area_arr,
                   const Real dt) {

#if AMREX_SPACEDIM == 1
  const int coord_type = geom.Coord();
#endif

  AMREX_PARALLEL_FOR_4D(bx, NUM_STATE, i, j, k, n,
  {

    flux(i,j,k,n) = dt * flux(i,j,k,n) * area_arr(i,j,k);
#if AMREX_SPACEDIM == 1
    // Correct the momentum flux with the grad p part.
    if (coord_type == 0 && n == UMX) {
      flux(i,j,k,n) += dt * area_arr(i,j,k) * qint(i,j,k,GDPRES);
    }
#endif
  });
}


#ifdef RADIATION
void
Castro::scale_rad_flux(const Box& bx,
                       Array4<Real> const rflux,
                       Array4<Real const> area_arr,
                       const Real dt) {

  AMREX_PARALLEL_FOR_4D(bx, Radiation::nGroups, i, j, k, g,
  {
    rflux(i,j,k,g) = dt * rflux(i,j,k,g) * area_arr(i,j,k);
  });
}
#endif

