#include <Castro.H>
#include <Castro_F.H>

#ifdef ROTATION
#include <Rotation.H>
#endif

#include <Castro_util.H>
#include <advection_util.H>

#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif

#ifdef RADIATION
#include <Radiation.H>
#include <fluxlimiter.H>
#include <rad_util.H>
#endif

#include <eos.H>

using namespace amrex;


void
Castro::ctoprim(const Box& bx,
                const Real time,
                Array4<Real const> const& uin,
#ifdef MHD
                Array4<Real const> const& Bx,
                Array4<Real const> const& By,
                Array4<Real const> const& Bz,
#endif
#ifdef RADIATION
                Array4<Real const> const& Erin,
                Array4<Real const> const& lam,
#endif
                Array4<Real> const& q_arr,
                Array4<Real> const& qaux_arr) {

#ifdef RADIATION
  int is_comoving = Radiation::comoving;
  int limiter = Radiation::limiter;
  int closure = Radiation::closure;
#endif

#ifdef ROTATION
  GeometryData geomdata = geom.data();
#endif

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {

#ifndef AMREX_USE_GPU
    if (uin(i,j,k,URHO) <= 0.0_rt) {
      std::cout << std::endl;
      std::cout << ">>> Error: advection_util_nd.F90::ctoprim " << i << " " << j << " " << k << std::endl;
      std::cout << ">>> ... negative density " << uin(i,j,k,URHO) << std::endl;
      amrex::Error("Error:: advection_util_nd.f90 :: ctoprim");
    } else if (uin(i,j,k,URHO) < castro::small_dens) {
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

#ifdef MHD
    q_arr(i,j,k,QMAGX) = 0.5_rt * (Bx(i+1,j,k) + Bx(i,j,k));
    q_arr(i,j,k,QMAGY) = 0.5_rt * (By(i,j+1,k) + By(i,j,k));
    q_arr(i,j,k,QMAGZ) = 0.5_rt * (Bz(i,j,k+1) + Bz(i,j,k));
#endif

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

    if ((uin(i,j,k,UEDEN) - kineng) > castro::dual_energy_eta1*uin(i,j,k,UEDEN)) {
      q_arr(i,j,k,QREINT) = (uin(i,j,k,UEDEN) - kineng) * rhoinv;
    } else {
      q_arr(i,j,k,QREINT) = uin(i,j,k,UEINT) * rhoinv;
    }

    // If we're advecting in the rotating reference frame,
    // then subtract off the rotation component here.

#ifdef ROTATION
    if (castro::do_rotation == 1 && castro::state_in_rotating_frame != 1) {
      GpuArray<Real, 3> vel;
      for (int n = 0; n < 3; n++) {
        vel[n] = uin(i,j,k,UMX+n) * rhoinv;
      }

      inertial_to_rotational_velocity(i, j, k, geomdata, time, vel);

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
      int n  = upassmap(ipassive);
      int iq = qpassmap(ipassive);
      q_arr(i,j,k,iq) = uin(i,j,k,n) * rhoinv;
    }

    // get gamc, p, T, c, csml using q state
    eos_rep_t eos_state;
    eos_state.T = q_arr(i,j,k,QTEMP);
    eos_state.rho = q_arr(i,j,k,QRHO);
    eos_state.e = q_arr(i,j,k,QREINT);
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n]  = q_arr(i,j,k,QFS+n);
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = q_arr(i,j,k,QFX+n);
    }
#endif

    eos(eos_input_re, eos_state);

    q_arr(i,j,k,QTEMP) = eos_state.T;
    q_arr(i,j,k,QREINT) = eos_state.e * q_arr(i,j,k,QRHO);
    q_arr(i,j,k,QPRES) = eos_state.p;
#ifdef TRUE_SDC
    q_arr(i,j,k,QGC) = eos_state.gam1;
#endif

#ifdef MHD
    q_arr(i,j,k,QPTOT) = q_arr(i,j,k,QPRES) +
      0.5_rt * (q_arr(i,j,k,QMAGX) * q_arr(i,j,k,QMAGX) +
                q_arr(i,j,k,QMAGY) * q_arr(i,j,k,QMAGY) +
                q_arr(i,j,k,QMAGZ) * q_arr(i,j,k,QMAGZ));
#endif

#ifdef RADIATION
    qaux_arr(i,j,k,QGAMCG) = eos_state.gam1;
    qaux_arr(i,j,k,QCG) = eos_state.cs;

    Real lams[NGROUPS];
    for (int g = 0; g < NGROUPS; g++) {
      lams[g] = lam(i,j,k,g);
    }
    Real qs[NQ];
    for (int n = 0; n < NQ; n++) {
      qs[n] = q_arr(i,j,k,n);
    }
    Real ptot;
    Real ctot;
    Real gamc_tot;
    compute_ptot_ctot(lams, qs,
                      is_comoving, limiter, closure,
                      qaux_arr(i,j,k,QCG),
                      ptot, ctot, gamc_tot);

    q_arr(i,j,k,QPTOT) = ptot;

    qaux_arr(i,j,k,QC) = ctot;
    qaux_arr(i,j,k,QGAMC) = gamc_tot;

    q_arr(i,j,k,QREITOT) = q_arr(i,j,k,QREINT);
    for (int g = 0; g < NGROUPS; g++) {
      qaux_arr(i,j,k,QLAMS+g) = lam(i,j,k,g);
      q_arr(i,j,k,QREITOT) += q_arr(i,j,k,QRAD+g);
    }

#else
    qaux_arr(i,j,k,QGAMC) = eos_state.gam1;
    qaux_arr(i,j,k,QC) = eos_state.cs;
#endif

  });
}


void
Castro::shock(const Box& bx,
              Array4<Real const> const& q_arr,
              Array4<Real> const& shk) {

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

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
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

#if AMREX_SPACEDIM <= 2
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
#endif

#if AMREX_SPACEDIM == 1
    } else if (coord_type == 2) {

      // 1-d spherical
      Real rc = (i + 0.5_rt) * dx[0];
      Real rm = (i - 1 + 0.5_rt) * dx[0];
      Real rp = (i + 1 + 0.5_rt) * dx[0];

      div_u += 0.5_rt * (rp * rp * q_arr(i+1,j,k,QU) - rm * rm * q_arr(i-1,j,k,QU)) / (rc * rc * dx[0]);
#endif

#ifndef AMREX_USE_GPU

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
             Array4<Real const> const& q_arr,
             Array4<Real> const& div) {
  // this computes the *node-centered* divergence

  const auto dx = geom.CellSizeArray();

#if AMREX_SPACEDIM <= 2
  const int coord_type = geom.Coord();
#endif

  const auto problo = geom.ProbLoArray();

  Real dxinv = 1.0_rt / dx[0];
#if AMREX_SPACEDIM >= 2
  Real dyinv = 1.0_rt / dx[1];
#else
  Real dyinv = 0.0_rt;
#endif
#if AMREX_SPACEDIM == 3
  Real dzinv = 1.0_rt / dx[2];
#endif

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
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
                 Array4<Real const> const& div,
                 Array4<Real const> const& uin,
                 Array4<Real> const& flux) {

  const auto dx = geom.CellSizeArray();

  Real diff_coeff = difmag;

  amrex::ParallelFor(bx, NUM_STATE,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
  {

    if (n == UTEMP) return;
#ifdef SHOCK_VAR
    if (n == USHK) return;
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
                     Array4<Real const> const& div,
                     Array4<Real const> const& Erin,
                     Array4<Real> const& radflux) {

  const auto dx = geom.CellSizeArray();

  Real diff_coeff = difmag;

  amrex::ParallelFor(bx, Radiation::nGroups,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
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
                                 Array4<Real> const& flux) {

  // Normalize the fluxes of the mass fractions so that
  // they sum to 0.  This is essentially the CMA procedure that is
  // defined in Plewa & Muller, 1999, A&A, 342, 179.

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {

    Real sum = 0.0_rt;

    for (int n = UFS; n < UFS+NumSpec; n++) {
      sum += flux(i,j,k,n);
    }

    Real fac = 1.0_rt;

    // We skip the normalization if the sum is zero or within epsilon.
    // There can be numerical problems here if the density flux is
    // approximately zero at the interface but not exactly, resulting in
    // division by a small number and/or resulting in one of the species
    // fluxes being negative because of roundoff error. There are also other
    // terms like artificial viscosity which can cause these problems.
    // So checking that sum is sufficiently large helps avoid this.

    if (std::abs(sum) > std::numeric_limits<Real>::epsilon() * std::abs(flux(i,j,k,URHO))) {
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
                   Array4<Real const> const& qint,
#endif
                   Array4<Real> const& flux,
                   Array4<Real const> const& area_arr,
                   const Real dt) {

#if AMREX_SPACEDIM == 1
  const int coord_type = geom.Coord();
#endif

  amrex::ParallelFor(bx, NUM_STATE,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
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
                       Array4<Real> const& rflux,
                       Array4<Real const> const& area_arr,
                       const Real dt) {

  amrex::ParallelFor(bx, Radiation::nGroups,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int g)
  {
    rflux(i,j,k,g) = dt * rflux(i,j,k,g) * area_arr(i,j,k);
  });
}
#endif



void
Castro::limit_hydro_fluxes_on_small_dens(const Box& bx,
                                         int idir,
                                         Array4<Real const> const& u,
                                         Array4<Real const> const& q,
                                         Array4<Real const> const& vol,
                                         Array4<Real> const& flux,
                                         Array4<Real const> const& area_arr,
                                         Real dt)
{

    // The following algorithm comes from Hu, Adams, and Shu (2013), JCP, 242, 169,
    // "Positivity-preserving method for high-order conservative schemes solving
    // compressible Euler equations." It has been modified to enforce not only positivity
    // but also the stronger requirement that rho > small_dens. We do not limit on pressure
    // (or, similarly, internal energy) because those cases are easily fixed by calls to
    // reset_internal_energy that enforce a thermodynamic floor. The density limiter, by
    // contrast, is very important because calls to enforce_minimum_density can yield
    // hydrodynamic states that are inconsistent (there is no clear strategy for what to do
    // when a density is negative).

    const Real density_floor_tolerance = 1.1_rt;

    // The density floor is the small density, modified by a small factor.
    // In practice numerical error can cause the density that is created
    // by this flux limiter to be slightly lower than the target density,
    // so we set the target to be slightly larger than the real density floor
    // to avoid density resets.

    Real density_floor = small_dens * density_floor_tolerance;

    // We apply this flux limiter on a per-edge basis. So we can guarantee
    // that any individual flux cannot cause a small density in one step,
    // but with the above floor we cannot guarantee that the sum of the
    // fluxes will enforce this constraint. The only way to guarantee that
    // is if the density floor is increased by a factor of the number of
    // edges, so that even if all edges are summed together, the density
    // will still be at the floor. So we multiply the floor by a factor of
    // 2 (two edges in each dimension) and a factor of AMREX_SPACEDIM.

    density_floor *= AMREX_SPACEDIM * 2;

    const Real* dx = geom.CellSize();

    Real dtdx = dt / dx[idir];
    Real lcfl = cfl;
    Real alpha = 1.0_rt / AMREX_SPACEDIM;

    auto coord = geom.Coord();
    GeometryData geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        // Grab the states on either side of the interface we are working with,
        // depending on which dimension we're currently calling this with.

        GpuArray<Real, NUM_STATE> uR;
        for (int n = 0; n < NUM_STATE; ++n) {
            uR[n] = u(i,j,k,n);
        }

        GpuArray<Real, NQ> qR;
        for (int n = 0; n < NQ; ++n) {
            qR[n] = q(i,j,k,n);
        }

        Real volR = vol(i,j,k);

        GpuArray<int, 3> idxR = {i,j,k};

        GpuArray<Real, NUM_STATE> uL;
        GpuArray<Real, NQ> qL;
        Real volL;
        GpuArray<int, 3> idxL;

        if (idir == 0) {

            for (int n = 0; n < NUM_STATE; ++n) {
                uL[n] = u(i-1,j,k,n);
            }

            for (int n = 0; n < NQ; ++n) {
                qL[n] = q(i-1,j,k,n);
            }

            volL = vol(i-1,j,k);

            idxL = {i-1,j,k};

        }
        else if (idir == 1) {

            for (int n = 0; n < NUM_STATE; ++n) {
                uL[n] = u(i,j-1,k,n);
            }

            for (int n = 0; n < NQ; ++n) {
                qL[n] = q(i,j-1,k,n);
            }

            volL = vol(i,j-1,k);

            idxL = {i,j-1,k};

        }
        else {

            for (int n = 0; n < NUM_STATE; ++n) {
                uL[n] = u(i,j,k-1,n);
            }

            for (int n = 0; n < NQ; ++n) {
                qL[n] = q(i,j,k-1,n);
            }

            volL = vol(i,j,k-1);

            idxL = {i,j,k-1};

        }

        // If an adjacent zone has a floor-violating density, set the flux to zero and move on.
        // At that point, the only thing to do is wait for a reset at a later point.

        if (uR[URHO] < density_floor || uL[URHO] < density_floor) {

            for (int n = 0; n < NUM_STATE; ++n) {
                flux(i,j,k,n) = 0.0_rt;
            }

            return;
        }

        // Construct cell-centered fluxes.

        GpuArray<Real, NUM_STATE> fluxL;
        dflux(uL, qL, idir, coord, geomdata, idxL, fluxL);

        GpuArray<Real, NUM_STATE> fluxR;
        dflux(uR, qR, idir, coord, geomdata, idxR, fluxR);

        // Construct the Lax-Friedrichs flux on the interface (Equation 12).
        // Note that we are using the information from Equation 9 to obtain the
        // effective maximum wave speed, (|u| + c)_max = CFL / lambda where
        // lambda = dt/(dx * alpha); alpha = 1 in 1D and may be chosen somewhat
        // freely in multi-D as long as alpha_x + alpha_y + alpha_z = 1.

        GpuArray<Real, NUM_STATE> fluxLF;
        for (int n = 0; n < NUM_STATE; ++n) {
            fluxLF[n] = 0.5_rt * (fluxL[n] + fluxR[n] + (lcfl / dtdx / alpha) * (uL[n] - uR[n]));
        }

        // Coefficients of fluxes on either side of the interface.

        Real flux_coefR = 2.0_rt * (dt / alpha) * area_arr(i,j,k) / volR;
        Real flux_coefL = 2.0_rt * (dt / alpha) * area_arr(i,j,k) / volL;

        // Obtain the one-sided update to the density, based on Hu et al., Eq. 11.
        // If we would violate the floor, then we need to limit the flux. Since the
        // flux adds to the density on one side and subtracts from the other, the floor
        // can only be violated in at most one direction, so we'll do an if-else test
        // below. This means that we can simplify the approach of Hu et al. -- whereas
        // they constructed two thetas for each interface (corresponding to either side)
        // we can complete the operation in one step with a single theta.

        Real drhoL = flux_coefL * flux(i,j,k,URHO);
        Real rhoL = uL[URHO] - drhoL;

        Real drhoR = flux_coefR * flux(i,j,k,URHO);
        Real rhoR = uR[URHO] + drhoR;

        Real theta = 1.0_rt;

        if (rhoL < density_floor) {

            // Obtain the final density corresponding to the LF flux.

            Real drhoLF = flux_coefL * fluxLF[URHO];
            Real rhoLF = uL[URHO] - drhoLF;

            // Solve for theta from (1 - theta) * rhoLF + theta * rho = density_floor.

            theta = amrex::min(theta, (density_floor - rhoLF) / (rhoL - rhoLF));

        }
        else if (rhoR < density_floor) {

            Real drhoLF = flux_coefR * fluxLF[URHO];
            Real rhoLF = uR[URHO] + drhoLF;

            theta = amrex::min(theta, (density_floor - rhoLF) / (rhoR - rhoLF));

        }

        // Limit theta to the valid range (this will deal with roundoff issues).

        theta = amrex::min(1.0_rt, amrex::max(theta, 0.0_rt));

        // Assemble the limited flux (Equation 16).

        for (int n = 0; n < NUM_STATE; ++n) {
            flux(i,j,k,n) = (1.0_rt - theta) * fluxLF[n] + theta * flux(i,j,k,n);
        }

        // Zero out fluxes for quantities that don't advect.

        flux(i,j,k,UTEMP) = 0.0_rt;
#ifdef SHOCK_VAR
        flux(i,j,k,USHK) = 0.0_rt;
#endif

        // Now, apply our requirement that the final flux cannot violate the density floor.

        drhoR = flux_coefR * flux(i,j,k,URHO);
        drhoL = flux_coefL * flux(i,j,k,URHO);

        if (uR[URHO] + drhoR < density_floor) {
            for (int n = 0; n < NUM_STATE; ++n) {
                flux(i,j,k,n) = flux(i,j,k,n) * std::abs((density_floor - uR[URHO]) / drhoR);
            }
        }
        else if (uL[URHO] - drhoL < density_floor) {
            for (int n = 0; n < NUM_STATE; ++n) {
                flux(i,j,k,n) = flux(i,j,k,n) * std::abs((density_floor - uL[URHO]) / drhoL);
            }
        }

    });

}



void
Castro::limit_hydro_fluxes_on_large_vel(const Box& bx,
                                        int idir,
                                        Array4<Real const> const& u,
                                        Array4<Real const> const& q,
                                        Array4<Real const> const& vol,
                                        Array4<Real> const& flux,
                                        Array4<Real const> const& area_arr,
                                        Real dt)
{

    // This limiter is similar to the density-based limiter above, but limits
    // on velocities that are too large instead. The comments are minimal since
    // the algorithm is effectively the same.

    if (castro::speed_limit <= 0.0_rt) return;

    const Real* dx = geom.CellSize();

    Real dtdx = dt / dx[idir];
    Real lcfl = cfl;
    Real alpha = 1.0_rt / AMREX_SPACEDIM;

    auto coord = geom.Coord();
    GeometryData geomdata = geom.data();

    Real lspeed_limit = speed_limit / (2 * AMREX_SPACEDIM);

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        GpuArray<Real, NUM_STATE> uR;
        for (int n = 0; n < NUM_STATE; ++n) {
            uR[n] = u(i,j,k,n);
        }

        GpuArray<Real, NQ> qR;
        for (int n = 0; n < NQ; ++n) {
            qR[n] = q(i,j,k,n);
        }

        Real volR = vol(i,j,k);

        GpuArray<int, 3> idxR = {i,j,k};

        GpuArray<Real, NUM_STATE> uL;
        GpuArray<Real, NQ> qL;
        Real volL;
        GpuArray<int, 3> idxL;

        if (idir == 0) {

            for (int n = 0; n < NUM_STATE; ++n) {
                uL[n] = u(i-1,j,k,n);
            }

            for (int n = 0; n < NQ; ++n) {
                qL[n] = q(i-1,j,k,n);
            }

            volL = vol(i-1,j,k);

            idxL = {i-1,j,k};

        }
        else if (idir == 1) {

            for (int n = 0; n < NUM_STATE; ++n) {
                uL[n] = u(i,j-1,k,n);
            }

            for (int n = 0; n < NQ; ++n) {
                qL[n] = q(i,j-1,k,n);
            }

            volL = vol(i,j-1,k);

            idxL = {i,j-1,k};

        }
        else {

            for (int n = 0; n < NUM_STATE; ++n) {
                uL[n] = u(i,j,k-1,n);
            }

            for (int n = 0; n < NQ; ++n) {
                qL[n] = q(i,j,k-1,n);
            }

            volL = vol(i,j,k-1);

            idxL = {i,j,k-1};

        }

        // Construct cell-centered fluxes.

         GpuArray<Real, NUM_STATE> fluxL;
         dflux(uL, qL, idir, coord, geomdata, idxL, fluxL);

         GpuArray<Real, NUM_STATE> fluxR;
         dflux(uR, qR, idir, coord, geomdata, idxR, fluxR);

         // Construct the Lax-Friedrichs flux on the interface.

         GpuArray<Real, NUM_STATE> fluxLF;
         for (int n = 0; n < NUM_STATE; ++n) {
             fluxLF[n] = 0.5_rt * (fluxL[n] + fluxR[n] + (lcfl / dtdx / alpha) * (uL[n] - uR[n]));
         }

         // Coefficients of fluxes on either side of the interface.

         Real flux_coefR = 2.0_rt * (dt / alpha) * area_arr(i,j,k) / volR;
         Real flux_coefL = 2.0_rt * (dt / alpha) * area_arr(i,j,k) / volL;

         Real theta = 1.0_rt;

         // Loop over all three momenta, and choose the strictest
         // limiter among them.

         for (int n = 0; n < 3; ++n) {

             int UMOM = UMX + n;

             // Obtain the one-sided update to the momentum.

             Real drhouL = flux_coefL * flux(i,j,k,UMOM);
             Real rhouL = std::abs(uL[UMOM] - drhouL);

             Real drhoL = flux_coefL * flux(i,j,k,URHO);
             Real rhoL = uL[URHO] - drhoL;

             Real drhouR = flux_coefR * flux(i,j,k,UMOM);
             Real rhouR = std::abs(uR[UMOM] + drhouR);

             Real drhoR = flux_coefR * flux(i,j,k,URHO);
             Real rhoR = uR[URHO] + drhoR;

             if (std::abs(rhouL) > rhoL * lspeed_limit) {

                 // Obtain the final density corresponding to the LF flux.

                 Real drhouLF = flux_coefL * fluxLF[UMOM];
                 Real rhouLF = std::abs(uL[UMOM] - drhouLF);

                 // Solve for theta from (1 - theta) * rhouLF + theta * rhou = rhoL * speed_limit.

                 theta = amrex::min(theta, std::abs(rhoL * lspeed_limit - rhouLF) / std::abs(rhouL - rhouLF));

             }
             else if (std::abs(rhouR) > rhoR * lspeed_limit) {

                 Real drhouLF = flux_coefR * fluxLF[UMOM];
                 Real rhouLF = std::abs(uR[UMOM] + drhouLF);

                 theta = amrex::min(theta, std::abs(rhoR * lspeed_limit - rhouLF) / std::abs(rhouR - rhouLF));

             }

         }

         // Limit theta to the valid range (this will deal with roundoff issues).

         theta = amrex::min(1.0_rt, amrex::max(theta, 0.0_rt));

         // Assemble the limited flux (Equation 16).

         for (int n = 0; n < NUM_STATE; ++n) {
             flux(i,j,k,n) = (1.0_rt - theta) * fluxLF[n] + theta * flux(i,j,k,n);
         }

         // Zero out fluxes for quantities that don't advect.

         flux(i,j,k,UTEMP) = 0.0_rt;
#ifdef SHOCK_VAR
         flux(i,j,k,USHK) = 0.0_rt;
#endif

    });

}


void
Castro::do_enforce_minimum_density(const Box& bx,
                                   Array4<Real> const& state_arr,
                                   const int verbose) {

#ifdef HYBRID_MOMENTUM
  GeometryData geomdata = geom.data();
#endif

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {

    if (state_arr(i,j,k,URHO) < small_dens) {

#ifndef AMREX_USE_GPU
      if (verbose > 0) {
        std::cout << " " << std::endl;
        if (state_arr(i,j,k,URHO) < 0.0_rt) {
          std::cout << ">>> RESETTING NEG.  DENSITY AT " << i << ", " << j << ", " << k << std::endl;
        }
        else if (state_arr(i,j,k,URHO) == 0.0_rt) {
          // If the density is *exactly* zero, that almost certainly means something has gone wrong,
          // like we failed to properly fill the state data on grid creation.
          amrex::Error("Density exactly zero at " + std::to_string(i) + ", " +
                                                    std::to_string(j) + ", " +
                                                    std::to_string(k));
        }
        else {
          std::cout << ">>> RESETTING SMALL DENSITY AT " << i << ", " << j << ", " << k << std::endl;
        }
        std::cout << ">>> FROM " << state_arr(i,j,k,URHO) << " TO " << small_dens << std::endl;
        std::cout << ">>> IN GRID " << bx << std::endl;
        std::cout << " " << std::endl;
      }
#endif

      for (int ipassive = 0; ipassive < npassive; ipassive++) {
        int n = upassmap(ipassive);
        state_arr(i,j,k,n) *= (small_dens / state_arr(i,j,k,URHO));
      }

      eos_re_t eos_state;
      eos_state.rho = small_dens;
      eos_state.T = small_temp;
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state_arr(i,j,k,UFS+n) / small_dens;
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = state_arr(i,j,k,UFX+n) / small_dens;
      }
#endif

      eos(eos_input_rt, eos_state);

      state_arr(i,j,k,URHO ) = eos_state.rho;
      state_arr(i,j,k,UTEMP) = eos_state.T;

      state_arr(i,j,k,UMX) = 0.0_rt;
      state_arr(i,j,k,UMY) = 0.0_rt;
      state_arr(i,j,k,UMZ) = 0.0_rt;

      state_arr(i,j,k,UEINT) = eos_state.rho * eos_state.e;
      state_arr(i,j,k,UEDEN) = state_arr(i,j,k,UEINT);

#ifdef HYBRID_MOMENTUM
      GpuArray<Real, 3> loc;

      position(i, j, k, geomdata, loc);

      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        loc[dir] -= problem::center[dir];
      }

      GpuArray<Real, 3> linear_mom;

      for (int dir = 0; dir < 3; ++dir) {
        linear_mom[dir] = state_arr(i,j,k,UMX+dir);
      }

      GpuArray<Real, 3> hybrid_mom;

      linear_to_hybrid(loc, linear_mom, hybrid_mom);

      for (int dir = 0; dir < 3; ++dir) {
        state_arr(i,j,k,UMR+dir) = hybrid_mom[dir];
      }
#endif
    }
  });
}
