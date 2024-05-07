#ifndef RADHYDRO_H
#define RADHYDRO_H

// do the advection in frequency space for the radiation energy.
// see Paper III, section 2.4

constexpr int rk_order = 3;
constexpr bool use_WENO = false;

constexpr Real cfl_fspace = 0.5_rt;

constexpr Real onethird = 1.0_rt/3.0_rt;
constexpr Real twothirds = 2.0_rt/3.0_rt;
constexpr Real onesixth = 1.0_rt/6.0_rt;

// RK5
constexpr Real B1 = 0.5_rt;
constexpr Real B2 = 1.0_rt/16.0_rt;
constexpr Real B3 = 0.5_rt;
constexpr Real B4 = 9.0_rt/16.0_rt;
constexpr Real B5 = 8.0_rt/7.0_rt;
constexpr Real B6 = 7.0_rt/90.0_rt;

constexpr Real C20 = 5.0_rt/8.0_rt;
constexpr Real C21 = 3.0_rt/8.0_rt;

constexpr Real C40 = 17.0_rt/8.0_rt;
constexpr Real C41 = 9.0_rt/8.0_rt;
constexpr Real C42 = -3.0_rt;
constexpr Real C43 = 0.75_rt;

constexpr Real C50 = -5.0_rt/21.0_rt;
constexpr Real C51 = 2.0_rt/7.0_rt;
constexpr Real C52 = 0.0_rt;
constexpr Real C53 = 4.0_rt;
constexpr Real C54 = -64.0_rt/21.0_rt;

constexpr Real C60 = -8.0_rt/27.0_rt;
constexpr Real C61 = -1.0_rt/5.0_rt;
constexpr Real C62 = 32.0_rt/45.0_rt;
constexpr Real C63 = -32.0_rt/45.0_rt;
constexpr Real C64 = 32.0_rt/27.0_rt;
constexpr Real C65 = 14.0_rt/45.0_rt;


AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
MC(const Real r) {
  return amrex::max(0.0_rt, amrex::min(2.0_rt * r, 0.5_rt * (1.0_rt + r), 2.0_rt));
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
weno5(const Real vm2, const Real vm1, const Real v, const Real vp1, const Real vp2) {

  constexpr Real epsw = 1.0e-6_rt;
  constexpr Real b1 = 13.0_rt/12.0_rt;
  constexpr Real b2 = 1.0_rt/6.0_rt;

  Real djm1 = vm2 - 2.0_rt * vm1 + v;
  Real ejm1 = vm2 - 4.0_rt * vm1 + 3.0_rt * v;
  Real dj = vm1 - 2.0_rt * v + vp1;
  Real ej = vm1 - vp1;
  Real djp1 = v - 2.0_rt * vp1 + vp2;
  Real ejp1 = 3.0_rt * v - 4.0_rt * vp1 + vp2;

  Real dis0 = b1 * djm1 * djm1 + 0.25_rt * ejm1 * ejm1 + epsw;
  Real dis1 = b1 * dj * dj + 0.25_rt * ej * ej + epsw;
  Real dis2 = b1 * djp1 * djp1 + 0.25_rt * ejp1 * ejp1 + epsw;

  Real q30 = 2.0_rt *vm2 - 7.0_rt * vm1 + 11.0_rt * v;
  Real q31 = -vm1 + 5.0_rt * v + 2.0_rt * vp1;
  Real q32 = 2.0_rt * v + 5.0_rt * vp1 - vp2;

  Real d01 = dis0 / dis1;
  Real d02 = dis0 / dis2;
  Real a1ba0 = 6.0_rt * d01 * d01;
  Real a2ba0 = 3.0_rt * d02 * d02;
  Real w0 = 1.0_rt / (1.0_rt + a1ba0 + a2ba0);
  Real w1 = a1ba0 * w0;
  Real w2 = 1.0_rt - w0 - w1;

  if (w0 <  1.0e-10_rt) {
    w0 = 0.0_rt;
  }
  if (w1 < 1.0e-10_rt) {
    w1 = 0.0_rt;
  }
  if (w2 < 1.0e-10_rt) {
    w2 = 0.e0_rt;
  }

  return b2 * (w0 * q30 + w1 * q31 + w2 * q32);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
dudt(Real* u, Real* a, const Real* dx, const int n, Real* dudt_tmp) {

  // compute the RHS for the g'th energy group out of n note, since
  // this is advection in frequency space, we will need to assume some
  // boundary conditions on the stencil used to reconstruct over nu

  Array1D<Real, 0, NGROUPS> f;

  // we need 2 ghost cells on each end
  Array1D<Real, -2, NGROUPS+1> ag;
  Array1D<Real, -2, NGROUPS+1> ug;
  Array1D<Real, -2, NGROUPS+1> fg;

  // ag and ug have 2 ghost cells on each end.  First fill the
  // interiors and then the ghost cells.
  for (int g = 0; g < n; g++) {
    ag(g) = a[g];
    ug(g) = u[g];
  }

  ag(-2) = -a[1];
  ag(-1) = -a[0];

  ag(n) = -ag(n-1);
  ag(n+1) = -ag(n-2);

  ug(-2) = u[1];
  ug(-1) = u[0];

  ug(n) = ug(n-1);
  ug(n+1) = ug(n-2);

  if (use_WENO) {

    for (int gg = 2; gg <= NGROUPS+1; gg++) {
      fg(gg) = ag(gg) * ug(gg);
      ag(gg) = std::abs(ag(gg));
    }

    f(0) = 0.0_rt;

    for (int i = 1; i < n; i++) {
      Real alpha = amrex::max(ag(i-3), ag(i-2), ag(i-1),
                              ag(i), ag(i+1), ag(i+2));
      Real fp[5];
      Real fm[5];
      for (int m = 0; m < 5; m++) {
        fp[m] = 0.5e0_rt * (fg(m-3) + alpha * ug(m-3));
        fm[m] = 0.5e0_rt * (fg(m-2) - alpha * ug(m-2));
      }
      Real fpw = weno5(fp[0], fp[1], fp[2], fp[3], fp[4]);
      Real fmw = weno5(fm[4], fm[3], fm[2], fm[1], fm[0]);
      f(i) = fpw + fmw;
    }
    f(n) = 0.0_rt;

  } else {

    f(0) = 0.0_rt;

    for (int i = 1; i < n; i++) {
      Real r = (ug(i-1) - ug(i-2)) / (ug(i) - ug(i-1) + 1.e-50_rt);
      Real ul = ug(i-1) + 0.5_rt * (ug(i) - ug(i-1)) * MC(r);

      r = (ag(i-1) - ag(i-2)) / (ag(i) - ag(i-1) + 1.e-50_rt);
      Real al = ag(i-1) + 0.5_rt * (ag(i) - ag(i-1)) * MC(r);

      Real fl = al*ul;

      r = (ug(i) - ug(i-1)) / (ug(i+1) - ug(i) + 1.e-50_rt);
      Real ur = ug(i) - 0.5_rt * (ug(i+1) - ug(i)) * MC(r);

      r = (ag(i) - ag(i-1)) / (ag(i+1) - ag(i) + 1.e-50_rt);
      Real ar = ag(i) - 0.5_rt * (ag(i+1) - ag(i)) * MC(r);

      Real fr = ar*ur;

      Real a_plus = amrex::max(0.0_rt, al, ar);
      Real a_minus = amrex::max(0.e0_rt, -al, -ar);
      f(i) = (a_plus*fl + a_minus*fr - a_plus*a_minus*(ur-ul)) /
        (a_plus + a_minus + 1.e-50_rt);
    }
    f(n) = 0.0_rt;

  }

  for (int i = 0; i < n; i++) {
    dudt_tmp[i] = (f(i) - f(i+1)) / dx[i];
  }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
update_one_species(const int n, Real* u, Real* a, const Real* dx, const Real tend, int& nstepmax) {

  Real dt = 1.e50_rt;
  Real acfl;

  int nstep;

  for (int i = 0; i < n; i++) {
    acfl = 1.e-50_rt + std::abs(a[i]);
    dt = amrex::min(dt, dx[i]/acfl*cfl_fspace);
  }

  if (dt >= tend) {
    nstep = 1;
    dt = tend;
  } else {
    nstep = std::ceil(tend/dt);
    dt = tend / static_cast<Real>(nstep);
  }

  Real u1[NGROUPS];
  Real u2[NGROUPS];
  Real u3[NGROUPS];
  Real u4[NGROUPS];
  Real u5[NGROUPS];
  Real dudt_tmp[NGROUPS];

  for (int istep = 1; istep <= nstep; istep++) {
    if (rk_order == 5) {
      // RK5
      dudt(u, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u1[g] = u[g] + B1 * dt * dudt_tmp[g];
      }
      dudt(u1, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u2[g] = (C20*u[g] + C21*u1[g]) + B2 * dt * dudt_tmp[g];
      }
      dudt(u2, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u3[g] = u[g] + B3 * dt * dudt_tmp[g];
      }
      dudt(u3, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u4[g] = (C40*u[g] + C41*u1[g] + C42*u2[g] + C43*u3[g]) + B4 * dt * dudt_tmp[g];
      }
      dudt(u4, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u5[g] = (C50*u[g] + C51*u1[g] + C52*u2[g] + C53*u3[g] + C54*u4[g]) + B5 * dt * dudt_tmp[g];
      }
      dudt(u5, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u[g] = (C60*u[g] + C61*u1[g] + C62*u2[g] + C63*u3[g] + C64*u4[g] + C65*u5[g]) + B6 * dt * dudt_tmp[g];
      }

    } else if (rk_order == 4) {
      // RK4
      dudt(u, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u1[g] = u[g] + 0.5_rt * dt * dudt_tmp[g];
      }
      dudt(u1, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u2[g] = u[g] + 0.5_rt * dt * dudt_tmp[g];
      }
      dudt(u2, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u3[g] = u[g] + dt * dudt_tmp[g];
      }
      dudt(u3, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u[g] = onethird * (u1[g] + 2.0_rt * u2[g] + u3[g] - u[g]) + onesixth * dt * dudt_tmp[g];
      }

    } else if (rk_order == 3) {
      // RK3
      dudt(u, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u1[g] = u[g] + dt * dudt_tmp[g];
      }
      dudt(u1, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u1[g] = 0.75e0_rt * u[g] + 0.25_rt * (u1[g] + dt * dudt_tmp[g]);
      }
      dudt(u1, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u[g] = onethird * u[g] + twothirds * (u1[g] + dt * dudt_tmp[g]);
      }

    } else {
      // first-order
      dudt(u, a, dx, n, dudt_tmp);
      for (int g = 0; g < n; g++) {
        u[g] += dt * dudt_tmp[g];
      }
    }
  }

  nstepmax = amrex::max(nstepmax, nstep);
}

#endif
