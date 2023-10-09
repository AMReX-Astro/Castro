#include <Castro.H>
#include <Castro_util.H>

#include <mhd_util.H>

using namespace amrex;

void
Castro::corner_couple(const Box& bx,
                      Array4<Real> const& qr_out,
                      Array4<Real> const& ql_out,
                      Array4<Real const> const& ur,
                      Array4<Real const> const& ul,
                      Array4<Real const> const& flxd2,
                      Array4<Real const> const& Ed1,
                      Array4<Real const> const& Ed3,
                      const int d1, const int d2, const int d3,
                      const Real dt) {

  // take conservative interface states ul and ur and update them
  // with with the transverse flux difference (corner coupling) to
  // produce ul_out and ur_out
  //
  // This implements MM step 3 of the CTU algorithm.
  // the normal direction (for the interface states) is d1
  // the transverse direction (for the flux difference) is d2

  GpuArray<Real, 3> dx;
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      dx[i] = geom.CellSizeArray()[i];
  }
  for (int i = AMREX_SPACEDIM; i < 3; ++i) {
      dx[i] = 0.0_rt;
  }

  // cl and cr are the offsets to the indices for the conserved state fluxes
  // they will be offset in d2 to capture the flux difference

  int cl[3] = {};
  int cr[3] = {};

  // b is for indexing into the magnetic field that is in the d1
  // (normal) direction

  int b[3] = {};

  // for indexing the electric field

  int err[3] = {};
  int elr[3] = {};
  int erl[3] = {};
  int ell[3] = {};

  // update the state on interface direction d1 with the input flux in direction d2

  // for the flux difference, F_r - F_l, we need to shift the indices in the first flux (F_r)
  // in d2 to get a difference across the interface.  We also need to shift by a zone in d1
  // for the left interface.  cr(:) and cl(:) will hold these shifts.

  // the first term of the flxd2 subtraction is shifted by 1 on the direction d2
  cr[d2] = 1;

  // for the normal B component
  b[d2] = 1;

  // err will capture the right state in both transverse directions
  // (e.g. Ez_{i+1/2,j+1/2,k})
  err[d2] = 1;
  err[d3] = 1;

  // elr will capture the right state in the second transverse direction
  // (e.g. Ez_{i-1/2,j+1/2,k})
  elr[d3] = 1;

  // erl will capture the right state in the first transverse direction
  // (e.g. Ez_{i+1/2,j-1/2,k})
  erl[d2] = 1;

  // ell is the lower-left E in both directions
  // (e.g. Ez_{i-1/2,j-1/2,k})

  int sgn = epsilon_ijk(d1, d2, d3);
  Real cdtdx = dt/(3.0_rt * dx[d1]);

  // UMAGD1 corresponds to d1, and UMAGD2 to d2, UMAGD3 to d3
  int UMAGD1 = UMAGX + d1;
  int UMAGD2 = UMAGX + d2;
  int UMAGD3 = UMAGX + d3;

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    // first the conserved state

    // right interface (e.g. U_{i-1/2,j,k,R} or the "-" state in MM notation)
    // MM Eq. 37

    Real utmp[NUM_STATE+3];
    Real qtmp[NQ];

    for (int n = 0; n < NUM_STATE; n++) {
      if (n != UTEMP) {
        utmp[n] = ur(i,j,k,n) -
          cdtdx*(flxd2(i+cr[0],j+cr[1],k+cr[2],n) - flxd2(i,j,k,n));

      } else {
        utmp[n] = 0.0_rt;
      }
    }

    // now magnetic fields

    // right state on the interface (e.g. B_{i-1/2,j,k,R} or `-` in MM notation)

    // d1 -- this is perpendicular to the face, MM Eq. 38
    // (note MM Eq. 38 has a sign error) e.g., for d1 = x and
    // d2 = y, this gets updated as
    //
    // Bx|y_{i-1/2,j,k,R} = Bx_{i-1/2,j,k) -
    //     dt/3dx (Ez_{i-1/2,j+1/2,k) - Ez_{i-1/2,j-1/2,k})
    //
    // we use b[] to captured the j+1/2 indexing into Ez

    utmp[UMAGD1] = ur(i,j,k,UMAGD1) - sgn * cdtdx *
      (Ed3(i+b[0],j+b[1],k+b[2]) - Ed3(i,j,k));

    // d3 -- this is in the plane of the face, MM Eq. 39
    // e.g.,g for d1 = x, and d2 = y, this gets updated as
    //
    // Bz|y_{i-1/2,j,k,R} = Bz_{i-1/2,j,k} + 1/2 dt/3dx
    //     (Ex_{i,j+1/2,k+1/2} - Ex_{i,j-1/2,k+1/2} +
    //      Ex_{i,j+1/2,k-1/2} - Ex_{i,j-1/2,k-1/2})
    //
    // we use err(:) for the first E term, elr(:) for the
    // second, and erl(:) for the third

    utmp[UMAGD3] = ur(i,j,k,UMAGD3) + sgn * 0.5_rt * cdtdx *
      ((Ed1(i+err[0],j+err[1],k+err[2]) - Ed1(i+elr[0],j+elr[1],k+elr[2])) +
       (Ed1(i+erl[0],j+erl[1],k+erl[2]) - Ed1(i,j,k)));

    // the component pointing in the transverse update direction, d2, is unchanged
    utmp[UMAGD2] = ur(i,j,k,UMAGD2);

    ConsToPrim(qtmp, utmp);

    for (int n = 0; n < NQ; n++) {
      qr_out(i,j,k,n) = qtmp[n];
    }
  });


  // left interface (e.g., U_{i-1/2,j,k,L} or the "+" state in MM notation)
  // note: this uses information one zone to the left in d1

  cl[d1] -= 1;
  cr[d1] -= 1;

  // The in-plane B component at B_{i-1/2,j,k,L} uses the information one zone to the left
  // in direction d1

  err[d1] -= 1;
  elr[d1] -= 1;
  erl[d1] -= 1;
  ell[d1] -= 1;

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    // conservative state

    Real utmp[NUM_STATE+3];
    Real qtmp[NQ];

    for (int n = 0; n < NUM_STATE; n++) {
      if (n != UTEMP) {
        utmp[n] = ul(i,j,k,n) -
          cdtdx*(flxd2(i+cr[0],j+cr[1],k+cr[2],n) -
                 flxd2(i+cl[0],j+cl[1],k+cl[2],n));

      } else {
        utmp[n] = 0.0_rt;
      }
    }

    // left state on the interface (e.g. B_{i-1/2,j,k,L} or `+` in MM notation)

    utmp[UMAGD1] = ul(i,j,k,UMAGD1) - sgn * cdtdx *
      (Ed3(i+b[0],j+b[1],k+b[2]) - Ed3(i,j,k));

    utmp[UMAGD3] = ul(i,j,k,UMAGD3) + sgn * 0.5_rt * cdtdx *
      ((Ed1(i+err[0],j+err[1],k+err[2]) - Ed1(i+elr[0],j+elr[1],k+elr[2])) +
       (Ed1(i+erl[0],j+erl[1],k+erl[2]) - Ed1(i+ell[0],j+ell[1],k+ell[2])));

    utmp[UMAGD2] = ul(i,j,k,UMAGD2);

    ConsToPrim(qtmp, utmp);

    for (int n = 0; n < NQ; n++) {
      ql_out(i,j,k,n) = qtmp[n];
    }
  });
}


void
Castro::half_step(const Box& bx,
                  Array4<Real> const& qr_out,
                  Array4<Real> const& ql_out,
                  Array4<Real const> const& ur,
                  Array4<Real const> const& ul,
                  Array4<Real const> const& flxd1,
                  Array4<Real const> const& flxd2,
                  Array4<Real const> const& Ed,
                  Array4<Real const> const& Ed1,
                  Array4<Real const> const& Ed2,
                  const int d, const int d1, const int d2, const Real dt) {

  // Final transverse flux corrections to the conservative state

  GpuArray<Real, 3> dx;
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      dx[i] = geom.CellSizeArray()[i];
  }
  for (int i = AMREX_SPACEDIM; i < 3; ++i) {
      dx[i] = 0.0_rt;
  }

  // c1l, c1r are for indexing flxd1 offsets, c2l, c2r are for flxd2

  int c1l[3] = {};
  int c1r[3] = {};
  int c2l[3] = {};
  int c2r[3] = {};

  int a1[3] = {};
  int a2[3] = {};

  int err[3] = {};
  int erl[3] = {};
  int elr[3] = {};
  int ell[3] = {};

  int e1rr[3] = {};
  int e1rl[3] = {};
  int e1lr[3] = {};

  int e2rr[3] = {};
  int e2rl[3] = {};
  int e2lr[3] = {};

  Real hdtdx = 0.5_rt * dt / dx[d];
  int sgn = -1 * epsilon_ijk(d, d1, d2);

  // UMAGD1 corresponds to d1, and UMAGD2 to d2, UMAGD3 to d3
  int UMAGD = UMAGX + d;
  int UMAGD1 = UMAGX + d1;
  int UMAGD2 = UMAGX + d2;


  c1r[d1] = 1;  // add +1 to the d1 direction in the first flxd1 term of the subtraction
  c2r[d2] = 1;  // add +1 to the d2 direction in the first flxd2 term of the subtraction

  // for Ed
  err[d1] = 1;
  err[d2] = 1;

  erl[d1] = 1;

  elr[d2] = 1;

  // for Ed1
  e1rr[d] = 1;
  e1rr[d2] = 1;

  e1lr[d2] = 1;

  e1rl[d] = 1;

  // for Ed2
  e2rr[d] = 1;
  e2rr[d1] = 1;

  e2lr[d1] = 1;

  e2rl[d] = 1;

  // for the normal component of B
  a1[d2] = 1;  // shift on first term of Ed1 subtraction, in d2 direction
  a2[d1] = 1;  // shift on first term of Ed2 subtraction, in d1 direction


  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    // first the conservative state

    // right interface (e.g. U_{i-1/2,j,k,R} or the "-" state in MM notation)
    // MM Eq. 44

    Real utmp[NUM_STATE+3];
    Real qtmp[NQ];

    for (int n = 0; n < NUM_STATE; n++) {
      if (n != UTEMP) {
        utmp[n] = ur(i,j,k,n) -
          hdtdx * (flxd1(i+c1r[0],j+c1r[1],k+c1r[2],n) - flxd1(i,j,k,n)) -
          hdtdx * (flxd2(i+c2r[0],j+c2r[1],k+c2r[2],n) - flxd2(i,j,k,n));

      } else {
        utmp[n] = 0.0_rt;
      }
    }

    // now the magnetic field

    // right state on the interface (e.g. B_{i-1/2,j,k,R} or `-` in MM notation)

    // Bd -- this is perpendicular to the face. Note MM eq.45
    // in Miniati has a sign error in the epsilon term

    utmp[UMAGD] = ur(i,j,k,UMAGD) - sgn * hdtdx *
      ((Ed1(i+a1[0],j+a1[1],k+a1[2]) - Ed1(i,j,k)) -
       (Ed2(i+a2[0],j+a2[1],k+a2[2]) - Ed2(i,j,k)));

    // Bd1 -- this is one of the components of B in the plane of the face d
    // Eq.46 in Miniati

    utmp[UMAGD1] = ur(i,j,k,UMAGD1) + sgn * 0.5_rt * hdtdx *
      ((Ed(i+err[0],j+err[1],k+err[2]) - Ed(i+erl[0],j+erl[1],k+erl[2])) +
       (Ed(i+elr[0],j+elr[1],k+elr[2]) - Ed(i+ell[0],j+ell[1],k+ell[2])) -
       (Ed2(i+e2rr[0],j+e2rr[1],k+e2rr[2]) - Ed2(i+e2lr[0],j+e2lr[1],k+e2lr[2])) -
       (Ed2(i+e2rl[0],j+e2rl[1],k+e2rl[2]) - Ed2(i+ell[0],j+ell[1],k+ell[2])));

    // Bd2 -- this is the other component of B in the plane of the face d
    // Eq. 46 in Miniati

    utmp[UMAGD2] = ur(i,j,k,UMAGD2) - sgn * 0.5_rt * hdtdx *
      ((Ed(i+err[0],j+err[1],k+err[2]) - Ed(i+elr[0],j+elr[1],k+elr[2])) +
       (Ed(i+erl[0],j+erl[1],k+erl[2]) - Ed(i+ell[0],j+ell[1],k+ell[2])) -
       (Ed1(i+e1rr[0],j+e1rr[1],k+e1rr[2]) - Ed1(i+e1lr[0],j+e1lr[1],k+e1lr[2])) -
       (Ed1(i+e1rl[0],j+e1rl[1],k+e1rl[2]) - Ed1(i+ell[0],j+ell[1],k+ell[2])));

    // convert to primitive

    ConsToPrim(qtmp, utmp);

    for (int n = 0; n < NQ; n++) {
      qr_out(i,j,k,n) = qtmp[n];
    }
  });


  // for the left state U components on the face dir, the flux
  // difference is in the zone to the left

    c1r[d] -= 1;
    c1l[d] -= 1;
    c2r[d] -= 1;
    c2l[d] -= 1;

    // left state on the interface (e.g., B_{i-1/2,j,k,L} or `+` in MM notation)

    // The in-plane B component at B_{i-1/2,j,k,L} uses the information one zone to the left
    // in direction d1

    err[d] -= 1;
    erl[d] -= 1;
    elr[d] -= 1;
    ell[d] -= 1;

    e1rr[d] -= 1;
    e1rl[d] -= 1;
    e1lr[d] -= 1;

    e2rr[d] -= 1;
    e2rl[d] -= 1;
    e2lr[d] -= 1;


  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    // left interface (e.g., U_{i+1/2,j,k,L} or the "+" state in MM notation)

    Real utmp[NUM_STATE+3];
    Real qtmp[NQ];

    // first the conserved state

    for (int n = 0; n < NUM_STATE; n++) {
      if (n != UTEMP) {
        utmp[n] = ul(i,j,k,n) -
          hdtdx * (flxd1(i+c1r[0],j+c1r[1],k+c1r[2],n) -
                   flxd1(i+c1l[0],j+c1l[1],k+c1l[2],n)) -
          hdtdx * (flxd2(i+c2r[0],j+c2r[1],k+c2r[2],n) -
                   flxd2(i+c2l[0],j+c2l[1],k+c2l[2],n));

      } else {
        utmp[n] = 0.0_rt;
      }
    }

    // now the B fields

    // Bd -- this is perpendicular to the face (MM Eq. 45 with sign fix)

    // this is the same face as the right state, so the update the identical
    utmp[UMAGD] = ul(i,j,k,UMAGD) - sgn * hdtdx *
      ((Ed1(i+a1[0],j+a1[1],k+a1[2]) - Ed1(i,j,k)) -
       (Ed2(i+a2[0],j+a2[1],k+a2[2]) - Ed2(i,j,k)));

    // Bd1 -- first component on face d, eq. 46 in Miniati

    utmp[UMAGD1] = ul(i,j,k,UMAGD1) + sgn * 0.5_rt * hdtdx *
      ((Ed(i+err[0],j+err[1],k+err[2]) - Ed(i+erl[0],j+erl[1],k+erl[2])) +
       (Ed(i+elr[0],j+elr[1],k+elr[2]) - Ed(i+ell[0],j+ell[1],k+ell[2])) -
       (Ed2(i+e2rr[0],j+e2rr[1],k+e2rr[2]) - Ed2(i+e2lr[0],j+e2lr[1],k+e2lr[2])) -
       (Ed2(i+e2rl[0],j+e2rl[1],k+e2rl[2]) - Ed2(i+ell[0],j+ell[1],k+ell[2])));

    // Bd2 -- second component on face d, eq. 46 in Miniati

    utmp[UMAGD2] = ul(i,j,k,UMAGD2) - sgn * 0.5_rt * hdtdx *
      ((Ed(i+err[0],j+err[1],k+err[2]) - Ed(i+elr[0],j+elr[1],k+elr[2])) +
       (Ed(i+erl[0],j+erl[1],k+erl[2]) - Ed(i+ell[0],j+ell[1],k+ell[2])) -
       (Ed1(i+e1rr[0],j+e1rr[1],k+e1rr[2]) - Ed1(i+e1lr[0],j+e1lr[1],k+e1lr[2])) -
       (Ed1(i+e1rl[0],j+e1rl[1],k+e1rl[2]) - Ed1(i+ell[0],j+ell[1],k+ell[2])));

    // convert to primitive

    ConsToPrim(qtmp, utmp);

    for (int n = 0; n < NQ; n++) {
      ql_out(i,j,k,n) = qtmp[n];
    }
  });

}
