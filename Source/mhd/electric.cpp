#include <Castro.H>

#include <mhd_util.H>

using namespace amrex;

void
Castro::electric_edge_x(const Box& bx,
                        Array4<Real const> const& q_arr,
                        Array4<Real> const& E,
                        Array4<Real const> const& flxy,
                        Array4<Real const> const& flxz) {

  // Compute Ex on an edge.  This will compute Ex(i, j-1/2, k-1/2)

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    Real q_zone[NQ];

    // Compute Ex(i, j-1/2, k-1/2) using MM Eq. 50

    // dEx/dy (Eq. 49), located at (i, j-3/4, k-1/2)

    // first compute dEx/dy_{i,j-3/4,k-1} using MM Eq. 49
    // note that the face value Ex_{i,j-1/2,k-1} = -F_{i,j-1/2,k-1}(Bz)
    // via Faraday's law (MM Eq. 15)

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j-1,k-1,n);
    }
    Real Ecen = 0.0_rt;
    electric(q_zone, Ecen, 0);
    Real a = 2.0_rt * (-flxy(i,j,k-1,UMAGZ) - Ecen);

    // now compute dEx/dy_{i,j-3/4,k}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j-1,k,n);
    }
    electric(q_zone, Ecen, 0);
    Real b = 2.0_rt *(-flxy(i,j,k,UMAGZ) - Ecen);

    // Upwind in the z direction to get dEx/dy i, j-3/4, k-1/2
    // using w_{i,j-1,k-1/2}
    // recall flxz(QRHO) = rho*w so sign(rho*w) = sign(w)

    Real d1 = 0.0;
    if (flxz(i,j-1,k,URHO) > 0.0_rt) {
      d1 = a;
    } else if (flxz(i,j-1,k,URHO) < 0.0_rt) {
      d1 = b;
    } else {
      d1 = 0.5_rt * (a + b);
    }

    // dEx/dy located at (i, j-1/4, k-1/2)

    // first compute dEx/dy_{i,j-1/4,k-1}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j,k-1,n);
    }
    electric(q_zone, Ecen, 0);
    a = 2.0_rt * (Ecen + flxy(i,j,k-1,UMAGZ));

    // now compute dEx/dy_{i,j-1/4,k}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j,k,n);
    }
    electric(q_zone, Ecen, 0);
    b = 2.0_rt * (Ecen + flxy(i,j,k,UMAGZ));

    // finally upwind in the z direction to get dEx/dy i, j-1/4, k-1/2
    // using w_{i,j,k-1/2}

    Real d2 = 0.0;
    if (flxz(i,j,k,URHO) > 0.0_rt) {
      d2 = a;
    } else if (flxz(i,j,k,URHO) < 0.0_rt) {
      d2 = b;
    } else {
      d2 = 0.5_rt * (a + b);
    }

    // Calculate the "second derivative" in the y direction for
    // d^2Ex/dy^2 i, j-1/2, k-1/2 (this is one of the terms in Eq. 50)
    // note: Stone 08 Eq. 79 has the signs backwards for this term.

    Real dd1 = 0.125_rt * (d1 - d2);


    // now dEx/dz located at (i, j-1/2, k-3/4)

    // first compute dEx/dz_{i,j-1,k-3/4}
    // note that the face value of Ex_{i,j-1,k-1/2} = F_{i,j-1,k-1/2}(By)

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j-1,k-1,n);
    }
    electric(q_zone, Ecen, 0);
    a = 2.0_rt * (flxz(i,j-1,k,UMAGY) - Ecen);

    // now compute dEx/dz_{i,j,k-3/4}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j,k-1,n);
    }
    electric(q_zone, Ecen, 0);
    b = 2.0_rt * (flxz(i,j,k,UMAGY) - Ecen);

    // upwind in the y direction to get dEx/dz i, j-1/2, k-3/4
    // using v_{i,j-1/2,k-1}

    if (flxy(i,j,k-1,URHO) > 0.0_rt) {
      d1 = a;
    } else if (flxy(i,j,k-1,URHO) < 0.0_rt) {
      d1 = b;
    } else {
      d1 = 0.5_rt * (a + b);
    }

    // dEx/dz located at (i, j-1/2, k-1/4)

    // first compute dEx/dz_{i,j-1,k-1/4}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j-1,k,n);
    }
    electric(q_zone, Ecen, 0);
    a = 2.0_rt * (Ecen - flxz(i,j-1,k,UMAGY));

    // now compute dEx/dz_{i,j,k-1/4}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j,k,n);
    }
    electric(q_zone, Ecen, 0);
    b = 2.0_rt * (Ecen - flxz(i,j,k,UMAGY));

    // upwind in the y direction to get dEx/dz i, j-1/2, k-1/4
    // using v_{i,j-1/2,k}

    if (flxy(i,j,k,URHO) > 0.0_rt) {
      d2 = a;
    } else if (flxy(i,j,k,URHO) < 0.0_rt) {
      d2 = b;
    } else {
      d2 = 0.5_rt * (a + b);
    }

    // calculate second derivative

    Real dd2 = 0.125_rt * (d1 - d2);

    // now the final Ex_{i,j-1/2,k-1/2}, using MM Eq. 50 (shifted to j-1/2, k-1/2)

    E(i,j,k) = 0.25_rt * (-flxy(i,j,k,UMAGZ) - flxy(i,j,k-1,UMAGZ) +
                          flxz(i,j-1,k,UMAGY) + flxz(i,j,k,UMAGY)) + dd1 + dd2;

  });
}

void
Castro::electric_edge_y(const Box& bx,
                        Array4<Real const> const& q_arr,
                        Array4<Real> const& E,
                        Array4<Real const> const& flxx,
                        Array4<Real const> const& flxz) {

  // Compute Ey on an edge.  This will compute Ey(i-1/2, j, k-1/2)

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    Real q_zone[NQ];

    // Compute Ey(i-1/2, j, k-1/2)

    // dEy/dz i-1/2, j, k-3/4

    // first compute dEy/dz_{i-1,j,k-3/4}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i-1,j,k-1,n);
    }
    Real Ecen = 0.0_rt;
    electric(q_zone, Ecen, 1);
    Real a = 2.0_rt * (-flxz(i-1,j,k,UMAGX) - Ecen);

    // now compute dEy/dz_{i,j,k-3/4}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j,k-1,n);
    }
    electric(q_zone, Ecen, 1);
    Real b = 2.0_rt * (-flxz(i,j,k,UMAGX) - Ecen);

    // upwind in the x direction to get dEy/dz i-1/2, j, k-3/4
    // using u_{i-1/2,j,k-1}

    Real d1 = 0.0;
    if (flxx(i,j,k-1,URHO) > 0.0_rt) {
      d1 = a;
    } else if (flxx(i,j,k-1,URHO) < 0.0_rt) {
      d1 = b;
    } else {
      d1 = 0.5_rt * (a + b);
    }

    // dEy/dz i-1/2, j, k-1/4

    // first compute dEy/dz_{i-1,j,k-1/4}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i-1,j,k,n);
    }
    electric(q_zone, Ecen, 1);
    a = 2.0_rt * (Ecen + flxz(i-1,j,k,UMAGX));

    // now compute dEy/dz_{i,j,k-1/4}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j,k,n);
    }
    electric(q_zone, Ecen, 1);
    b = 2.0_rt * (Ecen + flxz(i,j,k,UMAGX));

    // upwind in the x direction to get dEy/dz i-1/2, j, k-1/4
    // using u_{i-1/2.j,k}

    Real d2 = 0.0;
    if (flxx(i,j,k,URHO) > 0.0_rt) {
      d2 = a;
    } else if (flxx(i,j,k,URHO) < 0.0_rt) {
      d2 = b;
    } else {
      d2 = 0.5_rt * (a + b);
    }

    // calculate the "second derivative" in the y direction for
    // d^2Ey/dz^2 i-1/2, j, k-1/2

    Real dd1 = 0.125_rt * (d1 - d2);


    // dEy/dx i-3/4, j, k-1/2

    // first compute dEy/dz_{i-3/4,j,k-1}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i-1,j,k-1,n);
    }
    electric(q_zone, Ecen, 1);
    a = 2.0_rt * (flxx(i,j,k-1,UMAGZ) - Ecen);

    // next compute dEy/dz_{i-3/4,j,k}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i-1,j,k,n);
    }
    electric(q_zone, Ecen, 1);
    b = 2.0_rt * (flxx(i,j,k,UMAGZ) - Ecen);

    // upwind in the z direction to get dEy/dx i-3/4, j, k-1/2
    // using w_{i-1,j,k-1/2}

    if (flxz(i-1,j,k,URHO) > 0.0_rt) {
      d1 = a;
    } else if (flxz(i-1,j,k,URHO) < 0.0_rt) {
      d1 = b;
    } else {
      d1 = 0.5_rt * (a + b);
    }

    // dEy/dx i-1/4, j, k-1/2

    // first compute dEy/dx_{i-1/4,j,k-1}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j,k-1,n);
    }
    electric(q_zone, Ecen, 1);
    a = 2.0_rt * (Ecen - flxx(i,j,k-1,UMAGZ));

    // next compute dEy/dx_{i-1/4,j,k}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j,k,n);
    }
    electric(q_zone, Ecen, 1);
    b = 2.0_rt * (Ecen - flxx(i,j,k,UMAGZ));

    // upwind in the z direction for i-1/4, j, k-1/2
    // using w_{i,j,k-1/2}

    if (flxz(i,j,k,URHO) > 0.0_rt) {
      d2 = a;
    } else if (flxz(i,j,k,URHO) < 0.0_rt) {
      d2 = b;
    } else {
      d2 = 0.5_rt * (a + b);
    }

    // calculate second derivative

    Real dd2 = 0.125_rt * (d1 - d2);

    // now the final Ey_{i-1/2, j, k-1/2}

    E(i,j,k) = 0.25_rt * (-flxz(i,j,k,UMAGX) - flxz(i-1,j,k,UMAGX) +
                          flxx(i,j,k-1,UMAGZ) + flxx(i,j,k,UMAGZ)) + dd1 + dd2;

  });
}


void
Castro::electric_edge_z(const Box& bx,
                        Array4<Real const> const& q_arr,
                        Array4<Real> const& E,
                        Array4<Real const> const& flxx,
                        Array4<Real const> const& flxy) {

  // Compute Ez on an edge.  This will compute Ez(i-1/2, j-1/2, k)

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    Real q_zone[NQ];

    // Compute Ez(i-1/2, j-1/2, k)

    // dEz/dx i-3/4, j-1/2, k

    // first compute dEz/dx_{i-3/4,j-1,k}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i-1,j-1,k,n);
    }
    Real Ecen = 0.0_rt;
    electric(q_zone, Ecen, 2);
    Real a = 2.0_rt * (-flxx(i,j-1,k,UMAGY) - Ecen);

    //  next dEz/dx_{i-3/4,j,k}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i-1,j,k,n);
    }
    electric(q_zone, Ecen, 2);
    Real b = 2.0_rt * (-flxx(i,j,k,UMAGY) - Ecen);

    // upwind in the y direction to get dEz/dx i-3/4, j-1/2, k
    // using v_{i-1,j-1/2,k}

    Real d1 = 0.0_rt;
    if ( flxy(i-1,j,k,URHO) > 0.0_rt) {
      d1 = a;
    } else if (flxy(i-1,j,k,URHO) < 0.0_rt) {
      d1 = b;
    } else {
      d1 = 0.50_rt * (a + b);
    }

    // dEz/dx i-1/4, j-1/2, k

    // first compute dEz/dx_{i-1/4,j-1,k}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j-1,k,n);
    }
    electric(q_zone, Ecen, 2);
    a = 2.0_rt * (Ecen + flxx(i,j-1,k,UMAGY));

    // next dEz/dx_{i-1/4,j,k}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j,k,n);
    }
    electric(q_zone, Ecen, 2);
    b = 2.0_rt * (Ecen + flxx(i,j,k,UMAGY));

    // upwind in the y direction to get dEz/dx i-1/4, j-1/2, k
    // using v_{i,j-1/2,k}

    Real d2 = 0.0_rt;
    if (flxy(i,j,k,URHO) > 0.0_rt) {
      d2 = a;
    } else if (flxy(i,j,k,URHO) < 0.0_rt) {
      d2 = b;
    } else {
      d2 = 0.5_rt * (a + b);
    }

    // Calculate the "second derivative" in the x direction for
    // d^2Ez/dx^2 i-1/2, j-1/2, k

    Real dd1 = 0.125_rt * (d1 - d2);


    // dEz/dy i-1/2, j-3/4, k

    // first compute dEz/dy_{i-1,j-3/4,k}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i-1,j-1,k,n);
    }
    electric(q_zone, Ecen, 2);
    a = 2.0_rt * (flxy(i-1,j,k,UMAGX) - Ecen);

    // now compute dEz/dy_{i,j-3/4,k}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j-1,k,n);
    }
    electric(q_zone, Ecen, 2);
    b = 2.0_rt * (flxy(i,j,k,UMAGX) - Ecen);

    // upwind in the x direction to get dEz/dy i-1/2, j-3/4, k
    // using u_{i-1/2,j-1,k}

    if (flxx(i,j-1,k,URHO) > 0.0_rt) {
      d1 = a;
    } else if (flxx(i,j-1,k,URHO) < 0.0_rt) {
      d1 = b;
    } else {
      d1 = 0.5_rt * (a + b);
    }

    // dEz/dy i-1/2, j-1/4, k

    // first compute dEz/dy_{i-1,j-1/4,k}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i-1,j,k,n);
    }
    electric(q_zone, Ecen, 2);
    a = 2.0_rt * (Ecen - flxy(i-1,j,k,UMAGX));

    // now compute dEz/dy_{i,j-1/4,k}

    for (int n = 0; n < NQ; n++) {
      q_zone[n] = q_arr(i,j,k,n);
    }
    electric(q_zone, Ecen, 2);
    b = 2.0_rt * (Ecen - flxy(i,j,k,UMAGX));

    // Upwind in the x direction for i-1/2, j-1/4, k
    // using u_{i-1/2,j,k}

    if (flxx(i,j,k,URHO) > 0.0_rt) {
      d2 = a;
    } else if (flxx(i,j,k,URHO) < 0.0_rt) {
      d2 = b;
    } else {
      d2 = 0.5_rt * (a + b);
    }

    // calculate second derivative

    Real dd2 = 0.125_rt * (d1 - d2);

    // compute Ez i-1/2, j-1/2, k

    E(i,j,k) = 0.25_rt * (-flxx(i,j,k,UMAGY) - flxx(i,j-1,k,UMAGY) +
                          flxy(i-1,j,k,UMAGX) + flxy(i,j,k,UMAGX)) + dd1 + dd2;

  });
}

