#include "Castro.H"
#include "Castro_F.H"
#include "Castro_util.H"
#include "Castro_hydro_F.H"

#include "mhd_sound_speed.H"

using namespace amrex;

void
Castro::check_for_mhd_cfl_violation(const Box& bx,
                                    const Real dt,
                                    Array4<Real const> const& q_arr,
                                    Array4<Real const> const& qaux_arr) {

  auto dx = geom.CellSizeArray();

  Real dtdx = dt / dx[0];
  Real dtdy = dt / dx[1];
  Real dtdz = dt / dx[2];

  ReduceOps<ReduceOpMax> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  reduce_op.eval(bx, reduce_data,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
  {

    //sound speed for ideal mhd
    Real as2 = qaux_arr(i,j,k,QC) * qaux_arr(i,j,k,QC);
    Real rho_inv = 1.0_rt / q_arr(i,j,k,QRHO);

    Real ca = (q_arr(i,j,k,QMAGX) * q_arr(i,j,k,QMAGX) +
               q_arr(i,j,k,QMAGY) * q_arr(i,j,k,QMAGY) +
               q_arr(i,j,k,QMAGZ) * q_arr(i,j,k,QMAGZ)) * rho_inv;

    Real cx = 0.0;
    Real cad = q_arr(i,j,k,QMAGX) * q_arr(i,j,k,QMAGX) * rho_inv;
    eos_soundspeed_mhd(cx, as2, ca, cad);

    Real cy = 0.0;
    cad = q_arr(i,j,k,QMAGY) * q_arr(i,j,k,QMAGY) * rho_inv;
    eos_soundspeed_mhd(cy, as2, ca, cad);

    Real cz = 0.0;
    cad = q_arr(i,j,k,QMAGZ) * q_arr(i,j,k,QMAGZ) * rho_inv;
    eos_soundspeed_mhd(cz, as2, ca, cad);

    Real courx = (cx + std::abs(q_arr(i,j,k,QU))) * dtdx;
    Real coury = (cy + std::abs(q_arr(i,j,k,QV))) * dtdy;
    Real courz = (cz + std::abs(q_arr(i,j,k,QW))) * dtdz;

#ifndef AMREX_USE_CUDA
    if (verbose == 1) {

      if (courx > 1.0_rt) {
        std::cout << std::endl;
        std::cout << "Warning:: CFL violation in check_for_mhd_cfl_violation" << std::endl;
        std::cout << ">>> ... (u+c) * dt / dx > 1 " << courx << std::endl;
        std::cout << ">>> ... at cell i = " << i << " j = " << j << " k = " << k << std::endl;
        std::cout << ">>> ... u = " << q_arr(i,j,k,QU) << " c = " << cx << std::endl;
        std::cout << ">>> ... B = " << q_arr(i,j,k,QMAGX) << " " << q_arr(i,j,k,QMAGY) << " " << q_arr(i,j,k,QMAGZ) << std::endl;
        std::cout << ">>> ... density = " << q_arr(i,j,k,QRHO) << std::endl;
        std::cout << ">>> ... internal e = " << q_arr(i,j,k,QREINT) << std::endl;
        std::cout << ">>> ... pressure = " << q_arr(i,j,k,QPRES) << std::endl;
      }

      if (coury > 1.0_rt) {
        std::cout << std::endl;
        std::cout << "Warning:: CFL violation in check_for_mhd_cfl_violation" << std::endl;
        std::cout << ">>> ... (v+c) * dt / dx > 1 " << coury << std::endl;
        std::cout << ">>> ... at cell i = " << i << " j = " << j << " k = " << k << std::endl;
        std::cout << ">>> ... v = " << q_arr(i,j,k,QV) << " c = " << cy << std::endl;
        std::cout << ">>> ... B = " << q_arr(i,j,k,QMAGX) << " " << q_arr(i,j,k,QMAGY) << " " << q_arr(i,j,k,QMAGZ) << std::endl;
        std::cout << ">>> ... density = " << q_arr(i,j,k,QRHO) << std::endl;
        std::cout << ">>> ... internal e = " << q_arr(i,j,k,QREINT) << std::endl;
        std::cout << ">>> ... pressure = " << q_arr(i,j,k,QPRES) << std::endl;
      }

      if (courz > 1.0_rt) {
        std::cout << std::endl;
        std::cout << "Warning:: CFL violation in check_for_mhd_cfl_violation" << std::endl;
        std::cout << ">>> ... (w+c) * dt / dx > 1 " << courz << std::endl;
        std::cout << ">>> ... at cell i = " << i << " j = " << j << " k = " << k << std::endl;
        std::cout << ">>> ... w = " << q_arr(i,j,k,QW) << " c = " << cz << std::endl;
        std::cout << ">>> ... B = " << q_arr(i,j,k,QMAGX) << " " << q_arr(i,j,k,QMAGY) << " " << q_arr(i,j,k,QMAGZ) << std::endl;
        std::cout << ">>> ... density = " << q_arr(i,j,k,QRHO) << std::endl;
        std::cout << ">>> ... internal e = " << q_arr(i,j,k,QREINT) << std::endl;
        std::cout << ">>> ... pressure = " << q_arr(i,j,k,QPRES) << std::endl;

      }

    }
#endif

    return {amrex::max(courx, coury, courz)};

  });

  ReduceTuple hv = reduce_data.value();
  Real courno = amrex::get<0>(hv);

  if (courno > 1.0) {
    amrex::Print() << "WARNING -- EFFECTIVE CFL AT LEVEL " << level << " IS " << courno << std::endl << std::endl;

    cfl_violation = 1;
  }

}
