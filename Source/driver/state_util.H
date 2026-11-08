#ifndef CASTRO_STATE_UTIL_H
#define CASTRO_STATE_UTIL_H

#include <AMReX_Array4.H>
#include <Castro_util.H>
#include <castro_params.H>
#include <eos.H>

#ifdef HYBRID_MOMENTUM
#include <prob_parameters.H>
#include <hybrid.H>
#endif

namespace state_util {

// Shared cell operations for the Fab and MultiFab cleanup interfaces.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void reset_internal_energy(int i, int j, int k,
#ifdef MHD
                           amrex::Array4<amrex::Real const> const& Bx,
                           amrex::Array4<amrex::Real const> const& By,
                           amrex::Array4<amrex::Real const> const& Bz,
#endif
                           amrex::Array4<amrex::Real> const& u,
                           amrex::Real lsmall_temp, amrex::Real ldual_energy_eta2)
{

    auto U_cell = u.cellData(i, j, k);

    const Real rhoInv = 1.0_rt / U_cell[URHO];
    const Real Up = U_cell[UMX] * rhoInv;
    const Real Vp = U_cell[UMY] * rhoInv;
    const Real Wp = U_cell[UMZ] * rhoInv;
    const Real ke = 0.5_rt * (Up * Up + Vp * Vp + Wp * Wp);

    eos_re_t eos_state;

    eos_state.rho = U_cell[URHO];
    eos_state.T   = lsmall_temp;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = U_cell[UFS+n] * rhoInv;
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        eos_state.aux[n] = U_cell[UFX+n] * rhoInv;
    }
#endif

    eos(eos_input_rt, eos_state);

    const Real small_e = eos_state.e;

#ifdef MHD
    const Real bx_cell_c = 0.5_rt * (Bx(i,j,k) + Bx(i+1,j,k));
    const Real by_cell_c = 0.5_rt * (By(i,j,k) + By(i,j+1,k));
    const Real bz_cell_c = 0.5_rt * (Bz(i,j,k) + Bz(i,j,k+1));

    const Real B_ener = 0.5_rt * (bx_cell_c*bx_cell_c +
                                  by_cell_c*by_cell_c +
                                  bz_cell_c*bz_cell_c);
#else
    const Real B_ener = 0.0_rt;
#endif

    // Ensure the internal energy is at least as large as this minimum
    // from the EOS; the same holds true for the total energy.

    U_cell[UEINT] = amrex::max(U_cell[UEINT], U_cell[URHO] * small_e);
    U_cell[UEDEN] = amrex::max(U_cell[UEDEN], U_cell[URHO] * (small_e + ke) + B_ener);

    // Apply the dual energy criterion: get e from E if (E - K) > eta * E.

    const Real rho_eint = U_cell[UEDEN] - U_cell[URHO] * ke - B_ener;

    if (rho_eint > ldual_energy_eta2 * U_cell[UEDEN]) {
        U_cell[UEINT] = rho_eint;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void enforce_minimum_density(int i, int j, int k,
                             amrex::Array4<amrex::Real> const& state_arr,
                             int verbose_warnings,
                             amrex::GeometryData const& geomdata)
{
    amrex::ignore_unused(verbose_warnings, geomdata);

    if (state_arr(i,j,k,URHO) < castro::small_dens) {

#ifndef AMREX_USE_GPU
      if (verbose_warnings > 1 ||
          (verbose_warnings > 0 && state_arr(i,j,k,URHO) > castro::retry_small_density_cutoff)) {
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
        std::cout << ">>> FROM " << state_arr(i,j,k,URHO) << " TO " << castro::small_dens << std::endl;
        std::cout << " " << std::endl;
      }
#endif

      for (int ipassive = 0; ipassive < npassive; ipassive++) {
          const int n = upassmap(ipassive);
          state_arr(i,j,k,n) *= (castro::small_dens / state_arr(i,j,k,URHO));
      }

      eos_re_t eos_state;
      eos_state.rho = castro::small_dens;
      eos_state.T = castro::small_temp;
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state_arr(i,j,k,UFS+n) / castro::small_dens;
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = state_arr(i,j,k,UFX+n) / castro::small_dens;
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
}

}

#endif
