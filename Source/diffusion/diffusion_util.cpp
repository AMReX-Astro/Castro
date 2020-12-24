#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>

#include <state_indices.H>

#include <castro_params.H>
#include <eos.H>
#include <conductivity.H>

using namespace amrex;

void
fill_temp_cond(const Box& bx,
               Array4<Real const> const& U_arr,
               Array4<Real> const& coeff_arr) {

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    eos_t eos_state;
    eos_state.rho  = U_arr(i,j,k,URHO);
    Real rhoinv = 1.0_rt/eos_state.rho;

    eos_state.T = U_arr(i,j,k,UTEMP);   // needed as an initial guess
    eos_state.e = U_arr(i,j,k,UEINT) * rhoinv;
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = U_arr(i,j,k,UFS+n) * rhoinv;
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = U_arr(i,j,k,UFX+n) * rhoinv;
    }
#endif

    if (eos_state.e < 0.0_rt) {
      eos_state.T = castro::small_temp;
      eos(eos_input_rt, eos_state);
    } else {
      eos(eos_input_re, eos_state);
    }


    if (eos_state.rho > castro::diffuse_cutoff_density) {
      conductivity(eos_state);

      if (eos_state.rho < castro::diffuse_cutoff_density_hi) {
        Real multiplier = (eos_state.rho - castro::diffuse_cutoff_density) /
          (castro::diffuse_cutoff_density_hi - castro::diffuse_cutoff_density);
        eos_state.conductivity = eos_state.conductivity * multiplier;
      }
    } else {
      eos_state.conductivity = 0.0_rt;
    }
    coeff_arr(i,j,k) = castro::diffuse_cond_scale_fac * eos_state.conductivity;

  });
}


void
fill_temp_diff_coeff(const Box& bx,
                     Array4<Real const> const& U_arr,
                     Array4<Real> const& coeff_arr) {

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    eos_t eos_state;
    eos_state.rho  = U_arr(i,j,k,URHO);
    Real rhoinv = 1.0_rt/eos_state.rho;

    eos_state.T = U_arr(i,j,k,UTEMP);   // needed as an initial guess
    eos_state.e = U_arr(i,j,k,UEINT) * rhoinv;
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = U_arr(i,j,k,UFS+n) * rhoinv;
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = U_arr(i,j,k,UFX+n) * rhoinv;
    }
#endif

    if (eos_state.e < 0.0_rt) {
      eos_state.T = castro::small_temp;
      eos(eos_input_rt, eos_state);
    } else {
      eos(eos_input_re, eos_state);
    }


    if (eos_state.rho > castro::diffuse_cutoff_density) {
      conductivity(eos_state);

      if (eos_state.rho < castro::diffuse_cutoff_density_hi) {
        Real multiplier = (eos_state.rho - castro::diffuse_cutoff_density) /
          (castro::diffuse_cutoff_density_hi - castro::diffuse_cutoff_density);
        eos_state.conductivity = eos_state.conductivity * multiplier;
      }
    } else {
      eos_state.conductivity = 0.0_rt;
    }
    coeff_arr(i,j,k) = castro::diffuse_cond_scale_fac * eos_state.conductivity *
      rhoinv / eos_state.cv;

  });
}

