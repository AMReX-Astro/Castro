#ifndef problem_setup_H
#define problem_setup_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#include <model_parser.H>
#include <initial_model.H>
#include <ambient.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    // check to make sure that small_dens is less than low_density_cutoff
    // if not, funny things can happen above the atmosphere

    if (small_dens >= 0.99_rt * problem::low_density_cutoff) {
        amrex::Error("ERROR: small_dens should be set lower than low_density_cutoff");
    }

    // get the species indices

    bool species_defined = true;

    int ifuel1 = network_spec_index(fuel1_name);
    if (ifuel1 < 0) {
        species_defined = false;
    }

    int ifuel2;
    if (fuel2_name != "") {
        ifuel2 = network_spec_index(fuel2_name);
        if (ifuel2 < 0) {
            species_defined = false;
        }
    }

    int ifuel3;
    if (fuel3_name != "") {
        ifuel3 = network_spec_index(fuel3_name);
        if (ifuel3 < 0) {
            species_defined = false;
        }
    }

    int ifuel4;
    if (fuel4_name != "") {
        ifuel4 = network_spec_index(fuel4_name);
        if (ifuel4 < 0) {
            species_defined = false;
        }
    }

    int iash1 = network_species_index(ash1_name);
    if (iash1 < 0) {
        species_defined = false;
    }

    int iash2;
    if (ash2_name != "") {
        iash2 = network_spec_index(ash2_name);
        if (iash2 < 0) {
            species_defined = false;
        }
    }

    int iash3;
    if (ash3_name != "") {
        iash3 = network_spec_index(ash3_name);
        if (iash3 < 0) {
            species_defined = false;
        }
    }

    if (! species_defined) {
        std::cout << ifuel1 << " " << ifuel2 << " " << ifuel3 << " " << ifuel4 << std::endl;
        std::cout << iash1 << " " << iash2 << " "<< iash3 << std::endl;
        amrex::Error("ERROR: species not defined");
    }

    model_t model_params;

    // set the composition of the underlying star

    for (int n = 0; n < NumSpec; n++) {
        model_params.xn_star[n] = problem::smallx;
    }
    model_params.xn_star[iash1] = ash1_frac;
    if (ash2_name != "") {
        model_params.xn_star[iash2] = ash2_frac;
    }
    if (ash3_name != "") {
        model_params.xn_star[iash3] = ash3_frac;
    }

    // and the composition of the accreted layer

    for (int n = 0; n < NumSpec; n++) {
        model_params.xn_base[n] = problem::smallx;
    }
    model_params.xn_base[ifuel1] = fuel1_frac;
    if (fuel2_name != "") {
        model_params.xn_base[ifuel2] = fuel2_frac;
    }
    if (fuel3_name != "") {
        model_params.xn_base[ifuel3] = fuel3_frac;
    }
    if (fuel4_name != "") {
        model_params.xn_base[ifuel4] = fuel4_frac;
    }

    // check if they sum to 1

    Real sumX = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sumX += model_params.xn_star[n];
    }
    if (std::abs(sumX - 1.0_rt) > nspec * problem::smallx) {
        amrex::Error("ERROR: ash mass fractions don't sum to 1");
    }

    sumX = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sumX += model_params.xn_base[n];
    }
    if (std::abs(sumX - 1.0_rt) > nspec * problem::smallx) {
        amrex::Error("ERROR: fuel mass fractions don't sum to 1");
    }

    // we are going to generate an initial model from problo(2) to
    // probhi(2) with nx_model zones.  But to allow for a interpolated
    // lower boundary, we'll add 4 ghostcells to this, so we need to
    // compute dx

    int nx_model = static_cast<int>((probhi[AMREX_SPACEDIM-1] - problo[AMREX_SPACEDIM-1]) /
                                    problem::dx_model);

    int ng = 4;

    // now generate the initial models

    model_params.dens_base = problem::dens_base;
    model_params.T_star = problem::T_star;
    model_params.T_hi = problem::T_hi;
    model_params.T_lo = problem::T_lo;

    model_params.H_star = problem::H_star;
    model_params.atm_delta = problem::atm_delta;

    model_params.low_density_cutoff = problem::low_density_cutoff;

    generate_initial_model(nx_model + ng,
                           problo[AMREX_SPACEDIM-1] - ng * problem::dx_model,
                           probhi[AMREX_SPACEDIM-1],
                           model_params, 0);

    // now create a perturbed model -- we want the same base conditions
    // a hotter temperature

    model_params.T_hi = model_params.T_hi + problem::dtemp;

    generate_initial_model(nx_model + ng,
                           problo[AMREX_SPACEDIM-1] - ng * problem::dx_model,
                           probhi[AMREX_SPACEDIM-1],
                           model_params, 1);

    // set center

    for (int d = 0; d < AMREX_SPACEDIM; d++) {
        problem::center[d] = 0.5_rt * (problo[d] + probhi[d]);
    }

#if AMREX_SPACEDIM == 2
    // for axisymmetry, put the x-center on the x-axis
    // and the y-center at 0, so the height computation is okay
    problem::center[0] = 0.0_rt;
    problem::center[1] = 0.0_rt;
#endif

    // set the ambient state for the upper boundary condition

    ambient::ambient_state[URHO] = model::profile(0).state(model::npts-1, model::idens);
    ambient::ambient_state[UTEMP] = model::profile(0).state(model::npts-1, model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        ambient::ambient_state[UFS+n] =
            ambient::ambient_state[URHO] * model::profile(0).state(model::npts-1, model::ispec+n);
    }

    ambient::ambient_state[UMX] = 0.0_rt;
    ambient::ambient_state[UMY] = 0.0_rt;
    ambient::ambient_state[UMZ] = 0.0_rt;

    // make the ambient state thermodynamically consistent

    eos_t eos_state;
    eos_state.rho = ambient::ambient_state[URHO];
    eos_state.T = ambient::ambient_state[UTEMP];
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = ambient::ambient_state[UFS+n] / eos_state.rho;
    }

    eos(eos_input_rt, eos_state);

    ambient::ambient_state[UEINT] = eos_state.rho * eos_state.e;
    ambient::ambient_state[UEDEN] = eos_state.rho * eos_state.e;
}


subroutine ca_initdata(lo, hi, &
                       state, s_lo, s_hi, &
                       dx, problo) bind(C, name='ca_initdata')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: ZERO, HALF, ONE
#ifndef AMREX_USE_CUDA
  use castro_error_module, only: castro_error
#endif
  use probdata_module, only: x_half_width, x_half_max
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_rt, eos_input_tp, eos_input_rp
  use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use network, only: nspec
  use initial_model_module, only: gen_npts_model, gen_model_r, gen_model_state, &
                                  idens_model, itemp_model, ipres_model, ispec_model
  use interpolate_module, only: interpolate ! function
  use prob_params_module, only : center

  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
  real(rt), intent(in   ) :: dx(3), problo(3)

  real(rt) :: x, y, z, r, height
  integer :: i, j, k, n

  real(rt) :: temppres

  type (eos_t) :: eos_state
  real(rt) :: f

  !$gpu

  do k = lo(3), hi(3)
     z = problo(3) + (dble(k) + HALF) * dx(3)

     do j = lo(2), hi(2)
        y = problo(2) + (dble(j) + HALF) * dx(2) - center(2)

        do i = lo(1), hi(1)
           x = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

           ! lateral distance
           if (AMREX_SPACEDIM == 1) then
              r = 1.0_rt
              height = x
           else if (AMREX_SPACEDIM == 2) then
              r = x
              height = y
           else if (AMREX_SPACEDIM == 3) then
              r = sqrt(x**2 + y**2)
              height = z
           end if

           if (r < x_half_max) then
              f = 1.0_rt
           else if (r > x_half_max + x_half_width) then
              f = 0.0_rt
           else
              f = -(r - x_half_max)/x_half_width + ONE
           endif

           state(i,j,k,URHO)  = f * interpolate(height,gen_npts_model,gen_model_r(:,2), &
                                                gen_model_state(:,idens_model,2)) + &
                     (1.0_rt - f) * interpolate(height,gen_npts_model,gen_model_r(:,1), &
                                                gen_model_state(:,idens_model,1))

           state(i,j,k,UTEMP) = f * interpolate(height,gen_npts_model,gen_model_r(:,2), &
                                                gen_model_state(:,itemp_model,2)) + &
                     (1.0_rt - f) * interpolate(height,gen_npts_model,gen_model_r(:,1), &
                                                gen_model_state(:,itemp_model,1))

           temppres = f * interpolate(height,gen_npts_model,gen_model_r(:,2), &
                                             gen_model_state(:,ipres_model,2)) + &
                  (1.0_rt - f) * interpolate(height,gen_npts_model,gen_model_r(:,1), &
                                             gen_model_state(:,ipres_model,1))

           state(i,j,k,UFS:UFS-1+nspec) = ZERO

           do n = 1, nspec
              state(i,j,k,UFS-1+n) = f * interpolate(height,gen_npts_model,gen_model_r(:,2), &
                                                     gen_model_state(:,ispec_model-1+n,2)) + &
                          (1.0_rt - f) * interpolate(height,gen_npts_model,gen_model_r(:,1), &
                                                     gen_model_state(:,ispec_model-1+n,1))
           enddo

           eos_state%rho = state(i,j,k,URHO)
           eos_state%T = state(i,j,k,UTEMP)
           eos_state%p = temppres
           eos_state%xn(:) = state(i,j,k,UFS:UFS-1+nspec)

           call eos(eos_input_rp, eos_state)

           state(i,j,k,UTEMP) = eos_state % T
           state(i,j,k,UEINT) = eos_state % rho * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,UEINT)

           ! Initial velocities = 0
           state(i,j,k,UMX:UMZ) = 0.e0_rt

           ! convert to partial densities
           do n = 1, nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

        end do
     end do
  end do

end subroutine ca_initdata
