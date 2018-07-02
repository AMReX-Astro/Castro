module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  real(rt), save :: T_l, T_r, dens, cfrac, w_T, center_T, smallx, vel

  integer, save :: idir

  real(rt), save :: center(3)
      
  integer, save :: ihe4, ic12, io16
  real(rt), save, allocatable :: xn(:)

  logical, save :: fill_ambient_bc

  real(rt), save :: ambient_dens
  real(rt), save :: ambient_temp
  real(rt), save, allocatable :: ambient_comp(:)
  real(rt), save :: ambient_e_l, ambient_e_r

contains

  ! Given a zone state, fill it with ambient material.

  subroutine fill_ambient(state, x, time)

    use amrex_constants_module, only: ZERO, HALF
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEINT, UEDEN, UFS
#ifdef GRAVITY
    use meth_params_module, only: gravity_type, const_grav
#endif
    use network, only: nspec
    use prob_params_module, only: problo, probhi

    implicit none

    real(rt), intent(in   ) :: x, time
    real(rt), intent(inout) :: state(NVAR)

    real(rt) :: c_T

    real(rt) :: vel_g

    c_T = problo(1) + center_T * (probhi(1) - problo(1))

    state(URHO) = ambient_dens
    state(UFS:UFS-1+nspec) = ambient_dens * ambient_comp

    if (x < c_T) then
       state(UTEMP) = T_l
       state(UEINT) = state(URHO) * ambient_e_l
       state(UMX) = state(URHO) * vel
#ifdef GRAVITY
       if (gravity_type == "ConstantGrav") then
          state(UMX) = state(UMX) + state(URHO) * const_grav * time
       end if
#endif
    else
       state(UTEMP) = T_r
       state(UEINT) = state(URHO) * ambient_e_r
       state(UMX) = -state(URHO) * vel
#ifdef GRAVITY
       if (gravity_type == "ConstantGrav") then
          state(UMX) = state(UMX) + state(URHO) * const_grav * time
       end if
#endif
    end if

    state(UMY:UMZ) = ZERO
    state(UEDEN) = state(UEINT) + HALF * sum(state(UMX:UMZ)**2) / state(URHO)

  end subroutine fill_ambient

end module probdata_module
