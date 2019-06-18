module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  real(rt), allocatable, managed :: T_l, T_r, dens, cfrac, ofrac, w_T, center_T, smallx, vel

  integer,  allocatable, managed :: idir

  integer,  allocatable, managed :: ihe4, ic12, io16
  real(rt), allocatable, managed :: xn(:)

  logical,  allocatable, managed :: fill_ambient_bc

  real(rt), allocatable, managed :: ambient_dens
  real(rt), allocatable, managed :: ambient_temp
  real(rt), allocatable, managed :: ambient_comp(:)
  real(rt), allocatable, managed :: ambient_e_l, ambient_e_r

contains

  ! Given a zone state, fill it with ambient material.

  subroutine fill_ambient(state, s_lo, s_hi, i, j, k, x, time)

    use amrex_constants_module, only: ZERO, HALF
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEINT, UEDEN, UFS
#ifdef GRAVITY
    use meth_params_module, only: gravity_type_int, const_grav
#endif
    use network, only: nspec
    use prob_params_module, only: problo, probhi

    implicit none

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer,  intent(in   ) :: i, j, k
    real(rt), intent(in   ) :: x, time

    real(rt) :: c_T

    real(rt) :: vel_g

    !$gpu

    c_T = problo(1) + center_T * (probhi(1) - problo(1))

    state(i,j,k,URHO) = ambient_dens
    state(i,j,k,UFS:UFS-1+nspec) = ambient_dens * ambient_comp

    if (x < c_T) then
       state(i,j,k,UTEMP) = T_l
       state(i,j,k,UEINT) = state(i,j,k,URHO) * ambient_e_l
       state(i,j,k,UMX) = state(i,j,k,URHO) * vel
#ifdef GRAVITY
       if (gravity_type_int == 0) then
          state(i,j,k,UMX) = state(i,j,k,UMX) + state(i,j,k,URHO) * const_grav * time
       end if
#endif
    else
       state(i,j,k,UTEMP) = T_r
       state(i,j,k,UEINT) = state(i,j,k,URHO) * ambient_e_r
       state(i,j,k,UMX) = -state(i,j,k,URHO) * vel
#ifdef GRAVITY
       if (gravity_type_int == 0) then
          state(i,j,k,UMX) = state(i,j,k,UMX) + state(i,j,k,URHO) * const_grav * time
       end if
#endif
    end if

    state(i,j,k,UMY:UMZ) = ZERO
    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * sum(state(i,j,k,UMX:UMZ)**2) / state(i,j,k,URHO)

  end subroutine fill_ambient

end module probdata_module
