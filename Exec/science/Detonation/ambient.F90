module ambient_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  ! This is a state vector that contains "ambient" material
  ! that will be used for filling material in regions that
  ! are not of interest, like at the edges of the domain.

  real(rt), allocatable :: ambient_state(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: ambient_state
#endif

contains

  subroutine get_ambient_state(ambient_state_out, loc, time)

    use amrex_constants_module, only: ZERO, HALF
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEINT, UEDEN, UFS
#ifdef GRAVITY
    use meth_params_module, only: gravity_type_int, const_grav
#endif
    use network, only: nspec
    use prob_params_module, only: problo, probhi
    use probdata_module

    implicit none

    real(rt), intent(inout) :: ambient_state_out(NVAR)
    real(rt), intent(in   ) :: loc(3), time

    real(rt) :: c_T, vel_g

    !$gpu

    c_T = problo(1) + center_T * (probhi(1) - problo(1))

    ambient_state_out(URHO) = ambient_dens
    ambient_state_out(UFS:UFS-1+nspec) = ambient_dens * ambient_comp

    if (loc(1) < c_T) then
       ambient_state_out(UTEMP) = T_l
       ambient_state_out(UEINT) = ambient_state_out(URHO) * ambient_e_l
       ambient_state_out(UMX) = ambient_state_out(URHO) * vel
#ifdef GRAVITY
       if (gravity_type_int == 0) then
          ambient_state_out(UMX) = ambient_state_out(UMX) + ambient_state_out(URHO) * const_grav * time
       end if
#endif
    else
       ambient_state_out(UTEMP) = T_r
       ambient_state_out(UEINT) = ambient_state_out(URHO) * ambient_e_r
       ambient_state_out(UMX) = -ambient_state_out(URHO) * vel
#ifdef GRAVITY
       if (gravity_type_int == 0) then
          ambient_state_out(UMX) = ambient_state_out(UMX) + ambient_state_out(URHO) * const_grav * time
       end if
#endif
    end if

    ambient_state_out(UMY:UMZ) = ZERO
    ambient_state_out(UEDEN) = ambient_state_out(UEINT) + &
                               HALF * sum(ambient_state_out(UMX:UMZ)**2) / ambient_state_out(URHO)
    
  end subroutine get_ambient_state

end module ambient_module
