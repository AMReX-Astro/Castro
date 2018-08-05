! compute thermodynamic sources for the conserved hydro equations.
! At the moment, this is only the p div{U} term in the internal energy
! equation, and only for method-of-lines integration

module thermo_sources

  implicit none

contains

  subroutine ca_thermo_src(lo, hi, &
                           old_state, os_lo, os_hi, &
                           new_state, ns_lo, ns_hi, &
                           src, src_lo, src_hi, problo, dx, time, dt) bind(C, name="ca_thermo_src")

    use amrex_constants_module, only: ZERO, HALF, FOURTH
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UFS, UEINT
    use amrex_fort_module, only : rt => amrex_real
    use network, only : nspec
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos

    implicit none

    integer, intent(in) :: lo(3),hi(3)
    integer, intent(in) :: os_lo(3),os_hi(3)
    integer, intent(in) :: ns_lo(3),ns_hi(3)
    integer, intent(in) :: src_lo(3),src_hi(3)
    real(rt), intent(in) :: old_state(os_lo(1):os_hi(1),os_lo(2):os_hi(2),os_lo(3):os_hi(3),NVAR)
    real(rt), intent(in) :: new_state(ns_lo(1):ns_hi(1),ns_lo(2):ns_hi(2),ns_lo(3):ns_hi(3),NVAR)
    real(rt), intent(inout) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
    real(rt), intent(in) :: problo(3),dx(3)
    real(rt), intent(in), value :: time,dt

    integer :: i, j, k

    type(eos_t) :: eos_state_old, eos_state_new

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             src(i,j,k,UEINT) = FOURTH*((old_state(i+1,j,k,UMX)/old_state(i+1,j,k,URHO) + &
                                         new_state(i+1,j,k,UMX)/new_state(i+1,j,k,URHO)) - &
                                        (old_state(i-1,j,k,UMX)/old_state(i-1,j,k,URHO) + &
                                         new_state(i-1,j,k,UMX)/new_state(i-1,j,k,URHO)))/dx(1)
#if BL_SPACEDIM >= 2
             src(i,j,k,UEINT) = src(i,j,k,UEINT) + &
                                FOURTH*((old_state(i,j+1,k,UMY)/old_state(i,j+1,k,URHO) + &
                                         new_state(i,j+1,k,UMY)/new_state(i,j+1,k,URHO)) - &
                                        (old_state(i,j-1,k,UMY)/old_state(i,j-1,k,URHO) + &
                                         new_state(i,j-1,k,UMY)/new_state(i,j-1,k,URHO)))/dx(2)
#endif
#if BL_SPACEDIM == 3
             src(i,j,k,UEINT) = src(i,j,k,UEINT) + &
                                FOURTH*((old_state(i,j,k+1,UMZ)/old_state(i,j,k+1,URHO) + &
                                         new_state(i,j,k+1,UMZ)/new_state(i,j,k+1,URHO)) - &
                                        (old_state(i,j,k-1,UMZ)/old_state(i,j,k-1,URHO) + &
                                         new_state(i,j,k-1,UMZ)/new_state(i,j,k-1,URHO)))/dx(3)
#endif

             ! we now need the pressure -- we will assume that the
             ! temperature is consistent with the input state
             eos_state_old % rho = old_state(i,j,k,URHO)
             eos_state_old % T = old_state(i,j,k,UTEMP)
             eos_state_old % xn(:) = old_state(i,j,k,UFS:UFS-1+nspec)/old_state(i,j,k,URHO)

             call eos(eos_input_rt, eos_state_old)

             eos_state_new % rho = new_state(i,j,k,URHO)
             eos_state_new % T = new_state(i,j,k,UTEMP)
             eos_state_new % xn(:) = new_state(i,j,k,UFS:UFS-1+nspec)/new_state(i,j,k,URHO)

             call eos(eos_input_rt, eos_state_new)

             src(i,j,k,UEINT) = HALF*(eos_state_old % p + eos_state_new % p)*src(i,j,k,UEINT)

          enddo
       enddo
    enddo

  end subroutine ca_thermo_src

end module thermo_sources
