module thermo_sources
    ! Functions for implementing the source terms to the internal energy
    ! equation, for those integration methods where we treat it as a
    ! source and do not explicitly discretize it in the conservative
    ! update.

  implicit none

contains

  subroutine ca_thermo_src(lo, hi, &
                           old_state, os_lo, os_hi, &
                           new_state, ns_lo, ns_hi, &
                           src, src_lo, src_hi, problo, dx, time, dt) bind(C, name="ca_thermo_src")
   !  Compute thermodynamic sources for the internal energy equation.
   ! At the moment, this is only the -p div{U} term in the internal energy
   ! equation, and only for method-of-lines integration, including the new
   ! SDC method (the `-` is because it is on the RHS of the equation)

    use amrex_constants_module, only: ZERO, HALF, FOURTH
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UFS, UFX, UEINT, NSRC
    use amrex_fort_module, only : rt => amrex_real
    use network, only : nspec, naux
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use prob_params_module, only : coord_type

    implicit none

    integer, intent(in) :: lo(3),hi(3)   ! bounds of the box to operate on
    integer, intent(in) :: os_lo(3),os_hi(3)   ! bounds of the old_state array
    integer, intent(in) :: ns_lo(3),ns_hi(3)   ! bounds of the new_state array
    integer, intent(in) :: src_lo(3),src_hi(3)   ! bounds of the src array
    real(rt), intent(in) :: old_state(os_lo(1):os_hi(1),os_lo(2):os_hi(2),os_lo(3):os_hi(3),NVAR)   ! the old-time hydrodynamic conserved state
    real(rt), intent(in) :: new_state(ns_lo(1):ns_hi(1),ns_lo(2):ns_hi(2),ns_lo(3):ns_hi(3),NVAR)   ! the new-time hydrodynamic conserved state
    real(rt), intent(inout) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NSRC)   ! the source terms for the conserved variables
    real(rt), intent(in) :: problo(3)   ! physical coordinates of the lower left corner of the domain
    real(rt), intent(in) :: dx(3)   ! grid spacing
    real(rt), intent(in), value :: time   ! current simulation time
    real(rt), intent(in), value :: dt   ! current timestep

    integer :: i, j, k

    type(eos_t) :: eos_state_old, eos_state_new
    real(rt) :: r, rp, rm

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! radius for non-Cartesian
             rp = problo(1) + (dble(i) + 1.5_rt)*dx(1)
             rm = problo(1) + (dble(i) - HALF)*dx(1)
             r = HALF*(rm + rp)

             ! compute -div{U}
             if (coord_type == 0) then
                src(i,j,k,UEINT) = -FOURTH*((old_state(i+1,j,k,UMX)/old_state(i+1,j,k,URHO) + &
                                             new_state(i+1,j,k,UMX)/new_state(i+1,j,k,URHO)) - &
                                            (old_state(i-1,j,k,UMX)/old_state(i-1,j,k,URHO) + &
                                             new_state(i-1,j,k,UMX)/new_state(i-1,j,k,URHO)))/dx(1)

             else if (coord_type == 1) then
                ! axisymmetric
                src(i,j,k,UEINT) = -FOURTH*(rp*(old_state(i+1,j,k,UMX)/old_state(i+1,j,k,URHO) + &
                                                new_state(i+1,j,k,UMX)/new_state(i+1,j,k,URHO)) - &
                                            rm*(old_state(i-1,j,k,UMX)/old_state(i-1,j,k,URHO) + &
                                                new_state(i-1,j,k,UMX)/new_state(i-1,j,k,URHO)))/(r*dx(1))

             else if (coord_type == 2) then
                ! spherical
                src(i,j,k,UEINT) = -FOURTH*(rp**2*(old_state(i+1,j,k,UMX)/old_state(i+1,j,k,URHO) + &
                                                   new_state(i+1,j,k,UMX)/new_state(i+1,j,k,URHO)) - &
                                            rm**2*(old_state(i-1,j,k,UMX)/old_state(i-1,j,k,URHO) + &
                                                   new_state(i-1,j,k,UMX)/new_state(i-1,j,k,URHO)))/(r**2*dx(1))
             endif

#if BL_SPACEDIM >= 2
             src(i,j,k,UEINT) = src(i,j,k,UEINT) - &
                                FOURTH*((old_state(i,j+1,k,UMY)/old_state(i,j+1,k,URHO) + &
                                         new_state(i,j+1,k,UMY)/new_state(i,j+1,k,URHO)) - &
                                        (old_state(i,j-1,k,UMY)/old_state(i,j-1,k,URHO) + &
                                         new_state(i,j-1,k,UMY)/new_state(i,j-1,k,URHO)))/dx(2)
#endif
#if BL_SPACEDIM == 3
             src(i,j,k,UEINT) = src(i,j,k,UEINT) - &
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
             eos_state_old % aux(:) = old_state(i,j,k,UFX:UFX+naux-1)/old_state(i,j,k,URHO)

             call eos(eos_input_rt, eos_state_old)

             eos_state_new % rho = new_state(i,j,k,URHO)
             eos_state_new % T = new_state(i,j,k,UTEMP)
             eos_state_new % xn(:) = new_state(i,j,k,UFS:UFS-1+nspec)/new_state(i,j,k,URHO)
             eos_state_new % aux(:) = new_state(i,j,k,UFX:UFX+naux-1)/new_state(i,j,k,URHO)

             call eos(eos_input_rt, eos_state_new)

             ! final source term, -p div{U}
             src(i,j,k,UEINT) = HALF*(eos_state_old % p + eos_state_new % p)*src(i,j,k,UEINT)

          enddo
       enddo
    enddo

  end subroutine ca_thermo_src

end module thermo_sources
