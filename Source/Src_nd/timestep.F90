module timestep_module

  implicit none

  public
  
contains

  ! Courant-condition limited timestep
  
  subroutine ca_estdt(lo,hi,u,u_lo,u_hi,dx,dt) bind(C)

    use network, only: nspec, naux
    use eos_module
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEINT, UESGS, UTEMP, UFS, UFX, &
         allow_negative_energy
    use prob_params_module, only: dim
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_lo(3), u_hi(3)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    double precision :: dx(3), dt

    double precision :: rhoInv, ux, uy, uz, c, dt1, dt2, dt3
    double precision :: sqrtK, grid_scl, dt4
    integer          :: i, j, k

    type (eos_t) :: eos_state

    grid_scl = (dx(1)*dx(2)*dx(3))**THIRD

    if (allow_negative_energy .eq. 0) eos_state % reset = .true.

    ! Call EOS for the purpose of computing sound speed

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / u(i,j,k,URHO)

             eos_state % rho = u(i,j,k,URHO )
             eos_state % T   = u(i,j,k,UTEMP)
             eos_state % e   = u(i,j,k,UEINT) * rhoInv
             eos_state % xn  = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux = u(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos(eos_input_re, eos_state)

             ! Compute velocity and then calculate CFL timestep.     

             ux = u(i,j,k,UMX) * rhoInv
             uy = u(i,j,k,UMY) * rhoInv
             uz = u(i,j,k,UMZ) * rhoInv

             if (UESGS .gt. -1) &
                  sqrtK = dsqrt( rhoInv*u(i,j,k,UESGS) )

             c = eos_state % cs

             dt1 = dx(1)/(c + abs(ux))
             if (dim .ge. 2) then
                dt2 = dx(2)/(c + abs(uy))
             else
                dt2 = dt1
             endif
             if (dim .eq. 3) then
                dt3 = dx(3)/(c + abs(uz))
             else
                dt3 = dt1
             endif

             dt  = min(dt,dt1,dt2,dt3)

             ! Now let's check the diffusion terms for the SGS equations
             if (UESGS .gt. -1 .and. dim .eq. 3) then

                ! First for the term in the momentum equation
                ! This is actually dx^2 / ( 6 nu_sgs )
                ! Actually redundant as it takes the same form as below with different coeff
                ! dt4 = grid_scl / ( 0.42d0 * sqrtK )

                ! Now for the term in the K equation itself
                ! nu_sgs is 0.65
                ! That gives us 0.65*6 = 3.9
                ! Using 4.2 to be conservative (Mach1-256 broke during testing with 3.9)
                !               dt4 = grid_scl / ( 3.9d0 * sqrtK )
                dt4 = grid_scl / ( 4.2d0 * sqrtK )
                dt = min(dt,dt4)

             end if

          enddo
       enddo
    enddo

  end subroutine ca_estdt



  ! Reactions-limited timestep

#ifdef REACTIONS
  subroutine ca_estdt_burning(u, u_lo, u_hi, &
                              reactions, r_lo, r_hi, &
                              lo, hi, dx, dt_old, dt) bind(C)

    use bl_constants_module, only: ONE
    use network, only: nspec, naux
    use meth_params_module, only : NVAR, URHO, UEINT, UTEMP, UFS, UFX, dtnuc, dsnuc
    use prob_params_module, only : dim
    use actual_rhs_module, only: actual_rhs
    use eos_module
    use burner_module, only: ok_to_burn
    use burn_type_module
    use eos_type_module

    implicit none

    integer          :: u_lo(3), u_hi(3)
    integer          :: r_lo(3), r_hi(3)
    integer          :: lo(3), hi(3)
    double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    double precision :: reactions(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspec+2)
    double precision :: dx(3), dt, dt_old

    double precision :: dedt, dXdt(nspec)
    integer          :: i, j, k

    type (burn_t)    :: burn_state
    type (eos_t)     :: eos_state
    double precision :: rhoInv

    ! We want to limit the timestep so that it is not larger than
    ! dtnuc * (e / (de/dt)).  If the timestep factor dtnuc is
    ! equal to 1, this says that we don't want the
    ! internal energy to change by any more than its current
    ! magnitude in the next timestep. 
    !
    ! If dtnuc is less than one, it controls the fraction we will
    ! allow the internal energy to change in this timestep due to
    ! nuclear burning, provided that the last timestep's burning is a
    ! good estimate for the current timestep's burning.
    !
    ! We do the same thing for the species, using a timestep
    ! limiter dsnuc * (X_k / (dX_k/dt)).
    !
    ! To estimate de/dt and dX/dt, we are going to call the RHS of the
    ! burner given the current state data. We need to do an EOS
    ! call before we do the RHS call so that we have accurate
    ! values for the thermodynamic data like abar, zbar, etc.
    ! But we will call in (rho, T) mode, which is inexpensive.

    if (dtnuc > 1.d199 .and. dsnuc > 1.d199) return

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = ONE / u(i,j,k,URHO)

             burn_state % rho = u(i,j,k,URHO)
             burn_state % T   = u(i,j,k,UTEMP)
             burn_state % e   = u(i,j,k,UEINT) * rhoInv
             burn_state % xn  = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
             burn_state % aux = u(i,j,k,UFX:UFX+naux-1) * rhoInv

             if (.not. ok_to_burn(burn_state)) cycle

             call burn_to_eos(burn_state, eos_state)
             call eos(eos_input_rt, eos_state)
             call eos_to_burn(eos_state, burn_state)

             call actual_rhs(burn_state)

             dedt = max(abs(burn_state % ydot(net_ienuc)), 1.d-200)

             dt = min(dt, dtnuc * burn_state % e / dedt)

             dXdt = max(abs(burn_state % ydot(1:nspec) * aion), 1.d-200)

             dt = min(dt, dsnuc * minval(burn_state % xn / dXdt))

          enddo
       enddo
    enddo

  end subroutine ca_estdt_burning
#endif



  ! Diffusion-limited timestep

#ifdef DIFFUSION
  subroutine ca_estdt_diffusion(lo,hi,state,s_lo,s_hi,dx,dt) bind(C)

    use network, only: nspec, naux
    use eos_module
    use eos_type_module
    use meth_params_module, only: NVAR, URHO, UEINT, UTEMP, UFS, UFX, &
         diffuse_cutoff_density
    use prob_params_module, only: dim
    use bl_constants_module
    use conductivity_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: dx(3), dt

    double precision :: dt1, dt2, dt3, rho_inv
    integer          :: i, j, k, n
    double precision :: cond, D

    type (eos_t) :: eos_state

    ! dt < 0.5 dx**2 / D
    ! where D = k/(rho c_v), and k is the conductivity      
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,URHO) > diffuse_cutoff_density) then

                rho_inv = ONE/state(i,j,k,URHO)

                ! we need cv
                eos_state % rho = state(i,j,k,URHO )
                eos_state % T   = state(i,j,k,UTEMP)
                eos_state % e   = state(i,j,k,UEINT) * rho_inv

                eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rho_inv
                eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rho_inv

                call eos(eos_input_re, eos_state)

                ! we also need the conductivity
                call thermal_conductivity(eos_state, cond)

                ! maybe we should check (and take action) on negative cv here?
                D = cond*rho_inv/eos_state%cv

                dt1 = HALF*dx(1)**2/D

                if (dim >= 2) then
                   dt2 = HALF*dx(2)**2/D
                else
                   dt2 = dt1
                endif

                if (dim == 3) then
                   dt3 = HALF*dx(3)**2/D
                else
                   dt3 = dt1
                endif

                dt  = min(dt,dt1,dt2,dt3)

             endif

          enddo
       enddo
    enddo

  end subroutine ca_estdt_diffusion
#endif
  
end module timestep_module
