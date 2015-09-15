     subroutine ca_estdt(lo,hi,u,u_lo,u_hi,dx,dt)

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

     type (eos_t_3D) :: eos_state

     grid_scl = (dx(1)*dx(2)*dx(3))**THIRD

     call eos_allocate(eos_state, lo, hi)
     
     ! Call EOS for the purpose of computing sound speed

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              rhoInv = ONE / u(i,j,k,URHO)
              
              eos_state % rho(i,j,k)  = u(i,j,k,URHO )
              eos_state % T(i,j,k)    = u(i,j,k,UTEMP)
              eos_state % e(i,j,k)    = u(i,j,k,UEINT) * rhoInv
              eos_state % xn(i,j,k,:) = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
              eos_state % aux(i,j,k,1:naux) = u(i,j,k,UFX:UFX+naux-1) * rhoInv
           enddo
        enddo
     enddo

     if (allow_negative_energy .eq. 0) eos_state % reset = .true.
     
     call eos(eos_input_re, eos_state)

     ! Compute velocity and then calculate CFL timestep.     
     
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
     
              ux = u(i,j,k,UMX) * rhoInv
              uy = u(i,j,k,UMY) * rhoInv
              uz = u(i,j,k,UMZ) * rhoInv

              if (UESGS .gt. -1) &
                 sqrtK = dsqrt( rhoInv*u(i,j,k,UESGS) )

              c = eos_state % cs(i,j,k)
              
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

     call eos_deallocate(eos_state)
     
     end subroutine ca_estdt
