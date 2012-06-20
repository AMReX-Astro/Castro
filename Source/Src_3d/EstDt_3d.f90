     subroutine ca_estdt(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi,dx,dt)

     use network, only : nspec, naux
     use eos_module
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UESGS, UTEMP, UFS, &
                                    UFX, allow_negative_energy

     implicit none

     integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
     integer          :: lo(3), hi(3)
     double precision :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
     double precision :: dx(3), dt

     double precision :: p, e, gamc, c, T, dpdr, dpde, xn(nspec+naux)
     double precision :: rhoInv,ux,uy,uz,dt1,dt2,dt3
     double precision :: sqrtK,grid_scl,dt4
     integer          :: i,j,k
     integer          :: pt_index(3)

     double precision, parameter :: onethird = 1.d0/3.d0

     grid_scl = (dx(1)*dx(2)*dx(3))**onethird

     ! Translate to primitive variables, compute sound speed (call eos)
     !$OMP PARALLEL DO PRIVATE(i,j,k,rhoInv,ux,uy,uz,e,T,sqrtK,xn,pt_index,gamc,p,c,dpdr,dpde,dt1,dt2,dt3) REDUCTION(min:dt)
     do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               rhoInv = 1.d0/u(i,j,k,URHO)

               ux = u(i,j,k,UMX)*rhoInv
               uy = u(i,j,k,UMY)*rhoInv
               uz = u(i,j,k,UMZ)*rhoInv
               T  = u(i,j,k,UTEMP)

               ! Use internal energy for calculating dt 
               e  = u(i,j,k,UEINT)*rhoInv

               if (UESGS .gt. -1) &
                  sqrtK = dsqrt( rhoInv*u(i,j,k,UESGS) )

               xn(1:nspec)=u(i,j,k,UFS:UFS+nspec-1)*rhoInv

               if (naux > 0) &
                  xn(nspec+1:nspec+naux)=u(i,j,k,UFX:UFX+naux-1)*rhoInv

               ! Protect against negative e
               if (e .gt. 0.d0 .or. allow_negative_energy.eq.1) then
                  pt_index(1) = i
                  pt_index(2) = j
                  pt_index(3) = k
                  call eos_given_ReX(gamc,p,c,T,dpdr,dpde,u(i,j,k,URHO),e,xn,pt_index=pt_index)
               else
                  c = 0.d0
               end if

               dt1 = dx(1)/(c + abs(ux))
               dt2 = dx(2)/(c + abs(uy))
               dt3 = dx(3)/(c + abs(uz))
               dt = min(dt,dt1,dt2,dt3)

               ! Now let's check the diffusion terms for the SGS equations
               if (UESGS .gt. -1) then

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
     !$OMP END PARALLEL DO

     end subroutine ca_estdt
