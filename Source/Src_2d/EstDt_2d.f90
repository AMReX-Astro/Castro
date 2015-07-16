     subroutine ca_estdt(u,u_l1,u_l2,u_h1,u_h2,lo,hi,dx,dt)

     use network, only : nspec, naux
     use eos_module
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                    allow_negative_energy
     use bl_constants_module

     implicit none

     integer          :: u_l1,u_l2,u_h1,u_h2
     integer          :: lo(2), hi(2)
     double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)
     double precision :: dx(2),dt

     double precision :: rhoInv,c,ux,uy,dt1,dt2
     integer          :: i,j

     type (eos_t) :: eos_state

!    Translate to primitive variables, compute sound speed (call eos), get dtmax
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)

            rhoInv = ONE / u(i,j,URHO)

            ux = u(i,j,UMX) * rhoInv
            uy = u(i,j,UMY) * rhoInv

            ! Use internal energy for calculating dt 
            eos_state % e = u(i,j,UEINT) * rhoInv

            ! Protect against negative e
            if (eos_state % e .gt. ZERO .or. allow_negative_energy .eq. 1) then

               eos_state % rho = u(i,j,URHO)
               eos_state % T   = u(i,j,UTEMP)
               eos_state % xn  = u(i,j,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux = u(i,j,UFX:UFX+naux-1) * rhoInv

               call eos(eos_input_re, eos_state)
              
               c = eos_state % cs
            else
               c = ZERO
            end if

            dt1 = dx(1)/(c + abs(ux))
            dt2 = dx(2)/(c + abs(uy))
            dt = min(dt,dt1,dt2)
         enddo
      enddo

      end subroutine ca_estdt
