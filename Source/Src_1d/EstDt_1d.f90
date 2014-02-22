      subroutine ca_estdt(u,u_l1,u_h1,lo,hi,dx,dt)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer u_l1,u_h1
      integer lo(1), hi(1)
      double precision u(u_l1:u_h1,NVAR)
      double precision dx(1), dt

      double precision :: rhoInv,c,ux,dt1
      integer          :: i
      integer          :: pt_index(1)

      type (eos_t) :: eos_state

!     Translate to primitive variables, compute sound speed (call eos), get dtmax
      do i = lo(1),hi(1)

         rhoInv = ONE / u(i,URHO)

         ux = u(i,UMX) * rhoInv
         eos_state % T  = u(i,UTEMP)

         ! Use internal energy for calculating dt 
         eos_state % e  = u(i,UEINT) * rhoInv

         ! Protect against negative e
         if (eos_state % e .gt. ZERO .or. allow_negative_energy .eq. 1) then

            eos_state % rho = u(i,URHO)
            eos_state % T   = u(i,UTEMP)
            eos_state % xn  = u(i,UFS:UFS+nspec-1) * rhoInv
            eos_state % aux = u(i,UFX:UFX+naux-1) * rhoInv

            pt_index(1) = i

            call eos(eos_input_re, eos_state, pt_index = pt_index)

            c = eos_state % cs
         else
            c = ZERO
         end if

         dt1 = dx(1) /( c + abs(ux) )
         dt  = min(dt,dt1)

      enddo

      end subroutine ca_estdt
