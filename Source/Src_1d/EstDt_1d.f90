      subroutine ca_estdt(u,u_l1,u_h1,lo,hi,dx,dt)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      implicit none

      integer u_l1,u_h1
      integer lo(1), hi(1)
      double precision u(u_l1:u_h1,NVAR)
      double precision dx(1), dt

      double precision :: p, e, gamc, c, T, dpdr, dpde, xn(nspec+naux)
      double precision :: rhoInv,ux,dt1
      integer          :: i,n
      integer          :: pt_index(1)

!     Translate to primitive variables, compute sound speed (call eos), get dtmax
      do i = lo(1),hi(1)

         rhoInv = 1.d0/u(i,URHO)

         ux = u(i,UMX)*rhoInv
         T  = u(i,UTEMP)

         ! Use internal energy for calculating dt 
         e  = u(i,UEINT)*rhoInv

         xn(1:nspec)=u(i,UFS:UFS+nspec-1)*rhoInv

         if (naux > 0) &
            xn(nspec+1:nspec+naux)=u(i,UFX:UFX+naux-1)*rhoInv

         ! Protect against negative eint
         if (e .gt. 0.d0 .or. allow_negative_energy .eq. 1) then
            pt_index(1) = i
            call eos_given_ReX(gamc,p,c,T,dpdr,dpde,u(i,URHO),e,xn,pt_index=pt_index)
         else
            c = 0.d0
         end if

         dt1 = dx(1) /( c + abs(ux) )
         dt  = min(dt,dt1)

      enddo

      end subroutine ca_estdt
