     subroutine ca_estdt(u,u_l1,u_l2,u_h1,u_h2,lo,hi,dx,dt)

     use network, only : nspec, naux
     use eos_module
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                    allow_negative_energy

     implicit none

     integer          :: u_l1,u_l2,u_h1,u_h2
     integer          :: lo(2), hi(2)
     double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)
     double precision :: dx(2),dt

     double precision :: p, e, gamc, c, T, dpdr, dpde, xn(nspec+naux)
     double precision :: rhoInv,ux,uy,dt1,dt2
     integer          :: i,j
     integer          :: pt_index(2)

!    Translate to primitive variables, compute sound speed (call eos), get dtmax
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)

            rhoInv = 1.d0/u(i,j,URHO)

            ux = u(i,j,UMX)*rhoInv
            uy = u(i,j,UMY)*rhoInv
            T  = u(i,j,UTEMP)

            ! Use internal energy for calculating dt 
            e  = u(i,j,UEINT)*rhoInv

            xn(1:nspec)=u(i,j,UFS:UFS+nspec-1)*rhoInv

            if (naux > 0) &
               xn(nspec+1:nspec+naux)=u(i,j,UFX:UFX+naux-1)*rhoInv

            ! Protect against negative e
            if (e .gt. 0.d0 .or. allow_negative_energy.eq.1) then
               pt_index(1) = i
               pt_index(2) = j
               call eos_given_ReX(gamc, p, c, T, dpdr, dpde, u(i,j,URHO), e, xn, pt_index=pt_index)
            else
               c = 0.d0
            end if

            dt1 = dx(1)/(c + abs(ux))
            dt2 = dx(2)/(c + abs(uy))
            dt = min(dt,dt1,dt2)
         enddo
      enddo

      end subroutine ca_estdt
