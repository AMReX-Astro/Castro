! ::: 
! ::: ------------------------------------------------------------------
! ::: 

     subroutine ca_estdt_burning(u,u_l1,u_l2,u_h1,u_h2, &
                                 enuc,e_l1,e_l2,e_h1,e_h2, &
                                 omegadot,o_l1,o_l2,o_h1,o_h2,lo,hi,&
                                 dt_hydro, dt_nuc_out, initial_step)

     use network, only : nspec, naux
     use eos_module
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                    allow_negative_energy
     use eos_module 
     use burner_dt_module	
     
     implicit none

     integer          :: u_l1,u_l2,u_h1,u_h2
     integer          :: o_l1,o_l2,o_h1,o_h2
     integer          :: e_l1,e_l2,e_h1,e_h2
     integer          :: lo(2), hi(2)
     double precision ::        u(u_l1:u_h1,u_l2:u_h2,NVAR)
     double precision ::     enuc(e_l1:e_h1,e_l2:e_h2)
     double precision :: omegadot(o_l1:o_h1,o_l2:o_h2,nspec)
     double precision :: dx,dt,nucdt
     integer          :: initial_step

     ! Note that we don't change dt_hydro
     ! In this routine we define dt_nuc
     double precision :: dt_hydro, dt_nuc, dt_nuc_out

     double precision :: rho, T, xin(nspec), xout(nspec)
     double precision :: ein, eout
     double precision :: rhoInv,ux

     integer          :: i,j

!    NOTE:  enuc contains change in enuc / dt
!    NOTE:  omegadot_i contains only change in X_i, not scaled by dt.


     dt_nuc     = 1.d200
     dt_nuc_out = 1.d200

     if (initial_step .eq. 1) then

        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           rho          = u(i,j,URHO)
           T            = u(i,j,UTEMP)
           xin(1:nspec) = u(i,j,UFS:UFS+nspec-1) / rho

           call ca_init_nuclear_esdt(rho,T,xin,dt_nuc)

           dt_nuc_out = min(dt_nuc,dt_nuc_out)

        enddo
        enddo

     else 

        do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           rho          = u(i,j,URHO)
           T            = u(i,j,UTEMP)
           xin(1:nspec) = u(i,j,UFS:UFS+nspec-1) / rho
           ein          = u(i,j,UEINT) / rho

           eout = ein + enuc(i,j) * dt_hydro

           xout(1:nspec) = xin(1:nspec) + omegadot(i,j,1:nspec)

           call ca_nuclear_esdt(xin,xout,ein,eout,dt_hydro,dt_nuc)

           dt_nuc_out = min(dt_nuc,dt_nuc_out)

        enddo
        enddo

     end if

     end subroutine ca_estdt_burning
