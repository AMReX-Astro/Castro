! ::: 
! ::: ------------------------------------------------------------------
! ::: 

     subroutine ca_estdt_burning(u,u_l1,u_l2,u_h1,u_h2, u_h3,u_h3, &
                                 enuc,e_l1,e_l2,e_h1,e_h2,e_h3,e_h3,&
                                 omegadot,o_l1,o_l2,o_h1,o_h2,o_h3,o_h3,lo,hi, &
                                 dt_hydro, dt_nuc_out, initial_step)

     use network, only : nspec, naux
     use eos_module
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                    allow_negative_energy
     use eos_module 
     use burner_dt_module	

     implicit none

     integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
     integer          :: o_l1,o_l2,o_l3,o_h1,o_h2,o_h3
     integer          :: e_l1,e_l2,e_l3,e_h1,e_h2,e_h3
     integer          :: lo(3), hi(3)
     double precision ::        u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
     double precision ::     enuc(e_l1:e_h1,e_l2:e_h2,e_l3:e_h3)
     double precision :: omegadot(o_l1:o_h1,o_l2:o_h2,o_l3:o_h3,nspec)
     double precision :: dx,dt,nucdt
     integer          :: initial_step

     ! Note that we don't change dt_hydro
     ! In this routine we define dt_nuc
     double precision :: dt_hydro, dt_nuc, dt_nuc_out

     double precision :: rho, T, xin(nspec), xout(nspec)
     double precision :: ein, eout
     double precision :: rhoInv,ux

     integer          :: i,j,k

!    NOTE:  enuc contains change in enuc / dt
!    NOTE:  omegadot_i contains only change in X_i, not scaled by dt.

     dt_nuc     = 1.d200
     dt_nuc_out = 1.d200

     if (initial_step .eq. 1) then

        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           rho          = u(i,j,k,URHO)
           T            = u(i,j,k,UTEMP)
           xin(1:nspec) = u(i,j,k,UFS:UFS+nspec-1) / rho

           call ca_init_nuclear_esdt(rho,T,xin,dt_nuc)

           dt_nuc_out = min(dt_nuc,dt_nuc_out)

        enddo
        enddo
        enddo

     else

        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           rho          = u(i,j,k,URHO)
           T            = u(i,j,k,UTEMP)
           xin(1:nspec) = u(i,j,k,UFS:UFS+nspec-1) / rho
           ein          = u(i,j,k,UEINT) / rho

           eout = ein + enuc(i,j,k) * dt_hydro

           xout(1:nspec) = xin(1:nspec) + omegadot(i,j,k,1:nspec)

           call ca_nuclear_esdt(xin,xout,ein,eout,dt_hydro,dt_nuc)

           dt_nuc_out = min(dt_nuc,dt_nuc_out)

        enddo
        enddo
        enddo

     end if

     end subroutine ca_estdt_burning
