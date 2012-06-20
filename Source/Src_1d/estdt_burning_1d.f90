! ::: 
! ::: ------------------------------------------------------------------
! ::: 

     subroutine ca_estdt_burning(u,u_l1,u_h1, &
                                 enuc,e_l1,e_h1, &
                                 omegadot,o_l1,o_h1,lo,hi, &
                                 dt_hydro, dt_nuc_out, initial_step)

     use network, only : nspec, naux
     use eos_module
     use meth_params_module, only : NVAR, URHO, UMX, UEINT, UTEMP, UFS, UFX
     use burner_dt_module	

     implicit none

     integer          :: u_l1,u_h1
     integer          :: o_l1,o_h1
     integer          :: e_l1,e_h1
     integer          :: lo(1),hi(1)
     double precision :: u(u_l1:u_h1,NVAR)
     double precision :: enuc(e_l1:e_h1)
     double precision :: omegadot(o_l1:o_h1,nspec)
     integer          :: initial_step

     ! Note that we don't change dt_hydro
     ! In this routine we define dt_nuc
     double precision :: dt_hydro, dt_nuc, dt_nuc_out

     double precision :: rho, T, xin(nspec), xout(nspec)
     double precision :: ein, eout
     double precision :: rhoInv,ux

     integer          :: i

!    NOTE:  enuc contains change in enuc / dt
!    NOTE:  omegadot_i contains only change in X_i, not scaled by dt.

     dt_nuc     = 1.d200
     dt_nuc_out = 1.d200

     if (initial_step .eq. 1) then

        do i = lo(1),hi(1)
           rho          = u(i,URHO)
           T            = u(i,UTEMP)
           xin(1:nspec) = u(i,UFS:UFS+nspec-1) / rho

           call ca_init_nuclear_esdt(rho,T,xin,dt_nuc)

           dt_nuc_out = min(dt_nuc,dt_nuc_out)

        enddo

     else 

        do i = lo(1),hi(1)

           rho          = u(i,URHO)
           T            = u(i,UTEMP)
           xin(1:nspec) = u(i,UFS:UFS+nspec-1) / rho
           ein          = u(i,UEINT) / rho

           eout = ein + enuc(i) * dt_hydro
		
           xout(1:nspec) = xin(1:nspec) + omegadot(i,1:nspec)

           call ca_nuclear_esdt(xin,xout,ein,eout,dt_hydro,dt_nuc)

           dt_nuc_out = min(dt_nuc,dt_nuc_out)

        enddo
	
     end if
 
     end subroutine ca_estdt_burning
