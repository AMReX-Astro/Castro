  subroutine call_to_thornado(lo, hi, dt, &
                              S , s_lo, s_hi, ns , &
                              dS, d_lo, d_hi, nds, &
                              U_R_o, U_R_o_lo, U_R_o_hi, n_uro, &
                              U_R_n, U_R_n_lo, U_R_n_hi, n_urn, &
                              n_fluid_dof, n_moments, ng) &
                              bind(C, name="call_to_thornado")

    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEINT,UEDEN,UFX
    use ProgramHeaderModule, only : nE, nDOF, nNodesX, nNodesE
    use FluidFieldsModule, only : uCF, iCF_D, iCF_S1, iCF_S2, iCF_E, iCF_Ne
    use RadiationFieldsModule, only : nSpecies, uCR
    use TimeSteppingModule_Castro, only : Update_IMEX_PDARS
    use UnitsModule, only : Gram, Centimeter, Second

    implicit none
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) ::  s_lo(2),  s_hi(2)
    integer, intent(in) ::  d_lo(2),  d_hi(2)
    integer, intent(in) ::  U_R_o_lo(2),  U_R_o_hi(2)
    integer, intent(in) ::  U_R_n_lo(2),  U_R_n_hi(3)
    integer, intent(in) ::  ns, nds, n_uro, n_urn
    integer, intent(in) ::  n_fluid_dof, n_moments
    integer, intent(in) :: ng
    real(rt), intent(in) :: dt

    ! Here we expect  n_rad_comp = 20 x 16 x 6 x 4 (energy x dof x species x moments)

    ! Conserved fluid state (rho, rho u, rho v, rho w, rho E, rho e...)
    real(rt), intent(inout) ::  S(s_lo(1):s_hi(1),s_lo(2):s_hi(2),ns) 
    real(rt), intent(inout) :: dS(d_lo(1):d_hi(1),d_lo(2):d_hi(2),nds)

    ! Old and new radiation state
    real(rt), intent(inout) ::  U_R_o(U_R_o_lo(1): U_R_o_hi(1),  U_R_o_lo(2): U_R_o_hi(2), 0:n_uro-1) 
    real(rt), intent(inout) ::  U_R_n(U_R_n_lo(1): U_R_n_hi(1),  U_R_n_lo(2): U_R_n_hi(2), 0:n_urn-1) 

    ! Temporary variables
    integer  :: i,j,k
    integer  :: ic,jc,kc
    integer  :: ii,id,ie,im,is
    real(rt) :: conv_dens, conv_mom, conv_enr, conv_ne, conv_J, conv_H, testdt

    ! Sanity check on size of arrays
    ! Note that we have set ngrow_thornado = ngrow_state in Castro_setup.cpp
    if (U_R_o_lo(1) .ne. S_lo(1) .or.  U_R_o_hi(1) .ne. S_hi(1) .or. &
        U_R_o_lo(2) .ne. S_lo(2) .or.  U_R_o_hi(2) .ne. S_hi(2)) then
        print *,'INCONSISTENT ARRAY BOUNDS ON FLUID AND RADIATION VARS'
        print *,'U_R_lo: ', U_R_o_lo(:)
        print *,'U_R_hi: ', U_R_o_hi(:)
        print *,'  S_lo: ',   S_lo(:)
        print *,'  S_hi: ',   S_hi(:)
        stop
    endif
    ! End Sanity check 

    if (ng.ne. 2) &
      call amrex_abort("Need two ghost cells in call_to_thornado!")

    ! Zero out dS
    dS = 0.0e0

    conv_dens = Gram / Centimeter**3
    conv_mom  = Gram / Centimeter**2 / Second
    conv_enr  = Gram / Centimeter / Second**2
    conv_ne   = 1.d0 / Centimeter**3
    conv_J    = Gram/Second**2/Centimeter ! check that this is correct
    conv_H    = Gram/Second**3

    ! ************************************************************************************
    ! Copy from the Castro arrays into Thornado arrays from InitThornado_Patch
    ! ************************************************************************************
    do jc = lo(2)-ng,hi(2)+ng
    do ic = lo(1)-ng,hi(1)+ng

         ! The uCF array was allocated in CreateFluidFieldsConserved with 
         !     ALLOCATE( uCF &
         !      (1:nDOFX, &
         !       1-swX(1):nX(1)+swX(1), &
         !       1-swX(2):nX(2)+swX(2), &
         !       1-swX(3):nX(3)+swX(3), &
         !       1:nCF) )

         ! U_R_o spatial indices start at lo - (number of ghost zones)
         ! uCR spatial indices start at 1 - (number of ghost zones)
         i = ic - lo(1) + 1
         j = jc - lo(2) + 1
         k = 1

         ! Thornado uses units where c = G = k = 1, Meter = 1
         uCF(1:n_fluid_dof,i,j,k,iCF_D)  = S(ic,jc,URHO)  * conv_dens
         uCF(1:n_fluid_dof,i,j,k,iCF_S1) = S(ic,jc,UMX)   * conv_mom
         uCF(1:n_fluid_dof,i,j,k,iCF_S2) = S(ic,jc,UMY)   * conv_mom
         uCF(1:n_fluid_dof,i,j,k,iCF_E)  = S(ic,jc,UEDEN) * conv_enr
         uCF(1:n_fluid_dof,i,j,k,iCF_Ne) = S(ic,jc,UFX)   * conv_ne

         ! The uCF array was allocated in CreateRadiationdFields_Conserved with 
         ! ALLOCATE &
         !   ( uCR(1:nDOF, &
         !         1-swE:nE+swE, &
         !         1-swX(1):nX(1)+swX(1), &
         !         1-swX(2):nX(2)+swX(2), &
         !         1:nCR, 1:nSpecies) )

         do is = 1, nSpecies
         do im = 1, n_moments
         do ie = 1, nE
         do id = 1, nDOF
            ii   = (is-1)*(n_moments*nE*nDOF) + (im-1)*(nE*nDOF) + (ie-1)*nDOF + (id-1)
!            if (im .eq. 1) uCR(id,ie,i,j,k,im,is) = U_R_o(ic,jc,ii)*conv_J
            if (im .eq. 1) uCR(id,ie,i,j,k,im,is) = U_R_o(ic,jc,ii)
!            if (im   >  1) uCR(id,ie,i,j,k,im,is) = U_R_o(ic,jc,ii)*conv_H
            if (im   >  1) uCR(id,ie,i,j,k,im,is) = U_R_o(ic,jc,ii)
            write(*,*) 'jc,ic,is,im,ie,id,uCR=',jc,ic,is,im,ie,id,uCR(id,ie,i,j,k,im,is)
         end do
         end do
         end do
         end do

    end do
    end do

    ! ************************************************************************************
    ! Call the Fortran interface that lives in the thornado repo
    ! ************************************************************************************
    call Update_IMEX_PDARS(dt*Second, uCF, uCR)

    ! ************************************************************************************
    ! Copy back from the thornado arrays into Castro arrays
    ! ************************************************************************************
    do jc = lo(2),hi(2)
    do ic = lo(1),hi(1)

         ! uCR spatial indices start at 1 - ng
         ! U_R_n spatial indices start at lo
         i = ic - lo(1) + 1
         j = jc - lo(2) + 1
         k = 1

         ! We store dS as a source term which we can add to S outside of this routine
         ! uCF returned as a cell-averaged quantity so all components are the same,
         !  can just use the first component
         ! Update_IMEX_PC2 doesn't currently change the fluid density or momentum

!         dS(ic,jc,URHO ) = uCF(1,i,j,k,iCF_D)  / conv_dens - S(i,j,URHO)
!         dS(ic,jc,UMX  ) = uCF(1,i,j,k,iCF_S1) / conv_mom  - S(i,j,UMX)
!         dS(ic,jc,UMY  ) = uCF(1,i,j,k,iCF_S2) / conv_mom  - S(i,j,UMY)
          dS(ic,jc,UEDEN) = uCF(1,i,j,k,iCF_E)  / conv_enr  - S(i,j,UEDEN)
          dS(ic,jc,UEINT) = dS(ic,jc,UEDEN)     ! TRUE IFF NO MOMENTUM SOURCE TERMS
          dS(ic,jc,UFX  ) = uCF(1,i,j,k,iCF_Ne) / conv_ne   - S(i,j,UFX)

         do is = 1, nSpecies
         do im = 1, n_moments
         do ie = 1, nE
         do id = 1, nDOF
            ii   = (is-1)*(n_moments*nE*nDOF) + (im-1)*(nE*nDOF) + (ie-1)*nDOF + (id-1)

!            if (im .eq. 1) U_R_n(ic,jc,ii) = uCR(id,ie,i,j,k,im,is)/conv_J
!            if (im   >  1) U_R_n(ic,jc,ii) = uCR(id,ie,i,j,k,im,is)/conv_H
            if (im .eq. 1) U_R_n(ic,jc,ii) = uCR(id,ie,i,j,k,im,is)
            if (im   >  1) U_R_n(ic,jc,ii) = uCR(id,ie,i,j,k,im,is)

         end do
         end do
         end do
         end do

    end do
    end do

  end subroutine call_to_thornado
