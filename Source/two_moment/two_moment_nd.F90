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
    use FluidFieldsModule, only : uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
    use RadiationFieldsModule, only : nSpecies, uCR
    use TimeSteppingModule_Castro, only : Update_IMEX_PC2
    use UnitsModule, only : Gram, Centimeter, Second

    implicit none
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) ::  d_lo(3),  d_hi(3)
    integer, intent(in) ::  U_R_o_lo(3),  U_R_o_hi(3)
    integer, intent(in) ::  U_R_n_lo(3),  U_R_n_hi(3)
    integer, intent(in) ::  ns, nds, n_uro, n_urn
    integer, intent(in) ::  n_fluid_dof, n_moments
    integer, intent(in) :: ng
    real(rt), intent(in) :: dt

    ! Here we expect  n_rad_comp = 20 x 16 x 6 x 4 (energy x dof x species x moments)

    ! Conserved fluid state (rho, rho u, rho v, rho w, rho E, rho e...)
    real(rt), intent(inout) ::  S(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ns) 
    real(rt), intent(inout) :: dS(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nds)

    ! Old and new radiation state
    real(rt), intent(inout) ::  U_R_o(U_R_o_lo(1): U_R_o_hi(1),  U_R_o_lo(2): U_R_o_hi(2),   U_R_o_lo(3): U_R_o_hi(3), 0:n_uro-1) 
    real(rt), intent(inout) ::  U_R_n(U_R_n_lo(1): U_R_n_hi(1),  U_R_n_lo(2): U_R_n_hi(2),   U_R_n_lo(3): U_R_n_hi(3), 0:n_urn-1) 

    ! Temporary variables
    integer  :: i,j,k
    integer  :: ic,jc,kc
    integer  :: ii,id,ie,im,is
    real(rt) :: conv_dens, conv_mom, conv_enr, conv_ne

    ! Sanity check on size of arrays
    ! Note that we have set ngrow_thornado = ngrow_state in Castro_setup.cpp
    if (U_R_o_lo(1) .ne. S_lo(1) .or.  U_R_o_hi(1) .ne. S_hi(1) .or. &
        U_R_o_lo(2) .ne. S_lo(2) .or.  U_R_o_hi(2) .ne. S_hi(2) .or. &
        U_R_o_lo(3) .ne. S_lo(3) .or.  U_R_o_hi(3) .ne. S_hi(3)) then
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

    ! ************************************************************************************
    ! Copy from the Castro arrays into Thornado arrays from InitThornado_Patch
    ! ************************************************************************************
    do kc = lo(3)-ng,hi(3)+ng
    do jc = lo(2)-ng,hi(2)+ng
    do ic = lo(1)-ng,hi(1)+ng

         ! The uCF array was allocated in CreateFluidFieldsConserved with 
         !     ALLOCATE( uCF &
         !      (1:nDOFX, &
         !       1-swX(1):nX(1)+swX(1), &
         !       1-swX(2):nX(2)+swX(2), &
         !       1-swX(3):nX(3)+swX(3), &
         !       1:nCF) )

         i = ic+1
         j = jc+1
         k = kc+1

         ! Thornado uses units where c = G = k = 1, Meter = 1
         uCF(1:n_fluid_dof,i,j,k,iCF_D)  = S(ic,jc,kc,URHO)  * conv_dens
         uCF(1:n_fluid_dof,i,j,k,iCF_S1) = S(ic,jc,kc,UMX)   * conv_mom
         uCF(1:n_fluid_dof,i,j,k,iCF_S2) = S(ic,jc,kc,UMY)   * conv_mom
         uCF(1:n_fluid_dof,i,j,k,iCF_S3) = S(ic,jc,kc,UMZ)   * conv_mom
         uCF(1:n_fluid_dof,i,j,k,iCF_E)  = S(ic,jc,kc,UEDEN) * conv_enr
         uCF(1:n_fluid_dof,i,j,k,iCF_Ne) = S(ic,jc,kc,UFX)   * conv_ne

         ! The uCF array was allocated in CreatRadiationdFields_Conserved with 
         ! ALLOCATE &
         !   ( uCR(1:nDOF, &
         !         1-swE:nE+swE, &
         !         1-swX(1):nX(1)+swX(1), &
         !         1-swX(2):nX(2)+swX(2), &
         !         1-swX(3):nX(3)+swX(3), &
         !         1:nCR, 1:nSpecies) )

         do is = 1, nSpecies
         do im = 1, n_moments
         do ie = 1, nE
         do id = 1, nDOF
            ii   = (is-1)*(n_moments*nE*nDOF) + (im-1)*(nE*nDOF) + (ie-1)*nDOF + (id-1)
            uCR(id,ie,i,j,k,im,is) = U_R_o(ic,jc,kc,ii)
         end do
         end do
         end do
         end do

    end do
    end do
    end do

    ! ************************************************************************************
    ! Call the Fortran interface that lives in the thornado repo
    ! ************************************************************************************
    call Update_IMEX_PC2(dt, uCF, uCR)

    ! ************************************************************************************
    ! Copy back from the thornado arrays into Castro arrays
    ! ************************************************************************************
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)

         !! KS: if we fill UEINT, need the final fluid state from ComputeIncrement; is that uCF at this point?
         ! We store dS as a source term which we can add to S outside of this routine
         ! uCF returned as a cell-averaged quantity so all components are the same,
         !  can just use the first component
         ! Update_IMEX_PC2 doesn't currently change the fluid state
!         dS(i,j,k,URHO ) = uCF(1,i,j,k,iCF_D) - S(i,j,k,URHO)
!         dS(i,j,k,UMX  ) = uCF(1,i,j,k,iCF_S1) - S(i,j,k,UMX)
!         dS(i,j,k,UMY  ) = uCF(1,i,j,k,iCF_S2) - S(i,j,k,UMY)
!         dS(i,j,k,UMZ  ) = uCF(1,i,j,k,iCF_S3) - S(i,j,k,UMZ)
!         dS(i,j,k,UEDEN) = uCF(1,i,j,k,iCF_E) - S(i,j,k,UEDEN)
!         dS(i,j,k,UEINT) = ?
!         dS(i,j,k,UFX  ) = uCF(1,i,j,k,iCF_Ne) - S(i,j,k,UFX)

         do is = 1, nSpecies
         do im = 1, n_moments
         do ie = 1, nE
         do id = 1, nDOF
            ii   = (is-1)*(n_moments*nE*nDOF) + (im-1)*(nE*nDOF) + (ie-1)*nDOF + (id-1)
            U_R_n(i,j,k,ii) = uCR(id,ie,i,j,k,im,is) 
         end do
         end do
         end do
         end do

    end do
    end do
    end do

  end subroutine call_to_thornado
