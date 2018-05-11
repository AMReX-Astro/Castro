  subroutine call_to_thornado(lo, hi, dt, &
                              S, dS, s_lo, s_hi, n_fluid_comp, &
                              U_R_o, U_R_n, U_R_lo, U_R_hi, n_rad_comp, &
                              n_fluid_dof, n_energy, n_species, &
                              n_rad_dof, n_moments) &
                              bind(C, name="call_to_thornado")

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEINT,UEDEN,UFX
    use FluidFieldsModule, only : uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
    use RadiationFieldsModule, only : uCR
    use TimeSteppingModule_Castro, only : Update_IMEX_PC2

    implicit none
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) ::  U_R_lo(3),  U_R_hi(3)
    integer, intent(in) ::  n_fluid_comp, n_rad_comp
    integer, intent(in) ::  n_fluid_dof, n_energy, n_species, n_rad_dof, n_moments
    real(rt), intent(in) :: dt

    ! Here we expect  n_rad_comp = 20 x 16 x 6 x 4 (energy x dof x species x moments)

    ! Conserved fluid state (rho, rho u, rho v, rho w, rho E, rho e...)
    real(rt), intent(inout) ::  S(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),n_fluid_comp) 
    real(rt), intent(inout) :: dS(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),n_fluid_comp)

    ! Old and new radiation state
    real(rt), intent(inout) ::  U_R_o(U_R_lo(1): U_R_hi(1),  U_R_lo(2): U_R_hi(2),   U_R_lo(3): U_R_hi(3), n_rad_comp)
    real(rt), intent(inout) ::  U_R_n(U_R_lo(1): U_R_hi(1),  U_R_lo(2): U_R_hi(2),   U_R_lo(3): U_R_hi(3), n_rad_comp) 

    ! Temporary variables
    integer  :: ne,nn,ns,nm
    integer  :: i,j,k
    integer  :: ii,id,ie,im,is
    integer  :: ng

    ! Sanity check on size of arrays
    ! Note that we have set ngrow_thornado = ngrow_state in Castro_setup.cpp
    if (U_R_lo(1) .ne. S_lo(1) .or.  U_R_hi(1) .ne. S_hi(1) .or. &
        U_R_lo(2) .ne. S_lo(2) .or.  U_R_hi(2) .ne. S_hi(2) .or. &
        U_R_lo(3) .ne. S_lo(3) .or.  U_R_hi(3) .ne. S_hi(3)) then
        print *,'INCONSISTENT ARRAY BOUNDS ON FLUID AND RADIATION VARS'
        stop
    endif
    ! End Sanity check 

    !! KS: Do the corner ghost cells get passed too?
    ng = 2 ! 2 ghost zones for both fluid and radiation

    ! ************************************************************************************
    ! Copy from the Castro arrays into thornado arrays from InitThornado_Patch
    ! ************************************************************************************
    do k = lo(3)-ng,hi(3)+ng
    do j = lo(2)-ng,hi(2)+ng
    do i = lo(1)-ng,hi(1)+ng

         !! KS: need unit conversion from thornado variables to castro variables
         uCF(1:n_fluid_dof,i,j,k,iCF_D) = S(i,j,k,URHO)
         uCF(1:n_fluid_dof,i,j,k,iCF_S1) = S(i,j,k,UMX)
         uCF(1:n_fluid_dof,i,j,k,iCF_S2) = S(i,j,k,UMY)
         uCF(1:n_fluid_dof,i,j,k,iCF_S3) = S(i,j,k,UMZ)
         uCF(1:n_fluid_dof,i,j,k,iCF_E) = S(i,j,k,UEDEN)
         uCF(1:n_fluid_dof,i,j,k,iCF_Ne) = S(i,j,k,UFX) !! KS: Make sure that Ne is filled at some point and maintained

         do is = 1, n_species
         do im = 1, n_moments
         do ie = 1, n_energy
         do id = 1, n_rad_dof
            ii = (is-1)*(n_moments*n_energy*n_rad_dof) + (im-1)*(n_energy*n_rad_dof) + (ie-1)*n_rad_dof + (id-1)
            uCR(id,ie,i,j,k,im,is) = U_R_o(i,j,k,ii)
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

         !! KS: need unit conversion from thornado variables to castro variables
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

         do is = 1, n_species
         do im = 1, n_moments
         do ie = 1, n_energy
         do id = 1, n_rad_dof
            ii = (is-1)*(n_moments*n_energy*n_rad_dof) + (im-1)*(n_energy*n_rad_dof) + (ie-1)*n_rad_dof + (id-1)
            U_R_n(i,j,k,ii) = uCR(id,ie,i,j,k,im,is) 
         end do
         end do
         end do
         end do

    end do
    end do
    end do

  end subroutine call_to_thornado
