  subroutine call_to_thornado(lo, hi, dt, &
                              S , s_lo, s_hi, ns , &
                              dS, d_lo, d_hi, nds, &
                              U_R_o, U_R_o_lo, U_R_o_hi, n_uro, &
                              U_R_n, U_R_n_lo, U_R_n_hi, n_urn, &
                              n_fluid_dof, n_moments, ng) &
                              bind(C, name="call_to_thornado")

    use amrex_constants_module, only : fourth, half, zero
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEINT,UEDEN,UFX
    use ProgramHeaderModule, only : nE, nDOF, nNodesX, nNodesE
    use FluidFieldsModule, only : uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
    use RadiationFieldsModule, only : nSpecies, uCR
    use TimeSteppingModule_Castro, only : Update_IMEX_PDARS
    use UnitsModule, only : Gram, Centimeter, Second

    use ReferenceElementModuleX, only: NodesX_q, WeightsX_q

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
    integer  :: ii,jj
    integer  :: ic,jc,kc
    integer  :: id,ie,im,is,ind
    real(rt) :: x,y,z
    real(rt) :: conv_dens, conv_mom, conv_enr, conv_ne, conv_J, conv_H, testdt

    ! For interpolation
    real(rt) :: xslope(ns), yslope(ns), xyslope(ns), Sval(ns)
    real(rt) :: min_val_ll(ns), min_val_lh(ns), min_val_hl(ns), min_val_hh(ns)
    real(rt) :: max_val_ll(ns), max_val_lh(ns), max_val_hl(ns), max_val_hh(ns)

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

    ! print *,'NDOF ',n_fluid_dof
    ! print *,'SIZE OF 1 NODES ',size(NodesX_q,dim=1)
    ! print *,'SIZE OF 2 NODES ',size(NodesX_q,dim=2)
    ! print *,'IN X ',NodesX_q(1,:)
    ! print *,'IN Y ',NodesX_q(2,:)
    ! print *,'IN Z ',NodesX_q(3,:)

    ! ************************************************************************************
    ! Copy from the Castro "S" arrays into Thornado "uCF" arrays from InitThornado_Patch
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

         !   S spatial indices start at lo - (number of ghost zones)
         ! uCF spatial indices start at 1 - (number of ghost zones)
         i = ic - lo(1) + 1
         j = jc - lo(2) + 1
         k = 1

         ! Define the slopes as centered differences
          xslope(:) = (S(ic+1,jc,:) - S(ic-1,jc,:)) * half
          yslope(:) = (S(ic,jc+1,:) - S(ic,jc-1,:)) * half
         xyslope(:) = (S(ic+1,jc+1,:) - S(ic-1,jc+1,:) &
                      +S(ic-1,jc-1,:) - S(ic+1,jc-1,:)) * fourth

         ! Find the min and max of the current state over the cell and 8 closest neighbors
         min_val_ll(:) = S(ic,jc,:)
         max_val_ll(:) = S(ic,jc,:)
         min_val_lh(:) = S(ic,jc,:)
         max_val_lh(:) = S(ic,jc,:)
         min_val_hl(:) = S(ic,jc,:)
         max_val_hl(:) = S(ic,jc,:)
         min_val_hh(:) = S(ic,jc,:)
         max_val_hh(:) = S(ic,jc,:)

         do jj = jc-1,jc
         do ii = ic-1,ic
            min_val_ll(:) = min(min_val_ll(:), S(ii,jj,:))
            max_val_ll(:) = max(max_val_ll(:), S(ii,jj,:))
         end do
         end do

         do jj = jc  ,jc+1
         do ii = ic-1,ic
            min_val_lh(:) = min(min_val_lh(:), S(ii,jj,:))
            max_val_lh(:) = max(max_val_lh(:), S(ii,jj,:))
         end do
         end do

         do jj = jc-1,jc
         do ii = ic  ,ic+1
            min_val_hl(:) = min(min_val_hl(:), S(ii,jj,:))
            max_val_hl(:) = max(max_val_hl(:), S(ii,jj,:))
         end do
         end do

         do jj = jc  ,jc+1
         do ii = ic  ,ic+1
            min_val_hh(:) = min(min_val_hh(:), S(ii,jj,:))
            max_val_hh(:) = max(max_val_hh(:), S(ii,jj,:))
         end do
         end do

         do ind = 1, n_fluid_dof

            ! These are the locations of the DG nodes in the space [-.5:.5]
            x = NodesX_q(1,ind)
            y = NodesX_q(2,ind)
            z = NodesX_q(3,ind)

            ! Use the slopes to extrapolate from the center to the nodes
            Sval(:) = S(ic,jc,:) + x*xslope(:) + y*yslope(:) + x*y*xyslope(:)

            ! Make sure that the extrapolation onto node locations creates no new max or min
            if (x.gt.zero .and. y.gt.zero) then
               Sval(:) = min(Sval(:),max_val_hh(:))
               Sval(:) = max(Sval(:),min_val_hh(:))
            end if
            if (x.gt.zero .and. y.le.zero) then
               Sval(:) = min(Sval(:),max_val_hl(:))
               Sval(:) = max(Sval(:),min_val_hl(:))
            end if
            if (x.le.zero .and. y.gt.zero) then
               Sval(:) = min(Sval(:),max_val_lh(:))
               Sval(:) = max(Sval(:),min_val_lh(:))
            end if
            if (x.le.zero .and. y.le.zero) then
               Sval(:) = min(Sval(:),max_val_ll(:))
               Sval(:) = max(Sval(:),min_val_ll(:))
            end if

            ! Thornado uses units where c = G = k = 1, Meter = 1
            uCF(ind,i,j,k,iCF_D)  = Sval(URHO)  * conv_dens
            uCF(ind,i,j,k,iCF_S1) = Sval(UMX)   * conv_mom
            uCF(ind,i,j,k,iCF_S2) = Sval(UMY)   * conv_mom
            uCF(ind,i,j,k,iCF_S3) = Sval(UMZ)   * conv_mom
            uCF(ind,i,j,k,iCF_E)  = Sval(UEDEN) * conv_enr
            uCF(ind,i,j,k,iCF_Ne) = Sval(UFX)   * conv_ne

         end do

    end do
    end do

    ! ************************************************************************************
    ! Copy from the Castro U_R arrays into Thornado arrays from InitThornado_Patch
    ! ************************************************************************************
    do jc = lo(2)-ng,hi(2)+ng
    do ic = lo(1)-ng,hi(1)+ng

         ! The uCR array was allocated in CreateRadiationdFields_Conserved with 
         ! ALLOCATE &
         !   ( uCR(1:nDOF, &
         !         1-swE:nE+swE, &
         !         1-swX(1):nX(1)+swX(1), &
         !         1-swX(2):nX(2)+swX(2), &
         !         1:nCR, 1:nSpecies) )

         ! U_R_o spatial indices start at lo - (number of ghost zones)
         !   uCR spatial indices start at 1 - (number of ghost zones)
         i = ic - lo(1) + 1
         j = jc - lo(2) + 1
         k = 1

         do is = 1, nSpecies
         do im = 1, n_moments
         do ie = 1, nE
         do id = 1, nDOF
            ii   = (is-1)*(n_moments*nE*nDOF) + (im-1)*(nE*nDOF) + (ie-1)*nDOF + (id-1)
            if (im .eq. 1) uCR(id,ie,i,j,k,im,is) = U_R_o(ic,jc,ii)
            if (im   >  1) uCR(id,ie,i,j,k,im,is) = U_R_o(ic,jc,ii)
         end do
         end do
         end do
         end do

    end do
    end do

    ! ************************************************************************************
    ! Call the Fortran interface that lives in the thornado repo
    ! ************************************************************************************

!    call Update_IMEX_PDARS(dt*Second, uCF, uCR)

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
         ! We now use the weighting from thornado to convert the four node values back 
         !  into a single averaged value
         ! 
         ! Update_IMEX_PC2 doesn't currently change the fluid density or momentum
         ! 

          dS(ic,jc,:) = 0.d0
          do ind = 1, n_fluid_dof
!            dS(ic,jc,URHO ) = dS(ic,jc,URHO ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_D )
!            dS(ic,jc,UMX  ) = dS(ic,jc,UMX  ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_S1)
!            dS(ic,jc,UMY  ) = dS(ic,jc,UMY  ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_S2)
             dS(ic,jc,UEDEN) = dS(ic,jc,UEDEN) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_E )
             dS(ic,jc,UFX  ) = dS(ic,jc,UFX  ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_Ne)
          end do

!         dS(ic,jc,URHO ) = dS(ic,jc,URHO ) / conv_dens - S(ic,jc,URHO )
!         dS(ic,jc,UMX  ) = dS(ic,jc,UMX  ) / conv_mom  - S(ic,jc,UMX  )
!         dS(ic,jc,UMY  ) = dS(ic,jc,UMY  ) / conv_mom  - S(ic,jc,UMY  )
          dS(ic,jc,UEDEN) = dS(ic,jc,UEDEN) / conv_enr  - S(ic,jc,UEDEN)
          dS(ic,jc,UFX  ) = dS(ic,jc,UFX  ) / conv_ne   - S(ic,jc,UFX  )

          dS(ic,jc,UEINT) = dS(ic,jc,UEDEN)     ! TRUE IFF NO MOMENTUM SOURCE TERMS

         do is = 1, nSpecies
         do im = 1, n_moments
         do ie = 1, nE
         do id = 1, nDOF
            ii   = (is-1)*(n_moments*nE*nDOF) + (im-1)*(nE*nDOF) + (ie-1)*nDOF + (id-1)
            if (im .eq. 1) U_R_n(ic,jc,ii) = uCR(id,ie,i,j,k,im,is)
            if (im   >  1) U_R_n(ic,jc,ii) = uCR(id,ie,i,j,k,im,is)

         end do
         end do
         end do
         end do

    end do
    end do

  end subroutine call_to_thornado
