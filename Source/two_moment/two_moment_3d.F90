  subroutine call_to_thornado(lo, hi, dt, &
                              S , s_lo, s_hi, ns , &
                              dS, d_lo, d_hi, nds, &
                              U_R_o, U_R_o_lo, U_R_o_hi, n_uro, &
                              U_R_n, U_R_n_lo, U_R_n_hi, n_urn, &
                              n_fluid_dof, n_moments, ng) &
                              bind(C, name="call_to_thornado")

    use amrex_constants_module, only : fourth, half, zero, one, two
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
    real(rt), intent(inout) :: dS(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nds)

    ! Old and new radiation state
    real(rt), intent(inout) ::  U_R_o(U_R_o_lo(1): U_R_o_hi(1),  U_R_o_lo(2): U_R_o_hi(2),   U_R_o_lo(3): U_R_o_hi(3), 0:n_uro-1) 
    real(rt), intent(inout) ::  U_R_n(U_R_n_lo(1): U_R_n_hi(1),  U_R_n_lo(2): U_R_n_hi(2),   U_R_n_lo(3): U_R_n_hi(3), 0:n_urn-1) 

    ! Temporary variables
    integer  :: i,j,k,n
    integer  :: ic,jc,kc
    integer  :: ii,id,ie,im,is,ind
    real(rt) :: conv_dens, conv_mom, conv_enr, conv_ne, conv_J, conv_H, testdt

    ! For interpolation
    real(rt) :: x, y, z
    real(rt) :: xslope(ns), yslope(ns), zslope(ns)
    real(rt) :: Sval(ns)
    real(rt) :: dlft(ns), drgt(ns), dcen(ns), dlim, dsgn

    if (ng.ne. 2) &
      call amrex_abort("Need two ghost cells in call_to_thornado!")

    conv_dens = Gram / Centimeter**3
    conv_mom  = Gram / Centimeter**2 / Second
    conv_enr  = Gram / Centimeter / Second**2
    conv_ne   = 1.d0 / Centimeter**3
    conv_J    = Gram/Second**2/Centimeter
    conv_H    = Gram/Second**3

    ! ************************************************************************************
    ! Interpolate from the Castro "S" arrays into Thornado "uCF" arrays from InitThornado_Patch
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

         ! U_R_o spatial indices start at lo - (number of ghost zones)
         ! uCR spatial indices start at 1 - (number of ghost zones)
         i = ic - lo(1) + 1
         j = jc - lo(2) + 1
         k = kc - lo(3) + 1

        ! Define limited second-order slopes in x-direction
         dlft(:) = (S(ic  ,jc,kc,:) - S(ic-1,jc,kc,:)) * two
         drgt(:) = (S(ic+1,jc,kc,:) - S(ic  ,jc,kc,:)) * two
         dcen(:) = (S(ic+1,jc,kc,:) - S(ic-1,jc,kc,:)) * half

         do n = 1, ns
            dsgn = sign(one, dcen(n))
            xslope(n) = min( abs(dlft(n)), abs(drgt(n)) )
            if (dlft(n) * drgt(n) .ge. zero) then
               dlim = xslope(n)
            else
               dlim = zero
            endif
            xslope(n) = dsgn * min( dlim, abs(dcen(n)) )
         end do

        ! Define limited second-order slopes in y-direction
         dlft(:) = (S(ic,jc  ,kc,:) - S(ic,jc-1,kc,:)) * two
         drgt(:) = (S(ic,jc+1,kc,:) - S(ic,jc  ,kc,:)) * two
         dcen(:) = (S(ic,jc+1,kc,:) - S(ic,jc-1,kc,:)) * half

         do n = 1, ns
            dsgn = sign(one, dcen(n))
            yslope(n) = min( abs(dlft(n)), abs(drgt(n)) )
            if (dlft(n) * drgt(n) .ge. zero) then
               dlim = yslope(n)
            else
               dlim = zero
            endif
            yslope(n) = dsgn * min( dlim, abs(dcen(n)) )
         end do

        ! Define limited second-order slopes in z-direction
         dlft(:) = (S(ic,jc,kc  ,:) - S(ic,jc,kc-1,:)) * two
         drgt(:) = (S(ic,jc,kc+1,:) - S(ic,jc,kc  ,:)) * two
         dcen(:) = (S(ic,jc,kc+1,:) - S(ic,jc,kc-1,:)) * half

         do n = 1, ns
            dsgn = sign(one, dcen(n))
            zslope(n) = min( abs(dlft(n)), abs(drgt(n)) )
            if (dlft(n) * drgt(n) .ge. zero) then
               dlim = zslope(n)
            else
               dlim = zero
            endif
            zslope(n) = dsgn * min( dlim, abs(dcen(n)) )
         end do

         do ind = 1, n_fluid_dof

            ! These are the locations of the DG nodes in the space [-.5:.5]
            x = NodesX_q(1,ind)
            y = NodesX_q(2,ind)
            z = NodesX_q(3,ind)

            ! Use the slopes to extrapolate from the center to the nodes
            Sval(:) = S(ic,jc,kc,:) + x*xslope(:) + y*yslope(:) + z*zslope(:)

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
    end do

    ! ************************************************************************************
    ! Copy from the Castro U_R arrays into Thornado arrays from InitThornado_Patch
    ! ************************************************************************************
    do kc = lo(3)-ng,hi(3)+ng
    do jc = lo(2)-ng,hi(2)+ng
    do ic = lo(1)-ng,hi(1)+ng

         ! The uCR array was allocated in CreateRadiationdFields_Conserved with 
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
            if (im .eq. 1) uCR(id,ie,i,j,k,im,is) = U_R_o(ic,jc,kc,ii)
            if (im  >   1) uCR(id,ie,i,j,k,im,is) = U_R_o(ic,jc,kc,ii)
         end do
         end do
         end do
         end do

    end do
    end do
    end do

    ! ************************************************************************************
    ! Call the time stepper that lives in the thornado repo
    ! ************************************************************************************

    call Update_IMEX_PDARS(dt*Second, uCF, uCR)

    ! ************************************************************************************
    ! Copy back from the thornado arrays into Castro arrays
    ! ************************************************************************************
    do kc = lo(3),hi(3)
    do jc = lo(2),hi(2)
    do ic = lo(1),hi(1)

         ! uCR spatial indices start at 1 - ng
         ! U_R_n spatial indices start at lo
         i = ic - lo(1) + 1
         j = jc - lo(2) + 1
         k = kc - lo(3) + 1

         ! We store dS as a source term which we can add to S outside of this routine
         ! We now use the weighting from thornado to convert the four node values back
         !  into a single averaged value
         ! 
         ! Update_IMEX_PC2 doesn't currently change the fluid density or momentum
         ! 

         ! Zero out dS so we can accumulate weighted average in it
         dS(ic,jc,kc,:) = 0.d0

         do ind = 1, n_fluid_dof
             dS(ic,jc,kc,URHO ) = dS(ic,jc,kc,URHO ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_D ) 
             dS(ic,jc,kc,UMX  ) = dS(ic,jc,kc,UMX  ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_S1) 
             dS(ic,jc,kc,UMY  ) = dS(ic,jc,kc,UMY  ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_S2) 
             dS(ic,jc,kc,UMZ  ) = dS(ic,jc,kc,UMZ  ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_S3) 
             dS(ic,jc,kc,UEDEN) = dS(ic,jc,kc,UEDEN) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_E ) 
             dS(ic,jc,kc,UFX  ) = dS(ic,jc,kc,UFX  ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_Ne) 
         end do

!        dS(ic,jc,kc,URHO ) = dS(ic,jc,kc,URHO ) / conv_dens - S(ic,jc,kc,URHO )
!        dS(ic,jc,kc,UMX  ) = dS(ic,jc,kc,UMX  ) / conv_mom  - S(ic,jc,kc,UMX  )  
!        dS(ic,jc,kc,UMY  ) = dS(ic,jc,kc,UMY  ) / conv_mom  - S(ic,jc,kc,UMY  )  
!        dS(ic,jc,kc,UMZ  ) = dS(ic,jc,kc,UMZ  ) / conv_mom  - S(ic,jc,kc,UMZ  )  
         dS(ic,jc,kc,UEDEN) = dS(ic,jc,kc,UEDEN) / conv_enr  - S(ic,jc,kc,UEDEN)
         dS(ic,jc,kc,UFX  ) = dS(ic,jc,kc,UFX  ) / conv_ne   - S(ic,jc,kc,UFX  )  

         dS(ic,jc,kc,UEINT) = dS(ic,jc,kc,UEDEN)     ! TRUE IFF NO MOMENTUM SOURCE TERMS

         do is = 1, nSpecies
         do im = 1, n_moments
         do ie = 1, nE
         do id = 1, nDOF
            ii   = (is-1)*(n_moments*nE*nDOF) + (im-1)*(nE*nDOF) + (ie-1)*nDOF + (id-1)

            if (im .eq. 1) U_R_n(ic,jc,kc,ii) = uCR(id,ie,i,j,k,im,is)
            if (im   >  1) U_R_n(ic,jc,kc,ii) = uCR(id,ie,i,j,k,im,is)

         end do
         end do
         end do
         end do

    end do
    end do
    end do

  end subroutine call_to_thornado
