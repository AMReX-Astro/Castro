  subroutine compute_thornado_timestep(dx, dt) &
                              bind(C, name="compute_thornado_timestep")

    use amrex_fort_module, only : rt => amrex_real
    use TimeSteppingModule_Castro, only : ComputeTimeStep_TwoMoment

    real(rt), intent(in ) :: dx(3)
    real(rt), intent(out) :: dt

    real(rt) :: dX_CGS(3)

    dX_CGS(:) = dx(:)

    call ComputeTimeStep_TwoMoment( dX_CGS, dt )

  end subroutine compute_thornado_timestep

  subroutine call_to_thornado(lo, hi, dt, &
                              S , s_lo, s_hi, ns , &
                              dS, d_lo, d_hi, nds, &
                              U_R_o, U_R_o_lo, U_R_o_hi, n_uro, &
                              dR   ,    dr_lo,    dr_hi, n_urn, &
                              n_fluid_dof, n_moments, ng) &
                              bind(C, name="call_to_thornado")

    use amrex_constants_module, only : fourth, half, zero, one, two
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEINT,UEDEN,UFX,UFS
    use ProgramHeaderModule, only : nE, nDOF, nNodesX, nNodesE, swE
    use FluidFieldsModule, only : uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
    use FluidFieldsModule, only : CreateFluidFields, DestroyFluidFields
    use RadiationFieldsModule, only : CreateRadiationFields,DestroyRadiationFields,nSpecies, uCR
    use TimeSteppingModule_Castro, only : Update_IMEX_PDARS
    use UnitsModule, only : Gram, Centimeter, Second,  AtomicMassUnit, Erg

    use ReferenceElementModuleX, only: NodesX_q, WeightsX_q

    implicit none
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) ::  d_lo(3),  d_hi(3)
    integer, intent(in) ::  U_R_o_lo(3),  U_R_o_hi(3)
    integer, intent(in) ::     dr_lo(3),     dr_hi(3)
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
    real(rt), intent(inout) ::     dR(   dr_lo(1):    dr_hi(1),     dr_lo(2):    dr_hi(2),      dr_lo(3):    dr_hi(3), 0:n_urn-1)

    ! Temporary variables
    integer  :: i,j,k,n
    integer  :: ic,jc,kc
    integer  :: ii,id,ie,im,is,ind
    real(rt) :: conv_dens, conv_mom, conv_enr, conv_ne

    integer  :: nX(3)
    integer  :: swX(3)

    if (ng.ne. 2) &
      call amrex_abort("Need two ghost cells in call_to_thornado!")

    conv_dens = Gram / Centimeter**3
    conv_mom  = Gram / Centimeter**2 / Second
    conv_enr  = Gram / Centimeter / Second**2
    conv_enr  = Erg / Centimeter**3
    conv_ne   = 1.d0 / Centimeter**3

    nX(:)  = hi(:) - lo(:) + 1
    swX(:) = ng

    print *,'NX ', nx(:)

    call CreateFluidFields ( nX, swX )

    call CreateRadiationFields ( nX, swX, nE, swE, nSpecies_Option = nSpecies )

    ! ************************************************************************************
    ! Interpolate from the Castro "S" arrays into Thornado "uCF" arrays
    ! ************************************************************************************

    call interpolate_fluid (lo, hi, &
                            S , s_lo, s_hi, ns , &
                            n_fluid_dof, ng)

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

         ! U_R_o spatial indices start at lo - (number of ghost zones)
         !   uCR spatial indices start at 1 - (number of ghost zones)
         i = ic - lo(1) + 1
         j = jc - lo(2) + 1
         k = kc - lo(3) + 1

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

            if (im .eq. 1) dR(ic,jc,kc,ii) = uCR(id,ie,i,j,k,im,is) - U_R_o(ic,jc,kc,ii)
            if (im   >  1) dR(ic,jc,kc,ii) = uCR(id,ie,i,j,k,im,is) - U_R_o(ic,jc,kc,ii)

         end do
         end do
         end do
         end do

         ! Store electron molar fraction * density in the species
         ds(ic,jc,kc,UFS) = dS(ic,jc,kc,UFX) * AtomicMassUnit / Gram

    end do
    end do
    end do

    call DestroyFluidFields
    call DestroyRadiationFields

  end subroutine call_to_thornado
