  subroutine compute_thornado_timestep(dx, dt) &
                              bind(C, name="compute_thornado_timestep")

    use amrex_fort_module, only : rt => amrex_real
    use TimeSteppingModule_Castro, only : ComputeTimeStep_TwoMoment

    real(rt), intent(in ) :: dx(2)
    real(rt), intent(out) :: dt

    real(rt) :: dX_CGS(3)

    dX_CGS(1:2) = dx(1:2)
    dX_CGS(3)   = dX_CGS(2)

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
    use UnitsModule, only : Gram, Centimeter, Second, AtomicMassUnit, Erg

    use ReferenceElementModuleX, only: NodesX_q, WeightsX_q

    implicit none
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) ::  s_lo(2),  s_hi(2)
    integer, intent(in) ::  d_lo(2),  d_hi(2)
    integer, intent(in) ::  U_R_o_lo(2),  U_R_o_hi(2)
    integer, intent(in) ::     dr_lo(2),     dr_hi(2)
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
    real(rt), intent(inout) ::     dR(   dr_lo(1):    dr_hi(1),     dr_lo(2):    dr_hi(2), 0:n_urn-1) 

    ! Temporary variables
    integer  :: i,j,k,n
    integer  :: ic,jc
    integer  :: ii,id,ie,im,is,ind
    real(rt) :: x,y,z
    real(rt) :: conv_dens, conv_mom, conv_enr, conv_ne

    integer  :: nX(3)
    integer  :: swX(3)

    if (ng.ne. 2) &
      call amrex_abort("Need 2 ghost cells  n call_to_thornado!")

    conv_dens = Gram / Centimeter**3
    conv_mom  = Gram / Centimeter**2 / Second
    conv_enr  = Erg / Centimeter**3
    conv_ne   = 1.d0 / Centimeter**3

    nX(1) = hi(1) - lo(1) + 1
    nX(2) = hi(2) - lo(2) + 1
    nX(3) = 1

    swX(1) = ng
    swX(2) = ng
    swX(3) = 0

    call CreateFluidFields ( nX, swX, Verbose_Option = .FALSE. )

    call CreateRadiationFields ( nX, swX, nE, swE, nSpecies_Option = nSpecies, Verbose_Option = .FALSE. )

    ! ************************************************************************************
    ! Interpolate from the Castro "S" arrays into Thornado "uCF" arrays
    ! ************************************************************************************

    call interpolate_fluid (lo, hi, &
                            S , s_lo, s_hi, ns , &
                            n_fluid_dof, ng) 

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
    ! Call the time stepper that lives in the thornado repo
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
         ! We now use the weighting from thornado to convert the four node values back 
         !  into a single averaged value
         ! 
         ! Update_IMEX_PC2 doesn't currently change the fluid density or momentum
         ! 

         ! Zero out dS so we can accumulate weighted average in it
         dS(ic,jc,:) = 0.d0

         do ind = 1, n_fluid_dof
            dS(ic,jc,URHO ) = dS(ic,jc,URHO ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_D )
            dS(ic,jc,UMX  ) = dS(ic,jc,UMX  ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_S1)
            dS(ic,jc,UMY  ) = dS(ic,jc,UMY  ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_S2)
            dS(ic,jc,UEDEN) = dS(ic,jc,UEDEN) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_E )
            dS(ic,jc,UFX  ) = dS(ic,jc,UFX  ) + WeightsX_q(ind) * uCF(ind,i,j,k,iCF_Ne)
         end do

!        dS(ic,jc,URHO ) = dS(ic,jc,URHO ) / conv_dens - S(ic,jc,URHO )
!        dS(ic,jc,UMX  ) = dS(ic,jc,UMX  ) / conv_mom  - S(ic,jc,UMX  )
!        dS(ic,jc,UMY  ) = dS(ic,jc,UMY  ) / conv_mom  - S(ic,jc,UMY  )
         dS(ic,jc,UEDEN) = dS(ic,jc,UEDEN) / conv_enr  - S(ic,jc,UEDEN)
         dS(ic,jc,UFX  ) = dS(ic,jc,UFX  ) / conv_ne   - S(ic,jc,UFX  )

         dS(ic,jc,UEINT) = dS(ic,jc,UEDEN)     ! TRUE IFF NO MOMENTUM SOURCE TERMS

         do is = 1, nSpecies
         do im = 1, n_moments
         do ie = 1, nE
         do id = 1, nDOF

            ii   = (is-1)*(n_moments*nE*nDOF) + (im-1)*(nE*nDOF) + (ie-1)*nDOF + (id-1)
            if (im .eq. 1) dR(ic,jc,ii) = uCR(id,ie,i,j,k,im,is) - U_R_o(ic,jc,ii)
            if (im   >  1) dR(ic,jc,ii) = uCR(id,ie,i,j,k,im,is) - U_R_o(ic,jc,ii)

         end do
         end do
         end do
         end do

         ! Store electron molar fraction * density in the species
         ds(ic,jc,UFS) = dS(ic,jc,UFX) * AtomicMassUnit / Gram

    end do
    end do

    call DestroyFluidFields 
    call DestroyRadiationFields 

  end subroutine call_to_thornado
