  subroutine compute_thornado_timestep(dx, dt) &
                              bind(C, name="compute_thornado_timestep")

    use amrex_fort_module, only : rt => amrex_real
    use TimeSteppingModule_Castro, only : ComputeTimeStep_TwoMoment

    real(rt), intent(in ) :: dx(1)
    real(rt), intent(out) :: dt

    real(rt) :: dX_CGS(3)

    dX_CGS(1) = dx(1)
    dX_CGS(2:3)   = dX_CGS(1)

    call ComputeTimeStep_TwoMoment( dX_CGS, dt )

  end subroutine compute_thornado_timestep

  subroutine call_to_thornado(lo, hi, dt, &
                              S , s_lo, s_hi, ns , &
                              dS, d_lo, d_hi, nds, &
                              U_R_o, U_R_o_lo, U_R_o_hi, n_uro, &
                              dR   ,    dr_lo,    dr_hi, n_urn, &
                              n_moments, ng) &
                              bind(C, name="call_to_thornado")

    use amrex_constants_module, only : fourth, half, zero, one, two
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEINT,UEDEN,UFX,UFS
    use ProgramHeaderModule, only : nE, nNodesE, swE
    use FluidFieldsModule, only : uCF, nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
    use FluidFieldsModule, only : CreateFluidFields, DestroyFluidFields
    use GeometryFieldsModule, only : uGF
    use GeometryFieldsModuleE, only : uGE
    use RadiationFieldsModule, only : CreateRadiationFields,DestroyRadiationFields,nSpecies, uCR
    use TimeSteppingModule_Castro, only : Update_IMEX_PDARS
    use UnitsModule, only : Gram, Centimeter, Second, AtomicMassUnit, Erg

    use ReferenceElementModuleX, only: NodesX_q, WeightsX_q
    use SubcellReconstructionModule, only: ProjectionMatrix

    implicit none
    integer, intent(in) :: lo(1), hi(1)
    integer, intent(in) ::  s_lo(1),  s_hi(1)
    integer, intent(in) ::  d_lo(1),  d_hi(1)
    integer, intent(in) ::  U_R_o_lo(1),  U_R_o_hi(1)
    integer, intent(in) ::     dr_lo(1),     dr_hi(1)
    integer, intent(in) ::  ns, nds, n_uro, n_urn
    integer, intent(in) ::  n_moments
    integer, intent(in) :: ng
    real(rt), intent(in) :: dt

    ! Here we expect  n_rad_comp = 20 x 16 x 6 x 4 (energy x dof x species x moments)

    ! Conserved fluid state (rho, rho u, rho v, rho w, rho E, rho e...)
    real(rt), intent(inout) ::  S(s_lo(1):s_hi(1),ns)
    real(rt), intent(inout) :: dS(d_lo(1):d_hi(1),nds)

    ! Old and new radiation state
    real(rt), intent(inout) ::  U_R_o(U_R_o_lo(1): U_R_o_hi(1), 0:n_uro-1)
    real(rt), intent(inout) ::     dR(   dr_lo(1):    dr_hi(1), 0:n_urn-1)

    dS(:,:) = zero
    dR(:,:) = zero

  end subroutine call_to_thornado
