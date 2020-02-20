
  subroutine interpolate_fluid (lo, hi, &
                                S , s_lo, s_hi, ns , &
                                u0, u_lo, u_hi, nu, &
                                ng)

    use amrex_constants_module, only : fourth, half, zero, one, two
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEDEN,UFX
    use SubcellReconstructionModule, only: ReconstructionMatrix
    use FluidFieldsModule    , only : uCF, nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
    use EquationOfStateModule_TABLE, only : ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
                                            ComputeThermodynamicStates_Primitive_TABLE
    use UnitsModule, only : Gram, Centimeter, Second, AtomicMassUnit, Erg, Kelvin

    use ReferenceElementModuleX, only: NodesX_q

    implicit none
    integer, intent(in) :: lo(1), hi(1)
    integer, intent(in) :: s_lo(1),  s_hi(1)
    integer, intent(in) :: u_lo(3),  u_hi(3)
    integer, intent(in) :: ns, nu
    integer, intent(in) :: ng

    ! Conserved fluid state (rho, rho u, rho v, rho w, rho E, rho e...)
    real(rt), intent(inout) ::  S(s_lo(1):s_hi(1),ns)

    real(rt), intent(inout) :: u0(nu,u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nCF)

    u0(:,:,:,:,:) = zero

  end subroutine interpolate_fluid
