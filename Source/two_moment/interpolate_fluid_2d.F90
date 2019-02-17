
  subroutine interpolate_fluid (lo, hi, &
                                S , s_lo, s_hi, ns , &
                                u0, u_lo, u_hi, nu, &
                                ng) 

    use amrex_constants_module, only : fourth, half, zero, one, two
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : URHO,UMX,UMY,UMZ,UEDEN,UFX,UTEMP
    use FluidFieldsModule    , only : uCF, nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
    use RadiationFieldsModule, only : uCR
    use EquationOfStateModule_TABLE, only : ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
                                            ComputeThermodynamicStates_Primitive_TABLE
    use UnitsModule, only : Gram, Centimeter, Second, AtomicMassUnit, Erg, Kelvin

    use ReferenceElementModuleX, only: NodesX_q

    implicit none
    integer, intent(in) :: lo(2), hi(2)
    integer, intent(in) :: s_lo(2),  s_hi(2)
    integer, intent(in) :: u_lo(3),  u_hi(3)
    integer, intent(in) :: ns, nu
    integer, intent(in) :: ng

    ! Conserved fluid state (rho, rho u, rho v, rho w, rho E, rho e...)
    real(rt), intent(inout) ::  S(s_lo(1):s_hi(1),s_lo(2):s_hi(2),ns) 

    real(rt), intent(inout) :: u0(nu,u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nCF)

    ! Temporary variables
    integer  :: i,j,n
    integer  :: ic,jc,kc
    integer  :: ind
    real(rt) :: conv_dens, conv_mom, conv_enr, conv_ne

    ! For interpolation
    real(rt) :: xslope(ns), yslope(ns), xyslope(ns)
    real(rt) :: Sval(ns)
    real(rt) :: dlft(ns), drgt(ns), dcen(ns), dlim, dsgn

    real(rt) :: rho_in(1)
    real(rt) :: Ye_in(1)
    real(rt) :: E_in(1)
    real(rt) :: T_in(1)

    real(rt) :: T_out(1)
    real(rt) :: Ne_out(1)
    real(rt) :: Epervol_out(1)
    real(rt) :: Epermass_out(1)

    real(rt) :: fac(ns),fac_all
    real(rt) :: delta_S_slope(ns)
    real(rt) :: delta_S_val(ns)
    real(rt) :: S_ll(ns),S_lh(ns),S_hl(ns),S_hh(ns)

    real(rt) :: x(4)
    real(rt) :: y(4)
    real(rt) :: z(4)
    real(rt) :: xt, yt

    integer  :: interp_type


    ! No interpolation
    interp_type = 0 

    ! Use conserved variables for interpolation
    ! interp_type = 1 

    ! Use primitive variables (rho, E, Ye) for interpolation
    ! interp_type = 2 

    ! Use primitive variables (rho, T, Ye) for interpolation
    ! interp_type = 3 

    ! Use primitive variables (rho, T, Ye) for interpolation -- first fill the corners
    !     of the cell then use bilinear interpolation
    ! interp_type = 4

    if (s_lo(1) .gt. (lo(1)-ng-1) .or. s_hi(1) .lt. (hi(1)+ng+1) .or. &
        s_lo(2) .gt. (lo(2)-ng-1) .or. s_hi(2) .lt. (hi(2)+ng+1)) then
       print *,'WRONG NUMBER OF GHOST CELLS!'
       stop
    end if

    print *,"SIZE OF Castro   S   ", size(S,dim=1), size(S,dim=2), size(S,dim=3)
    print *,"SIZE OF Thornado uCF ", size(uCF,dim=1), size(uCF,dim=2), size(uCF,dim=3), &
                                     size(uCF,dim=4), size(uCF,dim=5)
    print *,"SIZE OF Thornado uCR ", size(uCR,dim=1), size(uCR,dim=2), size(uCR,dim=3), &
                                     size(uCR,dim=4), size(uCR,dim=5), size(uCR,dim=6), &
                                     size(uCR,dim=7)

    ! Conversion to thornado units
    conv_dens = Gram / Centimeter**3
    conv_mom  = Gram / Centimeter**2 / Second
    conv_enr  = Erg / Centimeter**3
    conv_ne   = 1.d0 / Centimeter**3

    ! These are the locations of the DG nodes in the space [-.5:.5]
    ! Note that in 2D with 2 nodes in each direction, they are ordered 
    ! ind = 1:  (Lo,Lo)
    ! ind = 2:  (Hi,Lo)
    ! ind = 3:  (Lo,Hi)
    ! ind = 4:  (Hi,Hi)
     do ind = 1, 4
        x(ind) = NodesX_q(1,ind)
        y(ind) = NodesX_q(2,ind)
        z(ind) = NodesX_q(3,ind)
     end do

     ! The uCF array was allocated in CreateFluidFieldsConserved with 
     !     ALLOCATE( uCF &
     !      (1:nDOFX, &
     !       1-swX(1):nX(1)+swX(1), &
     !       1-swX(2):nX(2)+swX(2), &
     !       1-swX(3):nX(3)+swX(3), &
     !       1:nCF) )

    ! ************************************************************************************
    ! No interpolation
    ! ************************************************************************************
    if (interp_type .eq. 0) then
       do jc = u_lo(2),u_hi(2)
       do ic = u_lo(1),u_hi(1)

          !   S spatial indices start at lo - (number of ghost zones)
          ! uCF spatial indices start at 1 - (number of ghost zones)
          i = lo(1) + 2*(ic-1)
          j = lo(2) + 2*(jc-1)

          ! In 2-d, kc = 1
          kc = 1

          ! Thornado uses units where c = G = k = 1, Meter = 1

          uCF(1,ic,jc,kc,iCF_D)  = S(i,j,URHO)  * conv_dens
          uCF(1,ic,jc,kc,iCF_S1) = S(i,j,UMX)   * conv_mom
          uCF(1,ic,jc,kc,iCF_S2) = S(i,j,UMY)   * conv_mom
          uCF(1,ic,jc,kc,iCF_S3) = S(i,j,UMZ)   * conv_mom
          uCF(1,ic,jc,kc,iCF_E)  = S(i,j,UEDEN) * conv_enr
          uCF(1,ic,jc,kc,iCF_Ne) = S(i,j,UFX)   * conv_ne

          uCF(2,ic,jc,kc,iCF_D)  = S(i+1,j,URHO)  * conv_dens
          uCF(2,ic,jc,kc,iCF_S1) = S(i+1,j,UMX)   * conv_mom
          uCF(2,ic,jc,kc,iCF_S2) = S(i+1,j,UMY)   * conv_mom
          uCF(2,ic,jc,kc,iCF_S3) = S(i+1,j,UMZ)   * conv_mom
          uCF(2,ic,jc,kc,iCF_E)  = S(i+1,j,UEDEN) * conv_enr
          uCF(2,ic,jc,kc,iCF_Ne) = S(i+1,j,UFX)   * conv_ne

          uCF(3,ic,jc,kc,iCF_D)  = S(i,j+1,URHO)  * conv_dens
          uCF(3,ic,jc,kc,iCF_S1) = S(i,j+1,UMX)   * conv_mom
          uCF(3,ic,jc,kc,iCF_S2) = S(i,j+1,UMY)   * conv_mom
          uCF(3,ic,jc,kc,iCF_S3) = S(i,j+1,UMZ)   * conv_mom
          uCF(3,ic,jc,kc,iCF_E)  = S(i,j+1,UEDEN) * conv_enr
          uCF(3,ic,jc,kc,iCF_Ne) = S(i,j+1,UFX)   * conv_ne

          uCF(4,ic,jc,kc,iCF_D)  = S(i+1,j+1,URHO)  * conv_dens
          uCF(4,ic,jc,kc,iCF_S1) = S(i+1,j+1,UMX)   * conv_mom
          uCF(4,ic,jc,kc,iCF_S2) = S(i+1,j+1,UMY)   * conv_mom
          uCF(4,ic,jc,kc,iCF_S3) = S(i+1,j+1,UMZ)   * conv_mom
          uCF(4,ic,jc,kc,iCF_E)  = S(i+1,j+1,UEDEN) * conv_enr
          uCF(4,ic,jc,kc,iCF_Ne) = S(i+1,j+1,UFX)   * conv_ne

          ! Make this copy so we can create dS on nodes, not after S has been
          !      averaged back to cell centers
          u0(1:4,ic,jc,kc,iCF_D ) = uCF(1:4,ic,jc,kc,iCF_D )
          u0(1:4,ic,jc,kc,iCF_S1) = uCF(1:4,ic,jc,kc,iCF_S1)
          u0(1:4,ic,jc,kc,iCF_S2) = uCF(1:4,ic,jc,kc,iCF_S2)
          u0(1:4,ic,jc,kc,iCF_S3) = uCF(1:4,ic,jc,kc,iCF_S3)
          u0(1:4,ic,jc,kc,iCF_E ) = uCF(1:4,ic,jc,kc,iCF_E )
          u0(1:4,ic,jc,kc,iCF_Ne) = uCF(1:4,ic,jc,kc,iCF_Ne)

       end do
       end do
    end if

  end subroutine interpolate_fluid
