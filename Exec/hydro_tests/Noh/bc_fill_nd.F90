module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

contains

  subroutine hypfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="hypfill")

    use amrex_constants_module, only: ONE, TWO
    use meth_params_module, only: NVAR, URHO, UTEMP, UMX, UMZ, UEDEN, UEINT, UFS
    use castro_util_module, only: position
    use eos_module, only: eos
    use extern_probin_module, only: eos_gamma
    use eos_type_module, only: eos_t, eos_input_rp
    use network, only: nspec

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer      :: i, j, k, n
    real(rt)     :: loc(3), vel(3), r
    type (eos_t) :: zone_state

    real(rt) :: pres_init = 1.0e-6_rt
    real(rt) :: rho_init = 1.0e0_rt

    ! Overwrite the outer boundary conditions

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (.not. (i .gt. domhi(1) .or. j .gt. domhi(2) .or. k .gt. domhi(3) )) cycle

             loc = position(i, j, k)

             r = sqrt( sum(loc**2) )

             zone_state % rho = rho_init * (ONE + time / r)**dble(AMREX_SPACEDIM-1)
             zone_state % P   = pres_init * (zone_state % rho / rho_init)**(ONE + eos_gamma)
             zone_state % xn  = ONE / nspec

             call eos(eos_input_rp, zone_state)

             ! Radial inflow with |v| = 1.

             vel(:) = -loc(:) / r

             adv(i,j,k,URHO)  = zone_state % rho
             adv(i,j,k,UTEMP) = zone_state % T
             adv(i,j,k,UEINT) = zone_state % e * zone_state % rho
             adv(i,j,k,UFS:UFS+nspec-1) = zone_state % xn(:) * zone_state % rho

             adv(i,j,k,UMX:UMZ) = adv(i,j,k,URHO) * vel(:)

             adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) + sum( adv(i,j,k,UMX:UMZ)**2 ) / ( TWO * adv(i,j,k,URHO) )

          enddo
       enddo
    enddo

  end subroutine hypfill



  subroutine denfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="denfill")

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

  end subroutine denfill

end module bc_fill_module
