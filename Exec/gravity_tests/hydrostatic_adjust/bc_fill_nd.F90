module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine hypfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="hypfill")

    use amrex_constants_module, only: ZERO, HALF
    use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UMY, UMZ, UTEMP, UEINT, UFS
    use probdata_module, only: hse_rho_top, hse_t_top, hse_X_top, &
                               hse_eint_top, hse_p_top
    use network, only: nspec

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2,NVAR)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer :: n, i, j, k
    real(rt) :: vel

    ! override the generic routine at the top physical boundary
    ! by resetting the velocity to zero there.

    if (AMREX_SPACEDIM == 1) then

       if (hi(1) > domhi(1)) then
          if (bc(1,2,UMX) == FOEXTRAP) then
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = domhi(1)+1, hi(1)

                      !adv(i,UMX) = adv(domhi(1),UMX)
                      vel = max(adv(i,j,k,UMX)/adv(i,j,k,URHO), ZERO)
                      adv(i,j,k,URHO)  = hse_rho_top
                      adv(i,j,k,UMX)   = adv(i,j,k,URHO)*vel
                      adv(i,j,k,UMY)   = ZERO
                      adv(i,j,k,UMZ)   = ZERO
                      adv(i,j,k,UTEMP) = hse_T_top
                      adv(i,j,k,UEINT) = hse_rho_top*hse_eint_top
                      adv(i,j,k,UEDEN) = hse_rho_top*hse_eint_top + &
                           HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                      adv(i,j,k,UFS:UFS+nspec-1) = hse_rho_top*hse_X_top(:)
                   end do
                end do
             end do
          end if
       end if

    else if (AMREX_SPACEDIM == 2) then

       if (hi(2) > domhi(2)) then
          if (bc(2,2,UMY) == FOEXTRAP) then
             do k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   do j = domhi(2)+1, hi(2)

                      !adv(i,UMX) = adv(domhi(1),UMX)
                      vel = max(adv(i,j,k,UMY)/adv(i,j,k,URHO), ZERO)
                      adv(i,j,k,URHO)  = hse_rho_top
                      adv(i,j,k,UMX)   = ZERO
                      adv(i,j,k,UMY)   = adv(i,j,k,URHO)*vel
                      adv(i,j,k,UMZ)   = ZERO
                      adv(i,j,k,UTEMP) = hse_T_top
                      adv(i,j,k,UEINT) = hse_rho_top*hse_eint_top
                      adv(i,j,k,UEDEN) = hse_rho_top*hse_eint_top + &
                           HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                      adv(i,j,k,UFS:UFS+nspec-1) = hse_rho_top*hse_X_top(:)
                   end do
                end do
             end do
          end if
       end if

    else

       if (hi(3) > domhi(3)) then
          if (bc(3,2,UMZ) == FOEXTRAP) then
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   do k = domhi(3)+1, hi(3)

                      !adv(i,UMX) = adv(domhi(1),UMX)
                      vel = max(adv(i,j,k,UMZ)/adv(i,j,k,URHO), ZERO)
                      adv(i,j,k,URHO)  = hse_rho_top
                      adv(i,j,k,UMX)   = ZERO
                      adv(i,j,k,UMY)   = ZERO
                      adv(i,j,k,UMZ)   = adv(i,j,k,URHO)*vel
                      adv(i,j,k,UTEMP) = hse_T_top
                      adv(i,j,k,UEINT) = hse_rho_top*hse_eint_top
                      adv(i,j,k,UEDEN) = hse_rho_top*hse_eint_top + &
                           HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                      adv(i,j,k,UFS:UFS+nspec-1) = hse_rho_top*hse_X_top(:)
                   end do
                end do
             end do
          end if
       end if

    end if

  end subroutine hypfill



  subroutine denfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="denfill")

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

  end subroutine denfill

end module bc_fill_module
