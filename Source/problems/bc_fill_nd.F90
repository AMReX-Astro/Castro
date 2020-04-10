module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine hypfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="hypfill")

    use amrex_filcc_module, only: amrex_filccn
    use meth_params_module, only: fill_ambient_bc
    use ambient_module, only: ambient_state
    use amrex_bc_types_module, only: amrex_bc_foextrap, amrex_bc_hoextrap

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2, NVAR)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer :: i, j, k

    !$gpu

    if (fill_ambient_bc == 1) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if ((i < domlo(1) .and. (bc(1,1,1) == amrex_bc_foextrap .or. bc(1,1,1) == amrex_bc_hoextrap))  &
                    .or. (i > domhi(1) .and. (bc(1,2,1) == amrex_bc_foextrap .or. bc(1,2,1) == amrex_bc_hoextrap))  &
#if AMREX_SPACEDIM >= 2
                    .or. (j < domlo(2) .and. (bc(2,1,1) == amrex_bc_foextrap .or. bc(2,1,1) == amrex_bc_hoextrap)) &
                    .or. (j > domhi(2) .and. (bc(2,2,1) == amrex_bc_foextrap .or. bc(2,2,1) == amrex_bc_hoextrap)) &
#endif
#if AMREX_SPACEDIM == 3
                    .or. (k < domlo(3) .and. (bc(3,1,1) == amrex_bc_foextrap .or. bc(3,1,1) == amrex_bc_hoextrap)) &
                    .or. (k > domhi(3) .and. (bc(3,2,1) == amrex_bc_foextrap .or. bc(3,2,1) == amrex_bc_hoextrap)) &
#endif
                    ) then
                   adv(i,j,k,:) = ambient_state(:)
                end if
             end do
          end do
       end do
    end if

  end subroutine hypfill



  subroutine denfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="denfill")

    use amrex_filcc_module, only: amrex_filccn
    use meth_params_module, only: fill_ambient_bc, URHO
    use ambient_module, only: ambient_state
    use amrex_bc_types_module, only: amrex_bc_foextrap, amrex_bc_hoextrap

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: i, j, k

    !$gpu

    if (fill_ambient_bc == 1) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if ((i < domlo(1) .and. (bc(1,1) == amrex_bc_foextrap .or. bc(1,1) == amrex_bc_hoextrap)) &
                    .or. (i > domhi(1) .and. (bc(1,2) == amrex_bc_foextrap .or. bc(1,2) == amrex_bc_hoextrap)) &
#if AMREX_SPACEDIM >= 2
                    .or. (j < domlo(2) .and. (bc(2,1) == amrex_bc_foextrap .or. bc(2,1) == amrex_bc_hoextrap)) &
                    .or. (j > domhi(2) .and. (bc(2,2) == amrex_bc_foextrap .or. bc(2,2) == amrex_bc_hoextrap)) &
#endif
#if AMREX_SPACEDIM == 3
                    .or. (k < domlo(3) .and. (bc(3,1) == amrex_bc_foextrap .or. bc(3,1) == amrex_bc_hoextrap)) &
                    .or. (k > domhi(3) .and. (bc(3,2) == amrex_bc_foextrap .or. bc(3,2) == amrex_bc_hoextrap)) &
#endif
                   ) then
                   adv(i,j,k) = ambient_state(URHO)
                end if
             end do
          end do
       end do
    end if

  end subroutine denfill
  
end module bc_fill_module
