module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine ambient_fill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, problo, dx, time, bc) bind(C, name="ambient_fill")

    use amrex_filcc_module, only: amrex_filccn
    use meth_params_module, only: fill_ambient_bc, URHO, NVAR
    use ambient_module, only: get_ambient_state
    use amrex_bc_types_module, only: amrex_bc_foextrap, amrex_bc_hoextrap

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2, NVAR)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ) :: problo(3), dx(3)
    real(rt), intent(in   ), value :: time

    integer :: i, j, k
    real(rt) :: ambient_state(NVAR), loc(3)

    !$gpu

    if (fill_ambient_bc == 1) then
       do k = lo(3), hi(3)
          loc(3) = problo(3) + dx(3) * (0.5_rt + k)
          do j = lo(2), hi(2)
             loc(2) = problo(2) + dx(2) * (0.5_rt + j)
             do i = lo(1), hi(1)
                loc(1) = problo(1) + dx(1) * (0.5_rt + i)
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
                   call get_ambient_state(ambient_state, loc, time)
                   adv(i,j,k,:) = ambient_state(:)
                end if
             end do
          end do
       end do
    end if

  end subroutine ambient_fill
  
end module bc_fill_module
