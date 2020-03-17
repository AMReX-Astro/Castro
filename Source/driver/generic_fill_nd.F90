module generic_fill_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  include 'AMReX_bc_types.fi'

contains

  subroutine generic_fill(lo, hi, state, s_lo, s_hi, numcomp, domlo, domhi, delta, xlo, bc) bind(C, name="generic_fill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2, numcomp)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),numcomp)
    integer,  intent(in   ), value :: numcomp

    integer :: d, n
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    do n = 1, numcomp

       ! Handle an external BC via extrapolation here -- users can override as needed
       bc_temp(:,:) = bc(:,:,n)

       do d = 1, AMREX_SPACEDIM
          if (bc(d,1,n) == EXT_DIR .and. s_lo(d) < domlo(d)) then
             bc_temp(d,1) = FOEXTRAP
          end if

          if (bc(d,2,n) == EXT_DIR .and. s_hi(d) > domhi(d)) then
             bc_temp(d,2) = FOEXTRAP
          end if
       end do

       call amrex_filccn(lo, hi, state(:,:,:,n), s_lo, s_hi, 1, domlo, domhi, delta, xlo, bc_temp)

    end do

  end subroutine generic_fill

end module generic_fill_module
