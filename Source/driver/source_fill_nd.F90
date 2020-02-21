module source_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NSRC

  implicit none

  include 'AMReX_bc_types.fi'

contains


  subroutine source_single_fill(lo, hi, state, s_lo, s_hi, domlo, domhi, delta, xlo, bc) bind(C, name="source_single_fill")
    ! Used for a source fill of any StateData.

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    integer :: d
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    ! handle an external BC via extrapolation here -- user's can override as needed
    bc_temp(:,:) = bc(:,:)

    do d = 1, AMREX_SPACEDIM
       if (bc(d,1) == EXT_DIR .and. s_lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       end if

       if (bc(d,2) == EXT_DIR .and. s_hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       end if
    end do

    call amrex_filccn(lo, hi, state, s_lo, s_hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine source_single_fill



  subroutine source_multi_fill(lo, hi, state, s_lo, s_hi, domlo, domhi, delta, xlo, bc) bind(C, name="source_multi_fill")

    use amrex_filcc_module, only: amrex_filccn
    
    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2, NSRC)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NSRC)

    integer :: d, n
    integer :: bc_temp(AMREX_SPACEDIM,2,NSRC)

    !$gpu

    ! handle an external BC via extrapolation here -- user's can override as needed
    bc_temp(:,:,:) = bc(:,:,:)

    do n = 1, NSRC
       do d = 1, AMREX_SPACEDIM
          if (bc(d,1,n) == EXT_DIR .and. s_lo(d) < domlo(d)) then
             bc_temp(d,1,n) = FOEXTRAP
          end if

          if (bc(d,2,n) == EXT_DIR .and. s_hi(d) > domhi(d)) then
             bc_temp(d,2,n) = FOEXTRAP
          end if
       end do
    end do

    call amrex_filccn(lo, hi, state, s_lo, s_hi, NSRC, domlo, domhi, delta, xlo, bc_temp)

  end subroutine source_multi_fill

end module source_fill_module
