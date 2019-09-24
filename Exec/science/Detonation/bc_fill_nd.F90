module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

contains

  subroutine hypfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="hypfill")

    use amrex_filcc_module, only: amrex_filccn
    use amrex_constants_module, only: HALF
    use meth_params_module, only: NVAR
    use probdata_module, only: fill_ambient_bc, fill_ambient
    use prob_params_module, only : problo
    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2, NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer  :: i, j, k, n
    real(rt) :: x

    !$gpu

    ! First, use the generic filling routines to make sure we have valid data everywhere
    ! on physical domain ghost cells.

    call amrex_filccn(lo, hi, adv, adv_lo, adv_hi, NVAR, domlo, domhi, delta, xlo, bc)

    ! Override the generic routine at the physical boundaries by
    ! setting the material to the ambient state. Note that we
    ! don't want to do this for interior/periodic boundaries,
    ! which have bc == 0, or for reflecting boundaries, which have
    ! bc == -1 or bc == 1.

    if (fill_ambient_bc) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                x = problo(1) + (dble(i) + HALF)*delta(1)

                if (AMREX_SPACEDIM .ge. 1) then
                   if (i .lt. domlo(1) .and. (bc(1,1,1) .ne. -1 .and. bc(1,1,1) .ne. 0 .and. bc(1,1,1) .ne. 1)) then
                      call fill_ambient(adv, adv_lo, adv_hi, i, j, k, x, time)
                   else if (i .gt. domhi(1) .and. (bc(1,2,1) .ne. -1 .and. bc(1,2,1) .ne. 0 .and. bc(1,2,1) .ne. 1)) then
                      call fill_ambient(adv, adv_lo, adv_hi, i, j, k, x, time)
                   end if
                end if

                if (AMREX_SPACEDIM .ge. 2) then
                   if (j .lt. domlo(2) .and. (bc(2,1,1) .ne. -1 .and. bc(2,1,1) .ne. 0 .and. bc(2,1,1) .ne. 1)) then
                      call fill_ambient(adv, adv_lo, adv_hi, i, j, k, x, time)
                   else if (j .gt. domhi(2) .and. (bc(2,2,1) .ne. -1 .and. bc(2,2,1) .ne. 0 .and. bc(2,2,1) .ne. 1)) then
                      call fill_ambient(adv, adv_lo, adv_hi, i, j, k, x, time)
                   end if
                end if
                
                if (AMREX_SPACEDIM .eq. 3) then
                   if (k .lt. domlo(3) .and. (bc(3,1,1) .ne. -1 .and. bc(3,1,1) .ne. 0 .and. bc(3,1,1) .ne. 1)) then
                      call fill_ambient(adv, adv_lo, adv_hi, i, j, k, x, time)
                   else if (k .gt. domhi(3) .and. (bc(3,2,1) .ne. -1 .and. bc(3,2,1) .ne. 0 .and. bc(3,2,1) .ne. 1)) then
                      call fill_ambient(adv, adv_lo, adv_hi, i, j, k, x, time)
                   endif
                endif

             enddo
          enddo
       enddo

    endif

  end subroutine hypfill



  subroutine denfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="denfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

    !$gpu

    call amrex_filccn(lo, hi, adv, adv_lo, adv_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine denfill



#ifdef GRAVITY
  subroutine gravxfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="gravxfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    !$gpu

    call amrex_filccn(lo, hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravxfill



  subroutine gravyfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="gravyfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    !$gpu

    call amrex_filccn(lo, hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravyfill



  subroutine gravzfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="gravzfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    !$gpu

    call amrex_filccn(lo, hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravzfill



  subroutine phigravfill(lo, hi, phi, phi_lo, phi_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="phigravfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

    !$gpu

    call amrex_filccn(lo, hi, phi, phi_lo, phi_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phigravfill
#endif

#ifdef REACTIONS
  subroutine reactfill(lo, hi, react, react_lo, react_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="reactfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: react_lo(3), react_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: react(react_lo(1):react_hi(1),react_lo(2):react_hi(2),react_lo(3):react_hi(3))

    !$gpu

    call amrex_filccn(lo, hi, react, react_lo, react_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine reactfill
#endif

end module bc_fill_module
