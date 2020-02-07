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

    call amrex_filccn(lo, hi, adv, adv_lo, adv_hi, NVAR, domlo, domhi, delta, xlo, bc)

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

    call amrex_filccn(lo, hi, adv, adv_lo, adv_hi, 1, domlo, domhi, delta, xlo, bc)

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


  
#ifdef GRAVITY
  subroutine phigravfill(lo, hi, phi, phi_lo, phi_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="phigravfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ), value :: time

    !$gpu

    call amrex_filccn(lo, hi, phi, phi_lo, phi_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phigravfill

  

  subroutine gravxfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="gravxfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: d
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    do d = 1, AMREX_SPACEDIM
       if (bc(d,1) == EXT_DIR .and. grav_lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       end if

       if (bc(d,2) == EXT_DIR .and. grav_hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       end if
    end do

    call amrex_filccn(lo, hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine gravxfill



  subroutine gravyfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="gravyfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: d
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    do d = 1, AMREX_SPACEDIM
       if (bc(d,1) == EXT_DIR .and. grav_lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       end if

       if (bc(d,2) == EXT_DIR .and. grav_hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       end if
    end do

    call amrex_filccn(lo, hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine gravyfill



  subroutine gravzfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="gravzfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: d
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    do d = 1, AMREX_SPACEDIM
       if (bc(d,1) == EXT_DIR .and. grav_lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       end if

       if (bc(d,2) == EXT_DIR .and. grav_hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       end if
    end do

    call amrex_filccn(lo, hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine gravzfill
#endif

  

#ifdef ROTATION
  subroutine phirotfill(lo, hi, phi, phi_lo, phi_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="phirotfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: d
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    do d = 1, AMREX_SPACEDIM
       if (bc(d,1) == EXT_DIR .and. phi_lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       end if

       if (bc(d,2) == EXT_DIR .and. phi_hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       end if
    end do

    call amrex_filccn(lo, hi, phi, phi_lo, phi_hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine phirotfill

  

  subroutine rotxfill(lo, hi, rot, rot_lo, rot_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="rotxfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: rot_lo(3), rot_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: d
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    do d = 1, AMREX_SPACEDIM
       if (bc(d,1) == EXT_DIR .and. rot_lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       end if

       if (bc(d,2) == EXT_DIR .and. rot_hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       end if
    end do

    call amrex_filccn(lo, hi, rot, rot_lo, rot_hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine rotxfill



  subroutine rotyfill(lo, hi, rot, rot_lo, rot_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="rotyfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: rot_lo(3), rot_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: d
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    do d = 1, AMREX_SPACEDIM
       if (bc(d,1) == EXT_DIR .and. rot_lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       end if

       if (bc(d,2) == EXT_DIR .and. rot_hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       end if
    end do

    call amrex_filccn(lo, hi, rot, rot_lo, rot_hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine rotyfill



  subroutine rotzfill(lo, hi, rot, rot_lo, rot_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="rotzfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: rot_lo(3), rot_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: d
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    do d = 1, AMREX_SPACEDIM
       if (bc(d,1) == EXT_DIR .and. rot_lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       end if

       if (bc(d,2) == EXT_DIR .and. rot_hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       end if
    end do

    call amrex_filccn(lo, hi, rot, rot_lo, rot_hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine rotzfill
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
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: react(react_lo(1):react_hi(1),react_lo(2):react_hi(2),react_lo(3):react_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: d
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    do d = 1, AMREX_SPACEDIM
       if (bc(d,1) == EXT_DIR .and. react_lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       end if

       if (bc(d,2) == EXT_DIR .and. react_hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       end if
    end do

    call amrex_filccn(lo, hi, react, react_lo, react_hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine reactfill
#endif



#ifdef RADIATION
  subroutine radfill(lo, hi, rad, rad_lo, rad_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="radfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: rad_lo(3), rad_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: rad(rad_lo(1):rad_hi(1),rad_lo(2):rad_hi(2),rad_lo(3):rad_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: d
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    ! handle an external BC via extrapolation here 
    bc_temp(:,:) = bc(:,:)

    do d = 1, AMREX_SPACEDIM
       if (bc(d,1) == EXT_DIR .and. rad_lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       end if

       if (bc(d,2) == EXT_DIR .and. rad_hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       end if
    end do

    call amrex_filccn(lo, hi, rad, rad_lo, rad_hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine radfill
#endif
  
end module bc_fill_module
