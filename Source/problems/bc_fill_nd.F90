module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine hypfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="hypfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2, NVAR)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    !$gpu

    call amrex_filccn(lo, hi, adv, adv_lo, adv_hi, NVAR, domlo, domhi, delta, xlo, bc)

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
