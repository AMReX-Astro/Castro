module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_filcc_module, only: filccn

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_hypfill")

    use meth_params_module, only: NVAR

    implicit none

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

    integer :: adv_lo(3), adv_hi(3)

    adv_lo = [adv_l1, adv_l2, adv_l3]
    adv_hi = [adv_h1, adv_h2, adv_h3]

    call filccn(adv_lo, adv_hi, adv, adv_lo, adv_hi, NVAR, domlo, domhi, delta, xlo, bc)

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_denfill")

    implicit none

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

    integer :: adv_lo(3), adv_hi(3)

    adv_lo = [adv_l1, adv_l2, adv_l3]
    adv_hi = [adv_h1, adv_h2, adv_h3]

    call filccn(adv_lo, adv_hi, adv, adv_lo, adv_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_denfill



#ifdef GRAVITY  
  subroutine ca_phigravfill(phi,phi_l1,phi_l2,phi_l3, &
                            phi_h1,phi_h2,phi_h3,domlo,domhi,delta,xlo,time,bc) &
                            bind(C, name="ca_phigravfill")

    
    implicit none

    integer,  intent(in   ) :: phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(in   ) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

    integer :: phi_lo(3), phi_hi(3)

    phi_lo = [phi_l1, phi_l2, phi_l3]
    phi_hi = [phi_h1, phi_h2, phi_h3]

    call filccn(phi_lo, phi_hi, phi, phi_lo, phi_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_phigravfill



  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravxfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(in   ) :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

    integer :: grav_lo(3), grav_hi(3)

    grav_lo = [grav_l1, grav_l2, grav_l3]
    grav_hi = [grav_h1, grav_h2, grav_h3]

    call filccn(grav_lo, grav_hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_gravxfill



  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravyfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(in   ) :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

    integer :: grav_lo(3), grav_hi(3)

    grav_lo = [grav_l1, grav_l2, grav_l3]
    grav_hi = [grav_h1, grav_h2, grav_h3]

    call filccn(grav_lo, grav_hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravzfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(in   ) :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

    integer :: grav_lo(3), grav_hi(3)

    grav_lo = [grav_l1, grav_l2, grav_l3]
    grav_hi = [grav_h1, grav_h2, grav_h3]

    call filccn(grav_lo, grav_hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_gravzfill
#endif



#ifdef ROTATION
  subroutine ca_phirotfill(phi,phi_l1,phi_l2,phi_l3, &
                           phi_h1,phi_h2,phi_h3,domlo,domhi,delta,xlo,time,bc) &
                           bind(C, name="ca_phirotfill")

    implicit none

    integer,  intent(in   ) :: phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(in   ) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

    integer :: phi_lo(3), phi_hi(3)

    phi_lo = [phi_l1, phi_l2, phi_l3]
    phi_hi = [phi_h1, phi_h2, phi_h3]

    call filccn(phi_lo, phi_hi, phi, phi_lo, phi_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_phirotfill



  subroutine ca_rotxfill(rot,rot_l1,rot_l2,rot_l3,rot_h1,rot_h2,rot_h3, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotxfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(in   ) :: rot(rot_l1:rot_h1,rot_l2:rot_h2,rot_l3:rot_h3)

    integer :: rot_lo(3), rot_hi(3)

    rot_lo = [rot_l1, rot_l2, rot_l3]
    rot_hi = [rot_h1, rot_h2, rot_h3]

    call filccn(rot_lo, rot_hi, rot, rot_lo, rot_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_rotxfill



  subroutine ca_rotyfill(rot,rot_l1,rot_l2,rot_l3,rot_h1,rot_h2,rot_h3, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotyfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(in   ) :: rot(rot_l1:rot_h1,rot_l2:rot_h2,rot_l3:rot_h3)

    integer :: rot_lo(3), rot_hi(3)

    rot_lo = [rot_l1, rot_l2, rot_l3]
    rot_hi = [rot_h1, rot_h2, rot_h3]

    call filccn(rot_lo, rot_hi, rot, rot_lo, rot_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_rotyfill



  subroutine ca_rotzfill(rot,rot_l1,rot_l2,rot_l3,rot_h1,rot_h2,rot_h3, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotzfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(in   ) :: rot(rot_l1:rot_h1,rot_l2:rot_h2,rot_l3:rot_h3)

    integer :: rot_lo(3), rot_hi(3)

    rot_lo = [rot_l1, rot_l2, rot_l3]
    rot_hi = [rot_h1, rot_h2, rot_h3]

    call filccn(rot_lo, rot_hi, rot, rot_lo, rot_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_rotzfill
#endif



#ifdef REACTIONS
  subroutine ca_reactfill(react,react_l1,react_l2,react_l3, &
                          react_h1,react_h2,react_h3,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_reactfill")

    implicit none

    integer,  intent(in   ) :: react_l1, react_l2, react_l3, react_h1, react_h2, react_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(in   ) :: react(react_l1:react_h1,react_l2:react_h2,react_l3:react_h3)

    integer :: react_lo(3), react_hi(3)

    react_lo = [react_l1, react_l2, react_l3]
    react_hi = [react_h1, react_h2, react_h3]

    call filccn(react_lo, react_hi, react, react_lo, react_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_reactfill
#endif



#ifdef RADIATION
  subroutine ca_radfill(rad,rad_l1,rad_l2,rad_l3, &
                        rad_h1,rad_h2,rad_h3,domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_radfill")

    implicit none

    integer,  intent(in   ) :: rad_l1, rad_l2, rad_l3, rad_h1, rad_h2, rad_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rad(rad_l1:rad_h1,rad_l2:rad_h2,rad_l3:rad_h3)

    integer :: rad_lo(3), rad_hi(3)

    rad_lo = [rad_l1, rad_l2, rad_l3]
    rad_hi = [rad_h1, rad_h2, rad_h3]

    call filccn(rad_lo, rad_hi, rad, rad_lo, rad_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_radfill
#endif

end module bc_fill_module
