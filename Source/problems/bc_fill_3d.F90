module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine hypfill(adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, &
                     domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

    integer  :: lo(3), hi(3)

    lo(1) = adv_l1
    lo(2) = adv_l2
    lo(3) = adv_l3
    hi(1) = adv_h1
    hi(2) = adv_h2
    hi(3) = adv_h3

    call amrex_filccn(lo, hi, adv, lo, hi, NVAR, domlo, domhi, delta, xlo, bc)

  end subroutine hypfill


  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")

    use meth_params_module, only: NVAR

    implicit none

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

    call hypfill(adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_hypfill


  subroutine denfill(adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, &
                     domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,1)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

    integer :: lo(3), hi(3)

    lo(1) = adv_l1
    lo(2) = adv_l2
    lo(3) = adv_l3
    hi(1) = adv_h1
    hi(2) = adv_h2
    hi(3) = adv_h3

    call amrex_filccn(lo, hi, adv, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine denfill


  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,1)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

    call denfill(adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_denfill


#ifdef GRAVITY
  subroutine phigravfill(phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                         domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

    integer :: lo(3), hi(3)

    lo(1) = phi_l1
    lo(2) = phi_l2
    lo(3) = phi_l3
    hi(1) = phi_h1
    hi(2) = phi_h2
    hi(3) = phi_h3

    call amrex_filccn(lo, hi, phi, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phigravfill


  subroutine ca_phigravfill(phi,phi_l1,phi_l2,phi_l3, &
                            phi_h1,phi_h2,phi_h3,domlo,domhi,delta,xlo,time,bc) &
                            bind(C, name="ca_phigravfill")

    
    implicit none

    integer,  intent(in   ) :: phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

    call phigravfill(phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                     domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_phigravfill


  subroutine gravxfill(grav, grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3, &
                       domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

    integer :: lo(3), hi(3)

    lo(1) = grav_l1
    lo(2) = grav_l2
    lo(3) = grav_l3
    hi(1) = grav_h1
    hi(2) = grav_h2
    hi(3) = grav_h3

    call amrex_filccn(lo, hi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravxfill


  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravxfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

    call gravxfill(grav, grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3, &
                   domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravxfill


  subroutine gravyfill(grav, grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3, &
                                    domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

    integer :: lo(3), hi(3)

    lo(1) = grav_l1
    lo(2) = grav_l2
    lo(3) = grav_l3
    hi(1) = grav_h1
    hi(2) = grav_h2
    hi(3) = grav_h3

    call amrex_filccn(lo, hi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravyfill


  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravyfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

    call gravyfill(grav, grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3, &
                   domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravyfill


  subroutine gravzfill(grav, grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3, &
                       domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

    integer :: lo(3), hi(3)

    lo(1) = grav_l1
    lo(2) = grav_l2
    lo(3) = grav_l3
    hi(1) = grav_h1
    hi(2) = grav_h2
    hi(3) = grav_h3

    call amrex_filccn(lo, hi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine gravzfill


  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravzfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

    call gravzfill(grav, grav_l1, grav_l2, grav_l3, grav_h1, grav_h2, grav_h3, &
                   domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravzfill
#endif



#ifdef ROTATION
  subroutine phirotfill(phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                        domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

    integer :: lo(3), hi(3)

    lo(1) = phi_l1
    lo(2) = phi_l2
    lo(3) = phi_l3
    hi(1) = phi_h1
    hi(2) = phi_h2
    hi(3) = phi_h3

    call amrex_filccn(lo, hi, phi, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phirotfill


  subroutine ca_phirotfill(phi,phi_l1,phi_l2,phi_l3, &
                           phi_h1,phi_h2,phi_h3,domlo,domhi,delta,xlo,time,bc) &
                           bind(C, name="ca_phirotfill")

    implicit none

    integer,  intent(in   ) :: phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

    call phirotfill(phi, phi_l1, phi_l2, phi_l3, phi_h1, phi_h2, phi_h3, &
                    domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_phirotfill


  subroutine rotxfill(rot, rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3, &
                      domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2,rot_l3:rot_h3)

    integer :: lo(3), hi(3)

    lo(1) = rot_l1
    lo(2) = rot_l2
    lo(3) = rot_l3
    hi(1) = rot_h1
    hi(2) = rot_h2
    hi(3) = rot_h3

    call amrex_filccn(lo, hi, rot, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine rotxfill


  subroutine ca_rotxfill(rot,rot_l1,rot_l2,rot_l3,rot_h1,rot_h2,rot_h3, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotxfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2,rot_l3:rot_h3)

    call rotxfill(rot, rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3, &
                  domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotxfill


  subroutine rotyfill(rot, rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3, &
                      domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2,rot_l3:rot_h3)

    integer :: lo(3), hi(3)

    lo(1) = rot_l1
    lo(2) = rot_l2
    lo(3) = rot_l3
    hi(1) = rot_h1
    hi(2) = rot_h2
    hi(3) = rot_h3

    call amrex_filccn(lo, hi, rot, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine rotyfill


  subroutine ca_rotyfill(rot,rot_l1,rot_l2,rot_l3,rot_h1,rot_h2,rot_h3, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotyfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2,rot_l3:rot_h3)

    call rotyfill(rot, rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3, &
                  domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotyfill


  subroutine rotzfill(rot, rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3, &
                      domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2,rot_l3:rot_h3)

    integer :: lo(3), hi(3)

    lo(1) = rot_l1
    lo(2) = rot_l2
    lo(3) = rot_l3
    hi(1) = rot_h1
    hi(2) = rot_h2
    hi(3) = rot_h3

    call amrex_filccn(lo, hi, rot, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine rotzfill


  subroutine ca_rotzfill(rot,rot_l1,rot_l2,rot_l3,rot_h1,rot_h2,rot_h3, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotzfill")

    implicit none

    integer,  intent(in   ) :: rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2,rot_l3:rot_h3)

    call rotzfill(rot, rot_l1, rot_l2, rot_l3, rot_h1, rot_h2, rot_h3, &
                  domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_rotzfill
#endif



#ifdef REACTIONS
  subroutine reactfill(react, react_l1, react_l2, react_l3, react_h1, react_h2, react_h3, &
                       domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: react_l1, react_l2, react_l3, react_h1, react_h2, react_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: react(react_l1:react_h1,react_l2:react_h2,react_l3:react_h3)

    integer :: lo(3), hi(3)

    lo(1) = react_l1
    lo(2) = react_l2
    lo(3) = react_l3
    hi(1) = react_h1
    hi(2) = react_h2
    hi(3) = react_h3

    call amrex_filccn(lo, hi, react, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine reactfill


  subroutine ca_reactfill(react,react_l1,react_l2,react_l3, &
                          react_h1,react_h2,react_h3,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_reactfill")

    implicit none

    integer,  intent(in   ) :: react_l1, react_l2, react_l3, react_h1, react_h2, react_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: react(react_l1:react_h1,react_l2:react_h2,react_l3:react_h3)

    call reactfill(react, react_l1, react_l2, react_l3, react_h1, react_h2, react_h3, &
                   domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_reactfill
#endif



#ifdef RADIATION
  subroutine radfill(rad, rad_l1, rad_l2, rad_l3, rad_h1, rad_h2, rad_h3, &
                     domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: rad_l1, rad_l2, rad_l3, rad_h1, rad_h2, rad_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rad(rad_l1:rad_h1,rad_l2:rad_h2,rad_l3:rad_h3)

    integer :: lo(3), hi(3)

    lo(1) = rad_l1
    lo(2) = rad_l2
    lo(3) = rad_l3
    hi(1) = rad_h1
    hi(2) = rad_h2
    hi(3) = rad_h3

    call amrex_filccn(lo, hi, rad, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine radfill


  subroutine ca_radfill(rad,rad_l1,rad_l2,rad_l3, &
                        rad_h1,rad_h2,rad_h3,domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_radfill")

    implicit none

    integer,  intent(in   ) :: rad_l1, rad_l2, rad_l3, rad_h1, rad_h2, rad_h3
    integer,  intent(in   ) :: bc(3,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: rad(rad_l1:rad_h1,rad_l2:rad_h2,rad_l3:rad_h3)

    call radfill(rad, rad_l1, rad_l2, rad_l3, rad_h1, rad_h2, rad_h3, &
                 domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_radfill
#endif

end module bc_fill_module
