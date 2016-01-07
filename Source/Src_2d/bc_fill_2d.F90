module bc_fill_module

  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C)

    use meth_params_module, only: NVAR

    implicit none

    include 'bc_types.fi'

    integer          :: adv_l1,adv_l2,adv_h1,adv_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer          :: n

    do n = 1,NVAR
       call filcc(adv(:,:,n),adv_l1,adv_l2,adv_h1,adv_h2, &
            domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'bc_types.fi'

    integer          :: adv_l1,adv_l2,adv_h1,adv_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2)

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_denfill



#ifdef GRAVITY  
  subroutine ca_phigravfill(phi,phi_l1,phi_l2, &
                            phi_h1,phi_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_h1,phi_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill

  
  
  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'bc_types.fi'

    integer          :: grav_l1,grav_l2,grav_h1,grav_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravxfill



  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'bc_types.fi'

    integer          :: grav_l1,grav_l2,grav_h1,grav_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'bc_types.fi'

    integer          :: grav_l1,grav_l2,grav_h1,grav_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravzfill
#endif



#ifdef ROTATION
  subroutine ca_phirotfill(phi,phi_l1,phi_l2, &
                           phi_h1,phi_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_h1,phi_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_phirotfill

  
  
  subroutine ca_rotxfill(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
                         domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'bc_types.fi'

    integer          :: rot_l1,rot_l2,rot_h1,rot_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: rot(rot_l1:rot_h1,rot_l2:rot_h2)

    call filcc(rot,rot_l1,rot_l2,rot_h1,rot_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_rotxfill



  subroutine ca_rotyfill(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
                         domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'bc_types.fi'

    integer          :: rot_l1,rot_l2,rot_h1,rot_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: rot(rot_l1:rot_h1,rot_l2:rot_h2)

    call filcc(rot,rot_l1,rot_l2,rot_h1,rot_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_rotyfill



  subroutine ca_rotzfill(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
                         domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'bc_types.fi'

    integer          :: rot_l1,rot_l2,rot_h1,rot_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: rot(rot_l1:rot_h1,rot_l2:rot_h2)

    call filcc(rot,rot_l1,rot_l2,rot_h1,rot_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_rotzfill
#endif
  


#ifdef REACTIONS  
  subroutine ca_reactfill(react,react_l1,react_l2, &
                          react_h1,react_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'bc_types.fi'

    integer          :: react_l1,react_l2,react_h1,react_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: react(react_l1:react_h1,react_l2:react_h2)

    call filcc(react,react_l1,react_l2,react_h1,react_h2, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_reactfill
#endif



#ifdef RADIATION
  subroutine ca_radfill(rad,rad_l1,rad_l2, &
       rad_h1,rad_h2,domlo,domhi,delta,xlo,time,bc) bind(C)

    implicit none

    include 'bc_types.fi'

    integer :: rad_l1,rad_l2,rad_h1,rad_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    double precision delta(2), xlo(2), time
    double precision rad(rad_l1:rad_h1,rad_l2:rad_h2)

    call filcc(rad,rad_l1,rad_l2,rad_h1,rad_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_radfill
#endif
  
end module bc_fill_module
