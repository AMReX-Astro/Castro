module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_h1, &
                        domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_hypfill")

    use meth_params_module, only: NVAR  

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: adv_l1,adv_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: adv(adv_l1:adv_h1,NVAR)

    integer          :: n

    do n = 1,NVAR
       call filcc(adv(:,n),adv_l1,adv_h1, &
                  domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_denfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: adv_l1,adv_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: adv(adv_l1:adv_h1)

    call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_denfill



#ifdef GRAVITY
  subroutine ca_phigravfill(phi,phi_l1,phi_h1, &
                            domlo,domhi,delta,xlo,time,bc) &
                            bind(C, name="ca_phigravfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: phi(phi_l1:phi_h1)

    call filcc(phi,phi_l1,phi_h1, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill



  subroutine ca_gravxfill(grav,grav_l1,grav_h1, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravxfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: grav_l1,grav_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: grav(grav_l1:grav_h1)

    call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravxfill



  subroutine ca_gravyfill(grav,grav_l1,grav_h1, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravyfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: grav_l1,grav_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: grav(grav_l1:grav_h1)

    call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav,grav_l1,grav_h1, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravzfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: grav_l1,grav_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: grav(grav_l1:grav_h1)

    call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravzfill
#endif



#ifdef ROTATION
  subroutine ca_phirotfill(phi,phi_l1,phi_h1, &
                           domlo,domhi,delta,xlo,time,bc) &
                           bind(C, name="ca_phirotfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: phi(phi_l1:phi_h1)

    call filcc(phi,phi_l1,phi_h1, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_phirotfill



  subroutine ca_rotxfill(rot,rot_l1,rot_h1, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotxfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: rot_l1,rot_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: rot(rot_l1:rot_h1)

    call filcc(rot,rot_l1,rot_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_rotxfill



  subroutine ca_rotyfill(rot,rot_l1,rot_h1, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotyfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: rot_l1,rot_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: rot(rot_l1:rot_h1)

    call filcc(rot,rot_l1,rot_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_rotyfill



  subroutine ca_rotzfill(rot,rot_l1,rot_h1, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotzfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: rot_l1,rot_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: rot(rot_l1:rot_h1)

    call filcc(rot,rot_l1,rot_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_rotzfill
#endif
  

  
#ifdef REACTIONS  
  subroutine ca_reactfill(react,react_l1,react_h1, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_reactfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: react_l1,react_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    real(rt)         :: delta(1), xlo(1), time
    real(rt)         :: react(react_l1:react_h1)

    call filcc(react,react_l1,react_h1, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_reactfill
#endif
  


#ifdef RADIATION
  subroutine ca_radfill(rad,rad_l1,rad_h1, &
                        domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_radfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer rad_l1,rad_h1
    integer bc(1,2,*)
    integer domlo(1), domhi(1)
    real(rt)         delta(1), xlo(1), time
    real(rt)         rad(rad_l1:rad_h1)

    call filcc(rad,rad_l1,rad_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_radfill
#endif
  
end module bc_fill_module
