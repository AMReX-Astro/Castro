module bc_fill_module

  implicit none

  public

contains

  ! All subroutines in this file must be threadsafe because they are called
  ! inside OpenMP parallel regions.
  
  subroutine ca_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_hypfill")

    use meth_params_module, only: NVAR
    use prob_params_module, only: dim

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: adv_lo(3),adv_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

    integer          :: n

    do n = 1,NVAR
       call filcc_nd(adv(:,:,:,n),adv_lo,adv_hi,domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_denfill")

    use prob_params_module, only: dim  

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: adv_lo(3),adv_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

    call filcc_nd(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_denfill


  
#ifdef GRAVITY
  subroutine ca_phigravfill(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_phigravfill")

    use prob_params_module, only: dim  

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_lo(3),phi_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

    call filcc_nd(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill
  

  
  subroutine ca_gravxfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_gravxfill")

    use prob_params_module, only: dim  

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: grav_lo(3),grav_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    call filcc_nd(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravxfill



  subroutine ca_gravyfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_gravyfill")

    use prob_params_module, only: dim  

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: grav_lo(3),grav_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    call filcc_nd(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_gravzfill")

    use prob_params_module, only: dim

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: grav_lo(3),grav_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    call filcc_nd(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravzfill
#endif

  

#ifdef ROTATION
  subroutine ca_phirotfill(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_phirotfill")

    use prob_params_module, only: dim  

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_lo(3),phi_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

    call filcc_nd(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_phirotfill

  

  subroutine ca_rotxfill(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_rotxfill")

    use prob_params_module, only: dim  

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: rot_lo(3),rot_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))

    call filcc_nd(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_rotxfill



  subroutine ca_rotyfill(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_rotyfill")

    use prob_params_module, only: dim  

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: rot_lo(3),rot_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))

    call filcc_nd(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_rotyfill



  subroutine ca_rotzfill(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_rotzfill")

    use prob_params_module, only: dim

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: rot_lo(3),rot_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3))

    call filcc_nd(rot,rot_lo,rot_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_rotzfill
#endif  


#ifdef REACTIONS  
  subroutine ca_reactfill(react,react_lo,react_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_reactfill")

    use prob_params_module, only: dim  

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: react_lo(3),react_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: react(react_lo(1):react_hi(1),react_lo(2):react_hi(2),react_lo(3):react_hi(3))

    call filcc_nd(react,react_lo,react_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_reactfill
#endif



#ifdef RADIATION
  subroutine ca_radfill(rad,rad_lo,rad_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_radfill")

    use prob_params_module, only: dim  

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: rad_lo(3),rad_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: rad(rad_lo(1):rad_hi(1),rad_lo(2):rad_hi(2),rad_lo(3):rad_hi(3))

    call filcc_nd(rad,rad_lo,rad_hi,domlo,domhi,delta,xlo,bc)

  end subroutine ca_radfill
#endif
  
end module bc_fill_module
