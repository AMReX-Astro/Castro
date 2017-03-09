module bc_fill_module

  use amrex_fort_module, only : rt => c_real
  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="ca_hypfill")

    use meth_params_module, only : NVAR
    use hse_bc_module

    use amrex_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer          :: adv_l1,adv_l2,adv_h1,adv_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer :: n

    do n=1,NVAR
       call filcc(adv(adv_l1,adv_l2,n),adv_l1,adv_l2,adv_h1,adv_h2, &
                  domlo,domhi,delta,xlo,bc(1,1,n))
    enddo

    ! YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
       call hse_bc_ylo(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                       domlo,domhi,delta,xlo,time,bc)
    end if

    ! YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
       call hse_bc_yhi(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                       domlo,domhi,delta,xlo,time,bc)
    end if

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) &
                        bind(C, name="ca_denfill")

    use probdata_module
    use hse_bc_module

    use amrex_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer          :: adv_l1,adv_l2,adv_h1,adv_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: adv(adv_l1:adv_h1,adv_l2:adv_h2)

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    !     YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
       call hse_bc_ylo(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                       domlo,domhi,delta,xlo,time,bc, density_only=.true.)
    end if

    !     YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
       call hse_bc_yhi(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                       domlo,domhi,delta,xlo,time,bc, density_only=.true.)
    end if

  end subroutine ca_denfill



  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravxfill")

    use probdata_module
    
    use amrex_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         grav(grav_l1:grav_h1,grav_l2:grav_h2)

    integer :: bc_temp(2,2)

    bc_temp(:,:) = bc(:,:,1)

    if ( bc(2,1,1).eq.EXT_DIR .and. grav_l2.lt.domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
               domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_gravxfill



  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravyfill")

    use probdata_module

    use amrex_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         grav(grav_l1:grav_h1,grav_l2:grav_h2)

    integer :: bc_temp(2,2)

    bc_temp(:,:) = bc(:,:,1)

    if ( bc(2,1,1).eq.EXT_DIR .and. grav_l2.lt.domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
               domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravzfill")

    use probdata_module
    
    use amrex_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         grav(grav_l1:grav_h1,grav_l2:grav_h2)

    integer :: bc_temp(2,2)

    bc_temp(:,:) = bc(:,:,1)

    if ( bc(2,1,1).eq.EXT_DIR .and. grav_l2.lt.domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
               domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_gravzfill



  subroutine ca_phigravfill(phi,phi_l1,phi_l2, &
                            phi_h1,phi_h2,domlo,domhi,delta,xlo,time,bc) &
                            bind(C, name="ca_phigravfill")

    use amrex_fort_module, only : rt => c_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_h1,phi_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    integer :: bc_temp(2,2)

    bc_temp(:,:) = bc(:,:,1)

    if ( bc(2,1,1).eq.EXT_DIR .and. phi_l2.lt.domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
               domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_phigravfill



  subroutine ca_rotxfill(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotxfill")

    use probdata_module
    
    use amrex_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: rot_l1,rot_l2,rot_h1,rot_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         rot(rot_l1:rot_h1,rot_l2:rot_h2)

    integer :: bc_temp(2,2)

    bc_temp(:,:) = bc(:,:,1)

    if ( bc(2,1,1).eq.EXT_DIR .and. rot_l2.lt.domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call filcc(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
               domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_rotxfill



  subroutine ca_rotyfill(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotyfill")

    use probdata_module

    use amrex_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: rot_l1,rot_l2,rot_h1,rot_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         rot(rot_l1:rot_h1,rot_l2:rot_h2)

    integer :: bc_temp(2,2)

    bc_temp(:,:) = bc(:,:,1)

    if ( bc(2,1,1).eq.EXT_DIR .and. rot_l2.lt.domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call filcc(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
               domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_rotyfill



  subroutine ca_rotzfill(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
                         domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="ca_rotzfill")

    use probdata_module
    
    use amrex_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: rot_l1,rot_l2,rot_h1,rot_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         rot(rot_l1:rot_h1,rot_l2:rot_h2)

    integer :: bc_temp(2,2)

    bc_temp(:,:) = bc(:,:,1)

    if ( bc(2,1,1).eq.EXT_DIR .and. rot_l2.lt.domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call filcc(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
               domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_rotzfill



  subroutine ca_phirotfill(phi,phi_l1,phi_l2, &
                           phi_h1,phi_h2,domlo,domhi,delta,xlo,time,bc) &
                           bind(C, name="ca_phirotfill")

    use amrex_fort_module, only : rt => c_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_h1,phi_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    integer :: bc_temp(2,2)

    bc_temp(:,:) = bc(:,:,1)

    if ( bc(2,1,1).eq.EXT_DIR .and. phi_l2.lt.domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
               domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_phirotfill



  subroutine ca_reactfill(react,react_l1,react_l2, &
                          react_h1,react_h2,domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_reactfill")

    use amrex_fort_module, only : rt => c_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: react_l1,react_l2,react_h1,react_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: react(react_l1:react_h1,react_l2:react_h2)

    integer :: bc_temp(2,2)

    bc_temp(:,:) = bc(:,:,1)

    if ( bc(2,1,1).eq.EXT_DIR .and. react_l2.lt.domlo(2)) then
       bc_temp(2,1) = FOEXTRAP
    endif

    call filcc(react,react_l1,react_l2,react_h1,react_h2, &
               domlo,domhi,delta,xlo,bc_temp)

  end subroutine ca_reactfill

end module bc_fill_module
