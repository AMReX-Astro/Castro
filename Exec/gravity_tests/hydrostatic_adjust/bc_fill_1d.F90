module bc_fill_module

  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_h1, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")

    use bl_constants_module
    use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UMY, UMZ, UTEMP, UEINT, UFS
    use probdata_module, only: hse_rho_top, hse_t_top, hse_X_top, &
         hse_eint_top, hse_p_top
    use network, only: nspec
    
    implicit none

    include 'bc_types.fi'

    integer          :: adv_l1,adv_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    double precision :: delta(1), xlo(1), time
    double precision :: adv(adv_l1:adv_h1,NVAR)

    integer n, i
    double precision :: vel

    ! call the generic ghostcell filling routine
    do n = 1,NVAR
       call filcc(adv(adv_l1,n), adv_l1,adv_h1, &
            domlo,domhi,delta,xlo,bc(1,1,n))
    enddo

    ! override the generic routine at the top physical boundary
    ! by resetting the velocity to zero there.
    if (adv_h1.gt.domhi(1)) then
       if (bc(1,2,UMX).eq.FOEXTRAP) then
          do i = domhi(1)+1,adv_h1
             !adv(i,UMX) = adv(domhi(1),UMX)
             vel = max(adv(i,UMX)/adv(i,URHO),0.d0)
             adv(i,URHO)  = hse_rho_top
             adv(i,UMX)   = adv(i,URHO)*vel
             adv(i,UMY)   = ZERO
             adv(i,UMZ)   = ZERO
             adv(i,UTEMP) = hse_T_top
             adv(i,UEINT) = hse_rho_top*hse_eint_top
             adv(i,UEDEN) = hse_rho_top*hse_eint_top + &
                  0.5*adv(i,UMX)**2/adv(i,URHO)
             adv(i,UFS:UFS+nspec-1) = hse_rho_top*hse_X_top(:)
          enddo
       end if
    end if

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_h1, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")
    
    implicit none

    integer          :: adv_l1,adv_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    double precision :: delta(1), xlo(1), time
    double precision :: adv(adv_l1:adv_h1)
    logical          :: rho_only

    call filcc(adv,adv_l1,adv_h1, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_denfill



#ifdef GRAVITY
  subroutine ca_gravxfill(grav,grav_l1,grav_h1, &
                          domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_gravxfill")

    use probdata_module
    
    implicit none

    integer          :: grav_l1,grav_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    double precision :: delta(1), xlo(1), time
    double precision :: grav(grav_l1:grav_h1)

    call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravxfill



  subroutine ca_gravyfill(grav,grav_l1,grav_h1, &
       domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_gravyfill")

    use probdata_module
    
    implicit none

    integer          :: grav_l1,grav_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    double precision :: delta(1), xlo(1), time
    double precision :: grav(grav_l1:grav_h1)

    call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav,grav_l1,grav_h1, &
                          domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_gravzfill")

    use probdata_module
    
    implicit none

    integer          :: grav_l1,grav_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    double precision :: delta(1), xlo(1), time
    double precision :: grav(grav_l1:grav_h1)

    call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravzfill



  subroutine ca_phigravfill(phi,phi_l1,phi_h1, &
                            domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_phigravfill")

    implicit none

    include 'bc_types.fi'

    integer          :: phi_l1,phi_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    double precision :: delta(1), xlo(1), time
    double precision :: phi(phi_l1:phi_h1)

    call filcc(phi,phi_l1,phi_h1, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill
#endif


#ifdef REACTIONS
  subroutine ca_reactfill(adv,adv_l1,adv_h1, &
                          domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_reactfill")

    implicit none

    integer          :: adv_l1,adv_h1
    integer          :: bc(1,2,*)
    integer          :: domlo(1), domhi(1)
    double precision :: delta(1), xlo(1), time
    double precision :: adv(adv_l1:adv_h1)
    logical          :: rho_only

    call filcc(adv,adv_l1,adv_h1, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_reactfill
#endif

end module bc_fill_module
