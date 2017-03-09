module bc_fill_module

  use bl_fort_module, only : rt => c_real
  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C)

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP
    use eos_module, only : gamma_const
    use probdata_module, only: p0_base

    use bl_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'
    
    integer :: adv_l1,adv_l2,adv_h1,adv_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer :: i,j,n
    real(rt)         :: y,pres

    do n = 1,NVAR
       call filcc(adv(adv_l1,adv_l2,n),adv_l1,adv_l2,adv_h1,adv_h2, &
                  domlo,domhi,delta,xlo,bc(1,1,n))
    enddo

    do n=1,NVAR
       !        XLO
       if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
          call bl_error('SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) ')
       end if

       !        XHI
       if ( bc(1,2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
          call bl_error('SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) ')
       end if

       !        YLO
       if ( bc(2,1,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
          do j=adv_l2,domlo(2)-1
             y = (j+0.5e0_rt)*delta(2)
             do i=adv_l1,adv_h1
                if (n .eq. URHO)  adv(i,j,n) = 1.0
                if (n .eq. UMX)   adv(i,j,n) = 0.0
                if (n .eq. UMY)   adv(i,j,n) = 0.0
                if (n .eq. UEDEN .or. n .eq. UEINT) then
                   pres = p0_base - y
                   adv(i,j,n) = pres / (gamma_const - 1.0e0_rt)
                end if

                if (n .eq. UFS)   adv(i,j,n) = 1.0
                if (n .eq. UTEMP) adv(i,j,n) = 0.0
             end do
          end do
       end if

       !        YHI
       if ( bc(2,2,n).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
          call bl_error('SHOULD NEVER GET HERE bc(2,2,n) .eq. EXT_DIR) ')
       end if

    end do

  end subroutine ca_hypfill

  ! ::: -----------------------------------------------------------

  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C)

    use bl_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'
    
    integer :: adv_l1,adv_l2,adv_h1,adv_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: adv(adv_l1:adv_h1,adv_l2:adv_h2)

    !     Note: this function should not be needed, technically, but is provided
    !     to filpatch because there are many times in the algorithm when just
    !     the density is needed.  We try to rig up the filling so that the same
    !     function is called here and in hypfill where all the states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    ! XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       call bl_error('SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) ')
    end if

    ! XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.lt.domhi(1)) then
       call bl_error('SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) ')
    end if

    ! YLO
    if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
       call bl_error('SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) ')
    end if

    ! YHI
    if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.gt.domhi(2)) then
       call bl_error('SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) ')
    end if

  end subroutine ca_denfill



  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    
    use bl_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravxfill



  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    
    use bl_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    
    use bl_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_h1,grav_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravzfill



  subroutine ca_phigravfill(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                            domlo,domhi,delta,xlo,time,bc) bind(C)

    use probdata_module
    
    use bl_fort_module, only : rt => c_real
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: phi_l1,phi_l2,phi_h1,phi_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2,domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill

end module bc_fill_module
