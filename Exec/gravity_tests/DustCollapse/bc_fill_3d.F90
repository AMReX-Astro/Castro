module bc_fill_module

  implicit none

  public

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                        domlo,domhi,delta,xlo,time,bc) bind(C,name="ca_hypfill")

    use bl_error_module
    use meth_params_module, only : NVAR,UMX,UMY,UMZ
    use prob_params_module, only : center

    implicit none
    
    include 'AMReX_bc_types.fi'
    
    integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
    integer bc(3,2,*)
    integer domlo(3), domhi(3)
    double precision delta(3), xlo(3), time
    double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

    integer          :: i, j, k, n
    integer          :: ic,jc,kc
    double precision :: x,y,z,r
    double precision :: xc,yc,zc,rc
    double precision :: mom,momc

    do n = 1,NVAR
       call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
            adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
            domlo,domhi,delta,xlo,bc(1,1,n))
    enddo

    ! Do this for all the variables, but we will overwrite the momenta below
    do n = 1,NVAR
       call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
            adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
            domlo,domhi,delta,xlo,bc(1,1,n))
    enddo

    if ( (bc(1,1,1).eq. EXT_DIR .or. bc(1,2,1).eq. EXT_DIR) .or.  &
         (bc(2,1,1).eq. EXT_DIR .or. bc(2,2,1).eq. EXT_DIR) .or. &
         (bc(3,1,1).eq. EXT_DIR .or. bc(3,2,1).eq. EXT_DIR) ) then
       call bl_error("NOT SET UP FOR EXT_DIR BCs IN HYPFILL")
    end if

    ! XLO
    if ( bc(1,1,1).eq. FOEXTRAP .and. adv_l1.lt.domlo(1)) then
       do k = adv_l3, adv_h3
          do j = adv_l2, adv_h2

             y = (dble(j) + 0.5d0) * delta(2) - center(2)
             z = (dble(k) + 0.5d0) * delta(3) - center(3)

             ic = domlo(1)
             xc = (dble(ic) + 0.5d0) * delta(1) - center(1)
             rc = sqrt(xc**2 + y**2 + z**2)

             momc = sqrt(adv(ic,j,k,UMX)**2 + adv(ic,j,k,UMY)**2 + adv(ic,j,k,UMZ)**2)

             do i = adv_l1, domlo(1)-1
                x = (dble(i) + 0.5d0) * delta(1) - center(1)
                r = sqrt(x**2 + y**2 + z**2)

                mom  = momc * (rc/r)**2

                ! Project along the normal
                adv(i,j,k,UMX) =  sign(1.d0,adv(ic,j,k,UMX)) * mom * (x/r)
                adv(i,j,k,UMY) =  sign(1.d0,adv(ic,j,k,UMY)) * mom * (y/r)
                adv(i,j,k,UMZ) =  sign(1.d0,adv(ic,j,k,UMZ)) * mom * (z/r)

             end do
          end do
       end do
    end if

    ! XHI
    if ( bc(1,2,1).eq. FOEXTRAP .and. adv_h1.gt.domhi(1)) then
       do k = adv_l3, adv_h3
          do j = adv_l2, adv_h2

             y = (dble(j) + 0.5d0) * delta(2) - center(2)
             z = (dble(k) + 0.5d0) * delta(3) - center(3)

             ic = domhi(1)
             xc = (dble(ic) + 0.5d0) * delta(1) - center(1)
             rc = sqrt(xc**2 + y**2 + z**2)

             momc = sqrt(adv(ic,j,k,UMX)**2 + adv(ic,j,k,UMY)**2 + adv(ic,j,k,UMZ)**2)

             do i = domhi(1)+1, adv_h1
                x = (dble(i) + 0.5d0) * delta(1) - center(1)
                r = sqrt(x**2 + y**2 + z**2)

                mom  = momc * (rc/r)**2

                ! Project along the normal
                adv(i,j,k,UMX) =  sign(1.d0,adv(ic,j,k,UMX)) * mom * (x/r)
                adv(i,j,k,UMY) =  sign(1.d0,adv(ic,j,k,UMY)) * mom * (y/r)
                adv(i,j,k,UMZ) =  sign(1.d0,adv(ic,j,k,UMZ)) * mom * (z/r)

             end do
          end do
       end do
    end if

    ! YLO
    if ( bc(2,1,1).eq. FOEXTRAP .and. adv_l2.lt.domlo(2)) then
       do k = adv_l3, adv_h3
          do i = adv_l1, adv_h1

             x = (dble(i) + 0.5d0) * delta(1) - center(1)
             z = (dble(k) + 0.5d0) * delta(3) - center(3)

             jc = domlo(2)
             yc = (dble(jc) + 0.5d0) * delta(2) - center(2)
             rc = sqrt(x**2 + yc**2 + z**2)

             momc = sqrt(adv(i,jc,k,UMX)**2 + adv(i,jc,k,UMY)**2 + adv(i,jc,k,UMZ)**2)

             do j = adv_l2, domlo(2)-1

                y = (dble(j) + 0.5d0) * delta(2) - center(2)
                r = sqrt(x**2 + y**2 + z**2)

                mom  = momc * (rc/r)**2

                ! Project along the normal
                adv(i,j,k,UMX) =  sign(1.d0,adv(i,jc,k,UMX)) * mom * (x/r)
                adv(i,j,k,UMY) =  sign(1.d0,adv(i,jc,k,UMY)) * mom * (y/r)
                adv(i,j,k,UMZ) =  sign(1.d0,adv(i,jc,k,UMZ)) * mom * (z/r)

             end do
          end do
       end do
    end if

    ! YHI
    if ( bc(2,2,1).eq. FOEXTRAP .and. adv_h2.gt.domhi(2)) then
       do k = adv_l3, adv_h3
          do i = adv_l1, adv_h1

             x = (dble(i) + 0.5d0) * delta(1) - center(1)
             z = (dble(k) + 0.5d0) * delta(3) - center(3)

             jc = domhi(2)
             yc = (dble(jc) + 0.5d0) * delta(2) - center(2)
             rc = sqrt(x**2 + yc**2 + z**2)

             momc = sqrt(adv(i,jc,k,UMX)**2 + adv(i,jc,k,UMY)**2 + adv(i,jc,k,UMZ)**2)

             do j = domhi(2)+1, adv_h2
                y = (dble(j) + 0.5d0) * delta(2) - center(2)
                r = sqrt(x**2 + y**2 + z**2)

                mom  = momc * (rc/r)**2

                ! Project along the normal
                adv(i,j,k,UMX) =  sign(1.d0,adv(i,jc,k,UMX)) * mom * (x/r)
                adv(i,j,k,UMY) =  sign(1.d0,adv(i,jc,k,UMY)) * mom * (y/r)
                adv(i,j,k,UMZ) =  sign(1.d0,adv(i,jc,k,UMZ)) * mom * (z/r)

             end do
          end do
       end do
    end if

    ! ZLO
    if ( bc(3,1,1).eq. FOEXTRAP .and. adv_l3.lt.domlo(3)) then
       do j = adv_l2, adv_h2
          do i = adv_l1, adv_h1

             x = (dble(i) + 0.5d0) * delta(1) - center(1)
             y = (dble(j) + 0.5d0) * delta(2) - center(2)

             kc = domlo(3)
             zc = (dble(kc) + 0.5d0) * delta(3) - center(3)
             rc = sqrt(x**2 + y**2 + zc**2)

             momc = sqrt(adv(i,j,kc,UMX)**2 + adv(i,j,kc,UMY)**2 + adv(i,j,kc,UMZ)**2)

             do k = adv_l3, domlo(3)-1
                z = (dble(k) + 0.5d0) * delta(3) - center(3)
                r = sqrt(x**2 + y**2 + z**2)

                mom  = momc * (rc/r)**2

                ! Project along the normal
                adv(i,j,k,UMX) =  sign(1.d0,adv(i,j,kc,UMX)) * mom * (x/r)
                adv(i,j,k,UMY) =  sign(1.d0,adv(i,j,kc,UMY)) * mom * (y/r)
                adv(i,j,k,UMZ) =  sign(1.d0,adv(i,j,kc,UMZ)) * mom * (z/r)

             end do
          end do
       end do
    end if

    ! ZHI
    if ( bc(3,2,1).eq. FOEXTRAP .and. adv_h3.gt.domhi(3)) then
       do j = adv_l2, adv_h2
          do i = adv_l1, adv_h1

             x = (dble(i) + 0.5d0) * delta(1) - center(1)
             y = (dble(j) + 0.5d0) * delta(2) - center(2)

             kc = domhi(3)
             zc = (dble(kc) + 0.5d0) * delta(3) - center(3)
             rc = sqrt(x**2 + y**2 + zc**2)

             momc = sqrt(adv(i,j,kc,UMX)**2 + adv(i,j,kc,UMY)**2 + adv(i,j,kc,UMZ)**2)

             do k = domhi(3)+1, adv_h3
                z = (dble(k) + 0.5d0) * delta(3) - center(3)
                r = sqrt(x**2 + y**2 + z**2)

                mom  = momc * (rc/r)**2

                ! Project along the normal
                adv(i,j,k,UMX) =  sign(1.d0,adv(i,j,kc,UMX)) * mom * (x/r)
                adv(i,j,k,UMY) =  sign(1.d0,adv(i,j,kc,UMY)) * mom * (y/r)
                adv(i,j,k,UMZ) =  sign(1.d0,adv(i,j,kc,UMZ)) * mom * (z/r)

             end do
          end do
       end do
    end if

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,delta,xlo,time,bc) bind(C,name="ca_denfill")
    
    implicit none
    
    include 'AMReX_bc_types.fi'
    
    integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
    integer bc(3,2,*)
    integer domlo(3), domhi(3)
    double precision delta(3), xlo(3), time
    double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)
    logical rho_only
    integer i,j,k

    call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,domlo,domhi,delta,xlo,bc)

  end subroutine ca_denfill



#ifdef GRAVITY
  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) bind(C,name="ca_gravxfill")

    use probdata_module
    
    implicit none
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    double precision delta(3), xlo(3), time
    double precision grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

    call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravxfill



  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) bind(C,name="ca_gravyfill")

    use probdata_module
    
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    double precision delta(3), xlo(3), time
    double precision grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

    call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravyfill



  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                          domlo,domhi,delta,xlo,time,bc) bind(C,name="ca_gravzfill")

    use probdata_module
    
    implicit none
    
    include 'AMReX_bc_types.fi'

    integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
    integer :: bc(3,2,*)
    integer :: domlo(3), domhi(3)
    double precision delta(3), xlo(3), time
    double precision grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

    call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,domlo,domhi,delta,xlo,bc)

  end subroutine ca_gravzfill



  subroutine ca_phigravfill(phi,phi_l1,phi_l2,phi_l3, &
                            phi_h1,phi_h2,phi_h3,domlo,domhi,delta,xlo,time,bc) bind(C,name="ca_phigravfill")

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
    integer          :: bc(3,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

    call filcc(phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
         domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill
#endif
  
end module bc_fill_module
