module bc_fill_module

  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only: dim

  implicit none

  include 'AMReX_bc_types.fi'


  public

contains

  subroutine ca_hypfill(adv, adv_lo, adv_hi, &
                        domlo, domhi, delta, xlo, time, bc) bind(C,name="ca_hypfill")

    use amrex_filcc_module, only: amrex_filccn
    use amrex_error_module
    use meth_params_module, only : NVAR,UMX,UMY,UMZ
    use prob_params_module, only : center

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(dim,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

    integer          :: i, j, k, n
    integer          :: ic,jc,kc
    real(rt)         :: x,y,z,r
    real(rt)         :: xc,yc,zc,rc
    real(rt)         :: mom,momc

    ! Do this for all the variables, but we will overwrite the momenta below
    call amrex_filccn(adv_lo, adv_hi, adv, adv_lo, adv_hi, NVAR, domlo, domhi, delta, xlo, bc)

#if AMREX_SPACEDIM == 3
    if ( (bc(1,1,1) == EXT_DIR .or. bc(1,2,1) == EXT_DIR) .or.  &
         (bc(2,1,1) == EXT_DIR .or. bc(2,2,1) == EXT_DIR) .or. &
         (bc(3,1,1) == EXT_DIR .or. bc(3,2,1) == EXT_DIR) ) then
       call amrex_error("NOT SET UP FOR EXT_DIR BCs IN HYPFILL")
    end if

    ! XLO
    if ( bc(1,1,1) == FOEXTRAP .and. adv_lo(1) < domlo(1)) then
       do k = adv_lo(3), adv_hi(3)
          do j = adv_lo(2), adv_hi(2)

             y = (dble(j) + HALF) * delta(2) - center(2)
             z = (dble(k) + HALF) * delta(3) - center(3)

             ic = domlo(1)
             xc = (dble(ic) + HALF) * delta(1) - center(1)
             rc = sqrt(xc**2 + y**2 + z**2)

             momc = sqrt(adv(ic,j,k,UMX)**2 + adv(ic,j,k,UMY)**2 + adv(ic,j,k,UMZ)**2)

             do i = adv_lo(1), domlo(1)-1
                x = (dble(i) + HALF) * delta(1) - center(1)
                r = sqrt(x**2 + y**2 + z**2)

                mom  = momc * (rc/r)**2

                ! Project along the normal
                adv(i,j,k,UMX) =  sign(1.e0_rt,adv(ic,j,k,UMX)) * mom * (x/r)
                adv(i,j,k,UMY) =  sign(1.e0_rt,adv(ic,j,k,UMY)) * mom * (y/r)
                adv(i,j,k,UMZ) =  sign(1.e0_rt,adv(ic,j,k,UMZ)) * mom * (z/r)

             end do
          end do
       end do
    end if

    ! XHI
    if ( bc(1,2,1) == FOEXTRAP .and. adv_hi(1) > domhi(1)) then
       do k = adv_lo(3), adv_hi(3)
          do j = adv_lo(2), adv_hi(2)

             y = (dble(j) + HALF) * delta(2) - center(2)
             z = (dble(k) + HALF) * delta(3) - center(3)

             ic = domhi(1)
             xc = (dble(ic) + HALF) * delta(1) - center(1)
             rc = sqrt(xc**2 + y**2 + z**2)

             momc = sqrt(adv(ic,j,k,UMX)**2 + adv(ic,j,k,UMY)**2 + adv(ic,j,k,UMZ)**2)

             do i = domhi(1)+1, adv_hi(1)
                x = (dble(i) + HALF) * delta(1) - center(1)
                r = sqrt(x**2 + y**2 + z**2)

                mom  = momc * (rc/r)**2

                ! Project along the normal
                adv(i,j,k,UMX) =  sign(1.e0_rt,adv(ic,j,k,UMX)) * mom * (x/r)
                adv(i,j,k,UMY) =  sign(1.e0_rt,adv(ic,j,k,UMY)) * mom * (y/r)
                adv(i,j,k,UMZ) =  sign(1.e0_rt,adv(ic,j,k,UMZ)) * mom * (z/r)

             end do
          end do
       end do
    end if

    ! YLO
    if ( bc(2,1,1) == FOEXTRAP .and. adv_lo(2) < domlo(2)) then
       do k = adv_lo(3), adv_hi(3)
          do i = adv_lo(1), adv_hi(1)

             x = (dble(i) + HALF) * delta(1) - center(1)
             z = (dble(k) + HALF) * delta(3) - center(3)

             jc = domlo(2)
             yc = (dble(jc) + HALF) * delta(2) - center(2)
             rc = sqrt(x**2 + yc**2 + z**2)

             momc = sqrt(adv(i,jc,k,UMX)**2 + adv(i,jc,k,UMY)**2 + adv(i,jc,k,UMZ)**2)

             do j = adv_lo(2), domlo(2)-1

                y = (dble(j) + HALF) * delta(2) - center(2)
                r = sqrt(x**2 + y**2 + z**2)

                mom  = momc * (rc/r)**2

                ! Project along the normal
                adv(i,j,k,UMX) =  sign(1.e0_rt,adv(i,jc,k,UMX)) * mom * (x/r)
                adv(i,j,k,UMY) =  sign(1.e0_rt,adv(i,jc,k,UMY)) * mom * (y/r)
                adv(i,j,k,UMZ) =  sign(1.e0_rt,adv(i,jc,k,UMZ)) * mom * (z/r)

             end do
          end do
       end do
    end if

    ! YHI
    if ( bc(2,2,1) == FOEXTRAP .and. adv_hi(2) > domhi(2)) then
       do k = adv_lo(3), adv_hi(3)
          do i = adv_lo(1), adv_hi(1)

             x = (dble(i) + HALF) * delta(1) - center(1)
             z = (dble(k) + HALF) * delta(3) - center(3)

             jc = domhi(2)
             yc = (dble(jc) + HALF) * delta(2) - center(2)
             rc = sqrt(x**2 + yc**2 + z**2)

             momc = sqrt(adv(i,jc,k,UMX)**2 + adv(i,jc,k,UMY)**2 + adv(i,jc,k,UMZ)**2)

             do j = domhi(2)+1, adv_hi(2)
                y = (dble(j) + HALF) * delta(2) - center(2)
                r = sqrt(x**2 + y**2 + z**2)

                mom  = momc * (rc/r)**2

                ! Project along the normal
                adv(i,j,k,UMX) =  sign(1.e0_rt,adv(i,jc,k,UMX)) * mom * (x/r)
                adv(i,j,k,UMY) =  sign(1.e0_rt,adv(i,jc,k,UMY)) * mom * (y/r)
                adv(i,j,k,UMZ) =  sign(1.e0_rt,adv(i,jc,k,UMZ)) * mom * (z/r)

             end do
          end do
       end do
    end if

    ! ZLO
    if ( bc(3,1,1) == FOEXTRAP .and. adv_lo(3) < domlo(3)) then
       do j = adv_lo(2), adv_hi(2)
          do i = adv_lo(1), adv_hi(1)

             x = (dble(i) + HALF) * delta(1) - center(1)
             y = (dble(j) + HALF) * delta(2) - center(2)

             kc = domlo(3)
             zc = (dble(kc) + HALF) * delta(3) - center(3)
             rc = sqrt(x**2 + y**2 + zc**2)

             momc = sqrt(adv(i,j,kc,UMX)**2 + adv(i,j,kc,UMY)**2 + adv(i,j,kc,UMZ)**2)

             do k = adv_lo(3), domlo(3)-1
                z = (dble(k) + HALF) * delta(3) - center(3)
                r = sqrt(x**2 + y**2 + z**2)

                mom  = momc * (rc/r)**2

                ! Project along the normal
                adv(i,j,k,UMX) =  sign(1.e0_rt,adv(i,j,kc,UMX)) * mom * (x/r)
                adv(i,j,k,UMY) =  sign(1.e0_rt,adv(i,j,kc,UMY)) * mom * (y/r)
                adv(i,j,k,UMZ) =  sign(1.e0_rt,adv(i,j,kc,UMZ)) * mom * (z/r)

             end do
          end do
       end do
    end if

    ! ZHI
    if ( bc(3,2,1) == FOEXTRAP .and. adv_hi(3) > domhi(3)) then
       do j = adv_lo(2), adv_hi(2)
          do i = adv_lo(1), adv_hi(1)

             x = (dble(i) + HALF) * delta(1) - center(1)
             y = (dble(j) + HALF) * delta(2) - center(2)

             kc = domhi(3)
             zc = (dble(kc) + HALF) * delta(3) - center(3)
             rc = sqrt(x**2 + y**2 + zc**2)

             momc = sqrt(adv(i,j,kc,UMX)**2 + adv(i,j,kc,UMY)**2 + adv(i,j,kc,UMZ)**2)

             do k = domhi(3)+1, adv_hi(3)
                z = (dble(k) + HALF) * delta(3) - center(3)
                r = sqrt(x**2 + y**2 + z**2)

                mom  = momc * (rc/r)**2

                ! Project along the normal
                adv(i,j,k,UMX) =  sign(1.e0_rt,adv(i,j,kc,UMX)) * mom * (x/r)
                adv(i,j,k,UMY) =  sign(1.e0_rt,adv(i,j,kc,UMY)) * mom * (y/r)
                adv(i,j,k,UMZ) =  sign(1.e0_rt,adv(i,j,kc,UMZ)) * mom * (z/r)

             end do
          end do
       end do
    end if
#endif

  end subroutine ca_hypfill



  subroutine ca_denfill(adv, adv_lo, adv_hi, &
                        domlo, domhi, delta, xlo, time, bc) bind(C,name="ca_denfill")

    use amrex_filcc_module, only: amrex_filccn
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

    call amrex_filccn(adv_lo, adv_hi, adv, adv_lo, adv_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_denfill



#ifdef GRAVITY
  subroutine ca_gravxfill(grav, grav_lo, grav_hi, &
                          domlo, domhi, delta, xlo, time, bc) bind(C,name="ca_gravxfill")

    use amrex_filcc_module, only: amrex_filccn
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    call amrex_filccn(grav_lo, grav_hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_gravxfill


  subroutine ca_gravyfill(grav, grav_lo, grav_hi, &
                          domlo, domhi, delta, xlo, time, bc) bind(C,name="ca_gravyfill")

    use amrex_filcc_module, only: amrex_filccn
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    call amrex_filccn(grav_lo, grav_hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_gravyfill


  subroutine ca_gravzfill(grav, grav_lo, grav_hi, &
                          domlo, domhi, delta, xlo, time, bc) bind(C,name="ca_gravzfill")

    use amrex_filcc_module, only: amrex_filccn
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))

    call amrex_filccn(grav_lo, grav_hi, grav, grav_lo, grav_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_gravzfill


  subroutine ca_phigravfill(phi, phi_lo, phi_hi, &
                            domlo, domhi, delta, xlo, time, bc) bind(C,name="ca_phigravfill")

    use amrex_filcc_module, only: amrex_filccn
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3), time
    real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

    call amrex_filccn(phi_lo, phi_hi, phi, phi_lo, phi_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine ca_phigravfill
#endif

end module bc_fill_module
