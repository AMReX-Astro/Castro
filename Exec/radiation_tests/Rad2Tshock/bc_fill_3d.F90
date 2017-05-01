module bc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  ! ::: -----------------------------------------------------------
  
  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")
 
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UFX, UTEMP
    use network, only : nspec, naux
    use eos_module
    use probdata_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'
    integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
    integer bc(3,2,*)
    integer domlo(3), domhi(3)
    real(rt)         delta(3), xlo(3), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

    integer i, j, k, n

    type(eos_t) :: eos_state
    real(rt)        , save :: eint0, etot0, eint1, etot1
    logical, save :: first_call = .true.

    if (first_call) then
       first_call = .false.

       eos_state % rho = rho0
       eos_state % T   =   T0
       eos_state % xn  = 1.e0_rt

       call eos(eos_input_rt, eos_state)

       eint0 = rho0 * eos_state % e
       etot0 = eint0 + 0.5*rho0*v0**2

       eos_state % rho = rho1
       eos_state % T   =   T1

       call eos(eos_input_rt, eos_state)

       eint1 = rho1 * eos_state % e
       etot1 = eint1 + 0.5*rho1*v1**2         
    end if

    do n = 1,NVAR
       call filcc(adv(adv_l1,adv_l2,adv_l3,n),adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                  domlo,domhi,delta,xlo,bc(1,1,n))
    enddo

    !     The strategy here is to set Dirichlet condition for inflow and
    !     outflow boundaries, and let the Riemann solver sort out the
    !     proper upwinding.  However, this decision makes this routine
    !     look somewhat non-orthodox, in that we need to set external
    !     values in either case....how do we know it's Outflow?  We have
    !     to assume that the setup routines converted Outflow to FOEXTRAP.

    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       do k = adv_l3, adv_h3
          do j = adv_l2, adv_h2
             do i = adv_l1, domlo(1)-1
                adv(i,j,k,URHO) = rho0
                adv(i,j,k,UMX) = rho0*v0
                adv(i,j,k,UMY) = 0.e0_rt
                adv(i,j,k,UMZ) = 0.e0_rt
                adv(i,j,k,UFS) = adv(i,j,k,URHO)
                if (naux > 0) then
                   adv(i,j,k,UFX) = adv(i,j,k,URHO)           
                end if
                adv(i,j,k,UEINT) = eint0
                adv(i,j,k,UEDEN) = etot0
                adv(i,j,k,UTEMP) = T0
             enddo
          enddo
       enddo
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
       do k = adv_l3, adv_h3
          do j = adv_l2, adv_h2
             do i = domhi(1)+1, adv_h1
                adv(i,j,k,URHO) = rho1
                adv(i,j,k,UMX) = rho1*v1
                adv(i,j,k,UMY) = 0.e0_rt
                adv(i,j,k,UMZ) = 0.e0_rt
                adv(i,j,k,UFS) = adv(i,j,k,URHO)
                if (naux > 0) then
                   adv(i,j,k,UFX) = adv(i,j,k,URHO)           
                end if
                adv(i,j,k,UEINT) = eint1
                adv(i,j,k,UEDEN) = etot1
                adv(i,j,k,UTEMP) = T1
             end do
          enddo
       enddo
    end if
       
  end subroutine ca_hypfill

  ! ::: 
  ! ::: -----------------------------------------------------------
  ! :::

  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")

    use probdata_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    include 'AMReX_bc_types.fi'
    integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
    integer bc(3,2,*)
    integer domlo(3), domhi(3)
    real(rt)         delta(3), xlo(3), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)
    integer i,j,k

    !     Note: this function should not be needed, technically, but is provided
    !     to filpatch because there are many times in the algorithm when just
    !     the density is needed.  We try to rig up the filling so that the same
    !     function is called here and in hypfill where all the states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,domlo,domhi,delta,xlo,bc)

    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       do k = adv_l3, adv_h3
          do j = adv_l2, adv_h2
             do i = adv_l1, domlo(1)-1
                adv(i,j,k) = rho0
             enddo
          enddo
       enddo
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
       do k = adv_l3, adv_h3
          do j = adv_l2, adv_h2
             do i = domhi(1)+1, adv_h1
                adv(i,j,k) = rho1
             enddo
          enddo
       enddo
    end if

  end subroutine ca_denfill

  ! ::: 
  ! ::: -----------------------------------------------------------
  ! :::

  subroutine ca_phigravfill(phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
                            domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_phigravfill")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
    integer          :: bc(3,2,*)
    integer          :: domlo(3), domhi(3)
    real(rt)         :: delta(3), xlo(3), time
    real(rt)         :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

    call filcc(phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill

  ! ::: 
  ! ::: -----------------------------------------------------------
  ! :::

  subroutine ca_radfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_radfill")

    use probdata_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    include 'AMReX_bc_types.fi'
    integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
    integer bc(3,2,*)
    integer domlo(3), domhi(3)
    real(rt)         delta(3), xlo(3), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)
    integer i,j,k

    !     Note: this function should not be needed, technically, but is provided
    !     to filpatch because there are many times in the algorithm when just
    !     the density is needed.  We try to rig up the filling so that the same
    !     function is called here and in hypfill where all the states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,domlo,domhi,delta,xlo,bc)

    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       do k = adv_l3, adv_h3
          do j = adv_l2, adv_h2
             do i = adv_l1, domlo(1)-1
                adv(i,j,k) = adv(domlo(1),j,k)
             enddo
          enddo
       enddo
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
       do k = adv_l3, adv_h3
          do j = adv_l2, adv_h2
             do i = domhi(1)+1, adv_h1
                adv(i,j,k) = adv(domhi(1),j,k)
             enddo
          enddo
       enddo
    end if

  end subroutine ca_radfill

end module bc_fill_module
