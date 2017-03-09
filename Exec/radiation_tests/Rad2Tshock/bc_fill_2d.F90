module bc_fill_module

  use amrex_fort_module, only : rt => c_real
  implicit none

  public

contains

  ! ::: -----------------------------------------------------------
  
  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")
 
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UFX, UTEMP
    use network, only : nspec, naux
    use eos_module
    use probdata_module

    use amrex_fort_module, only : rt => c_real
    implicit none

    include 'AMReX_bc_types.fi'
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer i, j, n

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
       call filcc(adv(adv_l1,adv_l2,n),adv_l1,adv_l2,adv_h1,adv_h2, &
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
       do j = adv_l2, adv_h2
          do i = adv_l1, domlo(1)-1
             adv(i,j,URHO) = rho0
             adv(i,j,UMX) = rho0*v0
             adv(i,j,UMY) = 0.e0_rt
             adv(i,j,UMZ) = 0.e0_rt
             adv(i,j,UFS) = adv(i,j,URHO)
             if (naux > 0) then
                adv(i,j,UFX) = adv(i,j,URHO)           
             end if
             adv(i,j,UEINT) = eint0
             adv(i,j,UEDEN) = etot0
             adv(i,j,UTEMP) = T0
          enddo
       enddo
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
       do j = adv_l2, adv_h2
          do i = domhi(1)+1, adv_h1
             adv(i,j,URHO) = rho1
             adv(i,j,UMX) = rho1*v1
             adv(i,j,UMY) = 0.e0_rt
             adv(i,j,UMZ) = 0.e0_rt
             adv(i,j,UFS) = adv(i,j,URHO)
             if (naux > 0) then
                adv(i,j,UFX) = adv(i,j,URHO)           
             end if
             adv(i,j,UEINT) = eint1
             adv(i,j,UEDEN) = etot1
             adv(i,j,UTEMP) = T1
          end do
       enddo
    end if
       
  end subroutine ca_hypfill

  ! ::: 
  ! ::: -----------------------------------------------------------
  ! :::

  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")

    use probdata_module

    use amrex_fort_module, only : rt => c_real
    implicit none
    include 'AMReX_bc_types.fi'
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2)
    integer i,j

    !     Note: this function should not be needed, technically, but is provided
    !     to filpatch because there are many times in the algorithm when just
    !     the density is needed.  We try to rig up the filling so that the same
    !     function is called here and in hypfill where all the states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       do j = adv_l2, adv_h2
          do i = adv_l1, domlo(1)-1
             adv(i,j) = rho0
          enddo
       enddo
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
       do j = adv_l2, adv_h2
          do i = domhi(1)+1, adv_h1
             adv(i,j) = rho1
          enddo
       enddo
    end if

  end subroutine ca_denfill

  ! ::: 
  ! ::: -----------------------------------------------------------
  ! :::

  subroutine ca_phigravfill(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                            domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_phigravfill")

    use amrex_fort_module, only : rt => c_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_h1,phi_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    real(rt)         :: delta(2), xlo(2), time
    real(rt)         :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
               domlo,domhi,delta,xlo,bc)

  end subroutine ca_phigravfill

  ! ::: 
  ! ::: -----------------------------------------------------------
  ! :::

  subroutine ca_radfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_radfill")

    use probdata_module

    use amrex_fort_module, only : rt => c_real
    implicit none
    include 'AMReX_bc_types.fi'
    integer adv_l1,adv_l2,adv_h1,adv_h2
    integer bc(2,2,*)
    integer domlo(2), domhi(2)
    real(rt)         delta(2), xlo(2), time
    real(rt)         adv(adv_l1:adv_h1,adv_l2:adv_h2)
    integer i,j

    !     Note: this function should not be needed, technically, but is provided
    !     to filpatch because there are many times in the algorithm when just
    !     the density is needed.  We try to rig up the filling so that the same
    !     function is called here and in hypfill where all the states are filled.

    call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

    !     XLO
    if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
       do j = adv_l2, adv_h2
          do i = adv_l1, domlo(1)-1
             adv(i,j) = adv(domlo(1),j)
          enddo
       enddo
    end if

    !     XHI
    if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
       do j = adv_l2, adv_h2
          do i = domhi(1)+1, adv_h1
             adv(i,j) = adv(domhi(1),j)
          enddo
       enddo
    end if

  end subroutine ca_radfill

end module bc_fill_module
