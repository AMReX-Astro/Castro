! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_h1, &
     domlo,domhi,delta,xlo,time,bc)
  
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UFX, UTEMP
  use network, only : nspec, naux
  use eos_module
  use probdata_module
  
  implicit none
  
  include 'bc_types.fi'
  integer adv_l1,adv_h1
  integer bc(1,2,*)
  integer domlo(1), domhi(1)
  double precision delta(1), xlo(1), time
  double precision adv(adv_l1:adv_h1,NVAR)
  
  integer i, n
  
  type(eos_t) :: eos_state
  double precision, save :: eint0, etot0, eint1, etot1
  logical, save :: first_call = .true.
  
  if (first_call) then
     first_call = .false.
     
     eos_state % rho = rho0
     eos_state % T   =   T0
     eos_state % xn  = 1.d0

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
     call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
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
     do i = adv_l1, domlo(1)-1
        adv(i,URHO) = rho0
        adv(i,UMX) = rho0*v0
        adv(i,UMY) = 0.d0
        adv(i,UMZ) = 0.d0
        adv(i,UFS) = adv(i,URHO)
        if (naux>0) then
           adv(i,UFX) = adv(i,URHO)           
        end if
        adv(i,UEINT) = eint0
        adv(i,UEDEN) = etot0
        adv(i,UTEMP) = T0
     end do
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     do i = domhi(1)+1, adv_h1
        adv(i,URHO) = rho1
        adv(i,UMX) = rho1*v1
        adv(i,UMY) = 0.d0
        adv(i,UMZ) = 0.d0
        adv(i,UFS) = adv(i,URHO)
        if (naux>0) then
           adv(i,UFX) = adv(i,URHO)           
        end if
        adv(i,UEINT) = eint1
        adv(i,UEDEN) = etot1
        adv(i,UTEMP) = T1
     end do
  end if
  
end subroutine ca_hypfill

! ::: 
! ::: -----------------------------------------------------------
! :::

subroutine ca_denfill(adv,adv_l1,adv_h1, &
     domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_h1
  integer bc(1,2,*)
  integer domlo(1), domhi(1)
  double precision delta(1), xlo(1), time
  double precision adv(adv_l1:adv_h1)
  integer i
  
  !     Note: this function should not be needed, technically, but is provided
  !     to filpatch because there are many times in the algorithm when just
  !     the density is needed.  We try to rig up the filling so that the same
  !     function is called here and in hypfill where all the states are filled.
  
  call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)
  
  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     do i = adv_l1, domlo(1)-1
        adv(i) = rho0
     end do
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     do i = domhi(1)+1, adv_h1
        adv(i) = rho1
     end do
  end if
  
end subroutine ca_denfill

! ::: 
! ::: -----------------------------------------------------------
! :::

subroutine ca_phigravfill(phi,phi_l1,phi_h1, &
                          domlo,domhi,delta,xlo,time,bc)

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

! ::: 
! ::: -----------------------------------------------------------
! :::

subroutine ca_radfill(adv,adv_l1,adv_h1, &
     domlo,domhi,delta,xlo,time,bc)
  
  use probdata_module
  
  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_h1
  integer bc(1,2,*)
  integer domlo(1), domhi(1)
  double precision delta(1), xlo(1), time
  double precision adv(adv_l1:adv_h1)
  integer i
  
  !     Note: this function should not be needed, technically, but is provided
  !     to filpatch because there are many times in the algorithm when just
  !     the density is needed.  We try to rig up the filling so that the same
  !     function is called here and in hypfill where all the states are filled.
  
  call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)
  
  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     do i = adv_l1, domlo(1)-1
        adv(i) = adv(domlo(1))
     end do
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     do i = domhi(1)+1, adv_h1
        adv(i) = adv(domhi(1))
     end do
  end if
  
end subroutine ca_radfill

