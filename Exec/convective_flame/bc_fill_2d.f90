! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,time,bc)

  use meth_params_module, only : NVAR
  use hse_bc_module

  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

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

! ::: -----------------------------------------------------------

subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
     domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  use interpolate_module
  use eos_module, only: gamma_const
  use hse_bc_module

  implicit none
  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

  integer i,j
  double precision y,H

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

! ::: -----------------------------------------------------------


subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
     domlo,domhi,delta,xlo,time,bc)

  ! note that this is used for both gravity and rotation.  For
  ! our lower boundary, where we have inflow, we'll implement
  ! our own zero-gradient (outflow) BC.

  use probdata_module
  implicit none
  include 'bc_types.fi'

  integer :: grav_l1,grav_l2,grav_h1,grav_h2
  integer :: bc(2,2,*)
  integer :: domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision grav(grav_l1:grav_h1,grav_l2:grav_h2)

  call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

end subroutine ca_gravxfill

! ::: -----------------------------------------------------------

subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
     domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  use meth_params_module, only: const_grav

  implicit none
  include 'bc_types.fi'

  integer :: grav_l1,grav_l2,grav_h1,grav_h2
  integer :: bc(2,2,*)
  integer :: domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision grav(grav_l1:grav_h1,grav_l2:grav_h2)

  integer :: i, j

  call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. grav_l2.lt.domlo(2)) then
     do j=domlo(2)-1,grav_l2,-1
        do i=grav_l1, grav_h1
           grav(i,j) = const_grav
        enddo
     enddo
     
  end if

  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. grav_h2.gt.domhi(2)) then
     do j=domhi(2)+1,grav_h2
        do i=grav_l1, grav_h1
           grav(i,j) = const_grav
        enddo
     enddo
     
  end if

end subroutine ca_gravyfill

! ::: -----------------------------------------------------------

subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
     domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  implicit none
  include 'bc_types.fi'

  integer :: grav_l1,grav_l2,grav_h1,grav_h2
  integer :: bc(2,2,*)
  integer :: domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision grav(grav_l1:grav_h1,grav_l2:grav_h2)

  call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

end subroutine ca_gravzfill

! ::: -----------------------------------------------------------

subroutine ca_phigravfill(phi,phi_l1,phi_l2, &
                          phi_h1,phi_h2,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: phi_l1,phi_l2,phi_h1,phi_h2
  integer          :: bc(2,2,*)
  integer          :: domlo(2), domhi(2)
  double precision :: delta(2), xlo(2), time
  double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

  call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_phigravfill


subroutine ca_reactfill(react,react_l1,react_l2, &
                        react_h1,react_h2,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: react_l1,react_l2,react_h1,react_h2
  integer          :: bc(2,2,*)
  integer          :: domlo(2), domhi(2)
  double precision :: delta(2), xlo(2), time
  double precision :: react(react_l1:react_h1,react_l2:react_h2)

  call filcc(react,react_l1,react_l2,react_h1,react_h2, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_reactfill
