! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                      domlo,domhi,delta,xlo,time,bc)

  use meth_params_module, only: NVAR
  use hse_bc_module

  implicit none
  include 'bc_types.fi'

  integer          :: adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

  integer :: n

  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
                adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

  ! ZLO
  if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     call hse_bc_zlo(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                     domlo,domhi,delta,xlo,time,bc)
  end if

  ! ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     call hse_bc_zhi(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                     domlo,domhi,delta,xlo,time,bc)
  end if

end subroutine ca_hypfill

! ::: -----------------------------------------------------------

subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                      domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  use hse_bc_module

  implicit none
  include 'bc_types.fi'

  integer          :: adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  ! Note: this function should not be needed, technically, but is
  ! provided to filpatch because there are many times in the algorithm
  ! when just the density is needed.  We try to rig up the filling so
  ! that the same function is called here and in hypfill where all the
  ! states are filled.

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
             domlo,domhi,delta,xlo,bc)

  !     ZLO
  if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     call hse_bc_zlo(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                     domlo,domhi,delta,xlo,time,bc, density_only=.true.)
  end if

  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_h3.gt.domhi(3)) then
     call hse_bc_zhi(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                     domlo,domhi,delta,xlo,time,bc, density_only=.true.)
  end if

end subroutine ca_denfill

! ::: -----------------------------------------------------------


subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                       domlo,domhi,delta,xlo,time,bc)

  implicit none
  include 'bc_types.fi'

  integer          :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

  integer :: bc_temp(3,2)

  bc_temp(:,:) = bc(:,:,1)

  if ( bc(3,1,1).eq.EXT_DIR .and. grav_l3.lt.domlo(3)) then
     bc_temp(3,1) = FOEXTRAP
  endif

  call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
             domlo,domhi,delta,xlo,bc_temp)

end subroutine ca_gravxfill


subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                       domlo,domhi,delta,xlo,time,bc)

  implicit none
  include 'bc_types.fi'

  integer          :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

  integer :: bc_temp(3,2)

  bc_temp(:,:) = bc(:,:,1)

  if ( bc(3,1,1).eq.EXT_DIR .and. grav_l3.lt.domlo(3)) then
     bc_temp(3,1) = FOEXTRAP
  endif

  call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
             domlo,domhi,delta,xlo,bc_temp)

end subroutine ca_gravyfill


subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                       domlo,domhi,delta,xlo,time,bc)

  implicit none
  include 'bc_types.fi'

  integer          :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

  integer :: bc_temp(3,2)

  bc_temp(:,:) = bc(:,:,1)

  if ( bc(3,1,1).eq.EXT_DIR .and. grav_l3.lt.domlo(3)) then
     bc_temp(3,1) = FOEXTRAP
  endif

  call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
             domlo,domhi,delta,xlo,bc_temp)

end subroutine ca_gravzfill


subroutine ca_reactfill(react,react_l1,react_l2,react_l3, &
                        react_h1,react_h2,react_h3,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: react_l1,react_l2,react_l3,react_h1,react_h2,react_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: react(react_l1:react_h1,react_l2:react_h2,react_l3:react_h3)

  integer :: bc_temp(3,2)

  bc_temp(:,:) = bc(:,:,1)

  if ( bc(3,1,1).eq.EXT_DIR .and. react_l3.lt.domlo(3)) then
     bc_temp(3,1) = FOEXTRAP
  endif

  call filcc(react,react_l1,react_l2,react_l3,react_h1,react_h2,react_h3, &
             domlo,domhi,delta,xlo,bc_temp)

end subroutine ca_reactfill



subroutine ca_phigravfill(phi,phi_l1,phi_l2,phi_l3, &
                          phi_h1,phi_h2,phi_h3,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

  integer :: bc_temp(3,2)

  bc_temp(:,:) = bc(:,:,1)

  if ( bc(3,1,1).eq.EXT_DIR .and. phi_l3.lt.domlo(3)) then
     bc_temp(3,1) = FOEXTRAP
  endif

  call filcc(phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
             domlo,domhi,delta,xlo,bc_temp)

end subroutine ca_phigravfill



subroutine ca_radfill(rad,rad_l1,rad_l2,rad_l3, &
                      rad_h1,rad_h2,rad_h3,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: rad_l1,rad_l2,rad_l3,rad_h1,rad_h2,rad_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: rad(rad_l1:rad_h1,rad_l2:rad_h2,rad_l3:rad_h3)

  integer :: bc_temp(3,2)

  bc_temp(:,:) = bc(:,:,1)

  if ( bc(3,1,1).eq.EXT_DIR .and. rad_l3.lt.domlo(3)) then
     bc_temp(3,1) = FOEXTRAP
  endif

  call filcc(rad,rad_l1,rad_l2,rad_l3,rad_h1,rad_h2,rad_h3, &
             domlo,domhi,delta,xlo,bc_temp)

end subroutine ca_radfill
