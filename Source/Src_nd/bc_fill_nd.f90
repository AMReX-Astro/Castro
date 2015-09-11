subroutine ca_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc)

  use meth_params_module, only: NVAR
  
  implicit none

  include 'bc_types.fi'

  integer          :: adv_lo(3),adv_hi(3)
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

  integer          :: n

  do n = 1,NVAR
     call filcc(adv(adv_lo(1),adv_lo(2),adv_lo(3),n), &
                adv_lo(1),adv_lo(2),adv_lo(3),adv_hi(1),adv_hi(2),adv_hi(3), &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

end subroutine ca_hypfill



subroutine ca_denfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: adv_lo(3),adv_hi(3)
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

  call filcc(adv,adv_lo(1),adv_lo(2),adv_lo(3),adv_hi(1),adv_hi(2),adv_hi(3), &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_denfill



subroutine ca_gravxfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: grav_lo(3),grav_hi(3)
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))
  
  call filcc(grav,grav_lo(1),grav_lo(2),grav_lo(3),grav_hi(1),grav_hi(2),grav_hi(3), &
       domlo,domhi,delta,xlo,bc)

end subroutine ca_gravxfill



subroutine ca_gravyfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: grav_lo(3),grav_hi(3)
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))
  
  call filcc(grav,grav_lo(1),grav_lo(2),grav_lo(3),grav_hi(1),grav_hi(2),grav_hi(3), &
       domlo,domhi,delta,xlo,bc)

end subroutine ca_gravyfill



subroutine ca_gravzfill(grav,grav_lo,grav_hi,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: grav_lo(3),grav_hi(3)
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))
  
  call filcc(grav,grav_lo(1),grav_lo(2),grav_lo(3),grav_hi(1),grav_hi(2),grav_hi(3), &
       domlo,domhi,delta,xlo,bc)

end subroutine ca_gravzfill



subroutine ca_reactfill(react,react_lo,react_hi,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: react_lo(3),react_hi(3)
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: react(react_lo(1):react_hi(1),react_lo(2):react_hi(2),react_lo(3):react_hi(3))

  call filcc(react,react_lo(1),react_lo(2),react_lo(3),react_hi(1),react_hi(2),react_hi(3), &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_reactfill



subroutine ca_phigravfill(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: phi_lo(3),phi_hi(3)
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

  call filcc(phi,phi_lo(1),phi_lo(2),phi_lo(3),phi_hi(1),phi_hi(2),phi_hi(3), &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_phigravfill



subroutine ca_radfill(rad,rad_lo,rad_hi,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: rad_lo(3),rad_hi(3)
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: rad(rad_lo(1):rad_hi(1),rad_lo(2):rad_hi(2),rad_lo(3):rad_hi(3))

  call filcc(rad,rad_lo(1),rad_lo(2),rad_lo(3),rad_hi(1),rad_hi(2),rad_hi(3), &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_radfill
