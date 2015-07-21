subroutine ca_hypfill(adv,adv_l1,adv_h1, &
                      domlo,domhi,delta,xlo,time,bc)

  use meth_params_module, only: NVAR  
  
  implicit none

  include 'bc_types.fi'

  integer          :: adv_l1,adv_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1,NVAR)

  integer          :: n

  do n = 1,NVAR
     call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

end subroutine ca_hypfill



subroutine ca_denfill(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: adv_l1,adv_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1)
  
  call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

end subroutine ca_denfill



subroutine ca_gravxfill(grav,grav_l1,grav_h1, &
                        domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: grav_l1,grav_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: grav(grav_l1:grav_h1)

  call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)

end subroutine ca_gravxfill



subroutine ca_reactfill(react,react_l1,react_h1, &
                        domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: react_l1,react_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: react(react_l1:react_h1)

  call filcc(react,react_l1,react_h1, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_reactfill
