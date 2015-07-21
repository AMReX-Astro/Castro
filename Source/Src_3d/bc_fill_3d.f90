subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                      adv_h3,domlo,domhi,delta,xlo,time,bc)

  use meth_params_module, only: NVAR
  
  implicit none

  include 'bc_types.fi'

  integer          :: adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

  integer          :: n

  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
                adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

end subroutine ca_hypfill



subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                      adv_h3,domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_denfill



subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                        domlo,domhi,delta,xlo,time,bc)

  implicit none

  include 'bc_types.fi'

  integer          :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
  integer          :: bc(3,2,*)
  integer          :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)
  
  call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
       domlo,domhi,delta,xlo,bc)

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
  
  call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
       domlo,domhi,delta,xlo,bc)

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
  
  call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
       domlo,domhi,delta,xlo,bc)

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

  call filcc(react,react_l1,react_l2,react_l3,react_h1,react_h2,react_h3, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_reactfill
