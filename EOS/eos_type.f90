module eos_type_module
  
  use bl_types
  use network

  implicit none

  type eos_t

     real (kind=dp_t) :: rho
     real (kind=dp_t) :: T
     real (kind=dp_t) :: xn(nspec)
     real (kind=dp_t) :: p
     real (kind=dp_t) :: h
     real (kind=dp_t) :: e
     real (kind=dp_t) :: cv
     real (kind=dp_t) :: cp
     real (kind=dp_t) :: xne
     real (kind=dp_t) :: eta
     real (kind=dp_t) :: pele
     real (kind=dp_t) :: dpdT
     real (kind=dp_t) :: dpdr
     real (kind=dp_t) :: dedT
     real (kind=dp_t) :: dedr
     real (kind=dp_t) :: dpdX(nspec)
     real (kind=dp_t) :: dhdX(nspec)
     real (kind=dp_t) :: gam1
     real (kind=dp_t) :: cs
     real (kind=dp_t) :: s
     real (kind=dp_t) :: dsdT
     real (kind=dp_t) :: dsdr         

  end type eos_t

end module eos_type_module
