module eos_type_module
  
  use bl_types
  use network

  implicit none

  real (kind=dp_t), parameter :: init_num = -1.0d200

  type eos_t

     real (kind=dp_t) :: rho         = init_num
     real (kind=dp_t) :: T           = init_num

     real (kind=dp_t) :: p           = init_num
     real (kind=dp_t) :: e           = init_num
     real (kind=dp_t) :: h           = init_num
     real (kind=dp_t) :: s           = init_num
     real (kind=dp_t) :: dpdT        
     real (kind=dp_t) :: dpdr        
     real (kind=dp_t) :: dedT        
     real (kind=dp_t) :: dedr        
     real (kind=dp_t) :: dhdT        
     real (kind=dp_t) :: dhdr        
     real (kind=dp_t) :: dsdT        
     real (kind=dp_t) :: dsdr        
     real (kind=dp_t) :: dpde        
     real (kind=dp_t) :: dpdr_e      

     real (kind=dp_t) :: xn(nspec)   
     real (kind=dp_t) :: cv          
     real (kind=dp_t) :: cp          
     real (kind=dp_t) :: xne         
     real (kind=dp_t) :: xnp         
     real (kind=dp_t) :: eta         
     real (kind=dp_t) :: pele        
     real (kind=dp_t) :: ppos        
     real (kind=dp_t) :: dedX(nspec) 
     real (kind=dp_t) :: dpdX(nspec) 
     real (kind=dp_t) :: dhdX(nspec) 
     real (kind=dp_t) :: gam1        
     real (kind=dp_t) :: cs          

     real (kind=dp_t) :: abar        
     real (kind=dp_t) :: zbar        
     real (kind=dp_t) :: dpa          
     real (kind=dp_t) :: dpz         
     real (kind=dp_t) :: dea         
     real (kind=dp_t) :: dez         


  end type eos_t

end module eos_type_module
