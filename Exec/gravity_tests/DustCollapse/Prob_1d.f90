subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_error_module
  use probdata_module
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rp
  use meth_params_module, only: small_temp
  use prob_params_module, only : center
  use network, only : nspec
  use amrex_fort_module, only : rt => amrex_real
  implicit none 

  integer :: init,namlen,untin,i,k
  integer :: name(namlen)

  real(rt)         :: center_x, center_y, center_z
  real(rt)         :: problo(1), probhi(1)

  type (eos_t) :: eos_state

  namelist /fortin/ &
       rho_0, r_0, p_0, rho_ambient, smooth_delta, &
       center_x, center_y, center_z

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
  integer, parameter :: maxlen = 127
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) call amrex_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
! set namelist defaults

  rho_0 = 1.e9_rt
  r_0 = 6.5e8_rt
  p_0 = 1.e10_rt
  rho_ambient = 1.e0_rt
  smooth_delta = 1.e-5_rt

!     Read namelists in probin file
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! in 1-d spherical, the lower domain boundary should be the origin
  center(1) = 0.0e0_rt

  xmin = problo(1)
  if (xmin /= 0.e0_rt) call amrex_error("ERROR: xmin should be 0!")


  xmax = probhi(1)

  ! set the composition to be uniform
  allocate(X_0(nspec))
  
  X_0(:) = 0.0
  X_0(1) = 1.0

  ! get the ambient temperature and sphere temperature, T_0

  eos_state % rho = rho_0
  eos_state % p   = p_0
  eos_state % xn  = x_0
  eos_state % T   = small_temp ! Initial guess for the EOS

  call eos(eos_input_rp, eos_state)

  T_0 = eos_state % T
  
  eos_state % rho = rho_ambient
  
  call eos(eos_input_rp, eos_state)

  T_ambient = eos_state % T

end subroutine amrex_probinit


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_l1,state_h1,delta,xlo,xhi)

  use probdata_module
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rp
  use network, only : nspec
  use interpolate_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS, small_temp
  use prob_params_module, only : center

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  
  integer          :: level, nscal
  integer          :: lo(1), hi(1)
  integer          :: state_l1,state_h1
  real(rt)         :: state(state_l1:state_h1,NVAR)
  real(rt)         :: time, delta(1)
  real(rt)         :: xlo(1), xhi(1)
  
  real(rt)         :: xl,xx,dist,pres,eint,temp,avg_rho, rho_n
  real(rt)         :: dx_sub
  integer          :: i,ii,n

  type (eos_t) :: eos_state

  integer, parameter :: nsub = 5

  dx_sub = delta(1)/dble(nsub)
  
  do i = lo(1), hi(1)
     xl = dble(i) * delta(1)
      
     avg_rho = 0.e0_rt

     do ii = 0, nsub-1

        xx = xl + (dble(ii) + 0.5e0_rt) * dx_sub
        
        dist = (xx-center(1))
        
        ! use a tanh profile to smooth the transition between rho_0                                    
        ! and rho_ambient                                                                              
        rho_n = rho_0 - (rho_0 - rho_ambient)* &
             0.5e0_rt * (1.e0_rt + tanh((dist - r_0)/smooth_delta))
                
        avg_rho = avg_rho + rho_n
        
     enddo

     state(i,URHO) = avg_rho/dble(nsub)

     eos_state % rho = state(i,URHO)
     eos_state % p   = p_0
     eos_state % T   = small_temp ! Initial guess for the EOS
     eos_state % xn  = X_0

     call eos(eos_input_rp, eos_state)

     temp = eos_state % T
     eint = eos_state % e
     
     state(i,UTEMP) = temp
     state(i,UMX) = 0.e0_rt
     state(i,UEDEN) = state(i,URHO) * eint
     state(i,UEINT) = state(i,URHO) * eint
     state(i,UFS:UFS+nspec-1) = state(i,URHO) * X_0(1:nspec)

  enddo
    
end subroutine ca_initdata

