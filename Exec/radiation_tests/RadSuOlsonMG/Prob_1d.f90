
subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use network, only : network_init
  use fundamental_constants_module, only : c_light, a_rad
  use amrex_fort_module, only : rt => c_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(1), probhi(1)
  
  integer untin,i
  
!  namelist /fortin/
  
  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer maxlen
  parameter (maxlen=256)
  character probin*(maxlen)
  
  call network_init()
  
  if (namlen .gt. maxlen) then
     write(6,*) 'probin file name too long'
     stop
  end if
  
  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
  
  ! set defaults
  x0 = 0.5e0_rt
  tau0 = 10.0e0_rt
  Q = 1.e0_rt
  Temp0 = 1.e6_rt
  kapbar = 1.e0_rt
  epsilon = 1.e0_rt
  p0 = 0.5e0_rt
  p1 = 0.5e0_rt

  xmin = problo(1)
  xmax = probhi(1)

  ! ! Read namelists
  ! untin = 9
  ! open(untin,file=probin(1:namlen),form='formatted',status='old')
  ! read(untin,fortin)
  ! close(unit=untin)
  
  t0 = tau0 / (epsilon*c_light*kapbar)
  x0 = x0 / kapbar
  qn(0) = a_rad * Temp0**4 * c_light * kapbar * Q * p0
  qn(1) = a_rad * Temp0**4 * c_light * kapbar * Q * p1

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
     state,state_l1,state_h1, &
     delta,xlo,xhi)
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UFX, UTEMP
  use network, only : nspec, naux
  
  use amrex_fort_module, only : rt => c_real
  implicit none
  
  integer :: level, nscal
  integer :: lo(1), hi(1)
  integer :: state_l1,state_h1
  real(rt)         :: state(state_l1:state_h1,NVAR)
  real(rt)         :: time, delta(1)
  real(rt)         :: xlo(1), xhi(1)
  
  integer :: i
  real(rt)         :: c_v, eint

  do i = lo(1), hi(1)
  
     state(i,URHO) = 1.0
     state(i,UTEMP) = 0.0
     state(i,UMX) = 0.0
        
     ! set the composition to be all in the first species
     state(i,UFS:UFS-1+nspec) = 0.e0_rt
     state(i,UFS  ) = state(i,URHO)
     if (naux > 0) then
        state(i,UFX) = state(i,URHO)
     end if
     
     state(i,UEINT) = 0.0
     state(i,UEDEN) = 0.0

  enddo
  
end subroutine ca_initdata


! ::: 
! ::: -----------------------------------------------------------
! :::
subroutine ca_initrad(level,time,lo,hi,nrad, &
     rad_state,rad_state_l1, &
     rad_state_h1, &
     delta,xlo,xhi)

  use probdata_module
  
  use amrex_fort_module, only : rt => c_real
  implicit none
  integer :: level, nrad
  integer :: lo(1), hi(1)
  integer :: rad_state_l1,rad_state_h1
  real(rt)         :: xlo(1), xhi(1), time, delta(1)
  real(rt)         ::  rad_state(rad_state_l1:rad_state_h1, nrad)

  ! local variables
  integer :: i,j
             
  do i = lo(1), hi(1)
     rad_state(i,:) = 0.0
  enddo
        
end subroutine ca_initrad


