
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use network, only : network_init
  use eos_module, only : gamma_const
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(1), probhi(1)
  
  integer untin,i
  
  namelist /fortin/ rho0, T0, v0, rho1, T1, v1

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
  
  ! set namelist defaults
  
  rho0 = 5.4588672836d-13
  T0 = 100.d0
  v0 = 235435.371882d0
  rho1 = 1.24793794736d-12
  T1 = 207.756999533d0
  v1 = 102986.727159d0
    
  xmin = problo(1)
  xmax = probhi(1)

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)
  
end subroutine PROBINIT

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
  use eos_module
  
  implicit none
  
  integer :: level, nscal
  integer :: lo(1), hi(1)
  integer :: state_l1,state_h1
  double precision :: state(state_l1:state_h1,NVAR)
  double precision :: time, delta(1)
  double precision :: xlo(1), xhi(1)
  
  integer :: i
  double precision :: xcell, rhoInv
  type(eos_t) :: eos_state

  do i = lo(1), hi(1)
  
     xcell = xlo(1) + delta(1) * (float(i-lo(1)) + 0.5d0)

     if (xcell < 0.d0) then
        state(i,URHO) = rho0
        
        ! set the composition to be all in the first species
        state(i,UFS:UFS-1+nspec) = 0.d0
        state(i,UFS  ) = state(i,URHO)

        state(i,UTEMP) = T0
        state(i,UMX) = rho0*v0
     else
        state(i,URHO) = rho1
        
        ! set the composition to be all in the first species
        state(i,UFS:UFS-1+nspec) = 0.d0
        state(i,UFS  ) = state(i,URHO)

        state(i,UTEMP) = T1
        state(i,UMX) = rho1*v1
     end if
 
     if (naux > 0) then
        state(i,UFX) = state(i,URHO)        
     end if

     eos_state % rho = state(i,URHO)
     eos_state % T   = state(i,UTEMP)

     rhoInv = 1.d0 / state(i,URHO)
     eos_state % xn  = state(i,UFS:UFS+nspec-1) * rhoInv
     eos_state % aux = state(i,UFX:UFX+naux-1) * rhoInv

     call eos(eos_input_rt, eos_state)
     
     state(i,UEINT) = state(i,URHO) * eos_state % e
     state(i,UEDEN) = state(i,UEINT) + &
          0.5*(state(i,UMX)**2)/state(i,URHO)                
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
  use fundamental_constants_module, only: a_rad
  use rad_params_module, only : xnu
  use blackbody_module, only : BGroup
  
  implicit none
  integer :: level, nrad
  integer :: lo(1), hi(1)
  integer :: rad_state_l1,rad_state_h1
  double precision :: xlo(1), xhi(1), time, delta(1)
  double precision ::  rad_state(rad_state_l1:rad_state_h1, 0:nrad-1)

  ! local variables
  integer :: i, igroup
  double precision xcell, t

  do i = lo(1), hi(1)

     xcell = xlo(1) + delta(1) * (float(i-lo(1)) + 0.5d0)
   
     if (xcell < 0.d0) then
        T = T0
     else
        T = T1
     end if
     
     if (nrad == 1) then
        rad_state(i,:) = a_rad*T**4
     else
        do igroup=0,nrad-1
           rad_state(i,igroup) = BGroup(T, xnu(igroup), xnu(igroup+1))
        end do
     end if

  enddo

end subroutine ca_initrad


