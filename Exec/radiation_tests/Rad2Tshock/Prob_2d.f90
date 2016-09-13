subroutine probinit(init, name, namlen, problo, probhi)

  use probdata_module
  use network, only : network_init

  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)
  
  integer untin,i
  
  ! build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter ::  maxlen = 256
  character probin*(maxlen)

  if (namlen > maxlen) call bl_error("probin file name too long")
    
  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  call network_init()
    
  ! read namelists -- this will override any defaults
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin,fortin)
  close(unit=untin)

  xmin = problo(1)
  xmax = probhi(1)

  ymin = problo(2)
  ymax = probhi(2)
  
end subroutine probinit

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
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UFX, UTEMP
  use network, only : nspec, naux
  use eos_module
  
  implicit none
  
  integer :: level, nscal
  integer :: lo(2), hi(2)
  integer :: state_l1,state_l2,state_h1,state_h2
  double precision :: state(state_l1:state_h1,state_l2:state_h2, NVAR)
  double precision :: time, delta(2)
  double precision :: xlo(2), xhi(2)
  
  integer :: i,j
  double precision :: xcell, rhoInv
  type(eos_t) :: eos_state

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
  
        xcell = xmin  + delta(1) * (dble(i) + 0.5d0)

        if (xcell < 0.d0) then
           state(i,j,URHO) = rho0
        
           ! set the composition to be all in the first species
           state(i,j,UFS:UFS-1+nspec) = 0.d0
           state(i,j,UFS  ) = state(i,j,URHO)

           state(i,j,UTEMP) = T0
           state(i,j,UMX) = rho0*v0
           state(i,j,UMY) = 0.0d0
        else
           state(i,j,URHO) = rho1
        
           ! set the composition to be all in the first species
           state(i,j,UFS:UFS-1+nspec) = 0.d0
           state(i,j,UFS  ) = state(i,j,URHO)

           state(i,j,UTEMP) = T1
           state(i,j,UMX) = rho1*v1
           state(i,j,UMY) = 0.0d0           
        end if
 
        if (naux > 0) then
           state(i,j,UFX) = state(i,j,URHO)        
        end if

        eos_state % rho = state(i,j,URHO)
        eos_state % T   = state(i,j,UTEMP)

        rhoInv = 1.d0 / state(i,j,URHO)
        eos_state % xn  = state(i,j,UFS:UFS+nspec-1) * rhoInv
        eos_state % aux = state(i,j,UFX:UFX+naux-1) * rhoInv

        call eos(eos_input_rt, eos_state)
     
        state(i,j,UEINT) = state(i,j,URHO) * eos_state % e
        state(i,j,UEDEN) = state(i,j,UEINT) + &
             0.5d0*(state(i,j,UMX)**2 + state(i,j,UMY)**2)/state(i,j,URHO)                
     enddo
  enddo
  
end subroutine ca_initdata


! ::: 
! ::: -----------------------------------------------------------
! :::
subroutine ca_initrad(level,time,lo,hi,nrad, &
                      rad_state, &
                      rad_state_l1,rad_state_l2, &
                      rad_state_h1,rad_state_h2, &
                      delta,xlo,xhi)

  use probdata_module
  use fundamental_constants_module, only: a_rad
  use rad_params_module, only : xnu
  use blackbody_module, only : BGroup
  
  implicit none
  integer :: level, nrad
  integer :: lo(2), hi(2)
  integer :: rad_state_l1,rad_state_l2,rad_state_h1,rad_state_h2
  double precision :: xlo(2), xhi(2), time, delta(2)
  double precision ::  rad_state(rad_state_l1:rad_state_h1,rad_state_l2:rad_state_h2, 0:nrad-1)

  ! local variables
  integer :: i, j, igroup
  double precision xcell, t

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        xcell = xmin + delta(1) * (dble(i) + 0.5d0)
   
        if (xcell < 0.d0) then
           T = T0
        else
           T = T1
        end if
     
        if (nrad == 1) then
           rad_state(i,j,:) = a_rad*T**4
        else
           do igroup=0,nrad-1
              rad_state(i,j,igroup) = BGroup(T, xnu(igroup), xnu(igroup+1))
           end do
        end if
        
     enddo
  enddo
  
end subroutine ca_initrad


