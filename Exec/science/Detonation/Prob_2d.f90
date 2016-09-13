subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use network   , only : network_species_index, nspec
  use bl_error_module

  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)
  
  integer untin,i

  namelist /fortin/ T_l, T_r, dens, cfrac, frac, idir, &
       denerr,  dengrad,  max_denerr_lev,  max_dengrad_lev, &
       velgrad,  max_velgrad_lev, &
       presserr, pressgrad,max_presserr_lev,max_pressgrad_lev

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
  integer, parameter :: maxlen = 256
  character probin*(maxlen)
  
  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
! set namelist defaults

  T_l = 1.d9
  T_r = 5.d7
  dens = 1.d8
  
  idir = 1                ! direction across which to jump
  frac = 0.5              ! fraction of the domain for the interface
  cfrac = 0.5

  denerr = 1.d20
  dengrad = 1.d20
  max_denerr_lev = -1
  max_dengrad_lev = -1
  
  presserr = 1.d20
  pressgrad = 1.d20
  max_presserr_lev = -1
  max_pressgrad_lev = -1
  
  velgrad = 1.d20
  max_velgrad_lev = -1
  
!     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)
  
  xmin = problo(1)
  xmax = probhi(1)

  ymin = problo(2)
  ymax = probhi(2)

  center(1) = frac*(problo(1)+probhi(1))
  center(2) = frac*(problo(2)+probhi(2))

  ! get the species indices
  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")

  if (ic12 < 0 .or. io16 < 0) then
     call bl_error("ERROR: species indices not found")
  endif


  ! make sure that the carbon fraction falls between 0 and 1
  if (cfrac > 1.d0 .or. cfrac < 0.d0) then
     call bl_error("ERROR: cfrac must fall between 0 and 1")
  endif

  ! set the default mass fractions
  allocate(xn(nspec))

  xn(:) = 0.d0
  xn(ic12) = cfrac
  xn(io16) = 1.d0 - cfrac


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
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use network, only: nspec
  use eos_module
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP

  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)
  double precision time, delta(2)
  double precision xlo(2), xhi(2)

  double precision xcen, ycen
  double precision p_temp, eint_temp
  integer i,j

  type (eos_t) :: eos_state

  do j = lo(2), hi(2)
     ycen = ymin + delta(2)*(dble(j) + 0.5d0)
         
     do i = lo(1), hi(1)
        xcen = xmin + delta(1)*(dble(i) + 0.5d0)
            
        state(i,j,URHO ) = dens

        if (xcen <= frac*(xmin + 0.5d0*(xmax-xmin))) then
           state(i,j,UTEMP) = T_l
        else
           state(i,j,UTEMP) = T_r
        endif

        state(i,j,UFS:UFS-1+nspec) = state(i,j,URHO)*xn(1:nspec)

        eos_state%rho = state(i,j,URHO)
        eos_state%T = state(i,j,UTEMP)
        eos_state%xn(:) = xn

        call eos(eos_input_rt, eos_state)

        state(i,j,UMX  ) = 0.d0
        state(i,j,UMY  ) = 0.d0
        state(i,j,UEDEN) = state(i,j,URHO)*eos_state%e
        state(i,j,UEINT) = state(i,j,URHO)*eos_state%e
        
     enddo
  enddo
  
end subroutine ca_initdata
