subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use network   , only : network_species_index, nspec
  use bl_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(1), probhi(1)

  integer untin,i

  namelist /fortin/ T_l, T_r, dens, cfrac, frac, idir, &
       denerr,  dengrad,  max_denerr_lev,  max_dengrad_lev, &
       velgrad,  max_velgrad_lev, &
       presserr,pressgrad,max_presserr_lev,max_pressgrad_lev, &
       temperr,tempgrad,max_temperr_lev,max_tempgrad_lev

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
  
  T_l = 1.e9_rt
  T_r = 5.e7_rt
  dens = 1.e8_rt

  idir = 1                ! direction across which to jump
  frac = 0.5              ! fraction of the domain for the interface
  cfrac = 0.5

  denerr = 1.e20_rt
  dengrad = 1.e20_rt
  max_denerr_lev = -1
  max_dengrad_lev = -1

  presserr = 1.e20_rt
  pressgrad = 1.e20_rt
  max_presserr_lev = -1
  max_pressgrad_lev = -1

  velgrad = 1.e20_rt
  max_velgrad_lev = -1
  
!     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  xmin = problo(1)
  xmax = probhi(1)

  center(1) = 0.5e0_rt*(problo(1)+probhi(1))

  ! get the species indices
  ihe4 = network_species_index("helium-4")
  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")

  if (ihe4 < 0 .or. ic12 < 0 .or. io16 < 0) then
     call bl_error("ERROR: species indices not found")
  endif


  ! make sure that the carbon fraction falls between 0 and 1
  if (cfrac > 1.e0_rt .or. cfrac < 0.e0_rt) then
     call bl_error("ERROR: cfrac must fall between 0 and 1")
  endif

  ! set the default mass fractions
  allocate(xn(nspec))

  xn(:) = 0.e0_rt
  xn(ic12) = cfrac
  xn(ihe4) = 1.e0_rt - cfrac
  
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

  use eos_module
  use network, only: nspec
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFS

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer level, nscal
  integer lo(1), hi(1)
  integer state_l1,state_h1
  real(rt)         state(state_l1:state_h1,NVAR)
  real(rt)         time, delta(1)
  real(rt)         xlo(1), xhi(1)
  
  real(rt)         xcen
  real(rt)         p_temp, eint_temp
  integer i

  real(rt)         :: L_x

  type (eos_t) :: eos_state

  L_x = xmax - xmin

  do i = lo(1), hi(1)
     xcen = xmin + delta(1)*(dble(i) + 0.5e0_rt)

     state(i,URHO ) = dens
            
     if (xcen <= xmin + frac*L_x) then
        state(i,UTEMP) = T_l
     else
        state(i,UTEMP) = T_r
     endif

     state(i,UFS:UFS-1+nspec) = state(i,URHO)*xn(1:nspec)

     eos_state%rho = state(i,URHO)
     eos_state%T = state(i,UTEMP)
     eos_state%xn(:) = xn

     call eos(eos_input_rt, eos_state)
    
     state(i,UMX  ) = 0.e0_rt
     state(i,UEDEN) = state(i,URHO)*eos_state%e  ! if vel /= 0, then KE needs to be added
     state(i,UEINT) = state(i,URHO)*eos_state%e
            
  enddo

end subroutine ca_initdata
