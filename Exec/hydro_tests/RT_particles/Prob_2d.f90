subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use eos_module, only : gamma_const
  use bl_error_module
  use bl_fort_module, only : rt => c_real
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  real(rt)         :: problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ frac, &
       rho_1, rho_2, p0_base

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call bl_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
  ! set namelist defaults here
  frac = 0.5e0_rt
  rho_1 = 1.0e0_rt
  rho_2 = 2.0e0_rt
  p0_base = 5.0e0_rt

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)


  ! set local variable defaults
  split(1) = frac*(problo(1)+probhi(1))
  split(2) = frac*(problo(2)+probhi(2))
  
  L_x = probhi(1) - problo(1)

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
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, &
       UEDEN, UEINT, UFS, UTEMP, small_temp
  use bl_constants_module, only: ZERO, HALF, M_PI
  use eos_module, only : gamma_const
  
  use bl_fort_module, only : rt => c_real
  implicit none
        
  integer :: level, nscal
  integer :: lo(2), hi(2)
  integer :: state_l1,state_l2,state_h1,state_h2
  real(rt)         :: xlo(2), xhi(2), time, delta(2)
  real(rt)         :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
  
  integer :: i,j
  real(rt)         :: x,y,pres,presmid,pertheight
  
  presmid  = p0_base - rho_1*split(2)
        
  state(:,:,UMX)   = ZERO
  state(:,:,UMY)   = ZERO
  state(:,:,UTEMP) = small_temp

  do j = lo(2), hi(2)
     y = (j+HALF)*delta(2)

     do i = lo(1), hi(1)
        
        if (y .lt. split(2)) then
           pres = p0_base - rho_1*y
           state(i,j,UEDEN) = pres / (gamma_const - 1.0e0_rt)
           state(i,j,UEINT) = pres / (gamma_const - 1.0e0_rt)
        else
           pres = presmid - rho_2*(y-split(2))
           state(i,j,UEDEN) = pres / (gamma_const - 1.0e0_rt)
           state(i,j,UEINT) = pres / (gamma_const - 1.0e0_rt)
        end if
        
     enddo
  enddo
        
  do j = lo(2), hi(2)
     y = (j+HALF)*delta(2)

     do i = lo(1), hi(1)
        x = (i+HALF)*delta(1)

        ! we explicitly make the perturbation symmetric here
        ! -- this prevents the RT from bending.
        pertheight = 0.01e0_rt*HALF*(cos(2.0e0_rt*M_PI*x/L_x) + &
                                  cos(2.0e0_rt*M_PI*(L_x-x)/L_x)) + 0.5e0_rt
        state(i,j,URHO) = rho_1 + ((rho_2-rho_1)/2.0e0_rt)* &
             (1+tanh((y-pertheight)/0.005e0_rt))
        state(i,j,UFS) = state(i,j,URHO)
        
     enddo
  enddo

end subroutine ca_initdata
