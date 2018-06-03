subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_error_module
  use amrex_constants_module
  use probdata_module
  use actual_eos_module, only : gamma_const

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  
  integer init, namlen
  integer name(namlen)
  real(rt)         problo(2), probhi(2)
  
  integer untin,i

  namelist /fortin/ p_ref, r_0, mach, ratio_c, r_circ

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) call amrex_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
  ! set namelist defaults

  ! These values are based on Lee & Koo, 1995 AIAA Journal, Figure 6

  ! Define reference pressure
  p_ref = 1.0e0_rt

  ! Define r_0
  r_0   = 0.25e0_rt

  ! Define rotating mach
  mach  = 0.0796e0_rt

  ! Define ratio_c = r_c/r_0
  ratio_c = 0.15e0_rt

  ! Define r_circ = circ/r_0*c_0
  r_circ  = 1.0e0_rt

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! Define rho_0
  rho_0 = p_ref**(1e0_rt/gamma_const)

  ! Define c_0
  c_0   = sqrt(gamma_const*p_ref/rho_0)

  ! Define r_c, radius of each vortex
  r_c   = ratio_c*r_0

  ! Define circ
  circ  = r_circ*r_0*c_0 !4e0_rt*M_PI*r_0*c_0*mach
 
  ! Center of first vortex
  x_c1  = HALF*probhi(1)
  y_c1  = HALF*probhi(2) + r_0
 
  ! Center of second vortex
  x_c2  = HALF*probhi(1)
  y_c2  = HALF*probhi(2) - r_0
  
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
  use amrex_constants_module
  use probdata_module
  use actual_eos_module, only : gamma_const
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS
  
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  
  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  real(rt)         xlo(2), xhi(2), time, delta(2)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,NVAR)
  
  real(rt)         rho, u,v
  real(rt)         x, y, r_1, r_2, vel_theta_1, vel_theta_2
  real(rt)         cos_theta_1, sin_theta_1, cos_theta_2, sin_theta_2
  integer i,j
  
  
  ! density
  rho = rho_0

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        x   = (dble(i)+HALF)*delta(1)
        y   = (dble(j)+HALF)*delta(2)

        r_1 = sqrt( (x-x_c1)**2 + (y-y_c1)**2 )
        r_2 = sqrt( (x-x_c2)**2 + (y-y_c2)**2 )

        vel_theta_1 = circ * r_1 / ( 2e0_rt * M_PI * (r_c**2 + r_1**2) )
        vel_theta_2 = circ * r_2 / ( 2e0_rt * M_PI * (r_c**2 + r_2**2) )

        sin_theta_1 = (y-y_c1) / r_1
        cos_theta_1 = (x-x_c1) / r_1

        sin_theta_2 = (y-y_c2) / r_2
        cos_theta_2 = (x-x_c2) / r_2

        u =   vel_theta_1 * sin_theta_1 + vel_theta_2 * sin_theta_2 
        v = - vel_theta_1 * cos_theta_1 - vel_theta_2 * cos_theta_2

        ! single species for all zones
        state(i,j,UFS) = 1.0e0_rt
           
        ! momentum field
        state(i,j,UMX) = rho * u
        state(i,j,UMY) = rho * v
        
        ! density       
        state(i,j,URHO) = rho
       
        ! internal energy
        state(i,j,UEINT) = p_ref / (gamma_const - 1.e0_rt)
        
        ! Total energy
        state(i,j,UEDEN) = state(i,j,UEINT) + HALF * &
             (state(i,j,UMX)**2 + state(i,j,UMY)**2) / state(i,j,URHO)
        
        ! Convert mass fractions to conserved quantity
        state(i,j,UFS) = state(i,j,UFS) * rho
        
        
     enddo
  enddo

end subroutine ca_initdata
