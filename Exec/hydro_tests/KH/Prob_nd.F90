subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

   use probdata_module
   use amrex_constants_module
   use castro_error_module
   use fundamental_constants_module
   use meth_params_module, only: small_temp, small_pres, small_dens
   
   use amrex_fort_module, only : rt => amrex_real
   implicit none

   integer :: init, namlen
   integer :: name(namlen)
   real(rt)         :: problo(3), probhi(3)

   integer :: untin
   integer :: i

   namelist /fortin/ &
        rho1, rho2, pressure, problem, bulk_velocity

   integer, parameter :: maxlen=127
   character :: probin*(maxlen)
   character :: model*(maxlen)
   integer :: ipp, ierr, ipp1

   ! Temporary storage variables in case we need to switch the primary and secondary.

   ! Build "probin" filename -- the name of file containing fortin namelist.
   if (namlen .gt. maxlen) then
      call castro_error("ERROR: probin file name too long")
   end if

   do i = 1, namlen
      probin(i:i) = char(name(i))
   end do

   ! Set namelist defaults

   problem = 2

   rho1 = 1.0
   rho2 = 2.0
   pressure = 2.5

   bulk_velocity = 0.0

   ! Read namelists -- override the defaults
   
   untin = 9 
   open(untin,file=probin(1:namlen),form='formatted',status='old')
   read(untin,fortin)
   close(unit=untin)

   ! Force a different pressure choice for problem 5

   if (problem .eq. 5) then
      pressure = 10.0
   endif

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
                       state,state_lo,state_hi, &
                       delta,xlo,xhi)

   use amrex_constants_module
   use castro_error_module

   use eos_module, only : eos
   use eos_type_module, only: eos_t, eos_input_rp
   use network, only : nspec

   use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
       UEDEN, UEINT, UFS, UFA
   use probdata_module
   use prob_params_module, only: problo, center, probhi
  
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_lo(3), state_hi(3)
  real(rt)         :: xlo(3), xhi(3), time, delta(3)
  real(rt)         :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  real(rt)         :: xx, yy, zz

  type (eos_t) :: eos_state

  integer :: i, j, k, n

  real(rt)         :: dens, velx, vely, velz
  real(rt)         :: w0, sigma, ramp, delta_y
  real(rt)         :: vel1, vel2
  real(rt)         :: y1, y2
  real(rt)         :: dye
  
  integer :: sine_n

  vel1 = -0.5
  vel2 =  0.5

  if (problem .eq. 1) then
     sine_n = 4
     w0 = 0.1
     sigma = 0.05 / 2**0.5
  else if (problem .eq. 2) then
     sine_n = 2
     w0 = 0.1
     delta_y = 0.05
  else if (problem .eq. 3) then
     sine_n = 4
     w0 = 0.01
     delta_y = 0.025
  else if (problem .eq. 4) then
     sine_n = 2
     w0 = 0.01
     delta_y = 0.025
  else if (problem .eq. 5) then
     sine_n = 2
     w0 = 0.01
     delta_y = 0.05
     sigma = 0.2     
     vel1 = ONE
     vel2 = ONE
  endif

  y1 = center(2) - (probhi(2) - problo(2)) * 0.25e0_rt
  y2 = center(2) + (probhi(2) - problo(2)) * 0.25e0_rt  

  velz = 0.0

  !$OMP PARALLEL DO PRIVATE(i, j, k, xx, yy, zz, dens, velx, vely, ramp, eos_state)
  do k = lo(3), hi(3)
     zz = xlo(3) + delta(3)*dble(k-lo(3)+HALF)

     do j = lo(2), hi(2)     
        yy = xlo(2) + delta(2)*dble(j-lo(2)+HALF)

        do i = lo(1), hi(1)   
           xx = xlo(1) + delta(1)*dble(i-lo(1)+HALF)

           ! Assume the initial y-velocity represents the bulk flow
           ! which will be perturbed in the following step

           vely = bulk_velocity
           dye = ZERO
           
           if (problem .eq. 1) then

              if (abs(yy - HALF * (y1 + y2)) < HALF * (y2 - y1)) then
                 dens = rho2
                 velx = 0.5
              else
                 dens = rho1
                 velx = -0.5
              endif

              vely = vely + w0 * sin(sine_n*M_PI*xx) * (exp(-(yy-y1)**2/(2*sigma**2)) + exp(-(yy-y2)**2/(2*sigma**2)))

           else if (problem .eq. 2) then

             ramp = ((ONE + exp(-TWO*(yy-y1)/delta_y))*(ONE + exp(TWO*(yy-y2)/delta_y)))**(-1)

             dens = rho1 + ramp * (rho2 - rho1)
             velx = vel1 + ramp * (vel2 - vel1)

             vely = vely + w0 * sin(sine_n*M_PI*xx)

           else if (problem .eq. 3 .or. problem .eq. 4) then

              if ( yy .lt. y1 ) then
                 dens = rho1 - (rho1 - rho2) / 2 * exp( (yy-y1) / delta_y )
                 velx = vel1 - (vel1 - vel2) / 2 * exp( (yy-y1) / delta_y )
              else if ( yy .le. HALF * (y1 + y2) ) then
                 dens = rho2 + (rho1 - rho2) / 2 * exp( (y1-yy) / delta_y )
                 velx = vel2 + (vel1 - vel2) / 2 * exp( (y1-yy) / delta_y )
              else if ( yy .lt. y2 ) then
                 dens = rho2 + (rho1 - rho2) / 2 * exp( (yy-y2) / delta_y )
                 velx = vel2 + (vel1 - vel2) / 2 * exp( (yy-y2) / delta_y )
              else
                 dens = rho1 - (rho1 - rho2) / 2 * exp( (y2-yy) / delta_y )
                 velx = vel1 - (vel1 - vel2) / 2 * exp( (y2-yy) / delta_y )
              endif

              vely = vely + w0 * sin(sine_n*M_PI*xx)

           else if (problem .eq. 5) then

              dens = rho1 + (rho2 - rho1) * HALF * (tanh( (yy - y1) / delta_y ) - tanh( (yy - y2) / delta_y ))
              velx = vel1 * (tanh( (yy - y1) / delta_y) - tanh( (yy - y2) / delta_y ) - ONE)
              vely = vely + w0 * sin(sine_n*M_PI*xx) * (exp(-(yy - y1)**2 / sigma**2) + exp(-(yy - y2)**2 / sigma**2))
              dye  = HALF * (tanh( (yy - y2) / delta_y) - tanh( (yy - y1) / delta_y ) + TWO)
              
           else

              call castro_error("Error: This problem choice is undefined.")

           endif

           state(i,j,k,URHO) = dens
           state(i,j,k,UMX)  = dens * velx
           state(i,j,k,UMY)  = dens * vely
           state(i,j,k,UMZ)  = dens * velz
           state(i,j,k,UFA)  = dye

           ! Establish the thermodynamic quantities

           state(i,j,k,UFS:UFS-1+nspec) = ONE / nspec

           eos_state % xn  = state(i,j,k,UFS:UFS-1+nspec)
           eos_state % rho = state(i,j,k,URHO)
           eos_state % p   = pressure

           call eos(eos_input_rp, eos_state)

           state(i,j,k,UTEMP) = eos_state % T

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * dens * (velx**2 + vely**2 + velz**2)

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO  

end subroutine ca_initdata

