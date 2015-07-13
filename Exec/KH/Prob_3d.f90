subroutine PROBINIT (init,name,namlen,problo,probhi)

   use probdata_module
   use bl_constants_module
   use fundamental_constants_module
   use meth_params_module, only: small_temp, small_pres, small_dens
   use eos_module

   implicit none

   integer :: init, namlen
   integer :: name(namlen)
   double precision :: problo(3), probhi(3)

   integer :: untin
   integer :: i

   namelist /fortin/ &
        rho1, rho2, pressure, problem, bulk_velocity

   integer, parameter :: maxlen=127
   character :: probin*(maxlen)
   character :: model*(maxlen)
   integer :: ipp, ierr, ipp1

   ! Temporary storage variables in case we need to switch the primary and secondary.

   integer :: ioproc

   ! For outputting -- determine if we are the IO processor
   call bl_pd_is_ioproc(ioproc)

   ! Build "probin" filename -- the name of file containing fortin namelist.
   if (namlen .gt. maxlen) then
      call bl_error("ERROR: probin file name too long")
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
                       state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                       delta,xlo,xhi)

  use probdata_module
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
       UEDEN, UEINT, UFS, UFA
  use network, only : nspec
  use bl_constants_module
  
  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision :: xlo(3), xhi(3), time, delta(3)
  double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

  double precision :: xx,yy,zz

  type (eos_t) :: eos_state

  integer :: i,j,k,n

  double precision :: dens, velx, vely, velz
  double precision :: w0, sigma, ramp, delta_y
  double precision :: vel1, vel2

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
  endif

  velz = 0.0
  
  !$OMP PARALLEL DO PRIVATE(i, j, k, xx, yy, zz, dens, velx, vely, velz, ramp, eos_state)
  do k = lo(3), hi(3)
     zz = xlo(3) + delta(3)*dble(k-lo(3)+HALF)

     do j = lo(2), hi(2)     
        yy = xlo(2) + delta(2)*dble(j-lo(2)+HALF)

        do i = lo(1), hi(1)   
           xx = xlo(1) + delta(1)*dble(i-lo(1)+HALF)

           ! Assume the initial y-velocity represents the bulk flow
           ! which will be perturbed in the following step

           vely = bulk_velocity

           if (problem .eq. 1) then

              if (abs(yy - 0.5) < 0.25) then
                 dens = rho2
                 velx = 0.5
              else
                 dens = rho1
                 velx = -0.5
              endif

              vely = vely + w0 * sin(sine_n*M_PI*xx) * (exp(-(yy-0.25)**2/(2*sigma**2)) + exp(-(yy-0.75)**2/(2*sigma**2)))

           else if (problem .eq. 2) then

             ramp = ((ONE + exp(-TWO*(yy-0.25)/delta_y))*(ONE + exp(TWO*(yy-0.75)/delta_y)))**(-1)

             dens = rho1 + ramp * (rho2 - rho1)
             velx = vel1 + ramp * (vel2 - vel1)

             vely = vely + w0 * sin(sine_n*M_PI*xx)

           else if (problem .eq. 3) then

              if ( yy .lt. 0.25 ) then
                 dens = rho1 - (rho1 - rho2) / 2 * exp( (yy-0.25) / delta_y )
                 velx = vel1 - (vel1 - vel2) / 2 * exp( (yy-0.25) / delta_y )
              else if ( yy .le. 0.50 ) then
                 dens = rho2 + (rho1 - rho2) / 2 * exp( (0.25-yy) / delta_y )
                 velx = vel2 + (vel1 - vel2) / 2 * exp( (0.25-yy) / delta_y )
              else if ( yy .lt. 0.75 ) then
                 dens = rho2 + (rho1 - rho2) / 2 * exp( (yy-0.75) / delta_y )
                 velx = vel2 + (vel1 - vel2) / 2 * exp( (yy-0.75) / delta_y )
              else
                 dens = rho1 - (rho1 - rho2) / 2 * exp( (0.75-yy) / delta_y )
                 velx = vel1 - (vel1 - vel2) / 2 * exp( (0.75-yy) / delta_y )
              endif

              vely = vely + w0 * sin(sine_n*M_PI*xx)

           else

              call bl_error("Error: This problem choice is undefined.")

           endif

           state(i,j,k,URHO) = dens
           state(i,j,k,UMX)  = dens * velx
           state(i,j,k,UMY)  = dens * vely
           state(i,j,k,UMZ)  = dens * velz

           ! Establish the thermodynamic quantities

           state(i,j,k,UFS:UFS-1+nspec) = 1.0d0 / nspec

           eos_state % xn  = state(i,j,k,UFS:UFS-1+nspec)
           eos_state % rho = state(i,j,k,URHO)
           eos_state % p   = pressure

           call eos(eos_input_rp, eos_state)

           state(i,j,k,UTEMP) = eos_state % T

           state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state % e
           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e

           state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + &
                HALF * dens * (velx**2 + vely**2 + velz**2)

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO  

end subroutine ca_initdata

! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                      domlo,domhi,delta,xlo,time,bc)

  use meth_params_module, only : NVAR

  implicit none
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

  integer n

  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
          adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
          domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

end subroutine ca_hypfill

! ::: -----------------------------------------------------------

subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                      domlo,domhi,delta,xlo,time,bc)

  implicit none
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,domlo,domhi,delta,xlo,bc)

end subroutine ca_denfill
