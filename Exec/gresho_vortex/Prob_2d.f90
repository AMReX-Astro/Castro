! This sets up the Gresho vortex problem as described in 
! Miczek, Roeple, and Edelmann 2015
!
! By choosing the reference pressure, p0, we can specify the
! Mach number

subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use prob_params_module, only : center
  use bl_constants_module
  use bl_error_module

  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(2), probhi(2)

  integer :: untin,i

  namelist /fortin/ p0, rho0, t_r, nsub

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) then
     call bl_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! Set namelist defaults
  p0 = 1.0
  rho0 = 1.0
  t_r = 1.0
  nsub = 4

  ! read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! problem center
  center(1) = (problo(1)+probhi(1))/2.d0
  center(2) = (problo(2)+probhi(2))/2.d0

  ! characteristic scales
  x_r = probhi(1) - problo(1)
  q_r = 0.4_dp_t*M_pi*x_r/t_r

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

  use probdata_module
  use eos_module, only : gamma_const
  use bl_constants_module
  use meth_params_module , only: NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS
  use prob_params_module, only : center, problo
  use network, only: nspec

  implicit none

  integer :: level, nscal
  integer :: lo(2), hi(2)
  integer :: state_l1,state_l2,state_h1,state_h2
  double precision :: xlo(2), xhi(2), time, delta(2)
  double precision :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

  double precision :: xl, yl, xx, yy
  double precision :: r
  double precision :: reint, p, u_phi, u_tot

  integer :: i,j,ii,jj

  do j = lo(2), hi(2)
     yl = problo(2) + delta(2)*dble(j)

     do i = lo(1), hi(1)
        xl = problo(1) + delta(1)*dble(i)

        reint = ZERO
        u_tot = ZERO

        do jj = 0, nsub-1
           yy = yl + delta(2)*dble(jj + HALF)/nsub

           do ii = 0, nsub-1
              xx = xl + delta(1)*dble(ii + HALF)/nsub

              r = sqrt((xx - center(1))**2 + (yy - center(2))**2)

              if (r < 0.2_dp_t) then
                 u_phi = FIVE*r
                 p = p0 + 12.5_dp_t*r**2

              else if (r < 0.4_dp_t) then
                 u_phi = TWO - FIVE*r
                 p = p0 + 12.5_dp_t*r**2 + FOUR*(ONE - FIVE*r - log(0.2_dp_t) + log(r))

              else
                 u_phi = ZERO
                 p = p0 - TWO + FOUR*log(TWO)
              endif

              u_tot = u_tot + u_phi
              reint = reint + p/(gamma_const - ONE)

           enddo
        enddo

        u_phi = u_tot/(nsub*nsub)
        reint = reint/(nsub*nsub)

        state(i,j,URHO) = rho0

        ! velocity is based on the reference velocity, q_r
        state(i,j,UMX) = rho0*q_r*u_phi*((xx-center(1))/r)  ! cos(phi) = x/r
        state(i,j,UMY) = rho0*q_r*u_phi*((yy-center(2))/r)  ! sin(phi) = y/r

        state(i,j,UEDEN) = reint +  &
             0.5d0*(state(i,j,UMX)**2/state(i,j,URHO) + &
                    state(i,j,UMY)**2/state(i,j,URHO))

        state(i,j,UEINT) = reint

        state(i,j,UFS:UFS-1+nspec) = ZERO
        state(i,j,UFS) = state(i,j,URHO)

     enddo
  enddo


end subroutine ca_initdata
