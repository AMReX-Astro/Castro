subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use castro_error_module
  use amrex_constants_module
  use probdata_module
  use actual_eos_module, only : gamma_const

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

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
subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       delta, xlo, xhi)

  use amrex_constants_module
  use probdata_module
  use actual_eos_module, only : gamma_const
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS
  use prob_params_module, only : problo
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: rho, u, v
  real(rt) :: x, y, r_1, r_2, vel_theta_1, vel_theta_2
  real(rt) :: cos_theta_1, sin_theta_1, cos_theta_2, sin_theta_2

  integer i, j, k


  ! density
  rho = rho_0

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        y = problo(2) + (dble(j)+HALF)*delta(2)

        do i = lo(1), hi(1)
           x = problo(1) + (dble(i)+HALF)*delta(1)

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
           state(i,j,k,UFS) = 1.0e0_rt

           ! momentum field
           state(i,j,k,UMX) = rho * u
           state(i,j,k,UMY) = rho * v
           state(i,j,k,UMZ) = 0.0e0_rt

           ! density
           state(i,j,k,URHO) = rho

           ! internal energy
           state(i,j,k,UEINT) = p_ref / (gamma_const - 1.e0_rt)

           ! Total energy
           state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * &
                sum(state(i,j,k,UMX:UMZ)**2) / state(i,j,k,URHO)

           ! Convert mass fractions to conserved quantity
           state(i,j,k,UFS) = state(i,j,k,UFS) * rho

        end do
     end do
  end do

end subroutine ca_initdata
