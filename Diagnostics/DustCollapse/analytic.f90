! Print out the analytic solution for the homologous dust collapse problem,
! from Colgate and White 1966, ApJ, 143, 626.
!
! This will yield r(t) as a function of t

module constants_module

  implicit none

  ! fundamental constants
  double precision, parameter :: Gconst = 6.67428d-8      ! cm^3/g/s^2
  double precision, parameter :: pi = 3.14159265358979323846d0
  double precision, parameter :: SMALL = 1.d-20

  ! user-defined problem parameters
  double precision, parameter :: rho_0 = 1.d9
  double precision, parameter :: r_0   = 6.5d8

  ! we try to do equal timesteps, but eventually, the radius of the
  ! sphere is changing so quickly that we need to adjust the size of
  ! the timestep.  The code will take nstep timesteps, with the
  ! initial timestep to be (t_f - t_i)/(nstep - 1).  Once the radius
  ! begins to change significantly, we will begin halving the timestep
  ! as needed.

  integer :: nstep = 100                         ! number of points to output
  double precision, parameter :: t_i = 1.d-3     ! initial output time
  double precision, parameter :: t_f = 0.1       ! final output time

end module constants_module


program analytic

  use constants_module

  implicit none

  ! local variables
  double precision :: r, r_guess, r_old, t
  integer :: i, n
  double precision :: dt
  double precision :: dr, f, dfdr

  logical :: converged
  double precision, parameter :: TOL = 1.d-6
  integer, parameter :: max_iter = 25

  ! compute the initial timestep
  dt = (t_f - t_i)/(nstep - 1)

  r_old = r_0

  t = t_i

  ! main timestepping loop
  do i = 1, nstep

     ! find the analytic solution via Newton iteration
     converged = .false.
     r_guess = 0.95*r_old

     do n = 1, max_iter
        dr = -f(r_guess,t)/dfdr(r_guess,t)

        ! we have lots of sqrt(1 - r/r_0), so r(t) has to be less than r_0
        if (r_guess + dr > r_0) then
           r_guess = 0.5d0*(r_guess + r_0)
        else
           r_guess = r_guess + dr
        endif

        if (abs(dr/r_guess) < TOL) then
           converged = .true.
           exit
        endif

     enddo

     if (.not. converged) then
        print *, "ERROR: Newton iterations failed to converge"
        stop
     endif

     ! we converged, set the new radius of the sphere
     r = r_guess

     print *, t, r

     ! toward the end, we need smaller timesteps, if the radius is
     ! changing significantly, halve the timestep
     if (abs(r - r_old)/r_old > 0.2) then
        dt = 0.5*dt
     endif

     ! store the old radius for the timestep check next cycle
     r_old = r

     ! update the desired output time
     t = t + dt

  enddo

end program analytic




function f(r,t) result (func)

  use constants_module

  implicit none

  double precision, intent(in) :: r, t
  double precision :: func

  ! the analytic solution from Colgate and White, Eq. 5 (with everything moved
  ! onto the LHS).  We use SMALL to prevent against zeros in some places.
  func = sqrt(8.d0*pi*Gconst*rho_0/3.d0)*t - &
         sqrt(1.d0 - r/r_0)*sqrt(r/r_0) - asin(sqrt(1.0 - r/r_0 + SMALL))

  return
end function f



function dfdr(r,t) result (dfuncdr)

  use constants_module

  implicit none

  double precision, intent(in) :: r, t
  double precision :: dfuncdr

  ! the derivative (wrt x) of the analytic function f, defined above.
  ! We use SMALL to prevent against zeros in some places.
  dfuncdr = (0.5d0/r_0)*(-sqrt(1.d0 - r/r_0 + SMALL)/sqrt(r/r_0) + &
                      sqrt(r/r_0)/sqrt(1 - r/r_0 + SMALL) + &
                      1.d0/(sqrt(r/r_0)*sqrt(1.d0 - r/r_0 + SMALL)))


  return
end function dfdr



