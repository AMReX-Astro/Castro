subroutine update_sponge_params(time) bind(C)

  use sponge_module
  use probdata_module, only: r_old_s
  use bl_error_module, only: bl_error

  implicit none    
  
  double precision, intent(in) :: time

  logical          :: converged
  double precision :: r, dr, r_guess, tol = 1.d-6
  integer          :: n, max_iter = 25
  double precision :: f, dfdr
  
  ! Find the analytic solution of the radius of the collapsing
  ! object via Newton-Raphson iteration.

  converged = .false.
  r_guess = 0.95*r_old_s

  do n = 1, max_iter
     dr = -f(r_guess,time) / dfdr(r_guess,time)

     ! We have lots of sqrt(1 - r/r_0), so r(t) has to be less than r_0.
     if (r_guess + dr > 6.5d8) then
        r_guess = 0.5d0*(r_guess + 6.5d8)
     else
        r_guess = r_guess + dr
     endif

     if (abs(dr/r_guess) < tol) then
        converged = .true.
        r = r_guess
        r_old_s = r_guess
        exit
     endif

  enddo
    
  if (.not. converged) then
     call bl_error("Newton iterations failed to converge in update_sponge_params.")
  endif
  
  sponge_lower_radius = r + 2.5d7
  sponge_upper_radius = r + 5.0d7
  sponge_timescale    = 1.0d-3
  
end subroutine update_sponge_params



! The analytic solution from Colgate and White, Eq. 5 (with everything moved
! onto the LHS).  We use 1.d-20 to prevent against zeros in some places.

double precision function f(r,t)

  implicit none

  double precision, intent(in) :: r, t

  f = sqrt(8.d0*3.14159265358979323846d0*6.67428d-8*1.d9/3.d0)*t - &
      sqrt(1.d0 - r/6.5d8)*sqrt(r/6.5d8) - asin(sqrt(1.0 - r/6.5d8 + 1.d-20))

end function f



! The derivative (wrt x) of the analytic function f, defined above.
! We use 1.d-20 to prevent against zeros in some places.

double precision function dfdr(r,t)

  implicit none

  double precision, intent(in) :: r, t
    
  dfdr = (0.5d0/6.5d8)*(-sqrt(1.d0 - r/6.5d8 + 1.d-20)/sqrt(r/6.5d8) + &
         sqrt(r/6.5d8)/sqrt(1 - r/6.5d8 + 1.d-20) + &
         1.d0/(sqrt(r/6.5d8)*sqrt(1.d0 - r/6.5d8 + 1.d-20)))

end function dfdr
