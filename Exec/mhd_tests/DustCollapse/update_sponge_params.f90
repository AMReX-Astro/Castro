subroutine update_sponge_params(time) bind(C)

  use sponge_module
  use probdata_module, only: r_old_s
  use castro_error_module, only: castro_error

  use amrex_fort_module, only : rt => amrex_real
  implicit none    
  
  real(rt)        , intent(in) :: time

  logical          :: converged
  real(rt)         :: r, dr, r_guess, tol = 1.e-6_rt
  integer          :: n, max_iter = 25
  real(rt)         :: f, dfdr
  
  ! Find the analytic solution of the radius of the collapsing
  ! object via Newton-Raphson iteration.

  converged = .false.
  r_guess = 0.95*r_old_s

  do n = 1, max_iter
     dr = -f(r_guess,time) / dfdr(r_guess,time)

     ! We have lots of sqrt(1 - r/r_0), so r(t) has to be less than r_0.
     if (r_guess + dr > 6.5e8_rt) then
        r_guess = 0.5e0_rt*(r_guess + 6.5e8_rt)
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
     call castro_error("Newton iterations failed to converge in update_sponge_params.")
  endif
  
  sponge_lower_radius = r + 2.5e7_rt
  sponge_upper_radius = r + 5.0e7_rt
  sponge_timescale    = 1.0e-3_rt
  
end subroutine update_sponge_params



! The analytic solution from Colgate and White, Eq. 5 (with everything moved
! onto the LHS).  We use 1.e-20_rt to prevent against zeros in some places.

double precision function f(r,t)

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , intent(in) :: r, t

  f = sqrt(8.e0_rt*3.14159265358979323846e0_rt*6.67428e-8_rt*1.e9_rt/3.e0_rt)*t - &
      sqrt(1.e0_rt - r/6.5e8_rt)*sqrt(r/6.5e8_rt) - asin(sqrt(1.0 - r/6.5e8_rt + 1.e-20_rt))

end function f



! The derivative (wrt x) of the analytic function f, defined above.
! We use 1.e-20_rt to prevent against zeros in some places.

double precision function dfdr(r,t)

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , intent(in) :: r, t
    
  dfdr = (0.5e0_rt/6.5e8_rt)*(-sqrt(1.e0_rt - r/6.5e8_rt + 1.e-20_rt)/sqrt(r/6.5e8_rt) + &
         sqrt(r/6.5e8_rt)/sqrt(1 - r/6.5e8_rt + 1.e-20_rt) + &
         1.e0_rt/(sqrt(r/6.5e8_rt)*sqrt(1.e0_rt - r/6.5e8_rt + 1.e-20_rt)))

end function dfdr
