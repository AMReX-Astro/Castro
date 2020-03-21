module rad_util_module

  use castro_error_module
  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real
  use rad_params_module, only : ngroups

  implicit none

contains

  function FLDlambda(r, limiter) result (lambda)

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    real(rt), intent(in) :: r
    integer,  intent(in) :: limiter
    real(rt) :: lambda

    !$gpu

    if (limiter .eq. 0) then
       ! no limiter
       lambda = 1.e0_rt/3.e0_rt

    else if (limiter < 10) then
       ! approximate LP
       lambda = (2.e0_rt + r) / (6.e0_rt + r * (3.e0_rt + r))

    else if (limiter < 20) then
       ! Bruenn
       lambda = 1.e0_rt / (3.e0_rt + r)

    else if (limiter < 30) then
       ! Larsen's square root
       lambda = 1.e0_rt / sqrt(9.e0_rt + r**2)

    else if (limiter < 40) then 
       ! Minerbo
       if (r .lt. 1.5e0_rt) then
          lambda = 2.e0_rt/(3.e0_rt + sqrt(9.e0_rt+12.e0_rt*r**2))
       else 
          lambda = 1.e0_rt/(1.e0_rt+r+sqrt(1.e0_rt+2.e0_rt*r))
       end if

    else
#ifndef AMREX_USE_GPU
       print *, "limiter = ", limiter
       call castro_error("Unknown limiter type")
#endif
    endif
  end function FLDlambda

end module rad_util_module
