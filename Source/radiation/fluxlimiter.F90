module fluxlimiter_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, allocatable, save :: limiter, closure

  real(rt), parameter :: OneThird = 1.e0_rt/3.e0_rt, TwoThirds=2.e0_rt/3.e0_rt, &
                         OneSixth=1.e0_rt/6.e0_rt, TwoNinths=2.e0_rt/9.e0_rt, &
                         FiveNinths=5.e0_rt/9.e0_rt

#ifdef AMREX_USE_CUDA
  attributes(managed) :: limiter, closure
#endif

contains

  subroutine init_fluxlimiter_module(limiter_in, closure_in)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: limiter_in, closure_in

    allocate(limiter)
    allocate(closure)

    limiter = limiter_in
    closure = closure_in

  end subroutine init_fluxlimiter_module

  function FLDlambda(r, limiter_in)

    use amrex_fort_module, only: rt => amrex_real

    integer,  intent(in), optional :: limiter_in
    real(rt), intent(in) :: r

    real(rt) :: FLDlambda
    integer :: limiter_local

    if (present(limiter_in))then
       limiter_local = limiter_in
    else
       limiter_local = limiter
    end if

    if (limiter_local .eq. 0) then ! no limiter
       FLDlambda = OneThird
    else if (limiter_local < 10) then  ! approximate LP, [123]
       FLDlambda = (2.e0_rt + r) / (6.e0_rt + r * (3.e0_rt + r))
    else if (limiter_local < 20) then  ! Bruenn, 1[123]
       FLDlambda = 1.e0_rt / (3.e0_rt + r)
    else if (limiter_local < 30) then  ! Larsen's square root, 2[123]
       FLDlambda = 1.e0_rt / sqrt(9.e0_rt + r**2)
    else if (limiter_local < 40) then  ! Minerbo, 3[123]
       if (r .lt. 1.5e0_rt) then
          FLDlambda = 2.e0_rt/(3.e0_rt + sqrt(9.e0_rt+12.e0_rt*r**2))
       else 
          FLDlambda = 1.e0_rt/(1.e0_rt+r+sqrt(1.e0_rt+2.e0_rt*r))
       end if
    else
       print *, "Unknown limiter ", limiter_local
       stop
    endif

  end function FLDlambda

  function Edd_factor(lambda) result(f)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: lambda
    real(rt) :: f

    !$gpu

    if (closure .eq. 0) then
       f = lambda
    else if (closure .eq. 1) then
       f = OneThird
    else if (closure .eq. 2) then
       f = 1.e0_rt - 2.e0_rt * lambda
    else if (closure .eq. 3) then ! lambda + (lambda*R)**2
       if (limiter .eq. 0) then ! no limiter
          f = OneThird
       else if (limiter < 10) then  ! approximate LP, [123]
          f = lambda + 0.25e0_rt*(max(0.0_rt, (1.e0_rt-3.e0_rt*lambda)) + &
              sqrt(max(0.0_rt, (1.e0_rt-3.e0_rt*lambda))*(1.e0_rt+5.e0_rt*lambda)))**2
       else if (limiter < 20) then  ! Bruenn, 1[123]
          f = 1.0e0_rt - 5.e0_rt*lambda + 9.e0_rt*lambda**2
       else if (limiter < 30) then  ! Larsen's square root, 2[123]
          f = 1.0e0_rt + lambda - 9.e0_rt*lambda**2
       else if (limiter < 40) then  ! Minerbo
          if (lambda .gt. TwoNinths) then
             f = OneThird
          else
             f = 1.e0_rt + 3.e0_rt*lambda - 2.e0_rt*sqrt(2.e0_rt*lambda)
          end if
       else
          print *, "Unknown limiter ", limiter
          stop
       endif
    else if (closure .eq. 4) then ! 1/3 + 2/3*(lambda*R)**2
       if (limiter .eq. 0) then ! no limiter
          f = OneThird
       else if (limiter < 10) then  ! approximate LP, [123]
          f = OneThird + OneSixth*(max(0.0_rt, (1.e0_rt-3.e0_rt*lambda)) + &
               sqrt(max(0.0_rt, (1.e0_rt-3.e0_rt*lambda))*(1.e0_rt+5.e0_rt*lambda)))**2
       else if (limiter < 20) then  ! Bruenn, 1[123]
          f = OneThird + TwoThirds*(1.0e0_rt - 6.e0_rt*lambda + 9.e0_rt*lambda**2)
       else if (limiter < 30) then  ! Larsen's square root, 2[123]
          f = OneThird + TwoThirds*(1.0e0_rt - 9.e0_rt*lambda**2)
       else if (limiter < 40) then  ! Minerbo
          if (lambda .gt. TwoNinths) then
             f = FiveNinths - TwoThirds*lambda
          else
             f = OneThird + TwoThirds*(1.e0_rt + 2.e0_rt*lambda - 2.e0_rt*sqrt(2.e0_rt*lambda))
          end if
       else
          print *, "Unknown limiter ", limiter
          stop
       endif
    end if

  end function Edd_factor

  function FLDalpha(lam) result(alpha)

    use amrex_fort_module, only: rt => amrex_real

    real(rt), intent(in) :: lam
    real(rt) :: alpha, R, omtl, cr

    omtl = max(0.e0_rt, 1.e0_rt-3.e0_rt*lam)

    if (limiter .eq. 0) then ! no limiter
       R = 0.e0_rt
    else if (limiter < 10) then  ! approximate LP, [123]
       R = (omtl + sqrt(omtl*(1.e0_rt+5.e0_rt*lam))) / (2.e0_rt*lam+1.e-50_rt)
    else if (limiter < 20) then  ! Bruenn, 1[123]
       R =omtl/(lam+1.e-50_rt)
    else if (limiter < 30) then  ! Larsen's square root, 2[123]
       R = sqrt(omtl*(1.e0_rt+3.e0_rt*lam)) / (lam+1.e-50_rt)
    else if (limiter < 40) then  ! Minerbo
       if (lam .gt. TwoNinths) then
          R = sqrt(omtl/3.e0_rt) / (lam+1.e-50_rt)
       else
          R = 1.e0_rt/(lam+1.e-50_rt) - sqrt(2.e0_rt/(lam+1.e-50_rt))
       end if
    else
       print *, "Unknown limiter ", limiter
       stop
    endif

    if (R .lt. 1.e-6_rt) then
       alpha = 0.25e0_rt
    else if ( R .gt. 300.e0_rt) then
       alpha = 0.5e0_rt
    else
       cr = cosh(R)
       alpha = cr*log(cr) / (2.e0_rt*R*sinh(R))
    end if

  end function FLDalpha

end module fluxlimiter_module

subroutine ca_initfluxlimiter(limiter, closure)  bind(C, name="ca_initfluxlimiter")

  use fluxlimiter_module, only: init_fluxlimiter_module
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in) :: limiter, closure

  call init_fluxlimiter_module(limiter, closure)

end subroutine ca_initfluxlimiter
