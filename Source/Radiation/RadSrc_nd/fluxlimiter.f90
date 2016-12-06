module fluxlimiter_module

  implicit none

  integer, save :: limiter = 0
  integer, save :: closure = 0

  double precision, parameter :: OneThird = 1.d0/3.d0, TwoThirds=2.d0/3.d0, &
       OneSixth=1.d0/6.d0, TwoNinths=2.d0/9.d0, FiveNinths=5.d0/9.d0

  contains

    subroutine init_fluxlimiter_module(limiter_in, closure_in)
      implicit none
      integer, intent(in) :: limiter_in, closure_in
      limiter = limiter_in
      closure = closure_in
    end subroutine init_fluxlimiter_module

    function FLDlambda(r, limiter_in)
      integer, intent(in), optional :: limiter_in
      double precision, intent(in) :: r
      double precision :: FLDlambda
      integer :: limiter_local
      if (present(limiter_in))then
         limiter_local = limiter_in
      else
         limiter_local = limiter
      end if
      if (limiter_local .eq. 0) then ! no limiter
         FLDlambda = OneThird
      else if (limiter_local < 10) then  ! approximate LP, [123]
         FLDlambda = (2.d0 + r) / (6.d0 + r * (3.d0 + r))
      else if (limiter_local < 20) then  ! Bruenn, 1[123]
         FLDlambda = 1.d0 / (3.d0 + r)
      else if (limiter_local < 30) then  ! Larsen's square root, 2[123]
         FLDlambda = 1.d0 / sqrt(9.d0 + r**2)
      else if (limiter_local < 40) then  ! Minerbo, 3[123]
         if (r .lt. 1.5d0) then
            FLDlambda = 2.d0/(3.d0 + sqrt(9.d0+12.d0*r**2))
         else 
            FLDlambda = 1.d0/(1.d0+r+sqrt(1.d0+2.d0*r))
         end if
      else
         print *, "Unknown limiter ", limiter_local
         stop
      endif
    end function FLDlambda

    function Edd_factor(lambda) result(f)
      double precision, intent(in) :: lambda
      double precision :: f
      if (closure .eq. 0) then
         f = lambda
      else if (closure .eq. 1) then
         f = OneThird
      else if (closure .eq. 2) then
         f = 1.d0 - 2.d0 * lambda
      else if (closure .eq. 3) then ! lambda + (lambda*R)**2
         if (limiter .eq. 0) then ! no limiter
            f = OneThird
         else if (limiter < 10) then  ! approximate LP, [123]
            f = lambda + 0.25d0*((1.d0-3.d0*lambda) + &
                 sqrt((1.d0-3.d0*lambda)*(1.d0+5.d0*lambda)))**2
         else if (limiter < 20) then  ! Bruenn, 1[123]
            f = 1.0d0 - 5.d0*lambda + 9.d0*lambda**2
         else if (limiter < 30) then  ! Larsen's square root, 2[123]
            f = 1.0d0 + lambda - 9.d0*lambda**2
         else if (limiter < 40) then  ! Minerbo
            if (lambda .gt. TwoNinths) then
               f = OneThird
            else
               f = 1.d0 + 3.d0*lambda - 2.d0*sqrt(2.d0*lambda)
            end if
         else
            print *, "Unknown limiter ", limiter
            stop
         endif
      else if (closure .eq. 4) then ! 1/3 + 2/3*(lambda*R)**2
         if (limiter .eq. 0) then ! no limiter
            f = OneThird
         else if (limiter < 10) then  ! approximate LP, [123]
            f = OneThird + OneSixth*((1.d0-3.d0*lambda) + &
                 sqrt((1.d0-3.d0*lambda)*(1.d0+5.d0*lambda)))**2
         else if (limiter < 20) then  ! Bruenn, 1[123]
            f = OneThird + TwoThirds*(1.0d0 - 6.d0*lambda + 9.d0*lambda**2)
         else if (limiter < 30) then  ! Larsen's square root, 2[123]
            f = OneThird + TwoThirds*(1.0d0 - 9.d0*lambda**2)
         else if (limiter < 40) then  ! Minerbo
            if (lambda .gt. TwoNinths) then
               f = FiveNinths - TwoThirds*lambda
            else
               f = OneThird + TwoThirds*(1.d0 + 2.d0*lambda - 2.d0*sqrt(2.d0*lambda))
            end if
         else
            print *, "Unknown limiter ", limiter
            stop
         endif
      end if
    end function Edd_factor

    function FLDalpha(lam) result(alpha)
      double precision, intent(in) :: lam
      double precision :: alpha, R, omtl, cr
      omtl = max(0.d0, 1.d0-3.d0*lam)
      if (limiter .eq. 0) then ! no limiter
         R = 0.d0
      else if (limiter < 10) then  ! approximate LP, [123]
         R = (omtl + sqrt(omtl*(1.d0+5.d0*lam))) / (2.d0*lam+1.d-50)
      else if (limiter < 20) then  ! Bruenn, 1[123]
         R =omtl/(lam+1.d-50)
      else if (limiter < 30) then  ! Larsen's square root, 2[123]
         R = sqrt(omtl*(1.d0+3.d0*lam)) / (lam+1.d-50)
      else if (limiter < 40) then  ! Minerbo
         if (lam .gt. TwoNinths) then
            R = sqrt(omtl/3.d0) / (lam+1.d-50)
         else
            R = 1.d0/(lam+1.d-50) - sqrt(2.d0/(lam+1.d-50))
         end if
      else
         print *, "Unknown limiter ", limiter
         stop
      endif
      
      if (R .lt. 1.d-6) then
         alpha = 0.25d0
      else if ( R .gt. 300.d0) then
         alpha = 0.5d0
      else
         cr = cosh(R)
         alpha = cr*log(cr) / (2.d0*R*sinh(R))
      end if
    end function FLDalpha

end module fluxlimiter_module

subroutine ca_initfluxlimiter(limiter, closure)  bind(C, name="ca_initfluxlimiter")
  use fluxlimiter_module, only : init_fluxlimiter_module
  implicit none
  integer, intent(in) :: limiter, closure
  call init_fluxlimiter_module(limiter, closure)
end subroutine ca_initfluxlimiter

