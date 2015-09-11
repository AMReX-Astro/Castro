! reference: SANDIA REPORT SAND2005-6988
!            Advances in Radiation Modeling in ALEGRA: ......
!            T. A. Brunner, T. Mehlhorn, R. McClarren, & C. J. Kurecka

module blackbody_module

  use fundamental_constants_module, only : a_rad, k_B, hplanck

  implicit none
  
  double precision, parameter :: pi = 3.1415926535897932384626d0
  double precision :: bk_const = a_rad * 15.d0/pi**4
  !                            = 8 * pi * k**4 / (c**3 * h**3)
  double precision :: magic_const = pi**4/15.d0
  double precision, parameter :: xmagic = 2.061981d0
  integer, parameter :: Ncoefs = 8
  double precision, dimension(0:Ncoefs) :: coefs =  &
       (/ 1.d0/60.d0, -1.d0/5040.d0, 1.d0/272160.d0, -1.d0/13305600.d0, &
       1.d0/622702080.d0, -691.d0/19615115520000.d0, 1.d0/1270312243200.d0,  &
       -3617.d0/202741834014720000.d0, 43867.d0/107290978560589824000.d0 /)
  double precision, parameter :: tol = 1.0d-10
  
  contains

    subroutine BdBdTIndefInteg(T, nu, B, dBdT)
      double precision, intent(in) :: T, nu
      double precision, intent(out) :: B, dBdT
      double precision :: x, integ, part
      
      x = hplanck*nu/(k_B*T)

      if (x .gt. xmagic) then
         integ = integlarge(x)
      else if (x .lt. 1.d-12) then
         integ = 0.d0
      else
         integ = integsmall(x)
      end if
      
      B = bk_const*T**4 * integ

      if (x .gt. 100.d0) then
         part = 0.d0
      else if (x .lt. 1.d-12) then
         part = 0.d0
      else
         part = x**4/(exp(x) - 1.d0)
      end if

      dBdT = bk_const*T**3 * (4.*integ - part)
    end subroutine BdBdTIndefInteg

    function BIndefInteg(T, nu)
      double precision BIndefInteg
      double precision, intent(in) :: T, nu
      double precision :: x, integ
      
      x = hplanck*nu/(k_B*T)

      if (x .gt. xmagic) then
         integ = integlarge(x)
      else
         integ = integsmall(x)
      end if
      
      BIndefInteg = bk_const*T**4 * integ
    end function BIndefInteg

    function BGroup(T, nu0, nu1)
      double precision BGroup
      double precision, intent(in) :: T, nu0, nu1
      BGroup = BIndefInteg(T,nu1) - BIndefInteg(T,nu0)
    end function BGroup
    
    function Li(n, z)
      integer, intent(in) :: n
      double precision, intent(in) :: z
      double precision :: Li, t
      integer :: k
      integer, parameter :: kmax = 18
      Li = z
      do k = 2, kmax
         t = z**k / k**n
         Li = Li + t
         if (t/Li .lt. tol) then
            exit
         end if
      end do
    end function Li

    function integlarge(x)
      double precision, intent(in) :: x
      double precision :: integlarge, z

      z = exp(-x)
      integlarge = magic_const & 
           - (x**3 * Li(1,z) + 3.d0 * x**2 * Li(2,z) &
           + 6.d0* x * Li(3,z) + 6.d0 * Li(4,z))
    end function integlarge

    function integsmall(x)
      double precision, intent(in) :: x
      double precision :: integsmall, t
      integer :: i
      double precision :: x2, x3, xfoo

      x2 = x**2
      x3 = x**3
      integsmall = x3/3.d0 - x2**2/8.d0
      xfoo = x3
      do i = 0, Ncoefs
         xfoo = xfoo * x2
         t = coefs(i) * xfoo
         integsmall = integsmall + t
         if (abs(t/integsmall) .lt. tol) then
            exit
         end if
      end do
    end function integsmall

end module blackbody_module
