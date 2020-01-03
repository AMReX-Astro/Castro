! reference: SANDIA REPORT SAND2005-6988
!            Advances in Radiation Modeling in ALEGRA: ......
!            T. A. Brunner, T. Mehlhorn, R. McClarren, & C. J. Kurecka

module blackbody_module

  use fundamental_constants_module, only: a_rad, k_B, hplanck
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), parameter :: pi = 3.1415926535897932384626e0_rt
  real(rt), parameter :: bk_const = a_rad * 15.e0_rt/pi**4
  real(rt), parameter :: magic_const = pi**4/15.e0_rt
  integer, parameter :: Ncoefs = 8
  real(rt), dimension(0:Ncoefs) :: coefs = &
                                           [1.e0_rt/60.e0_rt, &
                                           -1.e0_rt/5040.e0_rt, &
                                            1.e0_rt/272160.e0_rt, &
                                           -1.e0_rt/13305600.e0_rt, &
                                            1.e0_rt/622702080.e0_rt, &
                                           -691.e0_rt/19615115520000.e0_rt, &
                                            1.e0_rt/1270312243200.e0_rt, &
                                           -3617.e0_rt/202741834014720000.e0_rt, &
                                            43867.e0_rt/107290978560589824000.e0_rt]
  real(rt), parameter :: tol = 1.0e-10_rt

  real(rt), parameter :: xmagic = 2.061981e0_rt
  real(rt), parameter :: xsmall = 1.e-5_rt
  real(rt), parameter :: xlarge = 100.e0_rt

  private :: pi, bk_const, magic_const, Ncoefs, coefs, tol, xmagic, xsmall, xlarge

contains

  subroutine BdBdTIndefInteg(T, nu, B, dBdT)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: T, nu
    real(rt), intent(out) :: B, dBdT
    real(rt) :: x, integ, part

    x = hplanck*nu/(k_B*T)

    if (x .gt. xlarge) then
       B = a_rad * T**4
       dBdT = 4.e0_rt * a_rad * T**3
    else if ( x .lt. xsmall) then
       B = 0.e0_rt
       dBdT = 0.e0_rt
    else
       if (x .gt. xmagic) then
          integ = integlarge(x)
       else
          integ = integsmall(x)
       end if
       part = x**4/(exp(x) - 1.e0_rt)
       B = bk_const*T**4 * integ
       dBdT = bk_const*T**3 * (4.*integ - part)
    end if

  end subroutine BdBdTIndefInteg

  function BIndefInteg(T, nu)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: T, nu
    real(rt) :: BIndefInteg
    real(rt) :: x, integ

    x = hplanck*nu/(k_B*T)

    if (x .gt. xlarge) then
       BIndefInteg = a_rad * T**4
    else if ( x .lt. xsmall) then
       BIndefInteg = 0.e0_rt
    else
       if (x .gt. xmagic) then
          integ = integlarge(x)
       else
          integ = integsmall(x)
       end if
       BIndefInteg = bk_const*T**4 * integ
    end if

  end function BIndefInteg


  function BGroup(T, nu0, nu1)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: T, nu0, nu1
    real(rt) :: BGroup

    BGroup = BIndefInteg(T,nu1) - BIndefInteg(T,nu0)

  end function BGroup


  function Li(n, z)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: n
    real(rt), intent(in) :: z
    real(rt) :: Li, t
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

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: x
    real(rt) :: integlarge, z

    z = exp(-x)
    integlarge = magic_const &
                 - (x**3 * Li(1,z) + 3.e0_rt * x**2 * Li(2,z) &
                 + 6.e0_rt* x * Li(3,z) + 6.e0_rt * Li(4,z))

  end function integlarge


  function integsmall(x)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: x
    real(rt) :: integsmall, t
    integer :: i
    real(rt) :: x2, x3, xfoo

    x2 = x**2
    x3 = x**3
    integsmall = x3/3.e0_rt - x2**2/8.e0_rt
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
