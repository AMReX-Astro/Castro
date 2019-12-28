! Compute the Planck spectrum for a multigroup problem, to serve as
! initial conditions.

module constants_module

  use amrex_fort_module, only : rt => amrex_real

  ! fundamental constants
  real(rt), parameter :: h_planck = 6.62606896e-27_rt  ! erg s
  real(rt), parameter :: k_B = 1.3806504e-16_rt  ! erg / K
  real(rt), parameter :: c_light = 2.99792458e10_rt  ! cm / s
  real(rt), parameter :: ev2erg = 1.602176487e-12_rt ! erg/eV
  real(rt), parameter :: pi = 3.14159265358979323846e0_rt
  real(rt), parameter :: sigma_SB = 5.670400e-5_rt  ! erg/s/cm^2/K^4
  real(rt), parameter :: a_rad = 4.0e0_rt*sigma_SB/c_light

  ! problem parameters

  ! specify the energy groups.  Ultimately we will work in terms of
  ! frequency, but we set the limits in terms of energy.

  ! note: E_min and E_max MUST match the parameters used in CASTRO
  ! to define the group structure.
  integer, parameter :: ngroups_derive = 128          ! number of groups to compute
  real(rt), parameter :: E_min = 0.5e0_rt*ev2erg ! min group energy (erg)
  real(rt), parameter :: E_max = 306e3_rt*ev2erg ! max group energy (erg)

  ! instead of computing the group structure based on the above, read the group
  ! structure in from group_structure.dat
  logical, parameter :: use_group_file = .false.

  ! physical parameters
  real(rt), parameter :: E_rad = 1.e12_rt ! initial radiation energy density
  real(rt), parameter :: T_0 = (E_rad/a_rad)**0.25e0_rt


end module constants_module

function safe_print(x) result (sx)
  implicit none
  real(rt) :: x, sx
  sx = x
  if (x < 1.e-99_rt) sx = 0.e0_rt

  return
end function safe_print

function planck(nu,T) result (B)

  ! the Planck function for a Blackbody (actually, energy density,
  ! B = (4 pi / c) I, where I is the normal Planck function
  !
  ! nu = frequency (Hz)
  ! T  = temperature (K)

  ! see Swesty and Myra (2009), eq. 23, but note that we are working
  ! in Hz here, so we have units of erg / cm^3 / Hz, where they
  ! have units of erg / cm^3 / MeV.  As a result, we have one
  ! less factor of h_planck.
  use constants_module
  implicit none

  real(rt), intent(in) :: nu, T
  real(rt) :: B

  B = (8.e0_rt*pi*h_planck*nu**3/c_light**3)/(exp(h_planck*nu/(k_B*T)) - 1.e0_rt)


  return
end function planck


program analytic

  ! compute the analytic solution to the radiating sphere, as given
  ! by Swesty and Myra (2009), Eq. 76.

  use constants_module

  implicit none

  integer :: ngroups
  real(rt), allocatable :: nu_groups(:), dnu_groups(:), xnu(:)
  real(rt) :: alpha
  real(rt) :: E, F
  real(rt) :: F_radsphere, planck
  real(rt) :: P_simpsonrule
  real(rt) :: safe_print
  integer :: n

  character(len=256) :: header_line
  integer :: ipos

  if (use_group_file) then
     ! open the file
     open(unit=10,file="group_structure.dat")

     ! read in the number of groups
     read(10,fmt='(a256)') header_line
     ipos = index(header_line, "=") + 1
     read (header_line(ipos:),*) ngroups

     allocate (nu_groups(ngroups), dnu_groups(ngroups))

     ! skip the next header line
     read(10,*) header_line

     ! read in the group centers and weights (widths)
     do n = 1, ngroups
        read(10,*) nu_groups(n), dnu_groups(n)
     enddo

  else
     ngroups = ngroups_derive
     allocate(nu_groups(ngroups),dnu_groups(ngroups),xnu(ngroups+1))

     ! try to mimic the way the groups are setup in RadMultiGroup.cpp

     ! do geometrically spaced group BOUNDARIES
     xnu(1) = E_min/h_planck
     alpha = (E_max/E_min)**(1.e0_rt/(ngroups))

     do n = 2, ngroups+1
        xnu(n) = alpha*xnu(n-1)
     enddo

     do n = 1, ngroups
        nu_groups(n) = 0.5e0_rt*(xnu(n) + xnu(n+1))
     enddo

     do n = 1, ngroups
        if (n == 1) then
           dnu_groups(n) = 0.5e0_rt*(nu_groups(2) - nu_groups(1))

        else if (n < ngroups) then
           dnu_groups(n) = 0.5e0_rt*(nu_groups(n+1) - nu_groups(n-1))

        else
           dnu_groups(n) = 0.5e0_rt*(nu_groups(n) - nu_groups(n-1))
        endif

     enddo

  endif


  ! compute the analytic radiating sphere solution for each point
1000 format(1x, i3, 4(1x,g20.10))

  print *, "# group, group center (Hz), ambient BB spectrum x dnu (erg/cm^3)"
  do n = 1, ngroups
     P_simpsonrule = (dnu_groups(n)/6.0e0_rt)* &
          (     planck(xnu(n),T_0) + &
          4.0e0_rt*planck(nu_groups(n),T_0) + &
          planck(xnu(n+1),T_0))

     write(*,1000) n, nu_groups(n), dnu_groups(n), &
          safe_print(planck(nu_groups(n),T_0)*dnu_groups(n))

  enddo

end program analytic
