! Compute the analytic solution for the multigroup radiating sphere
! problem.  We follow the problem outline from Swesty and Myra (2009).
!
! The output from this is the radiation energy density as a function
! of energy group at a specified distance (r_obs) from the radiating
! sphere, at a time of t_obs.

module constants_module

  ! fundamental constants
  double precision, parameter :: h_planck =  6.62606957d-27  ! erg s
  double precision, parameter :: k_B = 1.3806488d-16  ! erg / K
  double precision, parameter :: c_light = 2.99792458d10  ! cm / s
  double precision, parameter :: ev2erg = 1.602176487d-12 ! erg/eV
  double precision, parameter :: pi = 3.14159265358979323846d0

  ! problem parameters

  ! we will specify the opacity as kappa = kapp_0 (nu_0 / nu)**3
  ! where nu_0 is a reference frequency
  ! 
  ! Swesty and Myra (2009) only say that the "opacity of the material
  ! is proportional to 1/E**3", but they don't give the proportionality
  ! constant.  The output depends critically on the constant.  Eric
  ! says that they used kappa = 1.e13 * (1.5e-6 MeV / E)**3
  !
  ! Converting 1.5e-6 MeV to a frequency gives nu_0 = 3.63e14 Hz
  double precision, parameter :: nu_0 = 3.6e14  ! ref. freq (Hz)
  double precision, parameter :: kappa_0 = 1.e13  ! scattering opacity (1/cm)

  ! specify the energy groups.  Ultimately we will work in terms of 
  ! frequency, but we set the limits in terms of energy.  
  integer, parameter :: ngroups_derive = 60          ! number of groups to compute
  double precision, parameter :: E_min = 0.5d0*ev2erg ! min group energy (erg)
  double precision, parameter :: E_max = 306d3*ev2erg ! max group energy (erg)

  ! instead of computing the group structure based on the above, read the group
  ! structure in from group_structure.dat
  logical, parameter :: use_group_file = .true.


  ! geometry parameters
  double precision, parameter :: R_sphere = 0.02d0 ! sphere radius (cm)
  double precision, parameter :: r_obs = 0.06d0  ! observer location (cm)

  ! physical parameters
  double precision, parameter :: T_0 = 50.d0*ev2erg/k_B ! ambient temp (K)
  double precision, parameter :: T_sphere = 1500.d0*ev2erg/k_B ! sphere temp (K
  
  ! output time (t_obs = 1.d-12 is a good value to try)
  double precision, parameter :: t_obs = 1.d-12  ! observation time (s)


end module constants_module

function safe_print(x) result (sx)
  implicit none
  double precision :: x, sx
  sx = x
  if (x < 1.d-99) sx = 0.d0

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

  double precision, intent(in) :: nu, T
  double precision :: B

  B = (8.d0*pi*h_planck*nu**3/c_light**3)/(exp(h_planck*nu/(k_B*T)) - 1.d0)


  return
end function planck


function F_radsphere(r,t,nu) result (F)

  ! the function F(r,t) as defined in Swesty and Myra
  !
  ! r      = position
  ! t      = time (s)
  ! nu     = frequency (Hz)

  use constants_module

  implicit none

  double precision, intent(in) :: r, t, nu

  double precision :: E, kappa, erfc_term1, erfc_term2, F

  double precision :: sferfc, sfexp, arg1

  ! to guard against over/underflows, we define some 'safe' functions
  ! (they come from Eric)
  sferfc(arg1) = erfc(max(-26.0d0,min(26.0d0,arg1)))
  sfexp(arg1) = exp(max(-650.0d0,min(650.0d0,arg1)))


  kappa = kappa_0*(nu_0/nu)**3


  erfc_term1 = sferfc(sqrt(3.d0*kappa/(4.d0*c_light*t))*(r - R_sphere) - &
                      sqrt(c_light*t*kappa))

  erfc_term2 = sferfc(sqrt(3.d0*kappa/(4.d0*c_light*t))*(r - R_sphere) + &
                    sqrt(c_light*t*kappa))


  F = 0.5*(sfexp(-sqrt(3.d0)*kappa*(r - R_sphere))*erfc_term1 + &
           sfexp( sqrt(3.d0)*kappa*(r - R_sphere))*erfc_term2)

  return

end function F_radsphere


program analytic

  ! compute the analytic solution to the radiating sphere, as given
  ! by Swesty and Myra (2009), Eq. 76.

  use constants_module

  implicit none

  integer :: ngroups
  double precision, allocatable :: nu_groups(:), dnu_groups(:), xnu(:)
  double precision :: alpha
  double precision :: E, F
  double precision :: F_radsphere, planck
  double precision :: safe_print
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
     alpha = (E_max/E_min)**(1.d0/(ngroups))

     do n = 2, ngroups+1
        xnu(n) = alpha*xnu(n-1)
     enddo

     do n = 1, ngroups
        nu_groups(n) = 0.5d0*(xnu(n) + xnu(n+1))
     enddo

     do n = 1, ngroups
        if (n == 1) then
           dnu_groups(n) = 0.5d0*(nu_groups(2) - nu_groups(1))

        else if (n < ngroups) then
           dnu_groups(n) = 0.5d0*(nu_groups(n+1) - nu_groups(n-1))

        else
           dnu_groups(n) = 0.5d0*(nu_groups(n) - nu_groups(n-1))
        endif

     enddo

  endif

  ! debug -- print out the group structure
  !do n = 1, ngroups
  !   print *, nu_groups(n), dnu_groups(n)
  !enddo

  ! compute the analytic radiating sphere solution for each point
1000 format(1x, i3, 4(1x,g20.10))

  print *, "# group, group center (Hz), analytic sol. (erg/cm^3/Hz), &
       &analytic sol. * group width (erg/cm^3), ambient BB spectrum (erg/cm^3/Hz)"
  do n = 1, ngroups
     F = F_radsphere(r_obs,t_obs,nu_groups(n))
     E = planck(nu_groups(n),T_0) + &
          (R_sphere/r_obs)*(planck(nu_groups(n),T_sphere) - &
                            planck(nu_groups(n),T_0))*F
     write(*,1000) n, nu_groups(n), &
          safe_print(E*dnu_groups(n)), safe_print(E)

  enddo     

end program analytic
  
