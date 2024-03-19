! write out the group dependent boundary conditions
!
! This routine takes as input the file "group_structures.dat" which is
! output by RadMultiGroup.cpp.  That file defines the group centers
! and weights (widths).
!
! NOTE: in CASTRO, for multigroup, the radiation energy is stored with
! units of erg / cm^3, so that the total radiation energy can be found
! simply by summing over the groups.  This means that we normalize by
! multiplying the Planck function in each group by the group's width.

module constants_module

  ! fundamental constants
  use amrex_fort_module, only : rt => amrex_real
  real(rt)        , parameter :: h_planck = 6.62606957e-27_rt   ! erg s
  real(rt)        , parameter :: k_B = 1.3806488e-16_rt  ! erg / K
  real(rt)        , parameter :: c_light = 2.99792458e10_rt  ! cm / s
  real(rt)        , parameter :: ev2erg = 1.602176487e-12_rt ! erg/eV
  real(rt)        , parameter :: pi = 3.14159265358979323846e0_rt

  ! physical parameters
  real(rt)        , parameter :: T_sphere = 1500.e0_rt*ev2erg/k_B ! sphere temp (K


end module constants_module


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
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , intent(in) :: nu, T
  real(rt)         :: B

  B = (8.e0_rt*pi*h_planck*nu**3/c_light**3)/(exp(h_planck*nu/(k_B*T)) - 1.e0_rt)

  return
end function planck


program bc

  ! compute the values of the blackbody function for the different
  ! energy groups.

  use constants_module

  implicit none

  integer :: ngroups
  real(rt)        , allocatable :: nu_groups(:), dnu_groups(:)
  real(rt)         :: planck

  integer :: n
  character(len=256) :: header_line
  integer :: ipos

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

  do n = 1, ngroups
     print *, n, nu_groups(n), planck(nu_groups(n), T_sphere)*dnu_groups(n)
  enddo

  print *, "radiation.lo_bcval0 = ", &
       ((planck(nu_groups(n), T_sphere)*dnu_groups(n)), n=1,ngroups)

end program bc

