
! This module stores physical constants and radiation group information
! used for multigroup photon and neutrino radiation diffusion.
! These parameters are initialized in ca_initgroups? to match the
! values used in the C++ radiation code.

module rad_params_module

  ! radiation energy group information

  use amrex_fort_module, only: rt => amrex_real
  use state_indices_module, only : ngroups

  implicit none

  integer, allocatable, save :: current_group, ng0, ng1
  real(rt), save, allocatable :: nugroup(:), dnugroup(:), xnu(:), dlognu(:), &
                                 erg2rhoYe(:), lognugroup(:)

  ! physical constants used for radiation
  real(rt), allocatable :: pi, clight, hplanck, kboltz, stefbol, arad, avogadro
  real(rt), allocatable :: Hz2MeV, mev2erg, tiny

  ! In our current solvers, E is stored in rad.  (In the past, J was stored.)
  ! So we use the following conversion factors to make sure the right variables are used
  real(rt), allocatable :: radtoE  !, radtoJ, Etorad, radfluxtoF
  real(rt), allocatable :: etafactor

  ! (yes, I know pi isn't a physical constant)
  ! (stefbol is derived from the other constants)
  ! (tiny a generic very small quantity without units, currently 1.e-50_rt)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: current_group, ng0, ng1
  attributes(managed) :: nugroup, dnugroup, xnu, dlognu, erg2rhoYe, lognugroup
  attributes(managed) :: pi, clight, hplanck, kboltz, stefbol, arad, avogadro
  attributes(managed) :: Hz2MeV, mev2erg, tiny
  attributes(managed) :: radtoE
  attributes(managed) :: etafactor
#endif

contains

  function get_ispec(g) result(ispec)

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: g
    integer :: ispec

    if (ng0 .eq. 0) then  ! photon
       ispec = 0
    else if (g < ng0) then
       ispec = 0
    else if (g < ng0+ng1) then
       ispec = 1
    else
       ispec = 2
    end if
  end function get_ispec

end module rad_params_module
