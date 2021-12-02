
! This module stores physical constants and radiation group information
! used for multigroup photon and neutrino radiation diffusion.
! These parameters are initialized in ca_initgroups? to match the
! values used in the C++ radiation code.

module rad_params_module

  ! radiation energy group information

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, parameter :: ngroups = NGROUPS

  integer, allocatable, save :: current_group, ng0, ng1
  real(rt), save, allocatable :: nugroup(:), dnugroup(:), xnu(:), dlognu(:), &
                                 erg2rhoYe(:), lognugroup(:)

end module rad_params_module
