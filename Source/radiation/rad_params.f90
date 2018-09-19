
! This module stores physical constants and radiation group information
! used for multigroup photon and neutrino radiation diffusion.
! These parameters are initialized in ca_initgroups? to match the
! values used in the C++ radiation code.

module rad_params_module

  ! radiation energy group information
  use amrex_fort_module, only : rt => amrex_real
  use state_sizes_module, only : ngroups

  integer         , save :: current_group, ng0, ng1, nnuspec
  integer, save :: nradspec = 1
  real(rt)        , save, allocatable :: nugroup(:), dnugroup(:), xnu(:), dlognu(:), &
       erg2rhoYe(:), lognugroup(:)

  ! physical constants used for radiation
  real(rt)        , save :: pi, clight, hplanck, kboltz, stefbol, arad, avogadro
  real(rt)        , save :: Hz2MeV, mev2erg, tiny

  ! In our current solvers, E is stored in rad.  (In the past, J was stored.)
  ! So we use the following conversion factors to make sure the right variables are used
  real(rt)        , save :: radtoE  !, radtoJ, Etorad, radfluxtoF
  real(rt)        , save :: etafactor

  ! (yes, I know pi isn't a physical constant)
  ! (stefbol is derived from the other constants)
  ! (tiny a generic very small quantity without units, currently 1.e-50_rt)

contains

  function get_ispec(g) result(ispec)
    use amrex_fort_module, only : rt => amrex_real
    integer, intent(in) :: g
    integer ispec

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
