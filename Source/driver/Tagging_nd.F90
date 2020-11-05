module tagging_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), allocatable ::    denerr,   dengrad, dengrad_rel
  real(rt), allocatable ::    enterr,   entgrad, entgrad_rel
  real(rt), allocatable ::    velerr,   velgrad, velgrad_rel
  real(rt), allocatable ::   temperr,  tempgrad, tempgrad_rel
  real(rt), allocatable ::  presserr, pressgrad, pressgrad_rel
  real(rt), allocatable ::    raderr,   radgrad, radgrad_rel
  real(rt), allocatable ::   enucerr

  integer, allocatable ::  max_denerr_lev,   max_dengrad_lev, max_dengrad_rel_lev
  integer, allocatable ::  max_velerr_lev,   max_velgrad_lev, max_velgrad_rel_lev
  integer, allocatable ::  max_temperr_lev,  max_tempgrad_lev, max_tempgrad_rel_lev
  integer, allocatable ::  max_presserr_lev, max_pressgrad_lev, max_pressgrad_rel_lev
  integer, allocatable ::  max_raderr_lev,   max_radgrad_lev, max_radgrad_rel_lev
  integer, allocatable ::  max_enucerr_lev

  ! limit the zone size based on how much the burning can change the
  ! internal energy of a zone. The zone size on the finest level must
  ! be smaller than dxnuc * c_s * (e/ \dot{e}) where c_s is the sound
  ! speed.  This ensures that the sound-crossing time is smaller than
  ! the nuclear energy injection timescale.
  real(rt), allocatable :: dxnuc_min

  ! Disable limiting based on dxnuc above this threshold. This allows
  !  zones that have already ignited or are about to ignite to be
  !  de-refined.
  real(rt), allocatable :: dxnuc_max

  ! Disable limiting based on dxnuc above this AMR level.
  integer, allocatable :: max_dxnuc_lev

  public

contains

  subroutine get_denerr_params(denerr_out, max_denerr_lev_out, &
                               dengrad_out, max_dengrad_lev_out, &
                               dengrad_rel_out, max_dengrad_rel_lev_out) &
                               bind(c, name='get_denerr_params')

    implicit none

    real(rt), intent(out) :: denerr_out, dengrad_out, dengrad_rel_out
    integer, intent(out) :: max_denerr_lev_out, max_dengrad_lev_out, max_dengrad_rel_lev_out

    denerr_out = denerr
    dengrad_out = dengrad
    dengrad_rel_out = dengrad_rel

    max_denerr_lev_out = max_denerr_lev
    max_dengrad_lev_out = max_dengrad_lev
    max_dengrad_rel_lev_out = max_dengrad_rel_lev

  end subroutine get_denerr_params

  subroutine get_velerr_params(velerr_out, max_velerr_lev_out, &
                               velgrad_out, max_velgrad_lev_out, &
                               velgrad_rel_out, max_velgrad_rel_lev_out) &
                               bind(c, name='get_velerr_params')

    implicit none

    real(rt), intent(out) :: velerr_out, velgrad_out, velgrad_rel_out
    integer, intent(out) :: max_velerr_lev_out, max_velgrad_lev_out, max_velgrad_rel_lev_out

    velerr_out = velerr
    velgrad_out = velgrad
    velgrad_rel_out = velgrad_rel

    max_velerr_lev_out = max_velerr_lev
    max_velgrad_lev_out = max_velgrad_lev
    max_velgrad_rel_lev_out = max_velgrad_rel_lev

  end subroutine get_velerr_params

  subroutine get_temperr_params(temperr_out, max_temperr_lev_out, &
                                tempgrad_out, max_tempgrad_lev_out, &
                                tempgrad_rel_out, max_tempgrad_rel_lev_out) &
                                bind(c, name='get_temperr_params')

    implicit none

    real(rt), intent(out) :: temperr_out, tempgrad_out, tempgrad_rel_out
    integer, intent(out) :: max_temperr_lev_out, max_tempgrad_lev_out, max_tempgrad_rel_lev_out

    temperr_out = temperr
    tempgrad_out = tempgrad
    tempgrad_rel_out = tempgrad_rel

    max_temperr_lev_out = max_temperr_lev
    max_tempgrad_lev_out = max_tempgrad_lev
    max_tempgrad_rel_lev_out = max_tempgrad_rel_lev

  end subroutine get_temperr_params

  subroutine get_presserr_params(presserr_out, max_presserr_lev_out, &
                                 pressgrad_out, max_pressgrad_lev_out, &
                                 pressgrad_rel_out, max_pressgrad_rel_lev_out) &
                                 bind(c, name='get_presserr_params')

    implicit none

    real(rt), intent(out) :: presserr_out, pressgrad_out, pressgrad_rel_out
    integer, intent(out) :: max_presserr_lev_out, max_pressgrad_lev_out, max_pressgrad_rel_lev_out

    presserr_out = presserr
    pressgrad_out = pressgrad
    pressgrad_rel_out = pressgrad_rel

    max_presserr_lev_out = max_presserr_lev
    max_pressgrad_lev_out = max_pressgrad_lev
    max_pressgrad_rel_lev_out = max_pressgrad_rel_lev

  end subroutine get_presserr_params

  subroutine get_raderr_params(raderr_out, max_raderr_lev_out, &
                               radgrad_out, max_radgrad_lev_out, &
                               radgrad_rel_out, max_radgrad_rel_lev_out) &
                               bind(c, name='get_raderr_params')

    implicit none

    real(rt), intent(out) :: raderr_out, radgrad_out, radgrad_rel_out
    integer, intent(out) :: max_raderr_lev_out, max_radgrad_lev_out, max_radgrad_rel_lev_out

    raderr_out = raderr
    radgrad_out = radgrad
    radgrad_rel_out = radgrad_rel

    max_raderr_lev_out = max_raderr_lev
    max_radgrad_lev_out = max_radgrad_lev
    max_radgrad_rel_lev_out = max_radgrad_rel_lev

  end subroutine get_raderr_params

  subroutine get_dxnuc_params(dxnuc_min_out, dxnuc_max_out, max_dxnuc_lev_out) &
                              bind(c, name='get_dxnuc_params')


    implicit none

    real(rt), intent(out) :: dxnuc_min_out, dxnuc_max_out
    integer, intent(out) :: max_dxnuc_lev_out

    dxnuc_min_out = dxnuc_min
    dxnuc_max_out = dxnuc_max

    max_dxnuc_lev_out = max_dxnuc_lev

  end subroutine get_dxnuc_params

  subroutine get_enuc_params(enucerr_out, max_enucerr_lev_out) &
                             bind(c, name='get_enuc_params')


    implicit none

    real(rt), intent(out) :: enucerr_out
    integer, intent(out) :: max_enucerr_lev_out

    enucerr_out = enucerr

    max_enucerr_lev_out = max_enucerr_lev

  end subroutine get_enuc_params

end module tagging_module
