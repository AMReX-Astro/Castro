subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: init, namlen
  integer          :: name(namlen)
  real(rt)         :: problo(3), probhi(3)

end subroutine amrex_probinit

#ifndef GPU_COMPATIBLE_PROBLEM

subroutine ca_initdata(level,time,lo,hi,nvar, &
                       state,state_lo,state_hi, &
                       dx,xlo,xhi)
    ! This routine is called at problem setup time and is used
    ! to initialize data on each grid.
    !
    ! .. note::
    !    all arrays have one cell of ghost zones surrounding
    !    the grid interior.  Values in these cells need not
    !    be set here.
    !
    ! INPUTS/OUTPUTS:
    !
    ! level     => amr level of grid
    ! time      => time at which to init data
    ! lo,hi     => index limits of grid interior (cell centered)
    ! nvar      => number of state components.
    ! state     <= scalar array
    ! dx        => cell size
    ! xlo, xhi  => physical locations of lower left and upper
    !              right hand corner of grid.  (does not include
    !              ghost region).

  use castro_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: level, nvar
  integer :: lo(3), hi(3)
  integer :: state_lo(3), state_hi(3)
  real(rt)         :: xlo(3), xhi(3), time, dx(3)
  real(rt)         :: state(state_lo(1):state_hi(1), &
                            state_lo(2):state_hi(2), &
                            state_lo(3):state_hi(3),nvar)

  ! This call does nothing by default; it should be copied to a problem directory
  ! and overwritten for the problem setup of interest.

end subroutine ca_initdata

#else

subroutine ca_initdata(lo, hi, &
                       state, state_lo, state_hi, &
                       dx, problo) bind(C, name='ca_initdata')

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR

  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: state_lo(3), state_hi(3)
  real(rt), intent(in   ) :: dx(3), problo(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  !$gpu

  ! This call does nothing by default; it should be copied to a problem directory
  ! and overwritten for the problem setup of interest.

end subroutine ca_initdata

#endif

#ifdef RADIATION
subroutine ca_initrad(level, time, lo, hi, nrad, &
                      rad_state, rad_state_lo, rad_state_hi, &
                      delta, xlo, xhi)

  use probdata_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: level, nrad
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: rad_state_lo(3), rad_state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: rad_state(rad_state_lo(1):rad_state_hi(1), &
                                       rad_state_lo(2):rad_state_hi(2), &
                                       rad_state_lo(3):rad_state_hi(3), 0:nrad-1)

  ! This call does nothing by default; it should be copied to a problem directory
  ! and overwritten for the problem setup of interest.

end subroutine ca_initrad
#endif
