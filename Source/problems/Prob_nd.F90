subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: init, namlen
  integer          :: name(namlen)
  real(rt)         :: problo(3), probhi(3)

end subroutine amrex_probinit


! ::: -----------------------------------------------------------

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
    !		   ghost region).

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

  ! Remove this call if you're defining your own problem; it is here to
  ! ensure that you cannot run CASTRO if you haven't got your own copy of this function.

  call castro_error("Prob_nd.f90 has not been defined for this problem!")

end subroutine ca_initdata
