
! :::
! ::: --------------------------------------------------------------------
! :::
subroutine ca_initmag(level, time, lo, hi, &
                      nbx, mag_x, bx_lo, bx_hi, &
                      nby, mag_y, by_lo, by_hi, &
                      nbz, mag_z, bz_lo, bz_hi, &
                      delta, xlo, xhi)

  use probdata_module
  use prob_params_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer :: level, nbx, nby, nbz
  integer :: lo(3), hi(3)
  integer :: bx_lo(3), bx_hi(3)
  integer :: by_lo(3), by_hi(3)
  integer :: bz_lo(3), bz_hi(3)
  real(rt) :: xlo(3), xhi(3), time, delta(3)

  real(rt) :: mag_x(bx_lo(1):bx_hi(1), bx_lo(2):bx_hi(2), bx_lo(3):bx_hi(3), nbx)
  real(rt) :: mag_y(by_lo(1):by_hi(1), by_lo(2):by_hi(2), by_lo(3):by_hi(3), nby)
  real(rt) :: mag_z(bz_lo(1):bz_hi(1), bz_lo(2):bz_hi(2), bz_lo(3):bz_hi(3), nbz)

  real(rt) :: xcen, ycen, zcen
  integer  :: i, j, k

  ! This call does nothing by default; it should be copied to a problem directory
  ! and overwritten for the problem setup of interest.

  ! Note that the B fields are nodal, so you will loop from lo to hi+1
  ! in the field direction.

end subroutine ca_initmag

