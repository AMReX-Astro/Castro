subroutine fdustcollapse1d(lo, hi, p, plo, phi, nc_p, nbins, dens, &
     imask, mask_size, r1, dens_comp, cnt) bind(C, name='fdustcollapse1d')

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  integer, intent(in), value :: nbins
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),1:nc_p)
  real(rt), intent(inout) :: dens(0:nbins-1)
  integer, intent(inout) :: imask(0:mask_size-1)
  integer, intent(in), value :: mask_size, r1, dens_comp
  integer, intent(inout) :: cnt

  integer :: i, j, k, index

  ! loop over all of the zones in the patch.  Here, we convert
  ! the cell-centered indices at the current level into the
  ! corresponding RANGE on the finest level, and test if we've
  ! stored data in any of those locations.  If we haven't then
  ! we store this level's data and mark that range as filled.
  j = lo(2)
  k = lo(3)

  do i = lo(1), hi(1)

     if ( any(imask(i*r1:(i+1)*r1-1) .eq. 1) ) then

        dens(cnt)   = p(i,j,k,dens_comp)

        imask(i*r1:(i+1)*r1-1) = 0

        cnt = cnt + 1

     end if

  enddo

end subroutine fdustcollapse1d
