! Process a 1-d sedov problem to produce rho, u, and p as a
! function of r, for comparison to the analytic solution.

subroutine fextract1d(lo, hi, p, plo, phi, nc_p, nbins, dens_bin, &
    vel_bin, pres_bin, e_bin, imask, mask_size, r1,&
    dens_comp, xmom_comp, pres_comp, rhoe_comp) bind(C, name='fextract1d')

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  integer, intent(in), value :: nbins
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),1:nc_p)
  real(rt), intent(inout) :: dens_bin(0:nbins-1)
  real(rt), intent(inout) :: vel_bin(0:nbins-1)
  real(rt), intent(inout) :: pres_bin(0:nbins-1)
  real(rt), intent(inout) :: e_bin(0:nbins-1)
  integer, intent(inout) :: imask(0:mask_size-1)
  integer, intent(in), value :: mask_size, r1, dens_comp, xmom_comp, pres_comp, rhoe_comp

  integer :: ii, index

  ! loop over all of the zones in the patch.  Here, we convert
  ! the cell-centered indices at the current level into the
  ! corresponding RANGE on the finest level, and test if we've
  ! stored data in any of those locations.  If we haven't then
  ! we store this level's data and mark that range as filled.
  do ii = plo(1), phi(1)

     if ( any(imask(ii*r1:(ii+1)*r1-1) .eq. 1) ) then

        index = ii * r1

        dens_bin(index:index+(r1-1)) = p(ii,1,1,dens_comp)

        vel_bin(index:index+(r1-1)) = &
             abs(p(ii,1,1,xmom_comp)) / p(ii,1,1,dens_comp)

        pres_bin(index:index+(r1-1)) = p(ii,1,1,pres_comp)

        e_bin(index:index+(r1-1)) = &
             p(ii,1,1,rhoe_comp)  / p(ii,1,1,dens_comp)

        imask(ii*r1:(ii+1)*r1-1) = .false.

     end if

  enddo

end subroutine fextract1d
