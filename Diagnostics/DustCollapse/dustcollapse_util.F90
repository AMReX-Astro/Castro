! Process a group of n-d plotfiles from the dustcollapse problem
! and output the position of the interface as a function of time.
!
! The initial dense sphere is assumed to be centered a r = 0 (x = 0).
! We take as default that it is centered vertically at y = 0, but
! this can be overridden with --yctr.
!
! The --profile option will write out the average density vs. radius
! profile to a file (plotfile name + '.profile')

subroutine fdustcollapse1d(lo, hi, p, plo, phi, nc_p, nbins, dens, &
     imask, mask_size, r1, dens_comp, cnt) bind(C, name='fdustcollapse1d')

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  integer, intent(in), value :: nbins
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),0:nc_p-1)
  real(rt), intent(inout) :: dens(0:nbins-1)
  integer, intent(inout) :: imask(0:mask_size-1)
  integer, intent(in), value :: mask_size, r1, dens_comp
  integer, intent(inout) :: cnt

  integer :: i, j, k

  ! loop over all of the zones in the patch.  Here, we convert
  ! the cell-centered indices at the current level into the
  ! corresponding RANGE on the finest level, and test if we've
  ! stored data in any of those locations.  If we haven't then
  ! we store this level's data and mark that range as filled.
  j = lo(2)
  k = lo(3)

  do i = lo(1), hi(1)

     if (any(imask(i*r1:(i+1)*r1-1) .eq. 1)) then

        dens(cnt) = p(i,j,k,dens_comp)

        imask(i*r1:(i+1)*r1-1) = 0

        cnt = cnt + 1

     end if

  enddo

end subroutine fdustcollapse1d

subroutine fdustcollapse2d(lo, hi, p, plo, phi, nc_p, nbins, dens, &
     volcount, imask, mask_size, r1, dx, dx_fine, yctr, dens_comp) &
     bind(C, name='fdustcollapse2d')

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  integer, intent(in), value :: nbins
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),0:nc_p-1)
  real(rt), intent(inout) :: dens(0:nbins-1)
  real(rt), intent(inout) :: volcount(0:nbins-1)
  integer, intent(inout) :: imask(0:mask_size-1,0:mask_size-1)
  integer, intent(in), value :: mask_size, r1, dens_comp
  real(rt), intent(in) :: dx(3)
  real(rt), intent(in), value :: yctr, dx_fine

  integer :: i, j, k, index
  real(rt) :: xx, xl, xr, yy, yl, yr, r_zone, vol

  ! loop over all of the zones in the patch.  Here, we convert
  ! the cell-centered indices at the current level into the
  ! corresponding RANGE on the finest level, and test if we've
  ! stored data in any of those locations.  If we haven't then
  ! we store this level's data and mark that range as filled.
  k = lo(3)

  do j = lo(2), hi(2)
     yy = (dble(j) + 0.5d0)*dx(2)
     yl = (dble(j))*dx(2)
     yr = (dble(j) + 1.d0)*dx(2)

     do i = lo(1), hi(1)
        xx = (dble(i) + 0.5d0)*dx(1)
        xl = (dble(i))*dx(1)
        xr = (dble(i) + 1.d0)*dx(1)

        if ( any(imask(i*r1:(i+1)*r1-1, &
             j*r1:(j+1)*r1-1) .eq. 1) ) then

           r_zone = sqrt((xx)**2 + (yy-yctr)**2)

           index = r_zone/dx_fine

           vol = (xr**2 - xl**2)*(yr - yl)

           ! weight the zone's data by its size
           dens(index) = dens(index) + p(i,j,k,dens_comp) * vol

           volcount(index) = volcount(index) + vol

           imask(i*r1:(i+1)*r1-1, &
                j*r1:(j+1)*r1-1) = 0

        end if

     enddo
  enddo

end subroutine fdustcollapse2d

subroutine fdustcollapse3d(lo, hi, p, plo, phi, nc_p, nbins, dens, &
     ncount, imask, mask_size, r1, dx, dx_fine, xctr, yctr, zctr, dens_comp) &
     bind(C, name='fdustcollapse3d')

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  integer, intent(in), value :: nbins
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),0:nc_p-1)
  real(rt), intent(inout) :: dens(0:nbins-1)
  real(rt), intent(inout) :: ncount(0:nbins-1)
  integer, intent(inout) :: imask(0:mask_size-1,0:mask_size-1,0:mask_size-1)
  integer, intent(in), value :: mask_size, r1, dens_comp
  real(rt), intent(in) :: dx(3)
  real(rt), intent(in), value :: xctr, yctr, zctr, dx_fine

  integer :: i, j, k, index
  real(rt) :: xx, yy, zz, r_zone

  ! loop over all of the zones in the patch.  Here, we convert
  ! the cell-centered indices at the current level into the
  ! corresponding RANGE on the finest level, and test if we've
  ! stored data in any of those locations.  If we haven't then
  ! we store this level's data and mark that range as filled.
  do k = lo(3), hi(3)
     zz = (dble(k) + 0.5d0)*dx(3)

     do j = lo(2), hi(2)
        yy = (dble(j) + 0.5d0)*dx(2)

        do i = lo(1), hi(1)
           xx = (dble(i) + 0.5d0)*dx(1)

           if ( any(imask(i*r1:(i+1)*r1-1, &
                j*r1:(j+1)*r1-1, &
                k*r1:(k+1)*r1-1) .eq. 1) ) then

              r_zone = sqrt((xx-xctr)**2 + (yy-yctr)**2 + (zz-zctr)**2)

              index = r_zone/dx_fine

              ! weight the zone's data by its size
              dens(index) = dens(index) + p(i,j,k,dens_comp)*r1**3

              ncount(index) = ncount(index) + r1**3

              imask(i*r1:(i+1)*r1-1, &
                   j*r1:(j+1)*r1-1, &
                   k*r1:(k+1)*r1-1) = 0

           end if

        enddo
     enddo
  enddo


end subroutine fdustcollapse3d
