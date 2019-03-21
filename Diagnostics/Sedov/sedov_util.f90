! Process a sedov problem to produce rho, u, p and e as a
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
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),0:nc_p-1)
  real(rt), intent(inout) :: dens_bin(0:nbins-1)
  real(rt), intent(inout) :: vel_bin(0:nbins-1)
  real(rt), intent(inout) :: pres_bin(0:nbins-1)
  real(rt), intent(inout) :: e_bin(0:nbins-1)
  integer, intent(inout) :: imask(0:mask_size-1)
  integer, intent(in), value :: mask_size, r1, dens_comp, xmom_comp, pres_comp, rhoe_comp

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

        index = i * r1

        dens_bin(index:index+(r1-1)) = p(i,j,k,dens_comp)

        vel_bin(index:index+(r1-1)) = &
             abs(p(i,j,k,xmom_comp)) / p(i,j,k,dens_comp)

        pres_bin(index:index+(r1-1)) = p(i,j,k,pres_comp)

        e_bin(index:index+(r1-1)) = &
             p(i,j,k,rhoe_comp)  / p(i,j,k,dens_comp)

        imask(i*r1:(i+1)*r1-1) = 0

     end if

  enddo

end subroutine fextract1d

subroutine fextract2d_cyl(lo, hi, p, plo, phi, nc_p, nbins, dens_bin, &
     vel_bin, pres_bin, e_bin, ncount, imask, mask_size, r1, &
     dens_comp, xmom_comp, ymom_comp, pres_comp, rhoe_comp, dx_fine, dx, &
     xctr, yctr) &
     bind(C, name='fextract2d_cyl')

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  integer, intent(in), value :: nbins
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),0:nc_p-1)
  real(rt), intent(inout) :: dens_bin(0:nbins-1)
  real(rt), intent(inout) :: vel_bin(0:nbins-1)
  real(rt), intent(inout) :: pres_bin(0:nbins-1)
  real(rt), intent(inout) :: e_bin(0:nbins-1)
  integer, intent(inout) :: ncount(0:nbins-1)
  integer, intent(in), value :: mask_size
  integer, intent(inout) :: imask(0:mask_size-1,0:mask_size-1)
  integer, intent(in), value :: r1, dens_comp, xmom_comp, ymom_comp, pres_comp, rhoe_comp
  real(rt), intent(in), value :: dx_fine, xctr, yctr
  real(rt), intent(in) :: dx(3)

  integer :: i, j, k, index
  real(rt) :: xx, yy, r_zone

  ! loop over all of the zones in the patch.  Here, we convert
  ! the cell-centered indices at the current level into the
  ! corresponding RANGE on the finest level, and test if we've
  ! stored data in any of those locations.  If we haven't then
  ! we store this level's data and mark that range as filled.
  k = lo(3)
  do j = lo(2), hi(2)
     yy = (j + HALF)*dx(2)

     do i = lo(1), hi(1)
        xx = (i + HALF)*dx(1)

        if ( any(imask(i*r1:(i+1)*r1-1, &
             j*r1:(j+1)*r1-1) .eq. 1) ) then

           r_zone = sqrt((xx-xctr)**2 + (yy-yctr)**2)

           index = r_zone/dx_fine

           dens_bin(index) = dens_bin(index) + &
                p(i,j,k,dens_comp)*r1**2

           vel_bin(index) = vel_bin(index) + &
                (sqrt(p(i,j,k,xmom_comp)**2 + &
                p(i,j,k,ymom_comp)**2)/ &
                p(i,j,k,dens_comp))*r1**2

           pres_bin(index) = pres_bin(index) + &
                p(i,j,k,pres_comp)*r1**2

           e_bin(index) = e_bin(index) + &
                (p(i,j,k,rhoe_comp)/p(i,j,k,dens_comp))*r1**2

           ncount(index) = ncount(index) + r1**2

           imask(i*r1:(i+1)*r1-1, &
                j*r1:(j+1)*r1-1) = 0

        end if

     enddo
  enddo

end subroutine fextract2d_cyl

subroutine fextract2d_sph(lo, hi, p, plo, phi, nc_p, nbins, dens_bin, &
     vel_bin, pres_bin, e_bin, volcount, imask, mask_size, r1,&
     dens_comp, xmom_comp, ymom_comp, pres_comp, rhoe_comp, dx_fine, dx, xctr, yctr) &
     bind(C, name='fextract2d_sph')

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  integer, intent(in), value :: nbins
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),0:nc_p-1)
  real(rt), intent(inout) :: dens_bin(0:nbins-1)
  real(rt), intent(inout) :: vel_bin(0:nbins-1)
  real(rt), intent(inout) :: pres_bin(0:nbins-1)
  real(rt), intent(inout) :: e_bin(0:nbins-1)
  real(rt), intent(inout) :: volcount(0:nbins-1)
  integer, intent(in), value :: mask_size
  integer, intent(inout) :: imask(0:mask_size-1,0:mask_size-1)
  integer, intent(in), value :: r1, dens_comp, xmom_comp, ymom_comp, pres_comp, rhoe_comp
  real(rt), intent(in), value :: dx_fine, xctr, yctr
  real(rt), intent(in) :: dx(3)

  integer :: i, j, k, index
  real(rt) :: xx, xl, xr, yy, yl, yr, r_zone, vol, vel

  ! write(*,*) "imask = ", imask

  ! loop over all of the zones in the patch.  Here, we convert
  ! the cell-centered indices at the current level into the
  ! corresponding RANGE on the finest level, and test if we've
  ! stored data in any of those locations.  If we haven't then
  ! we store this level's data and mark that range as filled.
  k = lo(3)
  do j = lo(2), hi(2)
     yy = (dble(j) + HALF)*dx(2)
     yl = (dble(j))*dx(2)
     yr = (dble(j) + ONE)*dx(2)

     do i = lo(1), hi(1)
        xx = (dble(i) + HALF)*dx(1)
        xl = (dble(i))*dx(1)
        xr = (dble(i) + ONE)*dx(1)

        if ( any(imask(i*r1:(i+1)*r1-1, &
             j*r1:(j+1)*r1-1) .eq. 1) ) then

           r_zone = sqrt((xx-xctr)**2 + (yy-yctr)**2)

           index = r_zone/dx_fine

           vol = (xr**2 - xl**2)*(yr - yl)

           ! weight the zone's data by its size
           dens_bin(index) = dens_bin(index) + &
                p(i,j,k,dens_comp) * vol

           vel = sqrt(p(i,j,k,xmom_comp)**2 + &
                p(i,j,k,ymom_comp)**2) / &
                p(i,j,k,dens_comp)
           vel_bin(index) = vel_bin(index) + vel * vol

           pres_bin(index) = pres_bin(index) + &
                p(i,j,k,pres_comp) * vol

           e_bin(index) = e_bin(index) + &
                (p(i,j,k,rhoe_comp) / p(i,j,k,dens_comp)) * vol

           volcount(index) = volcount(index) + vol

           imask(i*r1:(i+1)*r1-1, &
                j*r1:(j+1)*r1-1) = 0

        end if

     enddo
  enddo

end subroutine fextract2d_sph


subroutine fextract3d_cyl(lo, hi, p, plo, phi, nc_p, nbins, dens_bin, &
     vel_bin, pres_bin, ncount, imask, mask_size, r1, &
     dens_comp, xmom_comp, ymom_comp, zmom_comp, pres_comp, dx_fine, dx, &
     xctr, yctr) bind(C, name='fextract3d_cyl')

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  integer, intent(in), value :: nbins
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),0:nc_p-1)
  real(rt), intent(inout) :: dens_bin(0:nbins-1)
  real(rt), intent(inout) :: vel_bin(0:nbins-1)
  real(rt), intent(inout) :: pres_bin(0:nbins-1)
  integer, intent(inout) :: ncount(0:nbins-1)
  integer, intent(in), value :: mask_size
  integer, intent(inout) :: imask(0:mask_size-1,0:mask_size-1,0:mask_size-1)
  integer, intent(in), value :: r1, dens_comp, xmom_comp, ymom_comp, zmom_comp, pres_comp
  real(rt), intent(in), value :: dx_fine, xctr, yctr
  real(rt), intent(in) :: dx(3)

  integer :: i, j, k, index
  real(rt) :: xx, yy, zz, r_zone

  ! loop over all of the zones in the patch.  Here, we convert
  ! the cell-centered indices at the current level into the
  ! corresponding RANGE on the finest level, and test if we've
  ! stored data in any of those locations.  If we haven't then
  ! we store this level's data and mark that range as filled.
  do k = lo(3), hi(3)
     zz = (k + HALF)*dx(3)
     do j = lo(2), hi(2)
        yy = (j + HALF)*dx(2)
        do i = lo(1), hi(1)
           xx = (i + HALF)*dx(1)

           if ( any(imask(i*r1:(i+1)*r1-1, &
                j*r1:(j+1)*r1-1, &
                k*r1:(k+1)*r1-1) .eq. 1) ) then

              r_zone = sqrt((xx-xctr)**2 + (yy-yctr)**2)

              index = r_zone/dx_fine

              dens_bin(index) = dens_bin(index) + &
                   p(i,j,k,dens_comp)*r1**3

              vel_bin(index) = vel_bin(index) + &
                   (sqrt(p(i,j,k,xmom_comp)**2 + &
                   p(i,j,k,ymom_comp)**2 + &
                   p(i,j,k,zmom_comp)**2)/ &
                   p(i,j,k,dens_comp))*r1**3

              pres_bin(index) = pres_bin(index) + &
                   p(i,j,k,pres_comp)*r1**3

              ncount(index) = ncount(index) + r1**3

              imask(i*r1:(i+1)*r1-1, &
                   j*r1:(j+1)*r1-1, &
                   k*r1:(k+1)*r1-1) = 0

           end if

        enddo
     enddo
  enddo

end subroutine fextract3d_cyl

subroutine fextract3d_sph(lo, hi, p, plo, phi, nc_p, nbins, dens_bin, &
     vel_bin, pres_bin, e_bin, ncount, imask, mask_size, r1,&
     dens_comp, xmom_comp, ymom_comp, zmom_comp, pres_comp, rhoe_comp, dx_fine, dx, &
     xctr, yctr, zctr) &
     bind(C, name='fextract3d_sph')

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  integer, intent(in), value :: nbins
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),0:nc_p-1)
  real(rt), intent(inout) :: dens_bin(0:nbins-1)
  real(rt), intent(inout) :: vel_bin(0:nbins-1)
  real(rt), intent(inout) :: pres_bin(0:nbins-1)
  real(rt), intent(inout) :: e_bin(0:nbins-1)
  integer, intent(inout) :: ncount(0:nbins-1)
  integer, intent(in), value :: mask_size
  integer, intent(inout) :: imask(0:mask_size-1,0:mask_size-1,0:mask_size-1)
  integer, intent(in), value :: r1, dens_comp, xmom_comp, ymom_comp, zmom_comp, pres_comp, rhoe_comp
  real(rt), intent(in), value :: dx_fine, xctr, yctr, zctr
  real(rt), intent(in) :: dx(3)

  integer :: i, j, k, index
  real(rt) :: xx, yy, zz, r_zone, vol, vel

  ! write(*,*) "imask = ", imask

  ! loop over all of the zones in the patch.  Here, we convert
  ! the cell-centered indices at the current level into the
  ! corresponding RANGE on the finest level, and test if we've
  ! stored data in any of those locations.  If we haven't then
  ! we store this level's data and mark that range as filled.
  do k = lo(3), hi(3)
     zz = (dble(k) + HALF)*dx(3)
     do j = lo(2), hi(2)
        yy = (dble(j) + HALF)*dx(2)

        do i = lo(1), hi(1)
           xx = (dble(i) + HALF)*dx(1)

           if ( any(imask(i*r1:(i+1)*r1-1, &
                j*r1:(j+1)*r1-1, &
                k*r1:(k+1)*r1-1) .eq. 1) ) then

              r_zone = sqrt((xx-xctr)**2 + (yy-yctr)**2 + (zz-zctr)**2)

              index = r_zone/dx_fine

              ! write(*,*) p(i,1,1,dens_comp)

              ! weight the zone's data by its size
              dens_bin(index) = dens_bin(index) + &
                   p(i,j,k,dens_comp)*r1**3

              vel_bin(index) = vel_bin(index) + &
                   (sqrt(p(i,j,k,xmom_comp)**2 + &
                   p(i,j,k,ymom_comp)**2 + &
                   p(i,j,k,zmom_comp)**2)/ &
                   p(i,j,k,dens_comp))*r1**3

              pres_bin(index) = pres_bin(index) + &
                   p(i,j,k,pres_comp)*r1**3

              e_bin(index) = e_bin(index) + &
                   (p(i,j,k,rhoe_comp)/p(i,j,k,dens_comp))*r1**3

              ncount(index) = ncount(index) + r1**3

              imask(i*r1:(i+1)*r1-1, &
                   j*r1:(j+1)*r1-1, &
                   k*r1:(k+1)*r1-1) = 0

           end if

        enddo
     enddo
  enddo

end subroutine fextract3d_sph
