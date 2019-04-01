! Process a 2-d gaussian radiation pulse
subroutine fgaussian_pulse(lo, hi, p, plo, phi, nc_p, nbins, rad_bin, &
     ncount, imask, mask_size, r1,&
     rad_comp, dx, dx_fine, xctr, yctr) bind(C, name='fgaussian_pulse')

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  integer, intent(in), value :: nbins
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),0:nc_p-1)
  real(rt), intent(inout) :: rad_bin(0:nbins-1)
  integer, intent(inout) :: ncount(0:nbins-1)
  integer, intent(inout) :: imask(0:mask_size-1,0:mask_size-1)
  integer, intent(in), value :: mask_size, r1, rad_comp
  real(rt), intent(in), value :: dx_fine, xctr, yctr
  real(rt), intent(in) :: dx(3)

  integer :: ii, jj, k, index
  real(rt) :: xx, yy, r_zone

  ! loop over all of the zones in the patch.  Here, we convert
  ! the cell-centered indices at the current level into the
  ! corresponding RANGE on the finest level, and test if we've
  ! stored data in any of those locations.  If we haven't then
  ! we store this level's data and mark that range as filled.
  k = lo(3)
  do jj = lo(2), hi(2)
     yy = (jj + HALF)*dx(2)
     do ii = lo(1), hi(1)
        xx = (ii + HALF)*dx(1)

        if ( any(imask(ii*r1:(ii+1)*r1-1, &
             jj*r1:(jj+1)*r1-1) .eq. 1) ) then

           r_zone = sqrt((xx-xctr)**2 + (yy-yctr)**2)

           index = r_zone/dx_fine

           ! weight the zone's data by its size
           rad_bin(index) = rad_bin(index) + &
                p(ii,jj,k,rad_comp)*r1**2

           ncount(index) = ncount(index) + r1**2

           imask(ii*r1:(ii+1)*r1-1, &
                jj*r1:(jj+1)*r1-1) = 0

        end if

     enddo
  enddo

end subroutine fgaussian_pulse


! Process a 1-d sedov problem to produce rho, u, and p as a
! function of r, for comparison to the analytic solution.
subroutine flgt_frnt1d(lo, hi, p, plo, phi, nc_p, nbins, &
     dens_bin, vel_bin, pres_bin, rad_bin, &
     imask, mask_size, r1,&
     dens_comp, xmom_comp, pres_comp, &
     rad_comp, dx, dx_fine) bind(C, name='flgt_frnt1d')

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
  real(rt), intent(inout) :: rad_bin(0:nbins-1)
  integer, intent(inout) :: imask(0:mask_size-1)
  integer, intent(in), value :: mask_size, r1, dens_comp, xmom_comp, pres_comp, rad_comp
  real(rt), intent(in), value :: dx_fine
  real(rt), intent(in) :: dx(3)

  integer :: ii, j, k, index

  ! loop over all of the zones in the patch.  Here, we convert
  ! the cell-centered indices at the current level into the
  ! corresponding RANGE on the finest level, and test if we've
  ! stored data in any of those locations.  If we haven't then
  ! we store this level's data and mark that range as filled.
  k = lo(3)
  j = lo(2)

  do ii = lo(1), hi(1)

     if ( any(imask(ii*r1:(ii+1)*r1-1) .eq. 1) )then

        index = ii * r1

        dens_bin(index:index+(r1-1)) = p(ii,j,k,dens_comp)

        vel_bin(index:index+(r1-1)) = &
             abs(p(ii,j,k,xmom_comp)) / p(ii,j,k,dens_comp)

        pres_bin(index:index+(r1-1)) = p(ii,j,k,pres_comp)

        rad_bin(index:index+(r1-1)) = p(ii,j,k,rad_comp)

        imask(ii*r1:(ii+1)*r1-1) = 0

     end if

  enddo

end subroutine flgt_frnt1d

! extract a 1-d slice of the data (all variables or a single variable)
! along the specified coordinate direction from a plotfile.  The
! plotfile can be 1-, 2-, or 3-d.
!
! This routine is a generalized version is based on fextract3d, but geared
! toward the CASTRO radiating shock problem
!
! We read in all the variables, but only output a subset
subroutine fradshock(lo, hi, problo, probhi, p, plo, phi, nc_p, nbins, &
     vars_bin, imask, mask_size, r1, rr, dx, idir, cnt) bind(C, name='fradshock')

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  real(rt), intent(in) :: problo(3), probhi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  integer, intent(in), value :: nbins
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),0:nc_p-1)
  real(rt), intent(inout) :: vars_bin(0:nbins-1, 0:nc_p)
  integer, intent(inout) :: imask(0:mask_size-1)
  integer, intent(in), value :: mask_size, r1, rr, idir
  integer, intent(inout) :: cnt
  real(rt), intent(in) :: dx(3)

  integer :: ii, jj, kk
  real(rt) :: iloc, jloc, kloc, xmin, xmax, ymin, ymax, zmin, zmax

  iloc = (hi(1)-lo(1)+1)/2 + lo(1)

  xmin = problo(1)
  xmax = probhi(1)

#if (AMREX_SPACEDIM >= 2)
  jloc = (hi(2)-lo(2)+1)/2 + lo(2)
  ymin = problo(2)
  ymax = probhi(2)
  kloc = (hi(3)-lo(3)+1)/2 + lo(3)
  zmin = problo(3)
  zmax = probhi(3)
#endif

  select case (idir)

  case (1)

#if (AMREX_SPACEDIM == 1)
     jj = lo(2)
     kk = lo(3)

     do ii = lo(1), hi(1)
        if ( any(imask(ii*r1:(ii+1)*r1-1) .eq. 1) ) then
           vars_bin(cnt,0) = xmin + (ii + HALF)*dx(1)
           vars_bin(cnt,1:) = p(ii,jj,kk,:)

           imask(ii*r1:(ii+1)*r1-1) = 0

           cnt = cnt + 1

        end if
     end do

#else

     ! if the current patch stradles our slice, then get a data
     ! pointer to it
     if ( rr*jloc >= lo(2) .and. rr*jloc <= hi(2) .and. &
          ( (AMREX_SPACEDIM .eq. 2) .or. (rr*kloc >= lo(3) .and. rr*kloc <= hi(3)) ) ) then

        jj = jloc*rr

        if (AMREX_SPACEDIM .eq. 3) then
           kk = kloc*rr
        else
           kk = lo(3)
        end if

        ! loop over all of the zones in the slice direction.
        ! Here, we convert the cell-centered indices at the
        ! current level into the corresponding RANGE on the
        ! finest level, and test if we've stored data in any of
        ! those locations.  If we haven't then we store this
        ! level's data and mark that range as filled.
        do ii = lo(1), hi(1)
           if ( any(imask(ii*r1:(ii+1)*r1-1) .eq. 1 ) ) then

              vars_bin(cnt,0) = xmin + (ii + HALF)*dx(1)
              vars_bin(cnt,1:) = p(ii,jj,kk,:)

              imask(ii*r1:(ii+1)*r1-1) = 0
              cnt = cnt + 1
           end if
        end do

     end if

#endif

  case (2)

     ! if the current patch stradles our slice, then get a data
     ! pointer to it
     if ( rr*iloc >= lo(1) .and. rr*iloc <= hi(1) .and. &
          ( (AMREX_SPACEDIM .eq. 2) .or. (rr*kloc >= lo(3) .and. rr*kloc <= hi(3)) ) ) then

        ii = iloc*rr
        if (AMREX_SPACEDIM .eq. 3) then
           kk = kloc*rr
        else
           kk = lo(3)
        end if

        ! loop over all of the zones in the slice direction.
        ! Here, we convert the cell-centered indices at the
        ! current level into the corresponding RANGE on the
        ! finest level, and test if we've stored data in any of
        ! those locations.  If we haven't then we store this
        ! level's data and mark that range as filled.
        do jj = lo(2), hi(2)
           if ( any(imask(jj*r1:(jj+1)*r1-1) .eq. 1) ) then

              vars_bin(cnt,0) = ymin + (jj + HALF)*dx(2)
              vars_bin(cnt,1:) = p(ii,jj,kk,:)

              imask(jj*r1:(jj+1)*r1-1) = 0
              cnt = cnt + 1
           end if
        end do

     end if

  case (3)

     ! if the current patch stradles our slice, then get a data
     ! pointer to it
     if ( rr*iloc >= lo(1) .and. rr*iloc <= hi(1) .and. &
          rr*jloc >= lo(2) .and. rr*jloc <= hi(2)) then

        ii = iloc*rr
        jj = jloc*rr


        ! loop over all of the zones in the slice direction.
        ! Here, we convert the cell-centered indices at the
        ! current level into the corresponding RANGE on the
        ! finest level, and test if we've stored data in any of
        ! those locations.  If we haven't then we store this
        ! level's data and mark that range as filled.
        do kk = lo(3), hi(3)
           if ( any(imask(kk*r1:(kk+1)*r1-1) .eq. 1) ) then

              vars_bin(cnt,0) = zmin + (kk + HALF)*dx(3)
              vars_bin(cnt,1:) = p(ii,jj,kk,:)

              imask(kk*r1:(kk+1)*r1-1) = 0
              cnt = cnt + 1
           end if
        end do

     end if

  end select

end subroutine fradshock


! Analysis routine for the radiation source test.
!
! This problem is a thermal relaxiation problem.  The domain is
! completely uniform, so we just need to look at the state variables
! in a single zone.
!
! Take a list of files and print out (rho e) and the total radiation
! energy density in the first zone as a function of time.
subroutine fradsource(lo, hi, p, plo, phi, nc_p, &
     rhoe, rad, rhoe_comp, rad_comp) bind(C, name='fradsource')

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),0:nc_p-1)
  real(rt), intent(inout) :: rhoe, rad
  integer, intent(in), value :: rhoe_comp, rad_comp

  integer :: i, j, k

  k = lo(3)
  j = lo(2)
  i = lo(1)

  rhoe = p(i,j,k,rhoe_comp)
  rad = p(i,j,k,rad_comp)

end subroutine fradsource


subroutine fradsphere(lo, hi, problo, probhi, p, plo, phi, nc_p, nbins, &
     vars_bin, imask, mask_size, r1, dx, cnt) bind(C, name='fradsphere')

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  real(rt), intent(in) :: problo(3), probhi(3)
  integer, intent(in) :: plo(3), phi(3), nc_p
  integer, intent(in), value :: nbins
  real(rt), intent(in) :: p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3),0:nc_p-1)
  real(rt), intent(inout) :: vars_bin(0:nbins-1, 0:nc_p)
  integer, intent(inout) :: imask(0:mask_size-1)
  integer, intent(in), value :: mask_size, r1
  integer, intent(inout) :: cnt
  real(rt), intent(in) :: dx(3)

  integer :: i, j, k
  real(rt) :: rmin, rmax

  rmin = problo(1)
  rmax = probhi(1)

  j = lo(2)
  k = lo(3)
  do i = lo(1), hi(1)
     if ( any(imask(i*r1:(i+1)*r1-1) .eq. 1)) then

        vars_bin(cnt,0) = rmin + (i + HALF)*dx(1)
        vars_bin(cnt,1:) = p(i,j,k,:)

        imask(i*r1:(i+1)*r1-1) = 0

        cnt = cnt + 1
     end if
  end do

end subroutine fradsphere

! Analysis routine for RHD_shocktube
subroutine frhdshocktube(lo, hi, p, plo, phi, nc_p, nbins, &
     dens_bin, vel_bin, pres_bin, rad_bin, &
     dens_comp, xvel_comp, pres_comp, rad_comp, ngroups) &
     bind(C, name='frhdshocktube')

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
  real(rt), intent(inout) :: rad_bin(0:nbins-1)
  integer, intent(in), value :: dens_comp, xvel_comp, pres_comp, rad_comp, ngroups

  integer :: i, j, k, ig

  k = lo(3)
  j = lo(2)
  do i = lo(1), hi(1)

     dens_bin(i) = p(i,j,k,dens_comp)
     vel_bin(i) = p(i,j,k,xvel_comp)
     pres_bin(i) = p(i,j,k,pres_comp)

     do ig = 0, ngroups-1
         rad_bin(i) = rad_bin(i) + p(i,j,k,rad_comp+ig)
     enddo

  enddo

end subroutine frhdshocktube
