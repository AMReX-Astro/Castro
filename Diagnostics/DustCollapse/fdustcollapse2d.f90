! Process a group of 1-d plotfiles from the dustcollapse problem
! and output the position of the interface as a function of time.
!
! The initial dense sphere is assumed to be centered a r = 0 (x = 0).
! We take as default that it is centered vertically at y = 0, but
! this can be overridden with --yctr.
!
! The --profile option will write out the average density vs. radius
! profile to a file (plotfile name + '.profile')
program fdustcollapse2d

  use f2kcli
  use bl_space
  use bl_error_module
  use bl_IO_module
  use plotfile_module

  implicit none

  type(plotfile) pf
  integer :: unit, uno
  integer :: i, j, ii, jj
  integer :: f

  real(kind=dp_t) :: xx, yy, xl, yl, xr, yr

  integer :: rr, r1

  integer :: nbins
  real(kind=dp_t), allocatable :: r(:), rl(:)
  real(kind=dp_t) :: maxdist, x_maxdist, y_maxdist
  real(kind=dp_t) :: yctr

  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  real(kind=dp_t) :: dx_fine

  real(kind=dp_t) :: r_zone, vol
  real(kind=dp_t) :: max_dens
  real(kind=dp_t) :: rho_lo,rho_hi,x,r_interface

  integer :: index_r

  real(kind=dp_t), pointer :: p(:,:,:,:)

  real(kind=dp_t), allocatable :: volcount(:)
  real(kind=dp_t), allocatable :: dens_bin(:)

  integer :: dens_comp

  logical, allocatable :: imask(:,:)
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  integer :: narg, farg
  character(len=256) :: fname, outfile
  integer :: indslsh

  logical :: profile

  unit = unit_new()
  uno  = unit_new()

  ! set the defaults
  yctr = 0.d0

  profile = .false.

  narg = command_argument_count()

  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('--yctr')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) yctr

     case ('--profile')
        profile = .true.

     case default
        exit

     end select
     farg = farg + 1
  end do

  do f = farg, narg

     call get_command_argument(f, value = fname)

     call build(pf, fname, unit)

     if (pf%dim /= 2) then
        print *, 'ERROR: not a 2-d file'
        stop
     endif

     ! density index
     dens_comp = plotfile_var_index(pf, "density")

     if (dens_comp < 0) then 
        call bl_error("ERROR: variable(s) not defined")
     endif

     ! get the index bounds and dx for the coarse level.  Note, lo and hi are
     ! ZERO based indicies
     lo = lwb(plotfile_get_pd_box(pf, 1))
     hi = upb(plotfile_get_pd_box(pf, 1))

     dx = plotfile_get_dx(pf, 1)

     ! get the index bounds for the finest level
     flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
     fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

     ! compute the size of the radially-binned array -- we'll do it to
     ! the furtherest corner of the domain
     x_maxdist = max(abs(pf%phi(1)       ), abs(pf%plo(1)       ))
     y_maxdist = max(abs(pf%phi(2) - yctr), abs(pf%plo(2) - yctr))
  
     maxdist = sqrt(x_maxdist**2 + y_maxdist**2)

     dx_fine = minval(plotfile_get_dx(pf, pf%flevel))
     nbins = int(maxdist/dx_fine)

     allocate(r(0:nbins-1))
     allocate(rl(0:nbins))

     do i = 0, nbins-1
        r(i) = (dble(i) + 0.5d0)*dx_fine
        rl(i) = dble(i)*dx_fine
     enddo
     rl(nbins) = dble(nbins)*dx_fine


     ! imask will be set to false if we've already output the data.
     ! Note, imask is defined in terms of the finest level.  As we loop
     ! over levels, we will compare to the finest level index space to
     ! determine if we've already output here
     allocate(imask(flo(1):fhi(1),flo(2):fhi(2)))

     imask(:,:) = .true.

     ! allocate storage for the data 
     allocate(dens_bin(0:nbins-1))
     allocate(volcount(0:nbins-1))

     volcount(:) = 0.d0
     dens_bin(:) = 0.d0

     ! loop over the data, starting at the finest grid, and if we haven't
     ! already store data in that grid location (according to imask),
     ! store it.  

     ! r1 is the factor between the current level grid spacing and the
     ! FINEST level
     r1  = 1

     do i = pf%flevel, 1, -1

        ! rr is the factor between the COARSEST level grid spacing and
        ! the current level
        rr = product(pf%refrat(1:i-1,1))

        do j = 1, nboxes(pf, i)

           ! read in the data 1 patch at a time
           call fab_bind(pf, i, j)

           lo = lwb(get_box(pf, i, j))
           hi = upb(get_box(pf, i, j))

           ! get a pointer to the current patch
           p => dataptr(pf, i, j)

           ! loop over all of the zones in the patch.  Here, we convert
           ! the cell-centered indices at the current level into the
           ! corresponding RANGE on the finest level, and test if we've
           ! stored data in any of those locations.  If we haven't then
           ! we store this level's data and mark that range as filled.
           do jj = lbound(p,dim=2), ubound(p,dim=2)
              yy = (dble(jj) + 0.5d0)*dx(2)/rr
              yl = (dble(jj))*dx(2)/rr
              yr = (dble(jj) + 1.d0)*dx(2)/rr

              do ii = lbound(p,dim=1), ubound(p,dim=1)
                 xx = (dble(ii) + 0.5d0)*dx(1)/rr
                 xl = (dble(ii))*dx(1)/rr
                 xr = (dble(ii) + 1.d0)*dx(1)/rr

                 if ( any(imask(ii*r1:(ii+1)*r1-1, &
                                jj*r1:(jj+1)*r1-1) ) ) then
                 
                    r_zone = sqrt((xx)**2 + (yy-yctr)**2)
                    
                    index_r = r_zone/dx_fine

                    vol = (xr**2 - xl**2)*(yr - yl)

                    ! weight the zone's data by its size
                    dens_bin(index_r) = dens_bin(index_r) + &
                         p(ii,jj,1,dens_comp) * vol

                    volcount(index_r) = volcount(index_r) + vol

                    imask(ii*r1:(ii+1)*r1-1, &
                         jj*r1:(jj+1)*r1-1) = .false.
                 
                 end if

              end do
           enddo

           call fab_unbind(pf, i, j)
        end do

        ! adjust r1 for the next lowest level
        if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
     end do

     ! normalize
     do i = 0, nbins-1
        if (volcount(i) /= 0.d0) then
           dens_bin(i) = dens_bin(i)/volcount(i)
        endif
     enddo

     ! These are calculated analytically given initial density 1.e9 and the
     !   analytic expression for the radius as a function of time
     ! t = 0.00
     if (abs(pf%tm-0.00) .le. 1.e-8) then
         max_dens = 1.e9 

     ! t = 0.01
     else if (abs(pf%tm-0.01) .le. 1.e-8) then
         max_dens = 1.043345e9

     ! t = 0.02
     else if (abs(pf%tm-0.02) .le. 1.e-8) then
         max_dens = 1.192524e9

     ! t = 0.03
     else if (abs(pf%tm-0.03) .le. 1.e-8) then
         max_dens = 1.527201e9

     ! t = 0.04
     else if (abs(pf%tm-0.04) .le. 1.e-8) then
         max_dens = 2.312884e9

     ! t = 0.05
     else if (abs(pf%tm-0.05) .le. 1.e-8) then
         max_dens = 4.779133e9

     ! t = 0.06
     else if (abs(pf%tm-0.06) .le. 1.e-8) then
         max_dens = 24.472425e9

     ! t = 0.065
     else if (abs(pf%tm-0.065) .le. 1.e-8) then
         max_dens = 423.447291e9

     else 
         print *,'Dont know the maximum density at this time: ',pf%tm
         stop
     end if

     ! print *,'Using a threshold of 0.5d0*max_dens: ',0.5d0*max_dens

     ! loop over the solution, from r = 0 outward, and find the first
     ! place where the density drops below the threshold density
     index_r = -1
     do i = 0, nbins-1
        if (dens_bin(i) < 0.5d0*max_dens) then
           index_r = i
           exit
        endif
     enddo

     if (index_r < 0) then
        print *, 'ERROR: density never fell below threshold'
        stop
     else if (index_r < 2) then
        r_interface = r(index_r)
     else

        rho_lo = dens_bin(index_r  )
        rho_hi = dens_bin(index_r-1)

        x = ( (0.5d0 * max_dens) - rho_lo) / (rho_hi - rho_lo)

        r_interface = x * r(index_r-1) + (1.d0-x) * r(index_r)

     endif

     ! output
     print *, pf%tm, r_interface

     ! dump out the radial profile -- if desired
     if (profile) then

        ! get the basename of the file -- remove the trailing slash if present
        indslsh = index(fname, '/', back = .TRUE.)
        
        if ( indslsh /= 0 ) then
           outfile = trim(fname(:indslsh-1)) // ".profile"
        else
           outfile = trim(fname) // ".profile"
        end if

        open (unit=uno, file=outfile, status = 'replace')
        do i = 0, nbins-1
           write (uno,*) r(i), dens_bin(i)
        enddo

     endif

     ! clean-up
     deallocate(r)
     deallocate(rl)
     deallocate(imask)
     deallocate(dens_bin)
     deallocate(volcount)

     call destroy(pf)

  enddo

end program fdustcollapse2d
