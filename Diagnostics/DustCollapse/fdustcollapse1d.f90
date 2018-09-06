! Process a group of 1-d plotfiles from the dustcollapse problem
! and output the position of the interface as a function of time.

program fdustcollapse1d

  use f2kcli
  use amrex_error_module
  use plotfile_module
  use sort_d_module

  implicit none

  type(plotfile) pf
  integer :: unit
  integer :: i, j, ii
  integer :: f

  integer :: rr, r1

  integer :: cnt
  integer :: max_points

  real(rt) :: dx(MAX_SPACEDIM)

  integer :: index

  real(rt), pointer :: p(:,:,:,:)

  real(rt), allocatable :: rcoord(:), dens(:)
  integer, allocatable :: isort(:)

  integer :: dens_comp

  logical, allocatable :: imask(:)
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  real(rt) :: rmin
  real(rt) :: max_dens
  real(rt) :: rho_lo,rho_hi,x,r_interface

  integer :: narg, farg
  character(len=256) :: fname

  unit = unit_new()

  narg = command_argument_count()

  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case default
        exit

     end select
     farg = farg + 1
  end do

  print *, 'Using a density threshhold of half the maximum analytic density'

  do f = farg, narg

     call get_command_argument(f, value = fname)

     call build(pf, fname, unit)

     if (pf%dim /= 1) then
        print *, 'ERROR: not a 1-d file'
        stop
     endif

     rmin = pf%plo(1)

  
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

     ! compute the maximum number of zones, as if we were completely refined
     max_points = fhi(1) - flo(1) + 1


     ! imask will be set to false if we've already output the data.
     ! Note, imask is defined in terms of the finest level.  As we loop
     ! over levels, we will compare to the finest level index space to
     ! determine if we've already output here
     allocate(imask(flo(1):fhi(1)))

     imask(:) = .true.

     ! allocate storage for the data 
     allocate(rcoord(max_points))
     allocate(  dens(max_points))
     allocate( isort(max_points))

     rcoord(:) = 0.d0
     dens(:)   = 0.d0
     isort(:)  = 0.d0

     ! loop over the data, starting at the finest grid, and if we haven't
     ! already store data in that grid location (according to imask),
     ! store it. 

     cnt = 0


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
           do ii = lbound(p,dim=1), ubound(p,dim=1)

              if ( any(imask(ii*r1:(ii+1)*r1-1) ) ) then

                 cnt = cnt + 1

                 rcoord(cnt) = rmin + (ii + 0.5d0)*dx(1)/rr
                 dens(cnt)   = p(ii,1,1,dens_comp)

                 imask(ii*r1:(ii+1)*r1-1) = .false.
                 
              end if

           enddo

           call fab_unbind(pf, i, j)
        end do

        ! adjust r1 for the next lowest level
        if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
     end do

     ! sort the data based on the coordinates
     call sort(rcoord(1:cnt),isort(1:cnt))

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

     ! loop over the solution, from r = 0 outward, and find the first
     ! place where the density drops below the threshold density
     index = -1
     do i = 1, cnt
        if (dens(isort(i)) < 0.5d0*max_dens) then
           index = i
           exit
        endif
     enddo

     if (index < 0) then
        print *, 'ERROR: density never fell below threshold'
        stop
     else if (index < 2) then
        r_interface = rcoord(index)
     else

        rho_lo = dens(index  )
        rho_hi = dens(index-1)

        x = ( (0.5d0 * max_dens) - rho_lo) / (rho_hi - rho_lo)

        r_interface = x * rcoord(index-1) + (1.d0-x) * rcoord(index)

     endif

     ! output
     print *, pf%tm, r_interface

     ! clean-up
     deallocate(rcoord)
     deallocate(dens)
     deallocate(isort)
     deallocate(imask)

     call destroy(pf)
     
  enddo

end program fdustcollapse1d
