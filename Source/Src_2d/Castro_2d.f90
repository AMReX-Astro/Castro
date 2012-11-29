      subroutine ca_check_initial_species(lo,hi,&
                             state,state_l1,state_l2,state_h1,state_h2)

      use network           , only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(2), hi(2)
      integer          :: state_l1,state_l2,state_h1,state_h2
      double precision :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

      ! Local variables
      integer          :: i,j,n
      double precision :: sum

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)

         sum = 0.d0
         do n = 1, nspec
             sum = sum + state(i,j,UFS+n-1)
         end do
         if (abs(state(i,j,URHO)-sum).gt. 1.d-8 * state(i,j,URHO)) then
            print *,'Sum of (rho X)_i vs rho at (i,j): ',i,j,sum,state(i,j,URHO)
            call bl_error("Error:: Failed check of initial species summing to 1")
         end if

      enddo
      enddo

      end subroutine ca_check_initial_species

! :: ----------------------------------------------------------
! :: Volume-weight average the fine grid data onto the coarse
! :: grid.  Overlap is given in coarse grid coordinates.
! ::
! :: INPUTS / OUTPUTS:
! ::  crse      <=  coarse grid data
! ::  clo,chi    => index limits of crse array interior
! ::  ngc        => number of ghost cells in coarse array
! ::  nvar	 => number of components in arrays
! ::  fine       => fine grid data
! ::  flo,fhi    => index limits of fine array interior
! ::  ngf        => number of ghost cells in fine array
! ::  rfine      => (ignore) used in 2-D RZ calc
! ::  lo,hi      => index limits of overlap (crse grid)
! ::  lrat       => refinement ratio
! ::
! :: NOTE:
! ::  Assumes all data cell centered
! :: ----------------------------------------------------------
! ::
      subroutine ca_avgdown(crse,c_l1,c_l2,c_h1,c_h2,nvar, &
                            cv,cv_l1,cv_l2,cv_h1,cv_h2, &
                            fine,f_l1,f_l2,f_h1,f_h2, &
                            fv,fv_l1,fv_l2,fv_h1,fv_h2,lo,hi,lrat)
      implicit none
      integer c_l1,c_l2,c_h1,c_h2
      integer cv_l1,cv_l2,cv_h1,cv_h2
      integer f_l1,f_l2,f_h1,f_h2
      integer fv_l1,fv_l2,fv_h1,fv_h2
      integer lo(2), hi(2)
      integer nvar, lrat(2)
      double precision crse(c_l1:c_h1,c_l2:c_h2,nvar)
      double precision cv(cv_l1:cv_h1,cv_l2:cv_h2)
      double precision fine(f_l1:f_h1,f_l2:f_h2,nvar)
      double precision fv(fv_l1:fv_h1,fv_l2:fv_h2)

      integer i, j, n, ic, jc, ioff, joff
      integer lenx, leny, mxlen
      integer lratx, lraty

      lratx = lrat(1)
      lraty = lrat(2)
      lenx = hi(1)-lo(1)+1
      leny = hi(2)-lo(2)+1
      mxlen = max(lenx,leny)

      if (lenx .eq. mxlen) then
         do n = 1, nvar
 
!           Set coarse grid to zero on overlap
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,n) = 0.d0
               enddo
            enddo

!           Sum fine data
            do joff = 0, lraty-1
               do jc = lo(2), hi(2)
                  j = jc*lraty + joff
                  do ioff = 0, lratx-1
                     do ic = lo(1), hi(1)
                        i = ic*lratx + ioff
                        crse(ic,jc,n) = crse(ic,jc,n) + fv(i,j) * fine(i,j,n)
                     enddo
                  enddo
               enddo
            enddo

!           Divide out by volume weight
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,n) = crse(ic,jc,n) / cv(ic,jc)
               enddo
            enddo
            
         enddo

      else

         do n = 1, nvar

!           Set coarse grid to zero on overlap
            do ic = lo(1), hi(1)
               do jc = lo(2), hi(2)
                  crse(ic,jc,n) = 0.d0
               enddo
            enddo
 
!           Sum fine data
            do ioff = 0, lratx-1
               do ic = lo(1), hi(1)
                  i = ic*lratx + ioff
                  do joff = 0, lraty-1
                     do jc = lo(2), hi(2)
                        j = jc*lraty + joff
                        crse(ic,jc,n) = crse(ic,jc,n) + fv(i,j) * fine(i,j,n)
                     enddo
                  enddo
               enddo
            enddo
             
!           Divide out by volume weight
            do ic = lo(1), hi(1)
               do jc = lo(2), hi(2)
                  crse(ic,jc,n) = crse(ic,jc,n) / cv(ic,jc)
               enddo
            enddo
            
         enddo

      end if

      end subroutine ca_avgdown

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_compute_avgstate (lo,hi,dx,dr,nc,&
                                      state,s_l1,s_l2,s_h1,s_h2,radial_state, &
                                      vol,v_l1,v_l2,v_h1,v_h2,radial_vol, &
                                      problo,numpts_1d)

      use meth_params_module, only : URHO, UMX, UMY
      use probdata_module
      implicit none

      integer          :: lo(2),hi(2),nc
      double precision :: dx(2),dr,problo(2)

      integer          :: numpts_1d
      double precision :: radial_vol(0:numpts_1d-1)
      double precision :: radial_state(nc,0:numpts_1d-1)

      integer          :: s_l1,s_l2,s_h1,s_h2
      double precision :: state(s_l1:s_h1,s_l2:s_h2,nc)

      integer          :: v_l1,v_l2,v_h1,v_h2
      double precision :: vol(v_l1:v_h1,v_l2:v_h2)

      integer          :: i,j,n,index
      double precision :: x,y,r
      double precision :: x_mom,y_mom,radial_mom

      do j = lo(2), hi(2)
         y = problo(2) + (dble(j)+0.50d0) * dx(2) - center(2)
         do i = lo(1), hi(1)
            x = problo(1) + (dble(i)+0.50d0) * dx(1) - center(1)
            r = sqrt(x**2  + y**2)
            index = int(r/dr)
            if (index .gt. numpts_1d-1) then
              print *,'COMPUTE_AVGSTATE:INDEX TOO BIG ',index,' > ',numpts_1d-1
              print *,'AT (i,j) ',i,j
              print *,'R / DR IS ',r,dr
              call bl_error("Error:: Castro_2d.f90 :: ca_compute_avgstate")
            end if

            ! Store the radial component of the momentum in both the UMX and UMY components for now.
            x_mom = state(i,j,UMX)
            y_mom = state(i,j,UMY)
            radial_mom = x_mom * (x/r) + y_mom * (y/r)
            radial_state(UMX,index) = radial_state(UMX,index) + vol(i,j)*radial_mom
            radial_state(UMY,index) = radial_state(UMY,index) + vol(i,j)*radial_mom
           
            ! Store all the other variables as themselves
            radial_state(URHO,index) = radial_state(URHO,index) + vol(i,j)*state(i,j,URHO)
            do n = UMY+1,nc
               radial_state(n,index) = radial_state(n,index) + vol(i,j)*state(i,j,n)
            end do
            radial_vol(index) = radial_vol(index) + vol(i,j)
         enddo
      enddo

      end subroutine ca_compute_avgstate

! ::
! :: ----------------------------------------------------------
! ::

      subroutine normalize_species_fluxes(  &
                        flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                        flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                        lo,hi)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(2),hi(2)
      integer          :: flux1_l1,flux1_l2,flux1_h1,flux1_h2
      integer          :: flux2_l1,flux2_l2,flux2_h1,flux2_h2
      double precision :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
      double precision :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)

      ! Local variables
      integer          :: i,j,n
      double precision :: sum,fac

      do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               sum = 0.d0
               do n = UFS, UFS+nspec-1
                  sum = sum + flux1(i,j,n)
               end do
               if (sum .ne. 0.d0) then
                  fac = flux1(i,j,URHO) / sum
               else
                  fac = 1.d0
               end if
               do n = UFS, UFS+nspec-1
                  flux1(i,j,n) = flux1(i,j,n) * fac
               end do
            end do
      end do
      do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)
               sum = 0.d0
               do n = UFS, UFS+nspec-1
                  sum = sum + flux2(i,j,n)
               end do
               if (sum .ne. 0.d0) then
                  fac = flux2(i,j,URHO) / sum
               else
                  fac = 1.d0
               end if
               do n = UFS, UFS+nspec-1
                  flux2(i,j,n) = flux2(i,j,n) * fac
               end do
            end do
      end do

      end subroutine normalize_species_fluxes

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine enforce_minimum_density( uin,  uin_l1, uin_l2, uin_h1, uin_h2, &
                                          uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                                          lo, hi, mass_added, eint_added, eden_added, verbose)
      use network, only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UEINT, UEDEN, &
                                     UFS, UFX, UFA, small_dens, nadv

      implicit none

      integer          :: lo(2), hi(2), verbose
      integer          :: uin_l1,uin_l2,uin_h1,uin_h2
      integer          :: uout_l1,uout_l2,uout_h1,uout_h2
      double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
      double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
      double precision :: mass_added, eint_added, eden_added

      ! Local variables
      integer                       :: i,ii,j,jj,n
      double precision              :: min_dens
      double precision, allocatable :: fac(:,:)
      double precision              :: initial_mass, final_mass
      double precision              :: initial_eint, final_eint
      double precision              :: initial_eden, final_eden

      allocate(fac(lo(1):hi(1),lo(2):hi(2)))

      initial_mass = 0.d0
        final_mass = 0.d0
      initial_eint = 0.d0
        final_eint = 0.d0
      initial_eden = 0.d0
        final_eden = 0.d0

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)

            initial_mass = initial_mass + uout(i,j,URHO)
            initial_eint = initial_eint + uout(i,j,UEINT)
            initial_eden = initial_eden + uout(i,j,UEDEN)

            if (uout(i,j,URHO) .eq. 0.d0) then

               print *,'   '
               print *,'>>> Error: Castro_2d::enforce_minimum_density ',i,j
               print *,'>>> ... density exactly zero in grid ',lo(1),hi(1),lo(2),hi(2)
               print *,'    '
               call bl_error("Error:: Castro_2d.f90 :: enforce_minimum_density")

            else if (uout(i,j,URHO) < small_dens) then

               min_dens = uin(i,j,URHO)
               do jj = -1,1
               do ii = -1,1
                 min_dens = min(min_dens,uin(i+ii,j+jj,URHO))
                 if (uout(i+ii,j+jj,URHO) > small_dens) &
                   min_dens = min(min_dens,uout(i+ii,j+jj,URHO))
               end do
               end do

               if (verbose .gt. 0) then
                  if (uout(i,j,URHO) < 0.d0) then
                     print *,'   '
                     print *,'>>> Warning: Castro_2d::enforce_minimum_density ',i,j
                     print *,'>>> ... resetting negative density '
                     print *,'>>> ... from ',uout(i,j,URHO),' to ',min_dens
                     print *,'    '
                  else
                     print *,'   '
                     print *,'>>> Warning: Castro_2d::enforce_minimum_density ',i,j
                     print *,'>>> ... resetting small density '
                     print *,'>>> ... from ',uout(i,j,URHO),' to ',min_dens
                     print *,'    '
                  end if
               end if

               fac(i,j) = min_dens / uout(i,j,URHO)

            end if

         enddo
      enddo

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)

            if (uout(i,j,URHO) < small_dens) then

               uout(i,j,URHO ) = uout(i,j,URHO ) * fac(i,j)
               uout(i,j,UEINT) = uout(i,j,UEINT) * fac(i,j)
               uout(i,j,UEDEN) = uout(i,j,UEDEN) * fac(i,j)
               uout(i,j,UMX  ) = uout(i,j,UMX  ) * fac(i,j)
               uout(i,j,UMY  ) = uout(i,j,UMY  ) * fac(i,j)

               do n = UFS, UFS+nspec-1
                  uout(i,j,n) = uout(i,j,n) * fac(i,j)
               end do
               do n = UFX, UFX+naux-1
                  uout(i,j,n) = uout(i,j,n) * fac(i,j)
               end do
               do n = UFA, UFA+nadv-1
                  uout(i,j,n) = uout(i,j,n) * fac(i,j)
               end do

            end if

            final_mass = final_mass + uout(i,j,URHO)
            final_eint = final_eint + uout(i,j,UEINT)
            final_eden = final_eden + uout(i,j,UEDEN)

         enddo
      enddo

      mass_added = mass_added + (final_mass - initial_mass)
      eint_added = eint_added + (final_eint - initial_eint)
      eden_added = eden_added + (final_eden - initial_eden)

      end subroutine enforce_minimum_density

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(2), hi(2)
      integer          :: uout_l1,uout_l2,uout_h1,uout_h2
      double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)

      ! Local variables
      integer          :: i,j,n
      integer          :: int_dom_spec
      logical          :: any_negative
      double precision :: dom_spec,x,eps

      eps = -1.0d-16

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

         any_negative = .false.

         ! First deal with tiny undershoots by just setting them to zero
         do n = UFS, UFS+nspec-1
           if (uout(i,j,n) .lt. 0.d0) then
              x = uout(i,j,n)/uout(i,j,URHO)
              if (x .gt. eps) then
                 uout(i,j,n) = 0.d0
              else
                 any_negative = .true.
              end if
           end if
         end do

         ! We know there are one or more undershoots needing correction 
         if (any_negative) then

            ! Find the dominant species
            int_dom_spec  = UFS
            dom_spec      = uout(i,j,int_dom_spec)

            do n = UFS,UFS+nspec-1
              if (uout(i,j,n) .gt. dom_spec) then
                dom_spec = uout(i,j,n)
                int_dom_spec = n
              end if
            end do

           ! Now take care of undershoots greater in magnitude than 1e-16.
           do n = UFS, UFS+nspec-1

              if (uout(i,j,n) .lt. 0.d0) then

                 x = uout(i,j,n)/uout(i,j,URHO)

                 ! Here we only print the bigger negative values
                 if (x .lt. -1.d-2) then
                    print *,'At cell (i,j) = ',i,j
                    print *,'... Fixing negative species ',n           ,' with X = ',x
                    print *,'...   from dominant species ',int_dom_spec,' with X = ',&
                             uout(i,j,int_dom_spec) / uout(i,j,URHO)
                 end if

                 ! Take enough from the dominant species to fill the negative one.
                 uout(i,j,int_dom_spec) = uout(i,j,int_dom_spec) + uout(i,j,n)
   
                 ! Test that we didn't make the dominant species negative
                 if (uout(i,j,int_dom_spec) .lt. 0.d0) then 
                    print *,'Just made dominant species negative ',int_dom_spec,' at ',i,j
                    print *,'... We were fixing species ',n,' which had value ',x
                    print *,'... Dominant species became ',uout(i,j,int_dom_spec) / uout(i,j,URHO)
                    call bl_error("Error:: Castro_2d.f90 :: ca_enforce_nonnegative_species")
                 end if

                 ! Now the negative species to zero
                 uout(i,j,n) = 0.d0

              end if

           enddo
         end if

      enddo
      enddo

      end subroutine ca_enforce_nonnegative_species

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine normalize_new_species(u,u_l1,u_l2,u_h1,u_h2,lo,hi)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(2), hi(2)
      integer          :: u_l1,u_l2,u_h1,u_h2
      double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)

      ! Local variables
      integer          :: i,j,n
      double precision :: fac,sum

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            sum = 0.d0
            do n = UFS, UFS+nspec-1
               sum = sum + u(i,j,n)
            end do
            if (sum .ne. 0.d0) then
               fac = u(i,j,URHO) / sum
            else
               fac = 1.d0
            end if
            do n = UFS, UFS+nspec-1
               u(i,j,n) = u(i,j,n) * fac
            end do
         end do
      end do

      end subroutine normalize_new_species


! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine reset_internal_e(u,u_l1,u_l2,u_h1,u_h2,lo,hi,verbose)

      use eos_module
      use network, only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UFX, &
                                     small_temp, allow_negative_energy

      implicit none

      integer          :: lo(2), hi(2), verbose
      integer          :: u_l1,u_l2,u_h1,u_h2
      double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)

      ! Local variables
      integer          :: i,j
      integer          :: pt_index(2)
      double precision :: Up, Vp, ke
      double precision :: rho_eint, eint_new
      double precision :: x_in(1:nspec+naux), dummy_pres

      ! Reset internal energy if negative.
      if (allow_negative_energy .eq. 0) then

         do j = lo(2),hi(2)
         do i = lo(1),hi(1)

           Up = u(i,j,UMX) / u(i,j,URHO)
           Vp = u(i,j,UMY) / u(i,j,URHO)

           ke = 0.5d0 * u(i,j,URHO) * (Up**2 + Vp**2)

           rho_eint = u(i,j,UEDEN) - ke

           ! Reset (rho e) if e is greater than 0.01% of E.
           if (rho_eint .gt. 0.d0 .and. rho_eint / u(i,j,UEDEN) .gt. 1.d-4) then

               u(i,j,UEINT) = rho_eint

           ! If (e from E) < 0 or (e from E) < .0001*E but (e from e) > 0.
           else if (u(i,j,UEINT) .gt. 0.d0) then

               u(i,j,UEDEN) = u(i,j,UEINT) + ke

           ! If not resetting and little e is negative ...
           else if (u(i,j,UEINT) .le. 0.d0) then

              x_in(1:nspec) = u(i,j,UFS:UFS+nspec-1) / u(i,j,URHO)
              if (naux > 0) &
                x_in(nspec+1:nspec+naux)  = u(i,j,UFX:UFX+naux -1) / u(i,j,URHO)

              pt_index(1) = i
              pt_index(2) = j

              call eos_given_RTX(eint_new, dummy_pres, u(i,j,URHO), small_temp, x_in, pt_index=pt_index)
              if (verbose .gt. 0) then
                 print *,'   '
                 print *,'>>> Warning: Castro_2d::reset_internal_energy  ',i,j 
                 print *,'>>> ... resetting neg. e from EOS using small_temp'
                 print *,'>>> ... from ',u(i,j,UEINT)/u(i,j,URHO),' to ', eint_new
                 print *,'    '
              end if

              u(i,j,UEINT) = u(i,j,URHO) * eint_new
              u(i,j,UEDEN) = u(i,j,URHO) * eint_new + ke

           end if

         enddo
         enddo

      ! If (allow_negative_energy .eq. 1) then just reset (rho e) from (rho E)
      else

         do j = lo(2),hi(2)
         do i = lo(1),hi(1)

           Up = u(i,j,UMX) / u(i,j,URHO)
           Vp = u(i,j,UMY) / u(i,j,URHO)
           ke = 0.5d0 * u(i,j,URHO) * (Up**2 + Vp**2)

           u(i,j,UEINT) = u(i,j,UEDEN) - ke

         enddo
         enddo

      endif

    end subroutine reset_internal_e

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine get_center(center_out)

        use probdata_module, only : center

        implicit none

        double precision :: center_out(2)

        center_out(1:2) = center(1:2)

      end subroutine get_center

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine set_center(center_in)

        use probdata_module, only : center

        implicit none

        double precision :: center_in(2)

        center(1:2) = center_in(1:2)

      end subroutine set_center

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine find_center(data,new_center,icen,dx,problo)

        implicit none

        double precision :: data(-1:1,-1:1)
        double precision :: new_center(2)
        double precision :: dx(2),problo(2)
        double precision :: a,b,x,y,cen
        integer          :: icen(2)
        integer          :: i,j

        ! We do this to take care of precision issues
        cen = data(0,0)
        do j = -1,1
        do i = -1,1
           data(i,j) = data(i,j) - cen 
        end do
        end do

!       This puts the center at the cell center
        new_center(1) = problo(1) +  (icen(1)+0.5d0) * dx(1)
        new_center(2) = problo(2) +  (icen(2)+0.5d0) * dx(2)
   
        ! Fit parabola y = a x^2  + b x + c through three points
        ! a = 1/2 ( y_1 + y_-1)
        ! b = 1/2 ( y_1 - y_-1)
        ! x_vertex = -b / 2a

        ! ... in x-direction
        a = 0.5d0 * (data(1,0) + data(-1,0)) - data(0,0)
        b = 0.5d0 * (data(1,0) - data(-1,0)) - data(0,0)
        x = -b / (2.d0*a)
        new_center(1) = new_center(1) +  x*dx(1)

        ! ... in y-direction
        a = 0.5d0 * (data(0,1) + data(0,-1)) - data(0,0)
        b = 0.5d0 * (data(0,1) - data(0,-1)) - data(0,0)
        y = -b / (2.d0*a)
        new_center(2) = new_center(2) +  y*dx(2)

      end subroutine find_center

