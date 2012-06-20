      subroutine ca_check_initial_species(lo,hi,state,state_l1,state_h1)

      use network           , only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(1), hi(1)
      integer          :: state_l1,state_h1
      double precision :: state(state_l1:state_h1,NVAR)

      ! Local variables
      integer          :: i,n
      double precision :: sum

      do i = lo(1), hi(1)

         sum = 0.d0
         do n = 1, nspec
             sum = sum + state(i,UFS+n-1)
         end do
         if (abs(state(i,URHO)-sum).gt. 1.d-8 * state(i,URHO)) then
            print *,'Sum of (rho X)_n vs rho at (i): ',i,sum,state(i,URHO)
            call bl_error("Error:: Failed check of initial species summing to 1")
         end if

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
      subroutine ca_avgdown (crse,c_l1,c_h1,nvar, &
           cv,cv_l1,cv_h1, &
           fine,f_l1,f_h1, &
           fv,fv_l1,fv_h1,lo,hi,lrat)

      implicit none
      integer c_l1,c_h1
      integer cv_l1,cv_h1
      integer f_l1,f_h1
      integer fv_l1,fv_h1
      integer lo(1), hi(1)
      integer nvar, lrat(1)
      double precision crse(c_l1:c_h1,nvar)
      double precision cv(cv_l1:cv_h1)
      double precision fine(f_l1:f_h1,nvar)
      double precision fv(fv_l1:fv_h1)

      integer i, n, ic, ioff
      integer lratx
 
      lratx = lrat(1)
 
      do n = 1, nvar
 
!        Set coarse grid to zero on overlap
         do ic = lo(1), hi(1)
            crse(ic,n) = 0.d0
         enddo
 
 
!        Sum fine data
         do ioff = 0, lratx-1
            do ic = lo(1), hi(1)
               i = ic*lratx + ioff
               crse(ic,n) = crse(ic,n) + fv(i) * fine(i,n)
            enddo
         enddo
             
!        Divide out by volume weight
         do ic = lo(1), hi(1)
            crse(ic,n) = crse(ic,n) / cv(ic)
         enddo
            
      enddo

      end subroutine ca_avgdown

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(1),hi(1)
      integer          :: flux_l1,flux_h1
      double precision :: flux(flux_l1:flux_h1,NVAR)

      ! Local variables
      integer          :: i,n
      double precision :: sum,fac

      do i = lo(1),hi(1)+1
         sum = 0.d0
         do n = UFS, UFS+nspec-1
            sum = sum + flux(i,n)
         end do
         if (sum .ne. 0.d0) then
            fac = flux(i,URHO) / sum
         else
            fac = 1.d0
         end if
         do n = UFS, UFS+nspec-1
            flux(i,n) = flux(i,n) * fac
         end do
      end do

      end subroutine normalize_species_fluxes

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine enforce_minimum_density( uin,  uin_l1, uin_h1, &
                                          uout,uout_l1,uout_h1,&
                                          lo, hi, verbose)
      use network, only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UFS, UFX, UFA, small_dens, nadv

      implicit none

      integer          :: lo(1), hi(1), verbose
      integer          :: uin_l1 , uin_h1
      integer          :: uout_l1,uout_h1
      double precision :: uin(uin_l1:uin_h1,NVAR)
      double precision :: uout(uout_l1:uout_h1,NVAR)

      ! Local variables
      integer          :: i,ii,n
      double precision :: min_dens
      double precision, allocatable :: fac(:)

      allocate(fac(lo(1):hi(1)))

      do i = lo(1),hi(1)

         if (uout(i,URHO) .eq. 0.d0) then

            print *,'   '
            print *,'>>> Error: Castro_1d::enforce_minimum_density ',i
            print *,'>>> ... density exactly zero in grid ',lo(1),hi(1)
            print *,'    '
            call bl_error("Error:: Castro_1d.f90 :: enforce_minimum_density")

         else if (uout(i,URHO) < small_dens) then

            min_dens = uin(i,URHO)
            do ii = -1,1
              min_dens = min(min_dens,uin(i+ii,URHO))
              if (uout(i+ii,URHO) > small_dens) &
                min_dens = min(min_dens,uout(i+ii,URHO))
            end do

            if (verbose .gt. 0) then
               if (uout(i,URHO) < 0.d0) then
                  print *,'   '
                  print *,'>>> Warning: Castro_1d::enforce_minimum_density ',i
                  print *,'>>> ... resetting negative density '
                  print *,'>>> ... from ',uout(i,URHO),' to ',min_dens
                  print *,'    '
               else
                  print *,'   '
                  print *,'>>> Warning: Castro_1d::enforce_minimum_density ',i
                  print *,'>>> ... resetting small density '
                  print *,'>>> ... from ',uout(i,URHO),' to ',min_dens
                  print *,'    '
               end if
            end if

            fac(i) = min_dens / uout(i,URHO)

         end if

      enddo

      do i = lo(1),hi(1)

         if (uout(i,URHO) < small_dens) then

            uout(i,URHO ) = uout(i,URHO ) * fac(i)
            uout(i,UEDEN) = uout(i,UEDEN) * fac(i)
            uout(i,UEINT) = uout(i,UEINT) * fac(i)
            uout(i,UMX  ) = uout(i,UMX  ) * fac(i)

            do n = UFS, UFS+nspec-1
               uout(i,n) = uout(i,n) * fac(i)
            end do
            do n = UFX, UFX+naux-1
               uout(i,n) = uout(i,n) * fac(i)
            end do
            do n = UFA, UFA+nadv-1
               uout(i,n) = uout(i,n) * fac(i)
            end do

         end if

      enddo

      deallocate(fac)

      end subroutine enforce_minimum_density

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine ca_enforce_nonnegative_species(uout,uout_l1,uout_h1,lo,hi)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(1), hi(1)
      integer          :: uout_l1,uout_h1
      double precision :: uout(uout_l1:uout_h1,NVAR)

      ! Local variables
      integer          :: i,n,nn
      integer          :: int_dom_spec
      logical          :: any_negative
      double precision :: dom_spec,x,eps

      eps = -1.d-16

      do i = lo(1),hi(1)

         any_negative = .false.

         ! First deal with tiny undershoots by just setting them to zero
         do n = UFS, UFS+nspec-1
           if (uout(i,n) .lt. 0.d0) then
              x = uout(i,n)/uout(i,URHO)
              if (x .gt. eps) then
                 uout(i,n) = 0.d0
              else
                 any_negative = .true.
              end if
           end if
         end do

         ! We know there are one or more undershoots needing correction 
         if (any_negative) then

            ! Find the dominant species
            dom_spec = 0.d0
            int_dom_spec = 0
            do n = UFS,UFS+nspec-1
              if (uout(i,n) .gt. dom_spec) then
                dom_spec = uout(i,n)
                int_dom_spec = n
              end if
            end do

           ! Now take care of undershoots greater in magnitude than 1e-16.
           do n = UFS, UFS+nspec-1

              if (uout(i,n) .lt. 0.d0) then

                 x = uout(i,n)/uout(i,URHO)

                 ! Here we only print the bigger negative values
                 if (x .lt. -1.d-2) then
                    print *,'Correcting negative species   ',n
                    print *,'   at cell (i)                ',i
                    print *,'Negative (rho*X) is           ',uout(i,n)
                    print *,'Negative      X  is           ',x
                    print *,'Filling from dominant species ',int_dom_spec
                    print *,'  which had X =               ',&
                             uout(i,int_dom_spec) / uout(i,URHO)
                 end if

                 ! Take enough from the dominant species to fill the negative one.
                 uout(i,int_dom_spec) = uout(i,int_dom_spec) + uout(i,n)
   
                 ! Test that we didn't make the dominant species negative
                 if (uout(i,int_dom_spec) .lt. 0.d0) then 
                    print *,' Just made dominant species negative ',int_dom_spec,' at ',i
                    print *,'We were fixing species ',n,' which had value ',x
                    print *,'Dominant species became ',uout(i,int_dom_spec) / uout(i,URHO)
                    call bl_error("Error:: Castro_2d.f90 :: ca_enforce_nonnegative_species")
                 end if

                 ! Now set the negative species to zero
                 uout(i,n) = 0.d0

              end if

           enddo
         end if

      enddo

      end subroutine ca_enforce_nonnegative_species

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine normalize_new_species(u,u_l1,u_h1,lo,hi)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(1), hi(1)
      integer          :: u_l1,u_h1
      double precision :: u(u_l1:u_h1,NVAR)

      ! Local variables
      integer          :: i,n
      double precision :: fac,sum

      do i = lo(1),hi(1)
         sum = 0.d0
         do n = UFS, UFS+nspec-1
            sum = sum + u(i,n)
         end do
         if (sum .ne. 0.d0) then
            fac = u(i,URHO) / sum
         else
            fac = 1.d0
         end if
         do n = UFS, UFS+nspec-1
            u(i,n) = u(i,n) * fac
         end do
      end do

      end subroutine normalize_new_species

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine reset_internal_e(u,u_l1,u_h1,lo,hi,verbose)

      use eos_module
      use network, only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UFS, UFX, small_temp, allow_negative_energy

      implicit none

      integer          :: lo(1), hi(1), verbose
      integer          :: u_l1,u_h1
      double precision :: u(u_l1:u_h1,NVAR)

      ! Local variables
      integer          :: i
      integer          :: pt_index(1)
      double precision :: Up, ke, rho_eint, eint_new, x_in(1:nspec+naux), dummy_pres

      ! Reset internal energy if negative.
      if (allow_negative_energy .eq. 0) then
         do i = lo(1),hi(1)

            Up = u(i,UMX) / u(i,URHO)
            ke = 0.5d0 * (Up**2)
   
            rho_eint = u(i,UEDEN) - u(i,URHO) * ke
   
            ! Reset (rho e) if e is greater than 0.01% of E.
            if (rho_eint .gt. 0.d0 .and. rho_eint / u(i,UEDEN) .gt. 1.d-4) then
   
                u(i,UEINT) = rho_eint
   
            ! If (e from E) < 0 or (e from E) < .0001*E but (e from e) > 0.
            else if (u(i,UEINT) .gt. 0.d0) then

               u(i,UEDEN) = u(i,UEINT) + u(i,URHO) * ke

            ! If not resetting and little e is negative ...
            else if (u(i,UEINT) .le. 0.d0) then
   
               x_in(1:nspec) = u(i,UFS:UFS+nspec-1) / u(i,URHO)
               if (naux > 0) &
                 x_in(nspec+1:nspec+naux)  = u(i,UFX:UFX+naux -1) / u(i,URHO)
   
               pt_index(1) = i
               call eos_given_RTX(eint_new, dummy_pres, u(i,URHO), small_temp, x_in, pt_index=pt_index)
               if (verbose .gt. 0) then
                  print *,'   '
                  print *,'>>> Warning: Castro_1d::reset_internal_e ',i
                  print *,'>>> ... resetting neg. e from EOS using small_temp'
                  print *,'>>> ... from ',u(i,UEINT)/u(i,URHO),' to ', eint_new
                  print *,'    '
               end if
   
               u(i,UEINT) = u(i,URHO) *  eint_new
               u(i,UEDEN) = u(i,URHO) * (eint_new + ke)
   
            end if
         enddo

      ! If (allow_negative_energy .eq. 1) then just reset (rho e) from (rho E)
      else

         do i = lo(1),hi(1)

            Up = u(i,UMX) / u(i,URHO)
            ke = 0.5d0 * (Up**2)

            u(i,UEINT) = u(i,UEDEN) - u(i,URHO) * ke

         enddo

      endif
 
    end subroutine reset_internal_e

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine get_center(center_out)

        use probdata_module, only : center

        implicit none

        double precision, intent(inout) :: center_out(1)

        center_out(1) = center(1)

      end subroutine get_center

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine set_center(center_in)

        use probdata_module, only : center

        implicit none

        double precision :: center_in(1)

        center(1) = center_in(1)

      end subroutine set_center

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine find_center(data,new_center)

        implicit none

        double precision :: data(0:2)
        double precision :: new_center(1)
        integer          :: i

        ! In 1-D it only make sense to have the center at the origin
        new_center(1) = 0.d0 

      end subroutine find_center
