! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine ca_umdrv(is_finest_level,time,&
                          lo,hi,domlo,domhi,&
                          uin,uin_l1,uin_h1,&
                          uout,uout_l1,uout_h1,&
                          ugdnv,ugdnv_l1,ugdnv_h1,&
                          src,src_l1,src_h1, &
                          grav,gv_l1,gv_h1, &
                          delta,dt,&
                          flux,flux_l1,flux_h1,&
                          area,area_l1,area_h1,&
                          dloga,dloga_l1,dloga_h1,&
                          vol,vol_l1,vol_h1,courno,verbose,&
                          mass_added,eint_added,eden_added,&
                          E_added_flux, E_added_grav)

      use meth_params_module, only : QVAR, QU, NVAR, NHYP, URHO, &
                                     use_colglaz, do_sponge, &
                                     normalize_species, use_flattening
      use advection_module, only : umeth1d, ctoprim, consup, enforce_minimum_density, &
           normalize_new_species
      use CastroCG_module, only : ctoprimcg, umeth1dcg
      use sponge_module, only : sponge

      implicit none

      integer is_finest_level
      integer lo(1),hi(1),verbose
      integer domlo(1),domhi(1)
      integer uin_l1,uin_h1
      integer uout_l1,uout_h1
      integer ugdnv_l1,ugdnv_h1
      integer flux_l1,flux_h1
      integer area_l1,area_h1
      integer dloga_l1,dloga_h1
      integer vol_l1,vol_h1
      integer src_l1,src_h1
      integer gv_l1,gv_h1
      double precision   uin(  uin_l1:  uin_h1,NVAR)
      double precision  uout( uout_l1: uout_h1,NVAR)
      double precision ugdnv(ugdnv_l1:ugdnv_h1)
      double precision   src(  src_l1:  src_h1,NVAR)
      double precision  grav(   gv_l1:   gv_h1     )
      double precision  flux( flux_l1: flux_h1,NVAR)
      double precision  area( area_l1: area_h1     )
      double precision dloga(dloga_l1:dloga_h1     )
      double precision   vol(  vol_l1: vol_h1      )
      double precision delta(1),dt,time,courno

!     Automatic arrays for workspace
      double precision, allocatable:: q(:,:)
      double precision, allocatable:: gamc(:)
      double precision, allocatable:: flatn(:)
      double precision, allocatable:: c(:)
      double precision, allocatable:: csml(:)
      double precision, allocatable:: div(:)
      double precision, allocatable:: pgdnv(:)
      double precision, allocatable:: srcQ(:,:)
      double precision, allocatable:: pdivu(:)

      double precision :: dx,E_added_flux,E_added_grav
      double precision :: mass_added, eint_added, eden_added
      integer i,ngf

      allocate(     q(uin_l1:uin_h1,QVAR))
      allocate(     c(uin_l1:uin_h1))
      allocate(  gamc(uin_l1:uin_h1))
      allocate( flatn(uin_l1:uin_h1))
      allocate(  csml(uin_l1:uin_h1))

      allocate(  srcQ(src_l1:src_h1,QVAR))

      allocate(   div(lo(1):hi(1)+1))
      allocate( pdivu(lo(1):hi(1)  ))
      allocate( pgdnv(lo(1):hi(1)+1))

      dx = delta(1)

      ngf = 1

!     Translate to primitive variables, compute sound speeds
!     Note that (q,c,gamc,csml,flatn) are all dimensioned the same
!       and set to correspond to coordinates of (lo:hi)
   
      if (use_colglaz.eq.0) then

         call ctoprim(lo,hi,uin,uin_l1,uin_h1, &
                   q,c,gamc,csml,flatn,uin_l1,uin_h1, &
                   src,srcQ,src_l1,src_h1, &
                   courno,dx,dt,NHYP,ngf)

         call umeth1d(lo,hi,domlo,domhi, &
                   q,c,gamc,csml,flatn,uin_l1,uin_h1, &
                   srcQ, src_l1, src_h1, &
                   grav, gv_l1, gv_h1, &
                   lo(1),hi(1),dx,dt, &
                   flux,flux_l1,flux_h1, &
                   pgdnv,lo(1),hi(1)+1, &
                   ugdnv,ugdnv_l1,ugdnv_h1, &
                   dloga,dloga_l1,dloga_h1)

      else

          call ctoprimcg(lo,hi,uin,uin_l1,uin_h1, &
                    q,c,gamc,csml,flatn,uin_l1,uin_h1, &
                    src,srcQ,src_l1,src_h1, &
                    courno,dx,dt,NHYP,ngf,use_flattening)

          call umeth1dcg(q,c,gamc,csml,flatn,uin_l1,uin_h1, &
                    srcQ, src_l1, src_h1, &
                    grav, gv_l1, gv_h1, &
                    lo(1),hi(1),dx,dt, &
                    flux,flux_l1,flux_h1, & 
                    pgdnv,lo(1),hi(1)+1, &
                    dloga,dloga_l1,dloga_h1)

      endif

      ! Define p*divu
      do i = lo(1), hi(1)
         pdivu(i) = 0.5d0 * &
              (pgdnv(i+1)+pgdnv(i))*(ugdnv(i+1)*area(i+1)-ugdnv(i)*area(i)) / vol(i)
      end do

      ! Define divu on surroundingNodes(lo,hi)
      do i = lo(1),hi(1)+1
         div(i) = (q(i,QU)-q(i-1,QU)) / dx
      enddo

!     Conservative update
      call consup(uin,uin_l1,uin_h1, &
           uout,uout_l1,uout_h1, &
           pgdnv,lo(1),hi(1)+1, &
           src , src_l1, src_h1, &
           grav,  gv_l1,  gv_h1, &
           flux,flux_l1,flux_h1, &
           area,area_l1,area_h1, &
           vol , vol_l1, vol_h1, &
           div ,pdivu,lo,hi,dx,dt)

      ! Enforce the density >= small_dens.
      call enforce_minimum_density(uin,uin_l1,uin_h1,uout,uout_l1,uout_h1,lo,hi,&
                                   mass_added,eint_added,eden_added,verbose)

      ! Enforce that the species >= 0
      call ca_enforce_nonnegative_species(uout,uout_l1,uout_h1,lo,hi)

      ! Normalize the species
      if (normalize_species .eq. 1) &
         call normalize_new_species(uout,uout_l1,uout_h1,lo,hi)

      if (do_sponge .eq. 1) &
           call sponge(uout,uout_l1,uout_h1,lo,hi,time,dt,dx,domlo,domhi)

      deallocate(q,c,gamc,flatn,csml,srcQ,div,pdivu,pgdnv)

!     if ( (is_finest_level      .eq. 1) .and. &
!          (lo(1) .eq. 0               ) ) then
!        print *,'CEN0  ',time, uin(0,URHO)
!        print *,'CEN1  ',time, uin(1,URHO)
!        print *,'CEN2  ',time, uin(2,URHO)
!     end if

      end subroutine ca_umdrv


! ::: 
! ::: ------------------------------------------------------------------
! ::: 

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

      subroutine ca_enforce_nonnegative_species(uout,uout_l1,uout_h1,lo,hi)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(1), hi(1)
      integer          :: uout_l1,uout_h1
      double precision :: uout(uout_l1:uout_h1,NVAR)

      ! Local variables
      integer          :: i,n
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

        ! In 1-D it only make sense to have the center at the origin
        new_center(1) = 0.d0 

      end subroutine find_center

