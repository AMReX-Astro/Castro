! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_umdrv(is_finest_level,time,lo,hi,domlo,domhi, &
                          uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                          uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                          ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                          ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
                          src,src_l1,src_l2,src_h1,src_h2, &
                          grav,gv_lo,gv_hi, &
                          rot,rt_lo,rt_hi, &
                          delta,dt, &
                          flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                          flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                          area1,area1_l1,area1_l2,area1_h1,area1_h2, &
                          area2,area2_l1,area2_l2,area2_h1,area2_h2, &
                          dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                          vol,vol_l1,vol_l2,vol_h1,vol_h2,&
                          courno,verbose,mass_added,eint_added,eden_added,&
                          xmom_added_flux, ymom_added_flux, &
                          xmom_added_sponge, ymom_added_sponge, &
                          E_added_flux,E_added_sponge)

      use meth_params_module, only : QVAR, NVAR, NHYP, &
                                     do_sponge, normalize_species
      use advection_module, only : umeth2d, ctoprim, divu, consup, enforce_minimum_density, &
           normalize_new_species
      use sponge_module, only : sponge

      implicit none

      integer is_finest_level
      integer lo(2),hi(2),verbose
      integer domlo(2),domhi(2)
      integer uin_l1,uin_l2,uin_h1,uin_h2
      integer uout_l1,uout_l2,uout_h1,uout_h2
      integer ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2
      integer ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2
      integer flux1_l1,flux1_l2,flux1_h1,flux1_h2
      integer flux2_l1,flux2_l2,flux2_h1,flux2_h2
      integer area1_l1,area1_l2,area1_h1,area1_h2
      integer area2_l1,area2_l2,area2_h1,area2_h2
      integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
      integer vol_l1,vol_l2,vol_h1,vol_h2
      integer src_l1,src_l2,src_h1,src_h2
      integer gv_lo(3),gv_hi(3)
      integer rt_lo(3),rt_hi(3)

      double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
      double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
      double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)
      double precision src(src_l1:src_h1,src_l2:src_h2,NVAR)
      double precision grav(gv_lo(1):gv_hi(1),gv_lo(2):gv_hi(2),gv_lo(3):gv_hi(3),3)
      double precision  rot(rt_lo(1):rt_hi(1),rt_lo(2):rt_hi(2),rt_lo(3):rt_hi(3),3)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
      double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
      double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
      double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
      double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
      double precision delta(2),dt,time,courno
      double precision E_added_flux,E_added_sponge
      double precision xmom_added_flux, ymom_added_flux
      double precision xmom_added_sponge, ymom_added_sponge
      double precision mass_added,eint_added,eden_added

!     Automatic arrays for workspace
      double precision, allocatable:: q(:,:,:)
      double precision, allocatable:: gamc(:,:)
      double precision, allocatable:: flatn(:,:)
      double precision, allocatable:: c(:,:)
      double precision, allocatable:: csml(:,:)
      double precision, allocatable:: div(:,:)
      double precision, allocatable:: pgdx(:,:)
      double precision, allocatable:: pgdy(:,:)
      double precision, allocatable:: srcQ(:,:,:)
      double precision, allocatable:: pdivu(:,:)

      integer ngq,ngf
!     integer i_c,j_c

      double precision dx,dy

      integer q_l1, q_l2, q_h1, q_h2

      ngq = NHYP
      ngf = 1

      q_l1 = lo(1)-NHYP
      q_l2 = lo(2)-NHYP
      q_h1 = hi(1)+NHYP
      q_h2 = hi(2)+NHYP

      allocate(     q(q_l1:q_h1,q_l2:q_h2,QVAR))
      allocate(  gamc(q_l1:q_h1,q_l2:q_h2))
      allocate( flatn(q_l1:q_h1,q_l2:q_h2))
      allocate(     c(q_l1:q_h1,q_l2:q_h2))
      allocate(  csml(q_l1:q_h1,q_l2:q_h2))

      allocate(  srcQ(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,QVAR))

      allocate(   div(lo(1)  :hi(1)+1,lo(2)  :hi(2)+1))
      allocate( pdivu(lo(1)  :hi(1)  ,lo(2)  :hi(2)))
      allocate(  pgdx(lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1))
      allocate(  pgdy(lo(1)-1:hi(1)+1,lo(2)  :hi(2)+1))

      dx = delta(1)
      dy = delta(2)

!     Translate to primitive variables, compute sound speeds
!     Note that (q,c,gamc,csml,flatn) are all dimensioned the same
!       and set to correspond to coordinates of (lo:hi)
      call ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                   q,c,gamc,csml,flatn,q_l1,q_l2,q_h1,q_h2, &
                   src,src_l1,src_l2,src_h1,src_h2, &
                   srcQ,lo(1)-1,lo(2)-1,hi(1)+1,hi(2)+1, &
                   courno,dx,dy,dt,ngq,ngf)

!     Compute hyperbolic fluxes using unsplit Godunov
      call umeth2d(q,c,gamc,csml,flatn,q_l1,q_l2,q_h1,q_h2, &
                   srcQ, lo(1)-1,lo(2)-1,hi(1)+1,hi(2)+1, &
                   grav,gv_lo(1),gv_lo(2),gv_lo(3),gv_hi(1),gv_hi(2),gv_hi(3), &
                   rot, rt_lo(1),rt_lo(2),rt_lo(3),rt_hi(1),rt_hi(2),rt_hi(3), &
                   lo(1),lo(2),hi(1),hi(2),dx,dy,dt, &
                   flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                   flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                   pgdx, lo(1), lo(2)-1, hi(1)+1, hi(2)+1, &
                   pgdy, lo(1)-1, lo(2), hi(1)+1, hi(2)+1, &
                   ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                   ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
                   area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                   area2, area2_l1, area2_l2, area2_h1, area2_h2, &
                   pdivu, vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                   dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                   domlo, domhi)

      ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
      ! this is used for the artifical viscosity
      call divu(lo,hi,q,q_l1,q_l2,q_h1,q_h2, &
                delta,div,lo(1),lo(2),hi(1)+1,hi(2)+1)

!     Conservative update
      call consup(uin,    uin_l1,  uin_l2,  uin_h1,  uin_h2, &
                  uout,  uout_l1, uout_l2, uout_h1, uout_h2, &
                  pgdx,  lo(1), lo(2)-1, hi(1)+1, hi(2)+1, &
                  pgdy,lo(1)-1,   lo(2), hi(1)+1, hi(2)+1, &
                  src,    src_l1,  src_l2,  src_h1,  src_h2, &
                  flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                  flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                  area1,area1_l1,area1_l2,area1_h1,area1_h2, &
                  area2,area2_l1,area2_l2,area2_h1,area2_h2, &
                  vol,    vol_l1,  vol_l2,  vol_h1,  vol_h2, &
                  div,pdivu,lo,hi,dx,dy,dt,E_added_flux,&
                  xmom_added_flux,ymom_added_flux)

      ! Enforce the density >= small_dens.
      call enforce_minimum_density( uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                                   uout,uout_l1,uout_l2,uout_h1,uout_h2,&
                                   lo,hi,mass_added,eint_added,eden_added,verbose)

      ! Enforce the species >= 0
      call ca_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi)

      ! Normalize the species 
      if (normalize_species .eq. 1) &
         call normalize_new_species(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi)

      if (do_sponge .eq. 1) &
           call sponge(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi, &
                       time,dt, &
                       dx,dy,domlo,domhi, &
                       E_added_sponge,xmom_added_sponge,ymom_added_sponge)

      deallocate(q,gamc,flatn,c,csml,div,pgdx,pgdy,srcQ,pdivu)

      end subroutine ca_umdrv

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_check_initial_species(lo,hi,&
                             state,state_l1,state_l2,state_h1,state_h2)

      use network           , only : nspec
      use meth_params_module, only : NVAR, URHO, UFS
      use bl_constants_module

      implicit none

      integer          :: lo(2), hi(2)
      integer          :: state_l1,state_l2,state_h1,state_h2
      double precision :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

      ! Local variables
      integer          :: i,j,n
      double precision :: sum

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)

         sum = ZERO
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
! ::  nvar	 => number of components in arrays
! ::  fine       => fine grid data
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

      use bl_constants_module

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

      do n = 1, nvar
 
!           Set coarse grid to zero on overlap
         do jc = lo(2), hi(2)
            do ic = lo(1), hi(1)
               crse(ic,jc,n) = ZERO
            enddo
         enddo

!           Sum fine data
         do joff = 0, lrat(2)-1
            do jc = lo(2), hi(2)
               j = jc*lrat(2) + joff
               do ioff = 0, lrat(1)-1
                  do ic = lo(1), hi(1)
                     i = ic*lrat(1) + ioff
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

      end subroutine ca_avgdown

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_compute_avgstate (lo,hi,dx,dr,nc,&
                                      state,s_l1,s_l2,s_h1,s_h2,radial_state, &
                                      vol,v_l1,v_l2,v_h1,v_h2,radial_vol, &
                                      problo,numpts_1d)

      use meth_params_module, only: URHO, UMX, UMY
      use prob_params_module, only: center
      use bl_constants_module

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
         y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
         do i = lo(1), hi(1)
            x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
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

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS
      use bl_constants_module

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
           if (uout(i,j,n) .lt. ZERO) then
              x = uout(i,j,n)/uout(i,j,URHO)
              if (x .gt. eps) then
                 uout(i,j,n) = ZERO
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

              if (uout(i,j,n) .lt. ZERO) then

                 x = uout(i,j,n)/uout(i,j,URHO)

                 ! Here we only print the bigger negative values
                 if (x .lt. -1.d-2) then
                    print *,'At cell (i,j) = ',i,j
                    print *,'... Fixing negative species ',n-UFS+1           ,' with X = ',x
                    print *,'...   from dominant species ',int_dom_spec-UFS+1,' with X = ',&
                             uout(i,j,int_dom_spec) / uout(i,j,URHO)
                 end if

                 ! Take enough from the dominant species to fill the negative one.
                 uout(i,j,int_dom_spec) = uout(i,j,int_dom_spec) + uout(i,j,n)
   
                 ! Test that we didn't make the dominant species negative
                 if (uout(i,j,int_dom_spec) .lt. ZERO) then 
                    print *,'Just made dominant species negative ',int_dom_spec-UFS+1,' at ',i,j
                    print *,'... We were fixing species ',n-UFS+1,' which had value ',x
                    print *,'... Dominant species became ',uout(i,j,int_dom_spec) / uout(i,j,URHO)
                    call bl_error("Error:: Castro_2d.f90 :: ca_enforce_nonnegative_species")
                 end if

                 ! Now the negative species to zero
                 uout(i,j,n) = ZERO

              end if

           enddo
         end if

      enddo
      enddo

      end subroutine ca_enforce_nonnegative_species

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine get_center(center_out)

        use prob_params_module, only : center

        implicit none

        double precision :: center_out(2)

        center_out(1:2) = center(1:2)

      end subroutine get_center

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine set_center(center_in)

        use prob_params_module, only : center

        implicit none

        double precision :: center_in(2)

        center(1:2) = center_in(1:2)

      end subroutine set_center

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine find_center(data,new_center,icen,dx,problo)

        use bl_constants_module

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
        new_center(1) = problo(1) +  (icen(1)+HALF) * dx(1)
        new_center(2) = problo(2) +  (icen(2)+HALF) * dx(2)
   
        ! Fit parabola y = a x^2  + b x + c through three points
        ! a = 1/2 ( y_1 + y_-1)
        ! b = 1/2 ( y_1 - y_-1)
        ! x_vertex = -b / 2a

        ! ... in x-direction
        a = HALF * (data(1,0) + data(-1,0)) - data(0,0)
        b = HALF * (data(1,0) - data(-1,0)) - data(0,0)
        x = -b / (TWO*a)
        new_center(1) = new_center(1) +  x*dx(1)

        ! ... in y-direction
        a = HALF * (data(0,1) + data(0,-1)) - data(0,0)
        b = HALF * (data(0,1) - data(0,-1)) - data(0,0)
        y = -b / (TWO*a)
        new_center(2) = new_center(2) +  y*dx(2)

      end subroutine find_center

