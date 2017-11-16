! ::
! :: ----------------------------------------------------------
! ::
       subroutine pm_add_to_grav(point_mass,&
                                 phi,phi_lo,phi_hi, &
                                 grav,grav_lo,grav_hi, &
                                 problo,dx,lo,hi) &
                                 bind(C, name="pm_add_to_grav")

       use bl_constants_module         , only : HALF
       use fundamental_constants_module, only : Gconst
       use prob_params_module          , only : center

       use amrex_fort_module, only : rt => amrex_real
       implicit none
       integer         , intent(in   ) :: lo(3), hi(3)
       integer         , intent(in   ) :: grav_lo(3), grav_hi(3)
       integer         , intent(in   ) :: phi_lo(3), phi_hi(3)
       real(rt)        , intent(in   ) :: point_mass
       real(rt)        , intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
       real(rt)        , intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3),3)
       real(rt)        , intent(in   ) :: problo(3),dx(3)

       integer          :: i,j,k
       real(rt)         :: x,y,z,rsq,radial_force,rinv

!      This computes radial gravity due to a point mass at center().
       do k = lo(3), hi(3)
          z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
          do j = lo(2), hi(2)
             y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
             do i = lo(1), hi(1)
                x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)

                rsq = x*x + y*y + z*z
                radial_force = -Gconst * point_mass / rsq

                rinv = 1.e0_rt/sqrt(rsq)

                ! Note that grav may have more ghost zones than
                ! phi, so we need to check that we're doing
                ! valid indexing here.

                if ( i .ge. phi_lo(1) .and. i .le. phi_hi(1) .and. &
                     j .ge. phi_lo(2) .and. j .le. phi_hi(2) .and. &
                     k .ge. phi_lo(3) .and. k .le. phi_hi(3) ) then

                   phi(i,j,k) = phi(i,j,k) - Gconst * point_mass * rinv

                endif

                grav(i,j,k,1) = grav(i,j,k,1) + radial_force * (x*rinv)
                grav(i,j,k,2) = grav(i,j,k,2) + radial_force * (y*rinv)
                grav(i,j,k,3) = grav(i,j,k,3) + radial_force * (z*rinv)

             end do
          end do
       end do

       end subroutine pm_add_to_grav

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine pm_compute_delta_mass(delta_mass,lo,hi,&
                                       uin, uin_lo, uin_hi, &
                                       uout, uout_lo, uout_hi, &
                                       vol,  vol_lo, vol_hi, &
                                       problo,dx,time,dt) &
                                       bind(C, name="pm_compute_delta_mass")

      use meth_params_module, only : NVAR, URHO
      use prob_params_module, only : center

      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer, intent(in) :: lo(3),      hi(3)
      integer, intent(in) :: uin_lo(3),  uin_hi(3)
      integer, intent(in) :: uout_lo(3), uout_hi(3)
      integer, intent(in) :: vol_lo(3),  vol_hi(3)

      real(rt), intent(inout) ::   delta_mass
      real(rt), intent(in) ::   uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),  &
                                  uin_lo(3):uin_hi(3),NVAR)
      real(rt), intent(in) ::  uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2), &
                                  uout_lo(3):uout_hi(3),NVAR)
      real(rt), intent(in) ::   vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2), &
                                  vol_lo(3):vol_hi(3))
      real(rt), intent(in) :: problo(3),dx(3),time,dt

      real(rt)           :: eps
      integer            :: ii,icen,istart,iend
      integer            :: jj,jcen,jstart,jend
      integer            :: kk,kcen,kstart,kend
      integer, parameter :: box_size = 2

      ! This is just a small number to keep precision issues from making
      !   icen,jcen,kcen one cell too low.
      eps = 1.e-8_rt

      ! This should be the cell whose lower left corner is at "center"
      icen = floor( (center(1)-problo(1))/dx(1) + eps)
      jcen = floor( (center(2)-problo(2))/dx(2) + eps)
      kcen = floor( (center(3)-problo(3))/dx(3) + eps)

      ! Make sure we only count contributions from this grid
      istart = max(icen-box_size, lo(1))
      jstart = max(jcen-box_size, lo(2))
      kstart = max(kcen-box_size, lo(3))

      iend = min(icen+box_size-1, hi(1))
      jend = min(jcen+box_size-1, hi(2))
      kend = min(kcen+box_size-1, hi(3))

      do kk = kstart,kend
         do jj = jstart,jend
            do ii = istart,iend
               delta_mass = delta_mass + vol(ii,jj,kk) * (uout(ii,jj,kk,URHO)-uin(ii,jj,kk,URHO))
            end do
         end do
      end do

      end subroutine pm_compute_delta_mass

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine pm_fix_solution(lo,hi, &
                                 uin,uin_lo,uin_hi, &
                                 uout,uout_lo,uout_hi, &
                                 problo,dx,time,dt) &
                                 bind(C, name="pm_fix_solution")

      use meth_params_module, only : NVAR
      use prob_params_module, only : center

      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer, intent(in) :: lo(3),      hi(3)
      integer, intent(in) :: uin_lo(3),  uin_hi(3)
      integer, intent(in) :: uout_lo(3), uout_hi(3)

      real(rt), intent(in)::   uin(  uin_lo(1):uin_hi(1),    uin_lo(2):uin_hi(2),  &
                                    uin_lo(3):uin_hi(3),NVAR)
      real(rt), intent(inout) ::  uout( uout_lo(1):uout_hi(1),  uout_lo(2):uout_hi(2), &
                                    uout_lo(3):uout_hi(3),NVAR)
      real(rt), intent(in) :: problo(3),dx(3),time,dt

      real(rt)           :: eps
      integer            :: ii,icen,istart,iend
      integer            :: jj,jcen,jstart,jend
      integer            :: kk,kcen,kstart,kend
      integer, parameter :: box_size = 2

      ! This is just a small number to keep precision issues from making
      !   icen,jcen,kcen one cell too low.
      eps = 1.e-8_rt

      ! This should be the cell whose lower left corner is at "center"
      icen = floor( (center(1)-problo(1))/dx(1) + eps)
      jcen = floor( (center(2)-problo(2))/dx(2) + eps)
      kcen = floor( (center(3)-problo(3))/dx(3) + eps)

      ! Make sure we only count contributions from this grid
      istart = max(icen-box_size, lo(1))
      jstart = max(jcen-box_size, lo(2))
      kstart = max(kcen-box_size, lo(3))

      iend = min(icen+box_size-1, hi(1))
      jend = min(jcen+box_size-1, hi(2))
      kend = min(kcen+box_size-1, hi(3))

      do kk = kstart,kend
         do jj = jstart,jend
            do ii = istart,iend
               uout(ii,jj,kk,:) = uin(ii,jj,kk,:)
            end do
         end do
      end do

      end subroutine pm_fix_solution
