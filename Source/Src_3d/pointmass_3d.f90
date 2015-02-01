! ::
! :: ----------------------------------------------------------
! ::
       subroutine pm_add_to_grav(point_mass,&
                                 grav,grav_l1,grav_l2,grav_l3, &
                                      grav_h1,grav_h2,grav_h3, &
                                 problo,dx,lo,hi)

       use fundamental_constants_module, only : Gconst
       use prob_params_module          , only : center

       implicit none
       integer         , intent(in   ) :: lo(3), hi(3)
       integer         , intent(in   ) :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
       double precision, intent(in   ) :: point_mass
       double precision, intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3,3)
       double precision, intent(in   ) :: problo(3),dx(3)

       integer          :: i,j,k
       double precision :: x,y,z,r,rsq,radial_force

!      This computes radial gravity due to a point mass at center().
       do k = lo(3), hi(3)
          z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
          do j = lo(2), hi(2)
             y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
             do i = lo(1), hi(1)
                x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)

                rsq = x*x + y*y + z*z
                radial_force = -Gconst * point_mass / rsq

                r = sqrt(rsq)
                grav(i,j,k,1) = grav(i,j,k,1) + radial_force * (x/r)
                grav(i,j,k,2) = grav(i,j,k,2) + radial_force * (y/r)
                grav(i,j,k,3) = grav(i,j,k,3) + radial_force * (y/r)

             end do
          end do
       end do

       end subroutine pm_add_to_grav

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine pm_compute_delta_mass(delta_mass,lo,hi,&
             uin,  uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3, &
            uout, uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3, &
             vol,  vol_l1,  vol_l2,  vol_l3,  vol_h1,  vol_h2,  vol_h3, &
           problo,dx,time,dt) 

      use meth_params_module, only : NVAR, URHO
      use prob_params_module, only : center

      implicit none

      integer :: lo(3),hi(3)
      integer ::   uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
      integer ::  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      integer ::   vol_l1,  vol_l2,  vol_l3,  vol_h1,  vol_h2,  vol_h3

      double precision   ::   delta_mass
      double precision   ::   uin(  uin_l1:uin_h1,    uin_l2:uin_h2,  &
                                    uin_l3:uin_h3,NVAR)
      double precision   ::  uout( uout_l1:uout_h1,  uout_l2:uout_h2, &
                                    uin_l3:uin_h3,NVAR)
      double precision   ::   vol(  vol_l1:  vol_h1,  vol_l2:  vol_h2, &
                                    vol_l3:  vol_h3)
      double precision   :: problo(3),dx(3),time,dt

      double precision   :: eps
      integer            :: ii,icen,istart,iend 
      integer            :: jj,jcen,jstart,jend 
      integer            :: kk,kcen,kstart,kend 
      integer, parameter :: box_size = 2

      ! This is just a small number to keep precision issues from making
      !   icen,jcen,kcen one cell too low.
      eps = 1.d-8

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

      delta_mass = max(ZERO, delta_mass)

      end subroutine pm_compute_delta_mass

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine pm_fix_solution(lo,hi,&
             uin,  uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3, &
            uout, uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3, &
           problo,dx,time,dt) 

      use meth_params_module, only : NVAR
      use prob_params_module, only : center

      implicit none

      integer :: lo(3),hi(3)
      integer ::   uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
      integer ::  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3

      double precision   ::   uin(  uin_l1:uin_h1,    uin_l2:uin_h2,  &
                                    uin_l3:uin_h3,NVAR)
      double precision   ::  uout( uout_l1:uout_h1,  uout_l2:uout_h2, &
                                    uin_l3:uin_h3,NVAR)
      double precision   :: problo(3),dx(3),time,dt

      double precision   :: eps
      integer            :: ii,icen,istart,iend
      integer            :: jj,jcen,jstart,jend
      integer            :: kk,kcen,kstart,kend
      integer, parameter :: box_size = 2

      ! This is just a small number to keep precision issues from making
      !   icen,jcen,kcen one cell too low.
      eps = 1.d-8

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
