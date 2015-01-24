! ::
! :: ----------------------------------------------------------
! ::
       subroutine pm_add_to_grav(point_mass,&
                                 grav,grav_l1,grav_l2,grav_h1,grav_h2,&
                                 problo,dx,lo,hi)

       use fundamental_constants_module, only : Gconst
       use probdata_module             , only : center
       use bl_constants_module

       implicit none

       integer         , intent(in   ) :: lo(2), hi(2)
       integer         , intent(in   ) :: grav_l1,grav_l2,grav_h1,grav_h2
       double precision, intent(in   ) :: point_mass
       double precision, intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2,2)
       double precision, intent(in   ) :: problo(2),dx(2)

       integer          :: i,j
       double precision :: x,y,r,rsq,radial_force

       ! We must be in r-z coordinates with problo(1) = 0. and center(1) = 0.
       !    or this doesn't make sense

!      This computes radial gravity due to a point mass at center().
       do j = lo(2), hi(2)
          y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             x = problo(1) + abs(dble(i)+HALF) * dx(1) - center(1)

             rsq = x*x + y*y
             radial_force = -Gconst * point_mass / rsq

             r = sqrt(rsq)
             grav(i,j,1) = grav(i,j,1) + radial_force * (x/r)
             grav(i,j,2) = grav(i,j,2) + radial_force * (y/r)

          end do
       end do

       end subroutine pm_add_to_grav

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine pm_compute_delta_mass(delta_mass,lo,hi,&
             uin,  uin_l1,  uin_l2,  uin_h1,  uin_h2, &
            uout, uout_l1, uout_l2, uout_h1, uout_h2, &
             vol,  vol_l1,  vol_l2,  vol_h1,  vol_h2, &
           problo,dx,time,dt)

      use meth_params_module, only : NVAR, URHO
      use probdata_module   , only : center
      use bl_constants_module

      implicit none

      integer :: lo(2),hi(2)
      integer :: uin_l1,uin_l2,uin_h1,uin_h2
      integer :: uout_l1,uout_l2,uout_h1,uout_h2
      integer ::   vol_l1,  vol_l2,  vol_h1,  vol_h2

      double precision   ::   delta_mass
      double precision   ::   uin(  uin_l1:uin_h1,    uin_l2:uin_h2,  NVAR)
      double precision   ::  uout( uout_l1:uout_h1,  uout_l2:uout_h2, NVAR)
      double precision   ::   vol(  vol_l1:  vol_h1,  vol_l2:  vol_h2)
      double precision   :: problo(2),dx(2),time,dt

      double precision   :: eps
      integer            :: icen,istart,iend
      integer            :: jcen,jstart,jend
      integer            :: ii,jj,n
      integer, parameter :: box_size = 2

      ! We must be in r-z coordinates with problo(1) = 0. and center(1) = 0.
      !    or this doesn't make sense

      ! This is just a small number to keep precision issues from making
      !   icen,jcen,kcen one cell too low.
      eps = 1.d-8
 
      ! This should be the cell whose lower left corner is at "center"
      icen = floor( (center(1)-problo(1))/dx(1) + eps)
      jcen = floor( (center(2)-problo(2))/dx(2) + eps)
 
      ! Make sure we only count contributions from this grid
      istart = max(icen-box_size, lo(1))
      jstart = max(jcen-box_size, lo(2))
 
      iend = min(icen+box_size-1, hi(1))
      jend = min(jcen+box_size-1, hi(2))
 
      do jj = jstart,jend
      do ii = istart,iend
        delta_mass = delta_mass + vol(ii,jj) * (uout(ii,jj,URHO)-uin(ii,jj,URHO))
      end do
      end do

      end subroutine pm_compute_delta_mass

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine pm_fix_solution(lo,hi,&
             uin,  uin_l1,  uin_l2,  uin_h1,  uin_h2, &
            uout, uout_l1, uout_l2, uout_h1, uout_h2, &
           problo,dx,time,dt)

      use meth_params_module, only : NVAR
      use probdata_module   , only : center

      implicit none

      integer :: lo(2),hi(2)
      integer :: uin_l1,uin_l2,uin_h1,uin_h2
      integer :: uout_l1,uout_l2,uout_h1,uout_h2

      double precision   ::   uin(  uin_l1:uin_h1,    uin_l2:uin_h2,  NVAR)
      double precision   ::  uout( uout_l1:uout_h1,  uout_l2:uout_h2, NVAR)
      double precision   :: problo(2),dx(2),time,dt

      double precision   :: eps
      integer            :: ii,icen,istart,iend
      integer            :: jj,jcen,jstart,jend
      integer, parameter :: box_size = 2

      ! This is just a small number to keep precision issues from making
      !   icen,jcen,kcen one cell too low.
      eps = 1.d-8
 
      ! This should be the cell whose lower left corner is at "center"
      icen = floor( (center(1)-problo(1))/dx(1) + eps)
      jcen = floor( (center(2)-problo(2))/dx(2) + eps)
 
      ! Make sure we only count contributions from this grid
      istart = max(icen-box_size, lo(1))
      jstart = max(jcen-box_size, lo(2))
 
      iend = min(icen+box_size-1, hi(1))
      jend = min(jcen+box_size-1, hi(2))

      do jj = jstart,jend
      do ii = istart,iend
         uout(ii,jj,: ) = uin(ii,jj,:)
      end do
      end do

      end subroutine pm_fix_solution
