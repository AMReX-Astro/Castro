! ::
! :: ----------------------------------------------------------
! ::
       subroutine pm_add_to_grav(point_mass,grav,grav_l1,grav_h1,problo,dx)

       use fundamental_constants_module, only : Gconst
       use probdata_module             , only : center

       implicit none
       integer         , intent(in   ) :: grav_l1,grav_h1
       double precision, intent(in   ) :: point_mass
       double precision, intent(inout) :: grav(grav_l1  :grav_h1)
       double precision, intent(in   ) :: problo(1),dx(1)

       integer          :: i
       double precision :: x

       ! We must be in spherical coordinates with center(1) = 0.
       !    or this doesn't make sense

!      This computes radial gravity due to a point mass at the origin.
       do i = grav_l1, grav_h1

          x = problo(1) + (dble(i)+0.5d0) * dx(1)

          grav(i) = grav(i) - Gconst * point_mass / (x*x)

       end do

       end subroutine pm_add_to_grav
! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine pm_compute_delta_mass(delta_mass,lo,hi,&
             uin,  uin_l1,  uin_h1, &
            uout, uout_l1, uout_h1, &
             vol,  vol_l1,  vol_h1, &
           problo,dx,time,dt)

      use meth_params_module, only : NVAR, URHO, UMX

      implicit none

      integer :: lo(1),hi(1)
      integer ::   uin_l1,  uin_h1
      integer ::  uout_l1, uout_h1
      integer ::   vol_l1,  vol_h1

      double precision   :: delta_mass
      double precision   ::   uin( uin_l1: uin_h1, NVAR)
      double precision   ::  uout(uout_l1:uout_h1, NVAR)
      double precision   ::   vol( vol_l1: vol_h1)
      double precision   :: problo(1),dx(1),time,dt

      integer, parameter :: box_size = 2
      integer            :: ii

      ! Make sure we only count contributions from the grid nearest the origin
      if (lo(1) .eq. 0) then
          do ii = 0, box_size-1
             delta_mass = delta_mass + vol(ii) * (uout(ii,URHO)-uin(ii,URHO))
          end do
          delta_mass = max(0.d0, delta_mass)
      end if

      end subroutine pm_compute_delta_mass

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine pm_fix_solution(lo,hi,&
             uin,  uin_l1,  uin_h1, &
            uout, uout_l1, uout_h1, &
           problo,dx,time,dt)

      use meth_params_module, only : NVAR, UMX

      implicit none

      integer :: lo(1),hi(1)
      integer ::   uin_l1,  uin_h1
      integer ::  uout_l1, uout_h1

      double precision   ::   uin(  uin_l1:uin_h1,  NVAR)
      double precision   ::  uout( uout_l1:uout_h1, NVAR)
      double precision   :: problo(1),dx(1),time,dt

      integer, parameter :: box_size = 2
      integer            :: ii

      if (box_size .ne. 2) then
         print *,"We assume box_size = 2 below!" 
         stop
      end if

      do ii = 0,box_size-1
         uout(ii,:) = uout(box_size,:)
      end do

      ! Linear interpolation between the point (box_size) and the origin.
      uout(0,UMX) = 0.2d0 * uout(box_size,UMX)
      uout(1,UMX) = 0.6d0 * uout(box_size,UMX)

      end subroutine pm_fix_solution
