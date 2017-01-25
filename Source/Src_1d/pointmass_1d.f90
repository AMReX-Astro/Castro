! ::
! :: ----------------------------------------------------------
! ::
       subroutine pm_add_to_grav(point_mass,grav,grav_l1,grav_h1,problo,dx,lo,hi) &
            bind(C, name="pm_add_to_grav")

       use fundamental_constants_module, only : Gconst
       use bl_constants_module

       use bl_fort_module, only : rt => c_real
       implicit none

       integer         , intent(in   ) :: lo(1), hi(1)
       integer         , intent(in   ) :: grav_l1,grav_h1
       real(rt)        , intent(in   ) :: point_mass
       real(rt)        , intent(inout) :: grav(grav_l1  :grav_h1)
       real(rt)        , intent(in   ) :: problo(1),dx(1)

       integer          :: i
       real(rt)         :: x

       ! We must be in spherical coordinates with center(1) = 0.
       !    or this doesn't make sense

!      This computes radial gravity due to a point mass at the origin.
       do i = lo(1), hi(1)

          x = problo(1) + (dble(i)+HALF) * dx(1)

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
           problo,dx,time,dt) bind(C, name="pm_compute_delta_mass")

      use meth_params_module, only : NVAR, URHO, UMX
      use bl_constants_module

      use bl_fort_module, only : rt => c_real
      implicit none

      integer :: lo(1),hi(1)
      integer ::   uin_l1,  uin_h1
      integer ::  uout_l1, uout_h1
      integer ::   vol_l1,  vol_h1

      real(rt)           :: delta_mass
      real(rt)           ::   uin( uin_l1: uin_h1, NVAR)
      real(rt)           ::  uout(uout_l1:uout_h1, NVAR)
      real(rt)           ::   vol( vol_l1: vol_h1)
      real(rt)           :: problo(1),dx(1),time,dt

      integer, parameter :: box_size = 2
      integer            :: ii

      ! Make sure we only count contributions from the grid nearest the origin
      if (lo(1) .eq. 0) then
          do ii = 0, box_size-1
             delta_mass = delta_mass + vol(ii) * (uout(ii,URHO)-uin(ii,URHO))
          end do
      end if

      end subroutine pm_compute_delta_mass

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine pm_fix_solution(lo,hi,&
             uin,  uin_l1,  uin_h1, &
            uout, uout_l1, uout_h1, &
           problo,dx,time,dt) bind(C, name="pm_fix_solution")

      use meth_params_module, only : NVAR, UMX

      use bl_fort_module, only : rt => c_real
      implicit none

      integer :: lo(1),hi(1)
      integer ::   uin_l1,  uin_h1
      integer ::  uout_l1, uout_h1

      real(rt)           ::   uin(  uin_l1:uin_h1,  NVAR)
      real(rt)           ::  uout( uout_l1:uout_h1, NVAR)
      real(rt)           :: problo(1),dx(1),time,dt

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
      uout(0,UMX) = 0.2e0_rt * uout(box_size,UMX)
      uout(1,UMX) = 0.6e0_rt * uout(box_size,UMX)

      end subroutine pm_fix_solution
