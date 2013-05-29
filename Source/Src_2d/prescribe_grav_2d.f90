!     This is an example of how to specify a radial profile
!     Note that r_c and rho_c must be specified in probdata_module

!     function ca_prescribe_grav_gravityprofile(r) result(g)
!     use fundamental_constants_module, only : Gconst
!     use probdata_module
!     implicit none
!     double precision, intent(in) :: r
!     double precision             :: g
!         ! put function for g(r) in here
!         g =-rho_c*Gconst*2*M_PI*(1/sqrt(1+(r/r_c)**2)-atanh(1/sqrt(1+(r/r_c)**2)))
!     end function ca_prescribe_grav_gravityprofile

      subroutine ca_prescribe_grav (lo,hi,dx, &
                                    grav,g_l1,g_l2,g_h1,g_h2, &
                                    problo)

      use fundamental_constants_module, only : Gconst
      use probdata_module
      use bl_constants_module, only : M_PI

      implicit none
      integer          :: g_l1,g_l2,g_h1,g_h2
      integer          :: lo(2),hi(2)
      double precision :: grav(g_l1:g_h1,g_l2:g_h2,2)
      double precision :: dx(2)
      double precision :: problo(2)
  
      ! Local variables
!     integer          :: i,j
!     double precision :: ca_prescribe_grav_gravityprofile
!     double precision :: x,y
!     double precision :: r,maggrav

!     This is an example of how to use the radial profile above
!     do j = lo(2), hi(2)
!        y = problo(2) + (dble(j)+0.50d0) * dx(2) - center(2)

!        do i = lo(1), hi(1)
!           x = problo(1) + (dble(i)+0.50d0) * dx(1) - center(1)

!           r = sqrt(x**2+y**2)

!           maggrav = ca_prescribe_grav_gravityprofile(r)

!           grav(i,j,1) = maggrav*(x/r)
!           grav(i,j,2) = maggrav*(y/r)

!        enddo
!     enddo

      grav = 0.d0

      end subroutine ca_prescribe_grav



