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
                                    grav,g_l1,g_l2,g_l3,g_h1,g_h2,g_h3,&
                                    problo)

      use fundamental_constants_module, only : Gconst
      use probdata_module
      use bl_constants_module

      implicit none
      integer          :: g_l1,g_l2,g_l3,g_h1,g_h2,g_h3
      integer          :: lo(3),hi(3)
      double precision :: grav(g_l1:g_h1,g_l2:g_h2,g_l3:g_h3,3)
      double precision :: dx(3)
      double precision :: problo(3)

      ! Local variables
!     integer          :: i,j,k
!     double precision :: ca_prescribe_grav_gravityprofile
!     double precision :: x,y,z
!     double precision :: r,maggrav

!     This is an example of how to use the radial profile above
!     do k = lo(3), hi(3)
!        z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)

!        do j = lo(2), hi(2)
!           y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)

!           do i = lo(1), hi(1)
!              x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)

!              r = sqrt(x**2+y**2+z**2)

!              maggrav = ca_prescribe_grav_gravityprofile(r)

!              Put in angular dependence
!              grav(i,j,k,1) = maggrav* ...
!              grav(i,j,k,2) = maggrav* ...
!              grav(i,j,k,3) = maggrav* ...

!        enddo
!     enddo

      grav = ZERO

      end subroutine ca_prescribe_grav



