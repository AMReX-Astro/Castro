module prescribe_grav_module

  implicit none

contains
  
  subroutine ca_prescribe_grav (lo,hi,grav,g_lo,g_hi,dx) &
       bind(C, name="ca_prescribe_grav")

    use bl_constants_module, only: ZERO
    ! use fundamental_constants_module, only: Gconst, M_PI
    ! use prob_params_module, only: problo, center, dim

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: g_lo(3), g_hi(3)
    double precision :: grav(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),3)
    double precision :: dx(3)

    ! Local variables
    !     integer          :: i, j, k
    !     double precision :: x, y, z
    !     double precision :: r, maggrav
    !     double precision :: r_c, rho_c

    !     This is an example of how to specify a radial profile.
    !     Note that in this example r_c and rho_c could be saved 
    !     in another module instead, like a probdata_module.
    !     Note also that you'll have to be careful if you're working
    !     in fewer than three dimensions; you may want to set
    !     z = 0 for 2D and y = 0 for 1D.
    !
    !     r_c = 1.0d9
    !     rho_c = 1.0d8
    !
    !     do k = lo(3), hi(3)
    !        if (dim .eq. 3) then      
    !           z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
    !        else
    !           z = ZERO
    !        endif
    !
    !        do j = lo(2), hi(2)
    !           if (dim .ge. 2) then
    !              y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
    !           else
    !              y = ZERO
    !           endif
    !
    !           do i = lo(1), hi(1)
    !              x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
    !
    !              r = sqrt(x**2+y**2+z**2)
    !
    !              maggrav =-rho_c*Gconst*2*M_PI*(1/sqrt(1+(r/r_c)**2)-atanh(1/sqrt(1+(r/r_c)**2)))
    !
    !              Put in angular dependence
    !
    !              grav(i,j,k,1) = maggrav* ...
    !              grav(i,j,k,2) = maggrav* ...
    !              grav(i,j,k,3) = maggrav* ...
    !
    !           enddo
    !        enddo
    !     enddo

    grav(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

  end subroutine ca_prescribe_grav

end module prescribe_grav_module
