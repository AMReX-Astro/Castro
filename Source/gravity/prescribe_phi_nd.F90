module prescribe_phi_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains
  
  subroutine ca_prescribe_phi (lo,hi,phi,p_lo,p_hi,dx) &
       bind(C, name="ca_prescribe_phi")

    use bl_constants_module, only: ZERO
    ! use bl_constants_module, only: HALF
    ! use fundamental_constants_module, only: Gconst
    ! use prob_params_module, only: problo, center, dim

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: p_lo(3), p_hi(3)
    real(rt), intent(out) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    real(rt), intent(in) :: dx(3)

    ! Local variables
    !     integer          :: i, j, k
    !     real(rt)         :: x, y, z
    !     real(rt)         :: r, maggrav
    !     real(rt)         :: M_c

    !     This is an example of how to specify a radial profile.
    !     Note that in this example M_c could be saved 
    !     in another module instead, like a probdata_module.
    !     Note also that you'll have to be careful if you're working
    !     in fewer than three dimensions; you may want to set
    !     z = 0 for 2D and y = 0 for 1D.
    !
    !     M_c = 1.0e33_rt
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
    !              phi(i,j,k) = -M_c * Gconst / r
    !
    !           enddo
    !        enddo
    !     enddo

    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ZERO

  end subroutine ca_prescribe_phi

end module prescribe_phi_module
