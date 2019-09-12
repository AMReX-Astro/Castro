module prescribe_grav_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  pure function grav_zone(y) result (g)

    use amrex_constants_module, only: ZERO, HALF, ONE, M_PI
    use probdata_module, only : g0

    real(rt), intent(in) :: y
    real(rt) :: g
    real(rt) :: fg

    !$gpu

    if (y < 1.0625e0_rt * 4.e8_rt) then 
       fg = HALF * (ONE + sin(16.e0_rt * M_PI * (y/4.e8_rt - 1.03125e0_rt)))
    else if (y > 2.9375e0_rt * 4.e8_rt) then
       fg = HALF * (ONE - sin(16.e0_rt * M_PI * (y/4.e8_rt - 2.96875e0_rt)))
    else
       fg = ONE
    endif

    g = fg * g0 / ((y / 4.e8_rt)**1.25e0_rt)

  end function grav_zone

  subroutine ca_prescribe_grav (lo,hi,grav,g_lo,g_hi,dx) &
       bind(C, name="ca_prescribe_grav")

    use amrex_constants_module, only: ZERO, HALF, ONE, M_PI
    use model_parser_module
    use interpolate_module
    use prob_params_module, only: dim, center, problo
    use probdata_module, only : g0

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: g_lo(3), g_hi(3)
    real(rt), intent(out):: grav(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),dim)
    real(rt), intent(in) :: dx(3)

    ! Local variables
    integer          :: i, j, k
    real(rt)         :: y, fg

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          y = problo(2) + (dble(j)+HALF) * dx(2)

          do i = lo(1), hi(1)

            grav(i,j,k,1:dim) = ZERO

            grav(i,j,k,2) = grav_zone(y)

          enddo
       enddo
    enddo

  end subroutine ca_prescribe_grav

end module prescribe_grav_module
