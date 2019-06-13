module prescribe_grav_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine ca_prescribe_grav (lo,hi,grav,g_lo,g_hi,dx) &
       bind(C, name="ca_prescribe_grav")

    use amrex_constants_module, only: ZERO, HALF
    use model_parser_module
    use interpolate_module
    use prob_params_module, only: dim, center, problo

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: g_lo(3), g_hi(3)
    real(rt), intent(out):: grav(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),dim)
    real(rt), intent(in) :: dx(3)

    ! Local variables
    integer          :: i, j, k
    real(rt)         :: x, y, z
    real(rt)         :: r, maggrav
    real(rt)         :: normal(3)

    real(rt), parameter :: SMALL_GRAV = 1.e-6

    do k = lo(3), hi(3)
       z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)

       do j = lo(2), hi(2)
          y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)

          do i = lo(1), hi(1)
             x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)

             r = sqrt(x**2+y**2+z**2)

             maggrav = interpolate(r,npts_model,model_r,model_state(:,igrav_model))

             ! Put in angular dependence
             if (r > 0.0d0) then
                normal(1) = x / r
                normal(2) = y / r
                normal(3) = z / r
             endif

             grav(i,j,k,1:3) = -maggrav * normal(1:3)

          enddo
       enddo
    enddo

  end subroutine ca_prescribe_grav

end module prescribe_grav_module
