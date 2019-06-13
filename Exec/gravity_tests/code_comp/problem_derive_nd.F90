! problem-specific Fortran derive routines go here

subroutine derextheating(h, h_lo, h_hi, ncomp_h, &
                         s, s_lo, s_hi, ncomp_s, &
                         lo, hi, domlo, domhi, &
                         dx, time) &
                         bind(C, name='derextheating')

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: HALF
    use meth_params_module, only : UEDEN
    use prob_params_module, only : center
    use probdata_module, only : heating_factor

    implicit none

    integer,  intent(in   ) :: h_lo(3), h_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt), intent(inout) :: h(h_lo(1):h_hi(1),h_lo(2):h_hi(2),h_lo(3):h_hi(3),ncomp_h)
    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ncomp_s)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: ncomp_h ! == 1
    integer,  intent(in   ), value :: ncomp_s
    real(rt), intent(in   ), value :: time

    integer :: i, j, k

    !$ gpu

    do k = lo(3), hi(3)
       ! z = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          ! y = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             ! x = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             ! r = sqrt(x**2 + y**2 + z**2)

             ! h(i,j,k,1) = heating_factor * 6.7e5_rt * exp(-(r/3.2e10_rt)**2)
             h(i,j,k,1) = s(i,j,k,UEDEN)

          end do
       end do
    end do

end subroutine derextheating
