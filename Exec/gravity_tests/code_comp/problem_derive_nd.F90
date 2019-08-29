! problem-specific Fortran derive routines go here

subroutine derextheating(h, h_lo, h_hi, ncomp_h, &
                         s, s_lo, s_hi, ncomp_s, &
                         lo, hi, domlo, domhi, &
                         dx, time) &
                         bind(C, name='derextheating')

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, ONE, M_PI
    use prob_params_module, only : problo
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
    real(rt) :: fheat, y

    !$ gpu

    h(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:ncomp_h) = ZERO

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          y = problo(2) + (dble(j) + HALF) * dx(2) 
          if (y < 1.125e0_rt * 4.e8_rt) then 

            fheat = sin(8.e0_rt * M_PI * (y/ 4.e8_rt - ONE))

            do i = lo(1), hi(1)
    
               ! Source terms
               h(i,j,k,1) = heating_factor * fheat
    
            end do
          endif

       end do
    end do

end subroutine derextheating
