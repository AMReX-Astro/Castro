! compute the analytic solution -- this is a function of time and the
! problem parameters only

subroutine ca_deranalytic(a, a_lo, a_hi, ncomp_a, &
                          u, u_lo, u_hi, ncomp_u, &
                          lo, hi, domlo, domhi,&
                          dx, xlo, time, dt, bc, level, grid_no) bind(C, name="ca_deranalytic")

  use prob_params_module, only: center, problo
  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: HALF
  use prob_util_module, only: analytic

  implicit none

  integer,  intent(in   ) :: a_lo(3), a_hi(3), ncomp_a
  integer,  intent(in   ) :: u_lo(3), u_hi(3), ncomp_u
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: a(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3),ncomp_a)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
  integer,  intent(in   ) :: bc(3,2,ncomp_u), level, grid_no

  real(rt) :: r(3), temp
  integer :: i, j, k

  do k = lo(3), hi(3)
     r(3) = problo(3) + dx(3) * (dble(k) + HALF)

     do j = lo(2), hi(2)
        r(2) = problo(2) + dx(2) * (dble(j) + HALF)

        do i = lo(1), hi(1)
           r(1) = problo(1) + dx(1) * (dble(i) + HALF)

           call analytic(r, time, temp)
           a(i,j,k,1) = temp

        end do
     end do
  end do

end subroutine ca_deranalytic
