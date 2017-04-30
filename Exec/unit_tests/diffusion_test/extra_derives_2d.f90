! compute the analytic solution -- this is a function of time and the
! problem parameters only

subroutine ca_deranalytic(a,a_l1,a_l2,a_h1,a_h2,ncomp_a, &
                          u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
                          domhi,dx,xlo,time,dt,bc,level,grid_no) bind(C, name="ca_deranalytic")

  use prob_params_module, only: center, problo
  use amrex_fort_module, only : rt => amrex_real
  use bl_constants_module, only : HALF
  use prob_util_module, only: analytic

  implicit none

  integer, intent(in) :: a_l1, a_l2, a_h1, a_h2, ncomp_a
  integer, intent(in) :: u_l1, u_l2, u_h1, u_h2, ncomp_u
  integer, intent(in) :: lo(2), hi(2), domlo(2), domhi(2)
  real(rt), intent(inout) :: a(a_l1:a_h1,a_l2:a_h2,ncomp_a)
  real(rt), intent(in) :: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
  real(rt), intent(in) :: dx(2), xlo(2), time, dt
  integer, intent(in) :: bc(2,2,ncomp_u), level, grid_no

  real(rt) :: xc, yc, temp
  integer :: i, j

  do j = lo(2), hi(2)
     yc = problo(2) + dx(2)*(dble(j) + HALF)

     do i = lo(1), hi(1)
        xc = problo(1) + dx(1)*(dble(i) + HALF)

        call analytic(xc, yc, time, temp)
        a(i,j,1) = temp

     enddo
  enddo

end subroutine ca_deranalytic
