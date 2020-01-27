module rad_source_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine ca_rad_source(lo, hi, &
                           rhs, rhs_lo, rhs_hi, &
                           dx, dt, time, igroup) &
                           bind(C, name="ca_rad_source")

    use rad_params_module, only: ngroups, clight
    use prob_params_module, only: problo

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3) 
    integer,  intent(in   ) :: rhs_lo(3), rhs_hi(3)
    real(rt), intent(inout) :: rhs(rhs_lo(1):rhs_hi(1),rhs_lo(2):rhs_hi(2),rhs_lo(3):rhs_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, time
    integer,  intent(in   ), value :: igroup

    real(rt), parameter :: x0 = 0.5e0_rt
    real(rt), parameter :: t0 = 3.3356409519815202e-10_rt 
    real(rt), parameter :: qn = 1.134074546528399e20_rt 

    integer  :: i, j, k
    real(rt) :: x, y, z

    !$gpu

    do k = lo(3), hi(3)
       z = problo(3) + (dble(k) + 0.5_rt) * dx(3)
       do j = lo(2), hi(2)
          y = problo(2) + (dble(j) + 0.5_rt) * dx(2)
          do i = lo(1), hi(1)
             x = problo(1) + (dble(i) + 0.5_rt) * dx(1)

             if (time .le. t0 .and. abs(x) .le. x0) then
                rhs(i,j,k) = rhs(i,j,k) + qn ! (qn / dt) * dt
             end if

          end do
       end do
    end do

  end subroutine ca_rad_source

end module rad_source_module
