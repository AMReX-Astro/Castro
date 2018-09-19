! Some utility functions for initialization and analysis

module prob_util_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains

  subroutine analytic(xc, yc, time, temp)

    use prob_params_module, only : center, coord_type
    use probdata_module, only : T1, T2, diff_coeff, t_0, rho0

    implicit none

    real(rt), intent(in) :: xc, yc, time
    real(rt), intent(out) :: temp

    real(rt)         :: dist2

    dist2 = (xc - center(1))**2 + (yc - center(2))**2

    if (coord_type == 0) then
       temp = (T2 - T1)*t_0/(time + t_0) * &
            exp(-0.25_rt*dist2/(diff_coeff*(time + t_0)) ) + T1
    elseif (coord_type == 1) then
       ! this should be the 3-d solution
       temp = (T2 - T1)*(t_0/(time + t_0))**1.5_rt * &
            exp(-0.25_rt*dist2/(diff_coeff*(time + t_0)) ) + T1
    endif

  end subroutine analytic

end module prob_util_module
