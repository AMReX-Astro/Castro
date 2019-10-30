! Some utility functions for initialization and analysis

module prob_util_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine analytic(r, time, temp)

    use prob_params_module, only: center, coord_type, dim
    use probdata_module, only: T1, T2, diff_coeff, t_0, rho0

    implicit none

    real(rt), intent(in   ) :: r(3), time
    real(rt), intent(inout) :: temp

    real(rt) :: dist2, exponent

    !$gpu

    if (dim == 1 .and. coord_type == 2) then
       ! Handle spherical coordinates
       exponent = 3.0_rt / 2.0_rt
    else if (dim == 2 .and. coord_type == 1) then
       ! Handle cylindrical coordinates
       exponent = 3.0_rt / 2.0_rt
    else
       exponent = dim / 2.0_rt
    end if

    dist2 = sum((r - center)**2)

    temp = T1 + (T2 - T1) * (t_0 / (time + t_0))**exponent * &
                exp(-0.25_rt * dist2 / (diff_coeff * (time + t_0)))

  end subroutine analytic

end module prob_util_module
