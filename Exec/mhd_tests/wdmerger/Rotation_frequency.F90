module rotation_frequency_module

  implicit none

  private

  public get_omega, get_domegadt

contains

  function get_omega(time) result(omega)

    use meth_params_module, only: rot_period, rot_axis
    use bl_constants_module, only: ZERO, TWO, M_PI
    use fundamental_constants_module, only: Gconst

    implicit none

    double precision :: time
    double precision :: omega(3)

    ! If rot_period is less than zero, that means rotation is disabled, and so we should effectively
    ! shut off the source term by setting omega = 0. Note that by default rot_axis == 3 for Cartesian
    ! coordinates and rot_axis == 2 for cylindrical coordinates.

    omega = ZERO

    omega(rot_axis) = TWO * M_PI / rot_period

  end function get_omega



  function get_domegadt(time) result(domegadt)

    use bl_constants_module, only: ZERO

    implicit none

    double precision :: time
    double precision :: domegadt(3)

    domegadt = ZERO

  end function get_domegadt

end module rotation_frequency_module
