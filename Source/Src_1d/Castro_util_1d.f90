module castro_util_1d_module

  implicit none

  public

contains

  subroutine get_center(center_out) bind(C, name="get_center")

    use prob_params_module, only : center

    implicit none

    double precision, intent(inout) :: center_out(1)

    center_out(1) = center(1)

  end subroutine get_center



  subroutine set_center(center_in) bind(C, name="set_center")

    use prob_params_module, only : center

    implicit none

    double precision :: center_in(1)

    center(1) = center_in(1)

  end subroutine set_center


  subroutine find_center(data,new_center) bind(C, name="find_center")

    use bl_constants_module

    implicit none

    double precision :: data(0:2)
    double precision :: new_center(1)

    ! In 1-D it only make sense to have the center at the origin
    new_center(1) = ZERO 

  end subroutine find_center

end module castro_util_1d_module
