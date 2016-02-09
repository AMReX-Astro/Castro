module actual_rhs_module

  use burn_type_module

  implicit none

contains

  subroutine actual_rhs(state)

    implicit none

    type (burn_t) :: state

    ! Do nothing in this burner.

    state % ydot = ZERO

  end subroutine actual_rhs



  subroutine actual_jac(state)

    implicit none

    type (burn_t) :: state

    ! Do nothing in this burner.

    state % jac(:,:) = ZERO

  end subroutine actual_jac

end module actual_rhs_module
