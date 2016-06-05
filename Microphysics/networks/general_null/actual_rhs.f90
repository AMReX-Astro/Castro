module actual_rhs_module

  use burn_type_module

  implicit none

contains

  subroutine actual_rhs_init()

    implicit none

  end subroutine actual_rhs_init



  subroutine actual_rhs(state)

    implicit none

    type (burn_t) :: state

    ! Do nothing in this RHS.

    state % ydot = ZERO

  end subroutine actual_rhs



  subroutine actual_jac(state)

    implicit none

    type (burn_t) :: state

    ! Do nothing in this RHS.

    state % jac(:,:) = ZERO

  end subroutine actual_jac

end module actual_rhs_module
