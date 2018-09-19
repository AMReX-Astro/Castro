module actual_rhs_module

  implicit none

contains

  subroutine actual_rhs_init()

    implicit none

  end subroutine actual_rhs_init



  subroutine actual_rhs(state)

    !$acc routine seq

    use burn_type_module, only: burn_t
    use amrex_constants_module, only: ZERO

    implicit none

    type (burn_t) :: state

    ! Do nothing in this RHS.

    state % ydot = ZERO

  end subroutine actual_rhs



  subroutine actual_jac(state)

    !$acc routine seq

    use burn_type_module, only: burn_t
    use amrex_constants_module, only: ZERO

    implicit none

    type (burn_t) :: state

    ! Do nothing in this RHS.

    state % jac(:,:) = ZERO

  end subroutine actual_jac



  subroutine update_unevolved_species(state)

    !$acc routine seq

    use burn_type_module, only: burn_t

    implicit none

    type (burn_t)    :: state

  end subroutine update_unevolved_species

end module actual_rhs_module
