module castro_util_1d_module

  implicit none

  public

contains

  subroutine ca_check_initial_species(lo,hi,state,state_l1,state_h1) &
       bind(C, name="ca_check_initial_species")

    use network           , only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(1), hi(1)
    integer          :: state_l1,state_h1
    double precision :: state(state_l1:state_h1,NVAR)

    ! Local variables
    integer          :: i,n
    double precision :: sum

    do i = lo(1), hi(1)

       sum = ZERO
       do n = 1, nspec
          sum = sum + state(i,UFS+n-1)
       end do
       if (abs(state(i,URHO)-sum).gt. 1.d-8 * state(i,URHO)) then
          print *,'Sum of (rho X)_n vs rho at (i): ',i,sum,state(i,URHO)
          call bl_error("Error:: Failed check of initial species summing to 1")
       end if

    enddo

  end subroutine ca_check_initial_species



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
