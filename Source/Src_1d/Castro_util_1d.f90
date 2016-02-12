module castro_util_1d_module

  implicit none

  public

contains

  subroutine ca_check_initial_species(lo,hi,state,state_l1,state_h1) bind(C, name="ca_check_initial_species")

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



  subroutine ca_enforce_nonnegative_species(uout,uout_l1,uout_h1,lo,hi) bind(C, name="ca_enforce_nonnegative_species")

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(1), hi(1)
    integer          :: uout_l1,uout_h1
    double precision :: uout(uout_l1:uout_h1,NVAR)

    ! Local variables
    integer          :: i,n
    integer          :: int_dom_spec
    logical          :: any_negative
    double precision :: dom_spec,x,eps

    eps = -1.d-16

    do i = lo(1),hi(1)

       any_negative = .false.

       ! First deal with tiny undershoots by just setting them to zero
       do n = UFS, UFS+nspec-1
          if (uout(i,n) .lt. ZERO) then
             x = uout(i,n)/uout(i,URHO)
             if (x .gt. eps) then
                uout(i,n) = ZERO
             else
                any_negative = .true.
             end if
          end if
       end do

       ! We know there are one or more undershoots needing correction 
       if (any_negative) then

          ! Find the dominant species
          dom_spec = ZERO
          int_dom_spec = 0
          do n = UFS,UFS+nspec-1
             if (uout(i,n) .gt. dom_spec) then
                dom_spec = uout(i,n)
                int_dom_spec = n
             end if
          end do

          ! Now take care of undershoots greater in magnitude than 1e-16.
          do n = UFS, UFS+nspec-1

             if (uout(i,n) .lt. ZERO) then

                x = uout(i,n)/uout(i,URHO)

                ! Here we only print the bigger negative values
                if (x .lt. -1.d-2) then
                   print *,'Correcting negative species   ',n
                   print *,'   at cell (i)                ',i
                   print *,'Negative (rho*X) is           ',uout(i,n)
                   print *,'Negative      X  is           ',x
                   print *,'Filling from dominant species ',int_dom_spec
                   print *,'  which had X =               ',&
                        uout(i,int_dom_spec) / uout(i,URHO)
                end if

                ! Take enough from the dominant species to fill the negative one.
                uout(i,int_dom_spec) = uout(i,int_dom_spec) + uout(i,n)

                ! Test that we didn't make the dominant species negative
                if (uout(i,int_dom_spec) .lt. ZERO) then 
                   print *,' Just made dominant species negative ',int_dom_spec,' at ',i
                   print *,'We were fixing species ',n,' which had value ',x
                   print *,'Dominant species became ',uout(i,int_dom_spec) / uout(i,URHO)
                   call bl_error("Error:: Castro_2d.f90 :: ca_enforce_nonnegative_species")
                end if

                ! Now set the negative species to zero
                uout(i,n) = ZERO

             end if

          enddo
       end if

    enddo

  end subroutine ca_enforce_nonnegative_species



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
