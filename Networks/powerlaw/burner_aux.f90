module rpar_indices

  implicit none

  integer, save :: n_rpar_comps = 0

  integer, save :: irp_dens, irp_temp

contains

  function get_next_rpar_index(num) result (next)

    ! return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num
    integer :: num, next

    next = n_rpar_comps + 1
    n_rpar_comps = n_rpar_comps + num

    return
  end function get_next_rpar_index


  subroutine init_rpar_indices(nspec)

    integer, intent(in) :: nspec

    irp_dens  = get_next_rpar_index(1)
    irp_temp  = get_next_rpar_index(1)

  end subroutine init_rpar_indices

end module rpar_indices
