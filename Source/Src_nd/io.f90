module io_module

  implicit none

  public

contains

  subroutine flush_output()

    use iso_fortran_env, only: output_unit

    implicit none

!$omp critical(fortran_print)
    flush(output_unit)
!$omp end critical(fortran_print)

  end subroutine flush_output

end module io_module
