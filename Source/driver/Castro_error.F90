module castro_error_module

  use amrex_error_module

contains

  subroutine castro_error(message)
    ! report an error and abort

    character(len=*), intent(in) :: message
    call amrex_error(message)
  end subroutine castro_error

end module castro_error_module

