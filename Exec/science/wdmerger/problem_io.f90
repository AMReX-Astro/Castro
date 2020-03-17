module problem_io_module

  implicit none

  ! For determining if we are the I/O processor.
  
  logical, save :: ioproc

contains

  subroutine initialize_io(name, namlen)
    
    implicit none

    integer :: namlen, i
    integer :: name(namlen)

    integer :: is_ioprocessor

    ! Determine whether we are the I/O procoessor.
    
    call bl_pd_is_ioproc(is_ioprocessor)

    ioproc = is_ioprocessor .ne. 0

    ! Read in probdata.

    call probdata_init(name, namlen)

  end subroutine initialize_io
  
end module problem_io_module
