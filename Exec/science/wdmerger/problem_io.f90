module problem_io_module

  implicit none

  ! For determining if we are the I/O processor.
  
  logical, save :: ioproc

  ! Probin file

  character (len=:), allocatable, save :: probin  

contains

  subroutine initialize_io(name, namlen)
    
    implicit none

    integer :: namlen, i
    integer :: name(namlen)

    integer :: is_ioprocessor
    
    ! Build "probin" filename -- the name of the file containing the fortin namelist.
    
    allocate(character(len=namlen) :: probin)
    do i = 1, namlen
       probin(i:i) = char(name(i))
    enddo

    ! Determine whether we are the I/O procoessor.
    
    call bl_pd_is_ioproc(is_ioprocessor)

    ioproc = is_ioprocessor .ne. 0
    
  end subroutine initialize_io
  
end module problem_io_module
