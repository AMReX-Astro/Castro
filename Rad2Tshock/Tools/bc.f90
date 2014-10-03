program bc

  use blackbody_module, only : BGroup

  implicit none

  integer :: ngroups, igroup
  double precision, allocatable :: xnu(:)
  double precision :: T

  integer :: n
  character(len=256) :: line
  integer :: ipos

  ! open the file
  open(unit=10,file="group_structure.dat")

  ! read in the number of groups
  read(10,fmt='(a256)') line
  ipos = index(line, "=") + 1
  read (line(ipos:),*) ngroups

  allocate (xnu(0:ngroups))

  ! skip lines
  do n=1,ngroups+3
     read(10,fmt='(a256)') line
  end do

  do igroup = 0, ngroups
     read(10,fmt='(a256)') line
     ipos = index(line, "=") + 1
     read (line(ipos:),*) xnu(igroup)
  enddo

  print *, 'T = ?'
  read(*,*) T

  print *, "T = ", T
  print *, ((BGroup(T, xnu(igroup), xnu(igroup+1))), igroup=0,ngroups-1)

end program bc
