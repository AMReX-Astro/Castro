subroutine ca_tempdiffextrap(lo, hi, tdif, t_l1, t_l2, t_h1, t_h2)
  integer :: lo(2), hi(2)
  integer :: t_l1, t_l2, t_h1, t_h2
  double precision :: tdif(t_l1:t_h1,t_l2:t_h2)

  integer :: i,j

  ! left side
  i = lo(1)-1
  do j = lo(2), hi(2)
     tdif(i,j) = tdif(i+1,j)
  end do

  ! right side
  i = hi(1)+1
  do j = lo(2), hi(2)
     tdif(i,j) = tdif(i-1,j)
  end do

  ! bottom side
  j = lo(2)-1
  do i = lo(1), hi(1)
     tdif(i,j) = tdif(i,j+1)
  end do
  
  ! top side
  j = hi(2)+1
  do i = lo(1), hi(1)
     tdif(i,j) = tdif(i,j-1)
  end do
  
  ! corners
  i = lo(1)-1
  j = lo(2)-1
  tdif(i,j) = tdif(i+1,j+1)

  j = hi(2)+1
  tdif(i,j) = tdif(i+1,j-1)
  
  i = hi(1)+1
  tdif(i,j) = tdif(i-1,j-1)
  
  j = lo(2)-1
  tdif(i,j) = tdif(i-1,j+1)

end subroutine ca_tempdiffextrap
