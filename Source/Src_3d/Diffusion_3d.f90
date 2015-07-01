
subroutine ca_tempdiffextrap(lo, hi, &
     tdif, t_l1, t_l2, t_l3, t_h1, t_h2, t_h3)
  integer :: lo(3), hi(3)
  integer :: t_l1, t_l2, t_l3, t_h1, t_h2, t_h3
  double precision :: tdif(t_l1:t_h1,t_l2:t_h2,t_l3:t_h3)

  integer :: i,j,k

  ! left side
  i = lo(1)-1
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        tdif(i,j,k) = tdif(i+1,j,k)
     end do
  end do

  ! right side
  i = hi(1)+1
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        tdif(i,j,k) = tdif(i-1,j,k)
     end do
  end do
  
  ! bottom side
  j = lo(2)-1
  do k = lo(3), hi(3)
     do i = lo(1), hi(1)
        tdif(i,j,k) = tdif(i,j+1,k)
     end do
  end do
  
  ! top side
  j = hi(2)+1
  do k = lo(3), hi(3)
     do i = lo(1), hi(1)
        tdif(i,j,k) = tdif(i,j-1,k)
     end do
  end do
  
  ! down side
  k = lo(3)-1
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        tdif(i,j,k) = tdif(i,j,k+1)
     end do
  end do
  
  ! up side
  k = hi(3)+1
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        tdif(i,j,k) = tdif(i,j,k-1)
     end do
  end do
  
  ! k-edges
  i = lo(1)-1
  j = lo(2)-1
  do k = lo(3), hi(3)
     tdif(i,j,k) = tdif(i+1,j+1,k)
  end do
  
  i = lo(1)-1
  j = hi(2)+1
  do k = lo(3), hi(3)
     tdif(i,j,k) = tdif(i+1,j-1,k)
  end do
  
  i = hi(1)+1
  j = lo(2)-1
  do k = lo(3), hi(3)
     tdif(i,j,k) = tdif(i-1,j+1,k)
  end do
  
  i = hi(1)+1
  j = hi(2)+1
  do k = lo(3), hi(3)
     tdif(i,j,k) = tdif(i-1,j-1,k)
  end do
  
  ! j-edges
  i = lo(1)-1
  k = lo(3)-1
  do j = lo(2), hi(2)
     tdif(i,j,k) = tdif(i+1,j,k+1)
  end do
  
  i = lo(1)-1
  k = hi(3)+1
  do j = lo(2), hi(2)
     tdif(i,j,k) = tdif(i+1,j,k-1)
  end do
  
  i = hi(1)+1
  k = lo(3)-1
  do j = lo(2), hi(2)
     tdif(i,j,k) = tdif(i-1,j,k+1)
  end do
  
  i = hi(1)+1
  k = hi(3)+1
  do j = lo(2), hi(2)
     tdif(i,j,k) = tdif(i-1,j,k-1)
  end do
  
  ! i-edges
  j = lo(2)-1
  k = lo(3)-1
  do i = lo(1), hi(1)
     tdif(i,j,k) = tdif(i,j+1,k+1)
  end do
  
  j = lo(2)-1
  k = hi(3)+1
  do i = lo(1), hi(1)
     tdif(i,j,k) = tdif(i,j+1,k-1)
  end do
  
  j = hi(2)+1
  k = lo(3)-1
  do i = lo(1), hi(1)
     tdif(i,j,k) = tdif(i,j-1,k+1)
  end do
  
  j = hi(2)+1
  k = hi(3)+1
  do i = lo(1), hi(1)
     tdif(i,j,k) = tdif(i,j-1,k-1)
  end do
  
  ! corners
  i = lo(1)-1
  j = lo(2)-1
  k = lo(3)-1
  tdif(i,j,k) = tdif(i+1,j+1,k+1)

  i = lo(1)-1
  j = hi(2)+1
  k = lo(3)-1
  tdif(i,j,k) = tdif(i+1,j-1,k+1)
  
  i = hi(1)+1
  j = hi(2)+1
  k = lo(3)-1
  tdif(i,j,k) = tdif(i-1,j-1,k+1)
  
  i = hi(1)+1
  j = lo(2)-1
  k = lo(3)-1
  tdif(i,j,k) = tdif(i-1,j+1,k+1)
  
  i = lo(1)-1
  j = lo(2)-1
  k = hi(3)+1
  tdif(i,j,k) = tdif(i+1,j+1,k-1)
  
  i = lo(1)-1
  j = hi(2)+1
  k = hi(3)+1
  tdif(i,j,k) = tdif(i+1,j-1,k-1)
  
  i = hi(1)+1
  j = lo(2)-1
  k = hi(3)+1
  tdif(i,j,k) = tdif(i-1,j+1,k-1)
  
  i = hi(1)+1
  j = hi(2)+1
  k = hi(3)+1
  tdif(i,j,k) = tdif(i-1,j-1,k-1)

end subroutine ca_tempdiffextrap
