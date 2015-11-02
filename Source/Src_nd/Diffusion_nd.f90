
subroutine ca_tempdiffextrap(lo, hi, tdif, t_lo, t_hi)

  use prob_params_module, only: dg, dim
  
  implicit none
  
  integer :: lo(3), hi(3)
  integer :: t_lo(3), t_hi(3)
  double precision :: tdif(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))

  ! Local variables
  
  integer :: i, j, k

  ! left side
  i = lo(1)-1*dg(1)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        tdif(i,j,k) = tdif(i+1*dg(1),j,k)
     end do
  end do

  ! right side
  i = hi(1)+1*dg(1)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        tdif(i,j,k) = tdif(i-1*dg(1),j,k)
     end do
  end do

  ! bottom side
  j = lo(2)-1*dg(2)
  do k = lo(3), hi(3)
     do i = lo(1), hi(1)
        tdif(i,j,k) = tdif(i,j+1*dg(2),k)
     end do
  end do

  ! top side
  j = hi(2)+1*dg(2)
  do k = lo(3), hi(3)
     do i = lo(1), hi(1)
        tdif(i,j,k) = tdif(i,j-1*dg(2),k)
     end do
  end do

  ! down side
  k = lo(3)-1*dg(3)
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        tdif(i,j,k) = tdif(i,j,k+1*dg(3))
     end do
  end do

  ! up side
  k = hi(3)+1*dg(3)
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        tdif(i,j,k) = tdif(i,j,k-1*dg(3))
     end do
  end do

  ! k-edges
  i = lo(1)-1*dg(1)
  j = lo(2)-1*dg(2)
  do k = lo(3), hi(3)
     tdif(i,j,k) = tdif(i+1*dg(1),j+1*dg(2),k)
  end do

  i = lo(1)-1*dg(1)
  j = hi(2)+1*dg(2)
  do k = lo(3), hi(3)
     tdif(i,j,k) = tdif(i+1*dg(1),j-1*dg(2),k)
  end do

  i = hi(1)+1*dg(1)
  j = lo(2)-1*dg(2)
  do k = lo(3), hi(3)
     tdif(i,j,k) = tdif(i-1*dg(1),j+1*dg(2),k)
  end do

  i = hi(1)+1*dg(1)
  j = hi(2)+1*dg(2)
  do k = lo(3), hi(3)
     tdif(i,j,k) = tdif(i-1*dg(1),j-1*dg(2),k)
  end do

  ! j-edges
  i = lo(1)-1*dg(1)
  k = lo(3)-1*dg(3)
  do j = lo(2), hi(2)
     tdif(i,j,k) = tdif(i+1*dg(1),j,k+1*dg(3))
  end do

  i = lo(1)-1*dg(1)
  k = hi(3)+1*dg(3)
  do j = lo(2), hi(2)
     tdif(i,j,k) = tdif(i+1*dg(1),j,k-1*dg(3))
  end do

  i = hi(1)+1*dg(1)
  k = lo(3)-1*dg(3)
  do j = lo(2), hi(2)
     tdif(i,j,k) = tdif(i-1*dg(1),j,k+1*dg(3))
  end do

  i = hi(1)+1*dg(1)
  k = hi(3)+1*dg(3)
  do j = lo(2), hi(2)
     tdif(i,j,k) = tdif(i-1*dg(1),j,k-1*dg(3))
  end do

  ! i-edges
  j = lo(2)-1*dg(2)
  k = lo(3)-1*dg(3)
  do i = lo(1), hi(1)
     tdif(i,j,k) = tdif(i,j+1*dg(2),k+1*dg(3))
  end do
  
  j = lo(2)-1*dg(2)
  k = hi(3)+1*dg(3)
  do i = lo(1), hi(1)
     tdif(i,j,k) = tdif(i,j+1*dg(2),k-1*dg(3))
  end do
  
  j = hi(2)+1*dg(2)
  k = lo(3)-1*dg(3)
  do i = lo(1), hi(1)
     tdif(i,j,k) = tdif(i,j-1*dg(2),k+1*dg(3))
  end do
  
  j = hi(2)+1*dg(2)
  k = hi(3)+1*dg(3)
  do i = lo(1), hi(1)
     tdif(i,j,k) = tdif(i,j-1*dg(2),k-1*dg(3))
  end do
  
  ! corners
  i = lo(1)-1*dg(1)
  j = lo(2)-1*dg(2)
  k = lo(3)-1*dg(3)
  tdif(i,j,k) = tdif(i+1*dg(1),j+1*dg(2),k+1*dg(3))

  i = lo(1)-1*dg(1)
  j = hi(2)+1*dg(2)
  k = lo(3)-1*dg(3)
  tdif(i,j,k) = tdif(i+1*dg(1),j-1*dg(2),k+1*dg(3))
  
  i = hi(1)+1*dg(1)
  j = hi(2)+1*dg(2)
  k = lo(3)-1*dg(3)
  tdif(i,j,k) = tdif(i-1*dg(1),j-1*dg(2),k+1*dg(3))
  
  i = hi(1)+1*dg(1)
  j = lo(2)-1*dg(2)
  k = lo(3)-1*dg(3)
  tdif(i,j,k) = tdif(i-1*dg(1),j+1*dg(2),k+1*dg(3))
  
  i = lo(1)-1*dg(1)
  j = lo(2)-1*dg(2)
  k = hi(3)+1*dg(3)
  tdif(i,j,k) = tdif(i+1*dg(1),j+1*dg(2),k-1*dg(3))
  
  i = lo(1)-1*dg(1)
  j = hi(2)+1*dg(2)
  k = hi(3)+1*dg(3)
  tdif(i,j,k) = tdif(i+1*dg(1),j-1*dg(2),k-1*dg(3))
  
  i = hi(1)+1*dg(1)
  j = lo(2)-1*dg(2)
  k = hi(3)+1*dg(3)
  tdif(i,j,k) = tdif(i-1*dg(1),j+1*dg(2),k-1*dg(3))
  
  i = hi(1)+1*dg(1)
  j = hi(2)+1*dg(2)
  k = hi(3)+1*dg(3)
  tdif(i,j,k) = tdif(i-1*dg(1),j-1*dg(2),k-1*dg(3))

end subroutine ca_tempdiffextrap
