subroutine ca_initmag(level, time, lo, hi, &
                      nbx, mag_x, bx_lo, bx_hi, &
                      nby, mag_y, by_lo, by_hi, &
                      nbz, mag_z, bz_lo, bz_hi, &
                      delta, xlo, xhi)

  use probdata_module
  use prob_params_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer :: level, nbx, nby, nbz
  integer :: lo(3), hi(3)
  integer :: bx_lo(3), bx_hi(3)
  integer :: by_lo(3), by_hi(3)
  integer :: bz_lo(3), bz_hi(3)
  real(rt) :: xlo(3), xhi(3), time, delta(3)

  real(rt) :: mag_x(bx_lo(1):bx_hi(1), bx_lo(2):bx_hi(2), bx_lo(3):bx_hi(3), nbx)
  real(rt) :: mag_y(by_lo(1):by_hi(1), by_lo(2):by_hi(2), by_lo(3):by_hi(3), nby)
  real(rt) :: mag_z(bz_lo(1):bz_hi(1), bz_lo(2):bz_hi(2), bz_lo(3):bz_hi(3), nbz)

  real(rt) :: xcen, ycen, zcen
  integer  :: i, j, k

  print *, "Initializing magnetic field!!"

  if (idir .eq. 1) then
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
              xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)
              if (xcen <= split(1)+1) then
                 mag_x(i,j,k,1) = B_x_l
              else
                 mag_x(i,j,k,1) = B_x_r
              endif
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)
              xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)
              if (xcen <= split(1)) then
                 mag_y(i,j,k,1) = B_y_l
              else
                 mag_y(i,j,k,1) = B_y_r
              endif
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)+1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)
              if (xcen <= split(1)) then
                 mag_z(i,j,k,1) = B_z_l
              else
                 mag_z(i,j,k,1) = B_z_r
              endif
           enddo
        enddo
     enddo

  endif

  if (idir .eq. 2) then
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)
           do i = lo(1), hi(1)+1
              if (ycen <= split(2)) then
                 mag_x(i,j,k,1) = B_x_l
              else
                 mag_x(i,j,k,1) = B_x_r
              endif
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)+1
           ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)
           do i = lo(1), hi(1)
              if (ycen <= split(2)+1) then
                 mag_y(i,j,k,1) = B_y_l
              else
                 mag_y(i,j,k,1) = B_y_r
              endif
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)+1
        do j = lo(2), hi(2)
           ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)
           do i = lo(1), hi(1)
              if (ycen <= split(2)) then
                 mag_z(i,j,k,1) = B_z_l
              else
                 mag_z(i,j,k,1) = B_z_r
              endif
           enddo
        enddo
     enddo

  endif

  if (idir .eq. 3) then
     do k = lo(3), hi(3)
        zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
              if (zcen <= split(3)) then
                 mag_x(i,j,k,1) = B_x_l
              else
                 mag_x(i,j,k,1) = B_x_r
              endif
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)
        zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)
              if (zcen <= split(3)) then
                 mag_y(i,j,k,1) = B_y_l
              else
                 mag_y(i,j,k,1) = B_y_r
              endif
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)+1
        zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              if (zcen <= split(3)+1) then
                 mag_z(i,j,k,1) = B_z_l
              else
                 mag_z(i,j,k,1) = B_z_r
              endif
           enddo
        enddo
     enddo

  endif

end subroutine ca_initmag

