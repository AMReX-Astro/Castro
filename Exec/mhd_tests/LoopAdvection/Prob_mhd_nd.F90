subroutine ca_initmag(level, time, lo, hi, &
                      nbx, mag_x, bx_lo, bx_hi, &
                      nby, mag_y, by_lo, by_hi, &
                      nbz, mag_z, bz_lo, bz_hi, &
                      dx, xlo, xhi)

  use probdata_module
  use prob_params_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module, only : M_PI, TWO, FOUR


  implicit none

  integer :: level, nbx, nby, nbz
  integer :: lo(3), hi(3)
  integer :: bx_lo(3), bx_hi(3)
  integer :: by_lo(3), by_hi(3)
  integer :: bz_lo(3), bz_hi(3)
  real(rt) :: xlo(3), xhi(3), time, dx(3)

  real(rt) :: mag_x(bx_lo(1):bx_hi(1), bx_lo(2):bx_hi(2), bx_lo(3):bx_hi(3), nbx)
  real(rt) :: mag_y(by_lo(1):by_hi(1), by_lo(2):by_hi(2), by_lo(3):by_hi(3), nby)
  real(rt) :: mag_z(bz_lo(1):bz_hi(1), bz_lo(2):bz_hi(2), bz_lo(3):bz_hi(3), nbz)

  real(rt) :: A_z(lo(1)-1:hi(1)+2, lo(2)-1:hi(2)+2, lo(3):hi(3))
  real(rt) :: x, y, z, r
  integer  :: i, j, k


    !get the vector potential A
    A_z = 0.0e0_rt

    do k = lo(3), hi(3)
        do j = lo(2)-1, hi(2)+2
           y = problo(2) + dx(2)*(dble(j))
           do i = lo(1)-1, hi(1)+2
              x = problo(1) + dx(1)*(dble(i))
             
              r = sqrt(x**2+y**2)
              if ( r .lt. 0.3 ) then
                 A_z(i,j,k) = B_0 * (0.3 - r)
              endif 
           
           enddo
        enddo
     enddo



     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1

                 mag_x(i,j,k,1) = 0.5e0_rt*(A_z(i,j+1,k) - A_z(i,j-1,k))/dx(2)
           
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)

                 mag_y(i,j,k,1) = 0.5e0_rt*(A_z(i-1,j,k) - A_z(i+1,j,k))/dx(1)
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)+1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
                 mag_z(i,j,k,1) = 0.0e0_rt
           enddo
        enddo
     enddo

  end subroutine ca_initmag

