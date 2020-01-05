#include "AMReX_LO_BCTYPES.H"

module habec_nd_module

  ! habec is Hypre abec, where abec is the form of the linear equation
  ! we are solving:
  ! 
  ! alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine hmac(lo, hi, &
                  mat, m_lo, m_hi, &
                  a, a_lo, a_hi, &
                  alpha) &
                  bind(C, name="hmac")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    real(rt), intent(inout) :: mat(0:2*dim,m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(in   ) :: a(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
    real(rt), intent(in   ), value :: alpha

    integer :: i, j, k

    !$gpu

    if (alpha == 0.e0_rt) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = 0.e0_rt
             end do
          end do
       end do

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = alpha * a(i,j,k)
             end do
          end do
       end do

    end if

  end subroutine hmac



  subroutine hmbc(lo, hi, &
                  mat, m_lo, m_hi, &
                  b, b_lo, b_hi, &
                  beta, dx, n) &
                  bind(C, name="hmbc")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: b_lo(3), b_hi(3)
    real(rt), intent(inout) :: mat(0:2*dim,m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(inout) :: b(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: n
    real(rt), intent(in   ), value :: beta

    integer  :: i, j, k
    real(rt) :: fac

    !$gpu

    if (n == 0) then

       fac = beta / (dx(1)**2)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i+1,j,k))
                mat(1,i,j,k) = - fac * b(i,j,k)
                mat(2,i,j,k) = - fac * b(i+1,j,k)
             end do
          end do
       end do

    elseif (n == 1) then

       fac = beta / (dx(2)**2)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i,j+1,k))
                mat(3,i,j,k) = - fac * b(i,j,k)
                mat(4,i,j,k) = - fac * b(i,j+1,k)
             end do
          end do
       end do

    else

       fac = beta / (dx(3)**2)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i,j,k+1))
                mat(5,i,j,k) = - fac * b(i,j,k)
                mat(6,i,j,k) = - fac * b(i,j,k+1)
             end do
          end do
       end do

    end if

  end subroutine hmbc



  subroutine hma2c(lo, hi, &
                   mat, m_lo, m_hi, &
                   a2, a_lo, a_hi, &
                   alpha2, n) &
                   bind(C, name="hma2c")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    real(rt), intent(inout) :: mat(0:2*dim,m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(inout) :: a2(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
    integer,  intent(in   ), value :: n
    real(rt), intent(in   ), value :: alpha2

    integer  :: i, j, k
    real(rt) :: fac

    !$gpu

    fac = 0.25e0_rt * alpha2

    if (n == 0) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) + fac * (a2(i,j,k) + a2(i+1,j,k))
                mat(1,i,j,k) = mat(1,i,j,k) + fac * a2(i,j,k)
                mat(2,i,j,k) = mat(2,i,j,k) + fac * a2(i+1,j,k)
             end do
          end do
       end do

    else if (n == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) + fac * (a2(i,j,k) + a2(i,j+1,k))
                mat(3,i,j,k) = mat(3,i,j,k) + fac * a2(i,j,k)
                mat(4,i,j,k) = mat(4,i,j,k) + fac * a2(i,j+1,k)
             end do
          end do
       end do

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) + fac * (a2(i,j,k) + a2(i,j,k+1))
                mat(5,i,j,k) = mat(5,i,j,k) + fac * a2(i,j,k)
                mat(6,i,j,k) = mat(6,i,j,k) + fac * a2(i,j,k+1)
             end do
          end do
       end do

    end if

  end subroutine hma2c



  subroutine hmcc(lo, hi, &
                  mat, m_lo, m_hi, &
                  c, c_lo, c_hi, &
                  gamma, dx, n) &
                  bind(C, name="hmcc")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    real(rt), intent(inout) :: mat(0:2*dim,m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(in   ) :: c(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: gamma
    integer,  intent(in   ), value :: n

    integer  :: i, j, k
    real(rt) :: fac

    !$gpu

    if (n == 0) then

       fac = 0.5e0_rt * gamma / dx(1)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) - fac * (c(i,j,k) - c(i+1,j,k))
                mat(1,i,j,k) = mat(1,i,j,k) - fac * c(i,j,k)
                mat(2,i,j,k) = mat(2,i,j,k) + fac * c(i+1,j,k)
             end do
          end do
       end do

    else if (n == 1) then

       fac = 0.5e0_rt * gamma / dx(2)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) - fac * (c(i,j,k) - c(i,j+1,k))
                mat(3,i,j,k) = mat(3,i,j,k) - fac * c(i,j,k)
                mat(4,i,j,k) = mat(4,i,j,k) + fac * c(i,j+1,k)
             end do
          end do
       end do

    else

       fac = 0.5e0_rt * gamma / dx(3)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) - fac * (c(i,j,k) - c(i,j,k+1))
                mat(5,i,j,k) = mat(5,i,j,k) - fac * c(i,j,k)
                mat(6,i,j,k) = mat(6,i,j,k) + fac * c(i,j,k+1)
             end do
          end do
       end do

    end if

  end subroutine hmcc



  subroutine hmd1c(lo, hi, &
                   mat, m_lo, m_hi, &
                   d1, d_lo, d_hi, &
                   delta1, dx, n) &
                   bind(C, name="hmd1c")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    real(rt), intent(inout) :: mat(0:2*dim,m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(in   ) :: d1(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: delta1
    integer,  intent(in   ), value :: n

    integer  :: i, j, k
    real(rt) :: fac

    !$gpu

    if (n == 0) then

       fac = 0.5e0_rt * delta1 / dx(1)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(1,i,j,k) = mat(1,i,j,k) - fac * d1(i,j,k)
                mat(2,i,j,k) = mat(2,i,j,k) + fac * d1(i,j,k)
             end do
          end do
       end do

    else if (n == 1) then

       fac = 0.5e0_rt * delta1 / dx(2)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(3,i,j,k) = mat(3,i,j,k) - fac * d1(i,j,k)
                mat(4,i,j,k) = mat(4,i,j,k) + fac * d1(i,j,k)
             end do
          end do
       end do

    else

       fac = 0.5e0_rt * delta1 / dx(3)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(5,i,j,k) = mat(5,i,j,k) - fac * d1(i,j,k)
                mat(6,i,j,k) = mat(6,i,j,k) + fac * d1(i,j,k)
             end do
          end do
       end do

    end if

  end subroutine hmd1c



  subroutine hmd2c(lo, hi, &
                   mat, m_lo, m_hi, &
                   d2, d_lo, d_hi, &
                   delta2, dx, n) &
                   bind(C, name="hmd2c")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    real(rt), intent(inout) :: mat(0:2*dim,m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(in   ) :: d2(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: delta2
    integer,  intent(in   ), value :: n

    integer  :: i, j, k
    real(rt) :: fac

    !$gpu

    if (n == 0) then

       fac = 0.5e0_rt * delta2 / dx(1)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) + fac * (d2(i,j,k) - d2(i+1,j,k))
                mat(1,i,j,k) = mat(1,i,j,k) - fac * d2(i,j,k)
                mat(2,i,j,k) = mat(2,i,j,k) + fac * d2(i+1,j,k)
             end do
          end do
       end do

    else if (n == 1) then

       fac = 0.5e0_rt * delta2 / dx(2)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) + fac * (d2(i,j,k) - d2(i,j+1,k))
                mat(3,i,j,k) = mat(3,i,j,k) - fac * d2(i,j,k)
                mat(4,i,j,k) = mat(4,i,j,k) + fac * d2(i,j+1,k)
             end do
          end do
       end do

    else

       fac = 0.5e0_rt * delta2 / dx(3)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) + fac * (d2(i,j,k) - d2(i,j,k+1))
                mat(5,i,j,k) = mat(5,i,j,k) - fac * d2(i,j,k)
                mat(6,i,j,k) = mat(6,i,j,k) + fac * d2(i,j,k+1)
             end do
          end do
       end do

    end if

  end subroutine hmd2c



  subroutine hmmat(lo, hi, &
                   mat, mat_lo, mat_hi, &
                   cdir, bct, bho, bcl, &
                   mask, m_lo, m_hi, &
                   b, b_lo, b_hi, &
                   beta, dx) &
                   bind(C, name="hmmat")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: mat_lo(3), mat_hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: b_lo(3), b_hi(3)
    real(rt), intent(inout) :: mat(0:2*dim,mat_lo(1):mat_hi(1),mat_lo(2):mat_hi(2),mat_lo(3):mat_hi(3))
    integer,  intent(in   ) :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(in   ) :: b(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: cdir, bct, bho
    real(rt), intent(in   ), value :: bcl, beta

    integer  :: i, j, k
    real(rt) :: h, fac, bfm, bfv
    real(rt) :: bfm2, h2, th2
    logical  :: xlo, xhi, ylo, yhi, zlo, zhi

    !$gpu

    xlo = .false.
    ylo = .false.
    zlo = .false.

    xhi = .false.
    yhi = .false.
    zhi = .false.

    if (dim == 1) then

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 1) then
          xhi = .true.
          h = dx(1)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    else if (dim == 2) then

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 2) then
          xhi = .true.
          h = dx(1)
       else if (cdir == 1) then
          ylo = .true.
          h = dx(2)
       else if (cdir == 3) then
          yhi = .true.
          h = dx(2)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    else

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 3) then
          xhi = .true.
          h = dx(1)
       else if (cdir == 1) then
          ylo = .true.
          h = dx(2)
       else if (cdir == 4) then
          yhi = .true.
          h = dx(2)
       else if (cdir == 2) then
          zlo = .true.
          h = dx(3)
       else if (cdir == 5) then
          zhi = .true.
          h = dx(3)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    end if

    fac = beta / (h**2)

    if (bct == LO_DIRICHLET) then

       if (bho >= 1) then

          h2 = 0.5e0_rt * h
          th2 = 3.e0_rt * h2
          bfm = fac * (th2 - bcl) / (bcl + h2) - fac
          bfm2 = fac * (bcl - h2) / (bcl + th2)

       else

          bfv = (beta / h) / (0.5e0_rt * h + bcl)
          bfm = bfv - fac

       end if

    else if (bct == LO_NEUMANN) then

       bfm = -fac
       bfm2 = 0.e0_rt

#ifndef AMREX_USE_GPU
    else

       call castro_error("hmmat: unsupported boundary type")
#endif

    end if

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (i-1 .ge. m_lo(1) .and. i-1 .le. m_hi(1) .and. &
                 j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                 k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (xlo .and. mask(i-1,j,k) > 0) then

                   mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
                   mat(1,i,j,k) = 0.e0_rt
                   if (bho >= 1) then
                      mat(2,i,j,k) = mat(2,i,j,k) + bfm2 * b(i,j,k)
                   end if

                end if

             else if (i+1 .ge. m_lo(1) .and. i+1 .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (xhi .and. mask(i+1,j,k) > 0) then

                   mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i+1,j,k)
                   mat(2,i,j,k) = 0.e0_rt
                   if (bho >= 1) then
                      mat(1,i,j,k) = mat(1,i,j,k) + bfm2 * b(i+1,j,k)
                   end if

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j-1 .ge. m_lo(2) .and. j-1 .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (ylo .and. mask(i,j-1,k) > 0) then

                   mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
                   mat(3,i,j,k) = 0.e0_rt
                   if (bho >= 1) then
                      mat(4,i,j,k) = mat(4,i,j,k) + bfm2 * b(i,j,k)
                   end if

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j+1 .ge. m_lo(2) .and. j+1 .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (yhi .and. mask(i,j+1,k) > 0) then

                   mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j+1,k)
                   mat(4,i,j,k) = 0.e0_rt
                   if (bho >= 1) then
                      mat(3,i,j,k) = mat(3,i,j,k) + bfm2 * b(i,j+1,k)
                   end if

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k-1 .ge. m_lo(3) .and. k-1 .le. m_hi(3)) then

                if (zlo .and. mask(i,j,k-1) > 0) then

                   mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
                   mat(5,i,j,k) = 0.e0_rt
                   if (bho >= 1) then
                      mat(6,i,j,k) = mat(6,i,j,k) + bfm2 * b(i,j,k)
                   end if

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k+1 .ge. m_lo(3) .and. k+1 .le. m_hi(3)) then

                if (zhi .and. mask(i,j,k+1) > 0) then

                   mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k+1)
                   mat(6,i,j,k) = 0.e0_rt
                   if (bho >= 1) then
                      mat(5,i,j,k) = mat(5,i,j,k) + bfm2 * b(i,j,k+1)
                   end if

                end if

             end if

          end do
       end do
    end do

  end subroutine hmmat



  subroutine hmmat3(lo, hi, &
                    lo_x, hi_x, &
                    ori_lo, idir, &
                    mat, mat_lo, mat_hi, &
                    cdir, bctype, &
                    tf, tf_lo, tf_hi, &
                    bho, bcl, &
                    mask, m_lo, m_hi, &
                    b, b_lo, b_hi, &
                    beta, dx, c, &
                    spa, spa_lo, spa_hi) &
                    bind(C, name="hmmat3")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: mat_lo(3), mat_hi(3)
    integer,  intent(in   ) :: tf_lo(3), tf_hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: b_lo(3), b_hi(3)
    integer,  intent(in   ) :: spa_lo(3), spa_hi(3)
    real(rt), intent(inout) :: mat(0:2*dim,mat_lo(1):mat_hi(1),mat_lo(2):mat_hi(2),mat_lo(3):mat_hi(3))
    integer,  intent(in   ) :: tf(tf_lo(1):tf_hi(1),tf_lo(2):tf_hi(2),tf_lo(3):tf_hi(3))
    integer,  intent(in   ) :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(in   ) :: b(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    real(rt), intent(in   ) :: spa(spa_lo(1):spa_hi(1),spa_lo(2):spa_hi(2),spa_lo(3):spa_hi(3))
    integer,  intent(in   ), value :: cdir, bctype, bho, lo_x, hi_x, ori_lo, idir
    real(rt), intent(in   ), value :: bcl, beta, c
    real(rt), intent(in   ) :: dx(3)

    integer  :: i, j, k, bct
    real(rt) :: h, fac, bfm, bfv
    real(rt) :: bfm2, h2, th2
    real(rt) :: r
    logical  :: xlo, xhi, ylo, yhi, zlo, zhi

    !$gpu

    xlo = .false.
    ylo = .false.
    zlo = .false.

    xhi = .false.
    yhi = .false.
    zhi = .false.

    if (dim == 1) then

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 1) then
          xhi = .true.
          h = dx(1)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    else if (dim == 2) then

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 2) then
          xhi = .true.
          h = dx(1)
       else if (cdir == 1) then
          ylo = .true.
          h = dx(2)
       else if (cdir == 3) then
          yhi = .true.
          h = dx(2)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    else

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 3) then
          xhi = .true.
          h = dx(1)
       else if (cdir == 1) then
          ylo = .true.
          h = dx(2)
       else if (cdir == 4) then
          yhi = .true.
          h = dx(2)
       else if (cdir == 2) then
          zlo = .true.
          h = dx(3)
       else if (cdir == 5) then
          zhi = .true.
          h = dx(3)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    end if

    fac = beta / (h**2)

    if (bctype == LO_DIRICHLET) then

       if (bho >= 1) then

          h2 = 0.5e0_rt * h
          th2 = 3.e0_rt * h2
          bfm = fac * (th2 - bcl) / (bcl + h2) - fac
          bfm2 = fac * (bcl - h2) / (bcl + th2)

       else

          bfv = (beta / h) / (0.5e0_rt * h + bcl)
          bfm = bfv - fac

       end if

    else if (bctype == LO_NEUMANN) then

       bfm = -fac
       bfm2 = 0.e0_rt

#ifndef AMREX_USE_GPU
    else

       call castro_error("hmmat3: unsupported boundary type")
#endif
    end if

    ! The -fac * b(i,j,k) term applied to the matrix diagonal is the contribution
    ! from the interior stencil which must be removed at the boundary.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             call face_metric(i, j, k, lo_x, hi_x, dx(1), idir, ori_lo, r)

             if (i-1 .ge. m_lo(1) .and. i-1 .le. m_hi(1) .and. &
                 j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                 k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (xlo .and. mask(i-1,j,k) > 0) then

                   if (bctype == -1) then
                      bct = tf(i-1,j,k)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then

                      if (bho >= 1) then
                         h2 = 0.5e0_rt * h
                         th2 = 3.e0_rt * h2
                         bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                         bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                      else
                         bfv = (beta / h) / (0.5e0_rt * h + bcl)
                         bfm = bfv * b(i,j,k)
                      end if

                   else if (bct == LO_NEUMANN) then

                      bfm  = 0.e0_rt
                      bfm2 = 0.e0_rt

                   else if (bct == LO_MARSHAK) then

                      bfv = 2.e0_rt * c * beta * r / h

                      if (bho >= 1) then
                         bfm  =  0.375e0_rt * bfv
                         bfm2 = -0.125e0_rt * bfv
                      else
                         bfm = 0.25e0_rt * bfv
                      end if

                   else if (bct == LO_SANCHEZ_POMRANING) then

                      bfv = 2.e0_rt * c * beta * r / h

                      if (bho >= 1) then
                         bfm  =  1.5e0_rt * spa(i,j,k) * bfv
                         bfm2 = -0.5e0_rt * spa(i,j,k) * bfv
                      else
                         bfm = spa(i,j,k) * bfv
                      end if

#ifndef AMREX_USE_GPU
                   else

                      call castro_error("hmmat3: unsupported boundary type")
#endif
                   end if

                   mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
                   mat(1,i,j,k) = 0.e0_rt
                   if (bho >= 1) then
                      mat(2,i,j,k) = mat(2,i,j,k) + bfm2
                   end if

                end if

             else if (i+1 .ge. m_lo(1) .and. i+1 .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (xhi .and. mask(i+1,j,k) > 0) then

                   if (bctype == -1) then
                      bct = tf(i+1,j,k)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then

                      if (bho >= 1) then
                         h2 = 0.5e0_rt * h
                         th2 = 3.e0_rt * h2
                         bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i+1,j,k)
                         bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i+1,j,k)
                      else
                         bfv = (beta / h) / (0.5e0_rt * h + bcl)
                         bfm = bfv * b(i+1,j,k)
                      end if

                   else if (bct == LO_NEUMANN) then

                      bfm  = 0.e0_rt
                      bfm2 = 0.e0_rt

                   else if (bct == LO_MARSHAK) then

                      bfv = 2.e0_rt * c * beta * r / h

                      if (bho >= 1) then
                         bfm  =  0.375e0_rt * bfv
                         bfm2 = -0.125e0_rt * bfv
                      else
                         bfm = 0.25e0_rt * bfv
                      end if

                   else if (bct == LO_SANCHEZ_POMRANING) then

                      bfv = 2.e0_rt * c * beta * r / h

                      if (bho >= 1) then
                         bfm  =  1.5e0_rt * spa(i,j,k) * bfv
                         bfm2 = -0.5e0_rt * spa(i,j,k) * bfv
                      else
                         bfm = spa(i,j,k) * bfv
                      end if

#ifndef AMREX_USE_GPU
                   else

                      call castro_error("hmmat3: unsupported boundary type")
#endif
                   end if

                   mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i+1,j,k)
                   mat(2,i,j,k) = 0.e0_rt
                   if (bho >= 1) then
                      mat(1,i,j,k) = mat(1,i,j,k) + bfm2
                   end if

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j-1 .ge. m_lo(2) .and. j-1 .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (ylo .and. mask(i,j-1,k) > 0) then

                   if (bctype == -1) then
                      bct = tf(i,j-1,k)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then

                      if (bho >= 1) then
                         h2 = 0.5e0_rt * h
                         th2 = 3.e0_rt * h2
                         bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                         bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                      else
                         bfv = (beta / h) / (0.5e0_rt * h + bcl)
                         bfm = bfv * b(i,j,k)
                      end if

                   else if (bct == LO_NEUMANN) then

                      bfm  = 0.e0_rt
                      bfm2 = 0.e0_rt

                   else if (bct == LO_MARSHAK) then

                      bfv = 2.e0_rt * c * beta * r / h

                      if (bho >= 1) then
                         bfm  =  0.375e0_rt * bfv
                         bfm2 = -0.125e0_rt * bfv
                      else
                         bfm = 0.25e0_rt * bfv
                      end if

                   else if (bct == LO_SANCHEZ_POMRANING) then

                      bfv = 2.e0_rt * c * beta * r / h

                      if (bho >= 1) then
                         bfm  =  1.5e0_rt * spa(i,j,k) * bfv
                         bfm2 = -0.5e0_rt * spa(i,j,k) * bfv
                      else
                         bfm = spa(i,j,k) * bfv
                      end if

#ifndef AMREX_USE_GPU
                   else

                      call castro_error("hmmat3: unsupported boundary type")

#endif
                   end if

                   mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
                   mat(3,i,j,k) = 0.e0_rt
                   if (bho >= 1) then
                      mat(4,i,j,k) = mat(4,i,j,k) + bfm2
                   end if

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j+1 .ge. m_lo(2) .and. j+1 .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (yhi .and. mask(i,j+1,k) > 0) then

                   if (bctype == -1) then
                      bct = tf(i,j+1,k)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then

                      if (bho >= 1) then
                         h2 = 0.5e0_rt * h
                         th2 = 3.e0_rt * h2
                         bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j+1,k)
                         bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j+1,k)
                      else
                         bfv = (beta / h) / (0.5e0_rt * h + bcl)
                         bfm = bfv * b(i,j+1,k)
                      end if

                   else if (bct == LO_NEUMANN) then

                      bfm  = 0.e0_rt
                      bfm2 = 0.e0_rt

                   else if (bct == LO_MARSHAK) then

                      bfv = 2.e0_rt * c * beta * r / h

                      if (bho >= 1) then
                         bfm  =  0.375e0_rt * bfv
                         bfm2 = -0.125e0_rt * bfv
                      else
                         bfm = 0.25e0_rt * bfv
                      end if

                   else if (bct == LO_SANCHEZ_POMRANING) then

                      bfv = 2.e0_rt * c * beta * r / h

                      if (bho >= 1) then
                         bfm  =  1.5e0_rt * spa(i,j,k) * bfv
                         bfm2 = -0.5e0_rt * spa(i,j,k) * bfv
                      else
                         bfm = spa(i,j,k) * bfv
                      end if

#ifndef AMREX_USE_GPU
                   else

                      call castro_error("hmmat3: unsupported boundary type")
#endif
                   end if

                   mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j+1,k)
                   mat(4,i,j,k) = 0.e0_rt
                   if (bho >= 1) then
                      mat(3,i,j,k) = mat(3,i,j,k) + bfm2
                   end if

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k-1 .ge. m_lo(3) .and. k-1 .le. m_hi(3)) then

                if (zlo .and. mask(i,j,k-1) > 0) then

                   if (bctype == -1) then
                      bct = tf(i,j,k-1)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then

                      if (bho >= 1) then
                         h2 = 0.5e0_rt * h
                         th2 = 3.e0_rt * h2
                         bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                         bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                      else
                         bfv = (beta / h) / (0.5e0_rt * h + bcl)
                         bfm = bfv * b(i,j,k)
                      end if

                   else if (bct == LO_NEUMANN) then

                      bfm  = 0.e0_rt
                      bfm2 = 0.e0_rt

                   else if (bct == LO_MARSHAK) then

                      bfv = 2.e0_rt * c * beta * r / h

                      if (bho >= 1) then
                         bfm  =  0.375e0_rt * bfv
                         bfm2 = -0.125e0_rt * bfv
                      else
                         bfm = 0.25e0_rt * bfv
                      end if

                   else if (bct == LO_SANCHEZ_POMRANING) then

                      bfv = 2.e0_rt * c * beta * r / h

                      if (bho >= 1) then
                         bfm  =  1.5e0_rt * spa(i,j,k) * bfv
                         bfm2 = -0.5e0_rt * spa(i,j,k) * bfv
                      else
                         bfm = spa(i,j,k) * bfv
                      end if

#ifndef AMREX_USE_GPU
                   else

                      call castro_error("hmmat3: unsupported boundary type")
#endif
                   end if

                   mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
                   mat(5,i,j,k) = 0.e0_rt
                   if (bho >= 1) then
                      mat(6,i,j,k) = mat(6,i,j,k) + bfm2
                   end if

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k+1 .ge. m_lo(3) .and. k+1 .le. m_hi(3)) then

                if (zhi .and. mask(i,j,k+1) > 0) then

                   if (bctype == -1) then
                      bct = tf(i,j,k+1)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then

                      if (bho >= 1) then
                         h2 = 0.5e0_rt * h
                         th2 = 3.e0_rt * h2
                         bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k+1)
                         bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k+1)
                      else
                         bfv = (beta / h) / (0.5e0_rt * h + bcl)
                         bfm = bfv * b(i,j,k+1)
                      end if

                   else if (bct == LO_NEUMANN) then

                      bfm  = 0.e0_rt
                      bfm2 = 0.e0_rt

                   else if (bct == LO_MARSHAK) then

                      bfv = 2.e0_rt * c * beta * r / h

                      if (bho >= 1) then
                         bfm  =  0.375e0_rt * bfv
                         bfm2 = -0.125e0_rt * bfv
                      else
                         bfm = 0.25e0_rt * bfv
                      end if

                   else if (bct == LO_SANCHEZ_POMRANING) then

                      bfv = 2.e0_rt * c * beta * r / h

                      if (bho >= 1) then
                         bfm  =  1.5e0_rt * spa(i,j,k) * bfv
                         bfm2 = -0.5e0_rt * spa(i,j,k) * bfv
                      else
                         bfm = spa(i,j,k) * bfv
                      end if

#ifndef AMREX_USE_GPU
                   else

                      call castro_error("hmmat3: unsupported boundary type")
#endif
                   end if

                   mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k+1)
                   mat(6,i,j,k) = 0.e0_rt
                   if (bho >= 1) then
                      mat(5,i,j,k) = mat(5,i,j,k) + bfm2
                   end if

                end if

             end if

          end do
       end do
    end do

  end subroutine hmmat3




  subroutine hacoef(lo, hi, &
                    mat, m_lo, m_hi, &
                    a, a_lo, a_hi, &
                    alpha) &
                    bind(C, name="hacoef")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    real(rt), intent(inout) :: mat(0:dim,m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(in   ) :: a(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
    real(rt), intent(in   ), value :: alpha

    integer :: i, j, k

    !$gpu

    if (alpha == 0.e0_rt) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(dim,i,j,k) = 0.e0_rt
             end do
          end do
       end do

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(dim,i,j,k) = alpha * a(i,j,k)
             end do
          end do
       end do

    end if

  end subroutine hacoef



  subroutine hbcoef(lo, hi, &
                    mat, m_lo, m_hi, &
                    b, b_lo, b_hi, &
                    beta, dx, idir) &
                    bind(C, name="hbcoef")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: b_lo(3), b_hi(3)
    real(rt), intent(inout) :: mat(0:dim,m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(in   ) :: b(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: beta
    integer,  intent(in   ), value :: idir

    integer  :: i, j, k
    real(rt) :: fac

    !$gpu

    if (idir == 0) then

       fac = beta / (dx(1)**2)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = -fac * b(i,j,k)
                mat(dim,i,j,k) = mat(dim,i,j,k) + fac * (b(i,j,k) + b(i+1,j,k))
             end do
          end do
       end do

    elseif (idir == 1) then

       fac = beta / (dx(2)**2)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(1,i,j,k) = -fac * b(i,j,k)
                mat(dim,i,j,k) = mat(dim,i,j,k) + fac * (b(i,j,k) + b(i,j+1,k))
             end do
          end do
       end do

    else

       fac = beta / (dx(3)**2)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(2,i,j,k) = -fac * b(i,j,k)
                mat(dim,i,j,k) = mat(dim,i,j,k) + fac * (b(i,j,k) + b(i,j,k+1))
             end do
          end do
       end do

    end if

  end subroutine hbcoef



  subroutine hbmat(lo, hi, &
                   mat, a_lo, a_hi, &
                   cdir, bct, bcl, &
                   mask, m_lo, m_hi, &
                   b, b_lo, b_hi, &
                   beta, dx) &
                   bind(C, name="hbmat")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: b_lo(3), b_hi(3)
    real(rt), intent(inout) :: mat(0:dim,a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
    integer,  intent(in   ) :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(in   ) :: b(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: cdir, bct
    real(rt), intent(in   ), value :: bcl, beta

    integer  :: i, j, k
    real(rt) :: h, fac, bfm, bfv
    logical  :: xlo, xhi, ylo, yhi, zlo, zhi

    !$gpu

    xlo = .false.
    ylo = .false.
    zlo = .false.

    xhi = .false.
    yhi = .false.
    zhi = .false.

    if (dim == 1) then

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 1) then
          xhi = .true.
          h = dx(1)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    else if (dim == 2) then

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 2) then
          xhi = .true.
          h = dx(1)
       else if (cdir == 1) then
          ylo = .true.
          h = dx(2)
       else if (cdir == 3) then
          yhi = .true.
          h = dx(2)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    else

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 3) then
          xhi = .true.
          h = dx(1)
       else if (cdir == 1) then
          ylo = .true.
          h = dx(2)
       else if (cdir == 4) then
          yhi = .true.
          h = dx(2)
       else if (cdir == 2) then
          zlo = .true.
          h = dx(3)
       else if (cdir == 5) then
          zhi = .true.
          h = dx(3)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    end if

    fac = beta / (h**2)

    if (bct == LO_DIRICHLET) then
       bfv = fac * h / (0.5e0_rt * h + bcl)
       bfm = bfv - fac
    else if (bct == LO_NEUMANN) then
       bfv = beta / h
       bfm = -fac
#ifndef AMREX_USE_GPU
    else
       call castro_error("hbmat: unsupported boundary type")
#endif
    end if

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (i-1 .ge. m_lo(1) .and. i-1 .le. m_hi(1) .and. &
                 j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                 k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (xlo .and. mask(i-1,j,k) > 0) then

                   mat(dim,i,j,k) = mat(dim,i,j,k) + bfm * b(i,j,k)
                   mat(0,i,j,k) = 0.e0_rt

                end if

             else if (i+1 .ge. m_lo(1) .and. i+1 .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (xhi .and. mask(i+1,j,k) > 0) then

                   mat(dim,i,j,k) = mat(dim,i,j,k) + bfm * b(i+1,j,k)

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j-1 .ge. m_lo(2) .and. j-1 .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (ylo .and. mask(i,j-1,k) > 0) then

                   mat(dim,i,j,k) = mat(dim,i,j,k) + bfm * b(i,j,k)
                   mat(1,i,j,k) = 0.e0_rt

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j+1 .ge. m_lo(2) .and. j+1 .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (yhi .and. mask(i,j+1,k) > 0) then

                   mat(dim,i,j,k) = mat(dim,i,j,k) + bfm * b(i,j+1,k)

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k-1 .ge. m_lo(3) .and. k-1 .le. m_hi(3)) then

                if (zlo .and. mask(i,j,k-1) > 0) then

                   mat(dim,i,j,k) = mat(dim,i,j,k) + bfm * b(i,j,k)
                   mat(2,i,j,k) = 0.e0_rt

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k+1 .ge. m_lo(3) .and. k+1 .le. m_hi(3)) then

                if (zhi .and. mask(i,j,k+1) > 0) then

                   mat(3,i,j,k) = mat(3,i,j,k) + bfm * b(i,j,k+1)

                end if

             end if

          end do
       end do
    end do

  end subroutine hbmat



  subroutine hbmat3(lo, hi, &
                    lo_x, hi_x, &
                    ori_lo, idir, &
                    mat, a_lo, a_hi, &
                    cdir, bctype, &
                    tf, t_lo, t_hi, &
                    bcl, &
                    mask, m_lo, m_hi, &
                    b, b_lo, b_hi, &
                    beta, dx, c, &
                    spa, s_lo, s_hi) &
                    bind(C, name="hbmat3")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: b_lo(3), b_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: mat(0:dim,a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
    integer,  intent(in   ) :: tf(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    integer,  intent(in   ) :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(in   ) :: b(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    real(rt), intent(in   ) :: spa(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: cdir, bctype, lo_x, hi_x, ori_lo, idir
    real(rt), intent(in   ), value :: bcl, beta, c

    integer  :: i, j, k
    integer  :: bct
    real(rt) :: h, fac, bfm, bfv
    real(rt) :: r
    logical  :: xlo, xhi, ylo, yhi, zlo, zhi

    !$gpu

    xlo = .false.
    ylo = .false.
    zlo = .false.

    xhi = .false.
    yhi = .false.
    zhi = .false.

    if (dim == 1) then

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 1) then
          xhi = .true.
          h = dx(1)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    else if (dim == 2) then

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 2) then
          xhi = .true.
          h = dx(1)
       else if (cdir == 1) then
          ylo = .true.
          h = dx(2)
       else if (cdir == 3) then
          yhi = .true.
          h = dx(2)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    else

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 3) then
          xhi = .true.
          h = dx(1)
       else if (cdir == 1) then
          ylo = .true.
          h = dx(2)
       else if (cdir == 4) then
          yhi = .true.
          h = dx(2)
       else if (cdir == 2) then
          zlo = .true.
          h = dx(3)
       else if (cdir == 5) then
          zhi = .true.
          h = dx(3)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    end if

    fac = beta / (h**2)

    ! The -fac * b(i,j,k) term applied to the matrix diagonal is the contribution
    ! from the interior stencil which must be removed at the boundary.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             call face_metric(i, j, k, lo_x, hi_x, dx(1), idir, ori_lo, r)

             if (i-1 .ge. m_lo(1) .and. i-1 .le. m_hi(1) .and. &
                 j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                 k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (xlo .and. mask(i-1,j,k) > 0) then

                   if (bctype == -1) then
                      bct = tf(i-1,j,k)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then
                      bfv = fac * h / (0.5e0_rt * h + bcl)
                      bfm = bfv * b(i,j,k)
                   else if (bct == LO_NEUMANN) then
                      bfm = 0.e0_rt
                   else if (bct == LO_MARSHAK) then
                      bfv = 2.e0_rt * beta * r / h
                      bfm = 0.25e0_rt * c * bfv
                   else if (bct == LO_SANCHEZ_POMRANING) then
                      bfv = 2.e0_rt * beta * r / h
                      bfm = spa(i,j,k) * c * bfv
#ifndef AMREX_USE_GPU
                   else
                      call castro_error("hbmat3: unsupported boundary type")
#endif
                   end if

                   mat(dim,i,j,k) = mat(dim,i,j,k) + bfm - fac * b(i,j,k)
                   mat(0,i,j,k) = 0.e0_rt

                end if

             else if (i+1 .ge. m_lo(1) .and. i+1 .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (xhi .and. mask(i+1,j,k) > 0) then

                   if (bctype == -1) then
                      bct = tf(i+1,j,k)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then
                      bfv = fac * h / (0.5e0_rt * h + bcl)
                      bfm = bfv * b(i+1,j,k)
                   else if (bct == LO_NEUMANN) then
                      bfm = 0.e0_rt
                   else if (bct == LO_MARSHAK) then
                      bfv = 2.e0_rt * beta * r / h
                      bfm = 0.25e0_rt * c * bfv
                   else if (bct == LO_SANCHEZ_POMRANING) then
                      bfv = 2.e0_rt * beta * r / h
                      bfm = spa(i,j,k) * c * bfv
#ifndef AMREX_USE_GPU
                   else
                      call castro_error("hbmat3: unsupported boundary type")
#endif                      
                   end if

                   mat(dim,i,j,k) = mat(dim,i,j,k) + bfm - fac * b(i+1,j,k)

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j-1 .ge. m_lo(2) .and. j-1 .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (ylo .and. mask(i,j-1,k) > 0) then

                   if (bctype == -1) then
                      bct = tf(i,j-1,k)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then
                      bfv = fac * h / (0.5e0_rt * h + bcl)
                      bfm = bfv * b(i,j,k)
                   else if (bct == LO_NEUMANN) then
                      bfm = 0.e0_rt
                   else if (bct == LO_MARSHAK) then
                      bfv = 2.e0_rt * beta * r / h
                      bfm = 0.25e0_rt * c * bfv
                   else if (bct == LO_SANCHEZ_POMRANING) then
                      bfv = 2.e0_rt * beta * r / h
                      bfm = spa(i,j,k) * c * bfv
#ifndef AMREX_USE_GPU
                   else
                      call castro_error("hbmat3: unsupported boundary type")
#endif
                   end if

                   mat(dim,i,j,k) = mat(dim,i,j,k) + bfm - fac * b(i,j,k)
                   mat(1,i,j,k) = 0.e0_rt

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j+1 .ge. m_lo(2) .and. j+1 .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (yhi .and. mask(i,j+1,k) > 0) then

                   if (bctype == -1) then
                      bct = tf(i,j+1,k)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then
                      bfv = fac * h / (0.5e0_rt * h + bcl)
                      bfm = bfv * b(i,j+1,k)
                   else if (bct == LO_NEUMANN) then
                      bfm = 0.e0_rt
                   else if (bct == LO_MARSHAK) then
                      bfv = 2.e0_rt * beta * r / h
                      bfm = 0.25e0_rt * c * bfv
                   else if (bct == LO_SANCHEZ_POMRANING) then
                      bfv = 2.e0_rt * beta * r / h
                      bfm = spa(i,j,k) * c * bfv
#ifndef AMREX_USE_GPU
                   else
                      call castro_error("hbmat3: unsupported boundary type")
#endif
                   end if

                   mat(dim,i,j,k) = mat(dim,i,j,k) + bfm - fac * b(i,j+1,k)

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k-1 .ge. m_lo(3) .and. k-1 .le. m_hi(3)) then

                if (zlo .and. mask(i,j,k-1) > 0) then

                   if (bctype == -1) then
                      bct = tf(i,j,k-1)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then
                      bfv = fac * h / (0.5e0_rt * h + bcl)
                      bfm = bfv * b(i,j,k)
                   else if (bct == LO_NEUMANN) then
                      bfm = 0.e0_rt
                   else if (bct == LO_MARSHAK) then
                      bfv = 2.e0_rt * beta * r / h
                      bfm = 0.25e0_rt * c * bfv
                   else if (bct == LO_SANCHEZ_POMRANING) then
                      bfv = 2.e0_rt * beta * r / h
                      bfm = spa(i,j,k) * c * bfv
#ifndef AMREX_USE_GPU
                   else
                      call castro_error("hbmat3: unsupported boundary type")
#endif
                   end if

                   mat(dim,i,j,k) = mat(dim,i,j,k) + bfm - fac * b(i,j,k)
                   mat(2,i,j,k) = 0.e0_rt

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k+1 .ge. m_lo(3) .and. k+1 .le. m_hi(3)) then

                if (zhi .and. mask(i,j,k+1) > 0) then

                   if (bctype == -1) then
                      bct = tf(i,j,k+1)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then
                      bfv = fac * h / (0.5e0_rt * h + bcl)
                      bfm = bfv * b(i,j,k+1)
                   else if (bct == LO_NEUMANN) then
                      bfm = 0.e0_rt
                   else if (bct == LO_MARSHAK) then
                      bfv = 2.e0_rt * beta * r / h
                      bfm = 0.25e0_rt * c * bfv
                   else if (bct == LO_SANCHEZ_POMRANING) then
                      bfv = 2.e0_rt * beta * r / h
                      bfm = spa(i,j,k) * c * bfv
#ifndef AMREX_USE_GPU
                   else
                      call castro_error("hbmat3: unsupported boundary type")
#endif
                   end if

                   mat(dim,i,j,k) = mat(dim,i,j,k) + bfm - fac * b(i,j,k+1)

                end if

             end if

          end do
       end do
    end do

  end subroutine hbmat3



  subroutine hbvec(lo, hi, &
                   vec, v_lo, v_hi, &
                   cdir, bct, bho, bcl, &
                   bcval, c_lo, c_hi, &
                   mask, m_lo, m_hi, &
                   b, b_lo, b_hi, &
                   beta, dx) &
                   bind(C, name="hbvec")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: b_lo(3), b_hi(3)
    real(rt), intent(inout) :: vec(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    real(rt), intent(in   ) :: bcval(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    integer,  intent(in   ) :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(rt), intent(in   ) :: b(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: cdir, bct, bho
    real(rt), intent(in   ), value :: bcl, beta

    integer  :: i, j, k
    real(rt) :: h, bfv
    real(rt) :: h2, th2
    logical  :: xlo, xhi, ylo, yhi, zlo, zhi

    !$gpu

    xlo = .false.
    ylo = .false.
    zlo = .false.

    xhi = .false.
    yhi = .false.
    zhi = .false.

    if (dim == 1) then

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 1) then
          xhi = .true.
          h = dx(1)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    else if (dim == 2) then

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 2) then
          xhi = .true.
          h = dx(1)
       else if (cdir == 1) then
          ylo = .true.
          h = dx(2)
       else if (cdir == 3) then
          yhi = .true.
          h = dx(2)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    else

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 3) then
          xhi = .true.
          h = dx(1)
       else if (cdir == 1) then
          ylo = .true.
          h = dx(2)
       else if (cdir == 4) then
          yhi = .true.
          h = dx(2)
       else if (cdir == 2) then
          zlo = .true.
          h = dx(3)
       else if (cdir == 5) then
          zhi = .true.
          h = dx(3)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    end if

    if (bct == LO_DIRICHLET) then
       if (bho >= 1) then
          h2 = 0.5e0_rt * h
          th2 = 3.e0_rt * h2
          bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
       else
          bfv = (beta / h) / (0.5e0_rt * h + bcl)
       end if
    else if (bct == LO_NEUMANN) then
       bfv = beta / h
#ifndef AMREX_USE_GPU
    else
       call castro_error("hbvec: unsupported boundary type")
#endif
    end if

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (i-1 .ge. m_lo(1) .and. i-1 .le. m_hi(1) .and. &
                 j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                 k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (xlo .and. mask(i-1,j,k) > 0) then

                   vec(i,j,k) = vec(i,j,k) + bfv * b(i,j,k) * bcval(i-1,j,k)

                end if

             else if (i+1 .ge. m_lo(1) .and. i+1 .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then


                if (xhi .and. mask(i+1,j,k) > 0) then

                   vec(i,j,k) = vec(i,j,k) + bfv * b(i+1,j,k) * bcval(i+1,j,k)

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j-1 .ge. m_lo(2) .and. j-1 .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (ylo .and. mask(i,j-1,k) > 0) then

                   vec(i,j,k) = vec(i,j,k) + bfv * b(i,j,k) * bcval(i,j-1,k)

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j+1 .ge. m_lo(2) .and. j+1 .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (yhi .and. mask(i,j+1,k) > 0) then

                   vec(i,j,k) = vec(i,j,k) + bfv * b(i,j+1,k) * bcval(i,j+1,k)

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k-1 .ge. m_lo(3) .and. k-1 .le. m_hi(3)) then

                if (zlo .and. mask(i,j,k-1) > 0) then

                   vec(i,j,k) = vec(i,j,k) + bfv * b(i,j,k) * bcval(i,j,k-1)

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k+1 .ge. m_lo(3) .and. k+1 .le. m_hi(3)) then

                if (zhi .and. mask(i,j,k+1) > 0) then

                   vec(i,j,k) = vec(i,j,k) + bfv * b(i,j,k+1) * bcval(i,j,k+1)

                end if

             end if

          end do
       end do
    end do

  end subroutine hbvec



  subroutine hbvec3(lo, hi, &
                    lo_x, hi_x, &
                    ori_lo, idir, &
                    vec, v_lo, v_hi, &
                    cdir, bctype, &
                    tf, t_lo, t_hi, &
                    bho, bcl, &
                    bcval, c_lo, c_hi, &
                    mask, m_lo, m_hi, &
                    b, b_lo, b_hi, &
                    beta, dx) &
                    bind(C, name="hbvec3")

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ONE
    use prob_params_module, only: dim
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    integer,  intent(in   ) :: m_lo(3), m_hi(3)
    integer,  intent(in   ) :: b_lo(3), b_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    real(rt), intent(inout) :: vec(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    real(rt), intent(in   ) :: bcval(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(rt), intent(in   ) :: b(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    integer,  intent(in   ) :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    integer,  intent(in   ) :: tf(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: bcl, beta
    integer,  intent(in   ), value :: cdir, bctype, bho, lo_x, hi_x, ori_lo, idir

    integer  :: i, j, k
    integer  :: bct
    real(rt) :: h, bfv
    real(rt) :: h2, th2
    real(rt) :: r
    logical  :: xlo, xhi, ylo, yhi, zlo, zhi

    !$gpu

    xlo = .false.
    ylo = .false.
    zlo = .false.

    xhi = .false.
    yhi = .false.
    zhi = .false.

    if (dim == 1) then

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 1) then
          xhi = .true.
          h = dx(1)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    else if (dim == 2) then

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 2) then
          xhi = .true.
          h = dx(1)
       else if (cdir == 1) then
          ylo = .true.
          h = dx(2)
       else if (cdir == 3) then
          yhi = .true.
          h = dx(2)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    else

       if (cdir == 0) then
          xlo = .true.
          h = dx(1)
       else if (cdir == 3) then
          xhi = .true.
          h = dx(1)
       else if (cdir == 1) then
          ylo = .true.
          h = dx(2)
       else if (cdir == 4) then
          yhi = .true.
          h = dx(2)
       else if (cdir == 2) then
          zlo = .true.
          h = dx(3)
       else if (cdir == 5) then
          zhi = .true.
          h = dx(3)
#ifndef AMREX_USE_GPU
       else
          call castro_error("Unknown cdir")
#endif
       end if

    end if

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             call face_metric(i, j, k, lo_x, hi_x, dx(1), idir, ori_lo, r)

             if (i-1 .ge. m_lo(1) .and. i-1 .le. m_hi(1) .and. &
                 j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                 k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (xlo .and. mask(i-1,j,k) > 0) then

                   if (bctype == -1) then
                      bct = tf(i-1,j,k)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then
                      if (bho >= 1) then
                         h2 = 0.5e0_rt * h
                         th2 = 3.e0_rt * h2
                         bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                      else
                         bfv = (beta / h) / (0.5e0_rt * h + bcl)
                      end if
                      bfv = bfv * b(i,j,k)
                   else if (bct == LO_NEUMANN) then
                      bfv = beta * r / h
                   else if (bct == LO_MARSHAK .OR. bct == LO_SANCHEZ_POMRANING) then
                      bfv = 2.e0_rt * beta * r / h
#ifndef AMREX_USE_GPU
                   else
                      call castro_error("hbvec3: unsupported boundary type")
#endif
                   end if

                   vec(i,j,k) = vec(i,j,k) + bfv * bcval(i-1,j,k)

                end if

             else if (i+1 .ge. m_lo(1) .and. i+1 .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (xhi .and. mask(i+1,j,k) > 0) then

                   if (bctype == -1) then
                      bct = tf(i+1,j,k)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then
                      if (bho >= 1) then
                         h2 = 0.5e0_rt * h
                         th2 = 3.e0_rt * h2
                         bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                      else
                         bfv = (beta / h) / (0.5e0_rt * h + bcl)
                      end if
                      bfv = bfv * b(i+1,j,k)
                   else if (bct == LO_NEUMANN) then
                      bfv = beta * r / h
                   else if (bct == LO_MARSHAK .OR. bct == LO_SANCHEZ_POMRANING) then
                      bfv = 2.e0_rt * beta * r / h
#ifndef AMREX_USE_GPU
                   else
                      call castro_error("hbvec3: unsupported boundary type")
#endif
                   end if

                   vec(i,j,k) = vec(i,j,k) + bfv * bcval(i+1,j,k)

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j-1 .ge. m_lo(2) .and. j-1 .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (ylo .and. mask(i,j-1,k) > 0) then

                   if (bctype == -1) then
                      bct = tf(i,j-1,k)
                   else
                      bct = bctype
                   endif

                   if (bct == LO_DIRICHLET) then
                      if (bho >= 1) then
                         h2 = 0.5e0_rt * h
                         th2 = 3.e0_rt * h2
                         bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                      else
                         bfv = (beta / h) / (0.5e0_rt * h + bcl)
                      endif
                      bfv = bfv * b(i,j,k)
                   else if (bct == LO_NEUMANN) then
                      bfv = beta * r / h
                   else if (bct == LO_MARSHAK .OR. bct == LO_SANCHEZ_POMRANING) then
                      bfv = 2.e0_rt * beta * r / h
#ifndef AMREX_USE_GPU
                   else
                      call castro_error("hbvec3: unsupported boundary type")
#endif
                   end if

                   vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j-1,k)

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j+1 .ge. m_lo(2) .and. j+1 .le. m_hi(2) .and. &
                      k   .ge. m_lo(3) .and. k   .le. m_hi(3)) then

                if (yhi .and. mask(i,j+1,k) > 0) then

                   if (bctype == -1) then
                      bct = tf(i,j+1,k)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then
                      if (bho >= 1) then
                         h2 = 0.5e0_rt * h
                         th2 = 3.e0_rt * h2
                         bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                      else
                         bfv = (beta / h) / (0.5e0_rt * h + bcl)
                      end if
                      bfv = bfv * b(i,j+1,k)
                   else if (bct == LO_NEUMANN) then
                      bfv = beta * r / h
                   else if (bct == LO_MARSHAK .OR. bct == LO_SANCHEZ_POMRANING) then
                      bfv = 2.e0_rt * beta * r / h
#ifndef AMREX_USE_GPU
                   else
                      call castro_error("hbvec3: unsupported boundary type")
#endif
                   end if

                   vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j+1,k)

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k-1 .ge. m_lo(3) .and. k-1 .le. m_hi(3)) then

                if (zlo .and. mask(i,j,k-1) > 0) then

                   if (bctype == -1) then
                      bct = tf(i,j,k-1)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then
                      if (bho >= 1) then
                         h2 = 0.5e0_rt * h
                         th2 = 3.e0_rt * h2
                         bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                      else
                         bfv = (beta / h) / (0.5e0_rt * h + bcl)
                      end if
                      bfv = bfv * b(i,j,k)
                   else if (bct == LO_NEUMANN) then
                      bfv = beta * r / h
                   else if (bct == LO_MARSHAK .OR. bct == LO_SANCHEZ_POMRANING) then
                      bfv = 2.e0_rt * beta * r / h
#ifndef AMREX_USE_GPU
                   else
                      call castro_error("hbvec3: unsupported boundary type")
#endif
                   end if

                   vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j,k-1)

                end if

             else if (i   .ge. m_lo(1) .and. i   .le. m_hi(1) .and. &
                      j   .ge. m_lo(2) .and. j   .le. m_hi(2) .and. &
                      k+1 .ge. m_lo(3) .and. k+1 .le. m_hi(3)) then

                if (zhi .and. mask(i,j,k+1) > 0) then

                   if (bctype == -1) then
                      bct = tf(i,j,k+1)
                   else
                      bct = bctype
                   end if

                   if (bct == LO_DIRICHLET) then
                      if (bho >= 1) then
                         h2 = 0.5e0_rt * h
                         th2 = 3.e0_rt * h2
                         bfv = 2.e0_rt * beta / ((bcl + h2) * (bcl + th2))
                      else
                         bfv = (beta / h) / (0.5e0_rt * h + bcl)
                      end if
                      bfv = bfv * b(i,j,k+1)
                   else if (bct == LO_NEUMANN) then
                      bfv = beta * r / h
                   else if (bct == LO_MARSHAK .OR. bct == LO_SANCHEZ_POMRANING) then
                      bfv = 2.e0_rt * beta * r / h
#ifndef AMREX_USE_GPU
                   else
                      call castro_error("hbvec3: unsupported boundary type")
#endif
                   end if

                   vec(i,j,k) = vec(i,j,k) + bfv * bcval(i,j,k+1)

                end if

             end if

          end do
       end do
    end do

  end subroutine hbvec3



  subroutine set_abec_flux(lo, hi, &
                           dir, &
                           density, d_lo, d_hi, &
                           dcoef, c_lo, c_hi, &
                           beta, &
                           dx, &
                           flux, f_lo, f_hi) &
                           bind(C, name="set_abec_flux")

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: density(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    real(rt), intent(in   ) :: dcoef(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    integer,  intent(in   ), value :: dir
    real(rt), intent(in   ), value :: beta
    real(rt), intent(in   ) :: dx(3)

    integer  :: i, j, k
    real(rt) :: fac

    !$gpu

    if (dir == 0) then

       fac = -beta / dx(1)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux(i,j,k) = dcoef(i,j,k) * (density(i,j,k) - density(i-1,j,k)) * fac
             end do
          end do
       end do

    else if (dir == 1) then

       fac = -beta / dx(2)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux(i,j,k) = dcoef(i,j,k) * (density(i,j,k) - density(i,j-1,k)) * fac
             end do
          end do
       end do

    else

       fac = -beta / dx(3)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux(i,j,k) = dcoef(i,j,k) * (density(i,j,k) - density(i,j,k-1)) * fac
             end do
          end do
       end do

    end if

  end subroutine set_abec_flux



  subroutine face_metric(i, j, k, lo, hi, dx, ori_dir, ori_lo, r)

    use amrex_constants_module, only: ONE
    use prob_params_module, only: problo, coord_type
    use castro_util_module, only: position ! function

    implicit none

    integer,  intent(in   ) :: i, j, k
    integer,  intent(in   ) :: lo, hi
    integer,  intent(in   ) :: ori_dir, ori_lo
    real(rt), intent(in   ) :: dx
    real(rt), intent(inout) :: r

    real(rt) :: loc(3)

    !$gpu

    if (ori_dir == 1) then

       if (coord_type == 0) then

          r = ONE

       else

          if (ori_lo == 1) then

             r = problo(1) + lo * dx

          else

             r = problo(1) + (hi + 1) * dx

          end if

          if (coord_type == 2) then

             r = r * r

          end if

       end if

    else

       if (coord_type == 0) then

          r = ONE

       else

          ! This must be RZ since we only support spherical coordinates in 1D.

          loc = position(i, j, k)

          r = loc(1)

       end if

    end if

  end subroutine face_metric



  subroutine cell_center_metric(i, j, k, dx, r, s)

    use amrex_constants_module, only: ONE
    use prob_params_module, only: dim, coord_type
    use castro_util_module, only: position ! function

    implicit none

    integer,  intent(in   ) :: i, j, k
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: r, s

    real(rt) :: loc(3)
    real(rt) :: h1, h2, d1, d2
    integer  :: d

    !$gpu

    if (dim >= 2) then
       d = 2
    else
       d = 1
    end if

    if (coord_type == 0) then

       r = ONE
       s = ONE

    else if (coord_type == 1) then

       loc = position(i, j, k)

       r = loc(1)
       s = ONE

    else if (coord_type == 2) then

       loc = position(i, j, k)

       h1 = 0.5e0_rt * dx(1)
       d1 = 1.e0_rt / (3.e0_rt * dx(1))

       r = loc(1)
       r = d1 * ((r + h1)**3 - (r - h1)**3)

       h2 = 0.5e0_rt * dx(2)
       d2 = 1.e0_rt / dx(2)

       s = loc(d)
       s = d2 * (cos(s - h2) - cos(s + h2))

    end if

  end subroutine cell_center_metric



  subroutine edge_center_metric(i, j, k, idir, dx, r, s)

    use amrex_constants_module, only: ONE
    use prob_params_module, only: dim, coord_type
    use castro_util_module, only: position ! function

    implicit none

    integer,  intent(in   ) :: i, j, k, idir
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: r, s

    real(rt) :: loc(3)
    real(rt) :: h1, h2, d1, d2
    integer  :: d

    !$gpu

    if (dim >= 2) then
       d = 2
    else
       d = 1
    end if

    if (coord_type == 0) then

       r = ONE
       s = ONE

    else if (coord_type == 1) then

       if (idir == 1) then
          loc = position(i, j, k, ccx = .false.)
       else
          loc = position(i, j, k)
       end if

       r = loc(1)
       s = ONE

    else if (coord_type == 2) then

       if (idir == 1) then

          loc = position(i, j, k, ccx = .false.)
          r = loc(1)

          loc = position(i, j, k)
          s = loc(d)

          h2 = 0.5e0_rt * dx(2)
          d2 = 1.e0_rt / dx(2)

          r = r**2
          s = d2 * (cos(s - h2) - cos(s + h2))

       else

          loc = position(i, j, k)
          r = loc(1)

          if (d == 2) then
             loc = position(i, j, k, ccy = .false.)
          else
             loc = position(i, j, k, ccx = .false.)
          end if
          s = loc(d)

          h1 = 0.5e0_rt * dx(1)
          d1 = 1.e0_rt / (3.e0_rt * dx(1))

          r = d1 * ((r + h1)**3 - (r - h1)**3)
          s = sin(s)

       end if

    end if

  end subroutine edge_center_metric

end module habec_nd_module
