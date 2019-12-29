module habec_nd_module

  ! habec is Hypre abec, where abec is the form of the linear equation
  ! we are solving:
  ! 
  ! alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

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
