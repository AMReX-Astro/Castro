#include <AMReX_LO_BCTYPES.H>

module habec_nd_module

  ! habec is Hypre abec, where abec is the form of the linear equation
  ! we are solving:
  ! 
  ! alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

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

#if AMREX_SPACEDIM >= 2
       h2 = 0.5e0_rt * dx(2)
       d2 = 1.e0_rt / dx(2)

       s = loc(d)
       s = d2 * (cos(s - h2) - cos(s + h2))
#else
       s = ONE
#endif

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
          r = r**2

#if AMREX_SPACEDIM >= 2
          loc = position(i, j, k)
          s = loc(d)

          h2 = 0.5e0_rt * dx(2)
          d2 = 1.e0_rt / dx(2)

          s = d2 * (cos(s - h2) - cos(s + h2))
#else
          s = ONE
#endif

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
