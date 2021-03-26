module castro_util_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine ca_find_center(data,new_center,icen,dx,problo) &
       bind(C, name="ca_find_center")
    !
    ! .. note::
    !    Binds to C function ``ca_find_center``

    use amrex_constants_module, only: ZERO, HALF, TWO
    use prob_params_module, only: dg, dim
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(inout) :: data(-1:1,-1*dg(2):1*dg(2),-1*dg(3):1*dg(3))
    real(rt), intent(  out) :: new_center(3)
    real(rt), intent(in   ) :: dx(3),problo(3)

    real(rt) :: a,b,x,y,z,cen
    integer  :: icen(3)
    integer  :: i,j,k

    if (dim .eq. 1) then

       ! In 1-D it only make sense to have the center at the origin
       new_center = ZERO

    else if (dim .ge. 2) then

       ! We do this to take care of precision issues
       cen = data(0,0,0)
       do k = -1*dg(3),1*dg(3)
          do j = -1*dg(2),1*dg(2)
             do i = -1*dg(1),1*dg(1)
                data(i,j,k) = data(i,j,k) - cen
             end do
          end do
       end do

       ! This puts the "center" at the cell center
       new_center(1:dim) = problo(1:dim) +  (icen(1:dim)+HALF) * dx(1:dim)

       ! Fit parabola y = a x^2  + b x + c through three points
       ! a = 1/2 ( y_1 + y_-1)
       ! b = 1/2 ( y_1 - y_-1)
       ! x_vertex = -b / 2a

       ! ... in x-direction
       a = HALF * (data(1,0,0) + data(-1,0,0)) - data(0,0,0)
       b = HALF * (data(1,0,0) - data(-1,0,0)) - data(0,0,0)
       x = -b / (TWO*a)
       new_center(1) = new_center(1) +  x*dx(1)

       ! ... in y-direction
       a = HALF * (data(0,1,0) + data(0,-1,0)) - data(0,0,0)
       b = HALF * (data(0,1,0) - data(0,-1,0)) - data(0,0,0)
       y = -b / (TWO*a)
       new_center(2) = new_center(2) +  y*dx(2)

       if (dim .eq. 3) then

          ! ... in z-direction
          a = HALF * (data(0,0,1) + data(0,0,-1)) - data(0,0,0)
          b = HALF * (data(0,0,1) - data(0,0,-1)) - data(0,0,0)
          z = -b / (TWO*a)
          new_center(3) = new_center(3) +  z*dx(3)

       endif

    endif

  end subroutine ca_find_center

end module castro_util_module
