module electric_field

  use mhd_state_module
  implicit none

contains

  subroutine electric(Q, E, comp)

    ! this takes the cell-center primitive state, Q, and computes the cell-center
    ! electric field, E, using Faraday's law:
    ! E = -v X B

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NQ, QU,QV, QW, QMAGX, QMAGY, QMAGZ

    implicit none

    real(rt), intent(in)   :: Q(NQ)
    real(rt), intent(out)  ::E
    integer, intent(in)    ::comp

    if(comp == 1) then
       E = -Q(QV) * Q(QMAGZ) + Q(QW) * Q(QMAGY)
    else if (comp == 2) then
       E = -Q(QW) * Q(QMAGX) + Q(QU) * Q(QMAGZ)
    else if (comp == 3) then
       E = -Q(QU) * Q(QMAGY) + Q(QV) * Q(QMAGX)
    end if

  end subroutine electric

subroutine electric_edge_x(work_lo, work_hi, &
                           q, q_lo, q_hi, &
                           E, ex_lo, ex_hi, &
                           flxy, flxy_lo, flxy_hi, &
                           flxz, flxz_lo, flxz_hi) bind(C, name="electric_edge_x")

  ! Compute Ex on an edge.  This will compute Ex(i, j-1/2, k-1/2)

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module

  implicit none
  integer, intent(in) :: work_lo(3), work_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: ex_lo(3), ex_hi(3)

  integer, intent(in) :: flxy_lo(3), flxy_hi(3)
  integer, intent(in) :: flxz_lo(3), flxz_hi(3)

  real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

  real(rt), intent(in) :: flxy(flxy_lo(1):flxy_hi(1),flxy_lo(2):flxy_hi(2),flxy_lo(3):flxy_hi(3),NVAR+3)
  real(rt), intent(in) :: flxz(flxz_lo(1):flxz_hi(1),flxz_lo(2):flxz_hi(2),flxz_lo(3):flxz_hi(3),NVAR+3)

  real(rt), intent(out) :: E(ex_lo(1):ex_hi(1),ex_lo(2):ex_hi(2),ex_lo(3):ex_hi(3))

  real(rt) :: Ecen
  real(rt) :: a ,b ,d1 ,d2 ,dd1 ,dd2

  integer :: i ,j ,k

  do k = work_lo(3), work_hi(3)
     do j = work_lo(2), work_hi(2)
        do i = work_lo(1), work_hi(1)

           ! Compute Ex(i, j-1/2, k-1/2) using MM Eq. 50

           ! dEx/dy (Eq. 49), located at (i, j-3/4, k-1/2)

           ! first compute dEx/dy_{i,j-3/4,k-1} using MM Eq. 49
           ! note that the face value Ex_{i,j-1/2,k-1} = -F_{i,j-1/2,k-1}(Bz) 
           ! via Faraday's law (MM Eq. 15)
           call electric(q(i,j-1,k-1,:), Ecen,1)
           a = 2.d0*( -flxy(i,j,k-1,UMAGZ) - Ecen )

           ! now compute dEx/dy_{i,j-3/4,k}
           call electric(q(i,j-1,k,:), Ecen,1)
           b = 2.d0*( -flxy(i,j,k,UMAGZ) - Ecen )

           ! Upwind in the z direction to get dEx/dy i, j-3/4, k-1/2
           ! using w_{i,j-1,k-1/2}
           ! recall flxz(QRHO) = rho*w so sign(rho*w) = sign(w)
           if (flxz(i,j-1,k,URHO) .gt. 0.d0) then 
              d1 = a
           else if (flxz(i,j-1,k,URHO) .lt. 0.d0) then
              d1 = b
           else
              d1 = 0.5d0*(a + b)
           endif

           ! dEx/dy located at (i, j-1/4, k-1/2)

           ! first compute dEx/dy_{i,j-1/4,k-1}
           call electric(q(i,j,k-1,:), Ecen, 1)
           a = 2.d0 * (Ecen + flxy(i,j,k-1,UMAGZ) )

           ! now compute dEx/dy_{i,j-1/4,k}
           call electric(q(i,j,k,:), Ecen, 1)
           b = 2.d0 * (Ecen + flxy(i,j,k,UMAGZ) )

           ! finally upwind in the z direction to get dEx/dy i, j-1/4, k-1/2
           ! using w_{i,j,k-1/2}
           if (flxz(i,j,k,URHO) .gt. 0.d0) then
              d2 = a
           else if (flxz(i,j,k,URHO) .lt. 0.d0) then
              d2 = b
           else
              d2 = 0.5d0*(a + b)
           endif

           ! Calculate the "second derivative" in the y direction for
           ! d^2Ex/dy^2 i, j-1/2, k-1/2 (this is one of the terms in Eq. 50)
           ! note: Stone 08 Eq. 79 has the signs backwards for this term.
           dd1 = 0.125d0*(d1 - d2)


           ! now dEx/dz located at (i, j-1/2, k-3/4)

           ! first compute dEx/dz_{i,j-1,k-3/4}
           ! note that the face value of Ex_{i,j-1,k-1/2} = F_{i,j-1,k-1/2}(By)
           call electric(q(i,j-1,k-1,:), Ecen, 1)
           a = 2.d0 * (flxz(i,j-1,k,UMAGY) - Ecen )

           ! now compute dEx/dz_{i,j,k-3/4}
           call electric(q(i,j, k-1, :), Ecen, 1)
           b = 2.d0 * (flxz(i,j,k,UMAGY) - Ecen)

           ! upwind in the y direction to get dEx/dz i, j-1/2, k-3/4
           ! using v_{i,j-1/2,k-1}
           if (flxy(i,j,k-1,URHO) .gt. 0.d0) then
              d1 = a
           else if (flxy(i,j,k-1,URHO) .lt. 0.d0) then
              d1 = b
           else
              d1 = 0.5d0*(a + b)
           endif

           ! dEx/dz located at (i, j-1/2, k-1/4)

           ! first compute dEx/dz_{i,j-1,k-1/4}
           call electric(q(i,j-1,k,:), Ecen, 1)
           a = 2.d0 * (Ecen - flxz(i,j-1,k,UMAGY) )

           ! now compute dEx/dz_{i,j,k-1/4}
           call electric(q(i,j,k,:), Ecen, 1)
           b = 2.d0 * (Ecen - flxz(i,j,k,UMAGY) )

           ! upwind in the y direction to get dEx/dz i, j-1/2, k-1/4
           ! using v_{i,j-1/2,k}
           if (flxy(i,j,k,URHO) .gt. 0.d0) then
              d2 = a
           else if (flxy(i,j,k,URHO) .lt. 0.d0) then
              d2 = b
           else
              d2 = 0.5d0*(a + b)
           endif

           ! calculate second derivative
           dd2 = 0.125d0*(d1 - d2)

           ! now the final Ex_{i,j-1/2,k-1/2}, using MM Eq. 50 (shifted to j-1/2, k-1/2)

           E(i,j,k) = 0.25d0*(-flxy(i,j,k,UMAGZ) - flxy(i,j,k-1,UMAGZ) + &
                              flxz(i,j-1,k,UMAGY) + flxz(i,j,k,UMAGY)) + dd1 + dd2

        enddo
     enddo
  enddo

end subroutine electric_edge_x

subroutine electric_edge_y(work_lo, work_hi, &
                           q, q_lo, q_hi, &
                           E, ey_lo, ey_hi, &
                           flxx, flxx_lo, flxx_hi, &
                           flxz, flxz_lo, flxz_hi) bind(C, name="electric_edge_y")


  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module

  implicit none

  integer, intent(in) :: work_lo(3), work_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: ey_lo(3), ey_hi(3)
  integer, intent(in) :: flxx_lo(3), flxx_hi(3)
  integer, intent(in) :: flxz_lo(3), flxz_hi(3)

  real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

  real(rt), intent(in) :: flxx(flxx_lo(1):flxx_hi(1),flxx_lo(2):flxx_hi(2),flxx_lo(3):flxx_hi(3),NVAR+3)
  real(rt), intent(in) :: flxz(flxz_lo(1):flxz_hi(1),flxz_lo(2):flxz_hi(2),flxz_lo(3):flxz_hi(3),NVAR+3)

  real(rt), intent(out) :: E(ey_lo(1):ey_hi(1),ey_lo(2):ey_hi(2),ey_lo(3):ey_hi(3))

  real(rt)  :: Ecen
  real(rt)  :: a ,b ,d1 ,d2 ,dd1 ,dd2

  integer   :: i ,j ,k

  do k = work_lo(3), work_hi(3)
     do j = work_lo(2), work_hi(2)
        do i = work_lo(1), work_hi(1)

           ! Compute Ey(i-1/2, j, k-1/2)

           ! dEy/dz i-1/2, j, k-3/4

           ! first compute dEy/dz_{i-1,j,k-3/4}
           call electric(q(i-1,j,k-1,:), Ecen, 2)
           a = 2.d0 * (-flxz(i-1,j,k,UMAGX) - Ecen)

           ! now compute dEy/dz_{i,j,k-3/4}
           call electric(q(i,j,k-1,:), Ecen, 2)
           b = 2.d0 * (-flxz(i,j,k,UMAGX) - Ecen)

           ! upwind in the x direction to get dEy/dz i-1/2, j, k-3/4
           ! using u_{i-1/2,j,k-1}
           if (flxx(i,j,k-1,URHO) .gt. 0.d0) then
              d1 = a
           else if (flxx(i,j,k-1,URHO) .lt. 0.d0) then
              d1 = b
           else
              d1 = 0.5d0*(a + b)
           endif

           ! dEy/dz i-1/2, j, k-1/4

           ! first compute dEy/dz_{i-1,j,k-1/4}
           call electric(q(i-1,j,k,:), Ecen, 2)
           a = 2.d0 * (Ecen + flxz(i-1,j,k,UMAGX))

           ! now compute dEy/dz_{i,j,k-1/4}
           call electric(q(i,j,k,:), Ecen, 2)
           b = 2.d0 * (Ecen + flxz(i,j,k,UMAGX))

           ! upwind in the x direction to get dEy/dz i-1/2, j, k-1/4
           ! using u_{i-1/2.j,k}
           if (flxx(i,j,k,URHO) .gt. 0.d0) then
              d2 = a
           else if (flxx(i,j,k,URHO) .lt. 0.d0) then
              d2 = b
           else
              d2 = 0.5d0*(a+b)
           endif

           ! calculate the "second derivative" in the y direction for
           ! d^2Ey/dz^2 i-1/2, j, k-1/2
           dd1 = 0.125d0*(d1-d2)


           ! dEy/dx i-3/4, j, k-1/2

           ! first compute dEy/dz_{i-3/4,j,k-1}
           call electric(q(i-1,j,k-1,:), Ecen, 2)
           a = 2.d0*(flxx(i,j,k-1,UMAGZ) - Ecen)

           ! next compute dEy/dz_{i-3/4,j,k}
           call electric(q(i-1,j,k, :), Ecen, 2)
           b = 2.d0*(flxx(i,j,k,UMAGZ) - Ecen)

           ! upwind in the z direction to get dEy/dx i-3/4, j, k-1/2
           ! using w_{i-1,j,k-1/2}
           if (flxz(i-1,j,k,URHO) .gt. 0.d0) then
              d1 = a
           else if (flxz(i-1,j,k,URHO) .lt. 0.d0) then
              d1 = b
           else
              d1 = 0.5d0*(a + b)
           endif

           ! dEy/dx i-1/4, j, k-1/2

           ! first compute dEy/dx_{i-1/4,j,k-1}
           call electric(q(i,j,k-1,:), Ecen, 2)
           a = 2.d0*(Ecen - flxx(i,j,k-1,UMAGZ))

           ! next compute dEy/dx_{i-1/4,j,k}
           call electric(q(i,j,k,:), Ecen, 2)
           b = 2.d0*(Ecen - flxx(i,j,k,UMAGZ))

           ! upwind in the z direction for i-1/4, j, k-1/2
           ! using w_{i,j,k-1/2}
           if (flxz(i,j,k,URHO) .gt. 0.d0) then
              d2 = a
           else if (flxz(i,j,k,URHO) .lt. 0.d0) then
              d2 = b
           else
              d2 = 0.5d0*(a + b)
           endif

           ! calculate second derivative
           dd2 = 0.125d0*(d1-d2)

           ! now the final Ey_{i-1/2, j, k-1/2}
           E(i,j,k) = 0.25d0*(-flxz(i,j,k,UMAGX) - flxz(i-1,j,k,UMAGX) + &
                              flxx(i,j,k-1,UMAGZ) + flxx(i,j,k,UMAGZ)) + dd1 + dd2

        enddo
     enddo
  enddo

end subroutine electric_edge_y

subroutine electric_edge_z(work_lo, work_hi, &
                           q, q_lo, q_hi, &
                           E, ez_lo, ez_hi, &
                           flxx, flxx_lo, flxx_hi, &
                           flxy, flxy_lo, flxy_hi) bind(C, name="electric_edge_z")

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module

  implicit none
  integer, intent(in)   :: work_lo(3), work_hi(3)
  integer, intent(in)   :: q_lo(3), q_hi(3)
  integer, intent(in)   :: ez_lo(3), ez_hi(3)

  integer, intent(in)   :: flxx_lo(3), flxx_hi(3)
  integer, intent(in)   :: flxy_lo(3), flxy_hi(3)

  real(rt), intent(in)  :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

  real(rt), intent(in)  :: flxx(flxx_lo(1):flxx_hi(1),flxx_lo(2):flxx_hi(2),flxx_lo(3):flxx_hi(3),NVAR+3)
  real(rt), intent(in)  :: flxy(flxy_lo(1):flxy_hi(1),flxy_lo(2):flxy_hi(2),flxy_lo(3):flxy_hi(3),NVAR+3)

  real(rt), intent(out)   :: E(ez_lo(1):ez_hi(1),ez_lo(2):ez_hi(2),ez_lo(3):ez_hi(3))

  real(rt)  :: Ecen, u_face, v_face
  real(rt)  :: a ,b ,d1 ,d2 ,dd1 ,dd2
  integer   :: i ,j ,k

  do k = work_lo(3), work_hi(3)
     do j = work_lo(2), work_hi(2)
        do i = work_lo(1), work_hi(1)

           ! Compute Ez(i-1/2, j-1/2, k)


           ! dEz/dx i-3/4, j-1/2, k

           ! first compute dEz/dx_{i-3/4,j-1,k}
           call electric(q(i-1,j-1,k,:), Ecen, 3)
           a = 2.d0 * (-flxx(i,j-1,k,UMAGY) - Ecen)

           ! next dEz/dx_{i-3/4,j,k}
           call electric(q(i-1,j,k,:), Ecen, 3)
           b = 2.d0 * (-flxx(i,j,k,UMAGY) - Ecen)

           ! upwind in the y direction to get dEz/dx i-3/4, j-1/2, k
           ! using v_{i-1,j-1/2,k}
           if ( flxy(i-1,j,k,URHO) .gt. 0.d0) then
              d1 = a
           else if (flxy(i-1,j,k,URHO) .lt. 0.d0) then
              d1 = b
           else
              d1 = 0.5d0*(a + b)
           endif

           ! dEz/dx i-1/4, j-1/2, k

           ! first compute dEz/dx_{i-1/4,j-1,k}
           call electric(q(i,j-1,k,:), Ecen, 3)
           a = 2.d0 * (Ecen + flxx(i,j-1,k,UMAGY))

           ! next dEz/dx_{i-1/4,j,k}
           call electric(q(i,j,k,:), Ecen, 3)
           b = 2.d0 * (Ecen + flxx(i,j,k,UMAGY))

           ! upwind in the y direction to get dEz/dx i-1/4, j-1/2, k
           ! using v_{i,j-1/2,k}
           if (flxy(i,j,k,URHO) .gt. 0.d0) then
              d2 = a
           else if (flxy(i,j,k,URHO) .lt. 0.d0) then
              d2 = b
           else
              d2 = 0.5d0*(a+b)
           endif

           ! Calculate the "second derivative" in the x direction for
           ! d^2Ez/dx^2 i-1/2, j-1/2, k
           dd1 = 0.125d0*(d1-d2)


           ! dEz/dy i-1/2, j-3/4, k

           ! first compute dEz/dy_{i-1,j-3/4,k}
           call electric(q(i-1,j-1,k,:), Ecen, 3)
           a = 2.d0 * (flxy(i-1,j,k,UMAGX) - Ecen)

           ! now compute dEz/dy_{i,j-3/4,k}
           call electric(q(i,j-1, k, :), Ecen, 3)
           b = 2.d0 * (flxy(i,j,k,UMAGX) - Ecen)

           ! upwind in the x direction to get dEz/dy i-1/2, j-3/4, k
           ! using u_{i-1/2,j-1,k}
           if (flxx(i,j-1,k,URHO) .gt. 0.d0) then
              d1 = a
           else if (flxx(i,j-1,k,URHO) .lt. 0.d0) then
              d1 = b
           else
              d1 = 0.5d0*(a + b)
           endif

           ! dEz/dy i-1/2, j-1/4, k

           ! first compute dEz/dy_{i-1,j-1/4,k}
           call electric(q(i-1,j,k,:), Ecen, 3)
           a = 2.d0 * (Ecen - flxy(i-1,j,k,UMAGX))

           ! now compute dEz/dy_{i,j-1/4,k}
           call electric(q(i,j,k,:), Ecen, 3)
           b = 2.d0 * (Ecen - flxy(i,j,k,UMAGX)) 

           ! Upwind in the x direction for i-1/2, j-1/4, k
           ! using u_{i-1/2,j,k}
           if (flxx(i,j,k,URHO) .gt. 0.d0) then
              d2 = a
           else if (flxx(i,j,k,URHO) .lt. 0.d0) then
              d2 = b
           else
              d2 = 0.5d0*(a + b)
           endif

           ! calculate second derivative
           dd2 = 0.125d0*(d1-d2)

           ! compute Ez i-1/2, j-1/2, k
           E(i,j,k) = 0.25d0*(-flxx(i,j,k,UMAGY) - flxx(i,j-1,k,UMAGY) + &
                              flxy(i-1,j,k,UMAGX) + flxy(i,j,k,UMAGX)) + dd1 + dd2

        enddo
     enddo
  enddo


end subroutine electric_edge_z

end module electric_field
