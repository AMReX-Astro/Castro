module electric_field

implicit none

contains

subroutine electric(Q, E, comp) !Use ideal Ohm's Law

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module, only : NQ, QU,QV, QW, QMAGX, QMAGY, QMAGZ

  implicit none

  real(rt), intent(in)   ::Q(NQ)
  real(rt), intent(out)  ::E
  integer, intent(in)    ::comp

  !E = -v X B
  if(comp.eq. 1) then
     E       = -Q(QV)*Q(QMAGZ) + Q(QW)*Q(QMAGY)
  elseif(comp.eq. 2) then
     E       = -Q(QW)*Q(QMAGX) + Q(QU)*Q(QMAGZ)
  elseif(comp.eq. 3) then
     E       = -Q(QU)*Q(QMAGY) + Q(QV)*Q(QMAGX)
  else
  endif

end subroutine electric

subroutine electric_edge_x(work_lo, work_hi, &
                           q, q_lo, q_hi, &
                           E, ex_lo, ex_hi, &
                           flxy, flxy_lo, flxy_hi, &
                           flxz, flxz_lo, flxz_hi)


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
  integer :: UMAGX, UMAGY, UMAGZ

  UMAGX = NVAR+1
  UMAGY = NVAR+2
  UMAGZ = NVAR+3

  E = 0.d0
  do k = work_lo(3), work_hi(3)
     do j = work_lo(2), work_hi(2)
        do i = work_lo(1), work_hi(1)

           !===================================Ex i, j -1/2, k -1/2 ====================================

           !-------------------Y derivatives -----------------
           !dEx/dy i,j-3/4, k-1/2
           call electric(q(i,j-1,k-1,:),Ecen,1)
           a = 2.d0*( -flxy(i,j,k-1,UMAGZ) - Ecen )
           call electric(q(i,j-1,k,:),Ecen,1)
           b = 2.d0*( -flxy(i,j,k,UMAGZ) - Ecen )

           !Upwind in the z direction to get dEx/dy i, j-3/4, k-1/2
           if(flxz(i,j-1,k,URHO) .gt. 0.d0) then !recall flxz(QRHO) = rho*w so sign(rho*w) = sign(w)
              d1 = a
           else if(flxz(i,j-1,k,URHO) .lt. 0.d0) then
              d1 = b
           else
              d1 = 0.5d0*(a + b)
           endif

           !dEx/dy i,j-1/4, k-1/2
           call electric(q(i,j,k-1,:),Ecen, 1)
           a = 2.d0*( Ecen + flxy(i,j,k-1,UMAGZ) )
           call electric(q(i,j,k,:),Ecen, 1)
           b = 2.d0*( Ecen + flxy(i,j,k,UMAGZ) )

           !Upwind in the z direction to get dEx/dy i, j-1/4, k-1/2
           if(flxz(i,j,k,URHO) .gt. 0.d0) then
              d2 = a
           else if(flxz(i,j,k,URHO) .lt. 0.d0) then
              d2 = b
           else
              d2 = 0.5d0*(a+b)
           endif

           !Calculate the "second derivative" in the y direction for d^2Ex/dy^2 i, j-1/2, k-1/2
           dd1 = 0.125d0*(d1 - d2)

           !------------------ Z derivatives --------------------
           !dEx/dz i, j-1/2, k - 3/4
           call electric(q(i,j-1,k-1,:),Ecen,1)
           a = 2.d0*( flxz(i,j-1,k,UMAGY) - Ecen )
           call electric(q(i,j, k-1, :), Ecen, 1)
           b = 2.d0*( flxz(i,j,k,UMAGY) - Ecen )

           !upwind in the y direction to get dEx/dz i, j-1/2, k -3/4
           if(flxy(i,j,k-1,URHO).gt.0.d0) then
                    d1 = a
           elseif(flxy(i,j,k-1,URHO).lt.0.d0) then
                    d1 = b
           else
                    d1 = 0.5d0*(a + b)
           endif

           !dEx/dz i, j-1/2, k-1/4
           call electric(q(i,j-1,k,:), Ecen, 1)
           a = 2.d0*( Ecen - flxz(i,j-1,k,UMAGY) )
           call electric(q(i,j,k,:), Ecen, 1)
           b = 2.d0*( Ecen - flxz(i,j,k,UMAGY) )

           !Upwind in the y direction for i,j-1/2,k-1/2
           if(flxy(i,j,k,URHO).gt.0.d0) then
              d2 = a
           elseif(flxy(i,j,k,URHO).lt.0.d0) then
              d2 = b
           else
              d2 = 0.5d0*(a + b)
           endif

           !calculate second derivative
           dd2 = 0.125d0*(d1 - d2)

           !----------------------Prescribe Ex i, j -1/2, k -1/2 ------------------------
           E(i,j,k) = 0.25d0*(-flxy(i,j,k,UMAGZ) - flxy(i,j,k-1,UMAGZ) + flxz(i,j-1,k,UMAGY) + flxz(i,j,k,UMAGY)) + dd1 + dd2


        enddo
     enddo
  enddo

end subroutine electric_edge_x

subroutine electric_edge_y(work_lo, work_hi, &
                           q, q_lo, q_hi, &
                           E, ey_lo, ey_hi, &
                           flxx, flxx_lo, flxx_hi, &
                           flxz, flxz_lo, flxz_hi)


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
  integer :: UMAGX, UMAGY, UMAGZ

  UMAGX = NVAR+1
  UMAGY = NVAR+2
  UMAGZ = NVAR+3

  E = 0.d0
  do k = work_lo(3), work_hi(3)
     do j = work_lo(2), work_hi(2)
        do i = work_lo(1), work_hi(1)

           !====================================== Ey i-1/2, j, k-1/2 ============================================
           !-------------------Z derivatives -----------------
           !dEy/dz i-1/2,j, k-3/4
           call electric(q(i-1,j,k-1,:),Ecen,2)
           a = 2.d0*( -flxz(i-1,j,k,UMAGX) - Ecen )
           call electric(q(i,j,k-1,:),Ecen,2)
           b = 2.d0*( -flxz(i,j,k,UMAGX) - Ecen )

           !Upwind in the x direction to get dEy/dz i-1/2, j, k-3/4
           if(flxx(i,j,k-1,URHO) .gt. 0.d0) then
              d1 = a
           else if(flxx(i,j,k-1,URHO) .lt. 0.d0) then
              d1 = b
           else
              d1 = 0.5d0*(a + b)
           endif

           !dEy/dz i-1/2,j, k-1/4
           call electric(q(i-1,j,k,:),Ecen, 2)
           a = 2.d0*( Ecen + flxz(i-1,j,k,UMAGX) )
           call electric(q(i,j,k,:),Ecen, 2)
           b = 2.d0*( Ecen + flxz(i,j,k,UMAGX) )

           !Upwind in the x direction to get dEy/dz i-1/2, j, k-1/4
           if(flxx(i,j,k,URHO) .gt. 0.d0) then
              d2 = a
           else if(flxx(i,j,k,URHO) .lt. 0.d0) then
              d2 = b
           else
              d2 = 0.5d0*(a+b)
           endif

           !Calculate the "second derivative" in the y direction for d^2Ey/dz^2 i-1/2, j, k-1/2
           dd1 = 0.125d0*(d1-d2)

         !------------------ X derivatives --------------------
           !dEy/dx i -3/4, j, k - 1/2
           call electric(q(i-1,j,k-1,:),Ecen,2)
           a = 2.d0*( flxx(i,j,k-1,UMAGZ) - Ecen )
           call electric(q(i-1,j, k, :), Ecen, 2)
           b = 2.d0*( flxx(i,j,k,UMAGZ) - Ecen )

           !upwind in the z direction to get dEy/dx i-3/4, j, k -1/2
           if(flxz(i-1,j,k,URHO).gt.0.d0) then
                    d1 = a
           elseif(flxz(i-1,j,k,URHO).lt.0.d0) then
                    d1 = b
           else
                    d1 = 0.5d0*(a + b)
           endif

           !dEy/dx i-1/4, j, k-1/2
           call electric(q(i,j,k-1,:), Ecen, 2)
           a = 2.d0*( Ecen - flxx(i,j,k-1,UMAGZ) )
           call electric(q(i,j,k,:), Ecen, 2)
           b = 2.d0*( Ecen - flxx(i,j,k,UMAGZ) )

           !Upwind in the z direction for i-1/2,j,k-1/2
           if(flxz(i,j,k,URHO).gt.0.d0) then
              d2 = a
           elseif(flxz(i,j,k,URHO).lt.0.d0) then
              d2 = b
           else
              d2 = 0.5d0*(a + b)
           endif

           !calculate second derivative
           dd2 = 0.125d0*(d1-d2)

           !----------------------Prescribe Ex i-1/2, j , k -1/2 ------------------------
           E(i,j,k) = 0.25d0*(-flxz(i,j,k,UMAGX) - flxz(i-1,j,k,UMAGX) + flxx(i,j,k-1,UMAGZ) + flxx(i,j,k,UMAGZ)) + dd1 + dd2


        enddo
     enddo
  enddo

end subroutine electric_edge_y

subroutine electric_edge_z(work_lo, work_hi, &
                           q, q_lo, q_hi, &
                           E, ez_lo, ez_hi, &
                           flxx, flxx_lo, flxx_hi, &
                           flxy, flxy_lo, flxy_hi)

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
  integer   :: UMAGX, UMAGY, UMAGZ

  UMAGX = NVAR+1
  UMAGY = NVAR+2
  UMAGZ = NVAR+3

  E = 0.d0

  do k = work_lo(3), work_hi(3)
     do j = work_lo(2), work_hi(2)
        do i = work_lo(1), work_hi(1)

           !===================================Ez i- 1/2, j -1/2, k  ====================================

           !-------------------X derivatives -----------------
           !dEz/dx i-3/4 ,j-1/2, k
           call electric(q(i-1,j-1,k,:),Ecen,3)
           a = 2.d0*( -flxx(i,j-1,k,UMAGY) - Ecen )
           call electric(q(i-1,j,k,:),Ecen,3)
           b = 2.d0*( -flxx(i,j,k,UMAGY) - Ecen )

           !Upwind in the y direction to get dEz/dx i-3/4, j-1/2, k
           if(flxy(i-1,j,k,URHO) .gt. 0.d0) then
              d1 = a
           else if(flxy(i-1,j,k,URHO) .lt. 0.d0) then
              d1 = b
           else
              d1 = 0.5d0*(a + b)
           endif

           !dEz/dx i-1/4,j-1/2, k
           call electric(q(i,j-1,k,:),Ecen, 3)
           a = 2.d0*( Ecen + flxx(i,j-1,k,UMAGY) )
           call electric(q(i,j,k,:),Ecen, 3)
           b = 2.d0*( Ecen + flxx(i,j,k,UMAGY) )

           !Upwind in the y direction to get dEz/dx i-1/4, j-1/2, k
           if(flxy(i,j,k,URHO) .gt. 0.d0) then
              d2 = a
           else if(flxy(i,j,k,URHO) .lt. 0.d0) then
              d2 = b
           else
              d2 = 0.5d0*(a+b)
           endif

           !Calculate the "second derivative" in the x direction for d^2Ez/dx^2 i-1/2, j-1/2, k
           dd1 = 0.125d0*(d1-d2)

           !------------------ Y derivatives --------------------
           !dEz/dy i-1/2, j-3/4, k
           call electric(q(i-1,j-1,k,:),Ecen,3)
           a = 2.d0*( flxy(i-1,j,k,UMAGX) - Ecen )
           call electric(q(i,j-1, k, :), Ecen, 3)
           b = 2.d0*( flxy(i,j,k,UMAGX) - Ecen )

           !upwind in the x direction to get dEz/dy i-1/2, j-3/4, k
           if(flxx(i,j-1,k,URHO).gt.0.d0) then
                    d1 = a
           elseif(flxx(i,j-1,k,URHO).lt.0.d0) then
                    d1 = b
           else
                    d1 = 0.5d0*(a + b)
           endif

           !dEz/dy i-1/2, j-1/4, k
           call electric(q(i-1,j,k,:), Ecen, 3)
           a = 2.d0*( Ecen - flxy(i-1,j,k,UMAGX) )
           call electric(q(i,j,k,:), Ecen, 3)
           b = 2.d0*( Ecen - flxy(i,j,k,UMAGX) )

           !Upwind in the x direction for i-1/2,j-1/2,k
           if(flxx(i,j,k,URHO).gt.0.d0) then
              d2 = a
           elseif(flxx(i,j,k,URHO).lt.0.d0) then
              d2 = b
           else
              d2 = 0.5d0*(a + b)
           endif

           !calculate second derivative
           dd2 = 0.125d0*(d1-d2)

           !----------------------Prescribe Ez i -1/2, j -1/2, k ------------------------
           E(i,j,k) = 0.25d0*(-flxx(i,j,k,UMAGY) - flxx(i,j-1,k,UMAGY) + flxy(i-1,j,k,UMAGX) + flxy(i,j,k,UMAGX)) + dd1 + dd2



        enddo
     enddo
  enddo


end subroutine electric_edge_z

end module electric_field
