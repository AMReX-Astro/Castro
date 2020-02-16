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
                           q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                           E, ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                           flxy, flxy_l1 , flxy_l2 , flxy_l3 , flxy_h1 , flxy_h2 , flxy_h3, &
                           flxz, flxz_l1 , flxz_l2 , flxz_l3 , flxz_h1 , flxz_h2 , flxz_h3)


  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module

  implicit none
  integer, intent(in)   :: work_lo(3), work_hi(3)
  integer, intent(in)   :: q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
  integer, intent(in)   :: ex_l1,ex_l2,ex_l3,ex_h1,ex_h2, ex_h3

  integer, intent(in)   :: flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3
  integer, intent(in)   :: flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3

  real(rt), intent(in)  ::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ)

  real(rt), intent(in) :: flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,NVAR+3) 
  real(rt), intent(in) :: flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,NVAR+3)

  real(rt), intent(out) :: E(ex_l1:ex_h1,ex_l2:ex_h2,ex_l3:ex_h3)

  real(rt)   ::Ecen
  real(rt)   ::a ,b ,d1 ,d2 ,dd1 ,dd2 

  integer    ::i ,j ,k
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
           a = 2.d0*( -flxy(i,j,k-1,UMAGZ) - Ecen )   !(Ecen + flxy(i,j,k-1,QMAGZ))
           call electric(q(i,j-1,k,:),Ecen,1)
           b = 2.d0*( -flxy(i,j,k,UMAGZ) - Ecen )  !(Ecen + flxy(i,j,k,QMAGZ))
        
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
           a = 2.d0*( Ecen + flxy(i,j,k-1,UMAGZ) ) !(-flxy(i,j,k-1,QMAGZ) - Ecen)
           call electric(q(i,j,k,:),Ecen, 1)
           b = 2.d0*( Ecen + flxy(i,j,k,UMAGZ) )  !(-flxy(i,j,k,QMAGZ) - Ecen)

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
           a = 2.d0*( flxz(i,j-1,k,UMAGY) - Ecen )       !(Ecen - flxz(i,j-1,k,UMAGY))   
           call electric(q(i,j, k-1, :), Ecen, 1)
           b = 2.d0*( flxz(i,j,k,UMAGY) - Ecen )         !(Ecen - flxz(i,j,k,UMAGY))
           
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
           a = 2.d0*( Ecen - flxz(i,j-1,k,UMAGY) )      !(flxz(i,j-1,k,QMAGY) - Ecen)
           call electric(q(i,j,k,:), Ecen, 1)
           b = 2.d0*( Ecen - flxz(i,j,k,UMAGY) )  !(flxz(i,j,k,QMAGY) - Ecen)

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
                           q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                           E, ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                           flxx, flxx_l1 , flxx_l2 , flxx_l3 , flxx_h1 , flxx_h2 , flxx_h3, &
                           flxz, flxz_l1 , flxz_l2 , flxz_l3 , flxz_h1 , flxz_h2 , flxz_h3)


  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module

  implicit none

  integer, intent(in)   :: work_lo(3), work_hi(3)
  integer, intent(in)   :: q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
  integer, intent(in)   :: ey_l1,ey_l2,ey_l3,ey_h1,ey_h2, ey_h3
  integer, intent(in)   :: flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3
  integer, intent(in)   :: flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3

  real(rt), intent(in)  ::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ)

  real(rt), intent(in)  :: flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,NVAR+3)
  real(rt), intent(in)  :: flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,NVAR+3)

  real(rt), intent(out) :: E(ey_l1:ey_h1,ey_l2:ey_h2,ey_l3:ey_h3)

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
           a = 2.d0*( -flxz(i-1,j,k,UMAGX) - Ecen )   !(Ecen + flxy(i,j,k-1,QMAGZ))
           call electric(q(i,j,k-1,:),Ecen,2)
           b = 2.d0*( -flxz(i,j,k,UMAGX) - Ecen )  !(Ecen + flxy(i,j,k,QMAGZ))
        
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
           a = 2.d0*( Ecen + flxz(i-1,j,k,UMAGX) ) !(-flxy(i,j,k-1,QMAGZ) - Ecen)
           call electric(q(i,j,k,:),Ecen, 2)
           b = 2.d0*( Ecen + flxz(i,j,k,UMAGX) )  !(-flxy(i,j,k,QMAGZ) - Ecen)

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
           a = 2.d0*( flxx(i,j,k-1,UMAGZ) - Ecen )       !(Ecen - flxz(i,j-1,k,UMAGY))   
           call electric(q(i-1,j, k, :), Ecen, 2)
           b = 2.d0*( flxx(i,j,k,UMAGZ) - Ecen )         !(Ecen - flxz(i,j,k,UMAGY))
           
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
           a = 2.d0*( Ecen - flxx(i,j,k-1,UMAGZ) )      !(flxz(i,j-1,k,QMAGY) - Ecen)
           call electric(q(i,j,k,:), Ecen, 2)
           b = 2.d0*( Ecen - flxx(i,j,k,UMAGZ) )  !(flxz(i,j,k,QMAGY) - Ecen)

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
                           q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                           E, ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                           flxx, flxx_l1 , flxx_l2 , flxx_l3 , flxx_h1 , flxx_h2 , flxx_h3, &
                           flxy, flxy_l1 , flxy_l2 , flxy_l3 , flxy_h1 , flxy_h2 , flxy_h3)

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module

  implicit none
  integer, intent(in)   :: work_lo(3), work_hi(3)
  integer, intent(in)   :: q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
  integer, intent(in)   :: ez_l1,ez_l2,ez_l3,ez_h1,ez_h2, ez_h3

  integer, intent(in)   :: flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3
  integer, intent(in)   :: flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3

  real(rt), intent(in)  :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ)

  real(rt), intent(in)  :: flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,NVAR+3)
  real(rt), intent(in)  :: flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,NVAR+3)

  real(rt), intent(out)   :: E(ez_l1:ez_h1,ez_l2:ez_h2,ez_l3:ez_h3)

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
           a = 2.d0*( -flxx(i,j-1,k,UMAGY) - Ecen )   !(Ecen + flxy(i,j,k-1,QMAGZ))
           call electric(q(i-1,j,k,:),Ecen,3)
           b = 2.d0*( -flxx(i,j,k,UMAGY) - Ecen )  !(Ecen + flxy(i,j,k,QMAGZ))
        
           !Upwind in the y direction to get dEz/dx i-3/4, j-1/2, k
           if(flxy(i-1,j,k,URHO) .gt. 0.d0) then !recall flxz(QRHO) = rho*w so sign(rho*w) = sign(w)
              d1 = a
           else if(flxy(i-1,j,k,URHO) .lt. 0.d0) then
              d1 = b
           else
              d1 = 0.5d0*(a + b)
           endif

           !dEz/dx i-1/4,j-1/2, k
           call electric(q(i,j-1,k,:),Ecen, 3)
           a = 2.d0*( Ecen + flxx(i,j-1,k,UMAGY) ) !(-flxy(i,j,k-1,QMAGZ) - Ecen)
           call electric(q(i,j,k,:),Ecen, 3)
           b = 2.d0*( Ecen + flxx(i,j,k,UMAGY) )  !(-flxy(i,j,k,QMAGZ) - Ecen)

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
           a = 2.d0*( flxy(i-1,j,k,UMAGX) - Ecen )       !(Ecen - flxz(i,j-1,k,UMAGY))   
           call electric(q(i,j-1, k, :), Ecen, 3)
           b = 2.d0*( flxy(i,j,k,UMAGX) - Ecen )         !(Ecen - flxz(i,j,k,UMAGY))
           
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
           a = 2.d0*( Ecen - flxy(i-1,j,k,UMAGX) )      !(flxz(i,j-1,k,QMAGY) - Ecen)
           call electric(q(i,j,k,:), Ecen, 3)
           b = 2.d0*( Ecen - flxy(i,j,k,UMAGX) )  !(flxz(i,j,k,QMAGY) - Ecen)

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
