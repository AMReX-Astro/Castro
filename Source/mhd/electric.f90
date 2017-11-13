module electric_field

implicit none

contains

subroutine electric(Q, E, comp) !Use ideal Ohm's Law
use amrex_fort_module, only : rt => amrex_real
use meth_params_module, only : QVAR, QU,QV, QW, QMAGX, QMAGY, QMAGZ

implicit none

 real(rt), intent(in)	::Q(QVAR)
 real(rt), intent(out) 	::E
 integer, intent(in)    ::comp

		!E = -v X B
	if(comp.eq. 1) then
	E	= -Q(QV)*Q(QMAGZ) + Q(QW)*Q(QMAGY)
	elseif(comp.eq. 2) then
	E	= -Q(QW)*Q(QMAGX) + Q(QU)*Q(QMAGZ)
	elseif(comp.eq. 3) then
	E	= -Q(QU)*Q(QMAGY) + Q(QV)*Q(QMAGX)
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

	real(rt), intent(in)	::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)

        real(rt), intent(in) :: flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,QVAR) 
        real(rt), intent(in) :: flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,QVAR)

	real(rt), intent(out)	:: E(ex_l1:ex_h1,ex_l2:ex_h2,ex_l3:ex_h3)
	
	real(rt)				::Ecen
	real(rt)				::a ,b ,d1 ,d2 ,dd1 ,dd2 
	
	integer					::i ,j ,k	

	!E = 0.d0
	do k = work_lo(3), work_hi(3)
		do j = work_lo(2), work_hi(2)
			do i = work_lo(1), work_hi(1)

!============================================= Ex i, j -1/2, k -1/2 ===========================================================================
				call electric(q(i,j-1,k-1,:),Ecen,1)
				a = 2.d0*(flxz(i,j-1,k,QMAGY) - Ecen)

				call electric(q(i,j-1,k,:),Ecen,1)
				b = 2.d0*(Ecen - flxz(i,j-1,k,QMAGY))

				if(flxz(i,j-1,k,QRHO).gt. 0.d0) then
					d1 = a
				elseif(flxz(i,j-1,k,QRHO).lt. 0.d0) then
					d1 = b
				else
					d1 = 0.5d0*(a+b)
				endif
				a = b
				b = 2.d0*(flxz(i,j-1,k+1,QMAGY) - Ecen)
				
				if(q(i,j-1,k,QW).gt. 0.d0) then
					d2 = a
				elseif(q(i,j-1,k,QW).lt. 0.d0) then
					d2 = b
				else
					d2 = 0.5d0*(a+b)
				endif
				dd1 = 0.125d0*(d1 - d2)
			
				call electric(q(i,j-1,k-1,:),Ecen,1)
				a = 2.d0*(-flxy(i,j,k-1,QMAGZ) - Ecen)

				call electric(q(i,j,k-1,:),Ecen,1)
				b = 2.d0*(Ecen + flxy(i,j,k-1,QMAGZ))

				if(flxy(i,j,k-1,QRHO).gt. 0.d0) then
					d1 = a
				elseif(flxy(i,j,k-1,QRHO).lt. 0.d0) then
					d1 = b
				else
					d1 = 0.5d0*(a+b)
				endif

				a = b 
				b = 2.d0*(-flxy(i,j+1,k-1,QMAGZ) - Ecen) 

				if(q(i,j,k-1,QV).gt.0.d0) then
					d2 = a
				elseif(q(i,j,k-1,QV).lt. 0.d0) then
					d2 = b
				else 
					d2 = 0.5d0*(a+b)
				endif
				dd2 = 0.125*(d1 - d2)

				E(i,j,k) = 0.25d0*(-flxy(i,j,k-1,QMAGZ) - flxy(i,j,k,QMAGZ) + &
                                                    flxz(i,j-1,k,QMAGY) + flxz(i,j,k,QMAGY))! + dd1 + dd2
!				E(i,j,k) = 0.5d0*(-flxy(i,j,k-1,QMAGZ) + flxz(i,j-1,k,QMAGY))

			
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

	real(rt), intent(in)	::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)

        real(rt), intent(in) :: flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,QVAR)
        real(rt), intent(in) :: flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,QVAR)

	real(rt), intent(out)	:: E(ey_l1:ey_h1,ey_l2:ey_h2,ey_l3:ey_h3)
	
	real(rt)				::Ecen
	real(rt)				::a ,b ,d1 ,d2 ,dd1 ,dd2 
	
	integer					::i ,j ,k	

	E = 0.d0
	do k = work_lo(3), work_hi(3)
		do j = work_lo(2), work_hi(2)
			do i = work_lo(1), work_hi(1)

!============================================= Ey i-1/2, j, k-1/2 ===========================================================================

				call electric(q(i-1,j,k-1,:),Ecen,2)
				a = 2.d0*(-flxx(i,j,k-1,QMAGZ) - Ecen)

				call electric(q(i,j,k-1,:),Ecen,2)
				b = 2.d0*(Ecen - flxx(i,j,k-1,QMAGZ))

				if(flxx(i,j,k-1,QRHO).gt. 0.d0) then
					d1 = a
				elseif(flxx(i,j,k-1,QRHO).lt. 0.d0) then
					d1 = b
				else
					d1 = 0.5d0*(a+b)
				endif

				a = b
				b = 2.d0*(flxx(i+1,j,k-1,QMAGZ) - Ecen)

				if(q(i,j,k-1,QU).gt. 0.d0) then
					d2 = a
				elseif(q(i,j,k-1,QU).lt. 0.d0) then
					d2 = b
				else
					d2 = 0.5d0*(a+b)
				endif
				dd1 = 0.125d0*(d1 - d2)

				call electric(q(i-1,j,k-1,:),Ecen,2)
				a = 2.d0*(flxz(i-1,j,k,QMAGX) - Ecen)

    			        call electric(q(i-1,j,k,:),Ecen,2)
				b = 2.d0*(Ecen - flxz(i-1,j,k,QMAGX))

				if(flxz(i-1,j,k,QRHO).gt. 0.d0) then
					d1 = a
				elseif(flxz(i-1,j,k,QRHO).lt. 0.d0) then
					d1 = b
				else
					d1 = 0.5d0*(a+b)
				endif
				a = b 
				b = 2.d0*(flxz(i-1,j,k+1,QMAGX) - Ecen) 
				if(q(i-1,j,k,QW).gt.0.d0) then
					d2 = a
				elseif(q(i-1,j,k,QW).lt. 0.d0) then
					d2 = b
				else 
					d2 = 0.5d0*(a+b)
				endif
				dd2 = 0.125*(d1 - d2)

				E(i,j,k) = 0.25d0*( flxx(i,j,k-1,QMAGZ) + flxx(i,j,k,QMAGZ) &
                                                   -flxz(i-1,j,k,QMAGX) - flxz(i,j,k,QMAGX))!  + dd1 + dd2
!				E(i,j,k) = 0.5d0*( flxx(i,j,k-1,QMAGZ) -flxz(i-1,j,k,QMAGX))
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

        real(rt), intent(in)    ::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)

        real(rt), intent(in) :: flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,QVAR)
        real(rt), intent(in) :: flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,QVAR)

        real(rt), intent(out)   :: E(ez_l1:ez_h1,ez_l2:ez_h2,ez_l3:ez_h3)
	
	real(rt)		:: Ecen, u_face, v_face
	real(rt)		:: a ,b ,d1 ,d2 ,dd1 ,dd2 

	integer			:: i ,j ,k	

	!E = 0.d0

! ====  Ez i-1/2, j-1/2, k ====

	do k = work_lo(3), work_hi(3)
		do j = work_lo(2), work_hi(2)
			do i = work_lo(1), work_hi(1)

				call electric(q(i-1,j-1,k,:),Ecen,3)
				a = 2.d0*(-flxx(i,j-1,k,QMAGY) - Ecen)

				call electric(q(i,j-1,k,:),Ecen,3)
				b = 2.d0*(Ecen + flxx(i,j-1,k,QMAGY))

				u_face = flxx(i,j-1,k,QRHO)
				if (u_face.gt. 0.d0) then
					d1 = a
				elseif (u_face.lt. 0.d0) then
					d1 = b
				else 
					d1 = 0.5d0*(a+b)
				endif

				a = b
				b = 2.d0*(-flxx(i+1,j-1,k,QMAGY) - Ecen)

				u_face = q(i,j-1,k,QU)
				if (u_face.gt. 0.d0) then
					d2 = a
				elseif (u_face.lt. 0.d0) then
					d2 = b
				else
					d2 = 0.5d0*(a+b)
				endif
				dd1 = 0.125d0*(d1 - d2)
!			if((i.eq.3.and.j.eq.17.and.k.eq.1).or.(i.eq.4.and.j.eq.17.and.k.eq.1)) then
!				print *, "d1 = ", d1, "d2 = ", d2
!				print *, "u = ", u_face
!			endif

				call electric(q(i-1,j-1,k,:),Ecen,3)
				a = 2.d0*(flxy(i-1,j,k,QMAGX) - Ecen)

				call electric(q(i-1,j,k,:),Ecen,3)
				b = 2.d0*(Ecen - flxy(i-1,j,k,QMAGX))

				v_face = flxy(i-1,j,k,QRHO)
				if (v_face.gt. 0.d0) then
					d1 = a
				elseif (v_face.lt. 0.d0) then
					d1 = b
				else
					d1 = 0.5d0*(a+b)
				endif

				a = b 
				b = 2.d0*(flxy(i-1,j+1,k,QMAGX) - Ecen) 

				if(q(i-1,j,k,QV).gt.0.d0) then
					d2 = a
				elseif(q(i-1,j,k,QV).lt. 0.d0) then
					d2 = b
				else 
					d2 = 0.5d0*(a+b)
				endif
				dd2 = 0.125*(d1 - d2)

				E(i,j,k) = 0.25d0*(-flxx(i,j-1,k,QMAGY) - flxx(i,j,k,QMAGY) &
                                                   +flxy(i-1,j,k,QMAGX) + flxy(i,j,k,QMAGX))! + dd1 + dd2
!				 E(i,j,k) = 0.5d0*(-flxx(i,j-1,k,QMAGY) + flxy(i-1,j,k,QMAGX))
				!if((i.eq.3.and.j.eq.17.and.k.eq.1).or.(i.eq.4.and.j.eq.17.and.k.eq.1)) then
				!	print *, "Electric Z at ", i,j,k
				!	print *, E(i,j,k)
				!	print *,"Flux x", - flxx(i,j-1,k, QMAGY) , flxx(i,j,k,QMAGY)
				!	print *,"Flux y", flxy(i-1,j,k,QMAGX) , flxy(i,j,k,QMAGX)
				!	print *, "Q at", i, j, k, q(i,j,k,:)
				!	print *, "dd1 = ", dd1, "dd2 = ", dd2
				!	pause
				!endif
			enddo
		enddo
	enddo



	
end subroutine electric_edge_z

end module electric_field
