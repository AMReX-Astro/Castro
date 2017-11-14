module mhd_plm_module

!Module that gives a piecewise linear interpolation for the primitive variables 
!They are projected onto the characteristic variables for tracing. 

   use amrex_fort_module, only : rt => amrex_real
   use meth_params_module
   implicit none

   private vanleer, lvecx, lvecy, lvecz, rvecx, rvecy, rvecz, evals

   public plm 

contains 

  !
  ! characteristics based on u
  !
  !===========================================================================
  ! This is called from within threaded loops in advance_mhd_tile so *no* OMP here ...
  !===========================================================================

   subroutine plm(lo, hi, s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,&	
		  bx, bxl1, bxl2, bxl3, bxh1, bxh2, bxh3, &
		  by, byl1, byl2, byl3, byh1, byh2, byh3, &
		  bz, bzl1, bzl2, bzl3, bzh1, bzh2, bzh3, &
                  Ip,Im, ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dy,dz,dt)


    implicit none
    integer , intent(in   ) ::  	s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, lo(3), hi(3)
    integer , intent(in   ) ::  	ilo1,ilo2,ilo3,ihi1,ihi2,ihi3
    integer , intent(in   ) ::  	bxl1, bxl2, bxl3, bxh1, bxh2, bxh3	
    integer , intent(in   ) ::  	byl1, byl2, byl3, byh1, byh2, byh3
    integer , intent(in   ) ::  	bzl1, bzl2, bzl3, bzh1, bzh2, bzh3
 
    real(rt), intent(in   ) ::      s(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QVAR) !Primitive Vars
    real(rt), intent(in   ) ::      bx(bxl1:bxh1, bxl2:bxh2, bxl3:bxh3)	!Face Centered Magnetic Fields
    real(rt), intent(in   ) ::      by(byl1:byh1, byl2:byh2, byl3:byh3)
    real(rt), intent(in   ) ::      bz(bzl1:bzh1, bzl2:bzh2, bzl3:bzh3)

    real(rt), intent(out) :: 		Ip(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,QVAR,3)
    real(rt), intent(out) :: 		Im(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,QVAR,3)

    real(rt), intent(in   ) :: 		dx,dy,dz,dt

    real(rt) 				:: 	dQL(7), dQR(7), dW, dL, dR, leig(7,7), reig(7,7), lam(7), summ(7)
    real(rt)				:: 	temp(s_l1-1:s_h1+1,s_l2-1:s_h2+1,s_l3-1:s_h3+1,8), smhd(7)
    real(rt)				::      tbx(s_l1-1:s_h1+1,s_l2-1:s_h2+1,s_l3-1:s_h3+1)
    real(rt)				::      tby(s_l1-1:s_h1+1,s_l2-1:s_h2+1,s_l3-1:s_h3+1)
    real(rt)				::      tbz(s_l1-1:s_h1+1,s_l2-1:s_h2+1,s_l3-1:s_h3+1)
    real(rt)				:: 	dt_over_a
    integer          		        ::	ii,ibx,iby,ibz, i , j, k

    ibx = 6
    iby = 7
    ibz = 8

    dt_over_a = dt
!    Ip = 0.d0
!    Im = 0.d0

!------------------------workspace variables---------------------------------------------
		temp = 0.d0
		temp(:,:,:,1) = small_dens
		temp(s_l1: s_h1, s_l2: s_h2, s_l3: s_h3,1  :  5) = s(:,:,:,QRHO:QPRES) !Gas vars Cell Centered
		temp(s_l1: s_h1, s_l2: s_h2, s_l3: s_h3,ibx:ibz) = s(:,:,:,QMAGX:QMAGZ) !Mag vars Cell Centered
!-------------------- Fill Boundaries ---------------------------------------------------
		temp(s_l1-1,s_l2-1,s_l3-1,1:5) = s(s_l1,s_l2,s_l3,QRHO:QPRES)
		temp(s_l1-1,s_l2-1,s_l3-1,6:8) = s(s_l1,s_l2,s_l3,QMAGX:QMAGZ)
		temp(s_l1-1, s_l2: s_h2, s_l3: s_h3,1:5) = s(s_l1,:,:,QRHO:QPRES)
		temp(s_l1-1, s_l2: s_h2, s_l3: s_h3,6:8) = s(s_l1,:,:,QMAGX:QMAGZ)
		temp(s_l1:s_h1, s_l2-1, s_l3: s_h3,1:5) = s(:,s_l2,:,QRHO:QPRES)
		temp(s_l1:s_h1, s_l2-1, s_l3: s_h3,6:8) = s(:,s_l2,:,QMAGX:QMAGZ)
		temp(s_l1:s_h1, s_l2:s_h2, s_l3-1,1:5) = s(:,:,s_l3,QRHO:QPRES)
		temp(s_l1:s_h1, s_l2:s_h2, s_l3-1,6:8) = s(:,:,s_l3,QMAGX:QMAGZ)
		temp(s_h1+1, s_l2: s_h2, s_l3: s_h3,1:5) = s(s_h1,:,:,QRHO:QPRES)
		temp(s_h1+1, s_l2: s_h2, s_l3: s_h3,6:8) = s(s_h1,:,:,QMAGX:QMAGZ)
		temp(s_l1:s_h1, s_h2+1, s_l3: s_h3,1:5) = s(:,s_h2,:,QRHO:QPRES)
		temp(s_l1:s_h1, s_h2+1, s_l3: s_h3,6:8) = s(:,s_h2,:,QMAGX:QMAGZ)
		temp(s_l1:s_h1, s_l2:s_h2, s_h3+1,1:5) = s(:,:,s_h3,QRHO:QPRES)
		temp(s_l1:s_h1, s_l2:s_h2, s_h3+1,6:8) = s(:,:,s_h3,QMAGX:QMAGZ)
		temp(s_h1+1,s_h2+1,s_h3+1,1:5) = s(s_h1,s_h2,s_h3,QRHO:QPRES)
		temp(s_h1+1,s_h2+1,s_h3+1,6:8) = s(s_h1,s_h2,s_h3,QMAGX:QMAGZ)
! Temp face centered magnetic fields
	do k = s_l3-1, s_h3+1
		do j = s_l2-1, s_h2+1
			do i = s_l1-1, s_h1+1
	!---------------------------------------------- set up temporary bx ------------------------------------------------
				if(i.lt.bxl1) then 
					if(j.lt.bxl2) then
						if(k.lt.bxl3) then 
							tbx(i,j,k) = bx(bxl1,bxl2,bxl3)
						elseif(k.lt.bxh3) then
							tbx(i,j,k) = bx(bxl1, bxl2, k)
						else
							tbx(i,j,k) = bx(bxl1, bxl2, bxh3)
						endif
					elseif(j.lt.bxh2) then
						if(k.lt.bxl3) then 
							tbx(i,j,k) = bx(bxl1,j,bxl3)
						elseif(k.lt.bxh3) then
							tbx(i,j,k) = bx(bxl1, j, k)
						else
							tbx(i,j,k) = bx(bxl1, j, bxh3)
						endif
					else
						if(k.lt.bxl3) then 
							tbx(i,j,k) = bx(bxl1,bxh2,bxl3)
						elseif(k.lt.bxh3) then
							tbx(i,j,k) = bx(bxl1, bxh2, k)
						else
							tbx(i,j,k) = bx(bxl1, bxh2, bxh3)
						endif
					endif
				elseif(i.lt.bxh1) then 
					if(j.lt.bxl2) then
						if(k.lt.bxl3) then 
							tbx(i,j,k) = bx(i,bxl2,bxl3)
						elseif(k.lt.bxh3) then
							tbx(i,j,k) = bx(i, bxl2, k)
						else
							tbx(i,j,k) = bx(i, bxl2, bxh3)
						endif
					elseif(j.lt.bxh2) then
						if(k.lt.bxl3) then 
							tbx(i,j,k) = bx(i,j,bxl3)
						elseif(k.lt.bxh3) then
							tbx(i,j,k) = bx(i, j, k)
						else
							tbx(i,j,k) = bx(i, j, bxh3)
						endif
					else
						if(k.lt.bxl3) then 
							tbx(i,j,k) = bx(i,bxh2,bxl3)
						elseif(k.lt.bxh3) then
							tbx(i,j,k) = bx(i, bxh2, k)
						else
							tbx(i,j,k) = bx(i, bxh2, bxh3)
						endif
					endif
				else
					if(j.lt.bxl2) then
						if(k.lt.bxl3) then 
							tbx(i,j,k) = bx(bxh1,bxl2,bxl3)
						elseif(k.lt.bxh3) then
							tbx(i,j,k) = bx(bxh1, bxl2, k)
						else
							tbx(i,j,k) = bx(bxh1, bxl2, bxh3)
						endif
					elseif(j.lt.bxh2) then
						if(k.lt.bxl3) then 
							tbx(i,j,k) = bx(bxh1,j,bxl3)
						elseif(k.lt.bxh3) then
							tbx(i,j,k) = bx(bxh1, j, k)
						else
							tbx(i,j,k) = bx(bxh1, j, bxh3)
						endif
					else
						if(k.lt.bxl3) then 
							tbx(i,j,k) = bx(bxh1,bxh2,bxl3)
						elseif(k.lt.bxh3) then
							tbx(i,j,k) = bx(bxh1, bxh2, k)
						else
							tbx(i,j,k) = bx(bxh1, bxh2, bxh3)
						endif
					endif
				endif

	!---------------------------------------------- set up temporary by ------------------------------------------------
				if(i.lt.byl1) then 
					if(j.lt.byl2) then
						if(k.lt.byl3) then 
							tby(i,j,k) = by(byl1,byl2,byl3)
						elseif(k.lt.byh3) then
							tby(i,j,k) = by(byl1, byl2, k)
						else
							tby(i,j,k) = by(byl1, byl2, byh3)
						endif
					elseif(j.lt.byh2) then
						if(k.lt.byl3) then 
							tby(i,j,k) = by(byl1,j,byl3)
						elseif(k.lt.byh3) then
							tby(i,j,k) = by(byl1, j, k)
						else
							tby(i,j,k) = by(byl1, j, byh3)
						endif
					else
						if(k.lt.byl3) then 
							tby(i,j,k) = by(byl1,byh2,byl3)
						elseif(k.lt.byh3) then
							tby(i,j,k) = by(byl1, byh2, k)
						else
							tby(i,j,k) = by(byl1, byh2, byh3)
						endif
					endif
				elseif(i.lt.byh1) then 
					if(j.lt.byl2) then
						if(k.lt.byl3) then 
							tby(i,j,k) = by(i,byl2,byl3)
						elseif(k.lt.byh3) then
							tby(i,j,k) = by(i, byl2, k)
						else
							tby(i,j,k) = by(i, byl2, byh3)
						endif
					elseif(j.lt.byh2) then
						if(k.lt.byl3) then 
							tby(i,j,k) = by(i,j,byl3)
						elseif(k.lt.byh3) then
							tby(i,j,k) = by(i, j, k)
						else
							tby(i,j,k) = by(i, j, byh3)
						endif
					else
						if(k.lt.byl3) then 
							tby(i,j,k) = by(i,byh2,byl3)
						elseif(k.lt.byh3) then
							tby(i,j,k) = by(i, byh2, k)
						else
							tby(i,j,k) = by(i, byh2, byh3)
						endif
					endif
				else
					if(j.lt.byl2) then
						if(k.lt.byl3) then 
							tby(i,j,k) = by(byh1,byl2,byl3)
						elseif(k.lt.byh3) then
							tby(i,j,k) = by(byh1, byl2, k)
						else
							tby(i,j,k) = by(byh1, byl2, byh3)
						endif
					elseif(j.lt.byh2) then
						if(k.lt.byl3) then 
							tby(i,j,k) = by(byh1,j,byl3)
						elseif(k.lt.byh3) then
							tby(i,j,k) = by(byh1, j, k)
						else
							tby(i,j,k) = by(byh1, j, byh3)
						endif
					else
						if(k.lt.byl3) then 
							tby(i,j,k) = by(byh1,byh2,byl3)
						elseif(k.lt.byh3) then
							tby(i,j,k) = by(byh1, byh2, k)
						else
							tby(i,j,k) = by(byh1, byh2, byh3)
						endif
					endif
				endif

	!---------------------------------------------- set up temporary bz ------------------------------------------------
				if(i.lt.bzl1) then 
					if(j.lt.bzl2) then
						if(k.lt.bzl3) then 
							tbz(i,j,k) = bz(bzl1,bzl2,bzl3)
						elseif(k.lt.bzh3) then
							tbz(i,j,k) = bz(bzl1, bzl2, k)
						else
							tbz(i,j,k) = bz(bzl1, bzl2, bzh3)
						endif
					elseif(j.lt.bzh2) then
						if(k.lt.bzl3) then 
							tbz(i,j,k) = bz(bzl1,j,bzl3)
						elseif(k.lt.bzh3) then
							tbz(i,j,k) = bz(bzl1, j, k)
						else
							tbz(i,j,k) = bz(bzl1, j, bzh3)
						endif
					else
						if(k.lt.bzl3) then 
							tbz(i,j,k) = bz(bzl1,bzh2,bzl3)
						elseif(k.lt.bzh3) then
							tbz(i,j,k) = bz(bzl1, bzh2, k)
						else
							tbz(i,j,k) = bz(bzl1, bzh2, bzh3)
						endif
					endif
				elseif(i.lt.bzh1) then 
					if(j.lt.bzl2) then
						if(k.lt.bzl3) then 
							tbz(i,j,k) = bz(i,bzl2,bzl3)
						elseif(k.lt.bzh3) then
							tbz(i,j,k) = bz(i, bzl2, k)
						else
							tbz(i,j,k) = bz(i, bzl2, bzh3)
						endif
					elseif(j.lt.bzh2) then
						if(k.lt.bzl3) then 
							tbz(i,j,k) = bz(i,j,bzl3)
						elseif(k.lt.bzh3) then
							tbz(i,j,k) = bz(i, j, k)
						else
							tbz(i,j,k) = bz(i, j, bzh3)
						endif
					else
						if(k.lt.bzl3) then 
							tbz(i,j,k) = bz(i,bzh2,bzl3)
						elseif(k.lt.bzh3) then
							tbz(i,j,k) = bz(i, bzh2, k)
						else
							tbz(i,j,k) = bz(i, bzh2, bzh3)
						endif
					endif
				else
					if(j.lt.bzl2) then
						if(k.lt.bzl3) then 
							tbz(i,j,k) = bz(bzh1,bzl2,bzl3)
						elseif(k.lt.bzh3) then
							tbz(i,j,k) = bz(bzh1, bzl2, k)
						else
							tbz(i,j,k) = bz(bzh1, bzl2, bzh3)
						endif
					elseif(j.lt.bzh2) then
						if(k.lt.bzl3) then 
							tbz(i,j,k) = bz(bzh1,j,bzl3)
						elseif(k.lt.bzh3) then
							tbz(i,j,k) = bz(bzh1, j, k)
						else
							tbz(i,j,k) = bz(bzh1, j, bzh3)
						endif
					else
						if(k.lt.bzl3) then 
							tbz(i,j,k) = bz(bzh1,bzh2,bzl3)
						elseif(k.lt.bzh3) then
							tbz(i,j,k) = bz(bzh1, bzh2, k)
						else
							tbz(i,j,k) = bz(bzh1, bzh2, bzh3)
						endif
					endif
				endif
			enddo
		enddo
	enddo

	!============================================== PLM ================================================================
	do k = s_l3, s_h3
		do j = s_l2, s_h2
		do i = s_l1, s_h1
	
	!============================================ X Direction ==============================================
			summ = 0.d0
			smhd = 0.d0
			dQL = 0.d0
			dQR = 0.d0
			dW = 0.d0
			reig = 0.d0
			leig = 0.d0
			lam = 0.d0
			!Skip Bx
			dQL(1:5) = 	temp(i,j,k,1:ibx-1) - temp(i-1,j,k,1:ibx-1) !gas
			dQL(6:7) = 	temp(i,j,k,ibx+1:8) - temp(i-1,j,k,ibx+1:8)	!mag			
			dQR(1:5) = 	temp(i+1,j,k,1:ibx-1) - temp(i,j,k,1:ibx-1)
			dQR(6:7) = 	temp(i+1,j,k,ibx+1:8) - temp(i,j,k,ibx+1:8)				
!			do ii = 1,7
!				!call vanleer(dW(ii),dQL(ii),dQR(ii)) !!slope limiting
!				call minmod(dW(ii),dQL(ii),dQR(ii)) !!slope limiting
!			enddo
			call evals(lam, s(i,j,k,:), 1) !!X dir eigenvalues
			call lvecx(leig,s(i,j,k,:))    !!left eigenvectors
			call rvecx(reig,s(i,j,k,:))    !!right eigenvectors
	!MHD Source Terms 
			smhd(2) = temp(i,j,k,ibx)/temp(i,j,k,1)
			smhd(3) = temp(i,j,k,iby)/temp(i,j,k,1)
			smhd(4) = temp(i,j,k,ibz)/temp(i,j,k,1)
			smhd(5) = temp(i,j,k,ibx)*temp(i,j,k,2) + temp(i,j,k,iby)*temp(i,j,k,3) + temp(i,j,k,ibz)*temp(i,j,k,4)
			smhd(6) = temp(i,j,k,3)
			smhd(7) = temp(i,j,k,4)
			smhd 	= smhd*(tbx(i+1,j,k) - tbx(i,j,k))/dx !cross-talk of normal magnetic field direction
	!Interpolate
		!Plus
				!!Using HLLD so sum over all eigenvalues
			do ii = 1,7
			        dL = dot_product(leig(ii,:),dQL)
			        dR = dot_product(leig(ii,:),dQR)
			        call vanleer(dW,dL,dR)
					summ(:) = summ(:) + (1 - dt_over_a/dx)*lam(ii)*dW*reig(:,ii)
			enddo
			Ip(i,j,k,QRHO:QPRES,1) 	 = temp(i,j,k,1:ibx-1) + 0.5d0*summ(1:5) + 0.5d0*dt_over_a*smhd(1:5)
			Ip(i,j,k,QMAGX,1) 		 = temp(i+1,j,k,ibx) !! Bx stuff
			Ip(i,j,k,QMAGY:QMAGZ,1)  = temp(i,j,k,iby:ibz) + 0.5d0*summ(6:7) + 0.5d0*dt_over_a*smhd(6:7)
			Ip(i,j,k,QPRES,1)        = Ip(i,j,k,QPRES,1) + 0.5d0*dot_product(Ip(i,j,k,QMAGX:QMAGZ,1),Ip(i,j,k,QMAGX:QMAGZ,1))
		!Minus
			summ = 0.d0
			do ii = 1,7
			        dL = dot_product(leig(ii,:),dQL)
			        dR = dot_product(leig(ii,:),dQR)
			        call vanleer(dW,dL,dR)			
					summ(:) = summ(:) + (- 1 - dt_over_a/dx*lam(ii))*dW*reig(:,ii)
			enddo
			Im(i,j,k,QRHO:QPRES,1)	 = temp(i,j,k,1:ibx-1) +0.5d0*summ(1:5) + 0.5d0*dt_over_a*smhd(1:5)
			Im(i,j,k,QMAGX,1)		 = temp(i-1,j,k,ibx) !! Bx stuff
			Im(i,j,k,QMAGY:QMAGZ,1)  = temp(i,j,k,iby:ibz) +0.5d0*summ(6:7) + 0.5d0*dt_over_a*smhd(6:7)
			Im(i,j,k,QPRES,1)        = Im(i,j,k,QPRES,1) + 0.5d0*dot_product(Im(i,j,k,QMAGX:QMAGZ,1),Im(i,j,k,QMAGX:QMAGZ,1))


	
	!========================================= Y Direction ================================================				

			summ = 0.d0
			smhd = 0.d0
			dQL = 0.d0
			dQR = 0.d0
			dW = 0.d0
			reig = 0.d0
			leig = 0.d0
			lam = 0.d0
			!Skip By
			dQL(1:6) = 	temp(i,j,k,1:ibx) - temp(i,j-1,k,1:ibx) !gas + bx
			dQL(7) = 	temp(i,j,k,8) - temp(i,j-1,k,iby+1)		!bz			
			dQR(1:6) = 	temp(i,j+1,k,1:ibx) - temp(i,j,k,1:ibx)
			dQR(7) = 	temp(i,j+1,k,ibz) - temp(i,j,k,ibz)				
			!do ii = 1,7
			!	!call vanleer(dW(ii),dQL(ii),dQR(ii)) !!slope limiting
			!	call minmod(dW(ii),dQL(ii),dQR(ii)) !!slope limiting						
			!enddo
			call evals(lam, s(i,j,k,:), 2) !!Y dir eigenvalues
			call lvecy(leig,s(i,j,k,:))    !!left eigenvectors
			call rvecy(reig,s(i,j,k,:))    !!right eigenvectors
			!!Using HLLD so sum over all eigenvalues
			do ii = 1,7
			        dL = dot_product(leig(ii,:),dQL)
			        dR = dot_product(leig(ii,:),dQR)
			        call vanleer(dW,dL,dR)			
					summ(:) = summ(:) + (1 - dt_over_a/dx)*lam(ii)*dW*reig(:,ii)
			enddo
	!MHD Source Terms 
			smhd(2) = temp(i,j,k,ibx)/temp(i,j,k,1)
			smhd(3) = temp(i,j,k,iby)/temp(i,j,k,1)
			smhd(4) = temp(i,j,k,ibz)/temp(i,j,k,1)
			smhd(5) = temp(i,j,k,ibx)*temp(i,j,k,2) + temp(i,j,k,iby)*temp(i,j,k,3) + temp(i,j,k,ibz)*temp(i,j,k,4)
			smhd(6) = temp(i,j,k,2)
			smhd(7) = temp(i,j,k,4)
			smhd 	= smhd*(tby(i,j+1,k) - tby(i,j,k))/dy !cross-talk of normal magnetic field direction
	!Interpolate
			Ip(i,j,k,QRHO:QPRES,2) 	= temp(i,j,k,1:ibx-1) +0.5d0*summ(1:5) + 0.5d0*dt_over_a*smhd(1:5) !!GAS
			Ip(i,j,k,QMAGX,2) 		= temp(i,j,k,ibx) + 0.5d0*summ(6) + 0.5d0*dt_over_a*smhd(6)
			Ip(i,j,k,QMAGY,2) 		= temp(i,j+1,k,iby) !! By stuff
			Ip(i,j,k,QMAGZ,2)  		= temp(i,j,k,ibz) + 0.5d0*summ(7) + 0.5d0*dt_over_a*smhd(7)
			Ip(i,j,k,QPRES,2)       = Ip(i,j,k,QPRES,2) + 0.5d0*dot_product(Ip(i,j,k,QMAGX:QMAGZ,2),Ip(i,j,k,QMAGX:QMAGZ,2))

			summ = 0.d0
			do ii = 1,7
			        dL = dot_product(leig(ii,:),dQL)
			        dR = dot_product(leig(ii,:),dQR)
			        call vanleer(dW,dL,dR)			
					summ(:) = summ(:) + (- 1 - dt_over_a/dx*lam(ii))*dW*reig(:,ii)
			enddo
			Im(i,j,k,QRHO:QPRES,2)	= temp(i,j,k,1:ibx-1) + 0.5d0*summ(1:5) + 0.5d0*dt_over_a*smhd(1:5) !!GAS
			Im(i,j,k,QMAGX,2) 		= temp(i,j,k,ibx) + 0.5d0*summ(6) + 0.5d0*dt_over_a*smhd(6)
			Im(i,j,k,QMAGY,2)		= temp(i,j-1,k,iby) !! By stuff
			Im(i,j,k,QMAGZ,2) 		= temp(i,j,k,ibz) + 0.5d0*summ(7) + 0.5d0*dt_over_a*smhd(7)
			Im(i,j,k,QPRES,2)       = Im(i,j,k,QPRES,2) + 0.5d0*dot_product(Im(i,j,k,QMAGX:QMAGZ,2),Im(i,j,k,QMAGX:QMAGZ,2))
			
	!========================================= Z Direction ================================================				
			summ = 0.d0
			smhd = 0.d0
			dQL = 0.d0
			dQR = 0.d0
			dW = 0.d0
			reig = 0.d0
			leig = 0.d0
			lam = 0.d0
			!Skip Bz
			dQL(1:7) = 	temp(i,j,k,1:iby) - temp(i,j,k-1,1:iby) 
			dQR(1:7) = 	temp(i,j,k+1,1:iby) - temp(i,j,k,1:iby)
!			do ii = 1,7
				!call vanleer(dW(ii),dQL(ii),dQR(ii)) !!slope limiting
!				call minmod(dW(ii), dQL(ii), dQR(ii))
!			enddo
			call evals(lam, s(i,j,k,:), 3) !!Z dir eigenvalues
			call lvecz(leig,s(i,j,k,:))    !!left eigenvectors
			call rvecz(reig,s(i,j,k,:))    !!right eigenvectors

			!!Characteristic Tracing
			do ii = 1,7
			        dL = dot_product(leig(ii,:),dQL)
			        dR = dot_product(leig(ii,:),dQR)
			        call vanleer(dW,dL,dR)			
				summ(:) = summ(:) + (1 - dt_over_a/dx)*lam(ii)*dW*reig(:,ii)
			enddo
	!MHD Source Terms 
			smhd(2) = temp(i,j,k,ibx)/temp(i,j,k,1)
			smhd(3) = temp(i,j,k,iby)/temp(i,j,k,1)
			smhd(4) = temp(i,j,k,ibz)/temp(i,j,k,1)
			smhd(5) = temp(i,j,k,ibx)*temp(i,j,k,2) + temp(i,j,k,iby)*temp(i,j,k,3) + temp(i,j,k,ibz)*temp(i,j,k,4)
			smhd(6) = temp(i,j,k,2)
			smhd(7) = temp(i,j,k,3)
			smhd 	= smhd*(tbz(i,j,k+1) - tbz(i,j,k))/dz !cross-talk of normal magnetic field direction

	!Interpolate
			Ip(i,j,k,QRHO:QPRES,3) 	= temp(i,j,k,1:ibx-1) + 0.5d0*summ(1:5) + 0.5d0*dt_over_a*smhd(1:5) !!GAS
			Ip(i,j,k,QMAGX:QMAGY,3)	= temp(i,j,k,ibx:iby) + 0.5d0*summ(6:7) + 0.5d0*dt_over_a*smhd(6:7)
			Ip(i,j,k,QMAGZ,3) 		= temp(i,j,k+1,ibz) !! Bz stuff
			Ip(i,j,k,QPRES,3)       = Ip(i,j,k,QPRES,3) + 0.5d0*dot_product(Ip(i,j,k,QMAGX:QMAGZ,3),Ip(i,j,k,QMAGX:QMAGZ,3))
			summ = 0.d0
			do ii = 1,7
			        dL = dot_product(leig(ii,:),dQL)
			        dR = dot_product(leig(ii,:),dQR)
			        call vanleer(dW,dL,dR)			
				summ(:) = summ(:) + (- 1 - dt_over_a/dx*lam(ii))*dW*reig(:,ii)
			enddo
			Im(i,j,k,QRHO:QPRES,3)	= temp(i,j,k,1:ibx-1) + 0.5d0*summ(1:5) + 0.5d0*dt_over_a*smhd(1:5) !!GAS
			Im(i,j,k,QMAGX:QMAGY,3) = temp(i,j,k,ibx:iby) + 0.5d0*summ(6:7) + 0.5d0*dt_over_a*smhd(6:7)
			Im(i,j,k,QMAGZ,3)		= temp(i,j,k-1,ibz) !! Bz stuff
			Im(i,j,k,QPRES,3)       = Im(i,j,k,QPRES,3) + 0.5d0*dot_product(Im(i,j,k,QMAGX:QMAGZ,3),Im(i,j,k,QMAGX:QMAGZ,3))

		enddo
		enddo
	enddo

	Im(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QREINT,1)       = s(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QREINT)
	Im(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QREINT,2)       = s(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QREINT)
	Im(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QREINT,3)       = s(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QREINT)
	Ip(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QREINT,1)       = s(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QREINT)
	Ip(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QREINT,2)       = s(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QREINT)
	Ip(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QREINT,3)       = s(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,QREINT)

!Need to add source terms, heating cooling, gravity, etc.
	end subroutine plm
!======================================== Minmod TVD slope limiter =========================================
	subroutine minmod(dW, WR, WL)
	
	use amrex_fort_module, only : rt => amrex_real

	implicit none

	real(rt), intent(in) :: WR, WL
	real(rt), intent(out) :: dW
	dW = 0.d0

	if(abs(WR).lt.abs(WL).and.WR*WL.gt. 0.d0) then
		dW = WR
	elseif(abs(WL).lt.abs(WR).and.WR*WL.gt.0.d0) then
		dW = WL
	endif

	end subroutine


!========================================= VanLeer TVD slope limiter =======================================
	subroutine vanleer(dW, WR, WL) 
	
	use amrex_fort_module, only : rt => amrex_real
	
	implicit none
	
	real(rt), intent(in )	::  WR, WL
	real(rt), intent(out)	::  dW
	dW = 0.0d0 	
	
	if( WR*WL .gt. 0.0d0 ) then 
	dW = 2.0d0*WR*WL/(WR + WL)
	endif 

	end subroutine


!=========================================== Evals =========================================================

	subroutine evals(lam, Q, dir)

	use amrex_fort_module, only : rt => amrex_real
        use eos_module, only: eos
        use eos_type_module, only: eos_t, eos_input_rp
        use actual_network, only: nspec
	implicit none
	
	real(rt), intent(in)	:: Q(QVAR)
	real(rt), intent(out)	:: lam(7) !7 waves
	integer, intent(in)		:: dir !Choose direction, 1 for x, 2 for y, 3 for z

	!The characteristic speeds of the system 
	real(rt)				:: cfx, cfy, cfz, cax, cay, caz, csx, csy, csz, ca, as
        type(eos_t) :: eos_state

	!Speeeeeeeedssssss
        ! TODO: gamma_const -> gamma_1 from EOS
        eos_state % rho = Q(QRHO)
        eos_state % p   = Q(QPRES)
        eos_state % xn  = Q(QFS:QFS+nspec-1)

        call eos(eos_input_rp, eos_state)   

	as = eos_state % gam1 * (Q(QPRES) - 0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/Q(QRHO)
	!Alfven
	ca  = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	cax = (Q(QMAGX)**2)/Q(QRHO)
	cay = (Q(QMAGY)**2)/Q(QRHO)
	caz = (Q(QMAGZ)**2)/Q(QRHO)
	!Sloooooooooow
	csx = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))
	csy = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))
	csz = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))
	!Fassssst
	cfx = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))
	cfy = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))
	cfz = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))
	if(dir.eq.1) then	
		!Ax eigenvalues
		lam(1) = Q(QU) - sqrt(cfx)
		lam(2) = Q(QU) - sqrt(cax)
		lam(3) = Q(QU) - sqrt(csx)
		lam(4) = Q(QU)
		lam(5) = Q(QU) + sqrt(csx)
		lam(6) = Q(QU) + sqrt(cax)
		lam(7) = Q(QU) + sqrt(cfx)
	elseif(dir.eq.2) then
		!Ay eigenvalues
		lam(1) = Q(QV) - sqrt(cfy)
		lam(2) = Q(QV) - sqrt(cay)
		lam(3) = Q(QV) - sqrt(csy)
		lam(4) = Q(QV)
		lam(5) = Q(QV) + sqrt(csy)
		lam(6) = Q(QV) + sqrt(cay)
		lam(7) = Q(QV) + sqrt(cfy)
	else
		!Az eigenvalues
		lam(1) = Q(QW) - sqrt(cfz)
		lam(2) = Q(QW) - sqrt(caz)
		lam(3) = Q(QW) - sqrt(csz)
		lam(4) = Q(QW)
		lam(5) = Q(QW) + sqrt(csz)
		lam(6) = Q(QW) + sqrt(caz)
		lam(7) = Q(QW) + sqrt(cfz)
	endif
	end subroutine evals
	
!====================================== Left Eigenvectors ===============================================

!x direction
	subroutine lvecx(leig, Q) 
	use amrex_fort_module, only : rt => amrex_real
        use eos_module, only : eos
        use eos_type_module, only: eos_t, eos_input_rp
        use network, only : nspec
	
	implicit none
	
	!returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Ax
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::leig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfx, cax, csx, ca, as, S, N
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, bety, betz
        type (eos_t) :: eos_state

	!Speeeeeeeedssssss
        ! TODO: gamma_const -> gamma_1 from EOS
        eos_state % rho = Q(QRHO)
        eos_state % p   = Q(QPRES)
        eos_state % xn  = Q(QFS:QFS+nspec-1)

        call eos(eos_input_rp, eos_state)   

	as = eos_state % gam1 * (Q(QPRES) - 0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/Q(QRHO)
	!Alfven
	ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	cax = (Q(QMAGX)**2)/Q(QRHO)
	!Sloooooooooow
	csx = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))
	!Fassssst
	cfx = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))
	!useful constants
	alf = sqrt((as - csx)/(cfx - csx))
	als = sqrt((cfx - as)/(cfx - csx))
	if(cfx - as .lt. 0.d0) als = 0.d0
	if(as - csx .lt. 0.d0) alf = 0.d0
	if(abs(Q(QMAGY)).le. 1.d-14 .and.abs(Q(QMAGZ)).le. 1.d-14) then
		bety = 1.d0/sqrt(2.d0)
		betz = bety
	else
		bety = Q(QMAGY)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
		betz = Q(QMAGZ)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
	endif
	cff = sqrt(cfx)*alf
	css = sqrt(csx)*als
	S = sign(1.0d0, Q(QMAGX))
	Qf = sqrt(cfx)*alf*S
	Qs = sqrt(csx)*als*S
	N = 0.5d0/as
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))

		
	leig(1,:) = (/0.d0, -N*Cff	, N*Qs*bety		, N*Qs*betz		, N*alf/Q(QRHO)	, N*AAs*bety/Q(QRHO)			, N*AAs*betz/Q(QRHO)			/) !u - cf
	leig(2,:) = (/0.d0,  0.d0	, -0.5d0*betz	, 0.5d0*bety	, 0.d0			, -0.5d0*S*betz/(sqrt(Q(QRHO)))	, 0.5d0*bety*S/(sqrt(Q(QRHO)))	/) !u - cAx
	leig(3,:) = (/0.d0, -N*Css	, -N*Qf*bety	, -N*Qf*betz	, N*als/Q(QRHO)	, -N*AAf*bety/Q(QRHO)			, -N*AAf*betz/Q(QRHO)			/) !u - cs
	leig(4,:) = (/1.d0,  0.d0	,  0.d0			, 0.d0			, -1.d0/as		, 0.d0							, 0.d0							/) !u 
	leig(5,:) = (/0.d0,  N*Css	, N*Qf*bety		, N*Qf*betz		, N*als/Q(QRHO)	, -N*AAf*bety/Q(QRHO)			, -N*AAf*betz/Q(QRHO)			/) !u + cs
	leig(6,:) = (/0.d0,  0.d0	, 0.5d0*betz	, -0.5d0*bety	, 0.d0			, -0.5d0*betz*S/(sqrt(Q(QRHO)))	, 0.5d0*bety*S/(sqrt(Q(QRHO)))	/) !u + cAx
	leig(7,:) = (/0.d0, N*Cff	, -N*Qs*bety	, -N*Qs*betz	, N*alf/Q(QRHO)	, N*AAs*bety/Q(QRHO)			, N*AAs*betz/Q(QRHO)			/) !u + cf
	
	
	end subroutine lvecx

!y direction
	subroutine lvecy(leig, Q) 
	use amrex_fort_module, only : rt => amrex_real
	use eos_module, only : eos
        use eos_type_module, only: eos_t, eos_input_rp
        use network, only : nspec

	implicit none
	
	!returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Ay
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::leig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfy, cay, csy, ca, as, S, N
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, betz
        
        type (eos_t) :: eos_state

	!Speeeeeeeedssssss
        ! TODO: gamma_const -> gamma_1
        eos_state % rho = Q(QRHO)
        eos_state % p   = Q(QPRES)
        eos_state % xn  = Q(QFS:QFS+nspec-1)

        call eos(eos_input_rp, eos_state)   

	as = eos_state % gam1 * (Q(QPRES) - 0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/Q(QRHO)
	!Alfven
	ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	cay = (Q(QMAGY)**2)/Q(QRHO)
	!Sloooooooooow
	csy = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))
	!Fassssst
	cfy = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))
	!useful constants
	alf = sqrt((as - csy)/(cfy - csy))
	if(as - csy .lt. 0.d0) alf = 0.d0
	als = sqrt((cfy - as)/(cfy - csy))
	if(cfy - as .lt. 0.d0) als = 0.d0
	if(abs(Q(QMAGX)).le. 1.d-14 .and.abs(Q(QMAGZ)).le. 1.d-14) then
		betx = 1.d0/sqrt(2.d0)
		betz = betx
	else
		betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
		betz = Q(QMAGZ)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
	endif
	cff = sqrt(cfy)*alf
	css = sqrt(csy)*als
	S = sign(1.0d0, Q(QMAGY))
	Qf = sqrt(cfy)*alf*S
	Qs = sqrt(csy)*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
	N = 0.5d0/as
	
!Need to double check the rows
	leig(1,:) = (/0.d0, -N*Cff , N*Qs*betz		, N*Qs*betx 	, N*alf/Q(QRHO)	, N*AAs*betz/Q(QRHO)			, N*AAs*betx/Q(QRHO)			/) ! v - cf
	leig(2,:) = (/0.d0,  0.d0  , -0.5d0*betx	, 0.5d0*betz	, 0.d0			, -0.5d0*betx*S/(sqrt(Q(QRHO)))	, 0.5d0*betz*S/(sqrt(Q(QRHO)))	/) ! v - cAy
	leig(3,:) = (/0.d0, -N*Css , -N*Qf*betz 	, -N*Qf*betx	, N*als/Q(QRHO)	, -N*AAf*betz/Q(QRHO)			, -N*AAf*betx/Q(QRHO)			/) ! v - cs
	leig(4,:) = (/1.d0,  0.d0  ,  0.d0			, 0.d0			, -1.d0/as		, 0.d0							, 0.d0							/) ! v 
	leig(5,:) = (/0.d0,  N*Css , N*Qf*betz		, N*Qf*betx		, N*als/Q(QRHO)	, -N*AAf*betz/Q(QRHO)			, -N*AAf*betx/Q(QRHO)			/) ! v + cs
	leig(6,:) = (/0.d0,  0.d0  , 0.5d0*betx		, -0.5d0*betz	, 0.d0			, -0.5d0*betx*S/(sqrt(Q(QRHO)))	, 0.5d0*betz*S/(sqrt(Q(QRHO)))	/) ! v + cAy
	leig(7,:) = (/0.d0, N*Cff  , -N*Qs*betz		, -N*Qs*betx	, N*alf/Q(QRHO)	, N*AAs*betz/Q(QRHO)			, N*AAs*betx/Q(QRHO)			/) ! v + cf
	
	
	end subroutine lvecy

!z direction
	subroutine lvecz(leig, Q) 
	use amrex_fort_module, only : rt => amrex_real
	use eos_module, only : eos
        use eos_type_module, only: eos_t, eos_input_rp
        use network, only : nspec

	implicit none
	
	!returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Az
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::leig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfz, caz, csz, ca, as, S, N
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, bety
        
        type (eos_t) :: eos_state

	!Speeeeeeeedssssss
        ! TODO: gamma_const -> gamma_1

        eos_state % rho = Q(QRHO)
        eos_state % p   = Q(QPRES)
        eos_state % xn  = Q(QFS:QFS+nspec-1)

        call eos(eos_input_rp, eos_state)   

	as = eos_state % gam1 * (Q(QPRES) - 0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/Q(QRHO)
	!Alfven
	ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	caz = (Q(QMAGZ)**2)/Q(QRHO)
	!Sloooooooooow
	csz = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))
	!Fassssst
	cfz = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))
	!useful constants
	alf = sqrt((as - csz)/(cfz - csz))
	als = sqrt((cfz - as)/(cfz - csz))
	if(cfz - as .lt. 0.d0) als = 0.d0
	if(abs(Q(QMAGX)).le. 1.d-14 .and.abs(Q(QMAGY)).le. 1.d-14) then
		betx = 1.d0/sqrt(2.d0)
		bety = betx
	else
		betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
		bety = Q(QMAGY)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
	endif
	cff = sqrt(cfz)*alf
	css = sqrt(csz)*als
	S = sign(1.0d0, Q(QMAGZ))
	Qf = sqrt(cfz)*alf*S
	Qs = sqrt(csz)*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
	N = 0.5d0/as
	
!Need to double check the order
	leig(1,:) = (/0.d0, -N*Cff, N*Qs*betx	, N*Qs*bety		, N*alf/Q(QRHO) , N*AAs*betx/Q(QRHO)		    , N*AAs*bety/Q(QRHO)			/) !w - cf
	leig(2,:) = (/0.d0,  0.d0 , -0.5d0*bety , 0.5d0*betx	, 0.d0			, -0.5d0*S*bety/(sqrt(Q(QRHO))) , 0.5d0*betx*S/(sqrt(Q(QRHO)))	/) !w - cAz
	leig(3,:) = (/0.d0, -N*Css, -N*Qf*betx  , -N*Qf*bety	, N*als/Q(QRHO) , -N*AAf*betx/Q(QRHO)		    , -N*AAf*bety/Q(QRHO)			/) !w - cs
	leig(4,:) = (/1.d0,  0.d0 ,  0.d0	    , 0.d0			, -1.d0/as		, 0.d0						    , 0.d0							/) !w
 	leig(5,:) = (/0.d0,  N*Css, N*Qf*betx   , N*Qf*bety		, N*als/Q(QRHO) , -N*AAf*betx/Q(QRHO)		    , -N*AAf*bety/Q(QRHO)			/) !w + cs
	leig(6,:) = (/0.d0,  0.d0 , 0.5d0*bety  , -0.5d0*betx	, 0.d0			, -0.5d0*bety*S/(sqrt(Q(QRHO))) , 0.5d0*betx*S/(sqrt(Q(QRHO)))  /) !w + cAz
	leig(7,:) = (/0.d0, N*Cff, -N*Qs*betx   , -N*Qs*bety	, N*alf/Q(QRHO) , N*AAs*betx/Q(QRHO)			, N*AAs*bety/Q(QRHO)			/) !w + cf
	end subroutine lvecz

!====================================== Right Eigenvectors ===============================================
!x direction
	subroutine rvecx(reig, Q) 
	use amrex_fort_module, only : rt => amrex_real
	use eos_module, only : eos
        use eos_type_module, only: eos_t, eos_input_rp
        use network, only : nspec

	implicit none
	
	!returnes reig, where the cols are the right eigenvectors of the characteristic matrix Ax
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::reig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfx, cax, csx, ca, as, S
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, bety, betz
 
        type(eos_t) :: eos_state
        
	!Speeeeeeeedssssss
        ! TODO: gamma_const -> gamma_1
        eos_state % rho = Q(QRHO)
        eos_state % p   = Q(QPRES)
        eos_state % xn  = Q(QFS:QFS+nspec-1)

        call eos(eos_input_rp, eos_state)    

	as = eos_state % gam1 * (Q(QPRES) - 0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/Q(QRHO)
	!Alfven
	ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	cax = (Q(QMAGX)**2)/Q(QRHO)
	!Sloooooooooow
	csx = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))
	!Fassssst
	cfx = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))
	!useful constants
	alf = sqrt((as - csx)/(cfx - csx))
	als = sqrt((cfx - as)/(cfx - csx))
	if(cfx - as .lt. 0.d0) als = 0.d0
	if(as - csx .lt. 0.d0) alf = 0.d0
	if(abs(Q(QMAGY)).le. 1.d-14 .and.abs(Q(QMAGZ)).le. 1.d-14) then
		bety = 1.d0/sqrt(2.d0)
		betz = bety
	else
		bety = Q(QMAGY)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
		betz = Q(QMAGZ)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
	endif
	cff = sqrt(cfx)*alf
	css = sqrt(csx)*als
	S = sign(1.0d0, Q(QMAGX))
	Qf = sqrt(cfx)*alf*S
	Qs = sqrt(csx)*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
	
				!   u - cf 				u - Cax 				u - cs			u 		u + cs			u + Cax					u + cf
	reig(1,:) = (/	Q(QRHO)*alf		, 0.d0					, Q(QRHO)*als	, 1.d0	, Q(QRHO)*als	, 0.d0					, Q(QRHO)*alf	/)
	reig(2,:) = (/	-cff    	    , 0.d0  				, -css			, 0.d0  , css			, 0.d0  				, 	cff			/)
	reig(3,:) = (/	Qs*bety			, -betz					, -Qf*bety		, 0.d0  , Qf*bety		, betz 					, -Qs*bety		/)
	reig(4,:) = (/	Qs*betz			, bety 					, -Qf*betz		, 0.d0  , Qf*betz		, -bety 				, -Qs*betz		/)
	reig(5,:) = (/	Q(QRHO)*as*alf	, 0.d0	 				, Q(QRHO)*as*als, 0.d0  , Q(QRHO)*as*als, 0.d0  				, Q(QRHO)*as*alf/)
	reig(6,:) = (/	AAs*bety		, -betz*S*sqrt(Q(QRHO))	, -AAf*bety		, 0.d0  , -AAf*bety		, -betz*S*sqrt(Q(QRHO)) , AAs*bety		/)
	reig(7,:) = (/	AAs*betz		, bety*S*sqrt(Q(QRHO))	, -AAf*betz		, 0.d0	, -AAf*betz		, bety*S*sqrt(Q(QRHO))	, AAs*betz		/)
	
	
	end subroutine rvecx

!y direction
	subroutine rvecy(reig, Q) 
	use amrex_fort_module, only : rt => amrex_real
	use eos_module, only : eos
        use eos_type_module, only: eos_t, eos_input_rp
        use network, only : nspec

	implicit none
	
	!returnes reig, where the cols are the right eigenvectors of the characteristic matrix Ay
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::reig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfy, cay, csy, ca, as, S
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, betz
        
        type (eos_t) :: eos_state

	!Speeeeeeeedssssss
        ! TODO: gamma_const -> gamma_1
        eos_state % rho = Q(QRHO)
        eos_state % p   = Q(QPRES)
        eos_state % xn  = Q(QFS:QFS+nspec-1)

        call eos(eos_input_rp, eos_state)   

	as = eos_state % gam1 * (Q(QPRES) - 0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/Q(QRHO)
	!Alfven
	ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	cay = (Q(QMAGY)**2)/Q(QRHO)
	!Sloooooooooow
	csy = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))
	!Fassssst
	cfy = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))
	!useful constants
	alf = sqrt((as - csy)/(cfy - csy))
	if(as - csy .lt. 0.d0) alf = 0.d0
	als = sqrt((cfy - as)/(cfy - csy))
	if(cfy - as .lt. 0.d0) als = 0.d0
	if(abs(Q(QMAGX)).le. 1.d-14 .and.abs(Q(QMAGZ)).le. 1.d-14) then
		betx = 1.d0/sqrt(2.d0)
		betz = betx
	else
		betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
		betz = Q(QMAGZ)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
	endif
	cff = sqrt(cfy)*alf
	css = sqrt(csy)*als
	S = sign(1.0d0, Q(QMAGY))
	Qf = sqrt(cfy)*alf*S
	Qs = sqrt(csy)*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
	
				!   v - cf 				v - Cay 				v - cs			v 		v + cs			v + Cay					v + cf
	reig(1,:) = (/	Q(QRHO)*alf		, 0.d0					, Q(QRHO)*als	, 1.d0	, Q(QRHO)*als	, 0.d0					, Q(QRHO)*alf	/)
	reig(2,:) = (/	-cff    	    , 0.d0  				, -css			, 0.d0  , css			, 0.d0  				, 	cff			/)
	reig(3,:) = (/	Qs*betz			, -betx					, -Qf*betz		, 0.d0  , Qf*betz		, betx 					, -Qs*betz		/)
	reig(4,:) = (/	Qs*betx			, betz 					, -Qf*betx		, 0.d0  , Qf*betx		, -betz 				, -Qs*betx		/)
	reig(5,:) = (/	Q(QRHO)*as*alf	, 0.d0	 				, Q(QRHO)*as*als, 0.d0  , Q(QRHO)*as*als, 0.d0  				, Q(QRHO)*as*alf/)
	reig(6,:) = (/	AAs*betz		, -betx*S*sqrt(Q(QRHO))	, -AAf*betz		, 0.d0  , -AAf*betz		, -betx*S*sqrt(Q(QRHO)) , AAs*betz		/)
	reig(7,:) = (/	AAs*betx		, betz*S*sqrt(Q(QRHO))	, -AAf*betx		, 0.d0	, -AAf*betx		, betz*S*sqrt(Q(QRHO))	, AAs*betx		/)
	
	
	end subroutine rvecy

!z direction
	subroutine rvecz(reig, Q) 
	use amrex_fort_module, only : rt => amrex_real
        use eos_module, only : eos
        use eos_type_module, only: eos_t, eos_input_rp
        use network, only : nspec
	
	implicit none
	
	!returnes reig, where the cols are the right eigenvectors of the characteristic matrix Az
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::reig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfz, caz, csz, ca, as, S
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, bety

        type(eos_t) :: eos_state

	!Speeeeeeeedssssss
        ! TODO: gamma_const -> gamma_1
        eos_state % rho = Q(QRHO)
        eos_state % p   = Q(QPRES)
        eos_state % xn  = Q(QFS:QFS+nspec-1)

        call eos(eos_input_rp, eos_state)   

	as = eos_state % gam1 * (Q(QPRES) - 0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/Q(QRHO)
	!Alfven
	ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	caz = (Q(QMAGZ)**2)/Q(QRHO)
	!Sloooooooooow
	csz = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))
	!Fassssst
	cfz = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))
	!useful constants
	alf = sqrt((as - csz)/(cfz - csz))
	als = sqrt((cfz - as)/(cfz - csz))
	if(cfz - as .lt. 0.d0) als = 0.d0
	if(abs(Q(QMAGX)).le. 1.d-14 .and.abs(Q(QMAGY)).le. 1.d-14) then
		betx = 1.d0/sqrt(2.d0)
		bety = betx
	else
		betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
		bety = Q(QMAGY)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
	endif
	cff = sqrt(cfz)*alf
	css = sqrt(csz)*als
	S = sign(1.0d0, Q(QMAGZ))
	Qf = sqrt(cfz)*alf*S
	Qs = sqrt(csz)*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
				!   w - cf 				w - Caz 				w - cs			w 		w + cs			w + Caz					w + cf
	reig(1,:) = (/	Q(QRHO)*alf		, 0.d0					, Q(QRHO)*als	, 1.d0	, Q(QRHO)*als	, 0.d0					, Q(QRHO)*alf	/)
	reig(2,:) = (/	-cff    	    , 0.d0  				, -css			, 0.d0  , css			, 0.d0  				, 	cff			/)
	reig(3,:) = (/	Qs*betx			, -bety					, -Qf*betx		, 0.d0  , Qf*betx		, bety 					, -Qs*betx		/)
	reig(4,:) = (/	Qs*bety			, betx 					, -Qf*bety		, 0.d0  , Qf*bety		, -betx 				, -Qs*bety		/)
	reig(5,:) = (/	Q(QRHO)*as*alf	, 0.d0	 				, Q(QRHO)*as*als, 0.d0  , Q(QRHO)*as*als, 0.d0  				, Q(QRHO)*as*alf/)
	reig(6,:) = (/	AAs*betx		, -bety*S*sqrt(Q(QRHO))	, -AAf*betx		, 0.d0  , -AAf*betx		, -bety*S*sqrt(Q(QRHO)) , AAs*betx		/)
	reig(7,:) = (/	AAs*bety		, betx*S*sqrt(Q(QRHO))	, -AAf*bety		, 0.d0	, -AAf*bety		, betx*S*sqrt(Q(QRHO))	, AAs*bety		/)
	
	
	end subroutine rvecz
end module mhd_plm_module
