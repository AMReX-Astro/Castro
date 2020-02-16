module hll_solver

   implicit none
   public  hll

contains

subroutine hll(work_lo, work_hi, qm,qp,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3,dir)

  !Riemann solve:
  !Main assumption, the normal velocity/Mag field is constant in the Riemann fan, and is sM/Bn respectively. 
  !Total Pressure is constant throughout the Riemann fan, pst!

   use amrex_fort_module, only : rt => amrex_real
   use meth_params_module
   use eos_module, only : eos
   use eos_type_module, only: eos_t, eos_input_rp
   use network, only : nspec

   integer, intent(in)   :: q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
   integer, intent(in)   :: work_lo(3), work_hi(3)
   integer, intent(in)   :: flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
   integer, intent(in)   :: dir

   real(rt), intent(in)  :: qm(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
   real(rt), intent(in)  :: qp(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
   real(rt), intent(out) :: flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR)

   real(rt)       :: cfL2, cfR, sL, sR, sM, ssL, ssR, pst, caL, canL
   real(rt)       :: caR, canR, asL, asR, ptL, ptR, eL, eR
   real(rt)       :: QL(QVAR), QR(QVAR)
   real(rt)       :: uL(QVAR), uR(QVAR)
   real(rt)       :: FL(QVAR), FR(QVAR)
   real(rt)       :: uS(QVAR), FS(QVAR)

   integer           :: QVELN, QVELP1, QVELP2
   integer           :: QMAGN, QMAGP1, QMAGP2
   integer           :: UMN  , UMP1  , UMP2
   integer           :: i,j,k
   
   type (eos_t) :: eos_state

   character(len=10) :: choice
   
   if (dir .eq. 1) then
      QMAGN  = QMAGX
      QMAGP1 = QMAGY
      QMAGP2 = QMAGZ
      QVELN  = QU
      QVELP1 = QV
      QVELP2 = QW
      UMN    = UMX
      UMP1   = UMY
      UMP2   = UMZ
   else if (dir .eq. 2) then
      QMAGN  = QMAGY
      QMAGP1 = QMAGZ
      QMAGP2 = QMAGX
      QVELN  = QV
      QVELP1 = QW
      QVELP2 = QU
      UMN    = UMY
      UMP1   = UMZ
      UMP2   = UMX
   else if (dir .eq. 3) then
      QMAGN  = QMAGZ
      QMAGP1 = QMAGX
      QMAGP2 = QMAGY
      QVELN  = QW
      QVELP1 = QU
      QVELP2 = QV
      UMN    = UMZ
      UMP1   = UMX
      UMP2   = UMY
   end if

   do k = work_lo(3), work_hi(3)
   do j = work_lo(2), work_hi(2)
   do i = work_lo(1), work_hi(1)

      
      if (dir .eq. 1) then
         qL(:) = qp(i-1,j,k,:,dir)
      else if (dir .eq. 2) then
         qL(:) = qp(i,j-1,k,:,dir)
      else if (dir .eq. 3) then
         qL(:) = qp(i,j,k-1,:,dir)
      end if

      qR(:) = qm(i,j,k,:,dir)

      flx(i,j,k,:) = 0.d0  
      FL  = 0.d0; FR = 0.d0; uS = 0.d0; FS = 0.d0;

      call PToC(qL,uL)
      call PToC(qR,uR)
      ! Note this is actually (rho e)
      eos_state % rho = qL(QRHO)
      eos_state % p   = qL(QPRES)
      eos_state % xn  = qL(QFS:QFS+nspec-1) 

      call eos(eos_input_rp, eos_state)

      eL   = eos_state % rho * eos_state % e &
                        + 0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)) &
                            + 0.5d0*dot_product(qL(QU:QW),qL(QU:QW))*qL(QRHO)

      FL(URHO)  = qL(QRHO)*qL(QVELN)
      FL(UMN)   = qL(QRHO)*qL(QVELN)**2 + qL(QPRES) - qL(QMAGN)**2
      FL(UMP1)  = qL(QRHO)*qL(QVELN)*qL(QVELP1) - qL(QMAGN)*qL(QMAGP1)
      FL(UMP2)  = qL(QRHO)*qL(QVELN)*qL(QVELP2) - qL(QMAGN)*qL(QMAGP2)
      FL(UEDEN) = qL(QVELN)*(eL + qL(QPRES)) - qL(QMAGN)*dot_product(qL(QMAGX:QMAGZ),qL(QU:QW))
      FL(QMAGN) = 0.d0
      FL(QMAGP1) = qL(QVELN)*qL(QMAGP1) - qL(QVELP1)*qL(QMAGN) 
      FL(QMAGP2) = qL(QVELN)*qL(QMAGP2) - qL(QVELP2)*qL(QMAGN)

      ! Note this is actually (rho e)
      ! TODO: need to get rho e from the EOS using p, rho, x
      eos_state % rho = qR(QRHO)
      eos_state % p   = qR(QPRES)
      eos_state % xn  = qR(QFS:QFS+nspec-1)

      call eos(eos_input_rp, eos_state)
      eR   = eos_state % rho * eos_state % e &
                        + 0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)) &
                        + 0.5d0*dot_product(qR(QU:QW),qR(QU:QW))*qR(QRHO)

      FR(URHO)  = qR(QRHO)*qR(QVELN)
      FR(UMN)   = qR(QRHO)*qR(QVELN)**2 + qR(QPRES) - qR(QMAGN)**2
      FR(UMP1)  = qR(QRHO)*qR(QVELN)*qR(QVELP1) - qR(QMAGN)*qR(QMAGP1)
      FR(UMP2)  = qR(QRHO)*qR(QVELN)*qR(QVELP2) - qR(QMAGN)*qR(QMAGP2)
      FR(UEDEN) = qR(QVELN)*(eR + qR(QPRES)) - qR(QMAGN)*dot_product(qR(QMAGX:QMAGZ),qR(QU:QW))
      FR(QMAGN) = 0.d0
      FR(QMAGP1) = qR(QVELN)*qR(QMAGP1) - qR(QVELP1)*qR(QMAGN)  
      FR(QMAGP2) = qR(QVELN)*qR(QMAGP2) - qR(QVELP2)*qR(QMAGN)

        uS = sR*uR - sL*uL - FR + FL
        uS = uS/(sR - sL)

        !--------------------------------------------------------- Fluxes ----------------------------------------------------------------------

        FS  = sR*FL - sL*FR + sR*sL*(uR - uL)
        FS  = FS/(sR - sL)

        !Solve the RP
        if(sL .gt. 0.d0) then
           flx(i,j,k,:) = FL
           choice = "FL"
        elseif(sL .le. 0.d0 .and. sR .ge. 0.d0) then 
           flx(i,j,k,:) = FS
           choice = "FS"
        else 
           flx(i,j,k,:) = FR
           choice = "FR"
        endif

        !if(dir.eq.2.and.((i.eq.3.and.j.eq.16.and.k.eq.1).or.(i.eq.4.and.j.eq.16.and.k.eq.1))) then
        !       print*, "dir, i, j, k =", dir, i, j, k
        !       print*, "flux is ", choice, " = ", flx(i,j,k,UMX:UMZ), flx(i,j,k,QMAGX:QMAGZ)
        !       print*, "FL = ",  FL(UMX:UMZ), FL(QMAGX:QMAGZ)
        !       print*, "FsL = ", FsL(UMX:UMZ), FsL(QMAGX:QMAGZ)
        !       print*, "ssL = ", ssL, "sL = ", sL, "sM =", sM
        !       print*, "UsL = ", UsL(UMX:UMZ), UsL(QMAGX:QMAGZ) 
        !       print*, "UssL = ", UssL(UMX:UMZ), UssL(QMAGX:QMAGZ) 
        !       print*, "qL = ", qL(QMAGX:QMAGZ)
        !       print*, "qR = ", qR(QMAGX:QMAGZ)
        !       pause
        !endif
   end do
   end do
   end do

end subroutine hll

!================================================= Calculate the Conservative Variables ===============================================

subroutine PToC(q, u)

   use amrex_fort_module, only : rt => amrex_real
   use meth_params_module
   use eos_module, only : eos
   use eos_type_module, only: eos_t, eos_input_rp
   use actual_network, only : nspec

   implicit none

   real(rt), intent(in)  ::q(QVAR)
   real(rt), intent(out) ::u(QVAR)
   
   type(eos_t) :: eos_state

   u = 0.d0

   u(URHO)       = q(QRHO)
   u(UMX)        = q(QRHO)*q(QU)
   u(UMY)        = q(QRHO)*q(QV)
   u(UMZ)        = q(QRHO)*q(QW)
   ! TODO: need to get this (rho e) from the EOS using p, rho, X
   eos_state % rho = q(QRHO)
   eos_state % p   = q(QPRES)
   eos_state % xn  = q(QFS:QFS+nspec-1)

   call eos(eos_input_rp, eos_state)

   u(UEINT)       = eos_state % rho * eos_state % e
   u(UEDEN)       = u(UEINT)  + 0.5d0*q(QRHO)*dot_product(q(QU:QW),q(QU:QW)) &
                             + 0.5d0*(dot_product(q(QMAGX:QMAGZ),q(QMAGX:QMAGZ)))
   u(QMAGX:QMAGZ) = q(QMAGX:QMAGZ)

end subroutine PToC

end module hll_solver
