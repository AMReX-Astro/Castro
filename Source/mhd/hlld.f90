module hlld_solver

   implicit none
   public  hlld

contains

subroutine hlld(work_lo, work_hi, &
                qleft, ql_lo, ql_hi, &
                qright, qr_lo, qr_hi, &
                flx, flx_lo, flx_hi, &
                dir)

  !Riemann solve:
  !Main assumption, the normal velocity/Mag field is constant in the Riemann fan, and is sM/Bn respectively. 
  !Total Pressure is constant throughout the Riemann fan, pst!

   use amrex_fort_module, only : rt => amrex_real
   use meth_params_module
   use eos_module, only : eos
   use eos_type_module, only: eos_t, eos_input_rp
   use network, only : nspec

   integer, intent(in)   :: ql_lo(3), ql_hi(3)
   integer, intent(in)   :: qr_lo(3), qr_hi(3)
   integer, intent(in)   :: work_lo(3), work_hi(3)
   integer, intent(in)   :: flx_lo(3), flx_hi(3)
   integer, intent(in)   :: dir

   real(rt), intent(in)  :: qleft(ql_lo(1):ql_hi(1),ql_lo(2):ql_hi(2),ql_lo(3):ql_hi(3),NQ,3)
   real(rt), intent(in)  :: qright(qr_lo(1):qr_hi(1),qr_lo(2):qr_hi(2),qr_lo(3):qr_hi(3),NQ,3)
   real(rt), intent(out) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR+3)

   real(rt)       :: cfL2, cfR, sL, sR, sM, ssL, ssR, pst, caL, canL
   real(rt)       :: caR, canR, asL, asR, ptL, ptR, eL, eR, eintL, eintR
   real(rt)       :: QL(NQ), QR(NQ)
   real(rt)       :: FL(NVAR+3), FR(NVAR+3)
   real(rt)       :: uL(NVAR+3), uR(NVAR+3)
   real(rt)       :: UsL(NVAR+3), FsL(NVAR+3)
   real(rt)       :: UsR(NVAR+3), FsR(NVAR+3)
   real(rt)       :: UssL(NVAR+3), FssL(NVAR+3)
   real(rt)       :: UssR(NVAR+3), FssR(NVAR+3)
   real(rt) :: gam1_L, gam1_R

   integer           :: QVELN, QVELP1, QVELP2
   integer           :: QMAGN, QMAGP1, QMAGP2
   integer           :: UMN  , UMP1  , UMP2
   integer           :: UMAGN, UMAGP1, UMAGP2
   integer           :: UMAGX, UMAGY, UMAGZ
   integer           :: i,j,k
   character(len=10) :: choice

   type (eos_t) :: eos_state

   ! in the conserved part, let's have magx, magy and magz at end
   UMAGX = NVAR+1
   UMAGY = NVAR+2
   UMAGZ = NVAR+3
 
   ! `n` here is the normal
   ! `p` are the perpendicular
   
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
      UMAGN  = UMAGX    
      UMAGP1 = UMAGY
      UMAGP2 = UMAGZ
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
      UMAGN  = UMAGY
      UMAGP1 = UMAGZ
      UMAGP2 = UMAGX
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
      UMAGN  = UMAGZ
      UMAGP1 = UMAGX
      UMAGP2 = UMAGY
   end if

   do k = work_lo(3), work_hi(3)
   do j = work_lo(2), work_hi(2)
   do i = work_lo(1), work_hi(1)

      
      if (dir .eq. 1) then
         qL(:) = qleft(i-1,j,k,:,dir)
      else if (dir .eq. 2) then
         qL(:) = qleft(i,j-1,k,:,dir)
      else if (dir .eq. 3) then
         qL(:) = qleft(i,j,k-1,:,dir)
      end if

      qR(:) = qright(i,j,k,:,dir)

      flx(i,j,k,:) = 0.d0  
      FL  = 0.d0; FR = 0.d0; UsL = 0.d0; UsR = 0.d0; FsL = 0.d0; FsR = 0.d0; UssL = 0.d0; UssR = 0.d0; FssL = 0.d0; FssR = 0.d0

      !check to not go lower than small_dens and small_press
      qL(QRHO) = max(small_dens, qL(QRHO))
      qR(QRHO) = max(small_dens, qR(QRHO))
      qL(QPRES) = max(small_pres, qL(QPRES))
      qR(QPRES) = max(small_pres, qR(QPRES))


      call PToC(qL, uL)
      call PToC(qR, uR)
      
      eos_state % rho = qL(QRHO)
      eos_state % p   = qL(QPRES)
      eos_state % T   = 100.0 !dummy initial guess 
      eos_state % xn  = qL(QFS:QFS+nspec-1)

      ! Note this is actually (rho E)
      call eos(eos_input_rp, eos_state)

      eL   = eos_state % rho * eos_state % e & 
                        + 0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)) &
                            + 0.5d0*dot_product(qL(QU:QW),qL(QU:QW))*qL(QRHO)
      eintL = eos_state % e
      gam1_L = eos_state % gam1
!here use total p not just p_g eq.11 in Miniati
      FL(URHO)  = qL(QRHO)*qL(QVELN)
      FL(UMN)   = qL(QRHO)*qL(QVELN)**2 + (qL(QPRES)+ 0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)))&
                  - qL(QMAGN)**2
      FL(UMP1)  = qL(QRHO)*qL(QVELN)*qL(QVELP1) - qL(QMAGN)*qL(QMAGP1)
      FL(UMP2)  = qL(QRHO)*qL(QVELN)*qL(QVELP2) - qL(QMAGN)*qL(QMAGP2)
      FL(UEDEN) = qL(QVELN)*(eL + (qL(QPRES)+ 0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)))) &
                  - qL(QMAGN)*dot_product(qL(QMAGX:QMAGZ),qL(QU:QW))
      FL(UMAGN) = 0.d0
      FL(UMAGP1) = qL(QVELN)*qL(QMAGP1) - qL(QVELP1)*qL(QMAGN)  
      FL(UMAGP2) = qL(QVELN)*qL(QMAGP2) - qL(QVELP2)*qL(QMAGN) 
      FL(UFS:UFS+nspec-1) = qL(QRHO)*qL(QVELN)*qL(QFS:QFS+nspec-1) 
      FL(UEINT) = qL(QRHO)*qL(QVELN)*eintL

      eos_state % rho = qR(QRHO)
      eos_state % p   = qR(QPRES)
      eos_state % T   = 100.0 ! dummy initial guess?  
      eos_state % xn  = qR(QFS:QFS+nspec-1)

      call eos(eos_input_rp, eos_state)

      eR   = eos_state % rho * eos_state % e &
                        + 0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)) &
                        + 0.5d0*dot_product(qR(QU:QW),qR(QU:QW))*qR(QRHO)
      eintR = eos_state % e
      gam1_R = eos_state % gam1

      FR(URHO)  = qR(QRHO)*qR(QVELN)
      FR(UMN)   = qR(QRHO)*qR(QVELN)**2 + (qR(QPRES)+0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ))) &
                  - qR(QMAGN)**2
      FR(UMP1)  = qR(QRHO)*qR(QVELN)*qR(QVELP1) - qR(QMAGN)*qR(QMAGP1)
      FR(UMP2)  = qR(QRHO)*qR(QVELN)*qR(QVELP2) - qR(QMAGN)*qR(QMAGP2)
      FR(UEDEN) = qR(QVELN)*(eR + (qR(QPRES)+0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)))) &
                  - qR(QMAGN)*dot_product(qR(QMAGX:QMAGZ),qR(QU:QW))
      FR(UMAGN) = 0.d0
      FR(UMAGP1) = qR(QVELN)*qR(QMAGP1) - qR(QVELP1)*qR(QMAGN)   
      FR(UMAGP2) = qR(QVELN)*qR(QMAGP2) - qR(QVELP2)*qR(QMAGN) 
      FR(UFS:UFS+nspec-1) = qR(QRHO)*qR(QVELN)*qR(QFS:QFS+nspec-1) 
      FR(UEINT) = qR(QRHO)*qR(QVELN)*eintR

      ! From Miyoshi and Kusano paper eq.(3) 
        asL  = gam1_L * qL(QPRES)/qL(QRHO)
        asR  = gam1_R * qR(QPRES)/qR(QRHO)

        caL  = (qL(QMAGN)**2 + qL(QMAGP1)**2 + qL(QMAGP2)**2)/qL(QRHO) !Magnetic Speeds
        caR  = (qR(QMAGN)**2 + qR(QMAGP1)**2 + qR(QMAGP2)**2)/qR(QRHO)

        canL = (qL(QMAGN)**2)/qL(QRHO)
        canR = (qR(QMAGN)**2)/qR(QRHO)

        !Catch the fastest waves, brah
        !cf as in equation (3) of HLLD paper (Miyoshi and Kusano)
        cfL2  = sqrt(0.5d0*((asL + caL) + sqrt((asL + caL)**2 - 4.0d0*asL*canL)))
        cfR  = sqrt(0.5d0*((asR + caR) + sqrt((asR + caR)**2 - 4.0d0*asR*canR)))

        !Riemann Speeds
        !sL and sR, eq.(12) (Miyoshi and Kusano)
        sL   = min(qL(QVELN) - cfL2, qR(QVELN) - cfR)
        sR   = max(qL(QVELN) + cfL2, qR(QVELN) + cfR)
        !sM, eq.(38)
        sM   = ((sR - qR(QVELN))*qR(QRHO)*qR(QVELN) - (sL - qL(QVELN))*qL(QRHO)*qL(QVELN) - &
               (qR(QPRES)+0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ))) + &
               (qL(QPRES)+0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ))))/ &
               ((sR - qR(QVELN))*qR(QRHO) - (sL - qL(QVELN))*qL(QRHO))

        !Pressures in the Riemann Fan
        ptL  = qL(QPRES)+0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ))
        ptR  = qR(QPRES)+0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ))
        !pst, eq.(41)
        pst  = (sR - qR(QVELN))*qR(QRHO)*ptL - (sL - qL(QVELN))*qL(QRHO)*ptR + qL(QRHO)*qR(QRHO)*(sR - qR(QVELN))*(sL - qL(QVELN))*(qR(QVELN) - qL(QVELN))
        pst  = pst/((sR - qR(QVELN))*qR(QRHO) - (sL - qL(QVELN))*qL(QRHO))

        !------------------------------------------- Density * states-------------------------------------------------------------------------

        ! Density eq.(43)
        UsL(URHO) = qL(QRHO)*((sL - qL(QVELN))/(sL - sM))
        UsR(URHO) = qR(QRHO)*((sR - qR(QVELN))/(sR - sM))
        
        !------------------------------------------- Species * states ------------------------------------------------

        UsL(UFS:UFS+nspec-1) = qL(QFS:QFS+nspec-1)*UsL(URHO)
        UsR(UFS:UFS+nspec-1) = qR(QFS:QFS+nspec-1)*UsR(URHO)

        ! e   
        UsL(UEINT) = eintL * UsL(URHO)
        UsR(UEINT) = eintR * UsR(URHO)

        !--------------------------------------------------------- Vel * states----------------------------------------------------------------

        !Normal dir
        UsL(UMN)    = sM
        UsR(UMN)    = sM

        !Perpendicular dir
        if(abs(qL(QMAGN)*qL(QMAGP1)*(sM - qL(QVELN))).lt. 1d-14) then
                UsL(UMP1)    = qL(QVELP1)
        else
                UsL(UMP1)    = qL(QVELP1) - qL(QMAGN)*qL(QMAGP1)*((sM - qL(QVELN))/(qL(QRHO)*(sL - qL(QVELN))*(sL - sM) - qL(QMAGN)**2))
        endif
        if(abs(qR(QMAGN)*qR(QMAGP1)*(sM - qR(QVELN))).lt. 1d-14) then
                UsR(UMP1) = qR(QVELP1)
        else
                UsR(UMP1)    = qR(QVELP1) - qR(QMAGN)*qR(QMAGP1)*((sM - qR(QVELN))/(qR(QRHO)*(sR - qR(QVELN))*(sR - sM) - qR(QMAGN)**2))
        endif
        
        ! Second Perpendicular dir
        if(abs(qL(QMAGN)*qL(QMAGP2)*(sM - qL(QVELN))).le. 1d-14) then
                UsL(UMP2) = qL(QVELP2)
        else
                UsL(UMP2)    = qL(QVELP2) - qL(QMAGN)*qL(QMAGP2)*((sM - qL(QVELN))/(qL(QRHO)*(sL - qL(QVELN))*(sL - sM) - qL(QMAGN)**2))
        endif
        if(abs(qR(QMAGN)*qR(QMAGP2)*(sM - qR(QVELN))).le. 1d-14) then
                UsR(UMP2) = qR(QVELP2)
        else
                UsR(UMP2)    = qR(QVELP2) - qR(QMAGN)*qR(QMAGP2)*((sM - qR(QVELN))/(qR(QRHO)*(sR - qR(QVELN))*(sR - sM) - qR(QMAGN)**2))
        endif

        UsL(UMX:UMZ) = UsL(UMX:UMZ)*UsL(URHO)
        UsR(UMX:UMZ) = UsR(UMX:UMZ)*UsR(URHO)

        !--------------------------------------------------------- B * states ----------------------------------------------------------------------
        
        !Magnetic Fields
        
        !Normal dir
        UsL(UMAGN) = qL(QMAGN)
        UsR(UMAGN) = qR(QMAGN) 

        !Perpendicular dir
        if(abs(qL(QMAGP1)*(qL(QRHO)*(sL - qL(QVELN))**2 - qL(QMAGN)**2)).lt.1d-14) then
                UsL(UMAGP1) = qL(QMAGP1)
        else
                UsL(UMAGP1) = qL(QMAGP1)*(qL(QRHO)*(sL - qL(QVELN))**2 - qL(QMAGN)**2)/(qL(QRHO)*(sL - qL(QVELN))*(sL - sM) - qL(QMAGN)**2)
        endif
        if(abs(qR(QMAGP1)*(qR(QRHO)*(sR - qR(QVELN))**2 - qL(QMAGN)**2)).lt.1d-14) then 
                UsR(UMAGP1) = qR(QMAGP1)
        else
                UsR(UMAGP1) = qR(QMAGP1)*(qR(QRHO)*(sR - qR(QVELN))**2 - qR(QMAGN)**2)/(qR(QRHO)*(sR - qR(QVELN))*(sR - sM) - qR(QMAGN)**2)
        endif

        ! Second Perpendicular dir
        if(abs(qL(QMAGP2)*(qL(QRHO)*(sL - qL(QVELN))**2 - qL(QMAGN)**2)).lt. 1d-14) then
                UsL(UMAGP2) = qL(QMAGP2)
        else
                UsL(UMAGP2) = qL(QMAGP2)*(qL(QRHO)*(sL - qL(QVELN))**2 - qL(QMAGN)**2)/(qL(QRHO)*(sL - qL(QVELN))*(sL - sM) - qL(QMAGN)**2)
        endif
        if(abs(qR(QMAGP2)*(qR(QRHO)*(sR - qR(QVELN))**2 - qR(QMAGN)**2)).lt.1d-14) then 
                UsR(UMAGP2) = qR(QMAGP2)
        else
                UsR(UMAGP2) = qR(QMAGP2)*(qR(QRHO)*(sR - qR(QVELN))**2 - qR(QMAGN)**2)/(qR(QRHO)*(sR - qR(QVELN))*(sR - sM) - qR(QMAGN)**2)
        endif
        
        !Energy, eq.(48)
        UsL(UEDEN) = (sL - qL(QVELN))*eL - ptL*qL(QVELN) + pst*sM + &
                      qL(QMAGN)*(dot_product(qL(QU:QW),qL(QMAGX:QMAGZ)) &
                                 -dot_product(UsL(UMX:UMZ)/UsL(URHO),UsL(UMAGX:UMAGZ)))
        UsL(UEDEN) = UsL(UEDEN)/(sL - sM)
        UsR(UEDEN) = (sR - qR(QVELN))*eR - ptR*qR(QVELN) + pst*sM + &
                      qR(QMAGN)*( dot_product(qR(QU:QW),qR(QMAGX:QMAGZ)) &
                                 -dot_product(UsR(UMX:UMZ)/UsR(URHO),UsR(UMAGX:UMAGZ)))
        UsR(UEDEN) = UsR(UEDEN)/(sR - sM)

        !speeds, eq.(51)
        ssL = sM - abs(qL(QMAGN))/sqrt(UsL(URHO))
        ssR = sM + abs(qR(QMAGN))/sqrt(UsR(URHO))

        !----------------------------------------- ** states ------------------------------------------------------------------------------
        !Dens
        UssL(URHO)  = UsL(URHO)
        UssR(URHO)  = UsR(URHO)

        !species
        UssL(UFS:UFS+nspec-1) = UsL(UFS:UFS+nspec-1)
        UssR(UFS:UFS+nspec-1) = UsR(UFS:UFS+nspec-1)
       
        !e
        UssL(UEINT) = UsL(UEINT)
        UssR(UEINT) = UsR(UEINT)
        
        !Vel in normal direction
        UssL(UMN)    = sM
        UssR(UMN)    = sM

        !Vel in perpendicular directioni, eq.(59)
        UssL(UMP1)    = (sqrt(UsL(URHO))*UsL(UMP1)/UsL(URHO) + sqrt(UsR(URHO))*UsR(UMP1)/UsR(URHO)&
        + (UsR(UMAGP1) - UsL(UMAGP1))*sign(1.d0,qL(QMAGN)))/(sqrt(UsL(URHO)) + sqrt(UsR(URHO)))
        UssR(UMP1)    = UssL(UMP1)

        !Vel in second perpendicular direction, eq(60)
        UssL(UMP2)    = (sqrt(UsL(URHO))*UsL(UMP2)/UsL(URHO) + sqrt(UsR(URHO))*UsR(UMP2)/UsR(URHO)&
                                + (UsR(UMAGP2) - UsL(UMAGP2))*sign(1.d0,qL(QMAGN)))/(sqrt(UsL(URHO)) + sqrt(UsR(URHO)))

        UssR(UMP2)    = UssL(UMP2)

        UssL(UMX:UMZ) = UssL(UMX:UMZ)*UssL(URHO)
        UssR(UMX:UMZ) = UssR(UMX:UMZ)*UssR(URHO)
        
        !B in normal direction
        UssL(UMAGN) = UsL(UMAGN)
        UssR(UMAGN) = UsR(UMAGN)

        !B in perpendicular direction
        UssL(UMAGP1) = (sqrt(UsL(URHO))*UsR(UMAGP1) + sqrt(UsR(URHO))*UsL(UMAGP1) + sqrt(UsL(URHO)*UsR(URHO))& 
        *(UsR(UMP1)/UsR(URHO) - UsL(UMP1)/UsL(URHO))*sign(1.d0,qL(QMAGN)))/(sqrt(UsL(URHO)) + sqrt(UsR(URHO)))
        UssR(UMAGP1) = UssL(UMAGP1)

        !B in second perpendicular direction
        UssL(UMAGP2) = (sqrt(UsL(URHO))*UsR(UMAGP2) + sqrt(UsR(URHO))*UsL(UMAGP2) + sqrt(UsL(URHO)*UsR(URHO))&
          *(UsR(UMP2)/UsR(URHO) - UsL(UMP2)/UsL(URHO))*sign(1.d0,qL(QMAGN)))/(sqrt(UsL(URHO)) + sqrt(UsR(URHO)))
        UssR(UMAGP2) = UssL(UMAGP2)

        !Energy , eq.(63)
        UssL(UEDEN) = UsL(UEDEN) - sqrt(UsL(URHO))*(dot_product(UsL(UMX:UMZ)/UsL(URHO),UsL(UMAGX:UMAGZ)) &
        - dot_product(UssL(UMX:UMZ)/UssL(URHO),UssL(UMAGX:UMAGZ)))*sign(1.d0, qL(QMAGN))
        UssR(UEDEN) = UsR(UEDEN) + sqrt(UsR(QRHO))*(dot_product(UsR(UMX:UMZ)/UsR(URHO),UsR(UMAGX:UMAGZ)) &
        - dot_product(UssR(UMX:UMZ)/UssR(URHO),UssR(UMAGX:UMAGZ)))*sign(1.d0, qR(QMAGN))

        !--------------------------------------------------------- Fluxes ----------------------------------------------------------------------
       ! eq. (64)
        FsL  = FL + sL*UsL - sL*uL
        !eq. (65)
        FssL = FL + ssL*UssL - (ssL-sL)*UsL - sL*uL
        FsR  = FR + sR*UsR - sR*uR
        FssR = FR + ssR*UssR - (ssR-sR)*UsR - sR*uR

        !Solve the RP
        if(sL .gt. 0.d0) then
           flx(i,j,k,:) = FL
           choice = "FL"
        elseif(sL .le. 0.d0 .and. ssL .gt. 0.d0) then
           flx(i,j,k,:) = FsL
           choice = "FsL"
        elseif(ssL .le. 0.d0 .and. sM .gt. 0.d0) then
           flx(i,j,k,:) = FssL
           choice = "FssL"
        elseif(sM .le. 0.d0 .and. ssR .gt. 0.d0) then
           flx(i,j,k,:) = FssR
           choice = "FssR"
        elseif(ssR .le. 0.d0 .and. sR .gt. 0.d0) then
           flx(i,j,k,:) = FsR
           choice = "FsR"
        else 
           flx(i,j,k,:) = FR
           choice = "FR"
        endif


   end do
   end do
   end do

end subroutine hlld

!================================================= Calculate the Conservative Variables ===============================================

subroutine PToC(q, u)

   use amrex_fort_module, only : rt => amrex_real
   use meth_params_module
   use eos_module, only : eos
   use eos_type_module, only: eos_t, eos_input_rp
   use network, only : nspec

   implicit none

   real(rt), intent(in)  ::q(NQ)
   real(rt), intent(out) ::u(NVAR+3)
   
   type (eos_t) :: eos_state

   integer :: n

   u(:) = 0.d0 

   u(URHO)       = q(QRHO)
   u(UMX)        = q(QRHO)*q(QU)
   u(UMY)        = q(QRHO)*q(QV)
   u(UMZ)        = q(QRHO)*q(QW)

   eos_state % rho = q(QRHO)
   eos_state % p   = q(QPRES) 
   eos_state % xn  = q(QFS:QFS+nspec-1)
   eos_state % T   = 100.0   ! dummy initial guess

   call eos(eos_input_rp, eos_state)


  u(UEINT)       = eos_state % rho * eos_state % e
  u(UEDEN)       = u(UEINT)  + 0.5d0*q(QRHO)*dot_product(q(QU:QW),q(QU:QW)) &
                             + 0.5d0*(dot_product(q(QMAGX:QMAGZ),q(QMAGX:QMAGZ)))
  u(NVAR+1:NVAR+3) = q(QMAGX:QMAGZ)

  ! handle species too
  do n = 1, nspec
     u(UFS-1+n) = u(URHO)*q(QFS-1+n)
  enddo

end subroutine PToC

end module hlld_solver
