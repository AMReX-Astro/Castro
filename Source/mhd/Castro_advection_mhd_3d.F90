module mhd_advection_module

  implicit none

contains

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine ctoprim_mhd(lo, hi, &
                       uin, uin_lo, uin_hi, &
                       bcc, bcc_lo, bcc_hi, &
                       bx, bxin_lo, bxin_hi, &
                       by, byin_lo, byin_hi, &
                       bz, bzin_lo, bzin_hi, &
                       q, q_lo, q_hi, &
                       cx, cx_lo, cx_hi, &
                       cy, cy_lo, cy_hi, &
                       cz, cz_lo, cz_hi) bind(C, name="ctoprim_mhd")

  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec, naux
  !      use eos_params_module
  use eos_module
  use eos_type_module, only : eos_t, eos_input_re
  !      use flatten_module
  use amrex_constants_module
  use amrex_error_module
  use meth_params_module, only : URHO, UMX, UMY, UMZ, NVAR,&
       UEDEN, UEINT, UFA, UFS, UTEMP, &
       NQ, NSRC ,NQSRC, QRHO, QU, QV, QW, QC, NQAUX,&
       QREINT, QPRES, QFA, QFS, QTEMP, &
       QMAGX,  QMAGY, QMAGZ, &
       nadv, small_dens, small_pres, &
       npassive, upass_map, qpass_map

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: bcc_lo(3), bcc_hi(3)
  integer, intent(in) :: bxin_lo(3), bxin_hi(3)
  integer, intent(in) :: byin_lo(3), byin_hi(3)
  integer, intent(in) :: bzin_lo(3), bzin_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: cx_lo(3), cx_hi(3)
  integer, intent(in) :: cy_lo(3), cy_hi(3)
  integer, intent(in) :: cz_lo(3), cz_hi(3)

  real(rt), intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
  real(rt), intent(in) :: bx(bxin_lo(1):bxin_hi(1), bxin_lo(2):bxin_hi(2), bxin_lo(3):bxin_hi(3))
  real(rt), intent(in) :: by(byin_lo(1):byin_hi(1), byin_lo(2):byin_hi(2), byin_lo(3):byin_hi(3))
  real(rt), intent(in) :: bz(bzin_lo(1):bzin_hi(1), bzin_lo(2):bzin_hi(2), bzin_lo(3):bzin_hi(3))

  real(rt), intent(inout) :: bcc(bcc_lo(1):bcc_hi(1), bcc_lo(2):bcc_hi(2), bcc_lo(3):bcc_hi(3), 3)
  real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
  real(rt), intent(inout) :: cx(cx_lo(1):cx_hi(1), cx_lo(2):cx_hi(2), cx_lo(3):cx_hi(3))
  real(rt), intent(inout) :: cy(cy_lo(1):cy_hi(1), cy_lo(2):cy_hi(2), cy_lo(3):cx_hi(3))
  real(rt), intent(inout) :: cz(cz_lo(1):cz_hi(1), cz_lo(2):cz_hi(2), cz_lo(3):cx_hi(3))

  integer :: i, j, k
  integer :: n, nq2, ipassive

  real(rt) :: cad
  real(rt) :: rhoInv

  ! for the soundspeed
  real(rt) :: as, ca

  type(eos_t) :: eos_state


  !
  ! Make q (all but p), except put e in slot for rho.e, fix after eos call.
  ! The temperature is used as an initial guess for the eos call and will be overwritten.
  !
  ! Also calculate cell-centered magnetic field

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           q(i,j,k,QMAGX) = 0.5d0*(bx(i+1,j,k) + bx(i,j,k))
           bcc(i,j,k, 1) = q(i,j,k,QMAGX)

           q(i,j,k,QMAGY) = 0.5d0*(by(i,j+1,k) + by(i,j,k))
           bcc(i,j,k,2) = q(i,j,k,QMAGY)

           q(i,j,k,QMAGZ) = 0.5d0*(bz(i,j,k+1) + bz(i,j,k))
           bcc(i,j,k,3) = q(i,j,k,QMAGZ)

           if (uin(i,j,k,URHO) .le. ZERO) then
              !
              ! A critical region since we usually can't write from threads.
              !
              print *,'   '
              print *,'>>> Error: Castro_advection_mhd_3d::ctoprim ',i,j,k
              print *,'>>> ... negative density ',uin(i,j,k,URHO)
              call amrex_error("Error:: Castro_advection_mhd_3d.f90 :: ctoprim")
           elseif (uin(i,j,k,URHO) .le. small_dens) then 
              print *,'   '
              print *,'>>> Error: Castro_advection_mhd_3d::ctoprim ',i,j,k
              print *,'>>> ... small density ',uin(i,j,k,URHO)
              call amrex_error("Error:: Castro_advection_mhd_3d.f90 :: ctoprim") 
           endif

           rhoInv = ONE/uin(i,j,k,URHO)

           q(i,j,k,QRHO) = uin(i,j,k,URHO)
           q(i,j,k,QU)   = uin(i,j,k,UMX)*rhoInv
           q(i,j,k,QV)   = uin(i,j,k,UMY)*rhoInv
           q(i,j,k,QW)   = uin(i,j,k,UMZ)*rhoInv

           q(i,j,k,QFS:QFS+nspec-1) = uin(i,j,k,UFS:UFS+nspec-1)*rhoInv
           ! Convert "rho e" to "e"
           q(i,j,k,QREINT ) = uin(i,j,k,UEINT)*rhoInv

           ! Get p, T, c, using q state

           ! If necessary, reset the energy using small_temp ?
           ! should we use the dual energy formalism like in the regular hydro?

           if (q(i,j,k,QREINT) .lt. ZERO) then
              !
              ! A critical region since we usually can't write from threads.
              !
              print *,'   '
              print *,'>>> Error: Castro_advection_mhd_3d::ctoprim ',i,j,k
              print *,'>>> ... e is negative ',q(i,j,k,QREINT)
              print *,'    '
              call amrex_error("Error:: Castro_advection_mhd_3d.f90 :: ctoprim")
           end if

           ! Define the magneto-accoustic speed from the EOS
           ! Note: QREINT is currently just "e"

           eos_state % rho = q(i,j,k,QRHO)
           eos_state % e   = q(i,j,k, QREINT)
           eos_state % T   = uin(i,j,k,UTEMP)
           eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)

           call eos(eos_input_re, eos_state)

           q(i,j,k,QPRES) = eos_state % p

           q(i,j,k,QTEMP) = eos_state % T

           ! Convert "e" back to "rho e"
           q(i,j,k,QREINT) = q(i,j,k,QREINT)*q(i,j,k,QRHO)

           !sound speed for ideal mhd
           as = eos_state % gam1 * q(i,j,k,QPRES)/q(i,j,k,QRHO)
           ca = (q(i,j,k,QMAGX)**2 + q(i,j,k,QMAGY)**2 + q(i,j,k,QMAGZ)**2)/q(i,j,k,QRHO)

           cad = q(i,j,k,QMAGX)**2/q(i,j,k,QRHO)
           call eos_soundspeed_mhd(cx(i,j,k), as, ca, cad)

           cad = q(i,j,k,QMAGY)**2/q(i,j,k,QRHO)
           call eos_soundspeed_mhd(cy(i,j,k), as, ca, cad)

           cad = q(i,j,k,QMAGZ)**2/q(i,j,k,QRHO)
           call eos_soundspeed_mhd(cz(i,j,k), as, ca, cad)

        enddo
     enddo
  enddo

  ! Load advected quatities, c, into q, assuming they arrived in uin as rho.c
  do ipassive = 1, npassive
     n = upass_map(ipassive)
     nq2 = qpass_map(ipassive)
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              q(i,j,k,nq2) = uin(i,j,k,n)/q(i,j,k,QRHO)
           enddo
        enddo
     enddo
  enddo

end subroutine ctoprim_mhd

subroutine srctoprim_mhd(lo, hi, &
                         q, q_lo, q_hi, &
                         src, src_lo, src_hi, &
                         srcQ, srcq_lo, srcq_hi) bind(C, name="srctoprim_mhd")

  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec, naux
  !      use eos_params_module
  use eos_module
  use eos_type_module, only : eos_t, eos_input_rt
  !      use flatten_module
  use amrex_constants_module
  use amrex_error_module
  use meth_params_module, only : URHO, UMX, UMY, UMZ, NVAR,&
       UEDEN, UEINT, UFA, UFS, UTEMP, &
       NQ, NSRC ,NQSRC, QRHO, QU, QV, QW, QC, NQAUX,&
       QREINT, QPRES, QFA, QFS, QTEMP, &
       QMAGX,  QMAGY, QMAGZ, &
       nadv, small_dens, small_pres, &
       npassive, upass_map, qpass_map

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: src_lo(3), src_hi(3)
  integer, intent(in) :: srcq_lo(3), srcq_hi(3)

  real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
  real(rt), intent(in) :: src(src_lo(1):src_hi(1), src_lo(2):src_hi(2), src_lo(3):src_hi(3), NSRC)
  real(rt), intent(inout) :: srcQ(srcq_lo(1):srcq_hi(1), srcq_lo(2):srcq_hi(2), srcq_lo(3):srcq_hi(3), NQSRC)

  real(rt) :: dpdr, dpde

  integer :: i, j, k
  integer :: n, nq2, ipassive

  real(rt) :: rhoInv

  type(eos_t) :: eos_state


  ! NOTE - WE ASSUME HERE THAT src(i,j,k,URHO) = 0. --
  !        IF NOT THEN THE FORMULAE BELOW ARE INCOMPLETE.

  ! compute srcQ terms
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           rhoInv = ONE/q(i,j,k,QRHO)

           eos_state % rho = q(i,j,k,QRHO)
           eos_state % T   = q(i,j,k,UTEMP)
           eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)

           call eos(eos_input_rt, eos_state)

           dpdr = eos_state % dpdr_e
           dpde = eos_state % dpde

           srcQ(i,j,k,QRHO  ) = src(i,j,k,URHO)
           srcQ(i,j,k,QU    ) = (src(i,j,k,UMX) - q(i,j,k,QU) * srcQ(i,j,k,QRHO)) * rhoInv
           srcQ(i,j,k,QV    ) = (src(i,j,k,UMY) - q(i,j,k,QV) * srcQ(i,j,k,QRHO)) * rhoInv
           srcQ(i,j,k,QW    ) = (src(i,j,k,UMZ) - q(i,j,k,QW) * srcQ(i,j,k,QRHO)) * rhoInv
           srcQ(i,j,k,QREINT) = src(i,j,k,UEINT)
           srcQ(i,j,k,QPRES ) = dpde*(srcQ(i,j,k,QREINT) - &
                q(i,j,k,QREINT)*srcQ(i,j,k,QRHO)*rhoinv) * rhoInv + &
                dpdr * srcQ(i,j,k,QRHO)

           ! add ifdef PRIM_SPECIES_HAVE_SOURCES

           !            if (UFS .gt. 0) then
           !               do ispec = 1,nspec+naux
           !                  srcQ(i,j,k,QFS+ispec-1) = src(i,j,k,UFS+ispec-1)*rhoInv
           !               enddo
           !            end if ! UFS > 0

           !            do iadv = 1,nadv
           !               srcQ(i,j,k,QFA+iadv-1) = src(i,j,k,UFA+iadv-1)*rhoInv
           !            enddo

        enddo
     enddo
  enddo

end subroutine srctoprim_mhd


subroutine check_for_mhd_cfl_violation(lo, hi, &
                                       q, q_lo, q_hi, &
                                       cx, cx_lo, cx_hi, &
                                       cy, cy_lo, cy_hi, &
                                       cz, cz_lo, cz_hi, &
                                       courno, dx, dt) bind(C, name="check_for_mhd_cfl_violation")

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
  use amrex_error_module
  use meth_params_module, only : NQ, QRHO, QU, QV, QW, QC, NQAUX,&
       QREINT, QPRES, QFA, QFS, QTEMP, &
       QMAGX,  QMAGY, QMAGZ

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: cx_lo(3), cx_hi(3)
  integer, intent(in) :: cy_lo(3), cy_hi(3)
  integer, intent(in) :: cz_lo(3), cz_hi(3)

  real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
  real(rt), intent(in) :: cx(cx_lo(1):cx_hi(1), cx_lo(2):cx_hi(2), cx_lo(3):cx_hi(3))
  real(rt), intent(in) :: cy(cy_lo(1):cy_hi(1), cy_lo(2):cy_hi(2), cy_lo(3):cx_hi(3))
  real(rt), intent(in) :: cz(cz_lo(1):cz_hi(1), cz_lo(2):cz_hi(2), cz_lo(3):cx_hi(3))

  real(rt), intent(in) :: dx(AMREX_SPACEDIM)
  real(rt), intent(out) :: courno
  real(rt), intent(in), value :: dt

  integer :: i, j, k

  real(rt) :: courx, coury, courz, courmx, courmy, courmz
  real(rt) :: dtdx, dtdy, dtdz

  ! Compute running max of Courant number over grids
  courmx = courno
  courmy = courno
  courmz = courno

  dtdx = dt / dx(1)
  dtdy = dt / dx(2)
  dtdz = dt / dx(3)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           courx = ( cx(i,j,k)+abs(q(i,j,k,QU)) ) * dtdx
           coury = ( cy(i,j,k)+abs(q(i,j,k,QV)) ) * dtdy
           courz = ( cz(i,j,k)+abs(q(i,j,k,QW)) ) * dtdz

           courmx = max( courmx, courx )
           courmy = max( courmy, coury )
           courmz = max( courmz, courz )

           if (courx .gt. ONE) then
              !
              ! A critical region since we usually can't write from threads.
              !
              print *,'   '
              print *,'>>> ... (u+c) * a * dt / dx > 1 ', courx
              print *,'>>> ... at cell (i,j,k)   : ',i,j,k
              print *,'>>> ... u, c                ',q(i,j,k,QU), cx(i,j,k)
              print *,'>>> ... B                   ',q(i,j,k,QMAGX:QMAGZ)
              print *,'>>> ... Internal e          ',q(i,j,k,QREINT)
              print *,'>>> ... pressure            ',q(i,j,k,QPRES)
              print *,'>>> ... density             ',q(i,j,k,QRHO)
              call amrex_error("Error:: Castro_advection_mhd_3d.f90 :: CFL violation in x-dir in ctoprim")
           end if

           if (coury .gt. ONE) then
              !
              ! A critical region since we usually can't write from threads.
              !
              print *,'   '
              print *,'>>> ... (v+c) * a * dt / dx > 1 ', coury
              print *,'>>> ... at cell (i,j,k)   : ',i,j,k
              print *,'>>> ... v, c                ',q(i,j,k,QV), cy(i,j,k)
              print *,'>>> ... B                   ',q(i,j,k,QMAGX:QMAGZ)
              print *,'>>> ... pressure            ',q(i,j,k,QPRES)
              print *,'>>> ... density             ',q(i,j,k,QRHO)
              call amrex_error("Error:: Castro_advection_mhd_3d.f90 :: CFL violation in y-dir in ctoprim")
           end if

           if (courz .gt. ONE) then
              !
              ! A critical region since we usually can't write from threads.
              !
              print *,'   '
              print *,'>>> ... (w+c) * a * dt / dx > 1 ', courz
              print *,'>>> ... at cell (i,j,k)   : ',i,j,k
              print *,'>>> ... w, c                ',q(i,j,k,QW), cz(i,j,k)
              print *,'>>> ... B                   ',q(i,j,k,QMAGX:QMAGZ)
              print *,'>>> ... pressure            ',q(i,j,k,QPRES)
              print *,'>>> ... density             ',q(i,j,k,QRHO)
              call amrex_error("Error:: Castro_advection_mhd_3d.f90 :: CFL violation in z-dir in ctoprim")
           end if

        enddo
     enddo
  enddo
  courno = max( courmx, courmy, courmz )

end subroutine check_for_mhd_cfl_violation


! :::
! ::: ========================== Conservative Update ===============================================================
! :::

subroutine consup(lo, hi, &
                  uin, uin_lo, uin_hi, &
                  uout, uout_lo, uout_hi, &
                  bcc, bcc_lo, bcc_hi, &
                  fluxx, flux1_lo, flux1_hi, &
                  fluxy, flux2_lo, flux2_hi, &
                  fluxz, flux3_lo, flux3_hi, &
                  dx, dt) bind(C, name="consup")

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module, only : UMX,UMY,UMZ, NVAR, URHO, UEDEN, UEINT, UFS
  use network, only: nspec

  implicit none

  integer, intent(in)  :: uin_lo(3), uin_hi(3)
  integer, intent(in)  :: uout_lo(3), uout_hi(3)
  integer, intent(in)  :: bcc_lo(3), bcc_hi(3)
  integer, intent(in)  :: flux1_lo(3), flux1_hi(3)
  integer, intent(in)  :: flux2_lo(3), flux2_hi(3)
  integer, intent(in)  :: flux3_lo(3), flux3_hi(3)
  integer, intent(in)   :: lo(3), hi(3)

  real(rt), intent(in)  :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
  real(rt), intent(inout)  :: bcc(bcc_lo(1):bcc_hi(1), bcc_lo(2):bcc_hi(2), bcc_lo(3):bcc_hi(3), 3)
  real(rt), intent(in)  :: fluxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR+3)
  real(rt), intent(in)  :: fluxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR+3)
  real(rt), intent(in)  :: fluxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR+3)
  real(rt), intent(in)  :: dx(3)
  real(rt), intent(in), value :: dt
  real(rt), intent(out) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3),NVAR)

  integer :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           uout(i,j,k,URHO) = uout(i,j,k,URHO) &
                - 1.0/dx(1)*(fluxx(i+1,j,k,URHO) - fluxx(i,j,k,URHO)) &
                - 1.0/dx(2)*(fluxy(i,j+1,k,URHO) - fluxy(i,j,k,URHO)) &
                - 1.0/dx(3)*(fluxz(i,j,k+1,URHO) - fluxz(i,j,k,URHO))
           uout(i,j,k,UMX) = uout(i,j,k,UMX) &
                - 1.0/dx(1)*(fluxx(i+1,j,k,UMX) - fluxx(i,j,k,UMX)) &
                - 1.0/dx(2)*(fluxy(i,j+1,k,UMX) - fluxy(i,j,k,UMX)) &
                - 1.0/dx(3)*(fluxz(i,j,k+1,UMX) - fluxz(i,j,k,UMX))
           uout(i,j,k,UMY) = uout(i,j,k,UMY) &
                - 1.0/dx(1)*(fluxx(i+1,j,k,UMY) - fluxx(i,j,k,UMY)) &
                - 1.0/dx(2)*(fluxy(i,j+1,k,UMY) - fluxy(i,j,k,UMY)) &
                - 1.0/dx(3)*(fluxz(i,j,k+1,UMY) - fluxz(i,j,k,UMY))
           uout(i,j,k,UMZ) = uout(i,j,k,UMZ) &
                - 1.0/dx(1)*(fluxx(i+1,j,k,UMZ) - fluxx(i,j,k,UMZ)) &
                - 1.0/dx(2)*(fluxy(i,j+1,k,UMZ) - fluxy(i,j,k,UMZ)) &
                - 1.0/dx(3)*(fluxz(i,j,k+1,UMZ) - fluxz(i,j,k,UMZ))
           uout(i,j,k,UEDEN) = uout(i,j,k,UEDEN) &
                - 1.0/dx(1)*(fluxx(i+1,j,k,UEDEN) - fluxx(i,j,k,UEDEN)) &
                - 1.0/dx(2)*(fluxy(i,j+1,k,UEDEN) - fluxy(i,j,k,UEDEN)) &
                - 1.0/dx(3)*(fluxz(i,j,k+1,UEDEN) - fluxz(i,j,k,UEDEN))
           uout(i,j,k,UEINT) = uout(i,j,k,UEINT) &
                - 1.0/dx(1)*(fluxx(i+1,j,k,UEINT) - fluxx(i,j,k,UEINT)) &
                - 1.0/dx(2)*(fluxy(i,j+1,k,UEINT) - fluxy(i,j,k,UEINT)) &
                - 1.0/dx(3)*(fluxz(i,j,k+1,UEINT) - fluxz(i,j,k,UEINT))
           uout(i,j,k,UFS:UFS+nspec-1) = uout(i,j,k,UFS:UFS+nspec-1) &
                - 1.0/dx(1)*(fluxx(i+1,j,k,UFS:UFS+nspec-1) - fluxx(i,j,k,UFS:UFS+nspec-1)) &
                - 1.0/dx(2)*(fluxy(i,j+1,k,UFS:UFS+nspec-1) - fluxy(i,j,k,UFS:UFS+nspec-1)) &
                - 1.0/dx(3)*(fluxz(i,j,k+1,UFS:UFS+nspec-1) - fluxz(i,j,k,UFS:UFS+nspec-1))

           bcc(i,j,k,:) = bcc(i,j,k,:) - dt/dx(1)*(fluxx(i+1, j, k, NVAR+1:NVAR+3)- fluxx(i,j,k, NVAR+1:NVAR+3)) &
                                       - dt/dx(2)*(fluxy(i, j+1, k, NVAR+1:NVAR+3)- fluxy(i,j,k, NVAR+1:NVAR+3)) &
                                       - dt/dx(3)*(fluxz(i, j, k+1, NVAR+1:NVAR+3)- fluxz(i,j,k, NVAR+1:NVAR+3))
        enddo
     enddo
  enddo

end subroutine consup

! :::
! ::: ========================== Magnetic Update ===============================================================
! :::

subroutine magup(lo, hi, &
                 bxin, bxin_lo, bxin_hi, &
                 byin, byin_lo, byin_hi, &
                 bzin, bzin_lo, bzin_hi, &
                 bxout, bxout_lo, bxout_hi, &
                 byout, byout_lo, byout_hi, &
                 bzout, bzout_lo, bzout_hi, &
                 Ex, Ex_lo, Ex_hi, &
                 Ey, Ey_lo, Ey_hi, &
                 Ez, Ez_lo, Ez_hi, &
                 dx, dt) bind(C, name="magup")

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module!, only : QVAR, NVAR, UEINT

  implicit none

  integer, intent(in)   :: bxin_lo(3), bxin_hi(3)
  integer, intent(in)   :: byin_lo(3), byin_hi(3)
  integer, intent(in)   :: bzin_lo(3), bzin_hi(3)
  integer, intent(in)   :: bxout_lo(3), bxout_hi(3)
  integer, intent(in)   :: byout_lo(3), byout_hi(3)
  integer, intent(in)   :: bzout_lo(3), bzout_hi(3)
  integer, intent(in)   :: Ex_lo(3), Ex_hi(3)
  integer, intent(in)   :: Ey_lo(3), Ey_hi(3)
  integer, intent(in)   :: Ez_lo(3), Ez_hi(3)
  integer, intent(in)   :: lo(3), hi(3)

  real(rt), intent(in)  :: bxin(bxin_lo(1):bxin_hi(1), bxin_lo(2):bxin_hi(2), bxin_lo(3):bxin_hi(3))
  real(rt), intent(in)  :: byin(byin_lo(1):byin_hi(1), byin_lo(2):byin_hi(2), byin_lo(3):byin_hi(3))
  real(rt), intent(in)  :: bzin(bzin_lo(1):bzin_hi(1), bzin_lo(2):bzin_hi(2), bzin_lo(3):bzin_hi(3))

  real(rt), intent(in) ::  Ex(Ex_lo(1):Ex_hi(1), Ex_lo(2):Ex_hi(2), Ex_lo(3):Ex_hi(3))
  real(rt), intent(in) ::  Ey(Ey_lo(1):Ey_hi(1), Ey_lo(2):Ey_hi(2), Ey_lo(3):Ey_hi(3))
  real(rt), intent(in) ::  Ez(Ez_lo(1):Ez_hi(1), Ez_lo(2):Ez_hi(2), Ez_lo(3):Ez_hi(3))

  real(rt), intent(in)  :: dx(3)
  real(rt), intent(in), value :: dt


  real(rt), intent(out) :: bxout(bxout_lo(1):bxout_hi(1), bxout_lo(2):bxout_hi(2), bxout_lo(3):bxout_hi(3))
  real(rt), intent(out) :: byout(byout_lo(1):byout_hi(1), byout_lo(2):byout_hi(2), byout_lo(3):byout_hi(3))
  real(rt), intent(out) :: bzout(bzout_lo(1):bzout_hi(1), bzout_lo(2):bzout_hi(2), bzout_lo(3):bzout_hi(3))

  integer                                 :: i, j, k

  !***** TO DO ***** SOURCES
  !-------------------------------- bx --------------------------------------------------
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           bxout(i,j,k) = bxin(i,j,k) + dt/dx(1)*( (Ey(i,j,k+1) - Ey(i,j,k)) - (Ez(i,j+1,k) - Ez(i,j,k)) )

        enddo
     enddo
  enddo

  !------------------------------- by --------------------------------------------------
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
           byout(i,j,k) = byin(i,j,k) + dt/dx(2)*( (Ez(i+1,j,k) - Ez(i,j,k)) - (Ex(i,j,k+1) - Ex(i,j,k)) )
        enddo
     enddo
  enddo
 !------------------------------- bz --------------------------------------------------
  do k = lo(3), hi(3)+1
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           bzout(i,j,k) = bzin(i,j,k) + dt/dx(3)*( (Ex(i,j+1,k) - Ex(i,j,k)) - (Ey(i+1,j,k) - Ey(i,j,k)) )
        enddo
     enddo
  enddo

 end subroutine magup

end module mhd_advection_module
