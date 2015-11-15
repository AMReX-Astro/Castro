module riemann_util_module

  use bl_types
  use bl_constants_module

  implicit none

contains

  pure function bc_test(idir, i, j, domlo, domhi) result (f)
    
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall

    integer, intent(in) :: idir, i, j, domlo(*), domhi(*)
    integer :: f
    
    ! Enforce that fluxes through a symmetry plane or wall are hard zero.
    f = 1

    if (idir == 1) then
       if (i == domlo(1) .and. &
            (physbc_lo(1) == Symmetry .or. &
             physbc_lo(1) == SlipWall .or. &
             physbc_lo(1) == NoSlipWall) ) then
          f = 0
       endif

       if (i == domhi(1)+1 .and. &
            (physbc_hi(1) == Symmetry .or. &
             physbc_hi(1) == SlipWall .or. &
             physbc_hi(1) == NoSlipWall) ) then
          f = 0
       endif
    end if

    if (idir == 2) then
       if (j == domlo(2) .and. &
            (physbc_lo(2) == Symmetry .or. &
             physbc_lo(2) == SlipWall .or. &
             physbc_lo(2) == NoSlipWall) ) then
          f = 0
       endif

       if (j == domhi(2)+1 .and. &
            (physbc_hi(2) == Symmetry .or. &
             physbc_hi(2) == SlipWall .or. &
             physbc_hi(2) == NoSlipWall) ) then
          f = 0
       end if
    endif

  end function bc_test


  subroutine HLL(ql, qr, cl, cr, idir, ndim, f)

    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   npassive, upass_map, qpass_map

    double precision, intent(in) :: ql(QVAR), qr(QVAR), cl, cr
    double precision, intent(inout) :: f(NVAR)
    integer, intent(in) :: idir, ndim

    integer :: ivel, ivelt, iveltt, imom, imomt, imomtt
    double precision :: a1, a4, bd, bl, bm, bp, br
    double precision :: cavg, uavg
    double precision :: fl_tmp, fr_tmp
    double precision :: rhod, rhoEl, rhoEr, rhol_sqrt, rhor_sqrt
    integer :: n, nq

    integer :: ipassive

    double precision, parameter :: small = 1.d-10

    select case (idir)
    case (1)
       ivel = QU
       ivelt = QV
       iveltt = QW

       imom = UMX
       imomt = UMY
       imomtt = UMZ

    case (2)
       ivel = QV
       ivelt = QU
       iveltt = QW

       imom = UMY
       imomt = UMX
       imomtt = UMZ

    case (3)
       ivel = QW
       ivelt = QU
       iveltt = QV

       imom = UMZ
       imomt = UMX
       imomtt = UMY

    end select

    rhol_sqrt = sqrt(ql(QRHO))
    rhor_sqrt = sqrt(qr(QRHO))

    rhod = ONE/(rhol_sqrt + rhor_sqrt)



    ! compute the average sound speed. This uses an approximation from
    ! E88, eq. 5.6, 5.7 that assumes gamma falls between 1
    ! and 5/3
    cavg = sqrt( (rhol_sqrt*cl**2 + rhor_sqrt*cr**2)*rhod + &
         HALF*rhol_sqrt*rhor_sqrt*rhod**2*(qr(ivel) - ql(ivel))**2 )


    ! Roe eigenvalues (E91, eq. 5.3b)
    uavg = (rhol_sqrt*ql(ivel) + rhor_sqrt*qr(ivel))*rhod

    a1 = uavg - cavg
    a4 = uavg + cavg


    ! signal speeds (E91, eq. 4.5)
    bl = min(a1, ql(ivel) - cl)
    br = max(a4, qr(ivel) + cr)

    bm = min(ZERO, bl)
    bp = max(ZERO, br)

    bd = bp - bm

    if (abs(bd) < small*max(abs(bm),abs(bp))) return

    bd = ONE/bd


    ! compute the fluxes according to E91, eq. 4.4b -- note that the
    ! min/max above picks the correct flux if we are not in the star
    ! region

    ! density flux
    fl_tmp = ql(QRHO)*ql(ivel)
    fr_tmp = qr(QRHO)*qr(ivel)

    f(URHO) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO) - ql(QRHO))


    ! normal momentum flux.  Note for 1- and 2-d, we leave off the pressure
    ! term and handle that separately in the update, to accommodate different
    ! geometries
    if (ndim < 3) then
       fl_tmp = ql(QRHO)*ql(ivel)**2 
       fr_tmp = qr(QRHO)*qr(ivel)**2 
    else
       fl_tmp = ql(QRHO)*ql(ivel)**2 + ql(QPRES)
       fr_tmp = qr(QRHO)*qr(ivel)**2 + qr(QPRES)
    endif

    f(imom) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivel) - ql(QRHO)*ql(ivel))


    ! transverse momentum flux
    fl_tmp = ql(QRHO)*ql(ivel)*ql(ivelt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(ivelt)

    f(imomt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivelt) - ql(QRHO)*ql(ivelt))


    fl_tmp = ql(QRHO)*ql(ivel)*ql(iveltt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(iveltt)

    f(imomtt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(iveltt) - ql(QRHO)*ql(iveltt))


    ! total energy flux
    rhoEl = ql(QREINT) + HALF*ql(QRHO)*(ql(ivel)**2 + ql(ivelt)**2 + ql(iveltt)**2)
    fl_tmp = ql(ivel)*(rhoEl + ql(QPRES))

    rhoEr = qr(QREINT) + HALF*qr(QRHO)*(qr(ivel)**2 + qr(ivelt)**2 + qr(iveltt)**2)
    fr_tmp = qr(ivel)*(rhoEr + qr(QPRES))

    f(UEDEN) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(rhoEr - rhoEl)


    ! eint flux
    fl_tmp = ql(QREINT)*ql(ivel)
    fr_tmp = qr(QREINT)*qr(ivel)

    f(UEINT) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QREINT) - ql(QREINT))


    ! passively-advected scalar fluxes
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       fl_tmp = ql(QRHO)*ql(nq)*ql(ivel)
       fr_tmp = qr(QRHO)*qr(nq)*qr(ivel)

       f(n) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(nq) - ql(QRHO)*ql(nq))
    enddo

  end subroutine HLL

  
  pure subroutine cons_state(q, U)

    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, QREINT, &
         NVAR, URHO, UMX, UMY, UEDEN, UEINT, &
         npassive, upass_map, qpass_map

    real (kind=dp_t), intent(in)  :: q(QVAR)
    real (kind=dp_t), intent(out) :: U(NVAR)

    integer :: ipassive, n, nq

    U(URHO) = q(QRHO)
    U(UMX)  = q(QRHO)*q(QU)
    U(UMY)  = q(QRHO)*q(QV)

    U(UEDEN) = q(QREINT) + HALF*q(QRHO)*(q(QU)**2 + q(QV)**2 + q(QW)**2)
    U(UEINT) = q(QREINT)

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       U(n) = q(QRHO)*q(nq)
    enddo

  end subroutine cons_state


  pure subroutine HLLC_state(idir, S_k, S_c, q, U)

    use meth_params_module, only: QVAR, QRHO, QU, QV, QW, QREINT, QPRES, &
         NVAR, URHO, UMX, UMY, UEDEN, UEINT, &
         npassive, upass_map, qpass_map

    integer, intent(in) :: idir
    real (kind=dp_t), intent(in)  :: S_k, S_c
    real (kind=dp_t), intent(in)  :: q(QVAR)
    real (kind=dp_t), intent(out) :: U(NVAR)

    real (kind=dp_t) :: hllc_factor, u_k
    integer :: ipassive, n, nq

    if (idir == 1) then
       u_k = q(QU)
    else
       u_k = q(QV)
    endif

    hllc_factor = q(QRHO)*(S_k - u_k)/(S_k - S_c)
    U(URHO) = hllc_factor
    if (idir == 1) then
       U(UMX)  = hllc_factor*S_c
       U(UMY)  = hllc_factor*q(QV)
    else
       U(UMX)  = hllc_factor*q(QU)
       U(UMY)  = hllc_factor*S_c
    endif

    U(UEDEN) = hllc_factor*(q(QREINT)/q(QRHO) + HALF*(q(QU)**2 + q(QV)**2 + q(QW)**2) + &
                            (S_c - u_k)*(S_c + q(QPRES)/(q(QRHO)*(S_k - u_k))))
    U(UEINT) = hllc_factor*q(QREINT)/q(QRHO)

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       U(n) = hllc_factor*q(nq)
    enddo

  end subroutine HLLC_state

  
  pure subroutine compute_flux(idir, U, p, F)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UEDEN, UEINT, &
         npassive, upass_map

    integer, intent(in) :: idir
    real (kind=dp_t), intent(in) :: U(NVAR)
    real (kind=dp_t), intent(in) :: p
    real (kind=dp_t), intent(out) :: F(NVAR)

    integer :: ipassive, n
    real (kind=dp_t) :: u_flx

    if (idir == 1) then
       u_flx = U(UMX)/U(URHO)
    else
       u_flx = U(UMY)/U(URHO)
    endif

    F(URHO) = U(URHO)*u_flx

    F(UMX) = U(UMX)*u_flx
    F(UMY) = U(UMY)*u_flx

    F(UEINT) = U(UEINT)*u_flx
    F(UEDEN) = (U(UEDEN) + p)*u_flx

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       F(n) = U(n)*u_flx
    enddo

  end subroutine compute_flux

end module riemann_util_module
