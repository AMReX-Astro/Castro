module riemann_module

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NQ, NQAUX, NVAR, QRHO, QU, QV, QW, &
                                QPRES, QGAME, QREINT, QFS, &
                                QFX, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, &
                                UFS, UFX, &
                                NGDNV, GDRHO, GDPRES, GDGAME, &
                                QC, QGAMC, &
                                small_dens, small_temp, &
                                npassive, upass_map, qpass_map

  implicit none

  real(rt), parameter :: smallu = 1.e-12_rt

contains

  AMREX_DEVICE subroutine cmpflx(lo, hi, domlo, domhi, idir, &
                                 qm, qm_lo, qm_hi, &
                                 qp, qp_lo, qp_hi, &
                                 qint, qe_lo, qe_hi, &
                                 flx, flx_lo, flx_hi, &
                                 qaux, qa_lo, qa_hi)

    use network, only: nspec, naux
    use amrex_fort_module, only: rt => amrex_real
    use bl_constants_module, only: ZERO, HALF, ONE
    use prob_params_module, only: physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall

    integer,  intent(in   ) :: qm_lo(3), qm_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)
    integer,  intent(in   ) :: qe_lo(3), qe_hi(3)
    integer,  intent(in   ) :: flx_lo(3), flx_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ), value :: idir

    real(rt), intent(in   ) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ,3)
    real(rt), intent(in   ) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ,3)
    real(rt), intent(inout) :: qint(qe_lo(1):qe_hi(1),qe_lo(2):qe_hi(2),qe_lo(3):qe_hi(3),NGDNV)
    real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)
    real(rt), intent(in   ) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    ! local variables

    integer  :: i, j, k
    integer  :: is_shock
    real(rt) :: cl, cr

    integer :: n, nqp, ipassive

    real(rt) :: regdnv
    real(rt) :: rl, ul, v1l, v2l, pl, rel
    real(rt) :: rr, ur, v1r, v2r, pr, rer
    real(rt) :: wl, wr, rhoetot, scr
    real(rt) :: rstar, cstar, estar, pstar, ustar
    real(rt) :: ro, uo, po, reo, co, gamco, entho, drho
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall, smallc, gamcm, gamcp, cavg, qavg

    real(rt) :: u_adv

    integer :: iu, iv1, iv2, im1, im2, im3
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    real(rt) :: bnd_fac_x, bnd_fac_y, bnd_fac_z
    real(rt) :: wwinv, roinv, co2inv

    real(rt), parameter :: small = 1.e-8_rt
    real(rt), parameter :: small_pres = 1.e-200_rt

    if (idir .eq. 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
       im1 = UMX
       im2 = UMY
       im3 = UMZ

    else if (idir .eq. 2) then
       iu = QV
       iv1 = QU
       iv2 = QW
       im1 = UMY
       im2 = UMX
       im3 = UMZ

    else
       iu = QW
       iv1 = QU
       iv2 = QV
       im1 = UMZ
       im2 = UMX
       im3 = UMY
    end if

    special_bnd_lo = (physbc_lo(idir) .eq. Symmetry &
         .or.         physbc_lo(idir) .eq. SlipWall &
         .or.         physbc_lo(idir) .eq. NoSlipWall)
    special_bnd_hi = (physbc_hi(idir) .eq. Symmetry &
         .or.         physbc_hi(idir) .eq. SlipWall &
         .or.         physbc_hi(idir) .eq. NoSlipWall)

    if (idir .eq. 1) then
       special_bnd_lo_x = special_bnd_lo
       special_bnd_hi_x = special_bnd_hi
    else
       special_bnd_lo_x = .false.
       special_bnd_hi_x = .false.
    end if

    do k = lo(3), hi(3)

       bnd_fac_z = ONE
       if (idir.eq.3) then
          if ( k .eq. domlo(3)   .and. special_bnd_lo .or. &
               k .eq. domhi(3)+1 .and. special_bnd_hi ) then
             bnd_fac_z = ZERO
          end if
       end if

       do j = lo(2), hi(2)

          bnd_fac_y = ONE
          if (idir .eq. 2) then
             if ( j .eq. domlo(2)   .and. special_bnd_lo .or. &
                  j .eq. domhi(2)+1 .and. special_bnd_hi ) then
                bnd_fac_y = ZERO
             end if
          end if

          do i = lo(1), hi(1)

             if (idir == 1) then
                smallc = max( small, max( small*qaux(i,j,k,QC), small * qaux(i-1,j,k,QC)) )
                cavg   = HALF*( qaux(i,j,k,QC) + qaux(i-1,j,k,QC) )
                gamcm  = qaux(i-1,j,k,QGAMC)
                gamcp  = qaux(i,j,k,QGAMC)
             elseif (idir == 2) then
                smallc = max( small, max( small*qaux(i,j,k,QC), small * qaux(i,j-1,k,QC)) )
                cavg   = HALF*( qaux(i,j,k,QC) + qaux(i,j-1,k,QC) )
                gamcm  = qaux(i,j-1,k,QGAMC)
                gamcp  = qaux(i,j,k,QGAMC)
             else
                smallc = max( small, max( small*qaux(i,j,k,QC), small * qaux(i,j,k-1,QC)) )
                cavg   = HALF*( qaux(i,j,k,QC) + qaux(i,j,k-1,QC) )
                gamcm  = qaux(i,j,k-1,QGAMC)
                gamcp  = qaux(i,j,k,QGAMC)
             endif

             rl = max(qm(i,j,k,QRHO,idir), small_dens)

             ! pick left velocities based on direction
             ul  = qm(i,j,k,iu,idir)
             v1l = qm(i,j,k,iv1,idir)
             v2l = qm(i,j,k,iv2,idir)

             pl  = max(qm(i,j,k,QPRES ,idir), small_pres)
             rel =     qm(i,j,k,QREINT,idir)
             rr  = max(qp(i,j,k,QRHO,idir), small_dens)

             ! pick right velocities based on direction
             ur  = qp(i,j,k,iu,idir)
             v1r = qp(i,j,k,iv1,idir)
             v2r = qp(i,j,k,iv2,idir)

             pr  = max(qp(i,j,k,QPRES,idir), small_pres)
             rer =     qp(i,j,k,QREINT,idir)
             csmall = smallc
             wsmall = small_dens*csmall
             wl = max(wsmall,sqrt(abs(gamcm*pl*rl)))
             wr = max(wsmall,sqrt(abs(gamcp*pr*rr)))

             wwinv = ONE/(wl + wr)
             pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
             ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv

             pstar = max(pstar,small_pres)
             ! for symmetry preservation, if ustar is really small, then we
             ! set it to zero
             if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
                ustar = ZERO
             endif

             if (ustar > ZERO) then
                ro = rl
                uo = ul
                po = pl
                reo = rel
                gamco = gamcm
             else if (ustar < ZERO) then
                ro = rr
                uo = ur
                po = pr
                reo = rer
                gamco = gamcp
             else
                ro = HALF*(rl+rr)
                uo = HALF*(ul+ur)
                po = HALF*(pl+pr)
                reo = HALF*(rel+rer)
                gamco = HALF*(gamcm+gamcp)
             endif

             ro = max(small_dens,ro)

             roinv = ONE/ro

             co = sqrt(abs(gamco*po*roinv))
             co = max(csmall,co)
             co2inv = ONE/(co*co)

             drho = (pstar - po)*co2inv
             rstar = ro + drho
             rstar = max(small_dens,rstar)

             entho = (reo + po)*roinv*co2inv
             estar = reo + (pstar - po)*entho
             cstar = sqrt(abs(gamco*pstar/rstar))
             cstar = max(cstar,csmall)

             sgnm = sign(ONE,ustar)
             spout = co - sgnm*uo
             spin = cstar - sgnm*ustar
             ushock = HALF*(spin + spout)

             if (pstar-po > ZERO) then
                spin = ushock
                spout = ushock
             endif

             if (spout-spin == ZERO) then
                scr = small*cavg
             else
                scr = spout-spin
             endif

             frac = (ONE + (spout + spin)/scr)*HALF
             frac = max(ZERO,min(ONE,frac))

             if (ustar > ZERO) then
                qint(i,j,k,iv1) = v1l
                qint(i,j,k,iv2) = v2l
             else if (ustar < ZERO) then
                qint(i,j,k,iv1) = v1r
                qint(i,j,k,iv2) = v2r
             else
                qint(i,j,k,iv1) = HALF*(v1l+v1r)
                qint(i,j,k,iv2) = HALF*(v2l+v2r)
             endif
             qint(i,j,k,GDRHO) = frac*rstar + (ONE - frac)*ro
             qint(i,j,k,iu  ) = frac*ustar + (ONE - frac)*uo

             qint(i,j,k,GDPRES) = frac*pstar + (ONE - frac)*po
             regdnv = frac*estar + (ONE - frac)*reo
             if (spout < ZERO) then
                qint(i,j,k,GDRHO) = ro
                qint(i,j,k,iu  ) = uo
                qint(i,j,k,GDPRES) = po
                regdnv = reo
             endif

             if (spin >= ZERO) then
                qint(i,j,k,GDRHO) = rstar
                qint(i,j,k,iu  ) = ustar
                qint(i,j,k,GDPRES) = pstar
                regdnv = estar
             endif


             qint(i,j,k,GDGAME) = qint(i,j,k,GDPRES)/regdnv + ONE
             qint(i,j,k,GDPRES) = max(qint(i,j,k,GDPRES),small_pres)
             u_adv = qint(i,j,k,iu)

             ! Enforce that fluxes through a symmetry plane or wall are hard zero.
             if ( special_bnd_lo_x .and. i.eq.domlo(1) .or. &
                  special_bnd_hi_x .and. i.eq.domhi(1)+1 ) then
                bnd_fac_x = ZERO
             else
                bnd_fac_x = ONE
             end if
             u_adv = u_adv * bnd_fac_x*bnd_fac_y*bnd_fac_z


             ! Compute fluxes, order as conserved state (not q)
             flx(i,j,k,URHO) = qint(i,j,k,GDRHO)*u_adv

             flx(i,j,k,im1) = flx(i,j,k,URHO)*qint(i,j,k,iu ) + qint(i,j,k,GDPRES)
             flx(i,j,k,im2) = flx(i,j,k,URHO)*qint(i,j,k,iv1)
             flx(i,j,k,im3) = flx(i,j,k,URHO)*qint(i,j,k,iv2)

             rhoetot = regdnv + HALF*qint(i,j,k,GDRHO)*(qint(i,j,k,iu)**2 + qint(i,j,k,iv1)**2 + qint(i,j,k,iv2)**2)

             flx(i,j,k,UTEMP) = ZERO
             flx(i,j,k,UEDEN) = u_adv*(rhoetot + qint(i,j,k,GDPRES))
             flx(i,j,k,UEINT) = u_adv*regdnv

             ! passively advected quantities
             do ipassive = 1, npassive

                n  = upass_map(ipassive)
                nqp = qpass_map(ipassive)

                if (ustar > ZERO) then
                   flx(i,j,k,n) = flx(i,j,k,URHO)*qm(i,j,k,nqp,idir)

                else if (ustar < ZERO) then
                   flx(i,j,k,n) = flx(i,j,k,URHO)*qp(i,j,k,nqp,idir)

                else
                   qavg = HALF * (qm(i,j,k,nqp,idir) + qp(i,j,k,nqp,idir))
                   flx(i,j,k,n) = flx(i,j,k,URHO)*qavg
                endif

             end do

          end do
       end do
    end do

  end subroutine cmpflx

end module riemann_module
