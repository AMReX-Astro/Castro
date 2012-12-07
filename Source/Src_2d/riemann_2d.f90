module riemann_module

  implicit none

  private

  public riemannus

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine riemannus(ql, qr, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                           gamcl, gamcr, cav, smallc, gd_l1, gd_l2, gd_h1, gd_h2, &
                           uflx, uflx_l1, uflx_l2, uflx_h1, uflx_h2, &
                           pgdnv, pgd_l1, pgd_l2, pgd_h1, pgd_h2, &
                           ugdnv, ugd_l1, ugd_l2, ugd_h1, ugd_h2, &
                           idir, ilo1, ihi1, ilo2, ihi2)

      use network, only : nspec, naux
      use prob_params_module, only : physbc_lo,physbc_hi,Symmetry
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, QFA, QFS, QFX, &
                                     URHO, UMX, UMY, UEDEN, UEINT, UFA, UFS, UFX, nadv, &
                                     small_dens, small_pres

      implicit none

      double precision, parameter:: small = 1.d-8

      integer qpd_l1, qpd_l2, qpd_h1, qpd_h2
      integer gd_l1, gd_l2, gd_h1, gd_h2
      integer uflx_l1, uflx_l2, uflx_h1, uflx_h2
      integer pgd_l1, pgd_l2, pgd_h1, pgd_h2
      integer ugd_l1, ugd_l2, ugd_h1, ugd_h2
      integer idir, ilo1, ihi1, ilo2, ihi2
      integer ilo,ihi,jlo,jhi

      double precision ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision cav(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision smallc(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
      double precision pgdnv(pgd_l1:pgd_h1,pgd_l2:pgd_h2)
      double precision ugdnv(ugd_l1:ugd_h1,ugd_l2:ugd_h2)

      integer iadv, ispec, iaux, n, nq
      integer i, j

      double precision rgd, vgd, regd, ustar
      double precision rl, ul, vl, pl, rel
      double precision rr, ur, vr, pr, rer
      double precision wl, wr, rhoetot, scr
      double precision rstar, cstar, estar, pstar
      double precision ro, uo, po, reo, co, gamco, entho
      double precision sgnm, spin, spout, ushock, frac
      double precision wsmall, csmall,qavg

!************************************************************
!  set min/max based on normal direction
      if(idir.eq.1) then
         ilo = ilo1
         ihi = ihi1 + 1
         jlo = ilo2
         jhi = ihi2
      else
         ilo = ilo1
         ihi = ihi1
         jlo = ilo2
         jhi = ihi2+1
      endif

!     Solve Riemann Problem
!     NOTE: The calling routine will order velocity unknowns so that
!     for the purposes of this routine, the normal component is always
!     loaded in the QU slot.
      do j = jlo, jhi
         do i = ilo, ihi

            rl = ql(i,j,QRHO)

!  pick left velocities based on direction
            if(idir.eq.1) then
               ul = ql(i,j,QU)
               vl = ql(i,j,QV)
            else
               ul = ql(i,j,QV)
               vl = ql(i,j,QU)
            endif

            pl = ql(i,j,QPRES)
            rel = ql(i,j,QREINT)

            rr = qr(i,j,QRHO)

!  pick right velocities based on direction
            if(idir.eq.1) then
               ur = qr(i,j,QU)
               vr = qr(i,j,QV)
            else
               ur = qr(i,j,QV)
               vr = qr(i,j,QU)
            endif

            pr = qr(i,j,QPRES)
            rer = qr(i,j,QREINT)

            csmall = smallc(i,j)
            wsmall = small_dens*csmall
            wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
            wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

            pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
            pstar = max(pstar,small_pres)
            ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)

            if (ustar .gt. 0.d0) then
               ro = rl
               uo = ul
               po = pl
               reo = rel
               gamco = gamcl(i,j)
            else if (ustar .lt. 0.d0) then
               ro = rr
               uo = ur
               po = pr
               reo = rer
               gamco = gamcr(i,j)
            else
               ro = 0.5d0*(rl+rr)
               uo = 0.5d0*(ul+ur)
               po = 0.5d0*(pl+pr)
               reo = 0.5d0*(rel+rer)
               gamco = 0.5d0*(gamcl(i,j)+gamcr(i,j))               
            endif
            ro = max(small_dens,ro)

            co = sqrt(abs(gamco*po/ro))
            co = max(csmall,co)
            entho = (reo/ro + po/ro)/co**2
            rstar = ro + (pstar - po)/co**2
            rstar = max(small_dens,rstar)
            estar = reo + (pstar - po)*entho
            cstar = sqrt(abs(gamco*pstar/rstar))
            cstar = max(cstar,csmall)

            sgnm = sign(1.d0,ustar)
            spout = co - sgnm*uo
            spin = cstar - sgnm*ustar
            ushock = 0.5d0*(spin + spout)
            if (pstar-po .ge. 0.d0) then
               spin = ushock
               spout = ushock
            endif
            if (spout-spin .eq. 0.d0) then
               scr = small*cav(i,j)
            else
               scr = spout-spin
            endif
            frac = (1.d0 + (spout + spin)/scr)*0.5d0
            frac = max(0.d0,min(1.d0,frac))

            if (ustar .gt. 0.d0) then
               vgd = vl
            else if (ustar .lt. 0.d0) then
               vgd = vr
            else
               vgd = 0.5d0*(vl+vr)
            endif
            rgd = frac*rstar + (1.d0 - frac)*ro

            ugdnv(i,j) = frac*ustar + (1.d0 - frac)*uo
            pgdnv(i,j) = frac*pstar + (1.d0 - frac)*po

            regd = frac*estar + (1.d0 - frac)*reo
            if (spout .lt. 0.d0) then
               rgd = ro
               ugdnv(i,j) = uo
               pgdnv(i,j) = po
               regd = reo
            endif
            if (spin .ge. 0.d0) then
               rgd = rstar
               ugdnv(i,j) = ustar
               pgdnv(i,j) = pstar
               regd = estar
            endif

            ! Enforce that fluxes through a symmetry plane are hard zero.
            if (i.eq.0 .and. physbc_lo(1) .eq. Symmetry .and. idir .eq. 1) ugdnv(i,j) = 0.d0
            if (j.eq.0 .and. physbc_lo(2) .eq. Symmetry .and. idir .eq. 2) ugdnv(i,j) = 0.d0

            ! Compute fluxes, order as conserved state (not q)
            uflx(i,j,URHO) = rgd*ugdnv(i,j)
            if(idir.eq.1) then
               uflx(i,j,UMX) = uflx(i,j,URHO)*ugdnv(i,j)
               uflx(i,j,UMY) = uflx(i,j,URHO)*vgd
            else
               uflx(i,j,UMX) = uflx(i,j,URHO)*vgd
               uflx(i,j,UMY) = uflx(i,j,URHO)*ugdnv(i,j)
            endif

            rhoetot = regd + 0.5d0*rgd*(ugdnv(i,j)**2 + vgd**2)
            uflx(i,j,UEDEN) = ugdnv(i,j)*(rhoetot + pgdnv(i,j))
            uflx(i,j,UEINT) = ugdnv(i,j)*regd

            do iadv = 1, nadv
               n = UFA + iadv - 1
               nq = QFA + iadv - 1
               if (ustar .gt. 0.d0) then
                  uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
               else if (ustar .lt. 0.d0) then
                  uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
               else 
                  qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                  uflx(i,j,n) = uflx(i,j,URHO)*qavg
               endif
            enddo

            do ispec = 1, nspec
               n  = UFS + ispec - 1
               nq = QFS + ispec - 1
               if (ustar .gt. 0.d0) then
                  uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
               else if (ustar .lt. 0.d0) then
                  uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
               else 
                  qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                  uflx(i,j,n) = uflx(i,j,URHO)*qavg
               endif
            enddo

            do iaux = 1, naux
               n  = UFX + iaux - 1
               nq = QFX + iaux - 1
               if (ustar .gt. 0.d0) then
                  uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
               else if (ustar .lt. 0.d0) then
                  uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
               else 
                  qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                  uflx(i,j,n) = uflx(i,j,URHO)*qavg
               endif
            enddo

         enddo
      enddo
      end subroutine riemannus

end module riemann_module
