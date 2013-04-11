! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transx1(qym,qymo,qyp,qypo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fx,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                         ugdnvx,pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                         gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
                         cdtdx,ilo,ihi,jlo,jhi,kc,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, QESGS, QFA, QFS, QFX, &
                                     URHO, UMX, UMY, UMZ, UEDEN, UESGS, UFA, UFS, UFX, &
                                     nadv, small_pres

      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
      integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
      integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
      integer ilo,ihi,jlo,jhi,kc,k3d

      double precision  qym(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qyp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qymo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qypo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      double precision ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
      double precision cdtdx

      integer i, j
      integer n, nq
      integer iadv, ispec, iaux

      double precision rrnew, rr
      double precision rrry, rrly
      double precision rury, ruly
      double precision rvry, rvly
      double precision rwry, rwly
      double precision ekenry, ekenly
      double precision rery, rely
      double precision rrnewry, rrnewly
      double precision runewry, runewly
      double precision rvnewry, rvnewly
      double precision rwnewry, rwnewly
      double precision renewry, renewly
      double precision pnewry, pnewly
      double precision rhoekenry, rhoekenly
      double precision compn, compu, compsn, comps
      double precision pgp, pgm, ugp, ugm, dup, pav, du

      ! NOTE: it is better *not* to protect against small density in this routine

      ! Treat K as a passively advected quantity
      if (UESGS .gt. -1) then
         n  = UESGS
         nq = QESGS
         do j = jlo, jhi
            do i = ilo, ihi
    
               compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))
    
               rr = qyp(i,j,kc,QRHO)
               rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
               compu = rr*qyp(i,j,kc,nq) - compn
               qypo(i,j,kc,nq) = compu/rrnew
 
               rr = qym(i,j+1,kc,QRHO)
               rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
               compu = rr*qym(i,j+1,kc,nq) - compn
               qymo(i,j+1,kc,nq) = compu/rrnew
 
            enddo
         enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,nq,i,j,compn,rr,rrnew,compu) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = UFA + iadv - 1
         nq = QFA + iadv - 1
         do j = jlo, jhi 
            do i = ilo, ihi 

               compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

               rr = qyp(i,j,kc,QRHO)
               rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
               compu = rr*qyp(i,j,kc,nq) - compn
               qypo(i,j,kc,nq) = compu/rrnew

               rr = qym(i,j+1,kc,QRHO)
               rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
               compu = rr*qym(i,j+1,kc,nq) - compn
               qymo(i,j+1,kc,nq) = compu/rrnew

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,nq,i,j,compsn,rr,rrnew,comps) IF(nspec.gt.1)
      do ispec = 1, nspec
         n  = UFS + ispec - 1
         nq = QFS + ispec - 1
         do j = jlo, jhi 
            do i = ilo, ihi 

               compsn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

               rr = qyp(i,j,kc,QRHO)
               rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
               comps = rr*qyp(i,j,kc,nq) - compsn
               qypo(i,j,kc,nq) = comps/rrnew

               rr = qym(i,j+1,kc,QRHO)
               rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
               comps = rr*qym(i,j+1,kc,nq) - compsn
               qymo(i,j+1,kc,nq) = comps/rrnew

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,nq,i,j,compsn,rr,rrnew,comps) IF(naux.gt.1)
      do iaux = 1, naux
         n  = UFX + iaux - 1
         nq = QFX + iaux - 1
         do j = jlo, jhi 
            do i = ilo, ihi 

               compsn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

               rr = qyp(i,j,kc,QRHO)
               rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
               comps = rr*qyp(i,j,kc,nq) - compsn
               qypo(i,j,kc,nq) = comps/rrnew

               rr = qym(i,j+1,kc,QRHO)
               rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
               comps = rr*qym(i,j+1,kc,nq) - compsn
               qymo(i,j+1,kc,nq) = comps/rrnew

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,pgp,pgm,ugp,ugm,rrry,rury,rvry,rwry,ekenry,rery,rrly,ruly,rvly,rwly,ekenly,rely) &
      !$OMP PRIVATE(rrnewry,runewry,rvnewry,rwnewry,renewry,rrnewly,runewly,rvnewly,rwnewly,renewly,dup,pav,du,pnewry) &
      !$OMP PRIVATE(pnewly,rhoekenry,rhoekenly)
      do j = jlo, jhi 
         do i = ilo, ihi 

            pgp = pgdnvx(i+1,j,kc)
            pgm = pgdnvx(i,j,kc)
            ugp = ugdnvx(i+1,j,kc)
            ugm = ugdnvx(i,j,kc)

            ! Convert to conservation form
            rrry = qyp(i,j,kc,QRHO)
            rury = rrry*qyp(i,j,kc,QU)
            rvry = rrry*qyp(i,j,kc,QV)
            rwry = rrry*qyp(i,j,kc,QW)
            ekenry = 0.5d0*rrry*(qyp(i,j,kc,QU)**2 + qyp(i,j,kc,QV)**2 + qyp(i,j,kc,QW)**2)
            rery = qyp(i,j,kc,QREINT) + ekenry

            rrly = qym(i,j+1,kc,QRHO)
            ruly = rrly*qym(i,j+1,kc,QU)
            rvly = rrly*qym(i,j+1,kc,QV)
            rwly = rrly*qym(i,j+1,kc,QW)
            ekenly = 0.5d0*rrly* &
                 (qym(i,j+1,kc,QU)**2 + qym(i,j+1,kc,QV)**2 + qym(i,j+1,kc,QW)**2)
            rely = qym(i,j+1,kc,QREINT) + ekenly

            ! Add transverse predictor
            rrnewry = rrry - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
            runewry = rury - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
            rvnewry = rvry - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
            rwnewry = rwry - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
            renewry = rery - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))

            rrnewly = rrly - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
            runewly = ruly - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
            rvnewly = rvly - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
            rwnewly = rwly - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
            renewly = rely - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))

            dup = pgp*ugp - pgm*ugm
            pav = 0.5d0*(pgp+pgm)
            du = ugp-ugm

            pnewry = qyp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-1.d0))
            pnewly = qym(i,j+1,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-1.d0))

            ! Convert back to non-conservation form
            if (j.ge.jlo+1) then
               qypo(i,j,kc,QRHO) = rrnewry
               qypo(i,j,kc,QU) = runewry/qypo(i,j,kc,QRHO)
               qypo(i,j,kc,QV) = rvnewry/qypo(i,j,kc,QRHO)
               qypo(i,j,kc,QW) = rwnewry/qypo(i,j,kc,QRHO)
               rhoekenry = 0.5d0*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,kc,QRHO)
               qypo(i,j,kc,QREINT) = renewry - rhoekenry
               qypo(i,j,kc,QPRES) = max(pnewry,small_pres)
            end if

            if (j.le.jhi-1) then
               qymo(i,j+1,kc,QRHO) = rrnewly
               qymo(i,j+1,kc,QU) = runewly/qymo(i,j+1,kc,QRHO)
               qymo(i,j+1,kc,QV) = rvnewly/qymo(i,j+1,kc,QRHO)
               qymo(i,j+1,kc,QW) = rwnewly/qymo(i,j+1,kc,QRHO)
               rhoekenly = 0.5d0*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,kc,QRHO)
               qymo(i,j+1,kc,QREINT) = renewly - rhoekenly
               qymo(i,j+1,kc,QPRES) = max(pnewly,small_pres)
            end if

         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine transx1

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transx2(qzm,qzmo,qzp,qzpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fx,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                         ugdnvx,pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                         gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
                         cdtdx,ilo,ihi,jlo,jhi,kc,km,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, QESGS, QFA, QFS, QFX, &
                                     URHO, UMX, UMY, UMZ, UEDEN, UESGS, UFA, UFS, UFX, &
                                     nadv, small_pres

      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
      integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
      integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
      integer ilo,ihi,jlo,jhi,kc,km,k3d

      double precision  qzm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qzp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qzmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qzpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      double precision ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
      double precision cdtdx

      integer i, j
      integer n, nq
      integer iadv, ispec, iaux

      double precision rrnew, rr
      double precision rrrz, rrlz
      double precision rurz, rulz
      double precision rvrz, rvlz
      double precision rwrz, rwlz
      double precision ekenrz, ekenlz
      double precision rerz, relz
      double precision rrnewrz, rrnewlz
      double precision runewrz, runewlz
      double precision rvnewrz, rvnewlz
      double precision rwnewrz, rwnewlz
      double precision renewrz, renewlz
      double precision pnewrz, pnewlz
      double precision rhoekenrz, rhoekenlz
      double precision compn, compu, compsn, comps
      double precision pgp, pgm, ugp, ugm, dup, pav, du

      ! Treat K as a passively advected quantity
      if (UESGS .gt. -1) then
         n  = UESGS
         nq = QESGS
         do j = jlo, jhi
            do i = ilo, ihi
 
                compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))
 
                rr = qzp(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                compu = rr*qzp(i,j,kc,nq) - compn
                qzpo(i,j,kc,nq) = compu/rrnew
 
                compn = cdtdx*(fx(i+1,j,km,n) - fx(i,j,km,n))
 
                rr = qzm(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
                compu = rr*qzm(i,j,kc,nq) - compn
                qzmo(i,j,kc,nq) = compu/rrnew
 
             enddo
          enddo
       endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,nq,i,j,compn,rr,rrnew,compu) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = UFA + iadv - 1
         nq = QFA + iadv - 1
         do j = jlo, jhi 
            do i = ilo, ihi 

                compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

                rr = qzp(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                compu = rr*qzp(i,j,kc,nq) - compn
                qzpo(i,j,kc,nq) = compu/rrnew

                compn = cdtdx*(fx(i+1,j,km,n) - fx(i,j,km,n))
 
                rr = qzm(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
                compu = rr*qzm(i,j,kc,nq) - compn
                qzmo(i,j,kc,nq) = compu/rrnew

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DO PRIVATE(ispec,n,nq,i,j,compsn,rr,rrnew,comps) IF(nspec.gt.1)
       do ispec = 1, nspec
          n  = UFS + ispec - 1
          nq = QFS + ispec - 1
          do j = jlo, jhi 
             do i = ilo, ihi 

                compsn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

                rr = qzp(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                comps = rr*qzp(i,j,kc,nq) - compsn
                qzpo(i,j,kc,nq) = comps/rrnew

                compsn = cdtdx*(fx(i+1,j,km,n) - fx(i,j,km,n))

                rr = qzm(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
                comps = rr*qzm(i,j,kc,nq) - compsn
                qzmo(i,j,kc,nq) = comps/rrnew

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DO PRIVATE(iaux,n,nq,i,j,compsn,rr,rrnew,comps) IF(naux.gt.1)
       do iaux = 1, naux
          n  = UFX + iaux - 1
          nq = QFX + iaux - 1
          do j = jlo, jhi 
             do i = ilo, ihi 

                compsn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

                rr = qzp(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                comps = rr*qzp(i,j,kc,nq) - compsn
                qzpo(i,j,kc,nq) = comps/rrnew

                compsn = cdtdx*(fx(i+1,j,km,n) - fx(i,j,km,n))

                rr = qzm(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
                comps = rr*qzm(i,j,kc,nq) - compsn
                qzmo(i,j,kc,nq) = comps/rrnew

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DO PRIVATE(i,j,pgp,pgm,ugp,ugm,rrrz,rurz,rvrz,rwrz,ekenrz,rerz,rrlz,rulz,rvlz,rwlz,ekenlz) &
       !$OMP PRIVATE(relz,rrnewrz,runewrz,rvnewrz,rwnewrz,renewrz,rrnewlz,runewlz,rvnewlz,rwnewlz,renewlz,dup,pav) &
       !$OMP PRIVATE(du,pnewrz,pnewlz,rhoekenrz,rhoekenlz)
       do j = jlo, jhi 
          do i = ilo, ihi 

             pgp = pgdnvx(i+1,j,kc)
             pgm = pgdnvx(i,j,kc)
             ugp = ugdnvx(i+1,j,kc)
             ugm = ugdnvx(i,j,kc)

             dup = pgp*ugp - pgm*ugm
             pav = 0.5d0*(pgp+pgm)
             du = ugp-ugm

             ! Convert to conservation form
             rrrz = qzp(i,j,kc,QRHO)
             rurz = rrrz*qzp(i,j,kc,QU)
             rvrz = rrrz*qzp(i,j,kc,QV)
             rwrz = rrrz*qzp(i,j,kc,QW)
             ekenrz = 0.5d0*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 + qzp(i,j,kc,QW)**2)
             rerz = qzp(i,j,kc,QREINT) + ekenrz

             ! Add transverse predictor
             rrnewrz = rrrz - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             runewrz = rurz - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
             rvnewrz = rvrz - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
             rwnewrz = rwrz - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
             renewrz = rerz - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))

             pnewrz = qzp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-1.d0))

             ! Convert back to non-conservation form
             qzpo(i,j,kc,QRHO) = rrnewrz
             qzpo(i,j,kc,QU) = runewrz/qzpo(i,j,kc,QRHO)
             qzpo(i,j,kc,QV) = rvnewrz/qzpo(i,j,kc,QRHO)
             qzpo(i,j,kc,QW) = rwnewrz/qzpo(i,j,kc,QRHO)
             rhoekenrz = 0.5d0*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,kc,QRHO)
             qzpo(i,j,kc,QREINT) = renewrz - rhoekenrz
             qzpo(i,j,kc,QPRES) = max(pnewrz,small_pres)

             pgp = pgdnvx(i+1,j,km)
             pgm = pgdnvx(i,j,km)
             ugp = ugdnvx(i+1,j,km)
             ugm = ugdnvx(i,j,km)

             dup = pgp*ugp - pgm*ugm
             pav = 0.5d0*(pgp+pgm)
             du = ugp-ugm

             rrlz = qzm(i,j,kc,QRHO)
             rulz = rrlz*qzm(i,j,kc,QU)
             rvlz = rrlz*qzm(i,j,kc,QV)
             rwlz = rrlz*qzm(i,j,kc,QW)
             ekenlz = 0.5d0*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 + qzm(i,j,kc,QW)**2)
             relz = qzm(i,j,kc,QREINT) + ekenlz

             ! Add transverse predictor
             rrnewlz = rrlz - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
             runewlz = rulz - cdtdx*(fx(i+1,j,km,UMX) - fx(i,j,km,UMX))
             rvnewlz = rvlz - cdtdx*(fx(i+1,j,km,UMY) - fx(i,j,km,UMY))
             rwnewlz = rwlz - cdtdx*(fx(i+1,j,km,UMZ) - fx(i,j,km,UMZ))
             renewlz = relz - cdtdx*(fx(i+1,j,km,UEDEN) - fx(i,j,km,UEDEN))

             pnewlz = qzm(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d-1)-1.d0))

             qzmo(i,j,kc,QRHO) = rrnewlz
             qzmo(i,j,kc,QU) = runewlz/qzmo(i,j,kc,QRHO)
             qzmo(i,j,kc,QV) = rvnewlz/qzmo(i,j,kc,QRHO)
             qzmo(i,j,kc,QW) = rwnewlz/qzmo(i,j,kc,QRHO)
             rhoekenlz = 0.5d0*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,kc,QRHO)
             qzmo(i,j,kc,QREINT) = renewlz - rhoekenlz
             qzmo(i,j,kc,QPRES) = max(pnewlz,small_pres)

          enddo
       enddo
       !$OMP END PARALLEL DO

      end subroutine transx2

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transy1(qxm,qxmo,qxp,qxpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fy,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                         ugdnvy,pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                         gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
                         cdtdy,ilo,ihi,jlo,jhi,kc,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, QESGS, QFA, QFS, QFX, &
                                     URHO, UMX, UMY, UMZ, UEDEN, UESGS, UFA, UFS, UFX, &
                                     nadv, small_pres
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
      integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
      integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
      integer ilo,ihi,jlo,jhi,kc,k3d

      double precision  qxm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qxp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qxmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qxpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      double precision ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
      double precision cdtdy

      integer i, j
      integer n, nq
      integer iadv, ispec, iaux

      double precision rrnew, rr
      double precision compn, compu, compsn, comps
      double precision rrrx, rrlx
      double precision rurx, rulx
      double precision rvrx, rvlx
      double precision rwrx, rwlx
      double precision ekenrx, ekenlx
      double precision rerx, relx
      double precision rrnewrx, rrnewlx
      double precision runewrx, runewlx
      double precision rvnewrx, rvnewlx
      double precision rwnewrx, rwnewlx
      double precision renewrx, renewlx
      double precision pnewrx, pnewlx
      double precision rhoekenrx, rhoekenlx
      double precision pgp, pgm, ugp, ugm, dup, pav, du

      ! Treat K as a passively advected quantity
      if (UESGS .gt. -1) then
         n  = UESGS
         nq = QESGS
         do j = jlo, jhi
            do i = ilo, ihi
 
               compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))
 
               rr = qxp(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               compu = rr*qxp(i,j,kc,nq) - compn
               qxpo(i,j,kc,nq) = compu/rrnew
 
               rr = qxm(i+1,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               compu = rr*qxm(i+1,j,kc,nq) - compn
               qxmo(i+1,j,kc,nq) = compu/rrnew
 
            enddo
         enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,nq,i,j,compn,rr,rrnew,compu) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = UFA + iadv - 1
         nq = QFA + iadv - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

               rr = qxp(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               compu = rr*qxp(i,j,kc,nq) - compn
               qxpo(i,j,kc,nq) = compu/rrnew

               rr = qxm(i+1,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               compu = rr*qxm(i+1,j,kc,nq) - compn
               qxmo(i+1,j,kc,nq) = compu/rrnew

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,nq,i,j,compsn,rr,rrnew,comps) IF(nspec.gt.1)
      do ispec = 1, nspec 
         n  = UFS + ispec - 1
         nq = QFS + ispec - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               compsn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

               rr = qxp(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               comps = rr*qxp(i,j,kc,nq) - compsn
               qxpo(i,j,kc,nq) = comps/rrnew

               rr = qxm(i+1,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               comps = rr*qxm(i+1,j,kc,nq) - compsn
               qxmo(i+1,j,kc,nq) = comps/rrnew

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,nq,i,j,compsn,rr,rrnew,comps) IF(naux.gt.1)
      do iaux = 1, naux 
         n  = UFX + iaux - 1
         nq = QFX + iaux - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               compsn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

               rr = qxp(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               comps = rr*qxp(i,j,kc,nq) - compsn
               qxpo(i,j,kc,nq) = comps/rrnew

               rr = qxm(i+1,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               comps = rr*qxm(i+1,j,kc,nq) - compsn
               qxmo(i+1,j,kc,nq) = comps/rrnew

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,pgp,pgm,ugp,ugm,rrrx,rurx,rvrx,rwrx,ekenrx,rerx,rrlx,rulx,rvlx,rwlx,ekenlx,relx) &
      !$OMP PRIVATE(rrnewrx,runewrx,rvnewrx,rwnewrx,renewrx,rrnewlx,runewlx,rvnewlx,rwnewlx,renewlx,dup,pav,du,pnewrx) &
      !$OMP PRIVATE(pnewlx,rhoekenrx,rhoekenlx)
      do j = jlo, jhi
         do i = ilo, ihi

            pgp = pgdnvy(i,j+1,kc)
            pgm = pgdnvy(i,j,kc)
            ugp = ugdnvy(i,j+1,kc)
            ugm = ugdnvy(i,j,kc)

            ! Convert to conservation form
            rrrx = qxp(i,j,kc,QRHO)
            rurx = rrrx*qxp(i,j,kc,QU)
            rvrx = rrrx*qxp(i,j,kc,QV)
            rwrx = rrrx*qxp(i,j,kc,QW)
            ekenrx = 0.5d0*rrrx*(qxp(i,j,kc,QU)**2 + qxp(i,j,kc,QV)**2 &
                 + qxp(i,j,kc,QW)**2)
            rerx = qxp(i,j,kc,QREINT) + ekenrx

            rrlx = qxm(i+1,j,kc,QRHO)
            rulx = rrlx*qxm(i+1,j,kc,QU)
            rvlx = rrlx*qxm(i+1,j,kc,QV)
            rwlx = rrlx*qxm(i+1,j,kc,QW)
            ekenlx = 0.5d0*rrlx*(qxm(i+1,j,kc,QU)**2 + qxm(i+1,j,kc,QV)**2 &
                 + qxm(i+1,j,kc,QW)**2)
            relx = qxm(i+1,j,kc,QREINT) + ekenlx

            ! Add transverse predictor
            rrnewrx = rrrx - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
            runewrx = rurx - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
            rvnewrx = rvrx - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
            rwnewrx = rwrx - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
            renewrx = rerx - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))

            rrnewlx = rrlx - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
            runewlx = rulx - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
            rvnewlx = rvlx - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
            rwnewlx = rwlx - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
            renewlx = relx - cdtdy*(fy(i,j+1,kc,UEDEN)- fy(i,j,kc,UEDEN))

            dup = pgp*ugp - pgm*ugm
            pav = 0.5d0*(pgp+pgm)
            du = ugp-ugm

            ! Convert back to non-conservation form
            if (i.ge.ilo+1) then
               qxpo(i,j,kc,QRHO) = rrnewrx
               qxpo(i,j,kc,QU) = runewrx/qxpo(i,j,kc,QRHO)
               qxpo(i,j,kc,QV) = rvnewrx/qxpo(i,j,kc,QRHO)
               qxpo(i,j,kc,QW) = rwnewrx/qxpo(i,j,kc,QRHO)
               rhoekenrx = 0.5d0*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,kc,QRHO)
               qxpo(i,j,kc,QREINT)= renewrx - rhoekenrx

               pnewrx = qxp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - 1.d0))
               qxpo(i,j,kc,QPRES) = max(pnewrx,small_pres)
            end if

            if (i.le.ihi-1) then
               qxmo(i+1,j,kc,QRHO) = rrnewlx
               qxmo(i+1,j,kc,QU) = runewlx/qxmo(i+1,j,kc,QRHO)
               qxmo(i+1,j,kc,QV) = rvnewlx/qxmo(i+1,j,kc,QRHO)
               qxmo(i+1,j,kc,QW) = rwnewlx/qxmo(i+1,j,kc,QRHO)
               rhoekenlx = 0.5d0*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,kc,QRHO)
               qxmo(i+1,j,kc,QREINT)= renewlx - rhoekenlx

               pnewlx = qxm(i+1,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - 1.d0))
               qxmo(i+1,j,kc,QPRES) = max(pnewlx,small_pres)
            end if

         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine transy1

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transy2(qzm,qzmo,qzp,qzpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fy,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                         ugdnvy,pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                         gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
                         cdtdy,ilo,ihi,jlo,jhi,kc,km,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, QESGS, QFA, QFS, QFX, &
                                     URHO, UMX, UMY, UMZ, UEDEN, UESGS, UFA, UFS, UFX, &
                                     nadv, small_pres
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
      integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
      integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
      integer ilo,ihi,jlo,jhi,kc,km,k3d

      double precision  qzm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qzp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qzmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qzpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      double precision ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
      double precision cdtdy

      integer i, j
      integer n, nq
      integer iadv, ispec, iaux

      double precision rrnew, rr
      double precision compn, compu, compsn, comps
      double precision rrrz, rrlz
      double precision rurz, rulz
      double precision rvrz, rvlz
      double precision rwrz, rwlz
      double precision ekenrz, ekenlz
      double precision rerz, relz
      double precision rrnewrz, rrnewlz
      double precision runewrz, runewlz
      double precision rvnewrz, rvnewlz
      double precision rwnewrz, rwnewlz
      double precision renewrz, renewlz
      double precision pnewrz, pnewlz
      double precision rhoekenrz, rhoekenlz
      double precision pgp, pgm, ugp, ugm, dup, pav, du

      ! Treat K as a passively advected quantity
      if (UESGS .gt. -1) then
         n  = UESGS
         nq = QESGS
         do j = jlo, jhi
            do i = ilo, ihi
 
               compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))
 
               rr = qzp(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               compu = rr*qzp(i,j,kc,nq) - compn
               qzpo(i,j,kc,nq) = compu/rrnew
 
               compn = cdtdy*(fy(i,j+1,km,n) - fy(i,j,km,n))
 
               rr = qzm(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
               compu = rr*qzm(i,j,kc,nq) - compn
               qzmo(i,j,kc,nq) = compu/rrnew
 
            enddo
         enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,nq,i,j,compn,rr,rrnew,compu) IF(nadv.gt.1)
      do iadv = 1, nadv
         n  = UFA + iadv - 1
         nq = QFA + iadv - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

               rr = qzp(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               compu = rr*qzp(i,j,kc,nq) - compn
               qzpo(i,j,kc,nq) = compu/rrnew

               compn = cdtdy*(fy(i,j+1,km,n) - fy(i,j,km,n))

               rr = qzm(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
               compu = rr*qzm(i,j,kc,nq) - compn
               qzmo(i,j,kc,nq) = compu/rrnew

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,nq,i,j,compsn,rr,rrnew,comps) IF(nspec.gt.1)
      do ispec = 1, nspec 
         n  = UFS + ispec - 1
         nq = QFS + ispec - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               compsn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

               rr = qzp(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               comps = rr*qzp(i,j,kc,nq) - compsn
               qzpo(i,j,kc,nq) = comps/rrnew

               compsn = cdtdy*(fy(i,j+1,km,n) - fy(i,j,km,n))

               rr = qzm(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
               comps = rr*qzm(i,j,kc,nq) - compsn
               qzmo(i,j,kc,nq) = comps/rrnew

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,nq,i,j,compsn,rr,rrnew,comps) IF(naux.gt.1)
      do iaux = 1, naux 
         n  = UFX + iaux - 1
         nq = QFX + iaux - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               compsn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

               rr = qzp(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
               comps = rr*qzp(i,j,kc,nq) - compsn
               qzpo(i,j,kc,nq) = comps/rrnew

               compsn = cdtdy*(fy(i,j+1,km,n) - fy(i,j,km,n))

               rr = qzm(i,j,kc,QRHO)
               rrnew = rr - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
               comps = rr*qzm(i,j,kc,nq) - compsn
               qzmo(i,j,kc,nq) = comps/rrnew

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,pgp,pgm,ugp,ugm,rrrz,rurz,rvrz,rwrz,ekenrz,rerz,rrlz,rulz,rvlz,rwlz,ekenlz,relz) &
      !$OMP PRIVATE(rrnewrz,runewrz,rvnewrz,rwnewrz,renewrz,rrnewlz,runewlz,rvnewlz,rwnewlz,renewlz,dup,pav,du,pnewrz) &
      !$OMP PRIVATE(pnewlz,rhoekenrz,rhoekenlz)
      do j = jlo, jhi
         do i = ilo, ihi

            pgp = pgdnvy(i,j+1,kc)
            pgm = pgdnvy(i,j,kc)
            ugp = ugdnvy(i,j+1,kc)
            ugm = ugdnvy(i,j,kc)

            ! Convert to conservation form
            rrrz = qzp(i,j,kc,QRHO)
            rurz = rrrz*qzp(i,j,kc,QU)
            rvrz = rrrz*qzp(i,j,kc,QV)
            rwrz = rrrz*qzp(i,j,kc,QW)
            ekenrz = 0.5d0*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 &
                 + qzp(i,j,kc,QW)**2)
            rerz = qzp(i,j,kc,QREINT) + ekenrz

            ! Add transverse predictor
            rrnewrz = rrrz - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
            runewrz = rurz - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
            rvnewrz = rvrz - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
            rwnewrz = rwrz - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
            renewrz = rerz - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))

            dup = pgp*ugp - pgm*ugm
            pav = 0.5d0*(pgp+pgm)
            du = ugp-ugm

            pnewrz = qzp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - 1.d0))

            ! Convert back to non-conservation form
            qzpo(i,j,kc,QRHO) = rrnewrz
            qzpo(i,j,kc,QU) = runewrz/qzpo(i,j,kc,QRHO)
            qzpo(i,j,kc,QV) = rvnewrz/qzpo(i,j,kc,QRHO)
            qzpo(i,j,kc,QW) = rwnewrz/qzpo(i,j,kc,QRHO)
            rhoekenrz = 0.5d0*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,kc,QRHO)
            qzpo(i,j,kc,QREINT)= renewrz - rhoekenrz
            qzpo(i,j,kc,QPRES) = max(pnewrz,small_pres)

            pgp = pgdnvy(i,j+1,km)
            pgm = pgdnvy(i,j,km)
            ugp = ugdnvy(i,j+1,km)
            ugm = ugdnvy(i,j,km)

            rrlz = qzm(i,j,kc,QRHO)
            rulz = rrlz*qzm(i,j,kc,QU)
            rvlz = rrlz*qzm(i,j,kc,QV)
            rwlz = rrlz*qzm(i,j,kc,QW)
            ekenlz = 0.5d0*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 &
                 + qzm(i,j,kc,QW)**2)
            relz = qzm(i,j,kc,QREINT) + ekenlz

            ! Add transverse predictor
            rrnewlz = rrlz - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
            runewlz = rulz - cdtdy*(fy(i,j+1,km,UMX) - fy(i,j,km,UMX))
            rvnewlz = rvlz - cdtdy*(fy(i,j+1,km,UMY) - fy(i,j,km,UMY))
            rwnewlz = rwlz - cdtdy*(fy(i,j+1,km,UMZ) - fy(i,j,km,UMZ))
            renewlz = relz - cdtdy*(fy(i,j+1,km,UEDEN)- fy(i,j,km,UEDEN))

            dup = pgp*ugp - pgm*ugm
            pav = 0.5d0*(pgp+pgm)
            du = ugp-ugm

            pnewlz = qzm(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))

            qzmo(i,j,kc,QRHO) = rrnewlz
            qzmo(i,j,kc,QU) = runewlz/qzmo(i,j,kc,QRHO)
            qzmo(i,j,kc,QV) = rvnewlz/qzmo(i,j,kc,QRHO)
            qzmo(i,j,kc,QW) = rwnewlz/qzmo(i,j,kc,QRHO)
            rhoekenlz = 0.5d0*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,kc,QRHO)
            qzmo(i,j,kc,QREINT)= renewlz - rhoekenlz
            qzmo(i,j,kc,QPRES) = max(pnewlz,small_pres)

         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine transy2

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transz(qxm,qxmo,qxp,qxpo, &
                        qym,qymo,qyp,qypo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        fz,fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
                        ugdnvz,pgdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                        gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
                        cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, QESGS, QFA, QFS, QFX,&
                                     URHO, UMX, UMY, UMZ, UEDEN, UESGS, UFA, UFS, UFX, &
                                     nadv, small_pres
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
      integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
      integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
      integer ilo,ihi,jlo,jhi,km,kc,k3d

      double precision  qxm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qxp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qym(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qyp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qxmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qxpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qymo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qypo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fz(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
      double precision ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      double precision pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
      double precision cdtdz

      integer n, nq
      integer iadv, ispec, iaux
      integer i, j

      double precision rrnew, rr
      double precision compn, compu, compsn, comps
      double precision rrrx, rrry, rrlx, rrly
      double precision rurx, rury, rulx, ruly
      double precision rvrx, rvry, rvlx, rvly
      double precision rwrx, rwry, rwlx, rwly
      double precision ekenrx, ekenry, ekenlx, ekenly
      double precision rerx, rery, relx, rely
      double precision rrnewrx, rrnewry, rrnewlx, rrnewly
      double precision runewrx, runewry, runewlx, runewly
      double precision rvnewrx, rvnewry, rvnewlx, rvnewly
      double precision rwnewrx, rwnewry, rwnewlx, rwnewly
      double precision renewrx, renewry, renewlx, renewly
      double precision pnewrx, pnewry, pnewlx, pnewly
      double precision rhoekenrx, rhoekenry, rhoekenlx, rhoekenly
      double precision pgp, pgm, ugp, ugm, dup, pav, du

      ! Treat K as a passively advected quantity
      if (UESGS .gt. -1) then
         n  = UESGS
         nq = QESGS
         do j = jlo, jhi
            do i = ilo, ihi
 
                 compn = cdtdz*(fz(i,j,kc,n) - fz(i,j,km,n))
 
                 rr = qxp(i,j,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 compu = rr*qxp(i,j,km,nq) - compn
                 qxpo(i,j,km,nq) = compu/rrnew
 
                 rr = qyp(i,j,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 compu = rr*qyp(i,j,km,nq) - compn
                 qypo(i,j,km,nq) = compu/rrnew
 
                 rr = qxm(i+1,j,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 compu = rr*qxm(i+1,j,km,nq) - compn
                 qxmo(i+1,j,km,nq) = compu/rrnew
 
                 rr = qym(i,j+1,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 compu = rr*qym(i,j+1,km,nq) - compn
                 qymo(i,j+1,km,nq) = compu/rrnew
 
             enddo
          enddo
       endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,nq,i,j,compn,rr,rrnew,compu) IF(nadv.gt.1)
      do iadv = 1, nadv
          n = UFA + iadv - 1
          nq = QFA + iadv - 1

          do j = jlo, jhi 
              do i = ilo, ihi 

                 compn = cdtdz*(fz(i,j,kc,n) - fz(i,j,km,n))

                 rr = qxp(i,j,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 compu = rr*qxp(i,j,km,nq) - compn
                 qxpo(i,j,km,nq) = compu/rrnew

                 rr = qyp(i,j,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 compu = rr*qyp(i,j,km,nq) - compn
                 qypo(i,j,km,nq) = compu/rrnew

                 rr = qxm(i+1,j,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 compu = rr*qxm(i+1,j,km,nq) - compn
                 qxmo(i+1,j,km,nq) = compu/rrnew

                 rr = qym(i,j+1,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 compu = rr*qym(i,j+1,km,nq) - compn
                 qymo(i,j+1,km,nq) = compu/rrnew

              enddo
          enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,nq,i,j,compsn,rr,rrnew,comps) IF(nspec.gt.1)
      do ispec = 1, nspec 
          n = UFS + ispec - 1
          nq = QFS + ispec  - 1

          do j = jlo, jhi 
              do i = ilo, ihi 

                 compsn = cdtdz*(fz(i,j,kc,n) - fz(i,j,km,n))

                 rr = qxp(i,j,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 comps = rr*qxp(i,j,km,nq) - compsn
                 qxpo(i,j,km,nq) = comps/rrnew

                 rr = qyp(i,j,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 comps = rr*qyp(i,j,km,nq) - compsn
                 qypo(i,j,km,nq) = comps/rrnew

                 rr = qxm(i+1,j,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 comps = rr*qxm(i+1,j,km,nq) - compsn
                 qxmo(i+1,j,km,nq) = comps/rrnew

                 rr = qym(i,j+1,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 comps = rr*qym(i,j+1,km,nq) - compsn
                 qymo(i,j+1,km,nq) = comps/rrnew

              enddo
          enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,nq,i,j,compsn,rr,rrnew,comps) IF(naux.gt.1)
      do iaux = 1, naux 
          n  = UFX + iaux - 1
          nq = QFX + iaux  - 1

          do j = jlo, jhi 
              do i = ilo, ihi 

                 compsn = cdtdz*(fz(i,j,kc,n) - fz(i,j,km,n))

                 rr = qxp(i,j,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 comps = rr*qxp(i,j,km,nq) - compsn
                 qxpo(i,j,km,nq) = comps/rrnew

                 rr = qyp(i,j,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 comps = rr*qyp(i,j,km,nq) - compsn
                 qypo(i,j,km,nq) = comps/rrnew

                 rr = qxm(i+1,j,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 comps = rr*qxm(i+1,j,km,nq) - compsn
                 qxmo(i+1,j,km,nq) = comps/rrnew

                 rr = qym(i,j+1,km,QRHO)
                 rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                 comps = rr*qym(i,j+1,km,nq) - compsn
                 qymo(i,j+1,km,nq) = comps/rrnew

              enddo
          enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,pgp,pgm,ugp,ugm,rrrx,rurx,rvrx,rwrx,ekenrx,rerx,rrry,rury) &
      !$OMP PRIVATE(rvry,rwry,ekenry,rery,rrlx,rulx,rvlx,rwlx,ekenlx,relx,rrly,ruly,rvly,rwly,ekenly)&
      !$OMP PRIVATE(rely,rrnewrx,runewrx,rvnewrx,rwnewrx,renewrx,rrnewry,runewry,rvnewry,rwnewry)&
      !$OMP PRIVATE(renewry,rrnewlx,runewlx,rvnewlx,rwnewlx,renewlx,rrnewly,runewly,rvnewly,rwnewly)&
      !$OMP PRIVATE(renewly,dup,pav,du,pnewrx,pnewlx,pnewry,pnewly,rhoekenrx,rhoekenry,rhoekenlx,rhoekenly)
      do j = jlo, jhi 
          do i = ilo, ihi 

             pgp = pgdnvz(i,j,kc)
             pgm = pgdnvz(i,j,km)
             ugp = ugdnvz(i,j,kc)
             ugm = ugdnvz(i,j,km)

             ! Convert to conservation form
             rrrx = qxp(i,j,km,QRHO)
             rurx = rrrx*qxp(i,j,km,QU)
             rvrx = rrrx*qxp(i,j,km,QV)
             rwrx = rrrx*qxp(i,j,km,QW)
             ekenrx = 0.5d0*rrrx*(qxp(i,j,km,QU)**2 + qxp(i,j,km,QV)**2 &
                  + qxp(i,j,km,QW)**2)
             rerx = qxp(i,j,km,QREINT) + ekenrx

             rrry = qyp(i,j,km,QRHO)
             rury = rrry*qyp(i,j,km,QU)
             rvry = rrry*qyp(i,j,km,QV)
             rwry = rrry*qyp(i,j,km,QW)
             ekenry = 0.5d0*rrry*(qyp(i,j,km,QU)**2 + qyp(i,j,km,QV)**2 &
                  + qyp(i,j,km,QW)**2)
             rery = qyp(i,j,km,QREINT) + ekenry

             rrlx = qxm(i+1,j,km,QRHO)
             rulx = rrlx*qxm(i+1,j,km,QU)
             rvlx = rrlx*qxm(i+1,j,km,QV)
             rwlx = rrlx*qxm(i+1,j,km,QW)
             ekenlx = 0.5d0*rrlx*(qxm(i+1,j,km,QU)**2 + qxm(i+1,j,km,QV)**2 &
                  + qxm(i+1,j,km,QW)**2)
             relx = qxm(i+1,j,km,QREINT) + ekenlx

             rrly = qym(i,j+1,km,QRHO)
             ruly = rrly*qym(i,j+1,km,QU)
             rvly = rrly*qym(i,j+1,km,QV)
             rwly = rrly*qym(i,j+1,km,QW)
             ekenly = 0.5d0*rrly*(qym(i,j+1,km,QU)**2 + qym(i,j+1,km,QV)**2 &
                  + qym(i,j+1,km,QW)**2)
             rely = qym(i,j+1,km,QREINT) + ekenly

             ! Add transverse predictor
             rrnewrx = rrrx - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewrx = rurx - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewrx = rvrx - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewrx = rwrx - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewrx = rerx - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))

             rrnewry = rrry - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewry = rury - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewry = rvry - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewry = rwry - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewry = rery - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))

             rrnewlx = rrlx - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewlx = rulx - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewlx = rvlx - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewlx = rwlx - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewlx = relx - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))

             rrnewly = rrly - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewly = ruly - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewly = rvly - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewly = rwly - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewly = rely - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))

             dup = pgp*ugp - pgm*ugm
             pav = 0.5d0*(pgp+pgm)
             du = ugp-ugm

             ! Convert back to non-conservation form
             if (i.ge.ilo+1) then
                qxpo(i,j,km,QRHO) = rrnewrx
                qxpo(i,j,km,QU) = runewrx/qxpo(i,j,km,QRHO)
                qxpo(i,j,km,QV) = rvnewrx/qxpo(i,j,km,QRHO)
                qxpo(i,j,km,QW) = rwnewrx/qxpo(i,j,km,QRHO)
                rhoekenrx = 0.5d0*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,km,QRHO)
                qxpo(i,j,km,QREINT)= renewrx - rhoekenrx

                pnewrx = qxp(i,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))
                qxpo(i,j,km,QPRES) = max(pnewrx,small_pres)
             end if

             if (j.ge.jlo+1) then
                qypo(i,j,km,QRHO) = rrnewry
                qypo(i,j,km,QU) = runewry/qypo(i,j,km,QRHO)
                qypo(i,j,km,QV) = rvnewry/qypo(i,j,km,QRHO)
                qypo(i,j,km,QW) = rwnewry/qypo(i,j,km,QRHO)
                rhoekenry = 0.5d0*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,km,QRHO)
                qypo(i,j,km,QREINT)= renewry - rhoekenry

                pnewry = qyp(i,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))
                qypo(i,j,km,QPRES) = max(pnewry,small_pres)
             end if

             if (i.le.ihi-1) then
                qxmo(i+1,j,km,QRHO) = rrnewlx
                qxmo(i+1,j,km,QU) = runewlx/qxmo(i+1,j,km,QRHO)
                qxmo(i+1,j,km,QV) = rvnewlx/qxmo(i+1,j,km,QRHO)
                qxmo(i+1,j,km,QW) = rwnewlx/qxmo(i+1,j,km,QRHO)
                rhoekenlx = 0.5d0*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,km,QRHO)
                qxmo(i+1,j,km,QREINT)= renewlx - rhoekenlx

                pnewlx = qxm(i+1,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))
                qxmo(i+1,j,km,QPRES) = max(pnewlx,small_pres)
             end if

             if (j.le.jhi-1) then
                qymo(i,j+1,km,QRHO) = rrnewly
                qymo(i,j+1,km,QU) = runewly/qymo(i,j+1,km,QRHO)
                qymo(i,j+1,km,QV) = rvnewly/qymo(i,j+1,km,QRHO)
                qymo(i,j+1,km,QW) = rwnewly/qymo(i,j+1,km,QRHO)
                rhoekenly = 0.5d0*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,km,QRHO)
                qymo(i,j+1,km,QREINT)= renewly - rhoekenly

                pnewly = qym(i,j+1,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))
                qymo(i,j+1,km,QPRES) = max(pnewly,small_pres)
             end if

          enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine transz

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transxy(qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fxy,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                         fyx,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                         ugdnvx,pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                         ugdnvy,pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                         gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
                         srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                         grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                         hdt,cdtdx,cdtdy,ilo,ihi,jlo,jhi,kc,km,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, QESGS, QFA, QFS, QFX, &
                                     URHO, UMX, UMY, UMZ, UEDEN, UESGS, UFA, UFS, UFX, &
                                     nadv, small_pres
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
      integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
      integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
      integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
      integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
      integer ilo,ihi,jlo,jhi,km,kc,k3d

      double precision  qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fxy(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      double precision fyx(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      double precision ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)
      double precision hdt,cdtdx,cdtdy

      integer i, j
      integer n , nq
      integer iadv, ispec, iaux

      double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
      double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
      double precision rrnewr, runewr, rvnewr, rwnewr, renewr
      double precision rrnewl, runewl, rvnewl, rwnewl, renewl
      double precision pnewr, pnewl
      double precision pgxp, pgxm, ugxp, ugxm, duxp, pxav, dux, pxnew
      double precision pgyp, pgym, ugyp, ugym, duyp, pyav, duy, pynew
      double precision pgxpm, pgxmm, ugxpm, ugxmm, duxpm, pxavm, duxm, pxnewm
      double precision pgypm, pgymm, ugypm, ugymm, duypm, pyavm, duym, pynewm
      double precision compr, compl, compnr, compnl

     ! Treat K as a passively advected quantity
      if (UESGS .gt. -1) then
         n  = UESGS
         nq = QESGS
         do j = jlo, jhi
            do i = ilo, ihi
 
               rrr = qp(i,j,kc,QRHO)
               rrl = qm(i,j,kc,QRHO)
 
               compr = rrr*qp(i,j,kc,nq)
               compl = rrl*qm(i,j,kc,nq)
 
               rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                    - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
               rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                    - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))
 
               compnr = compr - cdtdx*(fxy(i+1,j,kc,n) - fxy(i,j,kc,n)) &
                    - cdtdy*(fyx(i,j+1,kc,n) - fyx(i,j,kc,n))
               compnl = compl - cdtdx*(fxy(i+1,j,km,n) - fxy(i,j,km,n)) &
                    - cdtdy*(fyx(i,j+1,km,n) - fyx(i,j,km,n))
 
               qpo(i,j,kc,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
               qmo(i,j,kc,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d-1,nq)
 
            enddo
         enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,nq,i,j,rrr,rrl,compr,compl,rrnewr,rrnewl,compnr,compnl) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = UFA + iadv - 1
         nq = QFA + iadv - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               rrr = qp(i,j,kc,QRHO)
               rrl = qm(i,j,kc,QRHO)

               compr = rrr*qp(i,j,kc,nq)
               compl = rrl*qm(i,j,kc,nq)

               rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                    - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
               rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                    - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))
                 
               compnr = compr - cdtdx*(fxy(i+1,j,kc,n) - fxy(i,j,kc,n)) &
                    - cdtdy*(fyx(i,j+1,kc,n) - fyx(i,j,kc,n))
               compnl = compl - cdtdx*(fxy(i+1,j,km,n) - fxy(i,j,km,n)) &
                    - cdtdy*(fyx(i,j+1,km,n) - fyx(i,j,km,n))

               qpo(i,j,kc,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
               qmo(i,j,kc,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d-1,nq)

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,nq,i,j,rrr,rrl,compr,compl,rrnewr,rrnewl,compnr,compnl) IF(nspec.gt.1)
      do ispec = 1, nspec
         n = UFS + ispec - 1
         nq = QFS + ispec - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               rrr = qp(i,j,kc,QRHO)
               rrl = qm(i,j,kc,QRHO)

               compr = rrr*qp(i,j,kc,nq)
               compl = rrl*qm(i,j,kc,nq)

               rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                    - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
               rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                    - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))

               compnr = compr - cdtdx*(fxy(i+1,j,kc,n) - fxy(i,j,kc,n)) &
                              - cdtdy*(fyx(i,j+1,kc,n) - fyx(i,j,kc,n))
               compnl = compl - cdtdx*(fxy(i+1,j,km,n) - fxy(i,j,km,n)) &
                              - cdtdy*(fyx(i,j+1,km,n) - fyx(i,j,km,n))

               qpo(i,j,kc,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
               qmo(i,j,kc,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d-1,nq)

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,nq,i,j,rrr,rrl,compr,compl,rrnewr,rrnewl,compnr,compnl) IF(naux.gt.1)
      do iaux = 1, naux
         n  = UFX + iaux - 1
         nq = QFX + iaux - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               rrr = qp(i,j,kc,QRHO)
               rrl = qm(i,j,kc,QRHO)

               compr = rrr*qp(i,j,kc,nq)
               compl = rrl*qm(i,j,kc,nq)

               rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                    - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
               rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                    - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))

               compnr = compr - cdtdx*(fxy(i+1,j,kc,n) - fxy(i,j,kc,n)) &
                              - cdtdy*(fyx(i,j+1,kc,n) - fyx(i,j,kc,n))
               compnl = compl - cdtdx*(fxy(i+1,j,km,n) - fxy(i,j,km,n)) &
                              - cdtdy*(fyx(i,j+1,km,n) - fyx(i,j,km,n))

               qpo(i,j,kc,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
               qmo(i,j,kc,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d-1,nq)

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,pgxp,pgxm,ugxp,ugxm,pgyp,pgym,ugyp,ugym,pgxpm,pgxmm,ugxpm)&
      !$OMP PRIVATE(ugxmm,pgypm,pgymm,ugypm,ugymm,rrr,rur,rvr,rwr,ekenr,rer,rrl,rul,rvl,rwl,ekenl,rel)&
      !$OMP PRIVATE(rrnewr,runewr,rvnewr,rwnewr,renewr,rrnewl,runewl,rvnewl,rwnewl,renewl,duxp,pxav)&
      !$OMP PRIVATE(dux,pxnew,duxpm,pxavm,duxm,pxnewm,duyp,pyav,duy,pynew,duypm,pyavm,duym,pynewm)&
      !$OMP PRIVATE(pnewr,pnewl,rhoekenr,rhoekenl)
      do j = jlo, jhi 
         do i = ilo, ihi 

            pgxp = pgdnvx(i+1,j,kc)
            pgxm = pgdnvx(i,j,kc)
            ugxp = ugdnvx(i+1,j,kc)
            ugxm = ugdnvx(i,j,kc)

            pgyp = pgdnvy(i,j+1,kc)
            pgym = pgdnvy(i,j,kc)
            ugyp = ugdnvy(i,j+1,kc)
            ugym = ugdnvy(i,j,kc)

            pgxpm = pgdnvx(i+1,j,km)
            pgxmm = pgdnvx(i,j,km)
            ugxpm = ugdnvx(i+1,j,km)
            ugxmm = ugdnvx(i,j,km)

            pgypm = pgdnvy(i,j+1,km)
            pgymm = pgdnvy(i,j,km)
            ugypm = ugdnvy(i,j+1,km)
            ugymm = ugdnvy(i,j,km)

            ! Convert to conservation form
            rrr = qp(i,j,kc,QRHO)
            rur = rrr*qp(i,j,kc,QU)
            rvr = rrr*qp(i,j,kc,QV)
            rwr = rrr*qp(i,j,kc,QW)
            ekenr = 0.5d0*rrr*(qp(i,j,kc,QU)**2 + qp(i,j,kc,QV)**2 + &
                 qp(i,j,kc,QW)**2)
            rer = qp(i,j,kc,QREINT) + ekenr

            rrl = qm(i,j,kc,QRHO)
            rul = rrl*qm(i,j,kc,QU)
            rvl = rrl*qm(i,j,kc,QV)
            rwl = rrl*qm(i,j,kc,QW)
            ekenl = 0.5d0*rrl*(qm(i,j,kc,QU)**2 + qm(i,j,kc,QV)**2 + &
                 qm(i,j,kc,QW)**2)
            rel = qm(i,j,kc,QREINT) + ekenl

            ! Add transverse predictor
            rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                         - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
            runewr = rur - cdtdx*(fxy(i+1,j,kc,UMX) - fxy(i,j,kc,UMX)) &
                         - cdtdy*(fyx(i,j+1,kc,UMX) - fyx(i,j,kc,UMX))
            rvnewr = rvr - cdtdx*(fxy(i+1,j,kc,UMY) - fxy(i,j,kc,UMY)) &
                         - cdtdy*(fyx(i,j+1,kc,UMY) - fyx(i,j,kc,UMY))
            rwnewr = rwr - cdtdx*(fxy(i+1,j,kc,UMZ) - fxy(i,j,kc,UMZ)) &
                         - cdtdy*(fyx(i,j+1,kc,UMZ) - fyx(i,j,kc,UMZ))
            renewr = rer - cdtdx*(fxy(i+1,j,kc,UEDEN) - fxy(i,j,kc,UEDEN)) &
                         - cdtdy*(fyx(i,j+1,kc,UEDEN) - fyx(i,j,kc,UEDEN))

            rhoekenr = 0.5d0*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr

            rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                         - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))
            runewl = rul - cdtdx*(fxy(i+1,j,km,UMX) - fxy(i,j,km,UMX)) &
                         - cdtdy*(fyx(i,j+1,km,UMX) - fyx(i,j,km,UMX))
            rvnewl = rvl - cdtdx*(fxy(i+1,j,km,UMY) - fxy(i,j,km,UMY)) &
                         - cdtdy*(fyx(i,j+1,km,UMY) - fyx(i,j,km,UMY))
            rwnewl = rwl - cdtdx*(fxy(i+1,j,km,UMZ) - fxy(i,j,km,UMZ)) &
                         - cdtdy*(fyx(i,j+1,km,UMZ) - fyx(i,j,km,UMZ))
            renewl = rel - cdtdx*(fxy(i+1,j,km,UEDEN) - fxy(i,j,km,UEDEN)) &
                         - cdtdy*(fyx(i,j+1,km,UEDEN) - fyx(i,j,km,UEDEN))

            rhoekenl = 0.5d0*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

            duxp = pgxp*ugxp - pgxm*ugxm
            pxav = 0.5d0*(pgxp+pgxm)
            dux = ugxp-ugxm
            pxnew = cdtdx*(duxp + pxav*dux*(gamc(i,j,k3d)-1.d0))

            duxpm = pgxpm*ugxpm - pgxmm*ugxmm
            pxavm = 0.5d0*(pgxpm+pgxmm)
            duxm = ugxpm-ugxmm
            pxnewm = cdtdx*(duxpm + pxavm*duxm*(gamc(i,j,k3d-1)-1.d0))

            duyp = pgyp*ugyp - pgym*ugym
            pyav = 0.5d0*(pgyp+pgym)
            duy = ugyp-ugym
            pynew = cdtdy*(duyp + pyav*duy*(gamc(i,j,k3d)-1.d0))

            duypm = pgypm*ugypm - pgymm*ugymm
            pyavm = 0.5d0*(pgypm+pgymm)
            duym = ugypm-ugymm
            pynewm = cdtdy*(duypm + pyavm*duym*(gamc(i,j,k3d-1)-1.d0))

            pnewr = qp(i,j,kc,QPRES) - pxnew - pynew
            pnewl = qm(i,j,kc,QPRES) - pxnewm - pynewm

            ! Convert back to non-conservation form
            qpo(i,j,kc,QRHO  ) = rrnewr        + hdt*srcQ(i,j,k3d,QRHO)
            qpo(i,j,kc,QU    ) = runewr/rrnewr + hdt*srcQ(i,j,k3d,QU)  + hdt*grav(i,j,k3d,1)
            qpo(i,j,kc,QV    ) = rvnewr/rrnewr + hdt*srcQ(i,j,k3d,QV)  + hdt*grav(i,j,k3d,2)
            qpo(i,j,kc,QW    ) = rwnewr/rrnewr + hdt*srcQ(i,j,k3d,QW)  + hdt*grav(i,j,k3d,3)
            qpo(i,j,kc,QREINT) = renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)

            qpo(i,j,kc,QPRES ) = pnewr         + hdt*srcQ(i,j,k3d,QPRES)
            qpo(i,j,kc,QPRES) = max(qpo(i,j,kc,QPRES),small_pres)

            qmo(i,j,kc,QRHO  ) = rrnewl        + hdt*srcQ(i,j,k3d-1,QRHO)
            qmo(i,j,kc,QU    ) = runewl/rrnewl + hdt*srcQ(i,j,k3d-1,QU) + hdt*grav(i,j,k3d-1,1)
            qmo(i,j,kc,QV    ) = rvnewl/rrnewl + hdt*srcQ(i,j,k3d-1,QV) + hdt*grav(i,j,k3d-1,2)
            qmo(i,j,kc,QW    ) = rwnewl/rrnewl + hdt*srcQ(i,j,k3d-1,QW) + hdt*grav(i,j,k3d-1,3)
            qmo(i,j,kc,QREINT) = renewl - rhoekenl + hdt*srcQ(i,j,k3d-1,QREINT)

            qmo(i,j,kc,QPRES ) = pnewl         + hdt*srcQ(i,j,k3d-1,QPRES)
            qmo(i,j,kc,QPRES) = max(qmo(i,j,kc,QPRES),small_pres)

         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine transxy

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transxz(qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fxz,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                         fzx,fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
                         ugdnvx,pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                         ugdnvz,pgdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                         gamc,gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3, &
                         srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
                         grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                         hdt,cdtdx,cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, QESGS, QFA, QFS, QFX, &
                                     URHO, UMX, UMY, UMZ, UEDEN, UESGS, UFA, UFS, UFX, &
                                     nadv, small_pres
      implicit none      

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
      integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
      integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
      integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
      integer gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
      integer ilo,ihi,jlo,jhi,km,kc,k3d

      double precision  qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fxz(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      double precision fzx(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
      double precision ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
      double precision ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      double precision pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2,gc_l3:gc_h3)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)
      double precision hdt,cdtdx,cdtdz

      integer i, j
      integer n, nq
      integer iadv, ispec, iaux

      double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
      double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
      double precision rrnewr, runewr, rvnewr, rwnewr, renewr
      double precision rrnewl, runewl, rvnewl, rwnewl, renewl
      double precision pnewr, pnewl
      double precision pgxp, pgxm, ugxp, ugxm, duxp, pxav, dux, pxnew
      double precision pgzp, pgzm, ugzp, ugzm, duzp, pzav, duz, pznew
      double precision compr, compl, compnr, compnl

      ! Treat K as a passively advected quantity
      if (UESGS .gt. -1) then
         n  = UESGS
         nq = QESGS
         do j = jlo, jhi
            do i = ilo, ihi
 
               rrr = qp(i,j,km,QRHO)
               rrl = qm(i,j+1,km,QRHO)
 
               compr = rrr*qp(i,j,km,nq)
               compl = rrl*qm(i,j+1,km,nq)
 
               rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                    - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
               rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                    - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
 
               compnr = compr - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                    - cdtdz*(fzx(i,j,kc,n) - fzx(i,j,km,n))
               compnl = compl - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                    - cdtdz*(fzx(i,j,kc,n) - fzx(i,j,km,n))
 
               qpo(i,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
               qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)
 
            enddo
         enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,nq,i,j,rrr,rrl,compr,compl,rrnewr,rrnewl,compnr,compnl) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = UFA + iadv - 1
         nq = QFA + iadv -1 

         do j = jlo, jhi 
            do i = ilo, ihi 

               rrr = qp(i,j,km,QRHO)
               rrl = qm(i,j+1,km,QRHO)

               compr = rrr*qp(i,j,km,nq)
               compl = rrl*qm(i,j+1,km,nq)

               rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                    - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
               rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                    - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))

               compnr = compr - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                    - cdtdz*(fzx(i,j,kc,n) - fzx(i,j,km,n))
               compnl = compl - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                    - cdtdz*(fzx(i,j,kc,n) - fzx(i,j,km,n))

               qpo(i,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
               qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

              enddo
          enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,nq,i,j,rrr,rrl,compr,compl,rrnewr,rrnewl,compnr,compnl) IF(nspec.gt.1)
      do ispec = 1, nspec
         n = UFS + ispec - 1
         nq = QFS + ispec - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               rrr = qp(i,j,km,QRHO)
               rrl = qm(i,j+1,km,QRHO)

               compr = rrr*qp(i,j,km,nq)
               compl = rrl*qm(i,j+1,km,nq)

               rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                    - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
               rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                    - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
                 
               compnr = compr - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                              - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))
               compnl = compl - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                              - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))

               qpo(i,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
               qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,nq,i,j,rrr,rrl,compr,compl,rrnewr,rrnewl,compnr,compnl) IF(naux.gt.1)
      do iaux = 1, naux
         n  = UFX + iaux - 1
         nq = QFX + iaux - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               rrr = qp(i,j,km,QRHO)
               rrl = qm(i,j+1,km,QRHO)

               compr = rrr*qp(i,j  ,km,nq)
               compl = rrl*qm(i,j+1,km,nq)

               rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                    - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
               rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                    - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
                 
               compnr = compr - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                              - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))
               compnl = compl - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                              - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))

               qpo(i,j  ,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
               qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,pgxp,pgxm,ugxp,ugxm,pgzp,pgzm,ugzp,ugzm,rrr,rur,rvr,rwr)&
      !$OMP PRIVATE(ekenr,rer,rrl,rul,rvl,rwl,ekenl,rel,rrnewr,runewr,rvnewr,rwnewr,renewr,rrnewl)&
      !$OMP PRIVATE(runewl,rvnewl,rwnewl,renewl,duxp,pxav,dux,pxnew,duzp,pzav,duz,pznew,pnewr,pnewl)&
      !$OMP PRIVATE(rhoekenr,rhoekenl)
      do j = jlo, jhi 
         do i = ilo, ihi 

            pgxp = pgdnvx(i+1,j,km)
            pgxm = pgdnvx(i,j,km)
            ugxp = ugdnvx(i+1,j,km)
            ugxm = ugdnvx(i,j,km)

            pgzp = pgdnvz(i,j,kc)
            pgzm = pgdnvz(i,j,km)
            ugzp = ugdnvz(i,j,kc)
            ugzm = ugdnvz(i,j,km)

            ! Convert to conservation form
            rrr = qp(i,j,km,QRHO)
            rur = rrr*qp(i,j,km,QU)
            rvr = rrr*qp(i,j,km,QV)
            rwr = rrr*qp(i,j,km,QW)
            ekenr = 0.5d0*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + qp(i,j,km,QW)**2)
            rer = qp(i,j,km,QREINT) + ekenr

            rrl = qm(i,j+1,km,QRHO)
            rul = rrl*qm(i,j+1,km,QU)
            rvl = rrl*qm(i,j+1,km,QV)
            rwl = rrl*qm(i,j+1,km,QW)
            ekenl = 0.5d0*rrl*(qm(i,j+1,km,QU)**2 + qm(i,j+1,km,QV)**2 + qm(i,j+1,km,QW)**2)
            rel = qm(i,j+1,km,QREINT) + ekenl

            ! Add transverse predictor
            rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                 - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
            runewr = rur - cdtdx*(fxz(i+1,j,km,UMX) - fxz(i,j,km,UMX)) &
                 - cdtdz*(fzx(i,j,kc,UMX) - fzx(i,j,km,UMX))
            rvnewr = rvr - cdtdx*(fxz(i+1,j,km,UMY) - fxz(i,j,km,UMY)) &
                 - cdtdz*(fzx(i,j,kc,UMY) - fzx(i,j,km,UMY))
            rwnewr = rwr - cdtdx*(fxz(i+1,j,km,UMZ) - fxz(i,j,km,UMZ)) &
                 - cdtdz*(fzx(i,j,kc,UMZ) - fzx(i,j,km,UMZ))
            renewr = rer - cdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
                 - cdtdz*(fzx(i,j,kc,UEDEN) - fzx(i,j,km,UEDEN))

            rhoekenr = 0.5d0*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr

            rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                 - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
            runewl = rul - cdtdx*(fxz(i+1,j,km,UMX) - fxz(i,j,km,UMX)) &
                 - cdtdz*(fzx(i,j,kc,UMX) - fzx(i,j,km,UMX))
            rvnewl = rvl - cdtdx*(fxz(i+1,j,km,UMY) - fxz(i,j,km,UMY)) &
                 - cdtdz*(fzx(i,j,kc,UMY) - fzx(i,j,km,UMY))
            rwnewl = rwl - cdtdx*(fxz(i+1,j,km,UMZ) - fxz(i,j,km,UMZ)) &
                 - cdtdz*(fzx(i,j,kc,UMZ) - fzx(i,j,km,UMZ))
            renewl = rel - cdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
                 - cdtdz*(fzx(i,j,kc,UEDEN) - fzx(i,j,km,UEDEN))

            rhoekenl = 0.5d0*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

            duxp = pgxp*ugxp - pgxm*ugxm
            pxav = 0.5d0*(pgxp+pgxm)
            dux = ugxp-ugxm
            pxnew = cdtdx*(duxp + pxav*dux*(gamc(i,j,k3d)-1.d0))

            duzp = pgzp*ugzp - pgzm*ugzm
            pzav = 0.5d0*(pgzp+pgzm)
            duz = ugzp-ugzm
            pznew = cdtdz*(duzp + pzav*duz*(gamc(i,j,k3d)-1.d0))

            pnewr = qp(i,j,km,QPRES) - pxnew - pznew
            pnewl = qm(i,j+1,km,QPRES) - pxnew - pznew

            ! Convert back to non-conservation form
            if (j.ge.jlo+1) then
               qpo(i,j,km,QRHO  ) = rrnewr        + hdt*srcQ(i,j,k3d,QRHO)
               qpo(i,j,km,QU    ) = runewr/rrnewr + hdt*srcQ(i,j,k3d,QU)  + hdt*grav(i,j,k3d,1)
               qpo(i,j,km,QV    ) = rvnewr/rrnewr + hdt*srcQ(i,j,k3d,QV)  + hdt*grav(i,j,k3d,2)
               qpo(i,j,km,QW    ) = rwnewr/rrnewr + hdt*srcQ(i,j,k3d,QW)  + hdt*grav(i,j,k3d,3)
               qpo(i,j,km,QREINT)= renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)
   
               qpo(i,j,km,QPRES ) = pnewr         + hdt*srcQ(i,j,k3d,QPRES)
               qpo(i,j,km,QPRES) = max(qpo(i,j,km,QPRES),small_pres)
            end if

            if (j.le.jhi-1) then
               qmo(i,j+1,km,QRHO  ) = rrnewl        + hdt*srcQ(i,j,k3d,QRHO)
               qmo(i,j+1,km,QU    ) = runewl/rrnewl + hdt*srcQ(i,j,k3d,QU) + hdt*grav(i,j,k3d,1)
               qmo(i,j+1,km,QV    ) = rvnewl/rrnewl + hdt*srcQ(i,j,k3d,QV) + hdt*grav(i,j,k3d,2)
               qmo(i,j+1,km,QW    ) = rwnewl/rrnewl + hdt*srcQ(i,j,k3d,QW) + hdt*grav(i,j,k3d,3)
               qmo(i,j+1,km,QREINT)= renewl - rhoekenl + hdt*srcQ(i,j,k3d,QREINT)
               qmo(i,j+1,km,QPRES ) = pnewl         + hdt*srcQ(i,j,k3d,QPRES)
               qmo(i,j+1,km,QPRES) = max(qmo(i,j+1,km,QPRES),small_pres)
            end if

         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine transxz

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transyz(qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         fyz,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                         fzy,fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
                         ugdnvy,pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                         ugdnvz,pgdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                         gamc,gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3, &
                         srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
                         grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                         hdt,cdtdy,cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                     QPRES, QREINT, QESGS, QFA, QFS, QFX, &
                                     URHO, UMX, UMY, UMZ, UEDEN, UESGS, UFA, UFS, UFX, &
                                     nadv, small_pres
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
      integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
      integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
      integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
      integer gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
      integer ilo,ihi,jlo,jhi,km,kc,k3d

      double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fyz(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      double precision fzy(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
      double precision ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
      double precision ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      double precision pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
      double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2,gc_l3:gc_h3)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)
      double precision hdt,cdtdy,cdtdz

      integer i, j
      integer n, nq
      integer iadv, ispec, iaux

      double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
      double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
      double precision rrnewr, runewr, rvnewr, rwnewr, renewr
      double precision rrnewl, runewl, rvnewl, rwnewl, renewl
      double precision pnewr, pnewl
      double precision pgyp, pgym, ugyp, ugym, duyp, pyav, duy, pynew
      double precision pgzp, pgzm, ugzp, ugzm, duzp, pzav, duz, pznew
      double precision compr, compl, compnr, compnl

      ! Treat K as a passively advected quantity
      if (UESGS .gt. -1) then
         n  = UESGS
         nq = QESGS
         do j = jlo, jhi
            do i = ilo, ihi
    
               rrr = qp(i,j,km,QRHO)
               rrl = qm(i+1,j,km,QRHO)
 
               compr = rrr*qp(i,j,km,nq)
               compl = rrl*qm(i+1,j,km,nq)
 
               rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                    - cdtdz*(fzy(i,j,kc,URHO) - fzy(i,j,km,URHO))
               rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                    - cdtdz*(fzy(i,j,kc,URHO) - fzy(i,j,km,URHO))
 
               compnr = compr - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                    - cdtdz*(fzy(i,j,kc,n) - fzy(i,j,km,n))
               compnl = compl - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                    - cdtdz*(fzy(i,j,kc,n) - fzy(i,j,km,n))
 
               qpo(i,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
               qmo(i+1,j,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)
 
            enddo
         enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,nq,i,j,rrr,rrl,compr,compl,rrnewr,rrnewl,compnr,compnl) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = UFA + iadv - 1
         nq = QFA + iadv - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               rrr = qp(i,j,km,QRHO)
               rrl = qm(i+1,j,km,QRHO)

               compr = rrr*qp(i,j,km,nq)
               compl = rrl*qm(i+1,j,km,nq)

               rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                    - cdtdz*(fzy(i,j,kc,URHO) - fzy(i,j,km,URHO))
               rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                    - cdtdz*(fzy(i,j,kc,URHO) - fzy(i,j,km,URHO))

               compnr = compr - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                    - cdtdz*(fzy(i,j,kc,n) - fzy(i,j,km,n))
               compnl = compl - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                    - cdtdz*(fzy(i,j,kc,n) - fzy(i,j,km,n))

               qpo(i,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
               qmo(i+1,j,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,nq,i,j,rrr,rrl,compr,compl,rrnewr,rrnewl,compnr,compnl) IF(nspec.gt.1)
      do ispec = 1, nspec
         n = UFS + ispec - 1
         nq = QFS + ispec - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               rrr = qp(i  ,j,km,QRHO)
               rrl = qm(i+1,j,km,QRHO)

               compr = rrr*qp(i  ,j,km,nq)
               compl = rrl*qm(i+1,j,km,nq)

               rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                            - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))
               rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                            - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))

               compnr = compr - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                              - cdtdz*(fzy(i,j  ,kc,n) - fzy(i,j,km,n))
               compnl = compl - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                              - cdtdz*(fzy(i,j  ,kc,n) - fzy(i,j,km,n))

               qpo(i  ,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
               qmo(i+1,j,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,nq,i,j,rrr,rrl,compr,compl,rrnewr,rrnewl,compnr,compnl) IF(naux.gt.1)
      do iaux = 1, naux
         n  = UFX + iaux - 1
         nq = QFX + iaux - 1

         do j = jlo, jhi 
            do i = ilo, ihi 

               rrr = qp(i,j,km,QRHO)
               rrl = qm(i+1,j,km,QRHO)

               compr = rrr*qp(i  ,j,km,nq)
               compl = rrl*qm(i+1,j,km,nq)

               rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                            - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))
               rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                            - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))

               compnr = compr - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                              - cdtdz*(fzy(i,j  ,kc,n) - fzy(i,j,km,n))
               compnl = compl - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                              - cdtdz*(fzy(i,j  ,kc,n) - fzy(i,j,km,n))

               qpo(i  ,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
               qmo(i+1,j,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,pgyp,pgym,ugyp,ugym,pgzp,pgzm,ugzp,ugzm,rrr,rur,rvr,rwr)&
      !$OMP PRIVATE(ekenr,rer,rrl,rul,rvl,rwl,ekenl,rel,rrnewr,runewr,rvnewr,rwnewr,renewr,rrnewl)&
      !$OMP PRIVATE(runewl,rvnewl,rwnewl,renewl,duyp,pyav,duy,pynew,duzp,pzav,duz,pznew,pnewr,pnewl)&
      !$OMP PRIVATE(rhoekenr,rhoekenl)
      do j = jlo, jhi 
         do i = ilo, ihi 

            pgyp = pgdnvy(i,j+1,km)
            pgym = pgdnvy(i,j,km)
            ugyp = ugdnvy(i,j+1,km)
            ugym = ugdnvy(i,j,km)

            pgzp = pgdnvz(i,j,kc)
            pgzm = pgdnvz(i,j,km)
            ugzp = ugdnvz(i,j,kc)
            ugzm = ugdnvz(i,j,km)

            ! Convert to conservation form
            rrr = qp(i,j,km,QRHO)
            rur = rrr*qp(i,j,km,QU)
            rvr = rrr*qp(i,j,km,QV)
            rwr = rrr*qp(i,j,km,QW)
            ekenr = 0.5d0*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + &
                 qp(i,j,km,QW)**2)
            rer = qp(i,j,km,QREINT) + ekenr

            rrl = qm(i+1,j,km,QRHO)
            rul = rrl*qm(i+1,j,km,QU)
            rvl = rrl*qm(i+1,j,km,QV)
            rwl = rrl*qm(i+1,j,km,QW)
            ekenl = 0.5d0*rrl*(qm(i+1,j,km,QU)**2 + qm(i+1,j,km,QV)**2 + &
                 qm(i+1,j,km,QW)**2)
            rel = qm(i+1,j,km,QREINT) + ekenl

            ! Add transverse predictor
            rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                 - cdtdz*(fzy(i,j,kc,URHO) - fzy(i,j,km,URHO))
            runewr = rur - cdtdy*(fyz(i,j+1,km,UMX) - fyz(i,j,km,UMX)) &
                 - cdtdz*(fzy(i,j,kc,UMX) - fzy(i,j,km,UMX))
            rvnewr = rvr - cdtdy*(fyz(i,j+1,km,UMY) - fyz(i,j,km,UMY)) &
                  - cdtdz*(fzy(i,j,kc,UMY) - fzy(i,j,km,UMY))
            rwnewr = rwr - cdtdy*(fyz(i,j+1,km,UMZ) - fyz(i,j,km,UMZ)) &
                 - cdtdz*(fzy(i,j,kc,UMZ) - fzy(i,j,km,UMZ))
            renewr = rer - cdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
                 - cdtdz*(fzy(i,j,kc,UEDEN) - fzy(i,j,km,UEDEN))

            rhoekenr = 0.5d0*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr

            rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                 - cdtdz*(fzy(i,j,kc,URHO) - fzy(i,j,km,URHO))
            runewl = rul - cdtdy*(fyz(i,j+1,km,UMX) - fyz(i,j,km,UMX)) &
                 - cdtdz*(fzy(i,j,kc,UMX) - fzy(i,j,km,UMX))
            rvnewl = rvl - cdtdy*(fyz(i,j+1,km,UMY) - fyz(i,j,km,UMY)) &
                 - cdtdz*(fzy(i,j,kc,UMY) - fzy(i,j,km,UMY))
            rwnewl = rwl - cdtdy*(fyz(i,j+1,km,UMZ) - fyz(i,j,km,UMZ)) &
                 - cdtdz*(fzy(i,j,kc,UMZ) - fzy(i,j,km,UMZ))
            renewl = rel - cdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
                 - cdtdz*(fzy(i,j,kc,UEDEN) - fzy(i,j,km,UEDEN))

            rhoekenl = 0.5d0*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

            duyp = pgyp*ugyp - pgym*ugym
            pyav = 0.5d0*(pgyp+pgym)
            duy = ugyp-ugym
            pynew = cdtdy*(duyp + pyav*duy*(gamc(i,j,k3d)-1.d0))

            duzp = pgzp*ugzp - pgzm*ugzm
            pzav = 0.5d0*(pgzp+pgzm)
            duz = ugzp-ugzm
            pznew = cdtdz*(duzp + pzav*duz*(gamc(i,j,k3d)-1.d0))

            pnewr = qp(i,j,km,QPRES) - pynew - pznew
            pnewl = qm(i+1,j,km,QPRES) - pynew - pznew

            ! Convert back to non-conservation form
           if (i.ge.ilo+1) then
               qpo(i,j,km,QRHO  ) = rrnewr        + hdt*srcQ(i,j,k3d,QRHO)
               qpo(i,j,km,QU    ) = runewr/rrnewr + hdt*srcQ(i,j,k3d,QU)  + hdt*grav(i,j,k3d,1)
               qpo(i,j,km,QV    ) = rvnewr/rrnewr + hdt*srcQ(i,j,k3d,QV)  + hdt*grav(i,j,k3d,2)
               qpo(i,j,km,QW    ) = rwnewr/rrnewr + hdt*srcQ(i,j,k3d,QW)  + hdt*grav(i,j,k3d,3)
               qpo(i,j,km,QREINT)= renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)

               qpo(i,j,km,QPRES ) = pnewr         + hdt*srcQ(i,j,k3d,QPRES)
               qpo(i,j,km,QPRES) = max(qpo(i,j,km,QPRES),small_pres)

           end if

           if (i.le.ihi-1) then
               qmo(i+1,j,km,QRHO   ) = rrnewl        + hdt*srcQ(i,j,k3d,QRHO)
               qmo(i+1,j,km,QU     ) = runewl/rrnewl + hdt*srcQ(i,j,k3d,QU)  + hdt*grav(i,j,k3d,1)
               qmo(i+1,j,km,QV     ) = rvnewl/rrnewl + hdt*srcQ(i,j,k3d,QV)  + hdt*grav(i,j,k3d,2)
               qmo(i+1,j,km,QW     ) = rwnewl/rrnewl + hdt*srcQ(i,j,k3d,QW)  + hdt*grav(i,j,k3d,3)
               qmo(i+1,j,km,QREINT ) = renewl - rhoekenl + hdt*srcQ(i,j,k3d,QREINT)
 
               qmo(i+1,j,km,QPRES  ) = pnewl         + hdt*srcQ(i,j,k3d,QPRES)
               qmo(i+1,j,km,QPRES  ) = max(qmo(i+1,j,km,QPRES),small_pres)
           end if

         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine transyz
