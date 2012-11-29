
! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine ca_umdrv(is_finest_level,time,&
                          lo,hi,domlo,domhi,&
                          uin,uin_l1,uin_h1,&
                          uout,uout_l1,uout_h1,&
                          ugdnv,ugdnv_l1,ugdnv_h1,&
                          src,src_l1,src_h1, &
                          grav,gv_l1,gv_h1, &
                          delta,dt,&
                          flux,flux_l1,flux_h1,&
                          area,area_l1,area_h1,&
                          dloga,dloga_l1,dloga_h1,&
                          vol,vol_l1,vol_h1,courno,verbose,&
                          mass_added,eint_added,eden_added,&
                          E_added_flux, E_added_grav)

      use meth_params_module, only : QVAR, QU, NVAR, NHYP, URHO, use_colglaz, do_sponge, &
                                     normalize_species

      implicit none

      integer is_finest_level
      integer lo(1),hi(1),verbose
      integer domlo(1),domhi(1)
      integer uin_l1,uin_h1
      integer uout_l1,uout_h1
      integer ugdnv_l1,ugdnv_h1
      integer flux_l1,flux_h1
      integer area_l1,area_h1
      integer dloga_l1,dloga_h1
      integer vol_l1,vol_h1
      integer src_l1,src_h1
      integer gv_l1,gv_h1
      double precision   uin(  uin_l1:  uin_h1,NVAR)
      double precision  uout( uout_l1: uout_h1,NVAR)
      double precision ugdnv(ugdnv_l1:ugdnv_h1)
      double precision   src(  src_l1:  src_h1,NVAR)
      double precision  grav(   gv_l1:   gv_h1     )
      double precision  flux( flux_l1: flux_h1,NVAR)
      double precision  area( area_l1: area_h1     )
      double precision dloga(dloga_l1:dloga_h1     )
      double precision   vol(  vol_l1: vol_h1      )
      double precision delta(1),dt,time,courno

!     Automatic arrays for workspace
      double precision, allocatable:: q(:,:)
      double precision, allocatable:: gamc(:)
      double precision, allocatable:: flatn(:)
      double precision, allocatable:: c(:)
      double precision, allocatable:: csml(:)
      double precision, allocatable:: div(:)
      double precision, allocatable:: pgdnv(:)
      double precision, allocatable:: srcQ(:,:)
      double precision, allocatable:: pdivu(:)

      double precision :: dx,E_added_flux,E_added_grav
      double precision :: mass_added, eint_added, eden_added
      integer i,ngf,iflaten

      allocate(     q(uin_l1:uin_h1,QVAR))
      allocate(     c(uin_l1:uin_h1))
      allocate(  gamc(uin_l1:uin_h1))
      allocate( flatn(uin_l1:uin_h1))
      allocate(  csml(uin_l1:uin_h1))

      allocate(  srcQ(src_l1:src_h1,QVAR))

      allocate(   div(lo(1):hi(1)+1))
      allocate( pdivu(lo(1):hi(1)  ))
      allocate( pgdnv(lo(1):hi(1)+1))

      dx = delta(1)

      ngf = 1
      iflaten = 1

!     Translate to primitive variables, compute sound speeds
!     Note that (q,c,gamc,csml,flatn) are all dimensioned the same
!       and set to correspond to coordinates of (lo:hi)
   
      if (use_colglaz.eq.0) then

         call ctoprim(lo,hi,uin,uin_l1,uin_h1, &
                   q,c,gamc,csml,flatn,uin_l1,uin_h1, &
                   src,srcQ,src_l1,src_h1, &
                   courno,dx,dt,NHYP,ngf,iflaten)

         call umeth1d(lo,hi,domlo,domhi, &
                   q,c,gamc,csml,flatn,uin_l1,uin_h1, &
                   srcQ, src_l1, src_h1, &
                   grav, gv_l1, gv_h1, &
                   lo(1),hi(1),dx,dt, &
                   flux,flux_l1,flux_h1, &
                   pgdnv,lo(1),hi(1)+1, &
                   ugdnv,ugdnv_l1,ugdnv_h1, &
                   dloga,dloga_l1,dloga_h1)

      else

          call ctoprimcg(lo,hi,uin,uin_l1,uin_h1, &
                    q,c,gamc,csml,flatn,uin_l1,uin_h1, &
                    src,srcQ,src_l1,src_h1, &
                    courno,dx,dt,NHYP,ngf,iflaten)

          call umeth1dcg(q,c,gamc,csml,flatn,uin_l1,uin_h1, &
                    srcQ, src_l1, src_h1, &
                    grav, gv_l1, gv_h1, &
                    lo(1),hi(1),dx,dt, &
                    flux,flux_l1,flux_h1, & 
                    pgdnv,lo(1),hi(1)+1, &
                    dloga,dloga_l1,dloga_h1)

      endif

      ! Define p*divu
      do i = lo(1), hi(1)
         pdivu(i) = 0.5d0 * &
              (pgdnv(i+1)+pgdnv(i))*(ugdnv(i+1)*area(i+1)-ugdnv(i)*area(i)) / vol(i)
      end do

      ! Define divu on surroundingNodes(lo,hi)
      do i = lo(1),hi(1)+1
         div(i) = (q(i,QU)-q(i-1,QU)) / dx
      enddo

!     Conservative update
      call consup(uin,uin_l1,uin_h1, &
           uout,uout_l1,uout_h1, &
           pgdnv,lo(1),hi(1)+1, &
           src , src_l1, src_h1, &
           grav,  gv_l1,  gv_h1, &
           flux,flux_l1,flux_h1, &
           area,area_l1,area_h1, &
           vol , vol_l1, vol_h1, &
           div ,pdivu,lo,hi,dx,dt)

      ! Enforce the density >= small_dens.
      call enforce_minimum_density(uin,uin_l1,uin_h1,uout,uout_l1,uout_h1,lo,hi,&
                                   mass_added,eint_added,eden_added,verbose)

      ! Enforce that the species >= 0
      call ca_enforce_nonnegative_species(uout,uout_l1,uout_h1,lo,hi)

      ! Normalize the species
      if (normalize_species .eq. 1) &
         call normalize_new_species(uout,uout_l1,uout_h1,lo,hi)

      if (do_sponge .eq. 1) &
           call sponge(uout,uout_l1,uout_h1,lo,hi,time,dt,dx,domlo,domhi)

      deallocate(q,c,gamc,flatn,csml,srcQ,div,pdivu,pgdnv)

!     if ( (is_finest_level      .eq. 1) .and. &
!          (lo(1) .eq. 0               ) ) then
!        print *,'CEN0  ',time, uin(0,URHO)
!        print *,'CEN1  ',time, uin(1,URHO)
!        print *,'CEN2  ',time, uin(2,URHO)
!     end if

      end subroutine ca_umdrv

! ::: ---------------------------------------------------------------
! ::: :: UMETH1D     Compute hyperbolic fluxes using unsplit second
! ::: ::               order Godunov integrator.
! ::: :: 
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: c           => (const)  sound speed
! ::: :: gamc        => (const)  sound speed gamma
! ::: :: csml        => (const)  local small c val
! ::: :: flatn       => (const)  flattening parameter
! ::: :: src         => (const)  source
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux      <=  (modify) flux in X direction on X edges
! ::: ----------------------------------------------------------------

      subroutine umeth1d(lo,hi,domlo,domhi, &
                         q,c,gamc,csml,flatn,qd_l1,qd_h1, &
                         srcQ,src_l1,src_h1, &
                         grav, gv_l1, gv_h1, &
                         ilo,ihi,dx,dt, &
                         flux ,   fd_l1,   fd_h1, &
                         pgdnv,pgdnv_l1,pgdnv_h1, &
                         ugdnv,ugdnv_l1,ugdnv_h1, &
                         dloga,dloga_l1,dloga_h1)

      use meth_params_module, only : QVAR, NVAR, ppm_type

      implicit none
      integer lo(1),hi(1)
      integer domlo(1),domhi(1)
      integer dloga_l1,dloga_h1
      integer qd_l1,qd_h1
      integer src_l1,src_h1
      integer fd_l1,fd_h1
      integer pgdnv_l1,pgdnv_h1
      integer ugdnv_l1,ugdnv_h1
      integer gv_l1,gv_h1
      integer ilo,ihi
      double precision dx, dt
      double precision     q(   qd_l1:qd_h1,QVAR)
      double precision  gamc(   qd_l1:qd_h1)
      double precision flatn(   qd_l1:qd_h1)
      double precision  csml(   qd_l1:qd_h1)
      double precision     c(   qd_l1:qd_h1)
      double precision  flux(fd_l1   :fd_h1,NVAR)
      double precision  srcQ(src_l1  :src_h1,NVAR)
      double precision  grav(gv_l1   :gv_h1)
      double precision pgdnv(pgdnv_l1:pgdnv_h1)
      double precision ugdnv(ugdnv_l1:ugdnv_h1)
      double precision dloga(dloga_l1:dloga_h1)

!     Left and right state arrays (edge centered, cell centered)
      double precision, allocatable:: dq(:,:),  qm(:,:),   qp(:,:)

!     Work arrays to hold 3 planes of riemann state and conservative fluxes
      allocate ( dq(ilo-1:ihi+1,QVAR))
      allocate ( qm(ilo-1:ihi+1,QVAR))
      allocate ( qp(ilo-1:ihi+1,QVAR))

!     Trace to edges w/o transverse flux correction terms
      if (ppm_type .gt. 0) then
         call trace_ppm(q,dq,c,flatn,gamc,qd_l1,qd_h1, &
                        dloga,dloga_l1,dloga_h1, &
                        srcQ,src_l1,src_h1, &
                        grav,gv_l1,gv_h1, &
                        qm,qp,ilo-1,ihi+1, &
                        ilo,ihi,domlo,domhi,dx,dt)
      else
         call trace(q,dq,c,flatn,qd_l1,qd_h1, &
                    dloga,dloga_l1,dloga_h1, &
                    srcQ,src_l1,src_h1, &
                    grav,gv_l1,gv_h1, &
                    qm,qp,ilo-1,ihi+1, &
                    ilo,ihi,domlo,domhi,dx,dt)
      end if

!     Solve Riemann problem, compute xflux from improved predicted states 
      call cmpflx(lo, hi, domlo, domhi, &
                  qm, qp, ilo-1,ihi+1, &
                  flux ,  fd_l1, fd_h1, &
                  pgdnv,pgdnv_l1,pgdnv_h1, &
                  ugdnv,ugdnv_l1,ugdnv_h1, &
                  gamc, csml,c,qd_l1,qd_h1,ilo,ihi)

      deallocate (dq,qm,qp)

      end subroutine umeth1d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ctoprim(lo,hi,uin,uin_l1,uin_h1, &
                         q,c,gamc,csml,flatn,q_l1,q_h1,&
                         src,srcQ,src_l1,src_h1, &
                         courno,dx,dt,ngp,ngf,iflaten)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFA, UFS, UFX, &
                                     QVAR, QRHO, QU, QREINT, QPRES, QTEMP, QFA, QFS, QFX, &
                                     nadv, small_temp, allow_negative_energy

      implicit none

      double precision, parameter:: small = 1.d-8

!     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
!     if iflaten=1.  Declared dimensions of q,c,gamc,csml,flatn are given
!     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
!     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
!     routine that computes flatn).  

      integer          :: lo(1), hi(1)
      integer          :: uin_l1,uin_h1
      integer          :: q_l1,q_h1
      integer          ::  src_l1,src_h1
      double precision ::   uin(uin_l1:uin_h1,NVAR)
      double precision ::     q(  q_l1:  q_h1,QVAR)
      double precision ::     c(  q_l1:  q_h1)
      double precision ::  gamc(  q_l1:  q_h1)
      double precision ::  csml(  q_l1:  q_h1)
      double precision :: flatn(  q_l1:  q_h1)
      double precision ::   src(src_l1:src_h1,NVAR)
      double precision ::  srcQ(src_l1:src_h1,QVAR)
      double precision :: dx, dt, courno
      integer iflaten

      integer          :: i
      integer          :: pt_index(1)
      integer          :: ngp, ngf, loq(1), hiq(1)
      integer          :: n, nq
      integer          :: iadv, ispec, iaux
      double precision :: courx, courmx

      double precision, allocatable :: dpdrho(:), dpde(:) !, dpdX_er(:,:)

      loq(1) = lo(1)-ngp
      hiq(1) = hi(1)+ngp

      allocate(dpdrho(q_l1:q_h1))
      allocate(dpde  (q_l1:q_h1))
!      allocate(dpdX_er(q_l1:q_h1,nspec))

!     Make q (all but p), except put e in slot for rho.e, fix after eos call
!     The temperature is used as an initial guess for the eos call and will be overwritten
      do i = loq(1),hiq(1)

         if (uin(i,URHO) .le. 0.d0) then
            print *,'   '
            print *,'>>> Error: Castro_1d::ctoprim ',i
            print *,'>>> ... negative density ',uin(i,URHO)
            print *,'    '
            call bl_error("Error:: Castro_1d.f90 :: ctoprim")
         end if

         q(i,QRHO) = uin(i,URHO)
         q(i,QU) = uin(i,UMX)/uin(i,URHO)
!        eken = 0.5d0*q(i,QU)**2
!        q(i,QREINT) = uin(i,UEDEN)/q(i,QRHO) - eken
         q(i,QREINT) = uin(i,UEINT)/q(i,QRHO)
         q(i,QTEMP ) = uin(i,UTEMP)
      enddo

!     Load advected quatities, c, into q, assuming they arrived in uin as rho.c
      do iadv = 1, nadv
         n  = UFA + iadv - 1
         nq = QFA + iadv - 1
         q(loq(1):hiq(1),nq) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
      enddo
      
!     Load species, c, into q, assuming they arrived in uin as rho.c
      do ispec = 1, nspec
         n  = UFS + ispec - 1
         nq = QFS + ispec - 1
         q(loq(1):hiq(1),nq) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
      enddo

!     Load auxiliary variables which are needed in the EOS
      do iaux = 1, naux
         n  = UFX + iaux - 1
         nq = QFX + iaux - 1
         q(loq(1):hiq(1),nq) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
      enddo

!     Get gamc, p, T, c, csml using q state
      do i = loq(1), hiq(1)

         ! If necessary, reset the energy using small_temp
         if (allow_negative_energy .eq. 0 .and. q(i,QREINT) .le. 0.d0) then
            q(i,QTEMP) = small_temp
            call eos_given_RTX(q(i,QREINT),q(i,QPRES),q(i,QRHO),q(i,QTEMP), q(i,QFS:))
            if (q(i,QREINT) .lt. 0.d0) then
               print *,'   '
               print *,'>>> Error: Castro_1d::ctoprim ',i
               print *,'>>> ... new e from eos_given_RTX call is negative ',q(i,QREINT)
               print *,'    '
               call bl_error("Error:: Castro_1d.f90 :: ctoprim")
            end if
         end if

         pt_index(1) = i
         call eos_given_ReX(gamc(i), q(i,QPRES), c(i), q(i,QTEMP), dpdrho(i), dpde(i), &
                            q(i,QRHO), q(i,QREINT), q(i,QFS:), pt_index=pt_index)!, &
!                            dpdX_er=dpdX_er(i,:))
         csml(i) = max(small, small * c(i))
      end do

!     Make this "rho e" instead of "e"
      do i = loq(1),hiq(1)
         q(i,QREINT ) = q(i,QREINT )*q(i,QRHO)
      enddo

      ! compute srcQ terms
      do i = lo(1)-1, hi(1)+1
         srcQ(i,QRHO   ) = src(i,URHO)
         srcQ(i,QU     ) = (src(i,UMX) - q(i,QU) * srcQ(i,QRHO)) / q(i,QRHO)
         srcQ(i,QREINT ) = src(i,UEDEN) - q(i,QU) * src(i,UMX) + 0.5d0 * q(i,QU)**2 * srcQ(i,QRHO)
         srcQ(i,QPRES  ) = dpde(i) * (srcQ(i,QREINT) - q(i,QREINT)*srcQ(i,QRHO)/q(i,QRHO)) / q(i,QRHO) + &
              dpdrho(i) * srcQ(i,QRHO)! - &
!              sum(dpdX_er(i,:)*(src(i,UFS:UFS+nspec-1) - q(i,QFS:QFS+nspec-1)*srcQ(i,QRHO)))/q(i,QRHO)

         do ispec=1,nspec
            srcQ(i,QFS+ispec-1) = ( src(i,UFS+ispec-1) - q(i,QFS+ispec-1) * srcQ(i,QRHO) ) / q(i,QRHO)
         end do

         do iaux=1,naux
            srcQ(i,QFX+iaux-1) = ( src(i,UFX+iaux-1) - q(i,QFX+iaux-1) * srcQ(i,QRHO) ) / q(i,QRHO)
         end do

         do iadv=1,nadv
            srcQ(i,QFA+iadv-1) = ( src(i,UFA+iadv-1) - q(i,QFA+iadv-1) * srcQ(i,QRHO) ) / q(i,QRHO)
         end do

      end do

!     Compute running max of Courant number over grids
      courmx = courno
      do i = lo(1),hi(1)

         courx  = ( c(i)+abs(q(i,QU)) ) * dt/dx
         courmx = max( courmx, courx )

         if (courx .gt. 1.d0) then
            print *,'   '
            call bl_warning("Warning:: Castro_1d.f90 :: CFL violation in ctoprim")
            print *,'>>> ... (u+c) * dt / dx > 1 ', courx
            print *,'>>> ... at cell (i)       : ',i
            print *,'>>> ... u, c                ',q(i,QU), c(i)
            print *,'>>> ... density             ',q(i,QRHO)
         end if
      enddo
      courno = courmx

!     Compute flattening coef for slope calculations
      if (iflaten.eq.1) then
         loq(1)=lo(1)-ngf
         hiq(1)=hi(1)+ngf
         call uflaten(loq,hiq, &
                      q(q_l1,QPRES), &
                      q(q_l1,QU), &
                      flatn,q_l1,q_h1)
      else
         flatn = 1.d0
      endif

      deallocate(dpdrho,dpde)

      end subroutine ctoprim

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
      subroutine trace(q,dq,c,flatn,qd_l1,qd_h1, &
                       dloga,dloga_l1,dloga_h1, &
                       srcQ,src_l1,src_h1,&
                       grav,gv_l1,gv_h1, &
                       qxm,qxp,qpd_l1,qpd_h1, &
                       ilo,ihi,domlo,domhi,dx,dt)

      use network, only : nspec, naux
      use meth_params_module, only : iorder, QVAR, QRHO, QU, QREINT, QPRES, QFA, QFS, QFX, & 
                                     nadv, small_dens, ppm_type, fix_mass_flux
      use prob_params_module, only : physbc_lo, physbc_hi, Outflow
      implicit none

      integer domlo(1),domhi(1)
      integer ilo,ihi
      integer    qd_l1,   qd_h1
      integer dloga_l1,dloga_h1
      integer   qpd_l1,  qpd_h1
      integer   src_l1,  src_h1
      integer    gv_l1,   gv_h1
      double precision dx, dt
      double precision     q( qd_l1: qd_h1,QVAR)
      double precision  srcQ(src_l1:src_h1,QVAR)
      double precision flatn(qd_l1:qd_h1)
      double precision     c(qd_l1:qd_h1)
      double precision dloga(dloga_l1:dloga_h1)

      double precision   dq( qpd_l1: qpd_h1,QVAR)
      double precision  qxm( qpd_l1: qpd_h1,QVAR)
      double precision  qxp( qpd_l1: qpd_h1,QVAR)
      double precision grav(  gv_l1:  gv_h1)

!     Local variables
      integer i
      integer n, iadv
      integer ns, ispec, iaux

      double precision hdt,dtdx
      double precision cc, csq, rho, u, p, rhoe
      double precision drho, du, dp, drhoe

      double precision enth, alpham, alphap, alpha0r, alpha0e
      double precision spminus, spplus, spzero
      double precision apright, amright, azrright, azeright
      double precision apleft, amleft, azrleft, azeleft
      double precision acmprght, acmpleft
      double precision ascmprght, ascmpleft
      double precision sourcr,sourcp,source,courn,eta,dlogatmp

      logical :: fix_mass_flux_lo, fix_mass_flux_hi

      fix_mass_flux_lo = &
           (fix_mass_flux .eq. 1) .and. (physbc_lo(1) .eq. Outflow) .and. (ilo .eq. domlo(1))
      fix_mass_flux_hi = &
           (fix_mass_flux .eq. 1) .and. (physbc_hi(1) .eq. Outflow) .and. (ihi .eq. domhi(1))

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in trace with ppm_type != 0'
        call bl_error("Error:: Castro_1d.f90 :: trace")
      end if

      hdt = 0.5d0 * dt
      dtdx = dt/dx

      ! Compute slopes
      if (iorder .eq. 1) then

         dq(ilo-1:ihi+1,1:QVAR) = 0.d0

      else

         call uslope(q,flatn,qd_l1,qd_h1, &
              dq,qpd_l1,qpd_h1, &
              ilo,ihi,QVAR)

         call pslope(q(:,QPRES),q(:,QRHO), &
              flatn      , qd_l1, qd_h1, &
              dq(:,QPRES),qpd_l1,qpd_h1, &
              grav       , gv_l1, gv_h1, &
              ilo,ihi,dx)

      endif

      ! Compute left and right traced states
      do i = ilo-1, ihi+1 

         cc = c(i)
         csq = cc**2
         rho = q(i,QRHO)
         u = q(i,QU)
         p = q(i,QPRES)
         rhoe = q(i,QREINT)
         enth = ( (rhoe+p)/rho )/csq

         drho  = dq(i,QRHO)
         du    = dq(i,QU)
         dp    = dq(i,QPRES) 
         drhoe = dq(i,QREINT)

         alpham = 0.5d0*(dp/(rho*cc) - du)*rho/cc
         alphap = 0.5d0*(dp/(rho*cc) + du)*rho/cc
         alpha0r = drho - dp/csq
         alpha0e = drhoe - dp*enth

         if (u-cc .gt. 0.d0) then
            spminus = -1.d0
         else
            spminus = (u-cc)*dtdx
         endif
         if (u+cc .gt. 0.d0) then
            spplus = -1.d0
         else
            spplus = (u+cc)*dtdx
         endif
         if (u .gt. 0.d0) then
            spzero = -1.d0
         else
            spzero = u*dtdx
         endif

         apright = 0.5d0*(-1.d0 - spplus )*alphap
         amright = 0.5d0*(-1.d0 - spminus)*alpham
         azrright = 0.5d0*(-1.d0 - spzero )*alpha0r
         azeright = 0.5d0*(-1.d0 - spzero )*alpha0e

         if (i .ge. ilo) then
            qxp(i,QRHO  ) = rho + apright + amright + azrright
            qxp(i,QRHO  ) = max(small_dens,qxp(i,QRHO))
            qxp(i,QU    ) = u + (apright - amright)*cc/rho
            qxp(i,QPRES ) = p + (apright + amright)*csq
            qxp(i,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

            ! add non-gravitational source term
            qxp(i  ,QRHO  ) = qxp(i,QRHO  ) + hdt*srcQ(i,QRHO)
            qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
            qxp(i  ,QU    ) = qxp(i,QU    ) + hdt*srcQ(i,QU)
            qxp(i  ,QREINT) = qxp(i,QREINT) + hdt*srcQ(i,QREINT)
            qxp(i  ,QPRES ) = qxp(i,QPRES ) + hdt*srcQ(i,QPRES)
            
            ! add gravitational source term
            qxp(i  ,QU) = qxp(i,QU) + hdt*grav(i)
         end if

         if (u-cc .ge. 0.d0) then
            spminus = (u-cc)*dtdx
         else
            spminus = 1.d0
         endif
         if (u+cc .ge. 0.d0) then
            spplus = (u+cc)*dtdx
         else
            spplus = 1.d0
         endif
         if (u .ge. 0.d0) then
            spzero = u*dtdx
         else
            spzero = 1.d0
         endif

         apleft = 0.5d0*(1.d0 - spplus )*alphap
         amleft = 0.5d0*(1.d0 - spminus)*alpham
         azrleft = 0.5d0*(1.d0 - spzero )*alpha0r
         azeleft = 0.5d0*(1.d0 - spzero )*alpha0e

         if (i .le. ihi) then
            qxm(i+1,QRHO  ) = rho + apleft + amleft + azrleft
            qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
            qxm(i+1,QU    ) = u + (apleft - amleft)*cc/rho
            qxm(i+1,QPRES ) = p + (apleft + amleft)*csq
            qxm(i+1,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

            ! add non-gravitational source term
            qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + hdt*srcQ(i,QRHO)
            qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
            qxm(i+1,QU    ) = qxm(i+1,QU    ) + hdt*srcQ(i,QU)
            qxm(i+1,QREINT) = qxm(i+1,QREINT) + hdt*srcQ(i,QREINT)
            qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + hdt*srcQ(i,QPRES)

            ! add gravitational source term
             qxm(i+1,QU) = qxm(i+1,QU) + hdt*grav(i)
         end if

         if(dloga(i).ne.0)then
            courn = dtdx*(cc+abs(u))
            eta = (1.d0-courn)/(cc*dt*abs(dloga(i)))
            dlogatmp = min(eta,1.d0)*dloga(i)
            sourcr = -0.5d0*dt*rho*dlogatmp*u
            sourcp = sourcr*csq
            source = sourcp*enth
            if (i .le. ihi) then
               qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + sourcr
               qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
               qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + sourcp
               qxm(i+1,QREINT) = qxm(i+1,QREINT) + source
            end if
            if (i .ge. ilo) then
               qxp(i  ,QRHO  ) = qxp(i  ,QRHO  ) + sourcr
               qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
               qxp(i  ,QPRES ) = qxp(i  ,QPRES ) + sourcp
               qxp(i  ,QREINT) = qxp(i  ,QREINT) + source
            end if
         endif
      enddo

      ! Enforce constant mass flux rate if specified
      if (fix_mass_flux_lo) then
          qxm(ilo,QRHO  ) = q(domlo(1)-1,QRHO)
          qxm(ilo,QU    ) = q(domlo(1)-1,QU  )
          qxm(ilo,QPRES ) = q(domlo(1)-1,QPRES)
          qxm(ilo,QREINT) = q(domlo(1)-1,QREINT)
      end if

      ! Enforce constant mass flux rate if specified
      if (fix_mass_flux_hi) then
          qxp(ihi+1,QRHO  ) = q(domhi(1)+1,QRHO)
          qxp(ihi+1,QU    ) = q(domhi(1)+1,QU  )
          qxp(ihi+1,QPRES ) = q(domhi(1)+1,QPRES)
          qxp(ihi+1,QREINT) = q(domhi(1)+1,QREINT)
      end if

      do iadv = 1, nadv
         n = QFA + iadv - 1

         ! Right state
         do i = ilo,ihi+1
            u = q(i,QU)
            if (u .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = u*dtdx
            endif
            acmprght = 0.5d0*(-1.d0 - spzero )*dq(i,n)
            qxp(i,n) = q(i,n) + acmprght
         enddo
         if (fix_mass_flux_hi) qxp(ihi+1,n) = q(ihi+1,n)

         ! Left state
         do i = ilo-1,ihi
            u = q(i,QU)
            if (u .ge. 0.d0) then
               spzero = u*dtdx
            else
               spzero = 1.d0
            endif
            acmpleft = 0.5d0*(1.d0 - spzero )*dq(i,n)
            qxm(i+1,n) = q(i,n) + acmpleft
         enddo
         if (fix_mass_flux_lo) qxm(ilo,n) = q(ilo-1,n)
      enddo

      do ispec = 1, nspec
         ns = QFS + ispec - 1

         ! Right state
         do i = ilo,ihi+1
            u = q(i,QU)
            if (u .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = u*dtdx
            endif
            ascmprght = 0.5d0*(-1.d0 - spzero )*dq(i,ns)
            qxp(i,ns) = q(i,ns) + ascmprght
            if (fix_mass_flux_hi .and. i.eq.domhi(1)+1) & 
               qxp(i,ns) = q(i,ns)
         enddo
         if (fix_mass_flux_hi) qxp(ihi+1,ns) = q(ihi+1,ns)

         ! Left state
         do i = ilo-1,ihi
            u = q(i,QU)
            if (u .ge. 0.d0) then
               spzero = u*dtdx
            else
               spzero = 1.d0
            endif
            ascmpleft = 0.5d0*(1.d0 - spzero )*dq(i,ns)
            qxm(i+1,ns) = q(i,ns) + ascmpleft
         enddo
         if (fix_mass_flux_lo) qxm(ilo,ns) = q(ilo-1,ns)
      enddo

      do iaux = 1, naux
         ns = QFX + iaux - 1

         ! Right state
         do i = ilo,ihi+1
            u = q(i,QU)
            if (u .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = u*dtdx
            endif
            ascmprght = 0.5d0*(-1.d0 - spzero )*dq(i,ns)
            qxp(i,ns) = q(i,ns) + ascmprght
         enddo
         if (fix_mass_flux_hi) qxp(ihi+1,ns) = q(ihi+1,ns)

         ! Left state
         do i = ilo-1,ihi
            u = q(i,QU)
            if (u .ge. 0.d0) then
               spzero = u*dtdx
            else
               spzero = 1.d0
            endif
            ascmpleft = 0.5d0*(1.d0 - spzero )*dq(i,ns)
            qxm(i+1,ns) = q(i,ns) + ascmpleft
         enddo
         if (fix_mass_flux_lo) qxm(ilo,ns) = q(ilo-1,ns)
      end do

     end subroutine trace

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
     
      subroutine consup( &
             uin,  uin_l1,  uin_h1, &
            uout, uout_l1 ,uout_h1, &
           pgdnv,pgdnv_l1,pgdnv_h1, &
             src,  src_l1,  src_h1, &
            grav, grav_l1, grav_h1, &
            flux, flux_l1, flux_h1, &
           area,area_l1,area_h1, &
           vol,vol_l1,vol_h1, &
           div,pdivu,lo,hi,dx,dt)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : difmag, NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFS, UFX, &
                                     normalize_species

      implicit none
      integer lo(1), hi(1)
      integer   uin_l1,  uin_h1
      integer  uout_l1, uout_h1
      integer pgdnv_l1,pgdnv_h1
      integer   src_l1,  src_h1
      integer  grav_l1, grav_h1
      integer  flux_l1, flux_h1
      integer  area_l1, area_h1
      integer   vol_l1,  vol_h1
      double precision   uin(uin_l1:uin_h1,NVAR)
      double precision  uout(uout_l1:uout_h1,NVAR)
      double precision pgdnv(pgdnv_l1:pgdnv_h1)
      double precision   src(  src_l1:  src_h1,NVAR)
      double precision  grav( grav_l1: grav_h1)
      double precision  flux( flux_l1: flux_h1,NVAR)
      double precision  area( area_l1: area_h1)
      double precision    vol(vol_l1:vol_h1)
      double precision    div(lo(1):hi(1)+1)
      double precision  pdivu(lo(1):hi(1)  )
      double precision dx, dt

      integer          :: i, n
      double precision :: div1, dpdx
      double precision :: SrU,Up,SrE

      ! Normalize the species fluxes
      if (normalize_species .eq. 1) &
         call normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)

      do n = 1, NVAR
         if ( n.eq.UTEMP ) then
            flux(:,n) = 0.d0
         else
            do i = lo(1),hi(1)+1
               div1 = difmag*min(0.d0,div(i))
               flux(i,n) = flux(i,n) &
                    + dx*div1*(uin(i,n) - uin(i-1,n))
               flux(i,n) = area(i) * flux(i,n) * dt
         enddo
         endif
      enddo

      do n = 1, NVAR
         if ( n.eq.UTEMP) then
            do i = lo(1),hi(1)
               uout(i,n) = uin(i,n)
            enddo
         else
            do i = lo(1),hi(1)
               uout(i,n) = uin(i,n) &
                    + ( flux(i,n) - flux(i+1,n) ) / vol(i) &
                    + dt * src(i,n)
            enddo
         end if
      enddo

      ! Add gradp term to momentum equation
      do i = lo(1),hi(1)

!        dpdx = 0.5d0 * (area(i)+area(i+1)) * (pgdnv(i+1)-pgdnv(i)) / vol(i)
         dpdx = ( pgdnv(i+1)-pgdnv(i) ) / dx
         uout(i,UMX)   = uout(i,UMX) - dt * dpdx

      enddo

      ! Add source term to (rho e)
      do i = lo(1),hi(1)
         uout(i,UEINT) = uout(i,UEINT)  - dt * pdivu(i)
      enddo

      ! Add gravitational source terms to momentum and energy equations 
      do i = lo(1),hi(1)

         Up  = uin(i,UMX) / uin(i,URHO)
         SrU = uin(i,URHO) * grav(i)

         ! This doesn't work
         ! SrE = SrU*(Up + SrU*dt/(2.d0*rho))

         ! This works 
         ! SrE = SrU*Up 

         SrU = uin(i,URHO) * grav(i)
         SrE = uin(i,UMX ) * grav(i)

         uout(i,UMX  ) = uout(i,UMX  ) + dt * SrU
         uout(i,UEDEN) = uout(i,UEDEN) + dt * SrE

      enddo

      do i = lo(1),hi(1)+1
           flux(i,UMX) = flux(i,UMX) + dt*area(i)*pgdnv(i)
      enddo

      end subroutine consup

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine cmpflx(lo,hi,domlo,domhi, &
                        qm,qp,qpd_l1,qpd_h1, &
                        flx,flx_l1,flx_h1, &
                        pgdnv,pg_l1,pg_h1, &
                        ugdnv,ug_l1,ug_h1, &
                        gamc,csml,c,qd_l1,qd_h1,ilo,ihi)

      use meth_params_module, only : QVAR, NVAR

      implicit none
      integer lo(1),hi(1)
      integer domlo(1),domhi(1)
      integer ilo,ihi
      integer qpd_l1,qpd_h1
      integer flx_l1, flx_h1
      integer  pg_l1, pg_h1
      integer  ug_l1, ug_h1
      integer  qd_l1,  qd_h1
      double precision    qm(qpd_l1:qpd_h1, QVAR)
      double precision    qp(qpd_l1:qpd_h1, QVAR)
      double precision   flx(flx_l1:flx_h1, NVAR)
      double precision pgdnv( pg_l1: pg_h1)
      double precision ugdnv( ug_l1: ug_h1)
      double precision  gamc( qd_l1: qd_h1)
      double precision     c( qd_l1: qd_h1)
      double precision  csml( qd_l1: qd_h1)

!     Local variables
      integer i
      double precision, allocatable :: smallc(:),cavg(:),gamcp(:), gamcm(:)
      
      allocate ( smallc(ilo:ihi+1) )
      allocate ( cavg(ilo:ihi+1) )
      allocate ( gamcp(ilo:ihi+1) )
      allocate ( gamcm(ilo:ihi+1) )

      do i = ilo, ihi+1 
          smallc(i) = max( csml(i), csml(i-1) )
          cavg(i) = 0.5d0*( c(i) + c(i-1) )
          gamcm(i) = gamc(i-1)
          gamcp(i) = gamc(i)
      enddo

!     Solve Riemann problem (gdnv state passed back, but only (u,p) saved)
      call riemannus(qm, qp,qpd_l1,qpd_h1, smallc, cavg, &
                     gamcm, gamcp, flx, flx_l1, flx_h1, &
                     pgdnv, pg_l1, pg_h1, &
                     ugdnv, ug_l1, ug_h1, ilo, ihi, domlo, domhi )

      deallocate (smallc,cavg,gamcp,gamcm)

      end subroutine cmpflx

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine riemannus(ql,qr,qpd_l1,qpd_h1,smallc,cav, &
                           gamcl,gamcr,uflx,uflx_l1,uflx_h1,&
                           pgdnv,pg_l1,pg_h1,ugdnv,ug_l1,ug_h1, &
                           ilo,ihi,domlo,domhi)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QPRES, QREINT, QFA, QFS, QFX, &
                                     URHO, UMX, UEDEN, UEINT, UFA, UFS, UFX, nadv, small_dens, small_pres, &
                                     fix_mass_flux
      use prob_params_module, only : physbc_lo, physbc_hi, Outflow, Symmetry

      implicit none

      double precision, parameter:: small = 1.d-8

      integer ilo,ihi
      integer domlo(1),domhi(1)
      integer  qpd_l1,  qpd_h1
      integer   pg_l1,   pg_h1
      integer   ug_l1,   ug_h1
      integer uflx_l1, uflx_h1
      double precision ql(qpd_l1:qpd_h1, QVAR)
      double precision qr(qpd_l1:qpd_h1, QVAR)
      double precision   cav(ilo:ihi+1), smallc(ilo:ihi+1)
      double precision gamcl(ilo:ihi+1), gamcr(ilo:ihi+1)
      double precision  uflx(uflx_l1:uflx_h1, NVAR)
      double precision pgdnv( pg_l1: pg_h1)
      double precision ugdnv( ug_l1: ug_h1)

      double precision rgdnv, regdnv, ustar
      double precision rl, ul, pl, rel
      double precision rr, ur, pr, rer
      double precision wl, wr, rhoetot, scr
      double precision rstar, cstar, estar, pstar
      double precision ro, uo, po, reo, co, gamco, entho
      double precision sgnm, spin, spout, ushock, frac

      double precision wsmall, csmall
      integer iadv, n, nq
      integer k,ispec, iaux
      logical :: fix_mass_flux_lo, fix_mass_flux_hi

!     Solve Riemann Problem

      fix_mass_flux_lo = (fix_mass_flux .eq. 1) .and. (physbc_lo(1) .eq. Outflow) .and. (ilo .eq. domlo(1))
      fix_mass_flux_hi = (fix_mass_flux .eq. 1) .and. (physbc_hi(1) .eq. Outflow) .and. (ihi .eq. domhi(1))

      do k = ilo, ihi+1
         rl  = ql(k,QRHO)
         ul  = ql(k,QU)
         pl  = ql(k,QPRES)
         rel = ql(k,QREINT)
         rr  = qr(k,QRHO)
         ur  = qr(k,QU)
         pr  = qr(k,QPRES)
         rer = qr(k,QREINT)

         csmall = smallc(k)
         wsmall = small_dens*csmall
         wl = max(wsmall,sqrt(abs(gamcl(k)*pl*rl)))
         wr = max(wsmall,sqrt(abs(gamcr(k)*pr*rr)))
         pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
         pstar = max(pstar,small_pres)
         ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)

         if (ustar .gt. 0.d0) then
            ro = rl
            uo = ul
            po = pl
            reo = rel
            gamco = gamcl(k)
         else if (ustar .lt. 0.d0) then
            ro = rr
            uo = ur
            po = pr
            reo = rer
            gamco = gamcr(k)
         else
            ro = 0.5d0*(rl+rr)
            uo = 0.5d0*(ul+ur)
            po = 0.5d0*(pl+pr)
            reo = 0.5d0*(rel+rer)
            gamco = 0.5d0*(gamcl(k)+gamcr(k))
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
            scr = small*cav(k)
         else
            scr = spout-spin
         endif
         frac = (1.d0 + (spout + spin)/scr)*0.5d0
         frac = max(0.d0,min(1.d0,frac))

         rgdnv = frac*rstar + (1.d0 - frac)*ro
         ugdnv(k) = frac*ustar + (1.d0 - frac)*uo
         pgdnv(k) = frac*pstar + (1.d0 - frac)*po
         regdnv = frac*estar + (1.d0 - frac)*reo

         if (spout .lt. 0.d0) then
            rgdnv = ro
            ugdnv(k) = uo
            pgdnv(k) = po
            regdnv = reo
         endif
         if (spin .ge. 0.d0) then
            rgdnv = rstar
            ugdnv(k) = ustar
            pgdnv(k) = pstar
            regdnv = estar
         endif

         if (k.eq.0 .and. physbc_lo(1) .eq. Symmetry) ugdnv(k) = 0.d0

         if (fix_mass_flux_lo .and. k.eq.domlo(1) .and. ugdnv(k) .ge. 0.d0) then
            rgdnv    = ql(k,QRHO)
            ugdnv(k) = ql(k,QU)
            regdnv   = ql(k,QREINT)
         end if
         if (fix_mass_flux_hi .and. k.eq.domhi(1)+1 .and. ugdnv(k) .le. 0.d0) then
            rgdnv    = qr(k,QRHO)
            ugdnv(k) = qr(k,QU)
            regdnv   = qr(k,QREINT)
         end if

         ! Compute fluxes, order as conserved state (not q)
         uflx(k,URHO) = rgdnv*ugdnv(k)
         uflx(k,UMX) = uflx(k,URHO)*ugdnv(k) 

         rhoetot = regdnv + 0.5d0*rgdnv*ugdnv(k)**2 
         uflx(k,UEDEN) = ugdnv(k)*(rhoetot + pgdnv(k))
         uflx(k,UEINT) = ugdnv(k)*regdnv

         do iadv = 1, nadv
            n = UFA + iadv - 1
            nq = QFA + iadv - 1
            if (ustar .ge. 0.d0) then
               uflx(k,n) = uflx(k,URHO)*ql(k,nq)
            else
               uflx(k,n) = uflx(k,URHO)*qr(k,nq)
            endif
         enddo

         do ispec = 1, nspec
            n  = UFS + ispec - 1
            nq = QFS + ispec - 1
            if (ustar .ge. 0.d0) then
               uflx(k,n) = uflx(k,URHO)*ql(k,nq)
            else
               uflx(k,n) = uflx(k,URHO)*qr(k,nq)
            endif
         enddo

         do iaux = 1, naux
            n  = UFX + iaux - 1
            nq = QFX + iaux - 1
            if (ustar .ge. 0.d0) then
               uflx(k,n) = uflx(k,URHO)*ql(k,nq)
            else
               uflx(k,n) = uflx(k,URHO)*qr(k,nq)
            endif
         enddo

      enddo
      end subroutine riemannus

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
      subroutine uslope(q,flatn,qd_l1,qd_h1,dq,qpd_l1,qpd_h1,ilo,ihi,nv)

      use meth_params_module, only : QPRES

      implicit none

      integer ilo, ihi, nv
      integer qd_l1,qd_h1,qpd_l1,qpd_h1
      double precision q(qd_l1:qd_h1, nv)
      double precision dq(qpd_l1:qpd_h1,nv)
      double precision flatn(qd_l1:qd_h1)

!     Local arrays 
      double precision, allocatable::dsgn(:),dlim(:),df(:),dcen(:)

      integer i, n
      double precision dlft, drgt, slop, dq1
      double precision four3rd, sixth

      four3rd = 4.d0/3.d0
      sixth = 1.d0/6.d0

      allocate (dsgn(ilo-2:ihi+2))
      allocate (dlim(ilo-2:ihi+2))
      allocate (  df(ilo-2:ihi+2))
      allocate (dcen(ilo-2:ihi+2))

      do n = 1, nv 
!       if (n .ne. QPRES) then

            ! first compute Fromm slopes
            do i = ilo-2, ihi+2 
                dlft = 2.d0*(q(i  ,n) - q(i-1,n))
                drgt = 2.d0*(q(i+1,n) - q(i  ,n))
                dcen(i) = .25d0 * (dlft+drgt)
                dsgn(i) = sign(1.d0, dcen(i))
                slop = min( abs(dlft), abs(drgt) )
!                dlim(i) = cvmgp( slop, 0.d0, dlft*drgt )
                if (dlft*drgt .ge. 0.d0) then
                   dlim(i) = slop
                else
                   dlim(i) = 0.d0
                endif
                df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
            enddo

            ! now limited fourth order slopes
            do i = ilo-1, ihi+1 
                dq1 = four3rd*dcen(i) - sixth*(df(i+1) + df(i-1))
                dq(i,n) = flatn(i)* &
                            dsgn(i)*min(dlim(i),abs(dq1))
            enddo
!       end if
      enddo

      deallocate (dsgn,dlim,df,dcen)

      end subroutine uslope

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
      subroutine pslope(p,rho,flatn,qd_l1,qd_h1,dp,qpd_l1,qpd_h1,grav,gv_l1,gv_h1,ilo,ihi,dx)

      implicit none

      integer ilo, ihi
      integer  qd_l1, qd_h1
      integer qpd_l1,qpd_h1
      integer  gv_l1, gv_h1
      double precision, intent(in   ) ::      p( qd_l1: qd_h1)
      double precision, intent(in   ) ::    rho( qd_l1: qd_h1)
      double precision, intent(in   ) ::  flatn( qd_l1: qd_h1)
      double precision, intent(  out) ::     dp(qpd_l1:qpd_h1)
      double precision, intent(in   ) ::   grav( gv_l1: gv_h1)
      double precision, intent(in   ) ::  dx

!     Local arrays
      double precision, allocatable::dsgn(:),dlim(:),df(:),dcen(:)

      integer i
      double precision dlft, drgt, dp1
      double precision four3rd, sixth

      four3rd = 4.d0/3.d0
      sixth = 1.d0/6.d0

      allocate (dsgn(ilo-2:ihi+2))
      allocate (dlim(ilo-2:ihi+2))
      allocate (  df(ilo-2:ihi+2))
      allocate (dcen(ilo-2:ihi+2))

      ! first compute Fromm slopes
      do i = ilo-2, ihi+2 
          dlft = p(i  ) - p(i-1)
          drgt = p(i+1) - p(i  )

          ! Here we subtract off (rho * grav) so as not to limit that part of the slope
          dlft = dlft - 0.25d0 * (rho(i)+rho(i-1))*(grav(i)+grav(i-1))*dx
          drgt = drgt - 0.25d0 * (rho(i)+rho(i+1))*(grav(i)+grav(i+1))*dx
!         dlft = dlft - rho(i)*grav(i)*dx
!         drgt = drgt - rho(i)*grav(i)*dx

          dcen(i) = 0.5d0*(dlft+drgt)
          dsgn(i) = sign(1.d0, dcen(i))

          if (dlft*drgt .ge. 0.d0) then
             dlim(i) = 2.d0 * min( abs(dlft), abs(drgt) )
          else
             dlim(i) = 0.d0
          endif
          df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
      enddo

      if (ilo .eq. 0) then
        df(-1) = -df(0)
        df(-2) = -df(1)
      end if

      ! now limited fourth order slopes
      do i = ilo-1, ihi+1 
          dp1 = four3rd*dcen(i) - sixth*(df(i+1) + df(i-1))
          dp(i) = flatn(i)* &
                      dsgn(i)*min(dlim(i),abs(dp1))
          dp(i) = dp(i) + rho(i)*grav(i)*dx
      enddo

      ! Here we are assuming a symmetry boundary condition
      if (ilo .eq. 0) dp(-1) = -dp(0)

      deallocate (dsgn,dlim,df,dcen)

      end subroutine pslope

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine uflaten(lo,hi,p,u,flatn,q_l1,q_h1)

      use meth_params_module, only : iorder, small_pres

      implicit none
      integer lo(1),hi(1)
      integer q_l1,q_h1
      double precision p(q_l1:q_h1)
      double precision u(q_l1:q_h1)
      double precision flatn(q_l1:q_h1)

!     Local arrays
      double precision, allocatable :: dp(:), z(:), chi(:)

      integer i, ishft
      double precision shktst, zcut1, zcut2, dzcut
      double precision denom, zeta, tst, tmp

!     Knobs for detection of strong shock
      data shktst /0.33d0/
      data zcut1 /0.75d0/
      data zcut2 /0.85d0/

      allocate(dp(lo(1)-1:hi(1)+1),z(lo(1)-1:hi(1)+1),chi(lo(1)-1:hi(1)+1))

      dzcut = 1.d0/(zcut2-zcut1)

      if (iorder .eq. 3) then
         do i = lo(1),hi(1) 
            flatn(i) = 1.d0
         enddo
         return
      endif

!     x-direction flattening coef
      do i = lo(1)-1,hi(1)+1
         denom = max(small_pres,abs(p(i+2)-p(i-2)))
         dp(i) = p(i+1) - p(i-1)
         zeta = abs(dp(i))/denom
         z(i) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )
         if (u(i-1)-u(i+1) .ge. 0.d0) then
            tst = 1.d0
         else
            tst = 0.d0
         endif
         tmp = min(p(i+1),p(i-1))
         if ((abs(dp(i))/tmp).gt.shktst) then
            chi(i) = tst
         else
            chi(i) = 0.d0
         endif
      enddo
      do i = lo(1),hi(1)
         if(dp(i).gt.0.d0)then
            ishft = 1
         else
            ishft = -1
         endif
         flatn(i) = 1.d0 - &
              max(chi(i-ishft)*z(i-ishft),chi(i)*z(i))
      enddo

      deallocate(dp,z,chi)
      
      end subroutine uflaten

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_corrgsrc(lo,hi, &
           gold,gold_l1,gold_h1, &
           gnew,gnew_l1,gnew_h1, &
           uold,uold_l1,uold_h1, &
           unew,unew_l1,unew_h1,dt)

      use meth_params_module, only : NVAR, URHO, UMX, UEDEN

      implicit none

      integer lo(1),hi(1)
      integer gold_l1,gold_h1
      integer gnew_l1,gnew_h1
      integer uold_l1,uold_h1
      integer unew_l1,unew_h1
      double precision   gold(gold_l1:gold_h1)
      double precision   gnew(gnew_l1:gnew_h1)
      double precision  uold(uold_l1:uold_h1,NVAR)
      double precision  unew(unew_l1:unew_h1,NVAR)
      double precision dt

      integer i
      double precision :: rhon, Upn
      double precision :: SrU_new, SrU_old
      double precision :: SrUcorr,SrEcorr
      double precision :: Upo

      do i = lo(1),hi(1)

         ! Define old source term
         SrU_old = uold(i,URHO) * gold(i)
            
         rhon = unew(i,URHO)
         Upn  = unew(i,UMX) / rhon
         Upo  = uold(i,UMX) / uold(i,URHO)
         
         ! Define new source term
         SrU_new = rhon * gnew(i)
         
         ! Define corrections to source terms
         SrUcorr = 0.5d0*(SrU_new - SrU_old)

         unew(i,UMX) = unew(i,UMX) + dt*SrUcorr

         ! This doesn't work
         ! SrEcorr = SrUcorr*(Upn + SrUcorr*dt/(2*rhon))

         ! This works
         SrEcorr = 0.5d0*(SrU_new * Upn - SrU_old * Upo)

         unew(i,UEDEN) = unew(i,UEDEN) + SrEcorr*dt

      enddo

      end subroutine ca_corrgsrc

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_syncgsrc(lo,hi, &
           gphi,gphi_l1,gphi_h1, &
           gdphi,gdphi_l1,gdphi_h1, &
           state,state_l1,state_h1, &
           dstate,dstate_l1,dstate_h1, &
           sync_src,src_l1,src_h1,dt)

      use meth_params_module, only : NVAR, URHO, UMX

      implicit none

      integer lo(1),hi(1)
      integer gphi_l1,gphi_h1
      integer gdphi_l1,gdphi_h1
      integer state_l1,state_h1
      integer dstate_l1,dstate_h1
      integer src_l1,src_h1
      double precision     gphi(gphi_l1:gphi_h1)
      double precision    gdphi(gdphi_l1:gdphi_h1)
      double precision    state(state_l1:state_h1,NVAR)
      double precision   dstate(dstate_l1:dstate_h1,1+1)
      double precision sync_src(src_l1:src_h1,1+1)
      double precision dt

!     Note that dstate is drho and drhoU, state is the entire state, and src
!     is S_rhoU and S_rhoE

      integer i
      double precision rho_pre, rhoU_pre
      double precision gx, dgx, SrU, SrE

      do i = lo(1),hi(1)
            
         rho_pre  = state(i,URHO) - dstate(i,1)
         rhoU_pre = state(i,UMX)  - dstate(i,2)
         
         gx  = gphi(i)
         dgx = gdphi(i)

         SrU = dstate(i,1)*gx + rho_pre*dgx

         SrE = SrU * (rhoU_pre + 0.5*SrU*dt)/rho_pre
         
         sync_src(i,1) = SrU
         sync_src(i,2) = SrE

      enddo
      end subroutine ca_syncgsrc
