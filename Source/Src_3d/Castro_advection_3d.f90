! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine ca_umdrv(is_finest_level,time,lo,hi,domlo,domhi, &
           uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
           uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
           ugdnvx_out,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
           ugdnvy_out,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
           ugdnvz_out,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
           src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
           grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
           delta,dt, &
           flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
           flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
           flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
           area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
           area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
           area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
           vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
           courno,verbose,E_added_flux,E_added_grav)

      use meth_params_module, only : QVAR, NVAR, NHYP, do_sponge, &
                                     normalize_species

      ! This is used for IsoTurb only
      ! use probdata_module   , only : radiative_cooling_type

      implicit none

      integer is_finest_level
      integer lo(3),hi(3),verbose
      integer domlo(3),domhi(3)
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
      integer ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3
      integer ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3
      integer ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
      integer area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3
      integer area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3
      integer area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3
      integer vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
      double precision   uin(  uin_l1:uin_h1,    uin_l2:uin_h2,     uin_l3:uin_h3,  NVAR)
      double precision  uout( uout_l1:uout_h1,  uout_l2:uout_h2,   uout_l3:uout_h3, NVAR)
      double precision ugdnvx_out(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3)
      double precision ugdnvy_out(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3)
      double precision ugdnvz_out(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3)
      double precision   src(  src_l1:src_h1,    src_l2:src_h2,     src_l3:src_h3,  NVAR)
      double precision  grav( gv_l1:gv_h1,  gv_l2:gv_h2,   gv_l3:gv_h3,    3)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,NVAR)
      double precision flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,NVAR)
      double precision area1(area1_l1:area1_h1,area1_l2:area1_h2, area1_l3:area1_h3)
      double precision area2(area2_l1:area2_h1,area2_l2:area2_h2, area2_l3:area2_h3)
      double precision area3(area3_l1:area3_h1,area3_l2:area3_h2, area3_l3:area3_h3)
      double precision vol(vol_l1:vol_h1,vol_l2:vol_h2, vol_l3:vol_h3)
      double precision delta(3),dt,time,courno,E_added_flux,E_added_grav

      ! Automatic arrays for workspace
      double precision, allocatable:: q(:,:,:,:)
      double precision, allocatable:: gamc(:,:,:)
      double precision, allocatable:: flatn(:,:,:)
      double precision, allocatable:: c(:,:,:)
      double precision, allocatable:: csml(:,:,:)
      double precision, allocatable:: div(:,:,:)
      double precision, allocatable:: pdivu(:,:,:)
      double precision, allocatable:: srcQ(:,:,:,:)

      double precision dx,dy,dz
      integer ngq,ngf,iflaten

      allocate(     q(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,QVAR))
      allocate(  gamc(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
      allocate( flatn(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
      allocate(     c(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
      allocate(  csml(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
      allocate(   div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1))

      allocate( pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

      allocate(  srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR))

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      ngq = NHYP
      ngf = 1
      iflaten = 1

      ! 1) Translate conserved variables (u) to primitive variables (q).
      ! 2) Compute sound speeds (c) and gamma (gamc).
      !    Note that (q,c,gamc,csml,flatn) are all dimensioned the same
      !    and set to correspond to coordinates of (lo:hi)
      ! 3) Translate source terms
      call ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                   q,c,gamc,csml,flatn,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                   src,srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                   courno,dx,dy,dz,dt,ngq,ngf,iflaten)


      ! Compute hyperbolic fluxes using unsplit Godunov
      call umeth3d(q,c,gamc,csml,flatn,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                   srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                   grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                   lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),dx,dy,dz,dt, &
                   flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                   flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                   flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                   ugdnvx_out,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
                   ugdnvy_out,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
                   ugdnvz_out,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
                   pdivu)

      ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
      call divu(lo,hi,q,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                dx,dy,dz,div,lo(1),lo(2),lo(3),hi(1)+1,hi(2)+1,hi(3)+1)

      ! Conservative update
      call consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                  uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                  src ,  src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3, &
                  flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                  flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                  flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                  area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
                  area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
                  area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
                  vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
                  div,pdivu,lo,hi,dx,dy,dz,dt,E_added_flux)

      ! Add the radiative cooling -- for SGS only.
      ! if (radiative_cooling_type.eq.2) then
      !    call post_step_radiative_cooling(lo,hi,dt, &
      !         uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3)
      ! endif

      ! Enforce the density >= small_dens.
      call enforce_minimum_density(uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                                   uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                                   lo,hi,verbose)

      ! Enforce species >= 0
      call ca_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
                                          uout_h1,uout_h2,uout_h3,lo,hi)
 
      ! Re-normalize the species
      if (normalize_species .eq. 1) then
         call normalize_new_species(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                                    lo,hi)
      end if

      call add_grav_source(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                           uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                           grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                           lo,hi,dt,E_added_grav)

      ! Impose sponge
      if (do_sponge .eq. 1) then
         call sponge(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3,lo,hi, &
                     time,dt, &
                     dx,dy,dz,domlo,domhi)
      end if

      deallocate(q,gamc,flatn,c,csml,div,srcQ,pdivu)

      end subroutine ca_umdrv


! ::: ---------------------------------------------------------------
! ::: :: UMETH3D     Compute hyperbolic fluxes using unsplit second
! ::: ::               order Godunov integrator.
! ::: :: 
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: c           => (const)  sound speed
! ::: :: gamc        => (const)  cound speed gamma
! ::: :: csml        => (const)  local small c val
! ::: :: flatn       => (const)  flattening parameter
! ::: :: src         => (const)  source
! ::: :: nx          => (const)  number of cells in X direction
! ::: :: ny          => (const)  number of cells in Y direction
! ::: :: nz          => (const)  number of cells in Z direction
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dy          => (const)  grid spacing in Y direction
! ::: :: dz          => (const)  grid spacing in Z direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux1      <=  (modify) flux in X direction on X edges
! ::: :: flux2      <=  (modify) flux in Y direction on Y edges
! ::: :: flux3      <=  (modify) flux in Z direction on Z edges
! L:: ----------------------------------------------------------------

      subroutine umeth3d(q, c, gamc, csml, flatn, qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3, &
                         srcQ, src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                         grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                         ilo1, ilo2, ilo3, ihi1, ihi2, ihi3, dx, dy, dz, dt, &
                         flux1, fd1_l1, fd1_l2, fd1_l3, fd1_h1, fd1_h2, fd1_h3, &
                         flux2, fd2_l1, fd2_l2, fd2_l3, fd2_h1, fd2_h2, fd2_h3, &
                         flux3, fd3_l1, fd3_l2, fd3_l3, fd3_h1, fd3_h2, fd3_h3, &
                         ugdnvx_out,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3, &
                         ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
                         ugdnvy_out,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3, &
                         ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
                         ugdnvz_out,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3, &
                         ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
                         pdivu)

      use meth_params_module, only : QVAR, NVAR, QPRES, QRHO, QU, ppm_type
      use trace_ppm_module
      use ppm_module

      implicit none

      integer qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3
      integer src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
      integer gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3
      integer ilo1, ilo2, ilo3, ihi1, ihi2, ihi3
      integer fd1_l1, fd1_l2, fd1_l3, fd1_h1, fd1_h2, fd1_h3
      integer fd2_l1, fd2_l2, fd2_l3, fd2_h1, fd2_h2, fd2_h3
      integer fd3_l1, fd3_l2, fd3_l3, fd3_h1, fd3_h2, fd3_h3
      integer ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3
      integer ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3
      integer ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3
      integer km,kc,kt,k3d,n
      integer i,j

      double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      double precision  gamc(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      double precision  csml(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      double precision  srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      double precision  grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)
      double precision flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,fd1_l3:fd1_h3,NVAR)
      double precision flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,fd2_l3:fd2_h3,NVAR)
      double precision flux3(fd3_l1:fd3_h1,fd3_l2:fd3_h2,fd3_l3:fd3_h3,NVAR)
      double precision ugdnvx_out(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3)
      double precision ugdnvy_out(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3)
      double precision ugdnvz_out(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3)
      double precision pdivu(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3)
      double precision dx, dy, dz, dt
      double precision dtdx, dtdy, dtdz, hdt
      double precision cdtdx, cdtdy, cdtdz
      double precision hdtdx, hdtdy, hdtdz

      ! Left and right state arrays (edge centered, cell centered)
      double precision, allocatable:: dqx(:,:,:,:), dqy(:,:,:,:), dqz(:,:,:,:)
      double precision, allocatable::qxm(:,:,:,:),qym(:,:,:,:), qzm(:,:,:,:)
      double precision, allocatable::qxp(:,:,:,:),qyp(:,:,:,:), qzp(:,:,:,:)

      double precision, allocatable::qmxy(:,:,:,:),qpxy(:,:,:,:)
      double precision, allocatable::qmxz(:,:,:,:),qpxz(:,:,:,:)

      double precision, allocatable::qmyx(:,:,:,:),qpyx(:,:,:,:)
      double precision, allocatable::qmyz(:,:,:,:),qpyz(:,:,:,:)

      double precision, allocatable::qmzx(:,:,:,:),qpzx(:,:,:,:)
      double precision, allocatable::qmzy(:,:,:,:),qpzy(:,:,:,:)

      double precision, allocatable::qxl(:,:,:,:),qxr(:,:,:,:)
      double precision, allocatable::qyl(:,:,:,:),qyr(:,:,:,:)
      double precision, allocatable::qzl(:,:,:,:),qzr(:,:,:,:)

      ! Work arrays to hold 3 planes of riemann state and conservative fluxes
      double precision, allocatable::   fx(:,:,:,:),  fy(:,:,:,:), fz(:,:,:,:)

      double precision, allocatable::fxy(:,:,:,:),fxz(:,:,:,:)
      double precision, allocatable::fyx(:,:,:,:),fyz(:,:,:,:)
      double precision, allocatable::fzx(:,:,:,:),fzy(:,:,:,:)

      double precision, allocatable:: pgdnvx(:,:,:), ugdnvx(:,:,:)
      double precision, allocatable:: pgdnvxf(:,:,:), ugdnvxf(:,:,:)
      double precision, allocatable:: pgdnvtmpx(:,:,:), ugdnvtmpx(:,:,:)

      double precision, allocatable:: pgdnvy(:,:,:), ugdnvy(:,:,:)
      double precision, allocatable:: pgdnvyf(:,:,:), ugdnvyf(:,:,:)
      double precision, allocatable:: pgdnvtmpy(:,:,:), ugdnvtmpy(:,:,:)

      double precision, allocatable:: pgdnvz(:,:,:), ugdnvz(:,:,:)
      double precision, allocatable:: pgdnvtmpz1(:,:,:), ugdnvtmpz1(:,:,:)
      double precision, allocatable:: pgdnvtmpz2(:,:,:), ugdnvtmpz2(:,:,:)

      double precision, allocatable:: pgdnvzf(:,:,:), ugdnvzf(:,:,:)

      double precision, allocatable:: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)

      allocate ( pgdnvx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( ugdnvx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( pgdnvxf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( ugdnvxf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( pgdnvtmpx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( ugdnvtmpx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))

      allocate ( pgdnvy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( ugdnvy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( pgdnvyf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( ugdnvyf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( pgdnvtmpy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( ugdnvtmpy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))

      allocate ( pgdnvz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( ugdnvz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( pgdnvtmpz1(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( ugdnvtmpz1(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( pgdnvtmpz2(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( ugdnvtmpz2(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( pgdnvzf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
      allocate ( ugdnvzf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))

      allocate ( dqx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( dqy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( dqz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

      allocate ( qxm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qxp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

      allocate ( qmxy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qpxy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

      allocate ( qmxz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qpxz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

      allocate ( qym(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qyp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

      allocate ( qmyx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qpyx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

      allocate ( qmyz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qpyz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

      allocate ( qzm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qzp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

      allocate ( qxl(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qxr(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qyl(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qyr(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qzl(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qzr(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

      allocate ( qmzx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qpzx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

      allocate ( qmzy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
      allocate ( qpzy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

      allocate ( fx(ilo1:ihi1+1,ilo2-1:ihi2+1,2,NVAR))
      allocate ( fy(ilo1-1:ihi1+1,ilo2:ihi2+1,2,NVAR))
      allocate ( fz(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,NVAR))

      allocate ( fxy(ilo1:ihi1+1,ilo2-1:ihi2+1,2,NVAR))
      allocate ( fxz(ilo1:ihi1+1,ilo2-1:ihi2+1,2,NVAR))

      allocate ( fyx(ilo1-1:ihi1+1,ilo2:ihi2+1,2,NVAR))
      allocate ( fyz(ilo1-1:ihi1+1,ilo2:ihi2+1,2,NVAR))

      allocate ( fzx(ilo1:ihi1,ilo2-1:ihi2+1,2,NVAR))
      allocate ( fzy(ilo1-1:ihi1+1,ilo2:ihi2,2,NVAR))

      ! x-index, y-index, z-index, dim, characteristics, variables
      allocate ( Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,3,QVAR))
      allocate ( Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,3,QVAR))

      ! Local constants
      dtdx = dt/dx
      dtdy = dt/dy
      dtdz = dt/dz
      hdt = 0.5d0*dt
      hdtdx = 0.5d0*dtdx
      hdtdy = 0.5d0*dtdy
      hdtdz = 0.5d0*dtdz
      cdtdx = dtdx/3.d0
      cdtdy = dtdy/3.d0
      cdtdz = dtdz/3.d0

      ! Initialize pdivu to zero
      pdivu(:,:,:) = 0.d0

      ! Initialize kc (current k-level) and km (previous k-level)
      kc = 1
      km = 2

      do k3d = ilo3-1, ihi3+1

         ! Swap pointers to levels
         kt = km
         km = kc
         kc = kt

         if (ppm_type .gt. 0) then

            do n=1,QVAR
               call ppm(q(:,:,:,n),qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        q(:,:,:,QU:),c,Ip(:,:,:,:,:,n),Im(:,:,:,:,:,n), &
                        ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc,n)
            end do

            ! Compute U_x and U_y at kc (k3d)
            call tracexy_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                             Ip,Im, &
                             qxm,qxp,qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d)

         else

            ! Compute all slopes at kc (k3d)
            call uslope(q,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        dqx,dqy,dqz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        ilo1,ilo2,ihi1,ihi2,kc,k3d,QVAR)
            
            call pslope(q(:,:,:,QPRES),q(:,:,:,QRHO), &
                        flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        dqx(:,:,:,QPRES),dqy(:,:,:,QPRES),dqz(:,:,:,QPRES), &
                        ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                        ilo1,ilo2,ihi1,ihi2,kc,k3d,dx,dy,dz)

            ! Compute U_x and U_y at kc (k3d)
            call tracexy(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         dqx,dqy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         qxm,qxp,qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d)

         end if

         ! Compute \tilde{F}^x at kc (k3d)
         call cmpflx(qxm,qxp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                     fx,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                     ugdnvx,pgdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                     gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     1,ilo1,ihi1+1,ilo2-1,ihi2+1,kc,kc,k3d)

         ! Compute \tilde{F}^y at kc (k3d)
         call cmpflx(qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                     fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                     ugdnvy,pgdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                     gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     2,ilo1-1,ihi1+1,ilo2,ihi2+1,kc,kc,k3d)

         ! Compute U'^y_x at kc (k3d)
         call transy1(qxm,qmxy,qxp,qpxy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                      ugdnvy,pgdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      cdtdy,ilo1-1,ihi1+1,ilo2,ihi2,kc,k3d)

         ! Compute U'^x_y at kc (k3d)
         call transx1(qym,qmyx,qyp,qpyx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      fx,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                      ugdnvx,pgdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      cdtdx,ilo1,ihi1,ilo2-1,ihi2+1,kc,k3d)

         ! Compute F^{x|y} at kc (k3d)
         call cmpflx(qmxy,qpxy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                     fxy,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                     ugdnvtmpx,pgdnvtmpx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                     gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     1,ilo1,ihi1+1,ilo2,ihi2,kc,kc,k3d)

         ! Compute F^{y|x} at kc (k3d)
         call cmpflx(qmyx,qpyx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                     fyx,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                     ugdnvtmpy,pgdnvtmpy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                     gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     2,ilo1,ihi1,ilo2,ihi2+1,kc,kc,k3d)

         if (k3d.ge.ilo3) then

            ! Compute U_z at kc (k3d)
            if (ppm_type .gt. 0) then
               call tracez_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                               Ip,Im, &
                               qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                               ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d)
            else
               call tracez(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           dqz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d)
            end if

            ! Compute \tilde{F}^z at kc (k3d)
            call cmpflx(qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        fz,ilo1-1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                        ugdnvz,pgdnvz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        3,ilo1-1,ihi1+1,ilo2-1,ihi2+1,kc,kc,k3d)

            ! Compute U'^y_z at kc (k3d)
            call transy2(qzm,qmzy,qzp,qpzy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                         ugdnvy,pgdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         cdtdy,ilo1-1,ihi1+1,ilo2,ihi2,kc,km,k3d)

            ! Compute U'^x_z at kc (k3d)
            call transx2(qzm,qmzx,qzp,qpzx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         fx,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                         ugdnvx,pgdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         cdtdx,ilo1,ihi1,ilo2-1,ihi2+1,kc,km,k3d)

            ! Compute F^{z|x} at kc (k3d)
            call cmpflx(qmzx,qpzx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        fzx,ilo1,ilo2-1,1,ihi1,ihi2+1,2, &
                        ugdnvtmpz1,pgdnvtmpz1,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        3,ilo1,ihi1,ilo2-1,ihi2+1,kc,kc,k3d)

            ! Compute F^{z|y} at kc (k3d)
            call cmpflx(qmzy,qpzy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        fzy,ilo1-1,ilo2,1,ihi1+1,ihi2,2, &
                        ugdnvtmpz2,pgdnvtmpz2,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        3,ilo1-1,ihi1+1,ilo2,ihi2,kc,kc,k3d)

            ! Compute U''_z at kc (k3d)
            call transxy(qzm,qzl,qzp,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         fxy,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                         fyx,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                         ugdnvtmpx,pgdnvtmpx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         ugdnvtmpy,pgdnvtmpy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                         grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3,&
                         hdt,hdtdx,hdtdy,ilo1,ihi1,ilo2,ihi2,kc,km,k3d)

            ! Compute F^z at kc (k3d) -- note that flux3 is indexed by k3d, not kc
            call cmpflx(qzl,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        flux3,fd3_l1,fd3_l2,fd3_l3,fd3_h1,fd3_h2,fd3_h3, &
                        ugdnvzf,pgdnvzf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        3,ilo1,ihi1,ilo2,ihi2,kc,k3d,k3d)

            do j=ilo2-1,ihi2+1
               do i=ilo1-1,ihi1+1
                  ugdnvz_out(i,j,k3d) = ugdnvzf(i,j,kc)
               end do
            end do

            if (k3d .ge. ilo3+1 .and. k3d .le. ihi3+1) then
               do j = ilo2,ihi2
                  do i = ilo1,ihi1
                     pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
                          0.5d0*(pgdnvzf(i,j,kc)+pgdnvzf(i,j,km)) * &
                                (ugdnvzf(i,j,kc)-ugdnvzf(i,j,km))/dz
                  end do
               end do
            end if

            if (k3d.gt.ilo3) then

               ! Compute U'^z_x and U'^z_y at km (k3d-1) -- note flux3 has physical index
               call transz(qxm,qmxz,qxp,qpxz, &
                           qym,qmyz,qyp,qpyz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           fz,ilo1-1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                           ugdnvz,pgdnvz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           cdtdz,ilo1-1,ihi1+1,ilo2-1,ihi2+1,km,kc,k3d)
         
               ! Compute F^{x|z} at km (k3d-1)
               call cmpflx(qmxz,qpxz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           fxz,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                           ugdnvx,pgdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           1,ilo1,ihi1+1,ilo2-1,ihi2+1,km,km,k3d-1)

               ! Compute F^{y|z} at km (k3d-1)
               call cmpflx(qmyz,qpyz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           fyz,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                           ugdnvy,pgdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           2,ilo1-1,ihi1+1,ilo2,ihi2+1,km,km,k3d-1)

               ! Compute U''_x at km (k3d-1)
               call transyz(qxm,qxl,qxp,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                            fyz,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                            fzy,ilo1-1,ilo2,1,ihi1+1,ihi2,2, &
                            ugdnvy,pgdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                            ugdnvtmpz2,pgdnvtmpz2,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                            gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                            srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                            grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3,&
                            hdt,hdtdy,hdtdz,ilo1-1,ihi1+1,ilo2,ihi2,km,kc,k3d-1)

               ! Compute U''_y at km (k3d-1)
               call transxz(qym,qyl,qyp,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                            fxz,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                            fzx,ilo1,ilo2-1,1,ihi1,ihi2+1,2, &
                            ugdnvx,pgdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                            ugdnvtmpz1,pgdnvtmpz1,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                            gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                            srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                            grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3,&
                            hdt,hdtdx,hdtdz,ilo1,ihi1,ilo2-1,ihi2+1,km,kc,k3d-1)

               ! Compute F^x at km (k3d-1)
               call cmpflx(qxl,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           flux1,fd1_l1,fd1_l2,fd1_l3,fd1_h1,fd1_h2,fd1_h3, &
                           ugdnvxf,pgdnvxf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           1,ilo1,ihi1+1,ilo2,ihi2,km,k3d-1,k3d-1)

               do j=ilo2-1,ihi2+1
                  do i=ilo1-1,ihi1+2
                     ugdnvx_out(i,j,k3d-1) = ugdnvxf(i,j,km)
                  end do
               end do

               ! Compute F^y at km (k3d-1)
               call cmpflx(qyl,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           flux2,fd2_l1,fd2_l2,fd2_l3,fd2_h1,fd2_h2,fd2_h3, &
                           ugdnvyf,pgdnvyf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           2,ilo1,ihi1,ilo2,ihi2+1,km,k3d-1,k3d-1)

               do j=ilo2-1,ihi2+2
                  do i=ilo1-1,ihi1+1
                     ugdnvy_out(i,j,k3d-1) = ugdnvyf(i,j,km)
                  end do
               end do

               do j = ilo2,ihi2
                  do i = ilo1,ihi1
                     pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
                          0.5d0*(pgdnvxf(i+1,j,km) + pgdnvxf(i,j,km)) *  &
                          (ugdnvxf(i+1,j,km)-ugdnvxf(i,j,km))/dx + &
                          0.5d0*(pgdnvyf(i,j+1,km) + pgdnvyf(i,j,km)) *  &
                          (ugdnvyf(i,j+1,km)-ugdnvyf(i,j,km))/dy
                  end do
               end do

            end if
         end if
      enddo

      ! Deallocate arrays
      deallocate(pgdnvx,ugdnvx)
      deallocate(pgdnvxf,ugdnvxf)
      deallocate(pgdnvtmpx,ugdnvtmpx)
      deallocate(pgdnvy,ugdnvy)
      deallocate(pgdnvyf,ugdnvyf)
      deallocate(pgdnvtmpy,ugdnvtmpy)
      deallocate(pgdnvz,ugdnvz)
      deallocate(pgdnvtmpz1,ugdnvtmpz1)
      deallocate(pgdnvtmpz2,ugdnvtmpz2)
      deallocate(pgdnvzf,ugdnvzf)
      deallocate(dqx,dqy,dqz)
      deallocate(qxm,qxp)
      deallocate(qmxy,qpxy)
      deallocate(qmxz,qpxz)
      deallocate(qym,qyp)
      deallocate(qmyx,qpyx)
      deallocate(qmyz,qpyz)
      deallocate(qzm,qzp)
      deallocate(qxl,qxr,qyl,qyr,qzl,qzr)
      deallocate(qmzx,qpzx)
      deallocate(qmzy,qpzy)
      deallocate(fx,fy,fz)
      deallocate(fxy,fxz)
      deallocate(fyx,fyz)
      deallocate(fzx,fzy)
      deallocate(Ip,Im)

      end subroutine umeth3d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ctoprim(lo,hi,          uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                         q,c,gamc,csml,flatn,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
                         src,srcQ,           src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                         courno,dx,dy,dz,dt,ngp,ngf,iflaten)
      !
      !     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
      !     if iflaten=1.  Declared dimensions of q,c,gamc,csml,flatn are given
      !     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
      !     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
      !     routine that computes flatn).  
      !
      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                     UEDEN, UEINT, UESGS, UTEMP, UFA, UFS, UFX, &
                                     QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QESGS, QPRES, QTEMP, QFA, QFS, QFX, &
                                     nadv, allow_negative_energy, small_temp
      implicit none

      double precision, parameter:: small = 1.d-8

      integer lo(3), hi(3)
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3

      double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
      double precision :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
      double precision :: c(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      double precision :: gamc(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      double precision :: csml(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      double precision :: flatn(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      double precision ::  src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
      double precision :: srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      double precision :: dx, dy, dz, dt, courno
      integer          :: iflaten

      double precision, allocatable:: dpdrho(:,:,:)
      double precision, allocatable:: dpde(:,:,:)
      double precision, allocatable:: dpdX_er(:,:,:,:)

      integer          :: i, j, k
      integer          :: pt_index(3)
      integer          :: ngp, ngf, loq(3), hiq(3)
      integer          :: n, nq
      integer          :: iadv, ispec, iaux
      double precision :: courx, coury, courz, courmx, courmy, courmz

      allocate( dpdrho(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3))
      allocate(   dpde(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3))
      allocate(dpdX_er(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,nspec))

      do i=1,3
         loq(i) = lo(i)-ngp
         hiq(i) = hi(i)+ngp
      enddo
      !
      ! Make q (all but p), except put e in slot for rho.e, fix after eos call.
      ! The temperature is used as an initial guess for the eos call and will be overwritten.
      !
      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = loq(3),hiq(3)
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)

               if (uin(i,j,k,URHO) .le. 0.d0) then
                  print *,'   '
                  print *,'>>> Error: Castro_advection_3d::ctoprim ',i,j,k
                  print *,'>>> ... negative density ',uin(i,j,k,URHO)
                  call bl_error("Error:: Castro_advection_3d.f90 :: ctoprim")
               end if

               q(i,j,k,QRHO) = uin(i,j,k,URHO)
               q(i,j,k,QU) = uin(i,j,k,UMX)/uin(i,j,k,URHO)
               q(i,j,k,QV) = uin(i,j,k,UMY)/uin(i,j,k,URHO)
               q(i,j,k,QW) = uin(i,j,k,UMZ)/uin(i,j,k,URHO)
               ! convert "rho e" to "e"
               q(i,j,k,QREINT ) = uin(i,j,k,UEINT)/q(i,j,k,QRHO)
               q(i,j,k,QTEMP  ) = uin(i,j,k,UTEMP)

               ! convert "rho K" to "K"
               if (QESGS .gt. -1) &
                  q(i,j,k,QESGS) = uin(i,j,k,UESGS)/q(i,j,k,QRHO)

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      ! Load advected quatities, c, into q, assuming they arrived in uin as rho.c
      !$OMP PARALLEL DO PRIVATE(iadv,n,nq,i,j,k) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = UFA + iadv - 1
         nq = QFA + iadv - 1
         do k = loq(3),hiq(3)
            do j = loq(2),hiq(2)
               do i = loq(1),hiq(1)
                  q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
               enddo
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Load chemical species, c, into q, assuming they arrived in uin as rho.c
      !$OMP PARALLEL DO PRIVATE(ispec,n,nq,i,j,k) IF(nspec.gt.1)
      do ispec = 1, nspec
         n  = UFS + ispec - 1
         nq = QFS + ispec - 1
         do k = loq(3),hiq(3)
            do j = loq(2),hiq(2)
               do i = loq(1),hiq(1)
                  q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
               enddo
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      
      ! Load auxiliary variables which are needed in the EOS
      !$OMP PARALLEL DO PRIVATE(iaux,n,nq,i,j,k) IF(naux.gt.1)
      do iaux = 1, naux
         n  = UFX + iaux - 1
         nq = QFX + iaux - 1
         do k = loq(3),hiq(3)
            do j = loq(2),hiq(2)
               do i = loq(1),hiq(1)
                  q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
               enddo
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      ! Get gamc, p, T, c, csml using q state
      !$OMP PARALLEL DO PRIVATE(i,j,k,pt_index)
      do k = loq(3), hiq(3)
         do j = loq(2), hiq(2)
            do i = loq(1), hiq(1)
 
               pt_index(1) = i
               pt_index(2) = j
               pt_index(3) = k

               ! If necessary, reset the energy using small_temp
               if ((allow_negative_energy .eq. 0) .and. (q(i,j,k,QREINT) .lt. 0)) then
                  q(i,j,k,QTEMP) = small_temp
                  call eos_given_RTX(q(i,j,k,QREINT),q(i,j,k,QPRES),q(i,j,k,QRHO), &
                                     q(i,j,k,QTEMP),q(i,j,k,QFS:),pt_index=pt_index)
                  if (q(i,j,k,QREINT) .lt. 0.d0) then
                     print *,'   '
                     print *,'>>> Error: Castro_advection_3d::ctoprim ',i,j,k
                     print *,'>>> ... new e from eos_given_RTX call is negative ' &
                          ,q(i,j,k,QREINT)
                     print *,'    '
                     call bl_error("Error:: Castro_advection_3d.f90 :: ctoprim")
                  end if
               end if

               call eos_given_ReX(gamc(i,j,k), q(i,j,k,QPRES), c(i,j,k), q(i,j,k,QTEMP), &
                                  dpdrho(i,j,k), dpde(i,j,k), &
                                  q(i,j,k,QRHO), q(i,j,k,QREINT), q(i,j,k,QFS:), &
                                  pt_index=pt_index)!, &
!                                  dpdX_er=dpdX_er(i,j,k,:))
               csml(i,j,k) = max(small, small * c(i,j,k))

               ! convert "e" back to "rho e"
               q(i,j,k,QREINT) = q(i,j,k,QREINT)*q(i,j,k,QRHO)

            end do
         end do
      end do
      !$OMP END PARALLEL DO

      ! compute srcQ terms
      !$OMP PARALLEL DO PRIVATE(i,j,k,ispec,iaux,iadv)
      do k = lo(3)-1, hi(3)+1
         do j = lo(2)-1, hi(2)+1
            do i = lo(1)-1, hi(1)+1

               srcQ(i,j,k,QRHO  ) = src(i,j,k,URHO)
               srcQ(i,j,k,QU    ) = (src(i,j,k,UMX) - q(i,j,k,QU) * srcQ(i,j,k,QRHO)) / q(i,j,k,QRHO)
               srcQ(i,j,k,QV    ) = (src(i,j,k,UMY) - q(i,j,k,QV) * srcQ(i,j,k,QRHO)) / q(i,j,k,QRHO)
               srcQ(i,j,k,QW    ) = (src(i,j,k,UMZ) - q(i,j,k,QW) * srcQ(i,j,k,QRHO)) / q(i,j,k,QRHO)
               srcQ(i,j,k,QREINT) = src(i,j,k,UEDEN) - q(i,j,k,QU)*src(i,j,k,UMX) &
                                                     - q(i,j,k,QV)*src(i,j,k,UMY) &
                                                     - q(i,j,k,QW)*src(i,j,k,UMZ) &
                           + 0.5d0 * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2) * srcQ(i,j,k,QRHO)

               srcQ(i,j,k,QPRES ) = dpde(i,j,k)*(srcQ(i,j,k,QREINT) - &
                                                 q(i,j,k,QREINT)*srcQ(i,j,k,QRHO)/q(i,j,k,QRHO)) /q(i,j,k,QRHO) + &
                                  dpdrho(i,j,k)*srcQ(i,j,k,QRHO)! + &
!                                    sum(dpdX_er(i,j,k,:)*(src(i,j,k,UFS:UFS+nspec-1) - &
!                                                          q(i,j,k,QFS:QFS+nspec-1)*srcQ(i,j,k,QRHO))) &
!                                    /q(i,j,k,QRHO)

               if (QESGS .gt. -1) &
                  srcQ(i,j,k,QESGS) = src(i,j,k,UESGS)/q(i,j,k,QRHO) - q(i,j,k,QESGS) * srcQ(i,j,k,QRHO)

               do ispec = 1,nspec
                  srcQ(i,j,k,QFS+ispec-1) = ( src(i,j,k,UFS+ispec-1) - q(i,j,k,QFS+ispec-1) * srcQ(i,j,k,QRHO) ) / &
                                             q(i,j,k,QRHO)
               enddo

               do iaux = 1,naux
                  srcQ(i,j,k,QFX+iaux-1) = ( src(i,j,k,UFX+iaux-1) - q(i,j,k,QFX+iaux-1) * srcQ(i,j,k,QRHO) ) / &
                                             q(i,j,k,QRHO)
               enddo

               do iadv = 1,nadv
                  srcQ(i,j,k,QFA+iadv-1) = ( src(i,j,k,UFA+iadv-1) - q(i,j,k,QFA+iadv-1) * srcQ(i,j,k,QRHO) ) / &
                                             q(i,j,k,QRHO)
               enddo

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      ! Compute running max of Courant number over grids
      courmx = courno
      courmy = courno
      courmz = courno
      !$OMP PARALLEL DO PRIVATE(i,j,k,courx,coury,courz) REDUCTION(max:courmx,courmy,courmz)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               courx = ( c(i,j,k)+abs(q(i,j,k,QU)) ) * dt/dx
               coury = ( c(i,j,k)+abs(q(i,j,k,QV)) ) * dt/dy
               courz = ( c(i,j,k)+abs(q(i,j,k,QW)) ) * dt/dz

               courmx = max( courmx, courx )
               courmy = max( courmy, coury )
               courmz = max( courmz, courz )

               if (courx .gt. 1.d0) then
                  print *,'   '
                  call bl_warning("Warning:: Castro_advection_3d.f90 :: CFL violation in ctoprim")
                  print *,'>>> ... (u+c) * dt / dx > 1 ', courx
                  print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                  print *,'>>> ... u, c                ',q(i,j,k,QU), c(i,j,k)
                  print *,'>>> ... density             ',q(i,j,k,QRHO)
               end if

               if (coury .gt. 1.d0) then
                  print *,'   '
                  call bl_warning("Warning:: Castro_advection_3d.f90 :: CFL violation in ctoprim")
                  print *,'>>> ... (v+c) * dt / dx > 1 ', coury
                  print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                  print *,'>>> ... v, c                ',q(i,j,k,QV), c(i,j,k)
                  print *,'>>> ... density             ',q(i,j,k,QRHO)
               end if

               if (courz .gt. 1.d0) then
                  print *,'   '
                  call bl_warning("Warning:: Castro_advection_3d.f90 :: CFL violation in ctoprim")
                  print *,'>>> ... (w+c) * dt / dx > 1 ', courz
                  print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                  print *,'>>> ... w, c                ',q(i,j,k,QW), c(i,j,k)
                  print *,'>>> ... density             ',q(i,j,k,QRHO)
               end if

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      courno = max( courmx, courmy, courmz )

      ! Compute flattening coef for slope calculations
      if(iflaten.eq.1)then
         do n=1,3
            loq(n)=lo(n)-ngf
            hiq(n)=hi(n)+ngf
         enddo
         call uflaten(loq,hiq, &
                      q(q_l1,q_l2,q_l3,QPRES), &
                      q(q_l1,q_l2,q_l3,QU), &
                      q(q_l1,q_l2,q_l3,QV), &
                      q(q_l1,q_l2,q_l3,QW), &
                      flatn,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3)
      else
         flatn = 1.d0
      endif

      deallocate(dpdrho,dpde)

      end subroutine ctoprim

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine tracexy(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         dqx,dqy,dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3, &
                         qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QESGS, QPRES, QFA, QFS, QFX, nadv, small_dens, &
                                     ppm_type
      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2
      integer kc,k3d

      double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

      double precision  dqx(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)
      double precision  dqy(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)

      double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision dx, dy, dt

      ! Local variables
      integer i, j, n
      integer iadv, ispec, iaux

      double precision dtdx, dtdy
      double precision cc, csq, rho, u, v, w, p, rhoe
      double precision drho, du, dv, dw, dp, drhoe

      double precision enth, alpham, alphap, alpha0r, alpha0e
      double precision alpha0u, alpha0v, alpha0w
      double precision spminus, spplus, spzero
      double precision apright, amright, azrright, azeright
      double precision azu1rght, azv1rght, azw1rght
      double precision apleft, amleft, azrleft, azeleft
      double precision azu1left, azv1left, azw1left
      double precision acmprght, acmpleft, acmpbot, acmptop
      double precision ascmprght, ascmpleft, ascmpbot, ascmptop

      dtdx = dt/dx
      dtdy = dt/dy

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in tracexy with ppm_type != 0'
        call bl_error("Error:: Castro_advection_3d.f90 :: tracexy")
      end if

      !!!!!!!!!!!!!!!
      ! NON-PPM CODE
      !!!!!!!!!!!!!!!
      
      ! Compute left and right traced states
      !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth,drho,du,dv,dw,dp,drhoe,alpham,alphap,alpha0r,alpha0e) &
      !$OMP PRIVATE(alpha0v,alpha0w,spminus,spplus,spzero,apright,amright,azrright,azeright,azv1rght,azw1rght,apleft) &
      !$OMP PRIVATE(amleft,azrleft,azeleft,azv1left,azw1left)
      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc**2
            rho = q(i,j,k3d,QRHO)
            u = q(i,j,k3d,QU)
            v = q(i,j,k3d,QV)
            w = q(i,j,k3d,QW)
            p = q(i,j,k3d,QPRES)
            rhoe = q(i,j,k3d,QREINT)
            enth = ( (rhoe+p)/rho )/csq

            drho = dqx(i,j,kc,QRHO)
            du = dqx(i,j,kc,QU)
            dv = dqx(i,j,kc,QV)
            dw = dqx(i,j,kc,QW)
            dp = dqx(i,j,kc,QPRES)
            drhoe = dqx(i,j,kc,QREINT)

            alpham = 0.5d0*(dp/(rho*cc) - du)*rho/cc
            alphap = 0.5d0*(dp/(rho*cc) + du)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0v = dv
            alpha0w = dw

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
            azv1rght = 0.5d0*(-1.d0 - spzero )*alpha0v
            azw1rght = 0.5d0*(-1.d0 - spzero )*alpha0w

            if (i .ge. ilo1) then
               qxp(i,j,kc,QRHO) = rho + apright + amright + azrright
               qxp(i,j,kc,QRHO) = max(small_dens,qxp(i,j,kc,QRHO))
               qxp(i,j,kc,QU) = u + (apright - amright)*cc/rho
               qxp(i,j,kc,QV) = v + azv1rght
               qxp(i,j,kc,QW) = w + azw1rght
               qxp(i,j,kc,QPRES) = p + (apright + amright)*csq
               qxp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
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
            azv1left = 0.5d0*(1.d0 - spzero )*alpha0v
            azw1left = 0.5d0*(1.d0 - spzero )*alpha0w

            if (i .le. ihi1) then
               qxm(i+1,j,kc,QRHO) = rho + apleft + amleft + azrleft
               qxm(i+1,j,kc,QRHO) = max(small_dens, qxm(i+1,j,kc,QRHO))
               qxm(i+1,j,kc,QU) = u + (apleft - amleft)*cc/rho
               qxm(i+1,j,kc,QV) = v + azv1left
               qxm(i+1,j,kc,QW) = w + azw1left
               qxm(i+1,j,kc,QPRES) = p + (apleft + amleft)*csq
               qxm(i+1,j,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
            endif

         enddo
      enddo
      !$OMP END PARALLEL DO

      ! Treat K as a passively advected quantity
      if (QESGS .gt. -1) then
         n = QESGS
         do j = ilo2-1, ihi2+1
            ! Right state
            do i = ilo1, ihi1+1
               u = q(i,j,k3d,QU)
               if (u .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = u*dtdx
               endif
               acmprght = 0.5d0*(-1.d0 - spzero )*dqx(i,j,kc,n)
               qxp(i,j,kc,n) = q(i,j,k3d,n) + acmprght
            enddo
 
            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,QU)
               if (u .ge. 0.d0) then
                  spzero = u*dtdx
               else
                  spzero = 1.d0
               endif
               acmpleft = 0.5d0*(1.d0 - spzero )*dqx(i,j,kc,n)
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + acmpleft
            enddo
         enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,i,j,u,spzero,acmprght,acmpleft) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = QFA + iadv - 1

         do j = ilo2-1, ihi2+1

            ! Right state
            do i = ilo1, ihi1+1
               u = q(i,j,k3d,QU)
               if (u .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = u*dtdx
               endif
               acmprght = 0.5d0*(-1.d0 - spzero )*dqx(i,j,kc,n)
               qxp(i,j,kc,n) = q(i,j,k3d,n) + acmprght
            enddo

            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,QU)
               if (u .ge. 0.d0) then
                  spzero = u*dtdx
               else
                  spzero = 1.d0
               endif
               acmpleft = 0.5d0*(1.d0 - spzero )*dqx(i,j,kc,n)
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + acmpleft
            enddo

         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,i,j,u,spzero,ascmprght,ascmpleft) IF(nspec.gt.1)
      do ispec = 1, nspec
         n = QFS + ispec - 1

         do j = ilo2-1, ihi2+1

            ! Right state
            do i = ilo1, ihi1+1
               u = q(i,j,k3d,QU)
               if (u .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = u*dtdx
               endif
               ascmprght = 0.5d0*(-1.d0 - spzero )*dqx(i,j,kc,n)
               qxp(i,j,kc,n) = q(i,j,k3d,n) + ascmprght
            enddo

            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,QU)
               if (u .ge. 0.d0) then
                  spzero = u*dtdx
               else
                  spzero = 1.d0
               endif
               ascmpleft = 0.5d0*(1.d0 - spzero )*dqx(i,j,kc,n)
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + ascmpleft
            enddo

         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,i,j,u,spzero,ascmprght,ascmpleft) IF(naux.gt.1)
      do iaux = 1, naux
         n = QFX + iaux - 1

         do j = ilo2-1, ihi2+1

            ! Right state
            do i = ilo1, ihi1+1
               u = q(i,j,k3d,QU)
               if (u .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = u*dtdx
               endif
               ascmprght = 0.5d0*(-1.d0 - spzero )*dqx(i,j,kc,n)
               qxp(i,j,kc,n) = q(i,j,k3d,n) + ascmprght
            enddo

            ! Left state
            do i = ilo1-1, ihi1
               u = q(i,j,k3d,QU)
               if (u .ge. 0.d0) then
                  spzero = u*dtdx
               else
                  spzero = 1.d0
               endif
               ascmpleft = 0.5d0*(1.d0 - spzero )*dqx(i,j,kc,n)
               qxm(i+1,j,kc,n) = q(i,j,k3d,n) + ascmpleft
            enddo

         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth,drho,du,dv,dw,dp,drhoe,alpham,alphap,alpha0r) &
      !$OMP PRIVATE(alpha0e,alpha0v,alpha0w,spminus,spplus,spzero,apright,amright,azrright,azeright,azv1rght) &
      !$OMP PRIVATE(azw1rght,apleft,amleft,azrleft,azeleft,azv1left,azw1left) &
      !$OMP PRIVATE(alpha0u,azu1rght,azu1left)
      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc**2
            rho = q(i,j,k3d,QRHO)
            u = q(i,j,k3d,QU)
            v = q(i,j,k3d,QV)
            w = q(i,j,k3d,QW)
            p = q(i,j,k3d,QPRES)
            rhoe = q(i,j,k3d,QREINT)
            enth = ( (rhoe+p)/rho )/csq

            drho = dqy(i,j,kc,QRHO)
            du = dqy(i,j,kc,QU)
            dv = dqy(i,j,kc,QV)
            dw = dqy(i,j,kc,QW)
            dp = dqy(i,j,kc,QPRES)
            drhoe = dqy(i,j,kc,QREINT)

            alpham = 0.5d0*(dp/(rho*cc) - dv)*rho/cc
            alphap = 0.5d0*(dp/(rho*cc) + dv)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0w = dw

            if (v-cc .gt. 0.d0) then
               spminus = -1.d0
            else
               spminus = (v-cc)*dtdy
            endif
            if (v+cc .gt. 0.d0) then
               spplus = -1.d0
            else
               spplus = (v+cc)*dtdy
            endif
            if (v .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = v*dtdy
            endif

            apright = 0.5d0*(-1.d0 - spplus )*alphap
            amright = 0.5d0*(-1.d0 - spminus)*alpham
            azrright = 0.5d0*(-1.d0 - spzero )*alpha0r
            azeright = 0.5d0*(-1.d0 - spzero )*alpha0e
            azu1rght = 0.5d0*(-1.d0 - spzero )*alpha0u
            azw1rght = 0.5d0*(-1.d0 - spzero )*alpha0w

            if (j .ge. ilo2) then
               qyp(i,j,kc,QRHO) = rho + apright + amright + azrright
               qyp(i,j,kc,QRHO) = max(small_dens, qyp(i,j,kc,QRHO))
               qyp(i,j,kc,QV) = v + (apright - amright)*cc/rho
               qyp(i,j,kc,QU) = u + azu1rght
               qyp(i,j,kc,QW) = w + azw1rght
               qyp(i,j,kc,QPRES) = p + (apright + amright)*csq
               qyp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
            end if

            if (v-cc .ge. 0.d0) then
               spminus = (v-cc)*dtdy
            else
               spminus = 1.d0
            endif
            if (v+cc .ge. 0.d0) then
               spplus = (v+cc)*dtdy
            else
               spplus = 1.d0
            endif
            if (v .ge. 0.d0) then
               spzero = v*dtdy
            else
               spzero = 1.d0
            endif

            apleft = 0.5d0*(1.d0 - spplus )*alphap
            amleft = 0.5d0*(1.d0 - spminus)*alpham
            azrleft = 0.5d0*(1.d0 - spzero )*alpha0r
            azeleft = 0.5d0*(1.d0 - spzero )*alpha0e
            azu1left = 0.5d0*(1.d0 - spzero )*alpha0u
            azw1left = 0.5d0*(1.d0 - spzero )*alpha0w

            if (j .le. ihi2) then
               qym(i,j+1,kc,QRHO) = rho + apleft + amleft + azrleft
               qym(i,j+1,kc,QRHO) = max(small_dens, qym(i,j+1,kc,QRHO))
               qym(i,j+1,kc,QV) = v + (apleft - amleft)*cc/rho
               qym(i,j+1,kc,QU) = u + azu1left
               qym(i,j+1,kc,QW) = w + azw1left
               qym(i,j+1,kc,QPRES) = p + (apleft + amleft)*csq
               qym(i,j+1,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
            endif

         enddo
      enddo
      !$OMP END PARALLEL DO

      ! Treat K as a passively advected quantity
      if (QESGS .gt. -1) then
         n = QESGS
         do i = ilo1-1, ihi1+1
 
            ! Top state
            do j = ilo2, ihi2+1
               v = q(i,j,k3d,QV)
               if (v .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = v*dtdy
               endif
               acmptop = 0.5d0*(-1.d0 - spzero )*dqy(i,j,kc,n)
               qyp(i,j,kc,n) = q(i,j,k3d,n) + acmptop
            enddo
 
            ! Bottom state
            do j = ilo2-1, ihi2
               v = q(i,j,k3d,QV)
               if (v .ge. 0.d0) then
                  spzero = v*dtdy
               else
                  spzero = 1.d0
               endif
               acmpbot = 0.5d0*(1.d0 - spzero )*dqy(i,j,kc,n)
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + acmpbot
            enddo
         enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,i,j,v,spzero,acmptop,acmpbot) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = QFA + iadv - 1

         do i = ilo1-1, ihi1+1

            ! Top state
            do j = ilo2, ihi2+1
               v = q(i,j,k3d,QV)
               if (v .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = v*dtdy
               endif
               acmptop = 0.5d0*(-1.d0 - spzero )*dqy(i,j,kc,n)
               qyp(i,j,kc,n) = q(i,j,k3d,n) + acmptop
            enddo

            ! Bottom state
            do j = ilo2-1, ihi2
               v = q(i,j,k3d,QV)
               if (v .ge. 0.d0) then
                  spzero = v*dtdy
               else
                  spzero = 1.d0
               endif
               acmpbot = 0.5d0*(1.d0 - spzero )*dqy(i,j,kc,n)
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + acmpbot
            enddo

         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,i,j,v,spzero,ascmptop,ascmpbot) IF(nspec.gt.1)
      do ispec = 1, nspec
         n = QFS + ispec - 1

         do i = ilo1-1, ihi1+1

            ! Top state
            do j = ilo2, ihi2+1
               v = q(i,j,k3d,QV)
               if (v .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = v*dtdy
               endif
               ascmptop = 0.5d0*(-1.d0 - spzero )*dqy(i,j,kc,n)
               qyp(i,j,kc,n) = q(i,j,k3d,n) + ascmptop
            enddo

            ! Bottom state
            do j = ilo2-1, ihi2
               v = q(i,j,k3d,QV)
               if (v .ge. 0.d0) then
                  spzero = v*dtdy
               else
                  spzero = 1.d0
               endif
               ascmpbot = 0.5d0*(1.d0 - spzero )*dqy(i,j,kc,n)
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + ascmpbot
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,i,j,v,spzero,ascmptop,ascmpbot) IF(naux.gt.1)
      do iaux = 1, naux
         n = QFX + iaux - 1

         do i = ilo1-1, ihi1+1

            ! Top state
            do j = ilo2, ihi2+1
               v = q(i,j,k3d,QV)
               if (v .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = v*dtdy
               endif
               ascmptop = 0.5d0*(-1.d0 - spzero )*dqy(i,j,kc,n)
               qyp(i,j,kc,n) = q(i,j,k3d,n) + ascmptop
            enddo

            ! Bottom state
            do j = ilo2-1, ihi2
               v = q(i,j,k3d,QV)
               if (v .ge. 0.d0) then
                  spzero = v*dtdy
               else
                  spzero = 1.d0
               endif
               ascmpbot = 0.5d0*(1.d0 - spzero )*dqy(i,j,kc,n)
               qym(i,j+1,kc,n) = q(i,j,k3d,n) + ascmpbot
            enddo

         enddo
      enddo
      !$OMP END PARALLEL DO

    end subroutine tracexy

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine tracez(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
           dqz,dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3, &
           qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
           ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QESGS, QPRES, QFA, QFS, QFX, nadv, small_dens, &
                                     ppm_type

      implicit none

      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer dq_l1,dq_l2,dq_l3,dq_h1,dq_h2,dq_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2
      integer km,kc,k3d

      double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

      double precision  dqz(dq_l1:dq_h1,dq_l2:dq_h2,dq_l3:dq_h3,QVAR)
      double precision qzm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision qzp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision dz, dt

      ! Local variables
      integer i, j
      integer n
      integer iadv, ispec, iaux

      double precision dtdz
      double precision cc, csq, rho, u, v, w, p, rhoe

      double precision drho, du, dv, dw, dp, drhoe
      double precision enth, alpham, alphap, alpha0r, alpha0e
      double precision alpha0u, alpha0v
      double precision spminus, spplus, spzero
      double precision apright, amright, azrright, azeright
      double precision azu1rght, azv1rght
      double precision apleft, amleft, azrleft, azeleft
      double precision azu1left, azv1left
      double precision acmpbot, acmptop
      double precision ascmpbot, ascmptop

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in tracez with ppm_type != 0'
        call bl_error("Error:: Castro_advection_3d.f90 :: tracez")
      end if

      dtdz = dt/dz
      
      !!!!!!!!!!!!!!!
      ! NON-PPM CODE
      !!!!!!!!!!!!!!!
      
      !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth,drho,du,dv,dw,dp,drhoe,alpham,alphap,alpha0r,alpha0e) &
      !$OMP PRIVATE(alpha0u,alpha0v,spminus,spplus,spzero,apright,amright,azrright,azeright,azu1rght,azv1rght,apleft) &
      !$OMP PRIVATE(amleft,azrleft,azeleft,azu1left,azv1left)
      do j = ilo2-1, ihi2+1
         do i = ilo1-1, ihi1+1

            cc = c(i,j,k3d)
            csq = cc**2
            rho = q(i,j,k3d,QRHO)
            u = q(i,j,k3d,QU)
            v = q(i,j,k3d,QV)
            w = q(i,j,k3d,QW)
            p = q(i,j,k3d,QPRES)
            rhoe = q(i,j,k3d,QREINT)
            enth = ( (rhoe+p)/rho )/csq

            drho = dqz(i,j,kc,QRHO)
            du = dqz(i,j,kc,QU)
            dv = dqz(i,j,kc,QV)
            dw = dqz(i,j,kc,QW)
            dp = dqz(i,j,kc,QPRES)
            drhoe = dqz(i,j,kc,QREINT)

            alpham = 0.5d0*(dp/(rho*cc) - dw)*rho/cc
            alphap = 0.5d0*(dp/(rho*cc) + dw)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0v = dv

            if (w-cc .gt. 0.d0) then
               spminus = -1.d0
            else
               spminus = (w-cc)*dtdz
            endif
            if (w+cc .gt. 0.d0) then
               spplus = -1.d0
            else
               spplus = (w+cc)*dtdz
            endif
            if (w .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = w*dtdz
            endif

            apright = 0.5d0*(-1.d0 - spplus )*alphap
            amright = 0.5d0*(-1.d0 - spminus)*alpham
            azrright = 0.5d0*(-1.d0 - spzero )*alpha0r
            azeright = 0.5d0*(-1.d0 - spzero )*alpha0e
            azu1rght = 0.5d0*(-1.d0 - spzero )*alpha0u
            azv1rght = 0.5d0*(-1.d0 - spzero )*alpha0v

            qzp(i,j,kc,QRHO) = rho + apright + amright + azrright
            qzp(i,j,kc,QRHO) = max(small_dens, qzp(i,j,kc,QRHO))
            qzp(i,j,kc,QW) = w + (apright - amright)*cc/rho
            qzp(i,j,kc,QU) = u + azu1rght
            qzp(i,j,kc,QV) = v + azv1rght
            qzp(i,j,kc,QPRES) = p + (apright + amright)*csq
            qzp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

            ! repeat above with km (k3d-1) to get qzm at kc
            cc = c(i,j,k3d-1)
            csq = cc**2
            rho = q(i,j,k3d-1,QRHO)
            u = q(i,j,k3d-1,QU)
            v = q(i,j,k3d-1,QV)
            w = q(i,j,k3d-1,QW)
            p = q(i,j,k3d-1,QPRES)
            rhoe = q(i,j,k3d-1,QREINT)
            enth = ( (rhoe+p)/rho )/csq

            drho = dqz(i,j,km,QRHO)
            du = dqz(i,j,km,QU)
            dv = dqz(i,j,km,QV)
            dw = dqz(i,j,km,QW)
            dp = dqz(i,j,km,QPRES)
            drhoe = dqz(i,j,km,QREINT)

            alpham = 0.5d0*(dp/(rho*cc) - dw)*rho/cc
            alphap = 0.5d0*(dp/(rho*cc) + dw)*rho/cc
            alpha0r = drho - dp/csq
            alpha0e = drhoe - dp*enth
            alpha0u = du
            alpha0v = dv

            if (w-cc .ge. 0.d0) then
               spminus = (w-cc)*dtdz
            else
               spminus = 1.d0
            endif
            if (w+cc .ge. 0.d0) then
               spplus = (w+cc)*dtdz
            else
               spplus = 1.d0
            endif
            if (w .ge. 0.d0) then
               spzero = w*dtdz
            else
               spzero = 1.d0
            endif

            apleft = 0.5d0*(1.d0 - spplus )*alphap
            amleft = 0.5d0*(1.d0 - spminus)*alpham
            azrleft = 0.5d0*(1.d0 - spzero )*alpha0r
            azeleft = 0.5d0*(1.d0 - spzero )*alpha0e
            azu1left = 0.5d0*(1.d0 - spzero )*alpha0u
            azv1left = 0.5d0*(1.d0 - spzero )*alpha0v

            qzm(i,j,kc,QRHO) = rho + apleft + amleft + azrleft
            qzm(i,j,kc,QRHO) = max(small_dens, qzm(i,j,kc,QRHO))
            qzm(i,j,kc,QW) = w + (apleft - amleft)*cc/rho
            qzm(i,j,kc,QU) = u + azu1left
            qzm(i,j,kc,QV) = v + azv1left
            qzm(i,j,kc,QPRES) = p + (apleft + amleft)*csq
            qzm(i,j,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

         enddo
      enddo
      !$OMP END PARALLEL DO

      ! Treat K as a passively advected quantity
      if (QESGS .gt. -1) then
         n = QESGS
         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1
 
               ! Top state
               w = q(i,j,k3d,QW)
               if (w .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = w*dtdz
               endif
               acmptop = 0.5d0*(-1.d0 - spzero )*dqz(i,j,kc,n)
               qzp(i,j,kc,n) = q(i,j,k3d,n) + acmptop
 
               ! Bottom state
               w = q(i,j,k3d-1,QW)
               if (w .ge. 0.d0) then
                  spzero = w*dtdz
               else
                  spzero = 1.d0
               endif
               acmpbot = 0.5d0*(1.d0 - spzero )*dqz(i,j,km,n)
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + acmpbot
 
            enddo
         enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iadv,n,i,j,w,acmptop,acmpbot,spzero) IF(nadv.gt.1)
      do iadv = 1, nadv
         n = QFA + iadv - 1

         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1

               ! Top state
               w = q(i,j,k3d,QW)
               if (w .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = w*dtdz
               endif
               acmptop = 0.5d0*(-1.d0 - spzero )*dqz(i,j,kc,n)
               qzp(i,j,kc,n) = q(i,j,k3d,n) + acmptop

               ! Bottom state
               w = q(i,j,k3d-1,QW)
               if (w .ge. 0.d0) then
                  spzero = w*dtdz
               else
                  spzero = 1.d0
               endif
               acmpbot = 0.5d0*(1.d0 - spzero )*dqz(i,j,km,n)
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + acmpbot
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ispec,n,i,j,w,ascmptop,ascmpbot,spzero) IF(nspec.gt.1)
      do ispec = 1, nspec
         n = QFS + ispec - 1

         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1

               ! Top state
               w = q(i,j,k3d,QW)
               if (w .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = w*dtdz
               endif
               ascmptop = 0.5d0*(-1.d0 - spzero )*dqz(i,j,kc,n)
               qzp(i,j,kc,n) = q(i,j,k3d,n) + ascmptop

               ! Bottom state
               w = q(i,j,k3d-1,QW)
               if (w .ge. 0.d0) then
                  spzero = w*dtdz
               else
                  spzero = 1.d0
               endif
               ascmpbot = 0.5d0*(1.d0 - spzero )*dqz(i,j,km,n)
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + ascmpbot
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iaux,n,i,j,w,ascmptop,ascmpbot,spzero) IF(naux.gt.1)
      do iaux = 1, naux
         n = QFX + iaux - 1

         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1

               ! Top state
               w = q(i,j,k3d,QW)
               if (w .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = w*dtdz
               endif
               ascmptop = 0.5d0*(-1.d0 - spzero )*dqz(i,j,kc,n)
               qzp(i,j,kc,n) = q(i,j,k3d,n) + ascmptop

               ! Bottom state
               w = q(i,j,k3d-1,QW)
               if (w .ge. 0.d0) then
                  spzero = w*dtdz
               else
                  spzero = 1.d0
               endif
               ascmpbot = 0.5d0*(1.d0 - spzero )*dqz(i,j,km,n)
               qzm(i,j,kc,n) = q(i,j,k3d-1,n) + ascmpbot
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

    end subroutine tracez

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

    subroutine consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                      uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                      src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                      flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                      flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                      flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                      area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
                      area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
                      area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
                      vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
                      div,pdivu,lo,hi,dx,dy,dz,dt,E_added_flux)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
           UEDEN, UEINT, UTEMP, normalize_species

      implicit none

      integer lo(3), hi(3)
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      integer   src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3 
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
      integer area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3
      integer area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3
      integer area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3
      integer vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3

      double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
      double precision   src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
      double precision flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
      double precision area1(area1_l1:area1_h1,area1_l2:area1_h2,area1_l3:area1_h3)
      double precision area2(area2_l1:area2_h1,area2_l2:area2_h2,area2_l3:area2_h3)
      double precision area3(area3_l1:area3_h1,area3_l2:area3_h2,area3_l3:area3_h3)
      double precision vol(vol_l1:vol_h1,vol_l2:vol_h2,vol_l3:vol_h3)
      double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      double precision pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      double precision dx, dy, dz, dt, E_added_flux

      double precision :: div1
      integer          :: i, j, k, n

      do n = 1, NVAR
         
         if ( n.eq.UTEMP ) then

            flux1(:,:,:,n) = 0.d0
            flux2(:,:,:,n) = 0.d0
            flux3(:,:,:,n) = 0.d0

         else

            !$OMP PARALLEL PRIVATE(i,j,k,div1)
            !$OMP DO
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)+1
                     div1 = .25d0*(div(i,j,k) + div(i,j+1,k) + div(i,j,k+1) + div(i,j+1,k+1))
                     div1 = difmag*min(0.d0,div1)
                     flux1(i,j,k,n) = flux1(i,j,k,n) + dx*div1*(uin(i,j,k,n)-uin(i-1,j,k,n))
                     flux1(i,j,k,n) = flux1(i,j,k,n) * area1(i,j,k) * dt
                  enddo
               enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP DO
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)+1
                  do i = lo(1),hi(1)
                     div1 = .25d0*(div(i,j,k) + div(i+1,j,k) + div(i,j,k+1) + div(i+1,j,k+1))
                     div1 = difmag*min(0.d0,div1)
                     flux2(i,j,k,n) = flux2(i,j,k,n) + dy*div1*(uin(i,j,k,n)-uin(i,j-1,k,n))
                     flux2(i,j,k,n) = flux2(i,j,k,n) * area2(i,j,k) * dt
                  enddo
               enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP DO
            do k = lo(3),hi(3)+1
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     div1 = .25d0*(div(i,j,k) + div(i+1,j,k) + div(i,j+1,k) + div(i+1,j+1,k))
                     div1 = difmag*min(0.d0,div1)
                     flux3(i,j,k,n) = flux3(i,j,k,n) + dz*div1*(uin(i,j,k,n)-uin(i,j,k-1,n))
                     flux3(i,j,k,n) = flux3(i,j,k,n) * area3(i,j,k) * dt
                  enddo
               enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

         endif

      enddo

      if (normalize_species .eq. 1) &
         call normalize_species_fluxes( &
                  flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                  flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                  flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                  lo,hi)

      do n = 1, NVAR

         ! pass temperature through
         if (n .eq. UTEMP) then
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     uout(i,j,k,n) = uin(i,j,k,n)
                  enddo
               enddo
            enddo
         else 
            ! update everything else with fluxes and source terms
            !$OMP PARALLEL DO PRIVATE(i,j,k)
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     uout(i,j,k,n) = uin(i,j,k,n) &
                          + ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) / vol(i,j,k) &
                          +   dt * src(i,j,k,n)
                     !
                     ! Add the source term to (rho e)
                     !
                     if (n .eq. UEINT) then
                        uout(i,j,k,n) = uout(i,j,k,n) - dt * pdivu(i,j,k)
                     else if (n .eq. UEDEN) then
                        E_added_flux = E_added_flux + &
                            ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) / vol(i,j,k) 
                     endif
                  enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         endif

      enddo

      end subroutine consup

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine cmpflx(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
                        ugdnv,pgdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                        gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        idir,ilo,ihi,jlo,jhi,kc,kflux,k3d)

      use meth_params_module, only : QVAR, NVAR

      implicit none

      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
      integer pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3
      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer idir,ilo,ihi,jlo,jhi
      integer i,j,kc,kflux,k3d

      double precision qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,NVAR)
      double precision ugdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
      double precision pgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
      double precision gamc(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      double precision csml(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      double precision    c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

      double precision, allocatable :: smallc(:,:),cavg(:,:)
      double precision, allocatable :: gamcm(:,:),gamcp(:,:)

      allocate ( smallc(ilo-1:ihi+1,jlo-1:jhi+1) )
      allocate (   cavg(ilo-1:ihi+1,jlo-1:jhi+1) )
      allocate (  gamcm(ilo-1:ihi+1,jlo-1:jhi+1) )
      allocate (  gamcp(ilo-1:ihi+1,jlo-1:jhi+1) )
      
      if(idir.eq.1) then
         do j = jlo, jhi
            do i = ilo, ihi
               smallc(i,j) = max( csml(i,j,k3d), csml(i-1,j,k3d) )
               cavg(i,j) = 0.5d0*( c(i,j,k3d) + c(i-1,j,k3d) )
               gamcm(i,j) = gamc(i-1,j,k3d)
               gamcp(i,j) = gamc(i,j,k3d)
            enddo
         enddo
      elseif(idir.eq.2) then
         do j = jlo, jhi
            do i = ilo, ihi
               smallc(i,j) = max( csml(i,j,k3d), csml(i,j-1,k3d) )
               cavg(i,j) = 0.5d0*( c(i,j,k3d) + c(i,j-1,k3d) )
               gamcm(i,j) = gamc(i,j-1,k3d)
               gamcp(i,j) = gamc(i,j,k3d)
            enddo
         enddo
      else
         do j = jlo, jhi
            do i = ilo, ihi
               smallc(i,j) = max( csml(i,j,k3d), csml(i,j,k3d-1) )
               cavg(i,j) = 0.5d0*( c(i,j,k3d) + c(i,j,k3d-1) )
               gamcm(i,j) = gamc(i,j,k3d-1)
               gamcp(i,j) = gamc(i,j,k3d)
            enddo
         enddo
      endif

      ! Solve Riemann problem
      call riemannus(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                     gamcm,gamcp,cavg,smallc,ilo-1,jlo-1,ihi+1,jhi+1, &
                     flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
                     ugdnv,pgdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                     idir,ilo,ihi,jlo,jhi,kc,kflux,k3d)

      deallocate(smallc,cavg,gamcm,gamcp)

      end subroutine cmpflx

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine riemannus(ql,qr,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                           gamcl,gamcr,cav,smallc,gd_l1,gd_l2,gd_h1,gd_h2, &
                           uflx,uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3, &
                           ugdnv,pgdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                           idir,ilo,ihi,jlo,jhi,kc,kflux,k3d)

      use network, only : nspec, naux
      use prob_params_module, only : physbc_lo,Symmetry
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QESGS, QFA, QFS, &
                                     QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, UFX, &
                                     nadv, small_dens, small_pres

      implicit none
      double precision, parameter:: small = 1.d-8
      double precision, parameter:: twothirds = 2.d0/3.d0

      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer gd_l1,gd_l2,gd_h1,gd_h2
      integer uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3
      integer pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3
      integer idir,ilo,ihi,jlo,jhi
      integer i,j,kc,kflux,k3d

      double precision ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision  gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision  gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision    cav(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision smallc(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,uflx_l3:uflx_h3,NVAR)
      double precision ugdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
      double precision pgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)

      integer n, nq
      integer iadv, ispec, iaux

      double precision rgdnv,v1gdnv,v2gdnv,regdnv,ustar
      double precision rl, ul, v1l, v2l, pl, rel
      double precision rr, ur, v1r, v2r, pr, rer
      double precision wl, wr, rhoetot, scr
      double precision rstar, cstar, estar, pstar
      double precision ro, uo, po, reo, co, gamco, entho
      double precision sgnm, spin, spout, ushock, frac
      double precision wsmall, csmall,qavg
      double precision rho_K_contrib

      !$OMP PARALLEL DO PRIVATE(i,j,rl,ul,v1l,v2l,pl,rel,rr,ur,v1r,v2r,pr,rer,csmall,wsmall,wl,wr,pstar,ustar,ro,uo) &
      !$OMP PRIVATE(po,reo,gamco,co,entho,rstar,estar,cstar,sgnm,spout,spin,ushock,scr,frac,v1gdnv,v2gdnv,rgdnv,regdnv) &
      !$OMP PRIVATE(rhoetot,iadv,n,nq,qavg,ispec,iaux)
      do j = jlo, jhi
         do i = ilo, ihi

            rl = max(ql(i,j,kc,QRHO),small_dens)

            ! pick left velocities based on direction
            if(idir.eq.1) then
               ul  = ql(i,j,kc,QU)
               v1l = ql(i,j,kc,QV)
               v2l = ql(i,j,kc,QW)
            elseif(idir.eq.2) then
               ul  = ql(i,j,kc,QV)
               v1l = ql(i,j,kc,QU)
               v2l = ql(i,j,kc,QW)
            else
               ul  = ql(i,j,kc,QW)
               v1l = ql(i,j,kc,QU)
               v2l = ql(i,j,kc,QV)
            endif

            pl  = max(ql(i,j,kc,QPRES ),small_pres)
            rel =     ql(i,j,kc,QREINT)

            rr = max(qr(i,j,kc,QRHO),small_dens)

            ! pick right velocities based on direction
            if(idir.eq.1) then
               ur  = qr(i,j,kc,QU)
               v1r = qr(i,j,kc,QV)
               v2r = qr(i,j,kc,QW)
            elseif(idir.eq.2) then
               ur  = qr(i,j,kc,QV)
               v1r = qr(i,j,kc,QU)
               v2r = qr(i,j,kc,QW)
            else
               ur  = qr(i,j,kc,QW)
               v1r = qr(i,j,kc,QU)
               v2r = qr(i,j,kc,QV)
            endif

            pr  = max(qr(i,j,kc,QPRES),small_pres)
            rer =     qr(i,j,kc,QREINT)

            csmall = smallc(i,j)
            wsmall = small_dens*csmall
            wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
            wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

            pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
            ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)
            pstar = max(pstar,small_pres)

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
               v1gdnv = v1l
               v2gdnv = v2l
            else if (ustar .lt. 0.d0) then
               v1gdnv = v1r
               v2gdnv = v2r
            else
               v1gdnv = 0.5d0*(v1l+v1r)
               v2gdnv = 0.5d0*(v2l+v2r)
            endif
            rgdnv = frac*rstar + (1.d0 - frac)*ro

            ugdnv(i,j,kc) = frac*ustar + (1.d0 - frac)*uo
            pgdnv(i,j,kc) = frac*pstar + (1.d0 - frac)*po

            regdnv = frac*estar + (1.d0 - frac)*reo
            if (spout .lt. 0.d0) then
               rgdnv = ro
               ugdnv(i,j,kc) = uo
               pgdnv(i,j,kc) = po
               regdnv = reo
            endif
            if (spin .ge. 0.d0) then
               rgdnv = rstar
               ugdnv(i,j,kc) = ustar
               pgdnv(i,j,kc) = pstar
               regdnv = estar
            endif

            pgdnv(i,j,kc) = max(pgdnv(i,j,kc),small_pres)

            ! Enforce that fluxes through a symmetry plane are hard zero.
            if (i    .eq.0 .and. physbc_lo(1) .eq. Symmetry .and. idir .eq. 1) &
                 ugdnv(i,j,kc) = 0.d0
            if (j    .eq.0 .and. physbc_lo(2) .eq. Symmetry .and. idir .eq. 2) &
                 ugdnv(i,j,kc) = 0.d0
            if (kflux.eq.0 .and. physbc_lo(3) .eq. Symmetry .and. idir .eq. 3) &
                 ugdnv(i,j,kc) = 0.d0

            ! Compute fluxes, order as conserved state (not q)
            uflx(i,j,kflux,URHO) = rgdnv*ugdnv(i,j,kc)

            if(idir.eq.1) then
               uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
               uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*v1gdnv
               uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*v2gdnv
            elseif(idir.eq.2) then
               uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*v1gdnv
               uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
               uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*v2gdnv
            else
               uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*v1gdnv
               uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*v2gdnv
               uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
            endif

            rhoetot = regdnv + 0.5d0*rgdnv*(ugdnv(i,j,kc)**2 + v1gdnv**2 + v2gdnv**2)

            uflx(i,j,kflux,UEDEN) = ugdnv(i,j,kc)*(rhoetot + pgdnv(i,j,kc))
            uflx(i,j,kflux,UEINT) = ugdnv(i,j,kc)*regdnv

            ! Treat K as a passively advected quantity but allow it to affect fluxes of (rho E) and momenta.
            if (UESGS .gt. -1) then
               n  = UESGS
               nq = QESGS
               if (ustar .gt. 0.d0) then
                  qavg = ql(i,j,kc,nq)
               else if (ustar .lt. 0.d0) then
                  qavg = qr(i,j,kc,nq)
               else
                  qavg = 0.5d0 * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
               endif
    
               uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
 
               rho_K_contrib =  twothirds * rgdnv * qavg
 
               if(idir.eq.1) then
                  uflx(i,j,kflux,UMX) = uflx(i,j,kflux,UMX) + rho_K_contrib
               elseif(idir.eq.2) then
                  uflx(i,j,kflux,UMY) = uflx(i,j,kflux,UMY) + rho_K_contrib
               elseif(idir.eq.3) then
                  uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,UMZ) + rho_K_contrib
               endif
 
               uflx(i,j,kflux,UEDEN) = uflx(i,j,kflux,UEDEN) + ugdnv(i,j,kc) * rho_K_contrib
            end if

            do iadv = 1, nadv
               n  = UFA + iadv - 1
               nq = QFA + iadv - 1
               if (ustar .gt. 0.d0) then
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
               else if (ustar .lt. 0.d0) then
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
               else
                  qavg = 0.5d0 * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
               endif
            enddo

            do ispec = 1, nspec
               n  = UFS + ispec - 1
               nq = QFS + ispec - 1
               if (ustar .gt. 0.d0) then
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
               else if (ustar .lt. 0.d0) then
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
               else
                  qavg = 0.5d0 * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
               endif
            enddo

            do iaux = 1, naux
               n  = UFX + iaux - 1
               nq = QFX + iaux - 1
               if (ustar .gt. 0.d0) then
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
               else if (ustar .lt. 0.d0) then
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
               else
                  qavg = 0.5d0 * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
               endif
            enddo
         
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine riemannus

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

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine uslope(q,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        dqx,dqy,dqz,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        ilo1,ilo2,ihi1,ihi2,kc,k3d,nv)

      use meth_params_module

      implicit none

      integer ilo,ihi
      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2,kc,k3d,nv

      double precision q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,nv)
      double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      double precision dqx(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)
      double precision dqy(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)
      double precision dqz(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)

      integer i, j, k, n

      double precision dlft, drgt, slop, dq1
      double precision dm, dp, dc, ds, sl, dl, dfm, dfp

      double precision, allocatable::dsgn(:,:),dlim(:,:),df(:,:),dcen(:,:)

      double precision, parameter :: four3rd = 4.d0/3.d0, sixth = 1.d0/6.d0

      ilo = MIN(ilo1,ilo2)
      ihi = MAX(ihi1,ihi2)

      allocate (dsgn(ilo-2:ihi+2,ilo-2:ihi+2))
      allocate (dlim(ilo-2:ihi+2,ilo-2:ihi+2))
      allocate (  df(ilo-2:ihi+2,ilo-2:ihi+2))
      allocate (dcen(ilo-2:ihi+2,ilo-2:ihi+2))

      if(iorder.eq.1) then

         do n = 1, nv
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1
                  dqx(i,j,kc,n) = 0.d0
                  dqy(i,j,kc,n) = 0.d0
                  dqz(i,j,kc,n) = 0.d0
               enddo
            enddo
         enddo

      else

         do n = 1, nv 

            ! Compute slopes in first coordinate direction
            !$OMP PARALLEL DO PRIVATE(i,j,dlft,drgt,slop,dq1)
            do j = ilo2-1, ihi2+1

               ! First compute Fromm slopes
               do i = ilo1-2, ihi1+2
                  dlft = 2.d0*(q(i ,j,k3d,n) - q(i-1,j,k3d,n))
                  drgt = 2.d0*(q(i+1,j,k3d,n) - q(i ,j,k3d,n))
                  dcen(i,j) = .25d0 * (dlft+drgt)
                  dsgn(i,j) = sign(1.d0, dcen(i,j))
                  slop = min( abs(dlft), abs(drgt) )
                  if (dlft*drgt .ge. 0.d0) then
                     dlim(i,j) = slop
                  else
                     dlim(i,j) = 0.d0
                  endif
                  df(i,j) = dsgn(i,j)*min( dlim(i,j), abs(dcen(i,j)) )
               enddo

               ! Now compute limited fourth order slopes
               do i = ilo1-1, ihi1+1
                  dq1       = four3rd*dcen(i,j) - sixth*(df(i+1,j) + df(i-1,j))
                  dqx(i,j,kc,n) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dq1))
               enddo

            enddo
            !$OMP END PARALLEL DO

            ! Compute slopes in second coordinate direction
            !$OMP PARALLEL DO PRIVATE(i,j,dlft,drgt,slop,dq1)
            do i = ilo1-1, ihi1+1
               ! First compute Fromm slopes for this column
               do j = ilo2-2, ihi2+2
                  dlft = 2.d0*(q(i,j ,k3d,n) - q(i,j-1,k3d,n))
                  drgt = 2.d0*(q(i,j+1,k3d,n) - q(i,j ,k3d,n))
                  dcen(i,j) = .25d0 * (dlft+drgt)
                  dsgn(i,j) = sign( 1.d0, dcen(i,j) )
                  slop = min( abs(dlft), abs(drgt) )
                  if (dlft*drgt .ge. 0.d0) then
                     dlim(i,j) = slop
                  else
                     dlim(i,j) = 0.d0
                  endif
                  df(i,j) = dsgn(i,j)*min( dlim(i,j),abs(dcen(i,j)) )
               enddo

               ! Now compute limited fourth order slopes
               do j = ilo2-1, ihi2+1
                  dq1 = four3rd*dcen(i,j) - sixth*( df(i,j+1) + df(i,j-1) )
                  dqy(i,j,kc,n) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dq1))
               enddo
            enddo
            !$OMP END PARALLEL DO

            ! Compute slopes in third coordinate direction
            !$OMP PARALLEL DO PRIVATE(i,j,k,dm,dp,dc,ds,sl,dl,dfm,dfp,dq1)
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1

                  ! Compute Fromm slope on slab below
                  k = k3d-1
                  dm = 2.d0*(q(i,j,k ,n) - q(i,j,k-1,n))
                  dp = 2.d0*(q(i,j,k+1,n) - q(i,j,k, n))
                  dc = .25d0*(dm+dp)
                  ds = sign( 1.d0, dc )
                  sl = min( abs(dm), abs(dp) )
                  if (dm*dp .ge. 0.d0) then
                     dl = sl
                  else
                     dl = 0.d0
                  endif
                  dfm = ds*min(dl,abs(dc))

                  ! Compute Fromm slope on slab above
                  k = k3d+1
                  dm = 2.d0*(q(i,j,k ,n) - q(i,j,k-1,n))
                  dp = 2.d0*(q(i,j,k+1,n) - q(i,j,k, n))
                  dc = .25d0*(dm+dp)
                  ds = sign( 1.d0, dc )
                  sl = min( abs(dm), abs(dp) )
                  if (dm*dp .ge. 0.d0) then
                     dl = sl
                  else
                     dl = 0.d0
                  endif
                  dfp = ds*min(dl,abs(dc))

                  ! Compute Fromm slope on current slab
                  k = k3d
                  dm = 2.d0*(q(i,j,k ,n) - q(i,j,k-1,n))
                  dp = 2.d0*(q(i,j,k+1,n) - q(i,j,k, n))
                  dc = .25d0*(dm+dp)
                  ds = sign( 1.d0, dc )
                  sl = min( abs(dm), abs(dp) )
                  if (dm*dp .ge. 0.d0) then
                     dl = sl
                  else
                     dl = 0.d0
                  endif

                  ! Now compute limited fourth order slopes
                  dq1 = four3rd*dc - sixth*( dfp + dfm )
                  dqz(i,j,kc,n) = flatn(i,j,k3d)*ds*min(dl,abs(dq1))
               enddo
            enddo
            !$OMP END PARALLEL DO
         enddo

      endif

      deallocate(dsgn,dlim,df,dcen)

      end subroutine uslope

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine pslope(p,rho,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        dpx,dpy,dpz,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                        ilo1,ilo2,ihi1,ihi2,kc,k3d,dx,dy,dz)
        
        use meth_params_module

        implicit none

        integer ilo,ihi
        integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
        integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
        integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
        integer ilo1,ilo2,ihi1,ihi2,kc,k3d

        double precision p  (qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
        double precision rho(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
        double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
        double precision dpx(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3)
        double precision dpy(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3)
        double precision dpz(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3)
        double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)
        double precision dx,dy,dz

        integer i, j, k

        double precision dlft, drgt, dp1
        double precision dm, dp, dc, dl, dfm, dfp, ds

        double precision, parameter :: four3rd = 4.d0/3.d0, sixth = 1.d0/6.d0

        !     Local arrays
        double precision, allocatable::dsgn(:,:),dlim(:,:),df(:,:),dcen(:,:)

        ilo = MIN(ilo1,ilo2)
        ihi = MAX(ihi1,ihi2)

        allocate (dsgn(ilo-2:ihi+2,ilo-2:ihi+2))
        allocate (dlim(ilo-2:ihi+2,ilo-2:ihi+2))
        allocate (  df(ilo-2:ihi+2,ilo-2:ihi+2))
        allocate (dcen(ilo-2:ihi+2,ilo-2:ihi+2))

        if(iorder.eq.1) then

           do j = ilo2-1, ihi2+1
              do i = ilo1-1, ihi1+1
                 dpx(i,j,kc) = 0.d0
                 dpy(i,j,kc) = 0.d0
                 dpz(i,j,kc) = 0.d0
              enddo
           enddo

        else
           ! Compute slopes in first coordinate direction
           !$OMP PARALLEL DO PRIVATE(i,j,dlft,drgt,dp1)
           do j = ilo2-1, ihi2+1

              ! First compute Fromm slopes
              do i = ilo1-2, ihi1+2

                 dlft = p(i  ,j,k3d) - p(i-1,j,k3d)
                 drgt = p(i+1,j,k3d) - p(i  ,j,k3d)

                 ! Subtract off (rho * grav) so as not to limit that part of the slope
                 dlft = dlft - 0.25d0 * &
                      (rho(i,j,k3d)+rho(i-1,j,k3d))*(grav(i,j,k3d,1)+grav(i-1,j,k3d,1))*dx
                 drgt = drgt - 0.25d0 * &
                      (rho(i,j,k3d)+rho(i+1,j,k3d))*(grav(i,j,k3d,1)+grav(i+1,j,k3d,1))*dx

                 dcen(i,j) = 0.5d0*(dlft+drgt)
                 dsgn(i,j) = sign(1.d0, dcen(i,j))
                 if (dlft*drgt .ge. 0.d0) then
                    dlim(i,j) = 2.d0 * min( abs(dlft), abs(drgt) )
                 else
                    dlim(i,j) = 0.d0
                 endif
                 df(i,j) = dsgn(i,j)*min( dlim(i,j), abs(dcen(i,j)) )
              enddo

              ! Now limited fourth order slopes
              do i = ilo1-1, ihi1+1
                 dp1         = four3rd*dcen(i,j) - sixth*(df(i+1,j) + df(i-1,j))
                 dpx(i,j,kc) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dp1))
                 dpx(i,j,kc) = dpx(i,j,kc) + rho(i,j,k3d)*grav(i,j,k3d,1)*dx
              enddo
           enddo
           !$OMP END PARALLEL DO

           ! Compute slopes in second coordinate direction
           !$OMP PARALLEL DO PRIVATE(i,j,dlft,drgt,dp1)
           do i = ilo1-1, ihi1+1

              ! First compute Fromm slopes
              do j = ilo2-2, ihi2+2
                 dlft = p(i,j  ,k3d) - p(i,j-1,k3d)
                 drgt = p(i,j+1,k3d) - p(i,j  ,k3d)

                 ! Subtract off (rho * grav) so as not to limit that part of the slope
                 dlft = dlft - 0.25d0 * &
                      (rho(i,j,k3d)+rho(i,j-1,k3d))*(grav(i,j,k3d,2)+grav(i,j-1,k3d,2))*dy
                 drgt = drgt - 0.25d0 * &
                      (rho(i,j,k3d)+rho(i,j+1,k3d))*(grav(i,j,k3d,2)+grav(i,j+1,k3d,2))*dy

                 dcen(i,j) = 0.5d0*(dlft+drgt)
                 dsgn(i,j) = sign( 1.d0, dcen(i,j) )
                 if (dlft*drgt .ge. 0.d0) then
                    dlim(i,j) = 2.d0 * min( abs(dlft), abs(drgt) )
                 else
                    dlim(i,j) = 0.d0
                 endif
                 df(i,j) = dsgn(i,j)*min( dlim(i,j),abs(dcen(i,j)) )
              enddo

              ! Now limited fourth order slopes
              do j = ilo2-1, ihi2+1
                 dp1 = four3rd*dcen(i,j) - sixth*( df(i,j+1) + df(i,j-1) )
                 dpy(i,j,kc) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dp1))
                 dpy(i,j,kc) = dpy(i,j,kc) + rho(i,j,k3d)*grav(i,j,k3d,2)*dy
              enddo
           enddo
           !$OMP END PARALLEL DO

           ! Compute slopes in third coordinate direction
           !$OMP PARALLEL DO PRIVATE(i,j,k,dm,dp,dc,ds,dl,dfm,dfp,dp1)
           do j = ilo2-1, ihi2+1
              do i = ilo1-1, ihi1+1

                 ! compute Fromm slopes on slab below
                 k = k3d-1
                 dm = p(i,j,k  ) - p(i,j,k-1)
                 dp = p(i,j,k+1) - p(i,j,k  )
                 dm = dm - 0.25d0 * (rho(i,j,k)+rho(i,j,k-1))* &
                      (grav(i,j,k,3)+grav(i,j,k-1,3))*dz
                 dp = dp - 0.25d0 * (rho(i,j,k)+rho(i,j,k+1))* &
                      (grav(i,j,k,3)+grav(i,j,k+1,3))*dz
                 dc = 0.5d0*(dm+dp)
                 ds = sign( 1.d0, dc )
                 if (dm*dp .ge. 0.d0) then
                    dl = 2.d0 * min( abs(dm), abs(dp) )
                 else
                    dl = 0.d0
                 endif
                 dfm = ds*min(dl,abs(dc))

                 ! compute Fromm slopes on slab above
                 k = k3d+1
                 dm = p(i,j,k  ) - p(i,j,k-1)
                 dp = p(i,j,k+1) - p(i,j,k  )
                 dm = dm - 0.25d0 * (rho(i,j,k)+rho(i,j,k-1))* &
                      (grav(i,j,k,3)+grav(i,j,k-1,3))*dz
                 dp = dp - 0.25d0 * (rho(i,j,k)+rho(i,j,k+1))* &
                      (grav(i,j,k,3)+grav(i,j,k+1,3))*dz
                 dc = 0.5d0*(dm+dp)
                 ds = sign( 1.d0, dc )
                 if (dm*dp .ge. 0.d0) then
                    dl = 2.d0 * min( abs(dm), abs(dp) )
                 else
                    dl = 0.d0
                 endif
                 dfp = ds*min(dl,abs(dc))

                 ! compute Fromm slopes on current slab
                 k = k3d
                 dm = p(i,j,k  ) - p(i,j,k-1)
                 dp = p(i,j,k+1) - p(i,j,k  )
                 dm = dm - 0.25d0 * (rho(i,j,k)+rho(i,j,k-1))* &
                      (grav(i,j,k,3)+grav(i,j,k-1,3))*dz
                 dp = dp - 0.25d0 * (rho(i,j,k)+rho(i,j,k+1))* &
                      (grav(i,j,k,3)+grav(i,j,k+1,3))*dz
                 dc = 0.5d0*(dm+dp)
                 ds = sign( 1.d0, dc )
                 if (dm*dp .ge. 0.d0) then
                    dl = 2.d0 * min( abs(dm), abs(dp) )
                 else
                    dl = 0.d0
                 endif

                 ! now limited fourth order slopes
                 dp1 = four3rd*dc - sixth*( dfp + dfm )
                 dpz(i,j,kc) = flatn(i,j,k3d)*ds*min(dl,abs(dp1))
                 dpz(i,j,kc) = dpz(i,j,kc) + rho(i,j,k3d)*grav(i,j,k3d,3)*dz
              enddo
           enddo
           !$OMP END PARALLEL DO

        endif

        deallocate(dsgn,dlim,df,dcen)

      end subroutine pslope

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine uflaten(lo,hi,p,u,v,w,flatn,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3)

      use meth_params_module, only : iorder, small_pres

      implicit none

      integer lo(3),hi(3)
      integer q_l1,q_l2,q_l3,q_h1,q_h2,q_h3

      double precision p(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      double precision u(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      double precision v(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      double precision w(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      double precision flatn(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)

      integer i, j, k, idx, ishft
      integer nx,ny,nz,nmax

      double precision denom, zeta, tst, tmp, ftmp

      ! Local arrays
      double precision, allocatable :: dp(:,:,:), z(:,:,:), chi(:,:,:)

      ! Knobs for detection of strong shock

      double precision, parameter :: shktst = 0.33d0, zcut1 = 0.75d0, zcut2 = 0.85d0, dzcut = 1.d0/(zcut2-zcut1)

      nx = hi(1)-lo(1)+3
      ny = hi(2)-lo(2)+3
      nz = hi(3)-lo(3)+3

      nmax = max(nx,ny,nz)

      if (iorder .eq. 3) then
         do k = lo(3),hi(3)
            do j = lo(2),hi(2) 
               do i = lo(1),hi(1) 
                  flatn(i,j,k) = 1.d0
               enddo
            enddo
         enddo
         return
      endif

      ! x-direction flattening coef
      allocate(dp (0:nmax-1,lo(2):hi(2),lo(3):hi(3)))
      allocate(z  (0:nmax-1,lo(2):hi(2),lo(3):hi(3)))
      allocate(chi(0:nmax-1,lo(2):hi(2),lo(3):hi(3)))
      !$OMP PARALLEL DO PRIVATE(i,j,k,idx,denom,zeta,tst,tmp,ishft)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2) 
            do i = lo(1)-1,hi(1)+1
               idx = i-lo(1)+1
               dp(idx,j,k) = p(i+1,j,k) - p(i-1,j,k)
               denom = max(small_pres,abs(p(i+2,j,k)-p(i-2,j,k)))
               zeta = abs(dp(idx,j,k))/denom
               z(idx,j,k) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )
               if (u(i-1,j,k)-u(i+1,j,k) .ge. 0.d0) then
                  tst = 1.d0
               else
                  tst = 0.d0
               endif
               tmp = min(p(i+1,j,k),p(i-1,j,k))
               if ((abs(dp(idx,j,k))/tmp).gt.shktst) then
                  chi(idx,j,k) = tst
               else
                  chi(idx,j,k) = 0.d0
               endif
            enddo
            do i = lo(1),hi(1)
               idx = i-lo(1)+1
               if(dp(idx,j,k).gt.0.d0)then
                  ishft = 1
               else
                  ishft = -1
               endif
               flatn(i,j,k) = 1.d0 - &
                    max(chi(idx-ishft,j,k)*z(idx-ishft,j,k),chi(idx,j,k)*z(idx,j,k))
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      deallocate(dp,z,chi)

      ! y-direction flattening coef
      allocate(dp (lo(1):hi(1),0:nmax-1,lo(3):hi(3)))
      allocate(z  (lo(1):hi(1),0:nmax-1,lo(3):hi(3)))
      allocate(chi(lo(1):hi(1),0:nmax-1,lo(3):hi(3)))
      !$OMP PARALLEL DO PRIVATE(i,j,k,idx,denom,zeta,tst,tmp,ishft,ftmp)
      do k = lo(3),hi(3)
         do i = lo(1),hi(1)
            do j = lo(2)-1,hi(2)+1
               idx = j-lo(2)+1
               dp(i,idx,k) = p(i,j+1,k) - p(i,j-1,k)
               denom = max(small_pres,abs(p(i,j+2,k)-p(i,j-2,k)))
               zeta = abs(dp(i,idx,k))/denom
               z(i,idx,k) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )
               if (v(i,j-1,k)-v(i,j+1,k) .ge. 0.d0) then
                  tst = 1.d0
               else
                  tst = 0.d0
               endif
               tmp = min(p(i,j+1,k),p(i,j-1,k))
               if ((abs(dp(i,idx,k))/tmp).gt.shktst) then
                  chi(i,idx,k) = tst
               else
                  chi(i,idx,k) = 0.d0
               endif
            enddo
            do j = lo(2),hi(2)
               idx = j-lo(2)+1
               if(dp(i,idx,k).gt.0.d0)then
                  ishft = 1
               else
                  ishft = -1
               endif
               ftmp = 1.d0 - &
                    max(chi(i,idx-ishft,k)*z(i,idx-ishft,k),chi(i,idx,k)*z(i,idx,k))
               flatn(i,j,k) = min( flatn(i,j,k), ftmp )
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      deallocate(dp,z,chi)

      ! z-direction flattening coef
      allocate(dp (lo(1):hi(1),lo(2):hi(2),0:nmax-1))
      allocate(z  (lo(1):hi(1),lo(2):hi(2),0:nmax-1))
      allocate(chi(lo(1):hi(1),lo(2):hi(2),0:nmax-1))
      !$OMP PARALLEL DO PRIVATE(i,j,k,idx,denom,zeta,tst,tmp,ishft,ftmp)
      do j = lo(2),hi(2) 
         do i = lo(1),hi(1)
            do k = lo(3)-1,hi(3)+1
               idx = k-lo(3)+1
               dp(i,j,idx) = p(i,j,k+1) - p(i,j,k-1)
               denom = max(small_pres,abs(p(i,j,k+2)-p(i,j,k-2)))
               zeta = abs(dp(i,j,idx))/denom
               z(i,j,idx) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )
               if (w(i,j,k-1)-w(i,j,k+1) .ge. 0.d0) then
                  tst = 1.d0
               else
                  tst = 0.d0
               endif
               tmp = min(p(i,j,k+1),p(i,j,k-1))
               if ((abs(dp(i,j,idx))/tmp).gt.shktst) then
                  chi(i,j,idx) = tst
               else
                  chi(i,j,idx) = 0.d0
               endif
            enddo
            do k = lo(3),hi(3)
               idx = k-lo(3)+1
               if(dp(i,j,idx).gt.0.d0)then
                  ishft = 1
               else
                  ishft = -1
               endif
               ftmp = 1.d0 - &
                    max(chi(i,j,idx-ishft)*z(i,j,idx-ishft),chi(i,j,idx)*z(i,j,idx))
               flatn(i,j,k) = min( flatn(i,j,k), ftmp )
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      deallocate(dp,z,chi)

      end subroutine uflaten

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine divu(lo,hi,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,dx,dy,dz, &
                      div,div_l1,div_l2,div_l3,div_h1,div_h2,div_h3)

      use meth_params_module, only : QU, QV, QW

      implicit none

      integer          :: lo(3),hi(3)
      integer          :: q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
      integer          :: div_l1,div_l2,div_l3,div_h1,div_h2,div_h3
      double precision :: dx, dy, dz
      double precision :: div(div_l1:div_h1,div_l2:div_h2,div_l3:div_h3)
      double precision :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,*)

      integer          :: i, j, k
      double precision :: ux, vy, wz

      !$OMP PARALLEL DO PRIVATE(i,j,k,ux,vy,wz)
      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)+1
               
               ux = .25d0*( &
                    + q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                    + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                    + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                    + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) )/ dx

               vy = .25d0*( &
                    + q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                    + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                    + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                    + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) )/ dy

               wz = .25d0*( &
                    + q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                    + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                    + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                    + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) )/ dz

               div(i,j,k) = ux + vy + wz

            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine divu

