module rad_advection_module

  use bl_constants_module

  implicit none

  private

  public umeth3d_rad, ctoprim_rad, consup_rad

contains

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine ctoprim_rad(lo,hi, &
                         uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                         Erin,Erin_l1,Erin_l2,Erin_l3,Erin_h1,Erin_h2,Erin_h3, &
                         lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
                         q,c,cg,gamc,gamcg,csml,flatn,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
                         src, src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                         srcQ, srQ_l1,srQ_l2,srQ_l3,srQ_h1,srQ_h2,srQ_h3, &
                         courno,dx,dy,dz,dt,ngp,ngf,iflaten)

    ! Will give primitive variables on lo-ngp:hi+ngp, and flatn on
    ! lo-ngf:hi+ngf if iflaten=1.  Declared dimensions of
    ! q,c,gamc,csml,flatn are given by DIMS(q).  This declared region
    ! is assumed to encompass lo-ngp:hi+ngp.  Also, uflaten call
    ! assumes ngp>=ngf+3 (ie, primitve data is used by the routine
    ! that computes flatn).

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, UFA, UFS, UFX, &
                                   QVAR, QRHO, QU, QV, QW, QGAME, &
                                   QREINT, QPRES, QTEMP, QFA, QFS, QFX, &
                                   nadv, allow_negative_energy, small_temp  
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, comoving, &
                                       flatten_pp_threshold, first_order_hydro
    use rad_params_module, only : ngroups
    use flatten_module, only : uflaten
    use fluxlimiter_module, only : Edd_factor

    implicit none

    double precision, parameter:: small = 1.d-8

    integer lo(3), hi(3)
    integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
    integer Erin_l1,Erin_l2,Erin_l3,Erin_h1,Erin_h2,Erin_h3
    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer srQ_l1,srQ_l2,srQ_l3,srQ_h1,srQ_h2,srQ_h3

    double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
    double precision Erin(Erin_l1:Erin_h1,Erin_l2:Erin_h2,Erin_l3:Erin_h3,0:ngroups-1)
    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)
    double precision ::     q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QRADVAR)
    double precision ::     c(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision ::    cg(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision :: gamc (q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision :: gamcg(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision :: csml (q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision :: flatn(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision ::  src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
    double precision :: srcQ(srQ_l1:srQ_h1,srQ_l2:srQ_h2,srQ_l3:srQ_h3,QVAR)
    double precision :: dx, dy, dz, dt, courno
    integer          :: ngp, ngf, iflaten

    ! Local variables

    double precision, allocatable:: dpdrho(:,:,:)
    double precision, allocatable:: dpde(:,:,:)
    double precision, allocatable:: flatg(:,:,:)

    integer          :: i, j, k, g
    integer          :: loq(3), hiq(3)
    integer          :: n, nq
    integer          :: iadv, ispec, iaux
    double precision :: courx, coury, courz, courmx, courmy, courmz

    double precision :: csrad2, prad, Eddf, gamr

    type(eos_t) :: eos_state

    allocate( dpdrho(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3))
    allocate(   dpde(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3))
    allocate(  flatg(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3))

    do i=1,3
       loq(i) = lo(i)-ngp
       hiq(i) = hi(i)+ngp
    enddo

    ! Make q (all but p), except put e in slot for rho.e, fix after eos call.
    ! The temperature is used as an initial guess for the eos call and will be overwritten.
    !
    do k = loq(3),hiq(3)
       do j = loq(2),hiq(2)
          do i = loq(1),hiq(1)

             if (uin(i,j,k,URHO) .le. 0.d0) then
                print *,'   '
                print *,'>>> Error: Castro_3d::ctoprim ',i,j,k
                print *,'>>> ... negative density ',uin(i,j,k,URHO)
                call bl_error("Error:: Castro_3d.f90 :: ctoprim")
             end if

             q(i,j,k,QRHO) = uin(i,j,k,URHO)
             q(i,j,k,QU) = uin(i,j,k,UMX)/uin(i,j,k,URHO)
             q(i,j,k,QV) = uin(i,j,k,UMY)/uin(i,j,k,URHO)
             q(i,j,k,QW) = uin(i,j,k,UMZ)/uin(i,j,k,URHO)
             ! convert "rho e" to "e"
             q(i,j,k,QREINT ) = uin(i,j,k,UEINT)/q(i,j,k,QRHO)
             q(i,j,k,QTEMP  ) = uin(i,j,k,UTEMP)
             q(i,j,k,qrad:qradhi) = Erin(i,j,k,:)
          enddo
       enddo
    enddo

    ! Load advected quatities, c, into q, assuming they arrived in uin as rho.c
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
    end do

    ! Load chemical species, c, into q, assuming they arrived in uin as rho.c
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

    ! Load auxiliary variables which are needed in the EOS
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

    ! Get gamc, p, T, c, csml using q state
    do k = loq(3), hiq(3)
       do j = loq(2), hiq(2)
          do i = loq(1), hiq(1)

             eos_state % rho = q(i,j,k,QRHO)
             eos_state % T   = q(i,j,k,QTEMP)
             eos_state % e   = q(i,j,k,QREINT)
             eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)
             eos_state % aux = q(i,j,k,QFX:QFX+naux-1)

             ! If necessary, reset the energy using small_temp
             if ((allow_negative_energy .eq. 0) .and. (q(i,j,k,QREINT) .lt. 0)) then
                q(i,j,k,QTEMP) = small_temp
                eos_state % T = q(i,j,k,QTEMP)
                call eos(eos_input_rt, eos_state)
                q(i,j,k,QPRES ) = eos_state % p
                q(i,j,k,QREINT) = eos_state % e
                if (q(i,j,k,QREINT) .lt. 0.d0) then
                   print *,'   '
                   print *,'>>> Error: ctoprim ',i,j,k
                   print *,'>>> ... new e from eos call is negative ' &
                        ,q(i,j,k,QREINT)
                   print *,'    '
                   call bl_error("Error:: ctoprim")
                end if
             end if

             call eos(eos_input_re, eos_state)

             q(i,j,k,QTEMP) = eos_state % T
             q(i,j,k,QPRES) = eos_state % p
             dpdrho(i,j,k)  = eos_state % dpdr_e
             dpde(i,j,k)    = eos_state % dpde
             gamcg(i,j,k)   = eos_state % gam1
             cg(i,j,k)      = eos_state % cs

             csrad2 = 0.d0
             prad = 0.d0
             do g=0, ngroups-1
                if (comoving) then
                   Eddf = Edd_factor(lam(i,j,k,g))
                   gamr = (3.d0-Eddf)/2.d0
                else
                   gamr = lam(i,j,k,g) + 1.d0
                end if
                prad = prad + lam(i,j,k,g)*q(i,j,k,qrad+g)
                csrad2 = csrad2 + gamr * (lam(i,j,k,g)*q(i,j,k,qrad+g)) / q(i,j,k,QRHO)
             end do

             q(i,j,k,qptot) = q(i,j,k,QPRES) + prad
             c(i,j,k) = cg(i,j,k)**2 + csrad2
             gamc(i,j,k) = c(i,j,k) * q(i,j,k,QRHO) / q(i,j,k,qptot)
             c(i,j,k) = sqrt(c(i,j,k))
             csml(i,j,k) = max(small, small * c(i,j,k))

             ! convert "e" back to "rho e"
             q(i,j,k,QREINT) = q(i,j,k,QREINT)*q(i,j,k,QRHO)
             q(i,j,k,qreitot) = q(i,j,k,QREINT) + sum(q(i,j,k,qrad:qradhi))

          end do
       end do
    end do

    ! compute srcQ terms
    do k = loq(3), hiq(3)
       do j = loq(2), hiq(2)
          do i = loq(1), hiq(1)

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
             !    sum(dpdX_er(i,j,k,:)*(src(i,j,k,UFS:UFS+nspec-1) - &
             !    q(i,j,k,QFS:QFS+nspec-1)*srcQ(i,j,k,QRHO))) &
             !    /q(i,j,k,QRHO)

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

    ! Compute running max of Courant number over grids
    courmx = courno
    courmy = courno
    courmz = courno
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
                call bl_warning("Warning:: Castro_3d.f90 :: CFL violation in ctoprim")
                print *,'>>> ... (u+c) * dt / dx > 1 ', courx
                print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                print *,'>>> ... u, c                ',q(i,j,k,QU), c(i,j,k)
                print *,'>>> ... density             ',q(i,j,k,QRHO)
             end if

             if (coury .gt. 1.d0) then
                print *,'   '
                call bl_warning("Warning:: Castro_3d.f90 :: CFL violation in ctoprim")
                print *,'>>> ... (v+c) * dt / dx > 1 ', coury
                print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                print *,'>>> ... v, c                ',q(i,j,k,QV), c(i,j,k)
                print *,'>>> ... density             ',q(i,j,k,QRHO)
             end if

             if (courz .gt. 1.d0) then
                print *,'   '
                call bl_warning("Warning:: Castro_3d.f90 :: CFL violation in ctoprim")
                print *,'>>> ... (w+c) * dt / dx > 1 ', courz
                print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                print *,'>>> ... w, c                ',q(i,j,k,QW), c(i,j,k)
                print *,'>>> ... density             ',q(i,j,k,QRHO)
             end if

          enddo
       enddo
    enddo

    courno = max( courmx, courmy, courmz )

    ! Compute flattening coef for slope calculations
    if (first_order_hydro) then
       flatn = 0.d0
    else if(iflaten.eq.1)then
       do n=1,3
          loq(n)=lo(n)-ngf
          hiq(n)=hi(n)+ngf
       enddo
       call uflaten(loq,hiq, &
            q(:,:,:,qptot), &
            q(:,:,:,QU), &
            q(:,:,:,QV), &
            q(:,:,:,QW), &
            flatn,(/q_l1,q_l2,q_l3/),(/q_h1,q_h2,q_h3/))
       call uflaten(loq,hiq, &
            q(:,:,:,qpres), &
            q(:,:,:,QU), &
            q(:,:,:,QV), &
            q(:,:,:,QW), &
            flatg,(/q_l1,q_l2,q_l3/),(/q_h1,q_h2,q_h3/))
       flatn = flatn * flatg

       if (flatten_pp_threshold > 0.d0) then
          call ppflaten(loq,hiq, &
               flatn, q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3)
       end if
    else
       flatn = 1.d0
    endif

    deallocate(dpdrho,dpde)

    q(:,:,:,QGAME) = 0.d0 ! QGAME is not used in radiation hydro. Setting it to 0 to mute valgrind.

  end subroutine ctoprim_rad


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

  subroutine umeth3d_rad(q, c,cg, gamc,gamcg, csml, flatn, qd_lo, qd_hi, &
                         lam,lam_lo,lam_hi, &
                         srcQ, src_lo, src_hi, &
                         lo, hi, dx, dy, dz, dt, &
                         flux1, fd1_lo, fd1_hi, &
                         flux2, fd2_lo, fd2_hi, &
                         flux3, fd3_lo, fd3_hi, &
                         rflux1,rfd1_lo, rfd1_hi, &
                         rflux2,rfd2_lo, rfd2_hi, &
                         rflux3,rfd3_lo, rfd3_hi, &
                         ugdnvx_out, ergdx_out, lmgdx_out, ugdnvx_lo, ugdnvx_hi, &
                         ugdnvy_out, ergdy_out, lmgdy_out, ugdnvy_lo, ugdnvy_hi, &
                         ugdnvz_out, ergdz_out, lmgdz_out, ugdnvz_lo, ugdnvz_hi, &
                         pdivu, uy_xfc, uz_xfc, ux_yfc, uz_yfc, ux_zfc, uy_zfc, domlo, domhi)

    ! TODO:
    ! dx needs to be a vector

    use meth_params_module, only : QVAR, NVAR, QU, ppm_type, hybrid_riemann, &
                                   GDPRES, GDU, GDV, GDW, GDERADS, GDLAMS, ngdnv
    use ppm_module
    use radhydro_params_module, only : QRADVAR
    use rad_params_module, only : ngroups
    use riemann_module, only : cmpflx, shock
    use trace_ppm_rad_module, only : tracexy_ppm_rad, tracez_ppm_rad
    use transverse_rad_module
    use mempool_module, only : bl_allocate, bl_deallocate

    implicit none

    integer :: lam_lo(3), lam_hi(3)
    integer :: src_lo(3), src_hi(3)
    integer :: rfd1_lo(3), rfd1_hi(3)
    integer :: rfd2_lo(3), rfd2_hi(3)
    integer :: rfd3_lo(3), rfd3_hi(3)
    integer :: qd_lo(3), qd_hi(3)
    integer :: lo(3), hi(3)

    integer :: fd1_lo(3), fd1_hi(3)
    integer :: fd2_lo(3), fd2_hi(3)
    integer :: fd3_lo(3), fd3_hi(3)
    integer :: ugdnvx_lo(3), ugdnvx_hi(3)
    integer :: ugdnvy_lo(3), ugdnvy_hi(3)
    integer :: ugdnvz_lo(3), ugdnvz_hi(3)

    double precision     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR)
    double precision     c(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision    cg(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision  gamc(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision gamcg(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision  csml(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    double precision flux1(fd1_lo(1):fd1_hi(1),fd1_lo(2):fd1_hi(2),fd1_lo(3):fd1_hi(3),NVAR)
    double precision flux2(fd2_lo(1):fd2_hi(1),fd2_lo(2):fd2_hi(2),fd2_lo(3):fd2_hi(3),NVAR)
    double precision flux3(fd3_lo(1):fd3_hi(1),fd3_lo(2):fd3_hi(2),fd3_lo(3):fd3_hi(3),NVAR)
    double precision ugdnvx_out(ugdnvx_lo(1):ugdnvx_hi(1),ugdnvx_lo(2):ugdnvx_hi(2),ugdnvx_lo(3):ugdnvx_hi(3))
    double precision ugdnvy_out(ugdnvy_lo(1):ugdnvy_hi(1),ugdnvy_lo(2):ugdnvy_hi(2),ugdnvy_lo(3):ugdnvy_hi(3))
    double precision ugdnvz_out(ugdnvz_lo(1):ugdnvz_hi(1),ugdnvz_lo(2):ugdnvz_hi(2),ugdnvz_lo(3):ugdnvz_hi(3))
    double precision pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision uy_xfc(ugdnvx_lo(1):ugdnvx_hi(1),ugdnvx_lo(2):ugdnvx_hi(2),ugdnvx_lo(3):ugdnvx_hi(3))
    double precision uz_xfc(ugdnvx_lo(1):ugdnvx_hi(1),ugdnvx_lo(2):ugdnvx_hi(2),ugdnvx_lo(3):ugdnvx_hi(3))
    double precision ux_yfc(ugdnvy_lo(1):ugdnvy_hi(1),ugdnvy_lo(2):ugdnvy_hi(2),ugdnvy_lo(3):ugdnvy_hi(3))
    double precision uz_yfc(ugdnvy_lo(1):ugdnvy_hi(1),ugdnvy_lo(2):ugdnvy_hi(2),ugdnvy_lo(3):ugdnvy_hi(3))
    double precision ux_zfc(ugdnvz_lo(1):ugdnvz_hi(1),ugdnvz_lo(2):ugdnvz_hi(2),ugdnvz_lo(3):ugdnvz_hi(3))
    double precision uy_zfc(ugdnvz_lo(1):ugdnvz_hi(1),ugdnvz_lo(2):ugdnvz_hi(2),ugdnvz_lo(3):ugdnvz_hi(3))
    double precision dx, dy, dz, dt

    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision ergdx_out(ugdnvx_lo(1):ugdnvx_hi(1),ugdnvx_lo(2):ugdnvx_hi(2),ugdnvx_lo(3):ugdnvx_hi(3),0:ngroups-1)
    double precision ergdy_out(ugdnvy_lo(1):ugdnvy_hi(1),ugdnvy_lo(2):ugdnvy_hi(2),ugdnvy_lo(3):ugdnvy_hi(3),0:ngroups-1)
    double precision ergdz_out(ugdnvz_lo(1):ugdnvz_hi(1),ugdnvz_lo(2):ugdnvz_hi(2),ugdnvz_lo(3):ugdnvz_hi(3),0:ngroups-1)
    double precision lmgdx_out(ugdnvx_lo(1):ugdnvx_hi(1),ugdnvx_lo(2):ugdnvx_hi(2),ugdnvx_lo(3):ugdnvx_hi(3),0:ngroups-1)
    double precision lmgdy_out(ugdnvy_lo(1):ugdnvy_hi(1),ugdnvy_lo(2):ugdnvy_hi(2),ugdnvy_lo(3):ugdnvy_hi(3),0:ngroups-1)
    double precision lmgdz_out(ugdnvz_lo(1):ugdnvz_hi(1),ugdnvz_lo(2):ugdnvz_hi(2),ugdnvz_lo(3):ugdnvz_hi(3),0:ngroups-1)
    double precision rflux1(rfd1_lo(1):rfd1_hi(1),rfd1_lo(2):rfd1_hi(2),rfd1_lo(3):rfd1_hi(3),0:ngroups-1)
    double precision rflux2(rfd2_lo(1):rfd2_hi(1),rfd2_lo(2):rfd2_hi(2),rfd2_lo(3):rfd2_hi(3),0:ngroups-1)
    double precision rflux3(rfd3_lo(1):rfd3_hi(1),rfd3_lo(2):rfd3_hi(2),rfd3_lo(3):rfd3_hi(3),0:ngroups-1)

    integer :: domlo(3), domhi(3)

    ! Local variables

    integer km,kc,kt,k3d,n
    integer i,j, g

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
    double precision, allocatable::  rfx(:,:,:,:), rfy(:,:,:,:),rfz(:,:,:,:)

    double precision, allocatable:: fxy(:,:,:,:), fxz(:,:,:,:)
    double precision, allocatable:: fyx(:,:,:,:), fyz(:,:,:,:)
    double precision, allocatable:: fzx(:,:,:,:), fzy(:,:,:,:)
    double precision, allocatable::rfxy(:,:,:,:),rfxz(:,:,:,:)
    double precision, allocatable::rfyx(:,:,:,:),rfyz(:,:,:,:)
    double precision, allocatable::rfzx(:,:,:,:),rfzy(:,:,:,:)

    double precision, allocatable:: qgdnvx(:,:,:,:), qgdnvxf(:,:,:,:), qgdnvtmpx(:,:,:,:)
    double precision, allocatable:: qgdnvy(:,:,:,:), qgdnvyf(:,:,:,:), qgdnvtmpy(:,:,:,:)
    double precision, allocatable:: qgdnvz(:,:,:,:), qgdnvtmpz1(:,:,:,:), qgdnvtmpz2(:,:,:,:), qgdnvzf(:,:,:,:)

    double precision, allocatable:: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)

    double precision, pointer :: shk(:,:,:)

    double precision :: pggdnvx, pggdnvy, pggdnvz

    integer :: qt_lo(3), qt_hi(3)
    integer :: It_lo(3), It_hi(3)
    integer :: shk_lo(3), shk_hi(3) 
    integer :: fx_lo(3), fx_hi(3)
    integer :: fy_lo(3), fy_hi(3)
    integer :: fz_lo(3), fz_hi(3)

    qt_lo = [lo(1)-1, lo(2)-1, 1]
    qt_hi = [hi(1)+2, hi(2)+2, 2]

    It_lo = [lo(1)-1, lo(2)-1, 1]
    It_hi = [hi(1)+1, hi(2)+1, 2]

    shk_lo(:) = lo(:) - 1
    shk_hi(:) = hi(:) + 1 

    fx_lo = [lo(1)    , lo(2) - 1, 1]
    fx_hi = [hi(1) + 1, hi(2) + 1, 2]

    fy_lo = [lo(1) - 1, lo(2)    , 1]
    fy_hi = [hi(1) + 1, hi(2) + 1, 2]

    fz_lo = [lo(1) - 1, lo(2) - 1, 1]
    fz_hi = [hi(1) + 1, hi(2) + 1, 2]

    allocate (    qgdnvx(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),ngdnv))
    allocate (   qgdnvxf(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),ngdnv))
    allocate ( qgdnvtmpx(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),ngdnv))

    allocate (    qgdnvy(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),ngdnv))
    allocate (   qgdnvyf(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),ngdnv))
    allocate ( qgdnvtmpy(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),ngdnv))

    allocate (     qgdnvz(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),ngdnv))
    allocate ( qgdnvtmpz1(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),ngdnv))
    allocate ( qgdnvtmpz2(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),ngdnv))
    allocate (    qgdnvzf(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),ngdnv))

    allocate ( dqx(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( dqy(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( dqz(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))

    allocate ( qxm(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qxp(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))

    allocate ( qmxy(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qpxy(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))

    allocate ( qmxz(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qpxz(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))

    allocate ( qym(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qyp(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))

    allocate ( qmyx(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qpyx(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))

    allocate ( qmyz(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qpyz(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))

    allocate ( qzm(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qzp(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))

    allocate ( qxl(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qxr(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qyl(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qyr(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qzl(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qzr(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))

    allocate ( qmzx(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qpzx(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))

    allocate ( qmzy(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))
    allocate ( qpzy(qt_lo(1):qt_hi(1),qt_lo(2):qt_hi(2),qt_lo(3):qt_hi(3),QRADVAR))

    allocate ( fx(fx_lo(1):fx_hi(1), fx_lo(2):fx_hi(2), fx_lo(3):fx_hi(3), NVAR))
    allocate (rfx(fx_lo(1):fx_hi(1), fx_lo(2):fx_hi(2), fx_lo(3):fx_hi(3), 0:ngroups-1))

    allocate ( fy(fy_lo(1):fy_hi(1), fy_lo(2):fy_hi(2), fy_lo(3):fy_hi(3), NVAR))
    allocate (rfy(fy_lo(1):fy_hi(1), fy_lo(2):fy_hi(2), fy_lo(3):fy_hi(3), 0:ngroups-1))

    allocate ( fz(fz_lo(1):fz_hi(1), fz_lo(2):fz_hi(2), fz_lo(3):fz_hi(3), NVAR))
    allocate (rfz(fz_lo(1):fz_hi(1), fz_lo(2):fz_hi(2), fz_lo(3):fz_hi(3), 0:ngroups-1))

    allocate ( fxy(fx_lo(1):fx_hi(1), fx_lo(2):fx_hi(2), fx_lo(3):fx_hi(3), NVAR))
    allocate (rfxy(fx_lo(1):fx_hi(1), fx_lo(2):fx_hi(2), fx_lo(3):fx_hi(3), 0:ngroups-1))
    allocate ( fxz(fx_lo(1):fx_hi(1), fx_lo(2):fx_hi(2), fx_lo(3):fx_hi(3), NVAR))
    allocate (rfxz(fx_lo(1):fx_hi(1), fx_lo(2):fx_hi(2), fx_lo(3):fx_hi(3), 0:ngroups-1))

    allocate ( fyx(fy_lo(1):fy_hi(1), fy_lo(2):fy_hi(2), fy_lo(3):fy_hi(3), NVAR))
    allocate (rfyx(fy_lo(1):fy_hi(1), fy_lo(2):fy_hi(2), fy_lo(3):fy_hi(3), 0:ngroups-1))
    allocate ( fyz(fy_lo(1):fy_hi(1), fy_lo(2):fy_hi(2), fy_lo(3):fy_hi(3), NVAR))
    allocate (rfyz(fy_lo(1):fy_hi(1), fy_lo(2):fy_hi(2), fy_lo(3):fy_hi(3), 0:ngroups-1))

    allocate ( fzx(fz_lo(1):fz_hi(1), fz_lo(2):fz_hi(2), fz_lo(3):fz_hi(3), NVAR))
    allocate (rfzx(fz_lo(1):fz_hi(1), fz_lo(2):fz_hi(2), fz_lo(3):fz_hi(3), 0:ngroups-1))
    allocate ( fzy(fz_lo(1):fz_hi(1), fz_lo(2):fz_hi(2), fz_lo(3):fz_hi(3), NVAR))
    allocate (rfzy(fz_lo(1):fz_hi(1), fz_lo(2):fz_hi(2), fz_lo(3):fz_hi(3), 0:ngroups-1))

    ! x-index, y-index, z-index, dim, characteristics, variables
    allocate ( Ip(It_lo(1):It_hi(1), It_lo(2):It_hi(2), It_lo(3):It_hi(3),3,3,QRADVAR))
    allocate ( Im(It_lo(1):It_hi(1), It_lo(2):It_hi(2), It_lo(3):It_hi(3),3,3,QRADVAR))

    ! for the hybrid Riemann solver
    call bl_allocate(shk, shk_lo(1), shk_hi(1), shk_lo(2), shk_hi(2), shk_lo(3), shk_hi(3))

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


    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q, qd_lo, qd_hi, shk, shk_lo, shk_hi, lo, hi, (/dx,dy,dz/))
    else
       shk(:,:,:) = ZERO
    endif


    ! Initialize kc (current k-level) and km (previous k-level)
    kc = 1
    km = 2

    do k3d = lo(3)-1, hi(3)+1

       ! Swap pointers to levels
       kt = km
       km = kc
       kc = kt

       if (ppm_type .le. 0) then
          call bl_error("ppm_type <=0 is not supported in umeth3d_rad")
       else

          do n=1,QRADVAR
             call ppm(q(:,:,:,n),qd_lo,qd_hi, &
                      q(:,:,:,QU:),c,qd_lo,qd_hi, &
                      flatn,qd_lo,qd_hi, &
                      Ip(:,:,:,:,:,n),Im(:,:,:,:,:,n), It_lo, It_hi, &
                      lo(1),lo(2),hi(1),hi(2),(/dx,dy,dz/),dt,k3d,kc)
          end do

          ! Compute U_x and U_y at kc (k3d)
          call tracexy_ppm_rad(lam, lam_lo, lam_hi, &
                               q, c, cg, flatn, qd_lo, qd_hi, &
                               Ip,Im, &
                               qxm, qxp, qym, qyp, qt_lo, qt_hi, &
                               lo(1),lo(2),hi(1),hi(2),dx,dy,dt,kc,k3d)

       end if

       ! Compute \tilde{F}^x at kc (k3d)
       call cmpflx(qxm, qxp, qt_lo, qt_hi, &
                   fx, fx_lo, fx_hi, &
                   qgdnvx, qt_lo, qt_hi, &
                   lam, lam_lo, lam_hi, &
                   rfx, fx_lo, fx_hi, &
                   gamcg, &
                   gamc, csml, c, qd_lo, qd_hi, &
                   shk, shk_lo, shk_hi, &
                   1, lo(1), hi(1)+1, lo(2)-1, hi(2)+1, kc, kc, k3d, domlo, domhi)

       ! Compute \tilde{F}^y at kc (k3d)
       call cmpflx(qym, qyp, qt_lo, qt_hi, &
                   fy, fy_lo, fy_hi, &
                   qgdnvy, qt_lo, qt_hi, &
                   lam, lam_lo, lam_hi, &
                   rfy, fy_lo, fy_hi, &
                   gamcg, &
                   gamc, csml, c, qd_lo, qd_hi, &
                   shk, shk_lo, shk_hi, &
                   2, lo(1)-1, hi(1)+1, lo(2), hi(2)+1, kc, kc, k3d, domlo, domhi)

       ! Compute U'^y_x at kc (k3d)
       call transy1_rad(lam, lam_lo, lam_hi, &
                        qxm, qmxy, qxp, qpxy, qt_lo, qt_hi, &
                        fy, rfy, fy_lo, fy_hi, &
                        qgdnvy, qt_lo, qt_hi, &
                        gamcg, qd_lo, qd_hi, &
                        cdtdy, lo(1)-1, hi(1)+1, lo(2), hi(2), kc, k3d)

       ! Compute U'^x_y at kc (k3d)
       call transx1_rad(lam, lam_lo, lam_hi, &
                        qym, qmyx, qyp, qpyx, qt_lo, qt_hi, &
                        fx, rfx, fx_lo, fx_hi, &
                        qgdnvx, qt_lo, qt_hi, &
                        gamcg, qd_lo, qd_hi, &
                        cdtdx, lo(1), hi(1), lo(2)-1, hi(2)+1, kc, k3d)

       ! Compute F^{x|y} at kc (k3d)
       call cmpflx(qmxy, qpxy, qt_lo, qt_hi, &
                   fxy, fx_lo, fx_hi, &
                   qgdnvtmpx, qt_lo, qt_hi, &
                   lam, lam_lo, lam_hi, &
                   rfxy, fx_lo, fx_hi, &
                   gamcg, &
                   gamc, csml, c, qd_lo, qd_hi, &
                   shk, shk_lo, shk_hi, &
                   1, lo(1), hi(1)+1, lo(2), hi(2), kc, kc, k3d, domlo, domhi)

       ! Compute F^{y|x} at kc (k3d)
       call cmpflx(qmyx, qpyx, qt_lo, qt_hi, &
                   fyx, fy_lo, fy_hi, &
                   qgdnvtmpy, qt_lo, qt_hi, &
                   lam, lam_lo, lam_hi, &
                   rfyx, fy_lo, fy_hi, &
                   gamcg, &
                   gamc, csml, c, qd_lo, qd_hi, &
                   shk, shk_lo, shk_hi, &
                   2, lo(1), hi(1), lo(2), hi(2)+1, kc, kc, k3d, domlo, domhi)

       if (k3d .ge. lo(3)) then

          ! Compute U_z at kc (k3d)
          call tracez_ppm_rad(lam, lam_lo, lam_hi, &
                              q, c, cg, flatn, qd_lo, qd_hi, &
                              Ip,Im, &
                              qzm, qzp, qt_lo, qt_hi, &
                              lo(1), lo(2), hi(1), hi(2), dz, dt, km, kc, k3d)

          ! Compute \tilde{F}^z at kc (k3d)
          call cmpflx(qzm, qzp, qt_lo, qt_hi, &
                      fz, fz_lo, fz_hi, &
                      qgdnvz, qt_lo, qt_hi, &
                      lam, lam_lo, lam_hi, &
                      rfz, fz_lo, fz_hi, &
                      gamcg, &
                      gamc, csml, c, qd_lo, qd_hi, &
                      shk, shk_lo, shk_hi, &
                      3, lo(1)-1, hi(1)+1, lo(2)-1, hi(2)+1, kc, kc, k3d, domlo, domhi)

          ! Compute U'^y_z at kc (k3d)
          call transy2_rad(lam, lam_lo, lam_hi, &
                           qzm, qmzy, qzp, qpzy, qt_lo, qt_hi, &
                           fy, rfy, fy_lo, fy_hi, &
                           qgdnvy, qt_lo, qt_hi, &
                           gamcg, qd_lo, qd_hi, &
                           cdtdy, lo(1)-1, hi(1)+1, lo(2), hi(2), kc, km, k3d)
          
          ! Compute U'^x_z at kc (k3d)
          call transx2_rad(lam, lam_lo, lam_hi, &
                           qzm, qmzx, qzp, qpzx, qt_lo, qt_hi, &
                           fx, rfx, fx_lo, fx_hi, &
                           qgdnvx, qt_lo, qt_hi, &
                           gamcg, qd_lo, qd_hi, &
                           cdtdx, lo(1), hi(1), lo(2)-1, hi(2)+1, kc, km, k3d)

          ! Compute F^{z|x} at kc (k3d)
          call cmpflx(qmzx, qpzx, qt_lo, qt_hi, &
                      fzx, fz_lo, fz_hi, &
                      qgdnvtmpz1, qt_lo, qt_hi, &
                      lam, lam_lo, lam_hi, &
                      rfzx, fz_lo, fz_hi, &
                      gamcg, &
                      gamc, csml, c, qd_lo, qd_hi, &
                      shk, shk_lo, shk_hi, &
                      3, lo(1), hi(1), lo(2)-1, hi(2)+1, kc, kc, k3d, domlo, domhi)

          ! Compute F^{z|y} at kc (k3d)
          call cmpflx(qmzy, qpzy, qt_lo, qt_hi, &
                      fzy, fz_lo, fz_hi, &
                      qgdnvtmpz2, qt_lo, qt_hi, &
                      lam, lam_lo, lam_hi, &
                      rfzy, fz_lo, fz_hi, &
                      gamcg, &
                      gamc, csml, c, qd_lo, qd_hi, &
                      shk, shk_lo, shk_hi, &
                      3, lo(1)-1, hi(1)+1, lo(2), hi(2), kc, kc, k3d, domlo, domhi)
          
          ! Compute U''_z at kc (k3d)
          call transxy_rad(lam, lam_lo, lam_hi, &
                           qzm, qzl, qzp, qzr, qt_lo, qt_hi, &
                           fxy, rfxy, fx_lo, fx_hi, &
                           fyx, rfyx, fy_lo, fy_hi, &
                           qgdnvtmpx, qt_lo, qt_hi, &
                           qgdnvtmpy, qt_lo, qt_hi, &
                           gamcg, qd_lo, qd_hi, &
                           srcQ, src_lo, src_hi, &
                           hdt, hdtdx, hdtdy, lo(1), hi(1), lo(2), hi(2), kc, km, k3d)

          ! Compute F^z at kc (k3d) -- note that flux3 is indexed by k3d, not kc
          call cmpflx(qzl, qzr, qt_lo, qt_hi, &
                      flux3, fd3_lo, fd3_hi, &
                      qgdnvzf, qt_lo, qt_hi, &
                      lam, lam_lo, lam_hi, &
                      rflux3, rfd3_lo, rfd3_hi, &
                      gamcg, &                    
                      gamc, csml, c, qd_lo, qd_hi, &
                      shk, shk_lo, shk_hi, &
                      3,lo(1),hi(1),lo(2),hi(2),kc,k3d,k3d,domlo,domhi)
          
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                ugdnvz_out(i,j,k3d) = qgdnvzf(i,j,kc,GDW)
                ux_zfc    (i,j,k3d) = qgdnvzf(i,j,kc,GDU)
                uy_zfc    (i,j,k3d) = qgdnvzf(i,j,kc,GDU)
             end do
          end do

          do g=0,ngroups-1
             do j=lo(2)-1,hi(2)+1
                do i=lo(1)-1,hi(1)+1
                   ergdz_out(i,j,k3d,g) = qgdnvzf(i,j,kc,GDERADS+g)
                   lmgdz_out(i,j,k3d,g) = qgdnvzf(i,j,kc,GDLAMS+g)
                end do
             end do
          end do

          if (k3d .ge. lo(3)+1 .and. k3d .le. hi(3)+1) then
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   pggdnvz = 0.5d0*(qgdnvzf(i,j,kc,GDPRES) + qgdnvzf(i,j,km,GDPRES)) 
                   pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
                        pggdnvz * (qgdnvzf(i,j,kc,GDW) - qgdnvzf(i,j,km,GDW))/dz
                end do
             end do
          end if

          if (k3d .gt. lo(3)) then

             ! Compute U'^z_x and U'^z_y at km (k3d-1) -- note flux3 has physical index
             call transz_rad(lam, lam_lo, lam_hi, &
                             qxm, qmxz, qxp, qpxz, qym, qmyz, qyp, qpyz, qt_lo, qt_hi, &
                             fz, rfz, fz_lo, fz_hi, &
                             qgdnvz, qt_lo, qt_hi, &
                             gamcg, qd_lo, qd_hi, &
                             cdtdz, lo(1)-1, hi(1)+1, lo(2)-1, hi(2)+1, km, kc, k3d)
             
             ! Compute F^{x|z} at km (k3d-1)
             call cmpflx(qmxz, qpxz, qt_lo, qt_hi, &
                         fxz, fx_lo, fx_hi, &
                         qgdnvx, qt_lo, qt_hi, &
                         lam, lam_lo, lam_hi, &
                         rfxz, fx_lo, fx_hi, &
                         gamcg, &
                         gamc, csml, c, qd_lo, qd_hi, &
                         shk, shk_lo, shk_hi, &
                         1, lo(1), hi(1)+1, lo(2)-1, hi(2)+1, km, km, k3d-1, domlo, domhi)
             
             ! Compute F^{y|z} at km (k3d-1)
             call cmpflx(qmyz, qpyz, qt_lo, qt_hi, &
                         fyz, fy_lo, fy_hi, &
                         qgdnvy, qt_lo, qt_hi, &
                         lam, lam_lo, lam_hi, &
                         rfyz, fy_lo, fy_hi, &
                         gamcg, &
                         gamc, csml, c, qd_lo, qd_hi, &
                         shk, shk_lo, shk_hi, &
                         2, lo(1)-1, hi(1)+1, lo(2), hi(2)+1, km, km, k3d-1, domlo, domhi)

             ! Compute U''_x at km (k3d-1)
             call transyz_rad(lam, lam_lo, lam_hi, &
                              qxm, qxl, qxp, qxr, qt_lo, qt_hi, &
                              fyz, rfyz, fy_lo, fy_hi, &
                              fzy, rfzy, fz_lo, fz_hi, & 
                              qgdnvy, qt_lo, qt_hi, &
                              qgdnvtmpz2, qt_lo, qt_hi, &
                              gamcg, qd_lo, qd_hi, &
                              srcQ, src_lo, src_hi, &
                              hdt, hdtdy, hdtdz, lo(1)-1, hi(1)+1, lo(2), hi(2), km, kc, k3d-1)
             
             ! Compute U''_y at km (k3d-1)
             call transxz_rad(lam, lam_lo, lam_hi, &
                              qym, qyl, qyp, qyr, qt_lo, qt_hi, &
                              fxz, rfxz, fx_lo, fx_hi, &
                              fzx, rfzx, fz_lo, fz_hi, &
                              qgdnvx, qt_lo, qt_hi, &
                              qgdnvtmpz1, qt_lo, qt_hi, &
                              gamcg, qd_lo, qd_hi, &
                              srcQ, src_lo, src_hi, &
                              hdt, hdtdx, hdtdz, lo(1), hi(1), lo(2)-1, hi(2)+1, km, kc, k3d-1)
             
             ! Compute F^x at km (k3d-1)
             call cmpflx(qxl, qxr, qt_lo, qt_hi, &
                         flux1, fd1_lo, fd1_hi, &
                         qgdnvxf, qt_lo, qt_hi, &
                         lam, lam_lo, lam_hi, &
                         rflux1, rfd1_lo, rfd1_hi, &
                         gamcg, &                           
                         gamc, csml, c, qd_lo, qd_hi, &
                         shk, shk_lo, shk_hi, &
                         1,lo(1),hi(1)+1,lo(2),hi(2),km,k3d-1,k3d-1,domlo,domhi)

             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+2
                   ugdnvx_out(i,j,k3d-1) = qgdnvxf(i,j,km,GDU)
                   uy_xfc    (i,j,k3d-1) = qgdnvxf(i,j,km,GDV)
                   uz_xfc    (i,j,k3d-1) = qgdnvxf(i,j,km,GDW)
                end do
             end do

             do g=0,ngroups-1
                do j = lo(2)-1, hi(2)+1
                   do i = lo(1)-1, hi(1)+2
                      ergdx_out(i,j,k3d-1,g) = qgdnvxf(i,j,km,GDERADS+g)
                      lmgdx_out(i,j,k3d-1,g) = qgdnvxf(i,j,km,GDLAMS+g)
                   end do
                end do
             end do

             ! Compute F^y at km (k3d-1)
             call cmpflx(qyl, qyr, qt_lo, qt_hi, &
                         flux2,  fd2_lo, fd2_hi, &
                         qgdnvyf, qt_lo, qt_hi, &
                         lam, lam_lo, lam_hi, &
                         rflux2, rfd2_lo, rfd2_hi, &
                         gamcg, &
                         gamc, csml, c, qd_lo, qd_hi, &
                         shk, shk_lo, shk_hi, &
                         2,lo(1),hi(1),lo(2),hi(2)+1,km,k3d-1,k3d-1,domlo,domhi)

             do j = lo(2)-1, hi(2)+2
                do i = lo(1)-1, hi(1)+1
                   ugdnvy_out(i,j,k3d-1) = qgdnvyf(i,j,km,GDV)
                   ux_yfc    (i,j,k3d-1) = qgdnvyf(i,j,km,GDU)
                   uz_yfc    (i,j,k3d-1) = qgdnvyf(i,j,km,GDW)
                end do
             end do

             do g=0,ngroups-1
                do j = lo(2)-1, hi(2)+2
                   do i = lo(1)-1, hi(1)+1
                      ergdy_out(i,j,k3d-1,g) = qgdnvyf(i,j,km,GDERADS+g)
                      lmgdy_out(i,j,k3d-1,g) = qgdnvyf(i,j,km,GDLAMS+g)
                   end do
                end do
             end do

             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   pggdnvx = 0.5d0*( qgdnvxf(i+1,j,km,GDPRES) + qgdnvxf(i,j,km,GDPRES)) 
                   pggdnvy = 0.5d0*( qgdnvyf(i,j+1,km,GDPRES) + qgdnvyf(i,j,km,GDPRES)) 
                   pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
                        pggdnvx * (qgdnvxf(i+1,j,km,GDU) - qgdnvxf(i,j,km,GDU))/dx + &
                        pggdnvy * (qgdnvyf(i,j+1,km,GDV) - qgdnvyf(i,j,km,GDV))/dy
                end do
             end do

          end if
       end if
    enddo

    ! Deallocate arrays
    deallocate(qgdnvx,qgdnvxf,qgdnvtmpx)
    deallocate(qgdnvy,qgdnvyf,qgdnvtmpy)
    deallocate(qgdnvz,qgdnvtmpz1,qgdnvtmpz2,qgdnvzf)

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

    deallocate( fx, fy, fz)
    deallocate(rfx,rfy,rfz)
    deallocate( fxy, fxz)
    deallocate(rfxy,rfxz)
    deallocate( fyx, fyz)
    deallocate(rfyx,rfyz)
    deallocate( fzx, fzy)
    deallocate(rfzx,rfzy)

    deallocate(Ip,Im)

    call bl_deallocate(shk)

  end subroutine umeth3d_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine consup_rad(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
       Erin,Erin_l1,Erin_l2,Erin_l3,Erin_h1,Erin_h2,Erin_h3, &
       Erout,Erout_l1,Erout_l2,Erout_l3,Erout_h1,Erout_h2,Erout_h3, &
       ugdx, ergdx, lmgdx, &
       ugdx_l1,ugdx_l2,ugdx_l3,ugdx_h1,ugdx_h2,ugdx_h3, &
       ugdy, ergdy, lmgdy, & 
       ugdy_l1,ugdy_l2,ugdy_l3,ugdy_h1,ugdy_h2,ugdy_h3, &
       ugdz, ergdz, lmgdz, &
       ugdz_l1,ugdz_l2,ugdz_l3,ugdz_h1,ugdz_h2,ugdz_h3, &
       src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
       flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
       flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
       flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
       radflux1,radflux1_l1,radflux1_l2,radflux1_l3,radflux1_h1,radflux1_h2,radflux1_h3, &
       radflux2,radflux2_l1,radflux2_l2,radflux2_l3,radflux2_h1,radflux2_h2,radflux2_h3, &
       radflux3,radflux3_l1,radflux3_l2,radflux3_l3,radflux3_h1,radflux3_h2,radflux3_h3, &
       area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
       area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
       area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
       vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
       div,pdivu, uy_xfc, uz_xfc, ux_yfc, uz_yfc, ux_zfc, uy_zfc, &
       lo,hi,dx,dy,dz,dt, nstep_fsp)

    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
         UEDEN, UEINT, UTEMP, normalize_species
    use rad_params_module, only : ngroups, nugroup, dlognu
    use radhydro_params_module, only : fspace_type, comoving
    use radhydro_nd_module, only : advect_in_fspace
    use fluxlimiter_module, only : Edd_factor
    use advection_util_module, only : normalize_species_fluxes

    implicit none

    integer nstep_fsp
    integer lo(3), hi(3)
    integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
    integer  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
    integer   src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3 
    integer Erout_l1,Erout_l2,Erout_l3,Erout_h1,Erout_h2,Erout_h3
    integer Erin_l1,Erin_l2,Erin_l3,Erin_h1,Erin_h2,Erin_h3
    integer ugdx_l1,ugdx_l2,ugdx_l3,ugdx_h1,ugdx_h2,ugdx_h3
    integer ugdy_l1,ugdy_l2,ugdy_l3,ugdy_h1,ugdy_h2,ugdy_h3
    integer ugdz_l1,ugdz_l2,ugdz_l3,ugdz_h1,ugdz_h2,ugdz_h3
    integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
    integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
    integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
    integer radflux1_l1,radflux1_l2,radflux1_l3,radflux1_h1,radflux1_h2,radflux1_h3
    integer radflux2_l1,radflux2_l2,radflux2_l3,radflux2_h1,radflux2_h2,radflux2_h3
    integer radflux3_l1,radflux3_l2,radflux3_l3,radflux3_h1,radflux3_h2,radflux3_h3
    integer area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3
    integer area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3
    integer area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3
    integer vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3

    double precision  Erin( Erin_l1: Erin_h1,  Erin_l2: Erin_h2,  Erin_l3: Erin_h3,0:ngroups-1)
    double precision Erout(Erout_l1:Erout_h1, Erout_l2:Erout_h2, Erout_l3:Erout_h3,0:ngroups-1)
    double precision  ugdx(ugdx_l1:ugdx_h1, ugdx_l2:ugdx_h2, ugdx_l3:ugdx_h3)
    double precision ergdx(ugdx_l1:ugdx_h1, ugdx_l2:ugdx_h2, ugdx_l3:ugdx_h3,0:ngroups-1)
    double precision lmgdx(ugdx_l1:ugdx_h1, ugdx_l2:ugdx_h2, ugdx_l3:ugdx_h3,0:ngroups-1)
    double precision  ugdy(ugdy_l1:ugdy_h1, ugdy_l2:ugdy_h2, ugdy_l3:ugdy_h3)
    double precision ergdy(ugdy_l1:ugdy_h1, ugdy_l2:ugdy_h2, ugdy_l3:ugdy_h3,0:ngroups-1)
    double precision lmgdy(ugdy_l1:ugdy_h1, ugdy_l2:ugdy_h2, ugdy_l3:ugdy_h3,0:ngroups-1)
    double precision  ugdz(ugdz_l1:ugdz_h1, ugdz_l2:ugdz_h2, ugdz_l3:ugdz_h3)
    double precision ergdz(ugdz_l1:ugdz_h1, ugdz_l2:ugdz_h2, ugdz_l3:ugdz_h3,0:ngroups-1)
    double precision lmgdz(ugdz_l1:ugdz_h1, ugdz_l2:ugdz_h2, ugdz_l3:ugdz_h3,0:ngroups-1)
    double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
    double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
    double precision   src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
    double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
    double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
    double precision flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
    double precision radflux1(radflux1_l1:radflux1_h1,radflux1_l2:radflux1_h2,radflux1_l3:radflux1_h3,0:ngroups-1)
    double precision radflux2(radflux2_l1:radflux2_h1,radflux2_l2:radflux2_h2,radflux2_l3:radflux2_h3,0:ngroups-1)
    double precision radflux3(radflux3_l1:radflux3_h1,radflux3_l2:radflux3_h2,radflux3_l3:radflux3_h3,0:ngroups-1)
    double precision area1(area1_l1:area1_h1,area1_l2:area1_h2,area1_l3:area1_h3)
    double precision area2(area2_l1:area2_h1,area2_l2:area2_h2,area2_l3:area2_h3)
    double precision area3(area3_l1:area3_h1,area3_l2:area3_h2,area3_l3:area3_h3)
    double precision vol(vol_l1:vol_h1,vol_l2:vol_h2,vol_l3:vol_h3)
    double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    double precision pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision uy_xfc(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2,ugdx_l3:ugdx_h3)
    double precision uz_xfc(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2,ugdx_l3:ugdx_h3)
    double precision ux_yfc(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2,ugdy_l3:ugdy_h3)
    double precision uz_yfc(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2,ugdy_l3:ugdy_h3)
    double precision ux_zfc(ugdz_l1:ugdz_h1,ugdz_l2:ugdz_h2,ugdz_l3:ugdz_h3)
    double precision uy_zfc(ugdz_l1:ugdz_h1,ugdz_l2:ugdz_h2,ugdz_l3:ugdz_h3)
    double precision dx, dy, dz, dt

    ! Local variables

    integer :: flux1_lo(3), flux1_hi(3)
    integer :: flux2_lo(3), flux2_hi(3)
    integer :: flux3_lo(3), flux3_hi(3)

    double precision :: div1
    double precision :: rho, Up, Vp, Wp
    double precision :: SrU, SrV, SrW, SrE
    integer          :: i, j, k, n, g

    double precision, dimension(0:ngroups-1) :: Erscale
    double precision, dimension(0:ngroups-1) :: ustar, af
    double precision :: Eddf, Eddfxm, Eddfxp, Eddfym, Eddfyp, Eddfzm, Eddfzp 
    double precision :: f1, f2, f1xm, f1xp, f1ym, f1yp, f1zm, f1zp
    double precision :: Gf1E(3)
    double precision :: ux, uy, uz, divu, lamc, Egdc
    double precision :: dudx(3), dudy(3), dudz(3), nhat(3), GnDotu(3), nnColonDotGu
    double precision :: dprdx, dprdy, dprdz, ek1, ek2, dek  

    if (ngroups .gt. 1) then
       if (fspace_type .eq. 1) then
          Erscale = dlognu
       else
          Erscale = nugroup*dlognu
       end if
    end if

    flux1_lo = (/ flux1_l1, flux1_l2, flux1_l3 /)
    flux1_hi = (/ flux1_h1, flux1_h2, flux1_h3 /)

    flux2_lo = (/ flux2_l1, flux2_l2, flux2_l3 /)
    flux2_hi = (/ flux2_h1, flux2_h2, flux2_h3 /)

    flux3_lo = (/ flux3_l1, flux3_l2, flux3_l3 /)
    flux3_hi = (/ flux3_h1, flux3_h2, flux3_h3 /)


    do n = 1, NVAR

       if ( n.eq.UTEMP ) then

          flux1(:,:,:,n) = 0.d0
          flux2(:,:,:,n) = 0.d0
          flux3(:,:,:,n) = 0.d0

       else

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

       endif

    enddo

    do g=0,ngroups-1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                div1 = .25d0*(div(i,j,k) + div(i,j+1,k) + div(i,j,k+1) + div(i,j+1,k+1))
                div1 = difmag*min(0.d0,div1)
                radflux1(i,j,k,g) = radflux1(i,j,k,g) + dx*div1*(Erin(i,j,k,g)-Erin(i-1,j,k,g))
                radflux1(i,j,k,g) = radflux1(i,j,k,g) * area1(i,j,k) * dt
             enddo
          enddo
       enddo
    enddo

    do g=0,ngroups-1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                div1 = .25d0*(div(i,j,k) + div(i+1,j,k) + div(i,j,k+1) + div(i+1,j,k+1))
                div1 = difmag*min(0.d0,div1)
                radflux2(i,j,k,g) = radflux2(i,j,k,g) + dy*div1*(Erin(i,j,k,g)-Erin(i,j-1,k,g))
                radflux2(i,j,k,g) = radflux2(i,j,k,g) * area2(i,j,k) * dt
             enddo
          enddo
       enddo
    enddo

    do g=0,ngroups-1
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                div1 = .25d0*(div(i,j,k) + div(i+1,j,k) + div(i,j+1,k) + div(i+1,j+1,k))
                div1 = difmag*min(0.d0,div1)
                radflux3(i,j,k,g) = radflux3(i,j,k,g) + dz*div1*(Erin(i,j,k,g)-Erin(i,j,k-1,g))
                radflux3(i,j,k,g) = radflux3(i,j,k,g) * area3(i,j,k) * dt
             enddo
          enddo
       enddo
    enddo

    if (normalize_species .eq. 1) &
         call normalize_species_fluxes(flux1,flux1_lo,flux1_hi, &
         flux2,flux2_lo,flux2_hi, &
         flux3,flux3_lo,flux3_hi, &
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
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   uout(i,j,k,n) = uin(i,j,k,n) &
                        + ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                        +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                        +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) / vol(i,j,k)
                   !
                   ! Add the source term to (rho e)
                   !
                   if (n .eq. UEINT) then
                      uout(i,j,k,n) = uout(i,j,k,n) - dt * pdivu(i,j,k)
                   endif
                enddo
             enddo
          enddo
       endif

    enddo

    ! update everything else with fluxes
    do g=0,ngroups-1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                Erout(i,j,k,g) = Erin(i,j,k,g) &
                     + ( radflux1(i,j,k,g) - radflux1(i+1,j,k,g) &
                     +   radflux2(i,j,k,g) - radflux2(i,j+1,k,g) &
                     +   radflux3(i,j,k,g) - radflux3(i,j,k+1,g)) / vol(i,j,k)
             enddo
          enddo
       enddo
    enddo

    ! add radiation force terms
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             dprdx = 0.d0
             dprdy = 0.d0
             dprdz = 0.d0
             do g=0,ngroups-1
                lamc = (lmgdx(i,j,k,g)+lmgdx(i+1,j,k,g) &
                     +  lmgdy(i,j,k,g)+lmgdy(i,j+1,k,g) &
                     +  lmgdz(i,j,k,g)+lmgdz(i,j,k+1,g) ) / 6.d0
                dprdx = dprdx + lamc*(ergdx(i+1,j,k,g)-ergdx(i,j,k,g))/dx
                dprdy = dprdy + lamc*(ergdy(i,j+1,k,g)-ergdy(i,j,k,g))/dy
                dprdz = dprdz + lamc*(ergdz(i,j,k+1,g)-ergdz(i,j,k,g))/dz
             end do

             ek1 = (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) &
                  / (2.d0*uout(i,j,k,URHO))

             uout(i,j,k,UMX) = uout(i,j,k,UMX) - dt * dprdx
             uout(i,j,k,UMY) = uout(i,j,k,UMY) - dt * dprdy
             uout(i,j,k,UMZ) = uout(i,j,k,UMZ) - dt * dprdz
             ek2 = (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) &
                  / (2.d0*uout(i,j,k,URHO))
             dek = ek2 - ek1

             uout(i,j,k,UEDEN) = uout(i,j,k,UEDEN) + dek
             if (.not. comoving) then ! mixed-frame (single group only)
                Erout(i,j,k,0) = Erout(i,j,k,0) - dek
             end if

          end do
       end do
    end do

    ! Add radiation source terms
    if (comoving) then
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                ux = 0.5d0*(ugdx(i,j,k) + ugdx(i+1,j,k))
                uy = 0.5d0*(ugdy(i,j,k) + ugdy(i,j+1,k))
                uz = 0.5d0*(ugdz(i,j,k) + ugdz(i,j,k+1))

                dudx(1) = (ugdx(i+1,j,k)-ugdx(i,j,k))/dx 
                dudx(2) = (uy_xfc(i+1,j,k)-uy_xfc(i,j,k))/dx 
                dudx(3) = (uz_xfc(i+1,j,k)-uz_xfc(i,j,k))/dx 

                dudy(1) = (ux_yfc(i,j+1,k)-ux_yfc(i,j,k))/dy 
                dudy(2) = (ugdy(i,j+1,k)-ugdy(i,j,k))/dy 
                dudy(3) = (uz_yfc(i,j+1,k)-uz_yfc(i,j,k))/dy 

                dudz(1) = (ux_zfc(i,j,k+1)-ux_zfc(i,j,k))/dz 
                dudz(2) = (uy_zfc(i,j,k+1)-uy_zfc(i,j,k))/dz 
                dudz(3) = (ugdz(i,j,k+1)-ugdz(i,j,k))/dz 

                divu = dudx(1) + dudy(2) + dudz(3)

                ! Note that for single group, fspace_type is always 1
                do g=0, ngroups-1

                   nhat(1) = (ergdx(i+1,j,k,g)-ergdx(i,j,k,g))/dx
                   nhat(2) = (ergdy(i,j+1,k,g)-ergdy(i,j,k,g))/dy
                   nhat(3) = (ergdz(i,j,k+1,g)-ergdz(i,j,k,g))/dz

                   GnDotu(1) = dot_product(nhat, dudx)
                   GnDotu(2) = dot_product(nhat, dudy)
                   GnDotu(3) = dot_product(nhat, dudz)

                   nnColonDotGu = dot_product(nhat, GnDotu) / (dot_product(nhat,nhat)+1.d-50)

                   lamc = (lmgdx(i,j,k,g)+lmgdx(i+1,j,k,g) &
                        +  lmgdy(i,j,k,g)+lmgdy(i,j+1,k,g) &
                        +  lmgdz(i,j,k,g)+lmgdz(i,j,k+1,g) ) / 6.d0
                   Eddf = Edd_factor(lamc)
                   f1 = (1.d0-Eddf)*0.5d0
                   f2 = (3.d0*Eddf-1.d0)*0.5d0
                   af(g) = -(f1*divu + f2*nnColonDotGu)

                   if (fspace_type .eq. 1) then
                      Eddfxp = Edd_factor(lmgdx(i+1,j  ,k  ,g))
                      Eddfxm = Edd_factor(lmgdx(i  ,j  ,k  ,g))
                      Eddfyp = Edd_factor(lmgdy(i  ,j+1,k  ,g))
                      Eddfym = Edd_factor(lmgdy(i  ,j  ,k  ,g)) 
                      Eddfzp = Edd_factor(lmgdz(i  ,j  ,k+1,g))
                      Eddfzm = Edd_factor(lmgdz(i  ,j  ,k  ,g)) 

                      f1xp = 0.5d0*(1.d0-Eddfxp)
                      f1xm = 0.5d0*(1.d0-Eddfxm)
                      f1yp = 0.5d0*(1.d0-Eddfyp)
                      f1ym = 0.5d0*(1.d0-Eddfym)
                      f1zp = 0.5d0*(1.d0-Eddfzp)
                      f1zm = 0.5d0*(1.d0-Eddfzm)

                      Gf1E(1) = (f1xp*ergdx(i+1,j,k,g) - f1xm*ergdx(i,j,k,g)) / dx
                      Gf1E(2) = (f1yp*ergdy(i,j+1,k,g) - f1ym*ergdy(i,j,k,g)) / dy
                      Gf1E(3) = (f1zp*ergdz(i,j,k+1,g) - f1zm*ergdz(i,j,k,g)) / dz

                      Egdc = (ergdx(i,j,k,g)+ergdx(i+1,j,k,g) &
                           +  ergdy(i,j,k,g)+ergdy(i,j+1,k,g) &
                           +  ergdz(i,j,k,g)+ergdz(i,j,k+1,g) ) / 6.d0

                      Erout(i,j,k,g) = Erout(i,j,k,g) + dt*(ux*Gf1E(1)+uy*Gf1E(2)+uz*Gf1E(3)) &
                           - dt*f2*Egdc*nnColonDotGu
                   end if

                end do

                if (ngroups.gt.1) then
                   ustar = Erout(i,j,k,:) / Erscale
                   call advect_in_fspace(ustar, af, dt, nstep_fsp)
                   Erout(i,j,k,:) = ustar * Erscale           
                end if
             enddo
          enddo
       enddo
    endif

  end subroutine consup_rad


  subroutine ppflaten(lof, hif, &
       flatn, q, q_l1,q_l2,q_l3, q_h1,q_h2,q_h3)
    use meth_params_module, only : QPRES, QU, QV, QW
    use radhydro_params_module, only : flatten_pp_threshold, QRADVAR, qptot
    implicit none
    integer, intent(in) :: lof(3), hif(3), q_l1, q_h1, q_l2, q_h2, q_l3, q_h3
    double precision, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QRADVAR)
    double precision, intent(inout) :: flatn(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)

    integer :: i,j,k

    do k=lof(3),hif(3)
       do j=lof(2),hif(2)
          do i=lof(1),hif(1)
             if ( q(i-1,j,k,QU)+q(i,j-1,k,QV)+q(i,j,k-1,QW) > &
                  q(i+1,j,k,QU)+q(i,j+1,k,QV)+q(i,j,k+1,QW) ) then
                if (q(i,j,k,QPRES) < flatten_pp_threshold* q(i,j,k,qptot)) then
                   flatn(i,j,k) = 0.d0
                end if
             end if
          end do
       end do
    end do

  end subroutine ppflaten

end module rad_advection_module
