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
             !                                    sum(dpdX_er(i,j,k,:)*(src(i,j,k,UFS:UFS+nspec-1) - &
             !                                                          q(i,j,k,QFS:QFS+nspec-1)*srcQ(i,j,k,QRHO))) &
             !                                    /q(i,j,k,QRHO)

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

  subroutine umeth3d_rad(q, c,cg, gamc,gamcg, csml, flatn, q_lo, q_hi, &
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
    
    use meth_params_module, only : QVAR, NVAR, QU, ppm_type, hybrid_riemann
    use ppm_module
    use radhydro_params_module, only : QRADVAR
    use rad_params_module, only : ngroups
    use riemann_module, only : cmpflx, shock
    use trace_ppm_rad_module, only : tracexy_ppm_rad, tracez_ppm_rad
    use mempool_module, only : bl_allocate, bl_deallocate

    implicit none

    integer :: lam_lo(3), lam_hi(3)
    integer :: src_lo(3), src_hi(3)
    integer :: rfd1_lo(3), rfd1_hi(3)
    integer :: rfd2_lo(3), rfd2_hi(3)
    integer :: rfd3_lo(3), rfd3_hi(3)
    integer :: q_lo(3), q_hi(3)
    integer :: lo(3), hi(3)

    integer :: fd1_lo(3), fd1_hi(3)
    integer :: fd2_lo(3), fd2_hi(3)
    integer :: fd3_lo(3), fd3_hi(3)
    integer :: ugdnvx_lo(3), ugdnvx_hi(3)
    integer :: ugdnvy_lo(3), ugdnvy_hi(3)
    integer :: ugdnvz_lo(3), ugdnvz_hi(3)

    double precision     q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QRADVAR)
    double precision     c(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision    cg(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision  gamc(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision gamcg(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision  csml(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision flatn(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
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

    double precision, allocatable:: qgdnvx(:,:,:), qgdnvxf(:,:,:), qgdnvtmpx(:,:,:)
    double precision, allocatable:: qgdnvy(:,:,:), qgdnvyf(:,:,:), qgdnvtmpy(:,:,:)
    double precision, allocatable:: qgdnvz(:,:,:), qgdnvtmpz1(:,:,:), qgdnvtmpz2(:,:,:), qgdnvzf(:,:,:)

    double precision, allocatable:: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)

    double precision, pointer :: shk(:,:,:)

    double precision :: pggdnvx, pggdnvy, pggdnvz

    integer :: qd_lo(3), qd_hi(3)
    integer :: It_lo(3), It_hi(3)
    integer :: shk_lo(3), shk_hi(3) 

    qd_lo = [lo(1)-1, lo(2)-1, 1]
    qd_hi = [hi(1)+2, hi(2)+2, 2]

    It_lo = [lo(1)-1, lo(2)-1, 1]
    It_hi = [hi(1)+1, hi(2)+1, 2]

    shk_lo(:) = lo(:) - 1
    shk_hi(:) = hi(:) + 1 

    allocate ( qgdnvx(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3)))
    allocate ( qgdnvxf(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3)))
    allocate ( qgdnvtmpx(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3)))

    allocate ( qgdnvy(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3)))
    allocate ( qgdnvyf(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3)))
    allocate ( qgdnvtmpy(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3)))

    allocate ( qgdnvz(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3)))
    allocate ( qgdnvtmpz1(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3)))
    allocate ( qgdnvtmpz2(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3)))
    allocate ( qgdnvzf(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3)))

    allocate ( dqx(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( dqy(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( dqz(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))

    allocate ( qxm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qxp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))

    allocate ( qmxy(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qpxy(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))

    allocate ( qmxz(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qpxz(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))

    allocate ( qym(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qyp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))

    allocate ( qmyx(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qpyx(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))

    allocate ( qmyz(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qpyz(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))

    allocate ( qzm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qzp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))

    allocate ( qxl(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qxr(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qyl(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qyr(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qzl(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qzr(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))

    allocate ( qmzx(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qpzx(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))

    allocate ( qmzy(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))
    allocate ( qpzy(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QRADVAR))

    allocate ( fx(ilo1:ihi1+1,ilo2-1:ihi2+1,2,NVAR))
    allocate (rfx(ilo1:ihi1+1,ilo2-1:ihi2+1,2,0:ngroups-1))
    allocate ( fy(ilo1-1:ihi1+1,ilo2:ihi2+1,2,NVAR))
    allocate (rfy(ilo1-1:ihi1+1,ilo2:ihi2+1,2,0:ngroups-1))
    allocate ( fz(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,NVAR))
    allocate (rfz(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,0:ngroups-1))

    allocate ( fxy(ilo1:ihi1+1,ilo2-1:ihi2+1,2,NVAR))
    allocate (rfxy(ilo1:ihi1+1,ilo2-1:ihi2+1,2,0:ngroups-1))
    allocate ( fxz(ilo1:ihi1+1,ilo2-1:ihi2+1,2,NVAR))
    allocate (rfxz(ilo1:ihi1+1,ilo2-1:ihi2+1,2,0:ngroups-1))

    allocate ( fyx(ilo1-1:ihi1+1,ilo2:ihi2+1,2,NVAR))
    allocate (rfyx(ilo1-1:ihi1+1,ilo2:ihi2+1,2,0:ngroups-1))
    allocate ( fyz(ilo1-1:ihi1+1,ilo2:ihi2+1,2,NVAR))
    allocate (rfyz(ilo1-1:ihi1+1,ilo2:ihi2+1,2,0:ngroups-1))

    allocate ( fzx(ilo1:ihi1,ilo2-1:ihi2+1,2,NVAR))
    allocate (rfzx(ilo1:ihi1,ilo2-1:ihi2+1,2,0:ngroups-1))
    allocate ( fzy(ilo1-1:ihi1+1,ilo2:ihi2,2,NVAR))
    allocate (rfzy(ilo1-1:ihi1+1,ilo2:ihi2,2,0:ngroups-1))

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
       call shock(q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
            shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
            ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dy,dz)
    else
       shk(:,:,:) = ZERO
    endif


    ! Initialize kc (current k-level) and km (previous k-level)
    kc = 1
    km = 2

    do k3d = ilo3-1, ihi3+1

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
                  ilo1,ilo2,ihi1,ihi2,(/dx,dy,dz/),dt,k3d,kc)
          end do

          ! Compute U_x and U_y at kc (k3d)
          call tracexy_ppm_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
               q,c,cg,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
               Ip,Im, &
               qxm,qxp,qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d)

       end if

       ! Compute \tilde{F}^x at kc (k3d)
       call cmpflx(qxm,qxp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
            fx,  ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
            ugdnvx, pgdnvx, gegdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
            lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
            rfx, ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
            v1gdnvtmp, v2gdnvtmp, ergdnvx, lmgdtmp, &
            gamcg, &
            gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
            shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &                 
            1,ilo1,ihi1+1,ilo2-1,ihi2+1,kc,kc,k3d,domlo,domhi)

       ! Compute \tilde{F}^y at kc (k3d)
       call cmpflx(qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
            fy,  ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
            ugdnvy, pgdnvy, gegdnvy, ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
            lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
            rfy, ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
            v1gdnvtmp, v2gdnvtmp, ergdnvy, lmgdtmp, &
            gamcg, &
            gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
            shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
            2,ilo1-1,ihi1+1,ilo2,ihi2+1,kc,kc,k3d,domlo,domhi)

       ! Compute U'^y_x at kc (k3d)
       call transy1_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
            qxm,qmxy,qxp,qpxy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
            fy, rfy, &
            ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
            ugdnvy, pgdnvy, ergdnvy, & 
            ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
            gamcg,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
            cdtdy,ilo1-1,ihi1+1,ilo2,ihi2,kc,k3d)

       ! Compute U'^x_y at kc (k3d)
       call transx1_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
            qym,qmyx,qyp,qpyx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
            fx, rfx, &
            ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
            ugdnvx, pgdnvx, ergdnvx, & 
            ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
            gamcg,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
            cdtdx,ilo1,ihi1,ilo2-1,ihi2+1,kc,k3d)

       ! Compute F^{x|y} at kc (k3d)
       call cmpflx(qmxy,qpxy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
            fxy,  ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
            ugdnvtmpx, pgdnvtmpx, gegdnvtmpx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
            lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
            rfxy, ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
            v1gdnvtmp, v2gdnvtmp, ergdnvtmpx, lmgdtmp, &
            gamcg, &
            gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
            shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
            1,ilo1,ihi1+1,ilo2,ihi2,kc,kc,k3d,domlo,domhi)

       ! Compute F^{y|x} at kc (k3d)
       call cmpflx(qmyx,qpyx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
            fyx,  ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
            ugdnvtmpy, pgdnvtmpy, gegdnvtmpy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
            lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
            rfyx, ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
            v1gdnvtmp, v2gdnvtmp, ergdnvtmpy, lmgdtmp, &
            gamcg, &
            gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
            shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
            2,ilo1,ihi1,ilo2,ihi2+1,kc,kc,k3d,domlo,domhi)

       if (k3d.ge.ilo3) then

          ! Compute U_z at kc (k3d)
          call tracez_ppm_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
               q,c,cg,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
               Ip,Im, &
               qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d)

          ! Compute \tilde{F}^z at kc (k3d)
          call cmpflx(qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               fz,  ilo1-1,ilo2-1,1,ihi1+1,ihi2+1,2, &
               ugdnvz, pgdnvz, gegdnvz, ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
               rfz, ilo1-1,ilo2-1,1,ihi1+1,ihi2+1,2, &
               v1gdnvtmp, v2gdnvtmp, ergdnvz, lmgdtmp, & 
               gamcg, &                    
               gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
               shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
               3,ilo1-1,ihi1+1,ilo2-1,ihi2+1,kc,kc,k3d,domlo,domhi)

          ! Compute U'^y_z at kc (k3d)
          call transy2_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
               qzm,qmzy,qzp,qpzy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               fy, rfy, &
               ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
               ugdnvy, pgdnvy, ergdnvy, & 
               ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               gamcg,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
               cdtdy,ilo1-1,ihi1+1,ilo2,ihi2,kc,km,k3d)

          ! Compute U'^x_z at kc (k3d)
          call transx2_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, & 
               qzm,qmzx,qzp,qpzx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               fx, rfx, &
               ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
               ugdnvx, pgdnvx, ergdnvx, & 
               ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               gamcg,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
               cdtdx,ilo1,ihi1,ilo2-1,ihi2+1,kc,km,k3d)

          ! Compute F^{z|x} at kc (k3d)
          call cmpflx(qmzx,qpzx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               fzx,  ilo1,ilo2-1,1,ihi1,ihi2+1,2, &
               ugdnvtmpz1, pgdnvtmpz1, gegdnvtmpz1, ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
               rfzx, ilo1,ilo2-1,1,ihi1,ihi2+1,2, &
               v1gdnvtmp, v2gdnvtmp, ergdnvtmpz1, lmgdtmp, &
               gamcg, &
               gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
               shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
               3,ilo1,ihi1,ilo2-1,ihi2+1,kc,kc,k3d,domlo,domhi)

          ! Compute F^{z|y} at kc (k3d)
          call cmpflx(qmzy,qpzy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               fzy,  ilo1-1,ilo2,1,ihi1+1,ihi2,2, &
               ugdnvtmpz2, pgdnvtmpz2, gegdnvtmpz2, ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
               rfzy, ilo1-1,ilo2,1,ihi1+1,ihi2,2, &
               v1gdnvtmp, v2gdnvtmp, ergdnvtmpz2, lmgdtmp, &                    
               gamcg, &
               gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
               shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
               3,ilo1-1,ihi1+1,ilo2,ihi2,kc,kc,k3d,domlo,domhi)

          ! Compute U''_z at kc (k3d)
          call transxy_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, & 
               qzm,qzl,qzp,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               fxy, rfxy, &
               ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
               fyx, rfyx, &
               ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
               ugdnvtmpx, pgdnvtmpx, ergdnvtmpx, &
               ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               ugdnvtmpy, pgdnvtmpy, ergdnvtmpy, &
               ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               gamcg,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
               srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
               hdt,hdtdx,hdtdy,ilo1,ihi1,ilo2,ihi2,kc,km,k3d)

          ! Compute F^z at kc (k3d) -- note that flux3 is indexed by k3d, not kc
          call cmpflx(qzl,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               flux3,   fd3_l1, fd3_l2, fd3_l3, fd3_h1, fd3_h2, fd3_h3, &
               ugdnvzf, pgdnvzf, gegdnvzf, ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
               lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
               rflux3, rfd3_l1,rfd3_l2,rfd3_l3,rfd3_h1,rfd3_h2,rfd3_h3, &
               v1gdnvtmp, v2gdnvtmp, ergdnvzf, lmgdtmp, &              
               gamcg, &                    
               gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
               shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
               3,ilo1,ihi1,ilo2,ihi2,kc,k3d,k3d,domlo,domhi)

          do j=ilo2-1,ihi2+1
             do i=ilo1-1,ihi1+1
                ugdnvz_out(i,j,k3d) =   ugdnvzf(i,j,kc)
                ux_zfc    (i,j,k3d) = v1gdnvtmp(i,j,kc)
                uy_zfc    (i,j,k3d) = v2gdnvtmp(i,j,kc)
             end do
          end do

          do g=0,ngroups-1
             do j=ilo2-1,ihi2+1
                do i=ilo1-1,ihi1+1
                   ergdz_out(i,j,k3d,g) = ergdnvzf(i,j,kc,g)
                   lmgdz_out(i,j,k3d,g) = lmgdtmp (i,j,kc,g)
                end do
             end do
          end do

          if (k3d .ge. ilo3+1 .and. k3d .le. ihi3+1) then
             do j = ilo2,ihi2
                do i = ilo1,ihi1
                   pggdnvz = 0.5d0*( pgdnvzf(i,j,kc) +  pgdnvzf(i,j,km)) 
                   pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
                        pggdnvz * (ugdnvzf(i,j,kc)-ugdnvzf(i,j,km))/dz
                end do
             end do
          end if

          if (k3d.gt.ilo3) then

             ! Compute U'^z_x and U'^z_y at km (k3d-1) -- note flux3 has physical index
             call transz_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
                  qxm,qmxz,qxp,qpxz,qym,qmyz,qyp,qpyz, &
                  ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  fz, rfz, &
                  ilo1-1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                  ugdnvz, pgdnvz, ergdnvz, &
                  ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  gamcg,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                  cdtdz,ilo1-1,ihi1+1,ilo2-1,ihi2+1,km,kc,k3d)

             ! Compute F^{x|z} at km (k3d-1)
             call cmpflx(qmxz,qpxz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  fxz,  ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                  ugdnvx, pgdnvx, gegdnvx, ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
                  rfxz, ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                  v1gdnvtmp, v2gdnvtmp, ergdnvx, lmgdtmp, &
                  gamcg, &
                  gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                  shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                  1,ilo1,ihi1+1,ilo2-1,ihi2+1,km,km,k3d-1,domlo,domhi)

             ! Compute F^{y|z} at km (k3d-1)
             call cmpflx(qmyz,qpyz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  fyz,  ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                  ugdnvy, pgdnvy, gegdnvy, ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
                  rfyz, ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                  v1gdnvtmp, v2gdnvtmp, ergdnvy, lmgdtmp, &
                  gamcg, &                       
                  gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                  shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &              
                  2,ilo1-1,ihi1+1,ilo2,ihi2+1,km,km,k3d-1,domlo,domhi)

             ! Compute U''_x at km (k3d-1)
             call transyz_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
                  qxm,qxl,qxp,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  fyz, rfyz, &
                  ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                  fzy, rfzy, & 
                  ilo1-1,ilo2,1,ihi1+1,ihi2,2, &
                  ugdnvy, pgdnvy, ergdnvy, &
                  ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  ugdnvtmpz2, pgdnvtmpz2, ergdnvtmpz2, &
                  ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  gamcg,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                  srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                  hdt,hdtdy,hdtdz,ilo1-1,ihi1+1,ilo2,ihi2,km,kc,k3d-1)

             ! Compute U''_y at km (k3d-1)
             call transxz_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
                  qym,qyl,qyp,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  fxz, rfxz, &
                  ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                  fzx, rfzx, &
                  ilo1,ilo2-1,1,ihi1,ihi2+1,2, &
                  ugdnvx, pgdnvx, ergdnvx, &
                  ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  ugdnvtmpz1, pgdnvtmpz1, ergdnvtmpz1, &
                  ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  gamcg,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                  srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                  hdt,hdtdx,hdtdz,ilo1,ihi1,ilo2-1,ihi2+1,km,kc,k3d-1)

             ! Compute F^x at km (k3d-1)
             call cmpflx(qxl, qxr, qd_lo, qd_hi, &
                  flux1, fd1_lo, fd1_hi, &
                  ugdnvxf, pgdnvxf, gegdnvxf, qd_lo, qd_hi, &
                  lam, lam_lo, lam_hi, &
                  rflux1, rfd1_lo, rfd1_hi, &
                  v1gdnvtmp, v2gdnvtmp, ergdnvxf, lmgdtmp, &
                  gamcg, &                           
                  gamc, csml, c, q_lo, q_hi, &
                  shk, shk_lo, shk_hi, &
                  1,lo(1),hi(1)+1,lo(2),hi(2),km,k3d-1,k3d-1,domlo,domhi)

             do j=ilo2-1,ihi2+1
                do i=ilo1-1,ihi1+2
                   ugdnvx_out(i,j,k3d-1) =   ugdnvxf(i,j,km)
                   uy_xfc    (i,j,k3d-1) = v1gdnvtmp(i,j,km)
                   uz_xfc    (i,j,k3d-1) = v2gdnvtmp(i,j,km)
                end do
             end do

             do g=0,ngroups-1
                do j=ilo2-1,ihi2+1
                   do i=ilo1-1,ihi1+2
                      ergdx_out(i,j,k3d-1,g) = ergdnvxf(i,j,km,g)
                      lmgdx_out(i,j,k3d-1,g) = lmgdtmp (i,j,km,g)
                   end do
                end do
             end do

             ! Compute F^y at km (k3d-1)
             call cmpflx(qyl,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  flux2,   fd2_l1, fd2_l2, fd2_l3, fd2_h1, fd2_h2, fd2_h3, &
                  ugdnvyf, pgdnvyf, gegdnvyf, ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                  lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
                  rflux2, rfd2_l1,rfd2_l2,rfd2_l3,rfd2_h1,rfd2_h2,rfd2_h3, &
                  v1gdnvtmp, v2gdnvtmp, ergdnvyf, lmgdtmp, &
                  gamcg, &
                  gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                  shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                  2,ilo1,ihi1,ilo2,ihi2+1,km,k3d-1,k3d-1,domlo,domhi)

             do j=ilo2-1,ihi2+2
                do i=ilo1-1,ihi1+1
                   ugdnvy_out(i,j,k3d-1) =   ugdnvyf(i,j,km)
                   ux_yfc    (i,j,k3d-1) = v1gdnvtmp(i,j,km)
                   uz_yfc    (i,j,k3d-1) = v2gdnvtmp(i,j,km)
                end do
             end do

             do g=0,ngroups-1
                do j=ilo2-1,ihi2+2
                   do i=ilo1-1,ihi1+1
                      ergdy_out(i,j,k3d-1,g) = ergdnvyf(i,j,km,g)
                      lmgdy_out(i,j,k3d-1,g) = lmgdtmp (i,j,km,g)
                   end do
                end do
             end do

             do j = ilo2,ihi2
                do i = ilo1,ihi1
                   pggdnvx = 0.5d0*( pgdnvxf(i+1,j,km) +  pgdnvxf(i,j,km)) 
                   pggdnvy = 0.5d0*( pgdnvyf(i,j+1,km) +  pgdnvyf(i,j,km)) 
                   pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
                        pggdnvx * (ugdnvxf(i+1,j,km)-ugdnvxf(i,j,km))/dx + &
                        pggdnvy * (ugdnvyf(i,j+1,km)-ugdnvyf(i,j,km))/dy
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

  subroutine transy1_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
       qxm,qxmo,qxp,qxpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
       fy, rfy, &
       fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
       ugdnvy, pgdnvy, ergdnvy, &
       pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
       gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
       cdtdy,ilo,ihi,jlo,jhi,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
         QPRES, QREINT, QFA, QFS, QFX, &
         URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX, &
         nadv
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
         fspace_type, comoving
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor

    implicit none

    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
    integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
    integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
    integer ilo,ihi,jlo,jhi,kc,k3d

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)
    double precision  qxm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  qxp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qxmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qxpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  fy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
    double precision rfy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,0:ngroups-1)
    double precision  ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision  pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision ergdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3,0:ngroups-1)
    double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
    double precision cdtdy

    ! Local variables

    integer i, j, g
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
    double precision ugp, ugm, dup, pav, du

    double precision :: pggp, pggm, dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    double precision eddf, f1, ugc

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

    do j = jlo, jhi
       do i = ilo, ihi

          lambda = lam(i,j,k3d,:)
          pggp  =  pgdnvy(i,j+1,kc)
          pggm  =  pgdnvy(i,j  ,kc)
          ugp  =  ugdnvy(i,j+1,kc)
          ugm  =  ugdnvy(i,j  ,kc)
          ugc = 0.5d0*(ugp+ugm)
          ergp = ergdnvy(i,j+1,kc,:)
          ergm = ergdnvy(i,j  ,kc,:)

          ! Convert to conservation form
          rrrx = qxp(i,j,kc,QRHO)
          rurx = rrrx*qxp(i,j,kc,QU)
          rvrx = rrrx*qxp(i,j,kc,QV)
          rwrx = rrrx*qxp(i,j,kc,QW)
          ekenrx = 0.5d0*rrrx*(qxp(i,j,kc,QU)**2 + qxp(i,j,kc,QV)**2 &
               + qxp(i,j,kc,QW)**2)
          rerx = qxp(i,j,kc,QREINT) + ekenrx
          err  = qxp(i,j,kc,qrad:qradhi)

          rrlx = qxm(i+1,j,kc,QRHO)
          rulx = rrlx*qxm(i+1,j,kc,QU)
          rvlx = rrlx*qxm(i+1,j,kc,QV)
          rwlx = rrlx*qxm(i+1,j,kc,QW)
          ekenlx = 0.5d0*rrlx*(qxm(i+1,j,kc,QU)**2 + qxm(i+1,j,kc,QV)**2 &
               + qxm(i+1,j,kc,QW)**2)
          relx = qxm(i+1,j,kc,QREINT) + ekenlx
          erl  = qxm(i+1,j,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewrx = rrrx - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
          runewrx = rurx - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
          rvnewrx = rvrx - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
          rwnewrx = rwrx - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdy*sum(lamge(:))
          rvnewrx = rvnewrx + dmom
          luge = ugc * lamge(:)
          dre = -cdtdy*sum(luge)
          renewrx = rerx - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))  &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if

          ernewr = err(:) - cdtdy*(rfy(i,j+1,kc,:) - rfy(i,j,kc,:)) &
               + der(:)

          rrnewlx = rrlx - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
          runewlx = rulx - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
          rvnewlx = rvlx - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
          rwnewlx = rwlx - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
          rvnewlx = rvnewlx + dmom
          renewlx = relx - cdtdy*(fy(i,j+1,kc,UEDEN)- fy(i,j,kc,UEDEN)) &
               + dre
          ernewl  = erl(:) - cdtdy*(rfy(i,j+1,kc,:) - rfy(i,j,kc,:)) &
               + der(:)

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          pnewrx = qxp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - 1.d0))
          pnewlx = qxm(i+1,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - 1.d0))

          ! Convert back to non-conservation form
          qxpo(i,j,kc,QRHO) = rrnewrx
          qxpo(i,j,kc,QU) = runewrx/qxpo(i,j,kc,QRHO)
          qxpo(i,j,kc,QV) = rvnewrx/qxpo(i,j,kc,QRHO)
          qxpo(i,j,kc,QW) = rwnewrx/qxpo(i,j,kc,QRHO)
          rhoekenrx = 0.5d0*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,kc,QRHO)
          qxpo(i,j,kc,QREINT)= renewrx - rhoekenrx
          qxpo(i,j,kc,QPRES) = pnewrx
          qxpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qxpo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qxpo(i,j,kc,QPRES)
          qxpo(i,j,kc,qreitot) = sum(qxpo(i,j,kc,qrad:qradhi)) + qxpo(i,j,kc,QREINT)

          qxmo(i+1,j,kc,QRHO) = rrnewlx
          qxmo(i+1,j,kc,QU) = runewlx/qxmo(i+1,j,kc,QRHO)
          qxmo(i+1,j,kc,QV) = rvnewlx/qxmo(i+1,j,kc,QRHO)
          qxmo(i+1,j,kc,QW) = rwnewlx/qxmo(i+1,j,kc,QRHO)
          rhoekenlx = 0.5d0*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,kc,QRHO)
          qxmo(i+1,j,kc,QREINT)= renewlx - rhoekenlx
          qxmo(i+1,j,kc,QPRES) = pnewlx
          qxmo(i+1,j,kc,qrad:qradhi) = ernewl(:)
          qxmo(i+1,j,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qxmo(i+1,j,kc,QPRES)
          qxmo(i+1,j,kc,qreitot) = sum(qxmo(i+1,j,kc,qrad:qradhi)) + qxmo(i+1,j,kc,QREINT)

       enddo
    enddo

  end subroutine transy1_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transx1_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
       qym,qymo,qyp,qypo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
       fx, rfx, &
       fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
       ugdnvx, pgdnvx, ergdnvx, & 
       pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
       gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
       cdtdx,ilo,ihi,jlo,jhi,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
         QPRES, QREINT, QFA, QFS, QFX, &
         URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX, &
         nadv
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
         fspace_type, comoving
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor

    implicit none

    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
    integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
    integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
    integer ilo,ihi,jlo,jhi,kc,k3d

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)
    double precision  qym(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  qyp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qymo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qypo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  fx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
    double precision rfx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,0:ngroups-1)
    double precision  ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision  pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision ergdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3,0:ngroups-1)
    double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
    double precision cdtdx

    ! Local variables

    integer i, j, g
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
    double precision ugp, ugm, dup, pav, du

    double precision :: pggp, pggm, dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    double precision eddf, f1, ugc

    ! NOTE: it is better *not* to protect against small density in this routine

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

    do j = jlo, jhi 
       do i = ilo, ihi 

          lambda = lam(i,j,k3d,:)
          pggp  =  pgdnvx(i+1,j,kc)
          pggm  =  pgdnvx(i  ,j,kc)
          ugp  =  ugdnvx(i+1,j,kc)
          ugm  =  ugdnvx(i  ,j,kc)
          ugc = 0.5d0*(ugp+ugm)
          ergp = ergdnvx(i+1,j,kc,:)
          ergm = ergdnvx(i  ,j,kc,:)

          ! Convert to conservation form
          rrry = qyp(i,j,kc,QRHO)
          rury = rrry*qyp(i,j,kc,QU)
          rvry = rrry*qyp(i,j,kc,QV)
          rwry = rrry*qyp(i,j,kc,QW)
          ekenry = 0.5d0*rrry*(qyp(i,j,kc,QU)**2 + qyp(i,j,kc,QV)**2 + qyp(i,j,kc,QW)**2)
          rery = qyp(i,j,kc,QREINT) + ekenry
          err  = qyp(i,j,kc,qrad:qradhi)

          rrly = qym(i,j+1,kc,QRHO)
          ruly = rrly*qym(i,j+1,kc,QU)
          rvly = rrly*qym(i,j+1,kc,QV)
          rwly = rrly*qym(i,j+1,kc,QW)
          ekenly = 0.5d0*rrly* &
               (qym(i,j+1,kc,QU)**2 + qym(i,j+1,kc,QV)**2 + qym(i,j+1,kc,QW)**2)
          rely = qym(i,j+1,kc,QREINT) + ekenly
          erl  = qym(i,j+1,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewry = rrry - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
          runewry = rury - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
          rvnewry = rvry - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
          rwnewry = rwry - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdx*sum(lamge(:))
          runewry = runewry + dmom
          luge = ugc * lamge(:)
          dre = -cdtdx*sum(luge)
          renewry = rery - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdx * luge
          end if

          ernewr  = err(:) - cdtdx*(rfx(i+1,j,kc,:) - rfx(i,j,kc,:)) &
               + der(:)

          rrnewly = rrly - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
          runewly = ruly - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
          rvnewly = rvly - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
          rwnewly = rwly - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
          runewly = runewly + dmom
          renewly = rely - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN)) &
               + dre
          ernewl  = erl(:) - cdtdx*(rfx(i+1,j,kc,:) - rfx(i,j,kc,:)) &
               + der(:)

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          pnewry = qyp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-1.d0))
          pnewly = qym(i,j+1,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-1.d0))

          ! Convert back to non-conservation form
          qypo(i,j,kc,QRHO) = rrnewry
          qypo(i,j,kc,QU) = runewry/qypo(i,j,kc,QRHO)
          qypo(i,j,kc,QV) = rvnewry/qypo(i,j,kc,QRHO)
          qypo(i,j,kc,QW) = rwnewry/qypo(i,j,kc,QRHO)
          rhoekenry = 0.5d0*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,kc,QRHO)
          qypo(i,j,kc,QREINT) = renewry - rhoekenry
          qypo(i,j,kc,QPRES) = pnewry
          qypo(i,j,kc,qrad:qradhi) = ernewr(:)
          qypo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qypo(i,j,kc,QPRES)
          qypo(i,j,kc,qreitot) = sum(qypo(i,j,kc,qrad:qradhi)) + qypo(i,j,kc,QREINT)

          qymo(i,j+1,kc,QRHO) = rrnewly
          qymo(i,j+1,kc,QU) = runewly/qymo(i,j+1,kc,QRHO)
          qymo(i,j+1,kc,QV) = rvnewly/qymo(i,j+1,kc,QRHO)
          qymo(i,j+1,kc,QW) = rwnewly/qymo(i,j+1,kc,QRHO)
          rhoekenly = 0.5d0*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,kc,QRHO)
          qymo(i,j+1,kc,QREINT) = renewly - rhoekenly
          qymo(i,j+1,kc,QPRES) = pnewly
          qymo(i,j+1,kc,qrad:qradhi) = ernewl(:)
          qymo(i,j+1,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qymo(i,j+1,kc,QPRES)
          qymo(i,j+1,kc,qreitot) = sum(qymo(i,j+1,kc,qrad:qradhi)) + qymo(i,j+1,kc,QREINT)

       enddo
    enddo

  end subroutine transx1_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transy2_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
       qzm,qzmo,qzp,qzpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
       fy, rfy, &
       fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
       ugdnvy, pgdnvy, ergdnvy, & 
       pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
       gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
       cdtdy,ilo,ihi,jlo,jhi,kc,km,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
         QPRES, QREINT, QFA, QFS, QFX, &
         URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX, &
         nadv
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
         fspace_type, comoving
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor

    implicit none

    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
    integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
    integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
    integer ilo,ihi,jlo,jhi,kc,km,k3d

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)
    double precision  qzm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  qzp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qzmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qzpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  fy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
    double precision rfy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,0:ngroups-1)
    double precision  ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision  pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision ergdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3,0:ngroups-1)
    double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
    double precision cdtdy

    ! Local variables

    integer i, j, g
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
    double precision ugp, ugm, dup, pav, du

    double precision :: pggp, pggm, dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    double precision eddf, f1, ugc

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

    do j = jlo, jhi
       do i = ilo, ihi

          lambda = lam(i,j,k3d,:)
          pggp  =  pgdnvy(i,j+1,kc)
          pggm  =  pgdnvy(i,j  ,kc)
          ugp  =  ugdnvy(i,j+1,kc)
          ugm  =  ugdnvy(i,j  ,kc)
          ugc = 0.5d0*(ugp+ugm)
          ergp = ergdnvy(i,j+1,kc,:)
          ergm = ergdnvy(i,j  ,kc,:)

          ! Convert to conservation form
          rrrz = qzp(i,j,kc,QRHO)
          rurz = rrrz*qzp(i,j,kc,QU)
          rvrz = rrrz*qzp(i,j,kc,QV)
          rwrz = rrrz*qzp(i,j,kc,QW)
          ekenrz = 0.5d0*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 &
               + qzp(i,j,kc,QW)**2)
          rerz = qzp(i,j,kc,QREINT) + ekenrz
          err  = qzp(i,j,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewrz = rrrz - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
          runewrz = rurz - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
          rvnewrz = rvrz - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
          rwnewrz = rwrz - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdy*sum(lamge(:))
          rvnewrz = rvnewrz + dmom
          luge = 0.5d0*(ugp+ugm) * lamge(:)
          dre = -cdtdy*sum(luge)
          renewrz = rerz - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if

          ernewr  = err(:) - cdtdy*(rfy(i,j+1,kc,:) - rfy(i,j,kc,:)) &
               + der(:)

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          pnewrz = qzp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - 1.d0))

          ! Convert back to non-conservation form
          qzpo(i,j,kc,QRHO) = rrnewrz
          qzpo(i,j,kc,QU) = runewrz/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QV) = rvnewrz/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QW) = rwnewrz/qzpo(i,j,kc,QRHO)
          rhoekenrz = 0.5d0*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QREINT)= renewrz - rhoekenrz
          qzpo(i,j,kc,QPRES) = pnewrz
          qzpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qzpo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,kc,QPRES)
          qzpo(i,j,kc,qreitot) = sum(qzpo(i,j,kc,qrad:qradhi)) + qzpo(i,j,kc,QREINT)

          lambda = lam(i,j,k3d-1,:)
          pggp  =  pgdnvy(i,j+1,km)
          pggm  =  pgdnvy(i,j  ,km)
          ugp  =  ugdnvy(i,j+1,km)
          ugm  =  ugdnvy(i,j  ,km)
          ugc = 0.5d0*(ugp+ugm)
          ergp = ergdnvy(i,j+1,km,:)
          ergm = ergdnvy(i,j  ,km,:)

          rrlz = qzm(i,j,kc,QRHO)
          rulz = rrlz*qzm(i,j,kc,QU)
          rvlz = rrlz*qzm(i,j,kc,QV)
          rwlz = rrlz*qzm(i,j,kc,QW)
          ekenlz = 0.5d0*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 &
               + qzm(i,j,kc,QW)**2)
          relz = qzm(i,j,kc,QREINT) + ekenlz
          erl  = qzm(i,j,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewlz = rrlz - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
          runewlz = rulz - cdtdy*(fy(i,j+1,km,UMX) - fy(i,j,km,UMX))
          rvnewlz = rvlz - cdtdy*(fy(i,j+1,km,UMY) - fy(i,j,km,UMY))
          rwnewlz = rwlz - cdtdy*(fy(i,j+1,km,UMZ) - fy(i,j,km,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdy*sum(lamge(:))
          rvnewlz = rvnewlz + dmom
          luge = 0.5d0*(ugp+ugm) * lamge(:)
          dre = -cdtdy*sum(luge)
          renewlz = relz - cdtdy*(fy(i,j+1,km,UEDEN)- fy(i,j,km,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if

          ernewl  = erl(:) - cdtdy*(rfy(i,j+1,km,:)- rfy(i,j,km,:)) &
               + der

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          pnewlz = qzm(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))

          qzmo(i,j,kc,QRHO) = rrnewlz
          qzmo(i,j,kc,QU) = runewlz/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QV) = rvnewlz/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QW) = rwnewlz/qzmo(i,j,kc,QRHO)
          rhoekenlz = 0.5d0*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QREINT)= renewlz - rhoekenlz
          qzmo(i,j,kc,QPRES) = pnewlz
          qzmo(i,j,kc,qrad:qradhi) = ernewl(:)
          qzmo(i,j,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,kc,QPRES)
          qzmo(i,j,kc,qreitot) = sum(qzmo(i,j,kc,qrad:qradhi)) + qzmo(i,j,kc,QREINT)

       enddo
    enddo

  end subroutine transy2_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transx2_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
       qzm,qzmo,qzp,qzpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
       fx, rfx, &
       fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
       ugdnvx, pgdnvx, ergdnvx, & 
       pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
       gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
       cdtdx,ilo,ihi,jlo,jhi,kc,km,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
         QPRES, QREINT, QFA, QFS, QFX, &
         URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX, &
         nadv
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
         fspace_type, comoving
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor

    implicit none

    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
    integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
    integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
    integer ilo,ihi,jlo,jhi,kc,km,k3d

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)
    double precision  qzm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  qzp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qzmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qzpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  fx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
    double precision rfx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,0:ngroups-1)
    double precision  ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision  pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision ergdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3,0:ngroups-1)
    double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
    double precision cdtdx

    ! Local variables

    integer i, j, g
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
    double precision ugp, ugm, dup, pav, du

    double precision :: pggp, pggm, dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    double precision eddf, f1, ugc

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

    do j = jlo, jhi 
       do i = ilo, ihi 

          lambda = lam(i,j,k3d,:)
          pggp  =  pgdnvx(i+1,j,kc)
          pggm  =  pgdnvx(i  ,j,kc)
          ugp  =  ugdnvx(i+1,j,kc)
          ugm  =  ugdnvx(i  ,j,kc)
          ugc = 0.5d0*(ugp+ugm)
          ergp = ergdnvx(i+1,j,kc,:)
          ergm = ergdnvx(i  ,j,kc,:)

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          ! Convert to conservation form
          rrrz = qzp(i,j,kc,QRHO)
          rurz = rrrz*qzp(i,j,kc,QU)
          rvrz = rrrz*qzp(i,j,kc,QV)
          rwrz = rrrz*qzp(i,j,kc,QW)
          ekenrz = 0.5d0*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 + qzp(i,j,kc,QW)**2)
          rerz = qzp(i,j,kc,QREINT) + ekenrz
          err  = qzp(i,j,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewrz = rrrz - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
          runewrz = rurz - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
          rvnewrz = rvrz - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
          rwnewrz = rwrz - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdx*sum(lamge(:))
          runewrz = runewrz + dmom
          luge = 0.5d0*(ugp+ugm) * lamge(:)
          dre = -cdtdx*sum(luge)
          renewrz = rerz - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdx * luge
          end if

          ernewr  = err(:) - cdtdx*(rfx(i+1,j,kc,:) - rfx(i,j,kc,:)) &
               + der(:)


          pnewrz = qzp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-1.d0))

          ! Convert back to non-conservation form
          qzpo(i,j,kc,QRHO) = rrnewrz
          qzpo(i,j,kc,QU) = runewrz/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QV) = rvnewrz/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QW) = rwnewrz/qzpo(i,j,kc,QRHO)
          rhoekenrz = 0.5d0*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QREINT) = renewrz - rhoekenrz
          qzpo(i,j,kc,QPRES) = pnewrz
          qzpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qzpo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,kc,QPRES)
          qzpo(i,j,kc,qreitot) = sum(qzpo(i,j,kc,qrad:qradhi)) + qzpo(i,j,kc,QREINT)

          lambda = lam(i,j,k3d-1,:)
          pggp  =  pgdnvx(i+1,j,km)
          pggm  =  pgdnvx(i  ,j,km)
          ugp  =  ugdnvx(i+1,j,km)
          ugm  =  ugdnvx(i  ,j,km)
          ugc = 0.5d0*(ugp+ugm)
          ergp = ergdnvx(i+1,j,km,:)
          ergm = ergdnvx(i  ,j,km,:)

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          rrlz = qzm(i,j,kc,QRHO)
          rulz = rrlz*qzm(i,j,kc,QU)
          rvlz = rrlz*qzm(i,j,kc,QV)
          rwlz = rrlz*qzm(i,j,kc,QW)
          ekenlz = 0.5d0*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 + qzm(i,j,kc,QW)**2)
          relz = qzm(i,j,kc,QREINT) + ekenlz
          erl  = qzm(i,j,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewlz = rrlz - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
          runewlz = rulz - cdtdx*(fx(i+1,j,km,UMX) - fx(i,j,km,UMX))
          rvnewlz = rvlz - cdtdx*(fx(i+1,j,km,UMY) - fx(i,j,km,UMY))
          rwnewlz = rwlz - cdtdx*(fx(i+1,j,km,UMZ) - fx(i,j,km,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdx*sum(lamge(:))
          runewlz = runewlz + dmom
          luge = 0.5d0*(ugp+ugm) * lamge(:)
          dre = -cdtdx*sum(luge)
          renewlz = relz - cdtdx*(fx(i+1,j,km,UEDEN) - fx(i,j,km,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdx * luge
          end if

          ernewl  = erl(:) - cdtdx*(rfx(i+1,j,km,:) - rfx(i,j,km,:)) &
               +der(:)

          pnewlz = qzm(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d-1)-1.d0))

          qzmo(i,j,kc,QRHO) = rrnewlz
          qzmo(i,j,kc,QU) = runewlz/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QV) = rvnewlz/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QW) = rwnewlz/qzmo(i,j,kc,QRHO)
          rhoekenlz = 0.5d0*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QREINT) = renewlz - rhoekenlz
          qzmo(i,j,kc,QPRES) = pnewlz
          qzmo(i,j,kc,qrad:qradhi) = ernewl(:)
          qzmo(i,j,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,kc,QPRES)
          qzmo(i,j,kc,qreitot) = sum(qzmo(i,j,kc,qrad:qradhi)) + qzmo(i,j,kc,QREINT)

       enddo
    enddo

  end subroutine transx2_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transxy_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
       qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
       fxy, rfxy, &
       fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
       fyx, rfyx, &
       fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
       ugdnvx, pgdnvx, ergdnvx, &
       pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
       ugdnvy, pgdnvy, ergdnvy, &
       pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
       gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
       srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
       hdt,cdtdx,cdtdy,ilo,ihi,jlo,jhi,kc,km,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
         QPRES, QREINT, QFA, QFS, QFX, &
         URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX, &
         nadv
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
         fspace_type, comoving
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor

    implicit none

    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
    integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
    integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
    integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
    integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer ilo,ihi,jlo,jhi,km,kc,k3d

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)
    double precision  qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  fxy(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
    double precision rfxy(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,0:ngroups-1)
    double precision  fyx(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
    double precision rfyx(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,0:ngroups-1)
    double precision  ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision  pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision ergdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3,0:ngroups-1)
    double precision  ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision  pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision ergdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3,0:ngroups-1)
    double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
    double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
    double precision hdt,cdtdx,cdtdy

    ! Local variables

    integer i, j, g
    integer n , nq
    integer iadv, ispec, iaux

    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision ugxp, ugxm, duxp, pxav, dux, pxnew
    double precision ugyp, ugym, duyp, pyav, duy, pynew
    double precision ugxpm, ugxmm, duxpm, pxavm, duxm, pxnewm
    double precision ugypm, ugymm, duypm, pyavm, duym, pynewm
    double precision compr, compl, compnr, compnl
    double precision rhotmp

    double precision :: dmx, dmy, dre, pggxp, pggyp, pggxm, pggym 
    double precision :: pggxpm, pggypm, pggxmm, pggymm 
    double precision, dimension(0:ngroups-1) :: der, lamc, lamm, lugex, lugey, lgex, lgey, &
         err, ernewr, erl, ernewl, ergxp, ergyp, ergxm, ergym, ergxpm, ergypm, ergxmm, ergymm
    double precision eddf, f1

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

    do j = jlo, jhi 
       do i = ilo, ihi 

          lamc = lam(i,j,k3d,:)

          pggxp  =  pgdnvx(i+1,j,kc)
          pggxm  =  pgdnvx(i  ,j,kc)
          ugxp  =  ugdnvx(i+1,j,kc)
          ugxm  =  ugdnvx(i  ,j,kc)
          ergxp = ergdnvx(i+1,j,kc,:)
          ergxm = ergdnvx(i  ,j,kc,:)

          pggyp  =  pgdnvy(i,j+1,kc)
          pggym  =  pgdnvy(i,j  ,kc)
          ugyp  =  ugdnvy(i,j+1,kc)
          ugym  =  ugdnvy(i,j  ,kc)
          ergyp = ergdnvy(i,j+1,kc,:)
          ergym = ergdnvy(i,j  ,kc,:)

          lamm = lam(i,j,k3d-1,:)

          pggxpm  =  pgdnvx(i+1,j,km)
          pggxmm  =  pgdnvx(i  ,j,km)
          ugxpm  =  ugdnvx(i+1,j,km)
          ugxmm  =  ugdnvx(i  ,j,km)
          ergxpm = ergdnvx(i+1,j,km,:)
          ergxmm = ergdnvx(i  ,j,km,:)

          pggypm  =  pgdnvy(i,j+1,km)
          pggymm  =  pgdnvy(i,j  ,km)
          ugypm  =  ugdnvy(i,j+1,km)
          ugymm  =  ugdnvy(i,j  ,km)
          ergypm = ergdnvy(i,j+1,km,:)
          ergymm = ergdnvy(i,j  ,km,:)

          ! Convert to conservation form
          rrr = qp(i,j,kc,QRHO)
          rur = rrr*qp(i,j,kc,QU)
          rvr = rrr*qp(i,j,kc,QV)
          rwr = rrr*qp(i,j,kc,QW)
          ekenr = 0.5d0*rrr*(qp(i,j,kc,QU)**2 + qp(i,j,kc,QV)**2 + &
               qp(i,j,kc,QW)**2)
          rer = qp(i,j,kc,QREINT) + ekenr
          err = qp(i,j,kc,qrad:qradhi)

          rrl = qm(i,j,kc,QRHO)
          rul = rrl*qm(i,j,kc,QU)
          rvl = rrl*qm(i,j,kc,QV)
          rwl = rrl*qm(i,j,kc,QW)
          ekenl = 0.5d0*rrl*(qm(i,j,kc,QU)**2 + qm(i,j,kc,QV)**2 + &
               qm(i,j,kc,QW)**2)
          rel = qm(i,j,kc,QREINT) + ekenl
          erl = qm(i,j,kc,qrad:qradhi)

          ! Add transverse predictor
          rrnewr = rrr - cdtdx*(fxy(i+1,j  ,kc,URHO) - fxy(i,j,kc,URHO)) &
               - cdtdy*(fyx(i  ,j+1,kc,URHO) - fyx(i,j,kc,URHO))
          runewr = rur - cdtdx*(fxy(i+1,j  ,kc,UMX ) - fxy(i,j,kc,UMX)) &
               - cdtdy*(fyx(i  ,j+1,kc,UMX ) - fyx(i,j,kc,UMX))
          rvnewr = rvr - cdtdx*(fxy(i+1,j  ,kc,UMY ) - fxy(i,j,kc,UMY)) &
               - cdtdy*(fyx(i  ,j+1,kc,UMY ) - fyx(i,j,kc,UMY))
          rwnewr = rwr - cdtdx*(fxy(i+1,j  ,kc,UMZ ) - fxy(i,j,kc,UMZ)) &
               - cdtdy*(fyx(i  ,j+1,kc,UMZ ) - fyx(i,j,kc,UMZ))
          lgex = lamc(:) * (ergxp(:)-ergxm(:))
          lgey = lamc(:) * (ergyp(:)-ergym(:))
          dmx = - cdtdx*sum(lgex)
          dmy = - cdtdy*sum(lgey)
          runewr = runewr + dmx
          rvnewr = rvnewr + dmy
          lugex = 0.5d0*(ugxp+ugxm) * lgex(:)
          lugey = 0.5d0*(ugyp+ugym) * lgey(:)
          dre = -cdtdx*sum(lugex) - cdtdy*sum(lugey)
          renewr = rer - cdtdx*(fxy(i+1,j  ,kc,UEDEN) - fxy(i,j,kc,UEDEN)) &
               - cdtdy*(fyx(i  ,j+1,kc,UEDEN) - fyx(i,j,kc,UEDEN))  &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamc(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdx*0.5d0*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                     +       cdtdy*0.5d0*(ugyp+ugym)*(ergyp(g)-ergym(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamc(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdx*0.5d0*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                     +       cdtdy*0.5d0*(ergyp(g)+ergym(g))*(ugym-ugyp) )
             end do
          else ! mixed frame
             der(:) = cdtdx * lugex + cdtdy * lugey
          end if

          ernewr = err(:) - cdtdx*(rfxy(i+1,j,kc,:) - rfxy(i,j,kc,:)) &
               &          - cdtdy*(rfyx(i,j+1,kc,:) - rfyx(i,j,kc,:))  &
               &          + der(:)

          rrnewl = rrl - cdtdx*(fxy(i+1,j  ,km,URHO) - fxy(i,j,km,URHO)) &
               - cdtdy*(fyx(i  ,j+1,km,URHO) - fyx(i,j,km,URHO))
          runewl = rul - cdtdx*(fxy(i+1,j  ,km,UMX ) - fxy(i,j,km,UMX)) &
               - cdtdy*(fyx(i  ,j+1,km,UMX ) - fyx(i,j,km,UMX))
          rvnewl = rvl - cdtdx*(fxy(i+1,j  ,km,UMY ) - fxy(i,j,km,UMY)) &
               - cdtdy*(fyx(i  ,j+1,km,UMY ) - fyx(i,j,km,UMY))
          rwnewl = rwl - cdtdx*(fxy(i+1,j  ,km,UMZ ) - fxy(i,j,km,UMZ)) &
               - cdtdy*(fyx(i  ,j+1,km,UMZ ) - fyx(i,j,km,UMZ))
          lgex = lamm(:) * (ergxpm(:)-ergxmm(:))
          lgey = lamm(:) * (ergypm(:)-ergymm(:))
          dmx = - cdtdx*sum(lgex)
          dmy = - cdtdy*sum(lgey)
          runewl = runewl + dmx
          rvnewl = rvnewl + dmy 
          lugex = 0.5d0*(ugxpm+ugxmm) * lgex(:)
          lugey = 0.5d0*(ugypm+ugymm) * lgey(:)
          dre = -cdtdx*sum(lugex) - cdtdy*sum(lugey)
          renewl = rel - cdtdx*(fxy(i+1,j  ,km,UEDEN) - fxy(i,j,km,UEDEN)) &
               - cdtdy*(fyx(i  ,j+1,km,UEDEN) - fyx(i,j,km,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamm(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdx*0.5d0*(ugxpm+ugxmm)*(ergxpm(g)-ergxmm(g)) &
                     +       cdtdy*0.5d0*(ugypm+ugymm)*(ergypm(g)-ergymm(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamm(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdx*0.5d0*(ergxpm(g)+ergxmm(g))*(ugxmm-ugxpm) &
                     +       cdtdy*0.5d0*(ergypm(g)+ergymm(g))*(ugymm-ugypm) )
             end do
          else ! mixed frame
             der(:) = cdtdx * lugex + cdtdy * lugey
          end if

          ernewl = erl(:) - cdtdx*(rfxy(i+1,j  ,km,:) - rfxy(i,j,km,:)) &
               &          - cdtdy*(rfyx(i  ,j+1,km,:) - rfyx(i,j,km,:)) &
               &          + der(:)

          duxp = pggxp*ugxp - pggxm*ugxm
          pxav = 0.5d0*(pggxp+pggxm)
          dux = ugxp-ugxm
          pxnew = cdtdx*(duxp + pxav*dux*(gamc(i,j,k3d)-1.d0))

          duxpm = pggxpm*ugxpm - pggxmm*ugxmm
          pxavm = 0.5d0*(pggxpm+pggxmm)
          duxm = ugxpm-ugxmm
          pxnewm = cdtdx*(duxpm + pxavm*duxm*(gamc(i,j,k3d-1)-1.d0))

          duyp = pggyp*ugyp - pggym*ugym
          pyav = 0.5d0*(pggyp+pggym)
          duy = ugyp-ugym
          pynew = cdtdy*(duyp + pyav*duy*(gamc(i,j,k3d)-1.d0))

          duypm = pggypm*ugypm - pggymm*ugymm
          pyavm = 0.5d0*(pggypm+pggymm)
          duym = ugypm-ugymm
          pynewm = cdtdy*(duypm + pyavm*duym*(gamc(i,j,k3d-1)-1.d0))

          pnewr = qp(i,j,kc,QPRES) - pxnew - pynew
          pnewl = qm(i,j,kc,QPRES) - pxnewm - pynewm

          ! Convert back to non-conservation form
          rhotmp = rrnewr
          qpo(i,j,kc,QRHO  ) = rhotmp        + hdt*srcQ(i,j,k3d,QRHO)
          qpo(i,j,kc,QU    ) = runewr/rhotmp + hdt*srcQ(i,j,k3d,QU)
          qpo(i,j,kc,QV    ) = rvnewr/rhotmp + hdt*srcQ(i,j,k3d,QV)
          qpo(i,j,kc,QW    ) = rwnewr/rhotmp + hdt*srcQ(i,j,k3d,QW)
          rhoekenr = 0.5d0*(runewr**2 + rvnewr**2 + rwnewr**2)/rhotmp
          qpo(i,j,kc,QREINT) = renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)
          qpo(i,j,kc,QPRES ) = pnewr             + hdt*srcQ(i,j,k3d,QPRES)
          qpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qpo(i,j,kc,qptot  ) = sum(lamc(:)*ernewr(:)) + qpo(i,j,kc,QPRES)
          qpo(i,j,kc,qreitot) = sum(qpo(i,j,kc,qrad:qradhi)) + qpo(i,j,kc,QREINT)

          rhotmp = rrnewl
          qmo(i,j,kc,QRHO  ) = rhotmp        + hdt*srcQ(i,j,k3d-1,QRHO)
          qmo(i,j,kc,QU    ) = runewl/rhotmp + hdt*srcQ(i,j,k3d-1,QU)
          qmo(i,j,kc,QV    ) = rvnewl/rhotmp + hdt*srcQ(i,j,k3d-1,QV)
          qmo(i,j,kc,QW    ) = rwnewl/rhotmp + hdt*srcQ(i,j,k3d-1,QW)
          rhoekenl = 0.5d0*(runewl**2 + rvnewl**2 + rwnewl**2)/rhotmp
          qmo(i,j,kc,QREINT) = renewl - rhoekenl + hdt*srcQ(i,j,k3d-1,QREINT)
          qmo(i,j,kc,QPRES ) = pnewl             + hdt*srcQ(i,j,k3d-1,QPRES)
          qmo(i,j,kc,qrad:qradhi) = ernewl(:)
          qmo(i,j,kc,qptot  ) = sum(lamm(:)*ernewl(:)) + qmo(i,j,kc,QPRES)
          qmo(i,j,kc,qreitot) = sum(qmo(i,j,kc,qrad:qradhi)) + qmo(i,j,kc,QREINT)

       enddo
    enddo

  end subroutine transxy_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transz_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
       qxm,qxmo,qxp,qxpo,qym,qymo,qyp,qypo, &
       qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
       fz, rfz, &
       fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
       ugdnvz, pgdnvz, ergdnvz, &
       pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
       gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
       cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
         QPRES, QREINT, QFA, QFS, QFX,&
         URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, UFX, &
         nadv
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
         fspace_type, comoving
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor

    implicit none

    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
    integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
    integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
    integer ilo,ihi,jlo,jhi,km,kc,k3d

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)
    double precision  qxm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  qxp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  qym(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  qyp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qxmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qxpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qymo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qypo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  fz(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
    double precision rfz(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,0:ngroups-1)
    double precision  ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision  pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision ergdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3,0:ngroups-1)
    double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
    double precision cdtdz

    ! Local variables

    integer n, nq
    integer iadv, ispec, iaux
    integer i, j, g

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
    double precision ugp, ugm, dup, pav, du

    double precision :: dmz, dre, pggp, pggm
    double precision, dimension(0:ngroups-1) :: der, lambda, luge, lamge, &
         ergp, errx, ernewrx, erry, ernewry, ergm, erlx, ernewlx, erly, ernewly
    double precision eddf, f1

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

    do j = jlo, jhi 
       do i = ilo, ihi 

          lambda = lam(i,j,k3d-1,:)

          pggp  =  pgdnvz(i,j,kc)
          pggm  =  pgdnvz(i,j,km)
          ugp  =  ugdnvz(i,j,kc)
          ugm  =  ugdnvz(i,j,km)
          ergp = ergdnvz(i,j,kc,:)
          ergm = ergdnvz(i,j,km,:)

          ! Convert to conservation form
          rrrx = qxp(i,j,km,QRHO)
          rurx = rrrx*qxp(i,j,km,QU)
          rvrx = rrrx*qxp(i,j,km,QV)
          rwrx = rrrx*qxp(i,j,km,QW)
          ekenrx = 0.5d0*rrrx*(qxp(i,j,km,QU)**2 + qxp(i,j,km,QV)**2 &
               + qxp(i,j,km,QW)**2)
          rerx = qxp(i,j,km,QREINT) + ekenrx
          errx = qxp(i,j,km,qrad:qradhi)

          rrry = qyp(i,j,km,QRHO)
          rury = rrry*qyp(i,j,km,QU)
          rvry = rrry*qyp(i,j,km,QV)
          rwry = rrry*qyp(i,j,km,QW)
          ekenry = 0.5d0*rrry*(qyp(i,j,km,QU)**2 + qyp(i,j,km,QV)**2 &
               + qyp(i,j,km,QW)**2)
          rery = qyp(i,j,km,QREINT) + ekenry
          erry = qyp(i,j,km,qrad:qradhi)

          rrlx = qxm(i+1,j,km,QRHO)
          rulx = rrlx*qxm(i+1,j,km,QU)
          rvlx = rrlx*qxm(i+1,j,km,QV)
          rwlx = rrlx*qxm(i+1,j,km,QW)
          ekenlx = 0.5d0*rrlx*(qxm(i+1,j,km,QU)**2 + qxm(i+1,j,km,QV)**2 &
               + qxm(i+1,j,km,QW)**2)
          relx = qxm(i+1,j,km,QREINT) + ekenlx
          erlx = qxm(i+1,j,km,qrad:qradhi)

          rrly = qym(i,j+1,km,QRHO)
          ruly = rrly*qym(i,j+1,km,QU)
          rvly = rrly*qym(i,j+1,km,QV)
          rwly = rrly*qym(i,j+1,km,QW)
          ekenly = 0.5d0*rrly*(qym(i,j+1,km,QU)**2 + qym(i,j+1,km,QV)**2 &
               + qym(i,j+1,km,QW)**2)
          rely = qym(i,j+1,km,QREINT) + ekenly
          erly = qym(i,j+1,km,qrad:qradhi)

          ! Add transverse predictor
          rrnewrx = rrrx - cdtdz*(fz(i,j,kc,URHO ) - fz(i,j,km,URHO))
          runewrx = rurx - cdtdz*(fz(i,j,kc,UMX  ) - fz(i,j,km,UMX))
          rvnewrx = rvrx - cdtdz*(fz(i,j,kc,UMY  ) - fz(i,j,km,UMY))
          rwnewrx = rwrx - cdtdz*(fz(i,j,kc,UMZ  ) - fz(i,j,km,UMZ))
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmz = - cdtdz*sum(lamge)
          rwnewrx = rwnewrx + dmz
          luge = 0.5d0*(ugp+ugm) * lamge(:)
          dre = -cdtdz*sum(luge)
          renewrx = rerx - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdz*0.5d0*(ugp+ugm)*(ergp(g)-ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdz*0.5d0*(ergp(g)+ergm(g))*(ugm-ugp)*f1
             end do
          else ! mixed frame
             der(:) = cdtdz * luge
          end if

          ernewrx = errx(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
               + der(:)

          rrnewry = rrry - cdtdz*(fz(i,j,kc,URHO ) - fz(i,j,km,URHO))
          runewry = rury - cdtdz*(fz(i,j,kc,UMX  ) - fz(i,j,km,UMX))
          rvnewry = rvry - cdtdz*(fz(i,j,kc,UMY  ) - fz(i,j,km,UMY))
          rwnewry = rwry - cdtdz*(fz(i,j,kc,UMZ  ) - fz(i,j,km,UMZ))
          rwnewry = rwnewry + dmz
          renewry = rery - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN)) &
               + dre
          ernewry = erry(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
               + der(:)

          rrnewlx = rrlx - cdtdz*(fz(i,j,kc,URHO ) - fz(i,j,km,URHO))
          runewlx = rulx - cdtdz*(fz(i,j,kc,UMX  ) - fz(i,j,km,UMX))
          rvnewlx = rvlx - cdtdz*(fz(i,j,kc,UMY  ) - fz(i,j,km,UMY))
          rwnewlx = rwlx - cdtdz*(fz(i,j,kc,UMZ  ) - fz(i,j,km,UMZ))
          rwnewlx = rwnewlx + dmz
          renewlx = relx - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN)) &
               + dre
          ernewlx = erlx(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
               + der(:)

          rrnewly = rrly - cdtdz*(fz(i,j,kc,URHO ) - fz(i,j,km,URHO))
          runewly = ruly - cdtdz*(fz(i,j,kc,UMX  ) - fz(i,j,km,UMX))
          rvnewly = rvly - cdtdz*(fz(i,j,kc,UMY  ) - fz(i,j,km,UMY))
          rwnewly = rwly - cdtdz*(fz(i,j,kc,UMZ  ) - fz(i,j,km,UMZ))
          rwnewly = rwnewly + dmz
          renewly = rely - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN)) &
               + dre
          ernewly = erly(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
               + der(:)

          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

          pnewrx = qxp(i  ,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))
          pnewlx = qxm(i+1,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))

          pnewry = qyp(i,j  ,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))
          pnewly = qym(i,j+1,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - 1.d0))

          ! Convert back to non-conservation form
          qxpo(i,j,km,QRHO) = rrnewrx
          qxpo(i,j,km,QU) = runewrx/qxpo(i,j,km,QRHO)
          qxpo(i,j,km,QV) = rvnewrx/qxpo(i,j,km,QRHO)
          qxpo(i,j,km,QW) = rwnewrx/qxpo(i,j,km,QRHO)
          rhoekenrx = 0.5d0*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,km,QRHO)
          qxpo(i,j,km,QREINT)= renewrx - rhoekenrx
          qxpo(i,j,km,QPRES) = pnewrx
          qxpo(i,j,km,qrad:qradhi) = ernewrx(:)
          qxpo(i,j,km,qptot  ) = sum(lambda(:)*ernewrx(:)) + qxpo(i,j,km,QPRES)
          qxpo(i,j,km,qreitot) = sum(qxpo(i,j,km,qrad:qradhi)) + qxpo(i,j,km,QREINT)

          qypo(i,j,km,QRHO) = rrnewry
          qypo(i,j,km,QU) = runewry/qypo(i,j,km,QRHO)
          qypo(i,j,km,QV) = rvnewry/qypo(i,j,km,QRHO)
          qypo(i,j,km,QW) = rwnewry/qypo(i,j,km,QRHO)
          rhoekenry = 0.5d0*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,km,QRHO)
          qypo(i,j,km,QREINT)= renewry - rhoekenry
          qypo(i,j,km,QPRES) = pnewry
          qypo(i,j,km,qrad:qradhi) = ernewry(:)
          qypo(i,j,km,qptot  ) = sum(lambda(:)*ernewry(:)) + qypo(i,j,km,QPRES)
          qypo(i,j,km,qreitot) = sum(qypo(i,j,km,qrad:qradhi)) + qypo(i,j,km,QREINT)

          qxmo(i+1,j,km,QRHO) = rrnewlx
          qxmo(i+1,j,km,QU) = runewlx/qxmo(i+1,j,km,QRHO)
          qxmo(i+1,j,km,QV) = rvnewlx/qxmo(i+1,j,km,QRHO)
          qxmo(i+1,j,km,QW) = rwnewlx/qxmo(i+1,j,km,QRHO)
          rhoekenlx = 0.5d0*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,km,QRHO)
          qxmo(i+1,j,km,QREINT)= renewlx - rhoekenlx
          qxmo(i+1,j,km,QPRES) = pnewlx
          qxmo(i+1,j,km,qrad:qradhi) = ernewlx(:)
          qxmo(i+1,j,km,qptot  ) = sum(lambda(:)*ernewlx(:)) + qxmo(i+1,j,km,QPRES)
          qxmo(i+1,j,km,qreitot) = sum(qxmo(i+1,j,km,qrad:qradhi)) + qxmo(i+1,j,km,QREINT)

          qymo(i,j+1,km,QRHO) = rrnewly
          qymo(i,j+1,km,QU) = runewly/qymo(i,j+1,km,QRHO)
          qymo(i,j+1,km,QV) = rvnewly/qymo(i,j+1,km,QRHO)
          qymo(i,j+1,km,QW) = rwnewly/qymo(i,j+1,km,QRHO)
          rhoekenly = 0.5d0*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,km,QRHO)
          qymo(i,j+1,km,QREINT)= renewly - rhoekenly
          qymo(i,j+1,km,QPRES) = pnewly
          qymo(i,j+1,km,qrad:qradhi) = ernewly(:)
          qymo(i,j+1,km,qptot  ) = sum(lambda(:)*ernewly(:)) + qymo(i,j+1,km,QPRES)
          qymo(i,j+1,km,qreitot) = sum(qymo(i,j+1,km,qrad:qradhi)) + qymo(i,j+1,km,QREINT)

       enddo
    enddo

  end subroutine transz_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transyz_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
       qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
       fyz, rfyz, &
       fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
       fzy, rfzy, &
       fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
       ugdnvy, pgdnvy, ergdnvy, &
       pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
       ugdnvz , pgdnvz, ergdnvz, &
       pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
       gamc,gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3, &
       srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
       hdt,cdtdy,cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
         QPRES, QREINT, QFA, QFS, QFX, &
         URHO, UMX, UMY, UMZ, UEDEN, UFA, UFS, UFX, &
         nadv
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
         fspace_type, comoving
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor

    implicit none

    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
    integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
    integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
    integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
    integer gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer ilo,ihi,jlo,jhi,km,kc,k3d

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)
    double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  fyz(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
    double precision rfyz(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,0:ngroups-1)
    double precision  fzy(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
    double precision rfzy(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,0:ngroups-1)
    double precision  ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision  pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision ergdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3,0:ngroups-1)
    double precision  ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision  pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision ergdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3,0:ngroups-1)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2,gc_l3:gc_h3)
    double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
    double precision hdt,cdtdy,cdtdz

    ! Local variables

    integer i, j, g
    integer n, nq
    integer iadv, ispec, iaux

    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision ugyp, ugym, duyp, pyav, duy, pynew
    double precision ugzp, ugzm, duzp, pzav, duz, pznew
    double precision compr, compl, compnr, compnl
    double precision rhotmp

    double precision :: dmy, dmz, dre, pggzp, pggzm, pggyp, pggym 
    double precision, dimension(0:ngroups-1) :: der, lambda, lugey, lugez, lgey, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergyp, ergzm, ergym
    double precision eddf, f1

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
                  - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))
             rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                  - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))

             compnr = compr - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                  - cdtdz*(fzy(i,j  ,kc,n) - fzy(i,j,km,n))
             compnl = compl - cdtdy*(fyz(i,j+1,km,n) - fyz(i,j,km,n)) &
                  - cdtdz*(fzy(i,j  ,kc,n) - fzy(i,j,km,n))

             qpo(i,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i+1,j,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

          enddo
       enddo
    enddo

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

    do j = jlo, jhi 
       do i = ilo, ihi 

          lambda = lam(i,j,k3d,:)

          pggyp  =  pgdnvy(i,j+1,km)
          pggym  =  pgdnvy(i,j  ,km)
          ugyp  =  ugdnvy(i,j+1,km)
          ugym  =  ugdnvy(i,j  ,km)
          ergyp = ergdnvy(i,j+1,km,:)
          ergym = ergdnvy(i,j  ,km,:)

          pggzp  =  pgdnvz(i,j,kc)
          pggzm  =  pgdnvz(i,j,km)
          ugzp  =  ugdnvz(i,j,kc)
          ugzm  =  ugdnvz(i,j,km)
          ergzp = ergdnvz(i,j,kc,:)
          ergzm = ergdnvz(i,j,km,:)

          ! Convert to conservation form
          rrr = qp(i,j,km,QRHO)
          rur = rrr*qp(i,j,km,QU)
          rvr = rrr*qp(i,j,km,QV)
          rwr = rrr*qp(i,j,km,QW)
          ekenr = 0.5d0*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + &
               qp(i,j,km,QW)**2)
          rer = qp(i,j,km,QREINT) + ekenr
          err = qp(i,j,km,qrad:qradhi)

          rrl = qm(i+1,j,km,QRHO)
          rul = rrl*qm(i+1,j,km,QU)
          rvl = rrl*qm(i+1,j,km,QV)
          rwl = rrl*qm(i+1,j,km,QW)
          ekenl = 0.5d0*rrl*(qm(i+1,j,km,QU)**2 + qm(i+1,j,km,QV)**2 + &
               qm(i+1,j,km,QW)**2)
          rel = qm(i+1,j,km,QREINT) + ekenl
          erl = qm(i+1,j,km,qrad:qradhi)

          ! Add transverse predictor
          rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO ) - fyz(i,j,km,URHO)) &
               - cdtdz*(fzy(i,j  ,kc,URHO ) - fzy(i,j,km,URHO))
          runewr = rur - cdtdy*(fyz(i,j+1,km,UMX  ) - fyz(i,j,km,UMX)) &
               - cdtdz*(fzy(i,j  ,kc,UMX  ) - fzy(i,j,km,UMX))
          rvnewr = rvr - cdtdy*(fyz(i,j+1,km,UMY  ) - fyz(i,j,km,UMY)) &
               - cdtdz*(fzy(i,j  ,kc,UMY  ) - fzy(i,j,km,UMY))
          rwnewr = rwr - cdtdy*(fyz(i,j+1,km,UMZ  ) - fyz(i,j,km,UMZ)) &
               - cdtdz*(fzy(i,j  ,kc,UMZ  ) - fzy(i,j,km,UMZ))
          lgey = lambda(:) * (ergyp(:)-ergym(:))
          lgez = lambda(:) * (ergzp(:)-ergzm(:))
          dmy = - cdtdy*sum(lgey)
          dmz = - cdtdz*sum(lgez)
          rvnewr = rvnewr + dmy 
          rwnewr = rwnewr + dmz
          lugey = 0.5d0*(ugyp+ugym) * lgey(:)
          lugez = 0.5d0*(ugzp+ugzm) * lgez(:)
          dre = -cdtdy*sum(lugey) - cdtdz*sum(lugez)
          renewr = rer - cdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
               - cdtdz*(fzy(i,j  ,kc,UEDEN) - fzy(i,j,km,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdy*0.5d0*(ugyp+ugym)*(ergyp(g)-ergym(g)) &
                     +       cdtdz*0.5d0*(ugzp+ugzm)*(ergzp(g)-ergzm(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdy*0.5d0*(ergyp(g)+ergym(g))*(ugym-ugyp) &
                     +       cdtdz*0.5d0*(ergzp(g)+ergzm(g))*(ugzm-ugzp) )
             end do
          else ! mixed frame
             der(:) = cdtdy*lugey + cdtdz*lugez
          end if

          ernewr = err(:) - cdtdy*(rfyz(i,j+1,km,:) - rfyz(i,j,km,:)) &
               &          - cdtdz*(rfzy(i,j  ,kc,:) - rfzy(i,j,km,:)) &
               &          + der(:)

          rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO ) - fyz(i,j,km,URHO)) &
               - cdtdz*(fzy(i,j  ,kc,URHO ) - fzy(i,j,km,URHO))
          runewl = rul - cdtdy*(fyz(i,j+1,km,UMX  ) - fyz(i,j,km,UMX)) &
               - cdtdz*(fzy(i,j  ,kc,UMX  ) - fzy(i,j,km,UMX))
          rvnewl = rvl - cdtdy*(fyz(i,j+1,km,UMY  ) - fyz(i,j,km,UMY)) &
               - cdtdz*(fzy(i,j  ,kc,UMY  ) - fzy(i,j,km,UMY))
          rwnewl = rwl - cdtdy*(fyz(i,j+1,km,UMZ  ) - fyz(i,j,km,UMZ)) &
               - cdtdz*(fzy(i,j  ,kc,UMZ  ) - fzy(i,j,km,UMZ))
          rvnewl = rvnewl + dmy
          rwnewl = rwnewl + dmz 
          renewl = rel - cdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
               - cdtdz*(fzy(i,j  ,kc,UEDEN) - fzy(i,j,km,UEDEN)) &
               + dre
          ernewl = erl(:) - cdtdy*(rfyz(i,j+1,km,:) - rfyz(i,j,km,:)) &
               &          - cdtdz*(rfzy(i,j  ,kc,:) - rfzy(i,j,km,:)) &
               &          + der(:)

          duyp = pggyp*ugyp - pggym*ugym
          pyav = 0.5d0*(pggyp+pggym)
          duy = ugyp-ugym
          pynew = cdtdy*(duyp + pyav*duy*(gamc(i,j,k3d)-1.d0))

          duzp = pggzp*ugzp - pggzm*ugzm
          pzav = 0.5d0*(pggzp+pggzm)
          duz = ugzp-ugzm
          pznew = cdtdz*(duzp + pzav*duz*(gamc(i,j,k3d)-1.d0))

          pnewr = qp(i,j,km,QPRES) - pynew - pznew
          pnewl = qm(i+1,j,km,QPRES) - pynew - pznew

          ! Convert back to non-conservation form
          rhotmp = rrnewr
          qpo(i,j,km,QRHO  ) = rhotmp        + hdt*srcQ(i,j,k3d,QRHO)
          qpo(i,j,km,QU    ) = runewr/rhotmp + hdt*srcQ(i,j,k3d,QU)
          qpo(i,j,km,QV    ) = rvnewr/rhotmp + hdt*srcQ(i,j,k3d,QV)
          qpo(i,j,km,QW    ) = rwnewr/rhotmp + hdt*srcQ(i,j,k3d,QW)
          rhoekenr = 0.5d0*(runewr**2 + rvnewr**2 + rwnewr**2)/rhotmp
          qpo(i,j,km,QREINT)= renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)
          qpo(i,j,km,QPRES ) = pnewr            + hdt*srcQ(i,j,k3d,QPRES)
          qpo(i,j,km,qrad:qradhi) = ernewr(:)
          qpo(i,j,km,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,km,QPRES)
          qpo(i,j,km,qreitot) = sum(qpo(i,j,km,qrad:qradhi)) + qpo(i,j,km,QREINT)

          rhotmp = rrnewl
          qmo(i+1,j,km,QRHO  ) = rhotmp        + hdt*srcQ(i,j,k3d,QRHO)
          qmo(i+1,j,km,QU    ) = runewl/rhotmp + hdt*srcQ(i,j,k3d,QU)
          qmo(i+1,j,km,QV    ) = rvnewl/rhotmp + hdt*srcQ(i,j,k3d,QV)
          qmo(i+1,j,km,QW    ) = rwnewl/rhotmp + hdt*srcQ(i,j,k3d,QW)
          rhoekenl = 0.5d0*(runewl**2 + rvnewl**2 + rwnewl**2)/rhotmp
          qmo(i+1,j,km,QREINT)= renewl - rhoekenl + hdt*srcQ(i,j,k3d,QREINT)
          qmo(i+1,j,km,QPRES ) = pnewl            + hdt*srcQ(i,j,k3d,QPRES)
          qmo(i+1,j,km,qrad:qradhi) = ernewl(:)
          qmo(i+1,j,km,qptot) = sum(lambda(:)*ernewl(:)) + qmo(i+1,j,km,QPRES)
          qmo(i+1,j,km,qreitot) = sum(qmo(i+1,j,km,qrad:qradhi)) + qmo(i+1,j,km,QREINT)
       enddo
    enddo

  end subroutine transyz_rad

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine transxz_rad(lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
       qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
       fxz, rfxz, &
       fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
       fzx, rfzx, &
       fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
       ugdnvx, pgdnvx, ergdnvx, &
       pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
       ugdnvz, pgdnvz, ergdnvz, & 
       pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
       gamc,gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3, &
       srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
       hdt,cdtdx,cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
         QPRES, QREINT, QFA, QFS, QFX, &
         URHO, UMX, UMY, UMZ, UEDEN, UFA, UFS, UFX, &
         nadv
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
         fspace_type, comoving
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor

    implicit none      

    integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
    integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
    integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
    integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
    integer gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer ilo,ihi,jlo,jhi,km,kc,k3d

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,lam_l3:lam_h3,0:ngroups-1)
    double precision  qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QRADVAR)
    double precision  fxz(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
    double precision rfxz(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,0:ngroups-1)
    double precision  fzx(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR) 
    double precision rfzx(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,0:ngroups-1)
    double precision  ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision  pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision ergdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3,0:ngroups-1)
    double precision  ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision  pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision ergdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3,0:ngroups-1)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2,gc_l3:gc_h3)
    double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
    double precision hdt,cdtdx,cdtdz

    ! Local variables

    integer i, j, g
    integer n, nq
    integer iadv, ispec, iaux

    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision ugxp, ugxm, duxp, pxav, dux, pxnew
    double precision ugzp, ugzm, duzp, pzav, duz, pznew
    double precision compr, compl, compnr, compnl
    double precision rhotmp

    double precision :: dmx, dmz, dre, pggxp, pggxm, pggzp, pggzm
    double precision, dimension(0:ngroups-1) :: der, lambda, lugex, lugez, lgex, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergxp, ergzm,  ergxm 
    double precision eddf, f1

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
                  - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))
             rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                  - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))

             compnr = compr - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                  - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))
             compnl = compl - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                  - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))

             qpo(i,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

          enddo
       enddo
    enddo

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
                  - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))
             rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                  - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))

             compnr = compr - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                  - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))
             compnl = compl - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                  - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))

             qpo(i,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

          enddo
       enddo
    enddo

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
                  - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))
             rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                  - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))

             compnr = compr - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                  - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))
             compnl = compl - cdtdx*(fxz(i+1,j,km,n) - fxz(i,j,km,n)) &
                  - cdtdz*(fzx(i  ,j,kc,n) - fzx(i,j,km,n))

             qpo(i,j  ,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)

          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          lambda = lam(i,j,k3d,:)

          pggxp  =  pgdnvx(i+1,j,km)
          pggxm  =  pgdnvx(i  ,j,km)
          ugxp  =  ugdnvx(i+1,j,km)
          ugxm  =  ugdnvx(i  ,j,km)
          ergxp = ergdnvx(i+1,j,km,:)
          ergxm = ergdnvx(i  ,j,km,:)

          pggzp  =  pgdnvz(i,j,kc)
          pggzm  =  pgdnvz(i,j,km)
          ugzp  =  ugdnvz(i,j,kc)
          ugzm  =  ugdnvz(i,j,km)
          ergzp = ergdnvz(i,j,kc,:)
          ergzm = ergdnvz(i,j,km,:)

          ! Convert to conservation form
          rrr = qp(i,j,km,QRHO)
          rur = rrr*qp(i,j,km,QU)
          rvr = rrr*qp(i,j,km,QV)
          rwr = rrr*qp(i,j,km,QW)
          ekenr = 0.5d0*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + qp(i,j,km,QW)**2)
          rer = qp(i,j,km,QREINT) + ekenr
          err = qp(i,j,km,qrad:qradhi)

          rrl = qm(i,j+1,km,QRHO)
          rul = rrl*qm(i,j+1,km,QU)
          rvl = rrl*qm(i,j+1,km,QV)
          rwl = rrl*qm(i,j+1,km,QW)
          ekenl = 0.5d0*rrl*(qm(i,j+1,km,QU)**2 + qm(i,j+1,km,QV)**2 + qm(i,j+1,km,QW)**2)
          rel = qm(i,j+1,km,QREINT) + ekenl
          erl = qm(i,j+1,km,qrad:qradhi)

          ! Add transverse predictor
          rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO ) - fxz(i,j,km,URHO)) &
               - cdtdz*(fzx(i  ,j,kc,URHO ) - fzx(i,j,km,URHO))
          runewr = rur - cdtdx*(fxz(i+1,j,km,UMX  ) - fxz(i,j,km,UMX)) &
               - cdtdz*(fzx(i  ,j,kc,UMX  ) - fzx(i,j,km,UMX))
          rvnewr = rvr - cdtdx*(fxz(i+1,j,km,UMY  ) - fxz(i,j,km,UMY)) &
               - cdtdz*(fzx(i  ,j,kc,UMY  ) - fzx(i,j,km,UMY))
          rwnewr = rwr - cdtdx*(fxz(i+1,j,km,UMZ  ) - fxz(i,j,km,UMZ)) &
               - cdtdz*(fzx(i  ,j,kc,UMZ  ) - fzx(i,j,km,UMZ))
          lgex = lambda(:) * (ergxp(:)-ergxm(:))
          lgez = lambda(:) * (ergzp(:)-ergzm(:))
          dmx = - cdtdx*sum(lgex)
          dmz = - cdtdz*sum(lgez)
          runewr = runewr + dmx
          rwnewr = rwnewr + dmz
          lugex = 0.5d0*(ugxp+ugxm) * lgex(:)
          lugez = 0.5d0*(ugzp+ugzm) * lgez(:)
          dre = -cdtdx * sum(lugex) - cdtdz * sum(lugez)
          renewr = rer - cdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
               - cdtdz*(fzx(i  ,j,kc,UEDEN) - fzx(i,j,km,UEDEN)) &
               + dre

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdx*0.5d0*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                     +       cdtdz*0.5d0*(ugzp+ugzm)*(ergzp(g)-ergzm(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = f1*(cdtdx*0.5d0*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                     +       cdtdz*0.5d0*(ergzp(g)+ergzm(g))*(ugzm-ugzp) )
             end do
          else ! mixed frame
             der(:) = cdtdx*lugex + cdtdz*lugez
          end if

          ernewr = err(:) - cdtdx*(rfxz(i+1,j,km,:) - rfxz(i,j,km,:)) &
               &          - cdtdz*(rfzx(i  ,j,kc,:) - rfzx(i,j,km,:)) &
               &          + der(:)

          rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO ) - fxz(i,j,km,URHO)) &
               - cdtdz*(fzx(i  ,j,kc,URHO ) - fzx(i,j,km,URHO))
          runewl = rul - cdtdx*(fxz(i+1,j,km,UMX  ) - fxz(i,j,km,UMX)) &
               - cdtdz*(fzx(i  ,j,kc,UMX  ) - fzx(i,j,km,UMX))
          rvnewl = rvl - cdtdx*(fxz(i+1,j,km,UMY  ) - fxz(i,j,km,UMY)) &
               - cdtdz*(fzx(i  ,j,kc,UMY  ) - fzx(i,j,km,UMY))
          rwnewl = rwl - cdtdx*(fxz(i+1,j,km,UMZ  ) - fxz(i,j,km,UMZ)) &
               - cdtdz*(fzx(i  ,j,kc,UMZ  ) - fzx(i,j,km,UMZ))
          runewl = runewl + dmx
          rwnewl = rwnewl + dmz
          renewl = rel - cdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
               - cdtdz*(fzx(i  ,j,kc,UEDEN) - fzx(i,j,km,UEDEN)) &
               + dre
          ernewl = erl(:) - cdtdx*(rfxz(i+1,j,km,:) - rfxz(i,j,km,:)) &
               &          - cdtdz*(rfzx(i  ,j,kc,:) - rfzx(i,j,km,:)) &
               &          + der(:)

          duxp = pggxp*ugxp - pggxm*ugxm
          pxav = 0.5d0*(pggxp+pggxm)
          dux = ugxp-ugxm
          pxnew = cdtdx*(duxp + pxav*dux*(gamc(i,j,k3d)-1.d0))

          duzp = pggzp*ugzp - pggzm*ugzm
          pzav = 0.5d0*(pggzp+pggzm)
          duz = ugzp-ugzm
          pznew = cdtdz*(duzp + pzav*duz*(gamc(i,j,k3d)-1.d0))

          pnewr = qp(i,j,km,QPRES) - pxnew - pznew
          pnewl = qm(i,j+1,km,QPRES) - pxnew - pznew

          ! Convert back to non-conservation form
          rhotmp = rrnewr
          qpo(i,j,km,QRHO  ) = rhotmp        + hdt*srcQ(i,j,k3d,QRHO)
          qpo(i,j,km,QU    ) = runewr/rhotmp + hdt*srcQ(i,j,k3d,QU)
          qpo(i,j,km,QV    ) = rvnewr/rhotmp + hdt*srcQ(i,j,k3d,QV)
          qpo(i,j,km,QW    ) = rwnewr/rhotmp + hdt*srcQ(i,j,k3d,QW)
          rhoekenr = 0.5d0*(runewr**2 + rvnewr**2 + rwnewr**2)/rhotmp
          qpo(i,j,km,QREINT)= renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)
          qpo(i,j,km,QPRES ) = pnewr            + hdt*srcQ(i,j,k3d,QPRES)
          qpo(i,j,km,qrad:qradhi) = ernewr(:)
          qpo(i,j,km,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,km,QPRES)
          qpo(i,j,km,qreitot) = sum(qpo(i,j,km,qrad:qradhi)) + qpo(i,j,km,QREINT)

          rhotmp = rrnewl
          qmo(i,j+1,km,QRHO  ) = rhotmp        + hdt*srcQ(i,j,k3d,QRHO)
          qmo(i,j+1,km,QU    ) = runewl/rhotmp + hdt*srcQ(i,j,k3d,QU)
          qmo(i,j+1,km,QV    ) = rvnewl/rhotmp + hdt*srcQ(i,j,k3d,QV)
          qmo(i,j+1,km,QW    ) = rwnewl/rhotmp + hdt*srcQ(i,j,k3d,QW)
          rhoekenl = 0.5d0*(runewl**2 + rvnewl**2 + rwnewl**2)/rhotmp
          qmo(i,j+1,km,QREINT)= renewl - rhoekenl + hdt*srcQ(i,j,k3d,QREINT)
          qmo(i,j+1,km,QPRES ) = pnewl            + hdt*srcQ(i,j,k3d,QPRES)
          qmo(i,j+1,km,qrad:qradhi) = ernewl(:)
          qmo(i,j+1,km,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j+1,km,QPRES)
          qmo(i,j+1,km,qreitot) = sum(qmo(i,j+1,km,qrad:qradhi)) + qmo(i,j+1,km,QREINT)

       enddo
    enddo

  end subroutine transxz_rad

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
