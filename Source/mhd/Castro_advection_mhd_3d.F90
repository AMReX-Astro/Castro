module mhd_advection_module

  implicit none

contains

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine ca_advance_mhd(time, lo, hi, &
                          uin, uin_lo, uin_hi, &
                          uout, uout_lo, uout_hi, &
                          bxin, bxin_lo, bxin_hi, &
                          byin, byin_lo, byin_hi, &
                          bzin, bzin_lo, bzin_hi, &
                          bxout, bxout_lo, bxout_hi, &
                          byout, byout_lo, byout_hi, &
                          bzout, bzout_lo, bzout_hi, &
                          src, src_lo, src_hi, &
                          delta, dt, &
                          flux1, flux1_lo, flux1_hi, &
                          flux2, flux2_lo, flux2_hi, &
                          flux3, flux3_lo, flux3_hi, &
                          Ex, ex_lo, ex_hi, &
                          Ey, ey_lo, ey_hi, &
                          Ez, ez_lo, ez_hi, &
                          courno, e_added, ke_added, print_fortran_warnings) &
                          bind(C, name="ca_advance_mhd")

  !--------------------- Dependencies ------------------------------------------------
  use amrex_fort_module, only : rt => amrex_real
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use ct_upwind, only : corner_transport
  use mhd_plm_module, only : plm
  use flatten_module_mhd, only : mhd_flatten
  use meth_params_module!, only : QVAR, NTHERM, NHYP, normalize_species, NVAR, URHO, UEDEN
  use prob_params_module, only : dg
  use amrex_constants_module
  use network, only: nspec

  implicit none

  !-------------------- Variables -----------------------------------------------------

  integer , intent(in   ) :: lo(3),hi(3), print_fortran_warnings
  integer , intent(in   ) :: uin_lo(3), uin_hi(3)
  integer , intent(in   ) :: uout_lo(3), uout_hi(3)
  integer , intent(in   ) :: bxin_lo(3), bxin_hi(3)
  integer , intent(in   ) :: byin_lo(3), byin_hi(3)
  integer , intent(in   ) :: bzin_lo(3), bzin_hi(3)
  integer , intent(in   ) :: bxout_lo(3), bxout_hi(3)
  integer , intent(in   ) :: byout_lo(3), byout_hi(3)
  integer , intent(in   ) :: bzout_lo(3), bzout_hi(3)
  integer , intent(in   ) :: src_lo(3), src_hi(3)
  integer , intent(in   ) :: flux1_lo(3), flux1_hi(3)
  integer , intent(in   ) :: flux2_lo(3), flux2_hi(3)
  integer , intent(in   ) :: flux3_lo(3), flux3_hi(3)
  integer , intent(in   ) :: ex_lo(3), ex_hi(3)
  integer , intent(in   ) :: ey_lo(3), ey_hi(3)
  integer , intent(in   ) :: ez_lo(3), ez_hi(3)

  real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
  real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)

  real(rt), intent(in   ) :: bxin(bxin_lo(1):bxin_hi(1), bxin_lo(2):bxin_hi(2), bxin_lo(3):bxin_hi(3))
  real(rt), intent(in   ) :: byin(byin_lo(1):byin_hi(1), byin_lo(2):byin_hi(2), byin_lo(3):byin_hi(3))
  real(rt), intent(in   ) :: bzin(bzin_lo(1):bzin_hi(1), bzin_lo(2):bzin_hi(2), bzin_lo(3):bzin_hi(3))

  real(rt), intent(inout) :: bxout(bxout_lo(1):bxout_hi(1), bxout_lo(2):bxout_hi(2), bxout_lo(3):bxout_hi(3))
  real(rt), intent(inout) :: byout(byout_lo(1):byout_hi(1), byout_lo(2):byout_hi(2), byout_lo(3):byout_hi(3))
  real(rt), intent(inout) :: bzout(bzout_lo(1):bzout_hi(1), bzout_lo(2):bzout_hi(2), bzout_lo(3):bzout_hi(3))

  real(rt), intent(in   ) :: src(src_lo(1):src_hi(1), src_lo(2):src_hi(2), src_lo(3):src_hi(3), NSRC)

  real(rt), intent(inout) :: flux1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3), NVAR)
  real(rt), intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
  real(rt), intent(inout) :: flux3(flux3_lo(1):flux3_hi(1), flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3), NVAR)

  real(rt), intent(inout) :: Ex(ex_lo(1):ex_hi(1), ex_lo(2):ex_hi(2), ex_lo(3):ex_hi(3))
  real(rt), intent(inout) :: Ey(ey_lo(1):ey_hi(1), ey_lo(2):ey_hi(2), ey_lo(3):ey_hi(3))
  real(rt), intent(inout) :: Ez(ez_lo(1):ez_hi(1), ez_lo(2):ez_hi(2), ez_lo(3):ez_hi(3))

  real(rt), intent(in   ) :: delta(3), time, dt, courno
  real(rt), intent(in   ) :: e_added, ke_added

  integer flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3
  integer flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3
  integer flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3

  integer extemp_l1,extemp_l2,extemp_l3,extemp_h1,extemp_h2,extemp_h3
  integer eytemp_l1,eytemp_l2,eytemp_l3,eytemp_h1,eytemp_h2,eytemp_h3
  integer eztemp_l1,eztemp_l2,eztemp_l3,eztemp_h1,eztemp_h2,eztemp_h3

  ! Automatic arrays for workspace
  real(rt), pointer :: q(:,:,:,:)
  real(rt), pointer :: bcc(:,:,:,:)
  real(rt), pointer :: cx(:,:,:)
  real(rt), pointer :: cy(:,:,:)
  real(rt), pointer :: cz(:,:,:)
  real(rt), pointer :: csml(:,:,:)
  real(rt), pointer :: flatn(:,:,:)
  real(rt), pointer :: srcQ(:,:,:,:)

  real(rt), allocatable :: flxx(:,:,:,:)
  real(rt), allocatable :: flxy(:,:,:,:)
  real(rt), allocatable :: flxz(:,:,:,:)

  real(rt), allocatable :: Extemp(:,:,:)
  real(rt), allocatable :: Eytemp(:,:,:)
  real(rt), allocatable :: Eztemp(:,:,:)

  real(rt), allocatable :: qp(:,:,:,:,:)
  real(rt), allocatable :: qm(:,:,:,:,:)

  real(rt) dx,dy,dz
  integer ngq,ngf
  integer q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  integer srcq_l1, srcq_l2, srcq_l3, srcq_h1, srcq_h2, srcq_h3
  integer   :: i,j,k

  ngq = NHYP
  ngf = 1
  q_l1 = lo(1)-NHYP
  q_l2 = lo(2)-NHYP
  q_l3 = lo(3)-NHYP
  q_h1 = hi(1)+NHYP
  q_h2 = hi(2)+NHYP
  q_h3 = hi(3)+NHYP

  srcq_l1 = lo(1)-1
  srcq_l2 = lo(2)-1
  srcq_l3 = lo(3)-1
  srcq_h1 = hi(1)+1
  srcq_h2 = hi(2)+1
  srcq_h3 = hi(3)+1


  call bl_allocate(     q, lo-NHYP, hi+NHYP, NQ)
  call bl_allocate(   bcc, lo-NHYP, hi+NHYP,  3  )
  call bl_allocate( flatn, lo-NHYP, hi+NHYP      )
  call bl_allocate(    cx, lo-NHYP, hi+NHYP      )
  call bl_allocate(    cy, lo-NHYP, hi+NHYP      )
  call bl_allocate(    cz, lo-NHYP, hi+NHYP      )
  call bl_allocate(  csml, lo-NHYP, hi+NHYP      )
  call bl_allocate(  srcQ, lo-1, hi+1, NQSRC)

  flxx_l1 = lo(1)-2
  flxx_l2 = lo(2)-3
  flxx_l3 = lo(3)-3
  flxx_h1 = hi(1)+3
  flxx_h2 = hi(2)+3
  flxx_h3 = hi(3)+3

  flxy_l1 = lo(1)-3
  flxy_l2 = lo(2)-2
  flxy_l3 = lo(3)-3
  flxy_h1 = hi(1)+3
  flxy_h2 = hi(2)+3
  flxy_h3 = hi(3)+3

  flxz_l1 = lo(1)-3
  flxz_l2 = lo(2)-3
  flxz_l3 = lo(3)-2
  flxz_h1 = hi(1)+3
  flxz_h2 = hi(2)+3
  flxz_h3 = hi(3)+3

  extemp_l1 = lo(1)-2
  extemp_l2 = lo(2)-2
  extemp_l3 = lo(3)-2
  extemp_h1 = hi(1)+2
  extemp_h2 = hi(2)+3
  extemp_h3 = hi(3)+3

  eytemp_l1 = lo(1)-2
  eytemp_l2 = lo(2)-2
  eytemp_l3 = lo(3)-2
  eytemp_h1 = hi(1)+3
  eytemp_h2 = hi(2)+2
  eytemp_h3 = hi(3)+3

  eztemp_l1 = lo(1)-2
  eztemp_l2 = lo(2)-2
  eztemp_l3 = lo(3)-2
  eztemp_h1 = hi(1)+3
  eztemp_h2 = hi(2)+3
  eztemp_h3 = hi(3)+2

  allocate(flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,NVAR+3))
  allocate(flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,NVAR+3))
  allocate(flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,NVAR+3))

  allocate(Extemp(extemp_l1:extemp_h1,extemp_l2:extemp_h2,extemp_l3:extemp_h3))
  allocate(Eytemp(eytemp_l1:eytemp_h1,eytemp_l2:eytemp_h2,eytemp_l3:eytemp_h3))
  allocate(Eztemp(eztemp_l1:eztemp_h1,eztemp_l2:eztemp_h2,eztemp_l3:eztemp_h3))

  allocate(  qp(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ, 3))
  allocate(  qm(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,NQ, 3))

  q = 0.d0

  dx = delta(1)
  dy = delta(2)
  dz = delta(3)

  !Step One, Calculate Primitives based on conservatives
  call ctoprim(lo,hi,uin,uin_lo,uin_hi, &
               bcc, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
               bxin, bxin_lo, bxin_hi, &
               byin, byin_lo, byin_hi, &
               bzin, bzin_lo, bzin_hi, &
               q , cx,cy,cz , csml, flatn,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
               src,  src_lo, src_hi, &
               srcQ, srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
               courno,dx,dy,dz,dt,ngq,ngf)

  flatn(:,:,:) = 0.0
  call mhd_flatten(lo-dg*ngf, hi+dg*ngf, &
                q, flatn,lo-NHYP, hi+NHYP)

  !Step Two, Interpolate Cell centered values to faces
  call plm(lo, hi, q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3,&
           flatn, &
           bxin, bxin_lo, bxin_hi, &
           byin, byin_lo, byin_hi, &
           bzin, bzin_lo, bzin_hi, &
           qp, qm, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
           srcQ, srcq_l1, srcq_l2, srcq_l3, srcq_h1, srcq_h2, srcq_h3, &
           dx, dy, dz, dt)

  flxx = 0.d0
  flxy = 0.d0
  flxz = 0.d0

  !Step Three, Corner Couple and find the correct fluxes + electric fields
  call corner_transport(q, qm, qp, q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3, &
                        flxx,flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
                        flxy,flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
                        flxz,flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
                        Extemp, extemp_l1,extemp_l2,extemp_l3,extemp_h1,extemp_h2,extemp_h3, &
                        Eytemp, eytemp_l1,eytemp_l2,eytemp_l3,eytemp_h1,eytemp_h2,eytemp_h3, &
                        Eztemp, eztemp_l1,eztemp_l2,eztemp_l3,eztemp_h1,eztemp_h2,eztemp_h3, &
                        dx , dy, dz, dt)

  !Step Four, Conservative update
  call consup(uin,  uin_lo, uin_hi, &
              uout, uout_lo, uout_hi, &
              bcc, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
              flxx,flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
              flxy,flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
              flxz,flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
              lo ,hi ,dx ,dy ,dz ,dt)

  bxout = 0.0
  byout = 0.0
  bzout = 0.0

  !Step Five Magnetic Update
  call magup(bxin, bxin_lo, bxin_hi, &
             byin, byin_lo, byin_hi, &
             bzin, bzin_lo, bzin_hi, &
             bxout, bxout_lo, bxout_hi, &
             byout, byout_lo, byout_hi, &
             bzout, bzout_lo, bzout_hi, &
             Extemp, extemp_l1,extemp_l2,extemp_l3,extemp_h1,extemp_h2,extemp_h3, &
             Eytemp, eytemp_l1,eytemp_l2,eytemp_l3,eytemp_h1,eytemp_h2,eytemp_h3, &
             Eztemp, eztemp_l1,eztemp_l2,eztemp_l3,eztemp_h1,eztemp_h2,eztemp_h3, &
             lo, hi, dx, dy, dz, dt)


  flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),URHO) = &
       flxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),URHO)
  flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),UMX) = &
       flxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),UMX)
  flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),UMY) = &
       flxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),UMY)
  flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),UMZ) = &
       flxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),UMZ)
  flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),UEDEN) = &
       flxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),UEDEN)
  flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),UFS:UFS+nspec-1) = &
       flxx(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3),UFS:UFS+nspec-1)

  flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),URHO) = &
       flxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),URHO)
  flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UMX) = &
       flxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UMX)
  flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UMY) = &
       flxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UMY)
  flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UMZ) = &
       flxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UMZ)
  flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UEDEN) = &
       flxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UEDEN)
  flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UFS:UFS+nspec-1) = &
       flxy(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3),UFS:UFS+nspec-1)


  flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),URHO) = &
       flxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),URHO)
  flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UMX) = &
       flxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UMX)
  flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UMY) = &
       flxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UMY)
  flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UMZ) = &
       flxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UMZ)
  flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UEDEN) = &
       flxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UEDEN)
  flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UFS:UFS+nspec-1) = &
       flxz(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3),UFS:UFS+nspec-1)


  Ex(ex_lo(1):ex_hi(1),ex_lo(2):ex_hi(2), ex_lo(3):ex_hi(3)) = Extemp(ex_lo(1):ex_hi(1),ex_lo(2):ex_hi(2),ex_lo(3):ex_hi(3))
  Ey(ey_lo(1):ey_hi(1),ey_lo(2):ey_hi(2), ey_lo(3):ey_hi(3)) = Eytemp(ey_lo(1):ey_hi(1),ey_lo(2):ey_hi(2),ey_lo(3):ey_hi(3))
  Ez(ez_lo(1):ez_hi(1),ez_lo(2):ez_hi(2), ez_lo(3):ez_hi(3)) = Eztemp(ez_lo(1):ez_hi(1),ez_lo(2):ez_hi(2),ez_lo(3):ez_hi(3))

  ! We are done with these here so can go ahead and free up the space
  call bl_deallocate(q)
  call bl_deallocate(bcc)
  call bl_deallocate(flatn)
  call bl_deallocate(cx)
  call bl_deallocate(cy)
  call bl_deallocate(cz)
  call bl_deallocate(csml)
  call bl_deallocate(srcQ)

  deallocate(qm)
  deallocate(qp)

  deallocate(flxx,flxy,flxz)
  deallocate(Extemp,Eytemp,Eztemp)



end subroutine ca_advance_mhd



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
                                 QREINT, QPRES, QFA, QFS, QTEMP, QDPDE, QDPDR,&
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
  integer, intent(in) :: src_lo(3), src_hi(3)
  integer, intent(in) :: srcq_lo(3), srcq_hi(3)

  real(rt), intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
  real(rt), intent(in) :: bx(bxin_lo(1):bxin_hi(1), bxin_lo(2):bxin_hi(2), bxin_lo(3):bxin_hi(3))
  real(rt), intent(in) :: by(byin_lo(1):byin_hi(1), byin_lo(2):byin_hi(2), byin_lo(3):byin_hi(3))
  real(rt), intent(in) :: bz(bzin_lo(1):bzin_hi(1), bzin_lo(2):bzin_hi(2), bzin_lo(3):bzin_hi(3))

  real(rt), intent(inout) :: bcc(bcc_lo(1):bcc_hi(1), bcc_lo(2):bcc_hi(2), bcc_lo(3):bcc_hi(3), 3)
  real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
  real(rt), intent(inout) :: cx(cx_lo(1):cx_hi(1), cx_lo(2):cx_hi(2), cx_lo(3):cx_hi(3))
  real(rt), intent(inout) :: cy(cy_lo(1):cy_hi(1), cy_lo(2):cy_hi(2), cy_lo(3):cx_hi(3))
  real(rt), intent(inout) :: cz(cz_lo(1):cz_hi(1), cz_lo(2):cz_hi(2), cz_lo(3):cx_hi(3))
  real(rt), intent(in) :: src(src_lo(1):src_hi(1), src_lo(2):src_hi(2), src_lo(3):src_hi(3), NSRC)
  real(rt), intent(inout) :: srcQ(srcq_lo(1):srcq_hi(1), srcq_lo(2):srcq_hi(2), srcq_lo(3):srcq_hi(3), NQSRC)

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
           end if

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
  use eos_type_module, only : eos_t, eos_input_re
  !      use flatten_module
  use amrex_constants_module
  use amrex_error_module
  use meth_params_module, only : URHO, UMX, UMY, UMZ, NVAR,&
                                 UEDEN, UEINT, UFA, UFS, UTEMP, &
                                 NQ, NSRC ,NQSRC, QRHO, QU, QV, QW, QC, NQAUX,&
                                 QREINT, QPRES, QFA, QFS, QTEMP, QDPDE, QDPDR,&
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
           eos_state % T   = uin(i,j,k,UTEMP)
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
                                 QREINT, QPRES, QFA, QFS, QTEMP, QDPDE, QDPDR,&
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

  integer :: i, j, k

  real(rt) :: courx, coury, courz, courmx, courmy, courmz
  real(rt) :: dtdx_old, dtdy_old, dtdz_old

  ! Compute running max of Courant number over grids
  courmx = courno
  courmy = courno
  courmz = courno

  dtdx = dt / dx
  dtdy = dt / dy
  dtdz = dt / dz

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
                  dx, dt)

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
  real(rt), intent(in)  :: dx(3), dt
  real(rt), intent(out) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3),NVAR)

  integer :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           uout(i,j,k,URHO) = uout(i,j,k,URHO) - 1.0/dx*(fluxx(i+1,j,k,URHO) - fluxx(i,j,k,URHO)) &
                - 1.0/dy*(fluxy(i,j+1,k,URHO) - fluxy(i,j,k,URHO)) &
                - 1.0/dz*(fluxz(i,j,k+1,URHO) - fluxz(i,j,k,URHO))
           uout(i,j,k,UMX) = uout(i,j,k,UMX) - 1.0/dx*(fluxx(i+1,j,k,UMX) - fluxx(i,j,k,UMX)) &
                - 1.0/dy*(fluxy(i,j+1,k,UMX) - fluxy(i,j,k,UMX)) &
                - 1.0/dz*(fluxz(i,j,k+1,UMX) - fluxz(i,j,k,UMX))
           uout(i,j,k,UMY) = uout(i,j,k,UMY) - 1.0/dx*(fluxx(i+1,j,k,UMY) - fluxx(i,j,k,UMY)) &
                - 1.0/dy*(fluxy(i,j+1,k,UMY) - fluxy(i,j,k,UMY)) &
                - 1.0/dz*(fluxz(i,j,k+1,UMY) - fluxz(i,j,k,UMY))
           uout(i,j,k,UMZ) = uout(i,j,k,UMZ) - 1.0/dx*(fluxx(i+1,j,k,UMZ) - fluxx(i,j,k,UMZ)) &
                - 1.0/dy*(fluxy(i,j+1,k,UMZ) - fluxy(i,j,k,UMZ)) &
                - 1.0/dz*(fluxz(i,j,k+1,UMZ) - fluxz(i,j,k,UMZ))
           uout(i,j,k,UEDEN) = uout(i,j,k,UEDEN) - 1.0/dx*(fluxx(i+1,j,k,UEDEN) - fluxx(i,j,k,UEDEN)) &
                - 1.0/dy*(fluxy(i,j+1,k,UEDEN) - fluxy(i,j,k,UEDEN)) &
                - 1.0/dz*(fluxz(i,j,k+1,UEDEN) - fluxz(i,j,k,UEDEN))
           uout(i,j,k,UEINT) = uout(i,j,k,UEINT) - 1.0/dx*(fluxx(i+1,j,k,UEINT) - fluxx(i,j,k,UEINT)) &
                - 1.0/dy*(fluxy(i,j+1,k,UEINT) - fluxy(i,j,k,UEINT)) &
                - 1.0/dz*(fluxz(i,j,k+1,UEINT) - fluxz(i,j,k,UEINT))
           uout(i,j,k,UFS:UFS+nspec-1) = uout(i,j,k,UFS:UFS+nspec-1) - 1.0/dx*(fluxx(i+1,j,k,UFS:UFS+nspec-1) - &
           fluxx(i,j,k,UFS:UFS+nspec-1)) - 1.0/dy*(fluxy(i,j+1,k,UFS:UFS+nspec-1) - fluxy(i,j,k,UFS:UFS+nspec-1)) &
                - 1.0/dz*(fluxz(i,j,k+1,UFS:UFS+nspec-1) - fluxz(i,j,k,UFS:UFS+nspec-1))

           bcc(i,j,k,:) = bcc(i,j,k,:) - dt/dx*(fluxx(i+1, j, k, NVAR+1:NVAR+3)- fluxx(i,j,k, NVAR+1:NVAR+3)) &
                                       - dt/dy*(fluxy(i, j+1, k, NVAR+1:NVAR+3)- fluxy(i,j,k, NVAR+1:NVAR+3)) &
                                       - dt/dz*(fluxz(i, j, k+1, NVAR+1:NVAR+3)- fluxz(i,j,k, NVAR+1:NVAR+3))
        enddo
     enddo
  enddo

end subroutine consup

! :::
! ::: ========================== Magnetic Update ===============================================================
! :::

subroutine magup(bxin, bxin_lo, bxin_hi, &
                 byin, byin_lo, byin_hi, &
                 bzin, bzin_lo, bzin_hi, &
                 bxout, bxout_lo, bxout_hi, &
                 byout, byout_lo, byout_hi, &
                 bzout, bzout_lo, bzout_hi, &
                 Ex,ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                 Ey,ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                 Ez,ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
                 lo, hi, dx, dy, dz, dt)

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module!, only : QVAR, NVAR, UEINT

  implicit none

  integer, intent(in)   :: bxin_lo(3), bxin_hi(3)
  integer, intent(in)   :: byin_lo(3), byin_hi(3)
  integer, intent(in)   :: bzin_lo(3), bzin_hi(3)
  integer, intent(in)   :: bxout_lo(3), bxout_hi(3)
  integer, intent(in)   :: byout_lo(3), byout_hi(3)
  integer, intent(in)   :: bzout_lo(3), bzout_hi(3)
  integer, intent(in)   ::  ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3
  integer, intent(in)   ::  ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3
  integer, intent(in)   ::  ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3
  integer, intent(in)   :: lo(3), hi(3)

  real(rt), intent(in)  :: bxin(bxin_lo(1):bxin_hi(1), bxin_lo(2):bxin_hi(2), bxin_lo(3):bxin_hi(3))
  real(rt), intent(in)  :: byin(byin_lo(1):byin_hi(1), byin_lo(2):byin_hi(2), byin_lo(3):byin_hi(3))
  real(rt), intent(in)  :: bzin(bzin_lo(1):bzin_hi(1), bzin_lo(2):bzin_hi(2), bzin_lo(3):bzin_hi(3))

  real(rt), intent(in) ::  Ex(ex_l1:ex_h1,ex_l2:ex_h2, ex_l3:ex_h3)
  real(rt), intent(in) ::  Ey(ey_l1:ey_h1,ey_l2:ey_h2, ey_l3:ey_h3)
  real(rt), intent(in) ::  Ez(ez_l1:ez_h1,ez_l2:ez_h2, ez_l3:ez_h3)

  real(rt), intent(in)  :: dx, dy, dz, dt


  real(rt), intent(out) :: bxout(bxout_lo(1):bxout_hi(1), bxout_lo(2):bxout_hi(2), bxout_lo(3):bxout_hi(3))
  real(rt), intent(out) :: byout(byout_lo(1):byout_hi(1), byout_lo(2):byout_hi(2), byout_lo(3):byout_hi(3))
  real(rt), intent(out) :: bzout(bzout_lo(1):bzout_hi(1), bzout_lo(2):bzout_hi(2), bzout_lo(3):bzout_hi(3))

  integer                                 :: i, j, k

  !***** TO DO ***** SOURCES
  !-------------------------------- bx --------------------------------------------------
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           bxout(i,j,k) = bxin(i,j,k) + dt/dx*( (Ey(i,j,k+1) - Ey(i,j,k)) - (Ez(i,j+1,k) - Ez(i,j,k)) )

        enddo
     enddo
  enddo

  !------------------------------- by --------------------------------------------------
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
           byout(i,j,k) = byin(i,j,k) + dt/dy*( (Ez(i+1,j,k) - Ez(i,j,k)) - (Ex(i,j,k+1) - Ex(i,j,k)) )
        enddo
     enddo
  enddo
 !------------------------------- bz --------------------------------------------------
  do k = lo(3), hi(3)+1
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           bzout(i,j,k) = bzin(i,j,k) + dt/dz*( (Ex(i,j+1,k) - Ex(i,j,k)) - (Ey(i+1,j,k) - Ey(i,j,k)) )
        enddo
     enddo
  enddo

 end subroutine magup
