
! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine ca_advance_mhd(time,lo,hi,&
           uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
           uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
	   bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
	   byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
	   bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
	   bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
	   byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
	   bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
           ugdnvx,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
           ugdnvy,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
           ugdnvz,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
           src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
           grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
           delta,dt, &
           flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
           flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
           flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
           Ex,ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
           Ey,ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
           Ez,ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
           courno,e_added,ke_added,print_fortran_warnings) &
           bind(C, name="ca_advance_mhd")

!--------------------- Dependencies ------------------------------------------------
      use amrex_fort_module, only : rt => amrex_real
      use mempool_module, only : bl_allocate, bl_deallocate
      use ct_upwind, only : corner_transport, checkisnan
      use mhd_plm_module, only : plm
      use meth_params_module!, only : QVAR, NTHERM, NHYP, normalize_species, NVAR, URHO, UEDEN
      !use enforce_module, only : enforce_nonnegative_species
      use bl_constants_module

      implicit none

!-------------------- Variables -----------------------------------------------------

      integer lo(3),hi(3),print_fortran_warnings
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
      integer bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3
      integer byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3
      integer bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3
      integer bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3
      integer byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3
      integer bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3
      integer ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3
      integer ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3
      integer ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
      integer ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3
      integer ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3
      integer ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3

      real(rt)  uin(uin_l1:uin_h1, uin_l2:uin_h2, uin_l3:uin_h3,  NVAR)
      real(rt)  uout(uout_l1:uout_h1, uout_l2:uout_h2, uout_l3:uout_h3, NVAR)
      real(rt)  bxin(bxin_l1:bxin_h1, bxin_l2:bxin_h2, bxin_l3:bxin_h3)
      real(rt)  bxout(bxout_l1:bxout_h1, bxout_l2:bxout_h2, bxout_l3:bxout_h3)
      real(rt)  byin(byin_l1:byin_h1, byin_l2:byin_h2, byin_l3:byin_h3)
      real(rt)  byout(byout_l1:byout_h1, byout_l2:byout_h2, byout_l3:byout_h3)
      real(rt)  bzin(bzin_l1:bzin_h1, bzin_l2:bzin_h2, bzin_l3:bzin_h3)
      real(rt)  bzout(bzout_l1:bzout_h1, bzout_l2:bzout_h2, bzout_l3:bzout_h3)
      real(rt)  src(src_l1:src_h1, src_l2:src_h2, src_l3:src_h3, NTHERM)
      real(rt) ugdnvx(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3)
      real(rt) ugdnvy(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3)
      real(rt) ugdnvz(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3)
      real(rt)  grav( gv_l1:gv_h1, gv_l2:gv_h2, gv_l3:gv_h3, 3)

      real(rt), intent(inout) ::  flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,NVAR)
      real(rt), intent(inout) ::  flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,NVAR)
      real(rt), intent(inout) ::  flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,NVAR)

      real(rt), intent(inout) ::  Ex(ex_l1:ex_h1,ex_l2:ex_h2,ex_l3:ex_h3)
      real(rt), intent(inout) ::  Ey(ey_l1:ey_h1,ey_l2:ey_h2,ey_l3:ey_h3)
      real(rt), intent(inout) ::  Ez(ez_l1:ez_h1,ez_l2:ez_h2,ez_l3:ez_h3)

      real(rt)  delta(3),dt,time,courno
      !real(rt)  a_old, a_new
      real(rt)  e_added,ke_added

      integer flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3
      integer flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3
      integer flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3

      integer extemp_l1,extemp_l2,extemp_l3,extemp_h1,extemp_h2,extemp_h3
      integer eytemp_l1,eytemp_l2,eytemp_l3,eytemp_h1,eytemp_h2,eytemp_h3
      integer eztemp_l1,eztemp_l2,eztemp_l3,eztemp_h1,eztemp_h2,eztemp_h3

      ! Automatic arrays for workspace
      real(rt), pointer :: q(:,:,:,:)
      real(rt), pointer :: cx(:,:,:)
      real(rt), pointer :: cy(:,:,:)
      real(rt), pointer :: cz(:,:,:)
      real(rt), pointer :: csml(:,:,:)
	  real(rt), pointer :: flatn(:,:,:)
 !     real(rt), pointer :: div(:,:,:)
 !     real(rt), pointer :: pdivu(:,:,:)
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
      integer 	:: i,j,k
    
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
		
	uout(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,:) = uin(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,:)
		
      call bl_allocate(     q, lo-NHYP, hi+NHYP, QVAR)
      call bl_allocate( flatn, lo-NHYP, hi+NHYP      )
      call bl_allocate(    cx, lo-NHYP, hi+NHYP      )
      call bl_allocate(    cy, lo-NHYP, hi+NHYP      )
      call bl_allocate(    cz, lo-NHYP, hi+NHYP      )
      call bl_allocate(  csml, lo-NHYP, hi+NHYP      )
      call bl_allocate(  srcQ, lo-1, hi+1, QVAR)

      flxx_l1 = lo(1)-3
      flxx_l2 = lo(2)-3
      flxx_l3 = lo(3)-3
      flxx_h1 = hi(1)+4
      flxx_h2 = hi(2)+3
      flxx_h3 = hi(3)+3

      flxy_l1 = lo(1)-3
      flxy_l2 = lo(2)-3
      flxy_l3 = lo(3)-3
      flxy_h1 = hi(1)+3
      flxy_h2 = hi(2)+4
      flxy_h3 = hi(3)+3

      flxz_l1 = lo(1)-3
      flxz_l2 = lo(2)-3
      flxz_l3 = lo(3)-3
      flxz_h1 = hi(1)+3
      flxz_h2 = hi(2)+3
      flxz_h3 = hi(3)+4

      extemp_l1 = lo(1)-3
      extemp_l2 = lo(2)-3
      extemp_l3 = lo(3)-3
      extemp_h1 = hi(1)+3
      extemp_h2 = hi(2)+4
      extemp_h3 = hi(3)+4

      eytemp_l1 = lo(1)-3
      eytemp_l2 = lo(2)-3
      eytemp_l3 = lo(3)-3
      eytemp_h1 = hi(1)+4
      eytemp_h2 = hi(2)+3
      eytemp_h3 = hi(3)+4

      eztemp_l1 = lo(1)-3
      eztemp_l2 = lo(2)-3
      eztemp_l3 = lo(3)-3
      eztemp_h1 = hi(1)+4
      eztemp_h2 = hi(2)+4
      eztemp_h3 = hi(3)+3

      allocate(flxx(flxx_l1:flxx_h1,flxx_l2:flxx_h2,flxx_l3:flxx_h3,QVAR))
      allocate(flxy(flxy_l1:flxy_h1,flxy_l2:flxy_h2,flxy_l3:flxy_h3,QVAR))
      allocate(flxz(flxz_l1:flxz_h1,flxz_l2:flxz_h2,flxz_l3:flxz_h3,QVAR))

      allocate(Extemp(extemp_l1:extemp_h1,extemp_l2:extemp_h2,extemp_l3:extemp_h3))
      allocate(Eytemp(eytemp_l1:eytemp_h1,eytemp_l2:eytemp_h2,eytemp_l3:eytemp_h3))
      allocate(Eztemp(eztemp_l1:eztemp_h1,eztemp_l2:eztemp_h2,eztemp_l3:eztemp_h3))

      allocate(  qp(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR, 3))
      allocate(  qm(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR, 3))

      q = 0.d0

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

!Step One, Calculate Primitives based on conservatives
    call ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3,&
	  	 bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
		 byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
		 bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
                 q , cx,cy,cz , csml, flatn,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
                 src,  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                 srcQ, srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                 grav,gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                 courno,dx,dy,dz,dt,ngq,ngf)

!Step Two, Interpolate Cell centered values to faces
	  call plm(lo, hi, q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3,&	
	  	   bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
		   byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
		   bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
           qp, qm, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, dx, dy, dz, dt)

!do i =1,3
!    qp(:,:,:,:,i) = q
!    qm(:,:,:,:,i) = q
!enddo    

flxx = 0.d0
flxy = 0.d0
flxz = 0.d0

!Step Three, Corner Couple and find the correct fluxes + electric fields
	  call corner_transport( q, qm, qp, q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3, &	
				flxx,flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
				flxy,flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
				flxz,flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
                                Extemp, extemp_l1,extemp_l2,extemp_l3,extemp_h1,extemp_h2,extemp_h3, &
                                Eytemp, eytemp_l1,eytemp_l2,eytemp_l3,eytemp_h1,eytemp_h2,eytemp_h3, &
                                Eztemp, eztemp_l1,eztemp_l2,eztemp_l3,eztemp_h1,eztemp_h2,eztemp_h3, &
                                dx , dy, dz, dt)
!Step Four, Conservative update
      call consup(uin,  uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                  uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                  src ,src_l1 ,src_l2 ,src_l3 ,src_h1 ,src_h2 ,src_h3, &
		  flxx,flxx_l1,flxx_l2,flxx_l3,flxx_h1,flxx_h2,flxx_h3, &
		  flxy,flxy_l1,flxy_l2,flxy_l3,flxy_h1,flxy_h2,flxy_h3, &
		  flxz,flxz_l1,flxz_l2,flxz_l3,flxz_h1,flxz_h2,flxz_h3, &
                  lo ,hi ,dx ,dy ,dz ,dt)

!Step Five Magnetic Update
     call magup(bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
		 byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
		 bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
		 bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
		 byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
		 bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
		 uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
		 src ,  src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3, &
         Extemp, extemp_l1,extemp_l2,extemp_l3,extemp_h1,extemp_h2,extemp_h3, &
         Eytemp, eytemp_l1,eytemp_l2,eytemp_l3,eytemp_h1,eytemp_h2,eytemp_h3, &
         Eztemp, eztemp_l1,eztemp_l2,eztemp_l3,eztemp_h1,eztemp_h2,eztemp_h3, &
		 lo, hi, dx, dy, dz, dt)

	  flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,URHO:UEDEN) = flxx(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,URHO:UEDEN)
	  flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,URHO:UEDEN) = flxy(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,URHO:UEDEN)
	  flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,URHO:UEDEN) = flxz(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,URHO:UEDEN)

	  Ex(ex_l1:ex_h1,ex_l2:ex_h2, ex_l3:ex_h3) = Extemp(ex_l1:ex_h1,ex_l2:ex_h2,ex_l3:ex_h3)
	  Ey(ey_l1:ey_h1,ey_l2:ey_h2, ey_l3:ey_h3) = Eytemp(ey_l1:ey_h1,ey_l2:ey_h2,ey_l3:ey_h3)
	  Ez(ez_l1:ez_h1,ez_l2:ez_h2, ez_l3:ez_h3) = Eztemp(ez_l1:ez_h1,ez_l2:ez_h2,ez_l3:ez_h3)

      ! We are done with these here so can go ahead and free up the space
      call bl_deallocate(q)
      call bl_deallocate(flatn)
      call bl_deallocate(cx)
      call bl_deallocate(cy)
      call bl_deallocate(cz)
      call bl_deallocate(csml)
!      call bl_deallocate(div)
      call bl_deallocate(srcQ)
!     call bl_deallocate(pdivu)

	  deallocate(qm)
	  deallocate(qp)

      deallocate(flxx,flxy,flxz)
      deallocate(Extemp,Eytemp,Eztemp)

      ! Enforce the density >= small_dens.  Make sure we do this immediately after consup.
      call enforce_minimum_density(uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                                        uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                                        lo,hi,print_fortran_warnings)
      
      if (do_grav .gt. 0)  then
          call add_grav_source(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                               uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                               grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                               lo,hi,dx,dy,dz,dt,e_added,ke_added)
      endif
      ! Enforce species >= 0
 !     call enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
  !                                     uout_h1,uout_h2,uout_h3,lo,hi,0)

      ! Re-normalize the species
   !   if (normalize_species .eq. 1) then
    !     call normalize_new_species(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
     !                               lo,hi)
     ! end if

end subroutine ca_advance_mhd



! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3,&
			 bx, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
			 by, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
			 bz, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
                         q,cx, cy, cz,csml,flatn,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
                         src,  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                         srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                         grav,gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                         courno,dx,dy,dz,dt,ngp,ngf)
      !
      !     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
      !     if use_flattening=1.  Declared dimensions of q,c,csml,flatn are given
      !     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
      !     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
      !     routine that computes flatn).
      !
      use amrex_fort_module, only : rt => amrex_real
      use network, only : nspec, naux
!      use eos_params_module
      use eos_module
      use eos_type_module, only : eos_t, eos_input_re 
!      use flatten_module
      use bl_constants_module
      use meth_params_module, only : NTHERM, URHO, UMX, UMY, UMZ, &
                                     UEDEN, UEINT, UFA, UFS, &
                                     QVAR, QRHO, QU, QV, QW, QC, &
                                     QREINT, QPRES, QFA, QFS, &
                                     QMAGX,  QMAGY, QMAGZ, &
                                     nadv, small_dens, small_pres, &
                                     use_flattening

      implicit none

      real(rt), parameter:: small = 1.d-8

      integer lo(3), hi(3)
      integer  uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
      integer bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3
      integer byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3
      integer bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3
      integer    q_l1,   q_l2,   q_l3,   q_h1,   q_h2,   q_h3
      integer   gv_l1,  gv_l2,  gv_l3,  gv_h1,  gv_h2,  gv_h3
      integer  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
      integer srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3

      real(rt) :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NTHERM)
      real(rt) :: bx(bxin_l1:bxin_h1, bxin_l2:bxin_h2, bxin_l3:bxin_h3)
      real(rt) :: by(byin_l1:byin_h1, byin_l2:byin_h2, byin_l3:byin_h3)
      real(rt) :: bz(bzin_l1:bzin_h1, bzin_l2:bzin_h2, bzin_l3:bzin_h3)

      real(rt) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR) !Contains Cell Centered Mag Field
      real(rt) :: cx(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: cy(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: cz(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: csml(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: flatn(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: src( src_l1: src_h1, src_l2: src_h2, src_l3: src_h3,NTHERM)
      real(rt) :: srcQ(srcq_l1:srcq_h1,srcq_l2:srcq_h2,srcq_l3:srcq_h3,QVAR)
      real(rt) :: grav( gv_l1: gv_h1, gv_l2: gv_h2, gv_l3: gv_h3,3)
      real(rt) :: dx, dy, dz, dt, courno
      real(rt) :: dpdr, dpde

      integer          :: i, j, k
      integer          :: ngp, ngf, loq(3), hiq(3)
      integer          :: n, nq
      integer          :: iadv, ispec
      real(rt) :: courx, coury, courz, courmx, courmy, courmz, cad
      real(rt) :: a_half, a_dot, rhoInv
      real(rt) :: dtdxaold, dtdyaold, dtdzaold, small_pres_over_dens
      
      type(eos_t) :: eos_state

      do i=1,3
         loq(i) = lo(i)-ngp
         hiq(i) = hi(i)+ngp
      enddo
      !
      ! Make q (all but p), except put e in slot for rho.e, fix after eos call.
      ! The temperature is used as an initial guess for the eos call and will be overwritten.
      !
	  !Calculate Cell Centered Magnetic Field x

      do k = loq(3),hiq(3)
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
		q(i,j,k,QMAGX) = 0.5d0*(bx(i+1,j,k) + bx(i,j,k))
	   end do
	 end do
      end do

      do k = loq(3),hiq(3)
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
		 q(i,j,k,QMAGY) = 0.5d0*(by(i,j+1,k) + by(i,j,k))
            end do
	end do
      end do

      do k = loq(3),hiq(3)
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
		 q(i,j,k,QMAGZ) = 0.5d0*(bz(i,j,k+1) + bz(i,j,k))
            end do
	end do
      end do
      do k = loq(3),hiq(3)
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
               if (uin(i,j,k,URHO) .le. ZERO) then
                  !
                  ! A critical region since we usually can't write from threads.
                  !
                  print *,'   '
                  print *,'>>> Error: Nyx_advection_3d::ctoprim ',i,j,k
                  print *,'>>> ... negative density ',uin(i,j,k,URHO)
                  call bl_error("Error:: Nyx_advection_3d.f90 :: ctoprim")
               end if

               rhoInv = ONE/uin(i,j,k,URHO)

               q(i,j,k,QRHO) = uin(i,j,k,URHO)
               q(i,j,k,QU)   = uin(i,j,k,UMX)*rhoInv
               q(i,j,k,QV)   = uin(i,j,k,UMY)*rhoInv
               q(i,j,k,QW)   = uin(i,j,k,UMZ)*rhoInv

               ! Convert "rho e" to "e"
               q(i,j,k,QREINT ) = uin(i,j,k,UEINT)*rhoInv
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
      enddo

      ! Load chemical species and aux. variables, c, into q, assuming they arrived in uin as rho.c
!      if (UFS .gt. 0) then
!         do ispec = 1, nspec+naux
!            n  = UFS + ispec - 1
!            nq = QFS + ispec - 1
!            do k = loq(3),hiq(3)
!               do j = loq(2),hiq(2)
!                  do i = loq(1),hiq(1)
!                     q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
!                  enddo
!               enddo
!            enddo
!         enddo
!      end if ! UFS > 0

      small_pres_over_dens = small_pres / small_dens


      ! Get p, T, c, csml using q state
      do k = loq(3), hiq(3)
         do j = loq(2), hiq(2)
            do i = loq(1), hiq(1)

               ! If necessary, reset the energy using small_temp
               if (q(i,j,k,QREINT) .lt. ZERO) then

!                 HACK HACK HACK 
!                 call nyx_eos_given_RT(q(i,j,k,QREINT),q(i,j,k,QPRES),q(i,j,k,QRHO), &
!                                       small_temp,diag_eos(i,j,k,NE_COMP),a_old)

                  if (q(i,j,k,QREINT) .lt. ZERO) then
                     !
                     ! A critical region since we usually can't write from threads.
                     !
                     print *,'   '
                     print *,'>>> Error: Nyx_advection_3d::ctoprim ',i,j,k
                     print *,'>>> ... new e from eos_given_RT call is negative ',q(i,j,k,QREINT)
                     print *,'    '
                     call bl_error("Error:: Nyx_advection_3d.f90 :: ctoprim")
                  end if
               end if

               ! Define the magneto-accoustic speed from the EOS
			   cad = q(i,j,k,QMAGX)!(q(i,j,k,QMAGX)**2)/q(i,j,k,QRHO)
               call eos_soundspeed_mhd(cx(i,j,k), q(i,j,k,QRHO), q(i,j,k,QREINT), &
					   q(i,j,k,QMAGX), q(i,j,k,QMAGY), q(i,j,k,QMAGZ), cad, &
                                           q(i,j,k,QFS:QFS+nspec-1))

			   cad = q(i,j,k,QMAGY)!(q(i,j,k,QMAGY)**2)/q(i,j,k,QRHO)
               call eos_soundspeed_mhd(cy(i,j,k), q(i,j,k,QRHO), q(i,j,k,QREINT), &
					   q(i,j,k,QMAGX), q(i,j,k,QMAGY), q(i,j,k,QMAGZ), cad, &
                                           q(i,j,k,QFS:QFS+nspec-1))

			   cad = q(i,j,k,QMAGZ)!(q(i,j,k,QMAGZ)**2)/q(i,j,k,QRHO)
               call eos_soundspeed_mhd(cz(i,j,k), q(i,j,k,QRHO), q(i,j,k,QREINT), &
					   q(i,j,k,QMAGX), q(i,j,k,QMAGY), q(i,j,k,QMAGZ), cad, &
                                           q(i,j,k,QFS:QFS+nspec-1))


               ! Convert "e" back to "rho e"
               q(i,j,k,QREINT) = q(i,j,k,QREINT)*q(i,j,k,QRHO)

               ! Pressure = (gamma - 1) * rho * e + 0.5 B dot B
               eos_state % rho = q(i, j, k,QRHO)
               eos_state % e = q(i,j,k,QREINT) / eos_state % rho
               eos_state % xn = q(i,j,k,QFS:QFS+nspec-1)

               call eos(eos_input_re, eos_state)

               ! Set csmal based on small_pres and small_dens 
               ! TODO: this is a small sound speed -- we should do this how we do in hydro
               csml(i,j,k) = max(small, small * eos_state % cs) !? 

               q(i,j,k,QPRES) = eos_state % p &
				+ 0.5d0*(q(i,j,k,QMAGX)**2 + q(i,j,k,QMAGY)**2 + q(i,j,k,QMAGZ)**2)
            end do
         end do
      end do

      !a_half = HALF * (a_old + a_new)
      !a_dot   = (a_new - a_old) / dt

      ! Make sure these are initialized to zero.
      srcQ = ZERO

      ! NOTE - WE ASSUME HERE THAT src(i,j,k,URHO) = 0. --
      !        IF NOT THEN THE FORMULAE BELOW ARE INCOMPLETE.

      ! compute srcQ terms
  !    do k = lo(3)-1, hi(3)+1
  !       do j = lo(2)-1, hi(2)+1
  !          do i = lo(1)-1, hi(1)+1

  !             rhoInv = ONE/q(i,j,k,QRHO)

  !            srcQ(i,j,k,QRHO  ) = src(i,j,k,URHO)
  !             srcQ(i,j,k,QU    ) = src(i,j,k,UMX) * rhoInv - a_dot * q(i,j,k,QU) + &
  !                                  grav(i,j,k,1)
  !             srcQ(i,j,k,QV    ) = src(i,j,k,UMY) * rhoInv - a_dot * q(i,j,k,QV) + &
  !                                  grav(i,j,k,2)
  !             srcQ(i,j,k,QW    ) = src(i,j,k,UMZ) * rhoInv - a_dot * q(i,j,k,QW) + &
  !                                  grav(i,j,k,3)
  !             srcQ(i,j,k,QREINT) = src(i,j,k,UEDEN) - q(i,j,k,QU)*src(i,j,k,UMX) - &
  !                                                     q(i,j,k,QV)*src(i,j,k,UMY) - &
  !                                                     q(i,j,k,QW)*src(i,j,k,UMZ) - &
  !                                                     a_dot * THREE * gamma_minus_1 * q(i,j,k,QREINT)

   !            dpde = gamma_minus_1 * q(i,j,k,QRHO)
   !            dpdr = gamma_minus_1 * q(i,j,k,QREINT)/q(i,j,k,QRHO)
   !            srcQ(i,j,k,QPRES ) = dpde * srcQ(i,j,k,QREINT) * rhoInv &
   !                               + dpdr * srcQ(i,j,k,QRHO)

   !            if (UFS .gt. 0) then
   !               do ispec = 1,nspec+naux
   !                  srcQ(i,j,k,QFS+ispec-1) = src(i,j,k,UFS+ispec-1)*rhoInv
   !               enddo
   !            end if ! UFS > 0

   !            do iadv = 1,nadv
   !               srcQ(i,j,k,QFA+iadv-1) = src(i,j,k,UFA+iadv-1)*rhoInv
   !            enddo

   !         enddo
   !      enddo
   !   enddo

      ! Compute running max of Courant number over grids
      courmx = courno
      courmy = courno
      courmz = courno

      dtdxaold = dt / dx
      dtdyaold = dt / dy
      dtdzaold = dt / dz

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               courx = ( cx(i,j,k)+abs(q(i,j,k,QU)) ) * dtdxaold
               coury = ( cy(i,j,k)+abs(q(i,j,k,QV)) ) * dtdyaold
               courz = ( cz(i,j,k)+abs(q(i,j,k,QW)) ) * dtdzaold

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
                  call bl_error("Error:: Nyx_advection_3d.f90 :: CFL violation in x-dir in ctoprim")
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
                  call bl_error("Error:: Nyx_advection_3d.f90 :: CFL violation in y-dir in ctoprim")
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
                  call bl_error("Error:: Nyx_advection_3d.f90 :: CFL violation in z-dir in ctoprim")
               end if

            enddo
         enddo
      enddo
      courno = max( courmx, courmy, courmz )
      end subroutine ctoprim
! :::
! ::: ========================== Conservative Update ===============================================================
! ::: 

	subroutine consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                          uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                          src ,  src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3, &
		          fluxx,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                          fluxy,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                          fluxz,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                          lo,hi,dx,dy,dz,dt)

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module, only : QVAR, UMX,UMY,UMZ, NVAR, URHO, UEDEN, UEINT

	implicit none

 	  integer,  intent(in)  :: uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
	  integer,  intent(in)  :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
      integer,  intent(in)  :: flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer,  intent(in)  :: flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer,  intent(in)  :: flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
	  integer,  intent(in)  :: src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3
	  integer, intent(in) 	:: lo(3), hi(3)

	  real(rt), intent(in)  :: uin(uin_l1:uin_h1, uin_l2:uin_h2, uin_l3:uin_h3, NVAR)
	  real(rt), intent(in)  :: src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3, NVAR)
	  real(rt), intent(in)  :: fluxx(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,QVAR)
	  real(rt), intent(in)  :: fluxy(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,QVAR)
	  real(rt), intent(in)  :: fluxz(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,QVAR)
	  real(rt), intent(in) 	:: dx,dy,dz,dt 
	  real(rt), intent(out) :: uout(uout_l1:uout_h1,uout_l2:uout_h2, uout_l3:uout_h3,NVAR)
	  real(rt)				:: u, v, w

	  integer 				:: i, j, k	
	  !****TO DO ******* SOURCES
		do k = lo(3), hi(3)
		do j = lo(2), hi(2)
		do i = lo(1), hi(1)
		   uout(i,j,k,URHO:UEDEN) = uin(i,j,k,URHO:UEDEN) - dt/dx*(fluxx(i+1,j,k,URHO:UEDEN) - fluxx(i,j,k,URHO:UEDEN)) &
		 						  - dt/dy*(fluxy(i,j+1,k,URHO:UEDEN) - fluxy(i,j,k,URHO:UEDEN)) &
		 						  - dt/dz*(fluxz(i,j,k+1,URHO:UEDEN) - fluxz(i,j,k,URHO:UEDEN)) !Add source terms later
		   u = uout(i,j,k,UMX)/uout(i,j,k,URHO)
		   v = uout(i,j,k,UMY)/uout(i,j,k,URHO)
   		   w = uout(i,j,k,UMZ)/uout(i,j,k,URHO)
		   uout(i,j,k,UEINT) = uout(i,j,k,UEDEN) - 0.5d0*uout(i,j,k,URHO)*(u**2 + v**2 + w**2)
     	   !if(uout(i,j,k,UEDEN).le.0.77d0) then
			!print*, uin(i,j,k,UEDEN), uout(i,j,k,UEDEN), "i j k = ", i, j, k
			!print*, "flux x = ", fluxx(i+1,j,k,UEDEN) , fluxx(i,j,k,UEDEN)
			!print*, "flux y = ", fluxy(i,j+1,k,UEDEN) , fluxy(i,j,k,UEDEN)
			!print*, "flux z = ", fluxz(i,j,k+1,UEDEN) , fluxz(i,j,k,UEDEN)
			!print*, " E = ", uout(i,j,k,UEDEN)
			!pause
!		   endif
		enddo
		enddo
		enddo
	end subroutine consup

! :::
! ::: ========================== Magnetic Update ===============================================================
! ::: 

	subroutine magup(bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
		         byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
        		 bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
        		 bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
        		 byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
        		 bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
        		 uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
        		 src ,  src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3, &
                         Ex,ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3, &
                         Ey,ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3, &
                         Ez,ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3, &
        		 lo, hi, dx, dy, dz, dt)

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module!, only : QVAR, NVAR, UEINT

	implicit none
	
	integer, intent(in)   :: bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3
	integer, intent(in)   :: byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3
	integer, intent(in)   :: bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3
	integer, intent(in)   :: bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3
	integer, intent(in)   :: byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3
	integer, intent(in)   :: bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3
	integer, intent(in)	  :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
	integer, intent(in)   :: src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3
    integer, intent(in)   ::  ex_l1,ex_l2,ex_l3,ex_h1,ex_h2,ex_h3
    integer, intent(in)   ::  ey_l1,ey_l2,ey_l3,ey_h1,ey_h2,ey_h3
    integer, intent(in)   ::  ez_l1,ez_l2,ez_l3,ez_h1,ez_h2,ez_h3
	integer, intent(in)   :: lo(3), hi(3)

	real(rt), intent(in)  :: bxin(bxin_l1:bxin_h1, bxin_l2:bxin_h2, bxin_l3:bxin_h3)
	real(rt), intent(in)  :: byin(byin_l1:byin_h1, byin_l2:byin_h2, byin_l3:byin_h3)
	real(rt), intent(in)  :: bzin(bzin_l1:bzin_h1, bzin_l2:bzin_h2, bzin_l3:bzin_h3)
	real(rt), intent(in)  :: src(src_l1:src_h1, src_l2:src_h2, src_l3:src_h3, QVAR)

    real(rt), intent(in) ::  Ex(ex_l1:ex_h1,ex_l2:ex_h2, ex_l3:ex_h3)
    real(rt), intent(in) ::  Ey(ey_l1:ey_h1,ey_l2:ey_h2, ey_l3:ey_h3)
    real(rt), intent(in) ::  Ez(ez_l1:ez_h1,ez_l2:ez_h2, ez_l3:ez_h3)

	real(rt), intent(in)  :: dx, dy, dz, dt

	real(rt), intent(inout) :: uout(uout_l1:uout_h1,uout_l2:uout_h2, uout_l3:uout_h3,NVAR)

	real(rt), intent(out) :: bxout(bxout_l1:bxout_h1, bxout_l2:bxout_h2, bxout_l3:bxout_h3)
	real(rt), intent(out) :: byout(byout_l1:byout_h1, byout_l2:byout_h2, byout_l3:byout_h3)
	real(rt), intent(out) :: bzout(bzout_l1:bzout_h1, bzout_l2:bzout_h2, bzout_l3:bzout_h3)

	real(rt)			  :: bx, by ,bz, e
	integer				  :: i, j, k
		
	!***** TO DO ***** SOURCES
	!-------------------------------- bx --------------------------------------------------
	do k = lo(3), hi(3)
	do j = lo(2), hi(2)
	do i = lo(1), hi(1)+1
		bxout(i,j,k) = bxin(i,j,k) + dt/dz*(Ey(i,j,k+1) - Ey(i,j,k)) - dt/dy*(Ez(i,j+1,k) - Ez(i,j,k))
	enddo
	enddo
	enddo

	!------------------------------- by --------------------------------------------------
	do k = lo(3), hi(3)
	do j = lo(2), hi(2)+1
	do i = lo(1), hi(1)
		byout(i,j,k) = byin(i,j,k) + dt/dx*(Ez(i+1,j,k) - Ez(i,j,k)) - dt/dz*(Ex(i,j,k+1) - Ex(i,j,k))
		!if(i.eq.3.and.j.eq.64.and.k.eq.2) then
		!	print *, "byout = ", byout(i,j,k), "at ", i, j ,k
		!	print *, "byin = ", byin(i,j,k)
		!	print *, "Ez =", Ez(i+1, j, k), Ez(i, j, k)
		!	print *, "Ex =", Ex(i, j, k+1), Ex(i, j, k)
		!	pause
		!endif
	enddo
	enddo
	enddo
!			print *, "byout = ", byout(2,17,1), "at ", 2, 17 ,1
!			print *, "byin = ", byin(2,17,1)
!			print *, "Ez =", Ez(3,17,1), Ez(2,17,1)
!			print *, "Ex =", Ex(2,17,2), Ex(2,17,1)
!			print *, "byout = ", byout(3,17,1), "at ", 3,17,1
!			print *, "byin = ", byin(3,17,1)
!			print *, "Ez =", Ez(4,17,1), Ez(3,17,1)
!			print *, "Ex =", Ex(3,17,2), Ex(3,17,1)
!			pause
	!------------------------------- bz --------------------------------------------------
	do k = lo(3), hi(3)+1
	do j = lo(2), hi(2)
	do i = lo(1), hi(1)
		bzout(i,j,k) = bzin(i,j,k) + dt/dy*(Ex(i,j+1,k) - Ex(i,j,k)) - dt/dx*(Ey(i+1,j,k) - Ey(i,j,k))
	enddo
	enddo
	enddo
	!-------------------------------- Internal Energy ----------------------------------------------------
	do k = lo(3), hi(3)
	do j = lo(2), hi(2)
	do i = lo(1), hi(1)
		bx = 0.5d0*(bxout(i+1,j,k)+bxout(i,j,k))
		by = 0.5d0*(byout(i,j+1,k)+byout(i,j,k))
		bz = 0.5d0*(bzout(i,j,k+1)+bzout(i,j,k))
		e = uout(i,j,k,UEINT)
		uout(i,j,k,UEINT) = e - 0.5d0*(bx**2 + by**2 + bz**2)
		if(uout(i,j,k,UEINT).le.0.d0) then
			print*, "e < 0 !", uout(i,j,k,UEINT)
			print*, "e before = ", e
			print*, "Total NRG = ", uout(i,j,k, UEDEN)
			print*, "bx = ", bx, "by = ", by, "bz = ", bz
			print*, "byout =", byout(i,j,k), byout(i,j+1,k)
			print*, "-1/2|B| = ", - 0.5d0*(bx**2 + by**2 + bz**2)
			print*, "i j k = ", i, j, k
			pause
		endif
	enddo
	enddo
	enddo			

	end subroutine magup
