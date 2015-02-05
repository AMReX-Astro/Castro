subroutine ca_umdrv_rad(is_finest_level,time,&
                        lo,hi,domlo,domhi,&
                        uin,uin_l1,uin_h1,&
                        uout,uout_l1,uout_h1,&
                        Erin,Erin_l1,Erin_h1,&
                        lam,lam_l1,lam_h1,&
                        Erout,Erout_l1,Erout_h1,&
                        ugdnv,ugdnv_l1,ugdnv_h1,&
                        src,src_l1,src_h1, &
                        grav,gv_l1,gv_h1, &
                        delta,dt,&
                        flux,flux_l1,flux_h1,&
                        radflux,radflux_l1,radflux_h1,&
                        area,area_l1,area_h1,&
                        dloga,dloga_l1,dloga_h1,&
                        vol,vol_l1,vol_h1,courno,verbose, &
                        nstep_fsp)

  use meth_params_module, only : QVAR, QU, NVAR, NHYP, do_sponge, normalize_species
  use rad_params_module, only : ngroups
  use radhydro_params_module, only : QRADVAR
  use advection_module, only : enforce_minimum_density, normalize_new_species
  use rad_advection_module, only : umeth1d_rad, ctoprim_rad, consup_rad
  use sponge_module, only : sponge

  implicit none

  integer nstep_fsp
  integer is_finest_level
  integer lo(1),hi(1),verbose
  integer domlo(1),domhi(1)
  integer uin_l1,uin_h1, Erin_l1, Erin_h1, lam_l1, lam_h1
  integer uout_l1,uout_h1, Erout_l1, Erout_h1
  integer ugdnv_l1,ugdnv_h1
  integer flux_l1,flux_h1, radflux_l1,radflux_h1
  integer area_l1,area_h1
  integer dloga_l1,dloga_h1
  integer vol_l1,vol_h1
  integer src_l1,src_h1
  integer gv_l1,gv_h1
  double precision   uin(  uin_l1:  uin_h1,NVAR)
  double precision  uout( uout_l1: uout_h1,NVAR)
  double precision  Erin( Erin_l1: Erin_h1, 0:ngroups-1)
  double precision  lam( lam_l1: lam_h1, 0:ngroups-1)
  double precision Erout(Erout_l1:Erout_h1, 0:ngroups-1)
  double precision ugdnv(ugdnv_l1:ugdnv_h1)
  double precision   src(  src_l1:  src_h1,NVAR)
  double precision  grav(   gv_l1:   gv_h1     )
  double precision  flux( flux_l1: flux_h1,NVAR)
  double precision radflux(radflux_l1: radflux_h1, 0:ngroups-1)
  double precision  area( area_l1: area_h1     )
  double precision dloga(dloga_l1:dloga_h1     )
  double precision   vol(  vol_l1: vol_h1      )
  double precision delta(1),dt,time,courno
  
  !     Automatic arrays for workspace
  double precision, allocatable:: q(:,:)
  double precision, allocatable:: gamc(:)
  double precision, allocatable:: gamcg(:)
  double precision, allocatable:: flatn(:)
  double precision, allocatable:: c(:)
  double precision, allocatable:: cg(:)
  double precision, allocatable:: csml(:)
  double precision, allocatable:: div(:)
  double precision, allocatable:: pgdnv(:)
  double precision, allocatable:: ergdnv(:,:)
  double precision, allocatable:: lamgdnv(:,:)
  double precision, allocatable:: srcQ(:,:)
  double precision, allocatable:: pdivu(:)
  
  double precision dx,mass_added,eint_added,eden_added
  integer i,ngf,ngq,iflaten
  integer q_l1, q_h1

  dx = delta(1)

  ngq = NHYP
  ngf = 1
  iflaten = 1

  q_l1 = lo(1)-NHYP
  q_h1 = hi(1)+NHYP

  allocate(     q(q_l1:q_h1,QRADVAR))  ! 
  allocate(    c (q_l1:q_h1))
  allocate(    cg(q_l1:q_h1))
  allocate( gamc (q_l1:q_h1))
  allocate( gamcg(q_l1:q_h1))
  allocate( flatn(q_l1:q_h1))
  allocate(  csml(q_l1:q_h1))
  
  allocate(  srcQ(lo(1)-1:hi(1)+1,QVAR))
  
  allocate(   div(lo(1):hi(1)+1))
  allocate( pdivu(lo(1):hi(1)  ))
  allocate( pgdnv(lo(1):hi(1)+1))
  allocate(ergdnv(lo(1):hi(1)+1, 0:ngroups-1))
  allocate(lamgdnv(lo(1):hi(1)+1, 0:ngroups-1))
  

  !     Translate to primitive variables, compute sound speeds
  !     Note that (q,c,gamc,csml,flatn) are all dimensioned the same
  !       and set to correspond to coordinates of (lo:hi)
  
  call ctoprim_rad(lo,hi,uin,uin_l1,uin_h1, &
       Erin, Erin_l1, Erin_h1, &
       lam, lam_l1, lam_h1, &
       q,c,cg,gamc,gamcg,csml,flatn,q_l1,q_h1, &
       src,src_l1,src_h1, &
       srcQ,lo(1)-1,hi(1)+1, &
       courno,dx,dt,NHYP,ngf,iflaten)

  call umeth1d_rad(lo,hi,domlo,domhi, &
       lam, lam_l1, lam_h1, &       
       q,c,cg,gamc,gamcg,csml,flatn,q_l1,q_h1, &
       srcQ,lo(1)-1,hi(1)+1, &
       grav, gv_l1, gv_h1, &
       lo(1),hi(1),dx,dt, &
       flux,flux_l1,flux_h1, &
       radflux,radflux_l1,radflux_h1, &
       pgdnv,lo(1),hi(1)+1, &
       ergdnv,lo(1),hi(1)+1, &
       lamgdnv,lo(1),hi(1)+1, &
       ugdnv,ugdnv_l1,ugdnv_h1, &
       dloga,dloga_l1,dloga_h1)

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
  call consup_rad(uin,uin_l1,uin_h1, &
       uout,uout_l1,uout_h1, &
       Erin,Erin_l1,Erin_h1, &
       Erout,Erout_l1,Erout_h1, &
       pgdnv,lo(1),hi(1)+1, &
       ergdnv,lo(1),hi(1)+1, &
       lamgdnv,lo(1),hi(1)+1, &
       ugdnv,ugdnv_l1,ugdnv_h1, &
       src , src_l1, src_h1, &
       grav,  gv_l1,  gv_h1, &
       flux,flux_l1,flux_h1, &
       radflux,radflux_l1,radflux_h1, &
       flatn,uin_l1,uin_h1, &
       area,area_l1,area_h1, &
       vol , vol_l1, vol_h1, &
       div ,pdivu,lo,hi,dx,dt, &
       nstep_fsp)

  ! Enforce the density >= small_dens.
  mass_added = 0.d0
  eint_added = 0.d0
  eden_added = 0.d0
  call enforce_minimum_density(uin,uin_l1,uin_h1,uout,uout_l1,uout_h1,lo,hi,&
       mass_added,eint_added,eden_added,verbose)
  
  ! Enforce that the species >= 0
  call ca_enforce_nonnegative_species(uout,uout_l1,uout_h1,lo,hi)
  
  ! Normalize the species
  if (normalize_species .eq. 1) &
       call normalize_new_species(uout,uout_l1,uout_h1,lo,hi)
  
  if (do_sponge .eq. 1) &
       call sponge(uout,uout_l1,uout_h1,lo,hi,time,dt,dx,domlo,domhi)
  
  deallocate(q,c,cg,gamc,gamcg,flatn,csml,srcQ,div,pdivu,pgdnv,ergdnv,lamgdnv)

end subroutine ca_umdrv_rad


! This subroutine cannot be tiled
subroutine ca_compute_lamborder(Er, Er_l1, Er_h1, &
     kap, kap_l1, kap_h1, &
     lam, lam_l1, lam_h1, &
     dx, ngrow, limiter, filter_T, S)

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDlambda
  use filter_module

  implicit none

  integer, intent(in) :: Er_l1, Er_h1, kap_l1, kap_h1, lam_l1, lam_h1
  integer, intent(in) :: ngrow, limiter, filter_T, S
  double precision, intent(in) :: dx
  double precision, intent(in) :: kap(kap_l1:kap_h1, 0:ngroups-1)
  double precision, intent(in) :: Er(Er_l1:Er_h1, 0:ngroups-1)
  double precision, intent(out) :: lam(lam_l1:lam_h1, 0:ngroups-1)

  integer :: i, reg_l1, reg_h1, g
  double precision :: r

  double precision, allocatable :: lamfil(:)

  lam = -1.d50

  reg_l1 = lam_l1 + ngrow
  reg_h1 = lam_h1 - ngrow

  if (filter_T .gt. 0) then
     allocate(lamfil(reg_l1:reg_h1))
  end if

  do g = 0, ngroups-1
     do i=lam_l1, lam_h1
        r = abs(Er(i+1,g) - Er(i-1,g)) / (2.*dx)
        r = r / (kap(i,g) * max(Er(i,g), 1.d-50))
        lam(i,g) = FLDlambda(r, limiter)
     end do

     if (Er(reg_l1-1,g) .eq. -1.d0) then
        r = abs(Er(reg_l1+1,g) - Er(reg_l1,g)) / dx
        r = r / (kap(reg_l1,g) * max(Er(reg_l1,g), 1.d-50))
        lam(reg_l1,g) = FLDlambda(r, limiter)
     end if
     
     if (Er(reg_h1+1,g) .eq. -1.d0) then
        r = abs(Er(reg_h1,g) - Er(reg_h1-1,g)) / dx
        r = r / (kap(reg_h1,g) * max(Er(reg_h1,g), 1.d-50))
        lam(reg_h1,g) = FLDlambda(r, limiter)
     end if

     ! filter
     if (filter_T .eq. 1) then

        do i=reg_l1, reg_h1
           lamfil(i) = ff1(0) * lam(i,g) &
                &    + ff1(1) * (lam(i-1,g)+lam(i+1,g))
           lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))
         end do
         
         if (Er(reg_l1-1,g) .eq. -1.d0) then
            i = reg_l1
            lamfil(i) = dot_product(ff1b, lam(i:i+1,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))
         end if
         
         if (Er(reg_h1+1,g) .eq. -1.d0) then
            i = reg_h1
            lamfil(i) = dot_product(ff1b(1:0:-1), lam(i-1:i,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))
         end if
         
         lam(reg_l1:reg_h1,g) = lamfil(reg_l1:reg_h1)
         
      else if (filter_T .eq. 2) then
         
         do i=reg_l1, reg_h1
            lamfil(i) = ff2(0,S) * lam(i,g) &
                 &    + ff2(1,S) * (lam(i-1,g)+lam(i+1,g)) &
                 &    + ff2(2,S) * (lam(i-2,g)+lam(i+2,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))
         end do
         
         if (Er(reg_l1-1,g) .eq. -1.d0) then
            i = reg_l1
            lamfil(i) = dot_product(ff2b0, lam(i:i+2,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))
            
            i = reg_l1 + 1
            lamfil(i) = dot_product(ff2b1, lam(i-1:i+2,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))
         end if

         if (Er(reg_h1+1,g) .eq. -1.d0) then
            i = reg_h1-1
            lamfil(i) = dot_product(ff2b1(2:-1:-1), lam(i-2:i+1,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))

            i = reg_h1
            lamfil(i) = dot_product(ff2b0(2:0:-1), lam(i-2:i,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))
         end if

         lam(reg_l1:reg_h1,g) = lamfil(reg_l1:reg_h1)

      else if (filter_T .eq. 3) then

         do i=reg_l1, reg_h1
            lamfil(i) = ff3(0,S) * lam(i,g) &
                 &    + ff3(1,S) * (lam(i-1,g)+lam(i+1,g)) &
                 &    + ff3(2,S) * (lam(i-2,g)+lam(i+2,g)) &
                 &    + ff3(3,S) * (lam(i-3,g)+lam(i+3,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))
         end do
         
         if (Er(reg_l1-1,g) .eq. -1.d0) then
            i = reg_l1
            lamfil(i) = dot_product(ff3b0, lam(i:i+3,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            
            
            i = reg_l1 + 1
            lamfil(i) = dot_product(ff3b1, lam(i-1:i+3,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            

            i = reg_l1 + 2
            lamfil(i) = dot_product(ff3b2, lam(i-2:i+3,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            
         end if

         if (Er(reg_h1+1,g) .eq. -1.d0) then
            i = reg_h1 - 2
            lamfil(i) = dot_product(ff3b2(3:-2:-1), lam(i-3:i+2,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            

            i = reg_h1 - 1
            lamfil(i) = dot_product(ff3b1(3:-1:-1), lam(i-3:i+1,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            

            i = reg_h1
            lamfil(i) = dot_product(ff3b0(3:0:-1), lam(i-3:i,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i))) 
         end if

         lam(reg_l1:reg_h1,g) = lamfil(reg_l1:reg_h1)

      else if (filter_T .eq. 4) then

         do i=reg_l1, reg_h1
            lamfil(i) = ff4(0,S) * lam(i,g) &
                 &    + ff4(1,S) * (lam(i-1,g)+lam(i+1,g)) &
                 &    + ff4(2,S) * (lam(i-2,g)+lam(i+2,g)) &
                 &    + ff4(3,S) * (lam(i-3,g)+lam(i+3,g)) &
                 &    + ff4(4,S) * (lam(i-4,g)+lam(i+4,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))
         end do
         
         if (Er(reg_l1-1,g) .eq. -1.d0) then
            i = reg_l1 
            lamfil(i) = dot_product(ff4b0, lam(i:i+4,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            
            
            i = reg_l1 + 1
            lamfil(i) = dot_product(ff4b1, lam(i-1:i+4,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            

            i = reg_l1 + 2
            lamfil(i) = dot_product(ff4b2, lam(i-2:i+4,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            

            i = reg_l1 + 3
            lamfil(i) = dot_product(ff4b3, lam(i-3:i+4,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            
         end if

         if (Er(reg_h1+1,g) .eq. -1.d0) then
            i = reg_h1 - 3
            lamfil(i) = dot_product(ff4b3(4:-3:-1), lam(i-4:i+3,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            

            i = reg_h1 - 2
            lamfil(i) = dot_product(ff4b2(4:-2:-1), lam(i-4:i+2,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            

            i = reg_h1 - 1
            lamfil(i) = dot_product(ff4b1(4:-1:-1), lam(i-4:i+1,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            

            i = reg_h1
            lamfil(i) = dot_product(ff4b0(4:0:-1), lam(i-4:i,g))
            lamfil(i) = min(1.d0/3.d0, max(1.d-25, lamfil(i)))            
         end if

         lam(reg_l1:reg_h1,g) = lamfil(reg_l1:reg_h1)

     end if

     ! boundary

     if (Er(reg_l1-1,g) .eq. -1.d0) then
        do i=lam_l1, reg_l1-1
           lam(i,g) = lam(reg_l1,g)
        end do
     end if
     
     if (Er(reg_h1+1,g) .eq. -1.d0) then
        do i=reg_h1+1, lam_h1
           lam(i,g) = lam(reg_h1,g)
        end do
     end if

  end do

  if (filter_T .gt. 0) then
     deallocate(lamfil)
  end if

end subroutine ca_compute_lamborder


subroutine ca_get_v_dcf( lo, hi, &
     er ,  er_l1,  er_h1, &
     s  ,   s_l1,   s_h1, &
     T  ,   T_l1,   T_h1, &
     c_v, c_v_l1, c_v_h1, &
     kr ,  kr_l1,  kr_h1, &
     kp ,  kp_l1,  kp_h1, &
     kp2, kp2_l1, kp2_h1, &
     dtemp, dtime, sigma, c, &
     v  ,   v_l1,   v_h1, &
     dcf, dcf_l1, dcf_h1)

  use meth_params_module, only : NVAR, URHO, UMX 

  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: er_l1,er_h1,s_l1,s_h1,T_l1,T_h1,c_v_l1,c_v_h1
  integer, intent(in) :: kr_l1,kr_h1,kp_l1,kp_h1,kp2_l1,kp2_h1
  integer, intent(in) :: v_l1,v_h1,dcf_l1,dcf_h1 
  double precision, intent(in)  ::  er( er_l1: er_h1)
  double precision, intent(in)  ::   s(  s_l1:  s_h1, NVAR)
  double precision, intent(in)  ::   T(  T_l1:  T_h1)
  double precision, intent(in)  :: c_v(c_v_l1:c_v_h1)
  double precision, intent(in ) ::  kr( kr_l1: kr_h1)
  double precision, intent(in ) ::  kp( kp_l1: kp_h1)
  double precision, intent(in ) :: kp2(kp2_l1:kp2_h1)
  double precision, intent(in) :: dtemp, dtime, sigma, c
  double precision, intent(out) ::   v(  v_l1:  v_h1)
  double precision, intent(out) :: dcf(dcf_l1:dcf_h1)

  integer :: i
  double precision :: etainv, fac0, fac2, alpha, frc

  fac0 = 4.d0 * sigma * dtime / dtemp
  fac2 = c * dtime / dtemp

  do i=lo(1),hi(1)
     v(i) = s(i,UMX)/s(i,URHO)

     alpha = fac0 * (kp2(i) * (T(i) + dtemp) ** 4    &
          -          kp (i) * (T(i)        ) ** 4)   &
          -  fac2 * (kp2(i) - kp(i)) * er(i)

     frc = s(i,URHO) * c_v(i) + 1.0d-50
     etainv = frc / (alpha + frc)

     dcf(i) = 2.d0 * etainv * (kp(i) / kr(i))
  end do

end subroutine ca_get_v_dcf


subroutine ca_compute_dcoefs( lo, hi, &
     d  ,   d_l1,   d_h1, &
     lam, lam_l1, lam_h1, &
     v ,    v_l1,   v_h1, &
     dcf, dcf_l1, dcf_h1, &
     r, idir)

  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: d_l1, d_h1, &
       & lam_l1, lam_h1, &
       &   v_l1,   v_h1, &
       & dcf_l1, dcf_h1, &
       idir

  double precision, intent(out) ::   d(  d_l1:  d_h1)
  double precision, intent(in)  :: lam(lam_l1:lam_h1)
  double precision, intent(in)  ::   v(  v_l1:  v_h1)
  double precision, intent(in)  :: dcf(dcf_l1:dcf_h1)
  double precision, intent(in)  ::   r( lo(1): hi(1))

  integer :: i

  do i = lo(1), hi(1)
     if (v(i-1) + v(i) .gt. 0.d0) then
        d(i) = dcf(i-1) * v(i-1) * lam(i)
     else if (v(i-1) + v(i) .lt. 0.d0) then
        d(i) = dcf(i) * v(i) * lam(i)
     else
        d(i) = 0.0
     end if
     d(i) = d(i) * r(i)
  end do

end subroutine ca_compute_dcoefs


subroutine ca_update_dcf(lo, hi, &
     dcf, dcf_l1, dcf_h1, &
     etainv, eti_l1, eti_h1, &
     kp, kp_l1, kp_h1, kr, kr_l1, kr_h1)

  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: dcf_l1, dcf_h1, eti_l1, eti_h1, &
       kp_l1, kp_h1, kr_l1, kr_h1
  double precision, intent(in) :: etainv(eti_l1:eti_h1)
  double precision, intent(in) :: kp(kp_l1:kp_h1)
  double precision, intent(in) :: kr(kr_l1:kr_h1)
  double precision, intent(out) :: dcf(dcf_l1:dcf_h1)

  integer :: i

  do i=lo(1),hi(1)
     dcf(i) = 2.d0 * etainv(i) * (kp(i)/kr(i))
  end do

end subroutine ca_update_dcf

subroutine ca_set_dterm_face( lo, hi, &
     Er, Er_l1, Er_h1, dc, dc_l1, dc_h1, &
     dtf, dtf_l1, dtf_h1, dx, idir)
  implicit none
  
  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: Er_l1, Er_h1, dc_l1, dc_h1, dtf_l1, dtf_h1, idir
  double precision, intent(in) :: dx
  double precision, intent(in) :: Er(Er_l1:Er_h1)
  double precision, intent(in) :: dc(dc_l1:dc_h1)
  double precision, intent(out) :: dtf(dtf_l1:dtf_h1)
  integer :: i

  do i=lo(1),hi(1)
     dtf(i) = (Er(i) - Er(i-1)) / dx * dc(i)
  end do

end subroutine ca_set_dterm_face

subroutine ca_face2center( lo, hi, &
  foox, foox_l1, foox_h1, &
  fooc, fooc_l1, fooc_h1)

  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: foox_l1, foox_h1
  integer, intent(in) :: fooc_l1, fooc_h1
  double precision, intent(in)  :: foox(foox_l1:foox_h1)
  double precision, intent(out) :: fooc(fooc_l1:fooc_h1)

  integer :: i

  do i=lo(1), hi(1)
     fooc(i) = (foox(i) + foox(i+1))/2.d0
  end do

end subroutine ca_face2center

! no tiling
subroutine ca_correct_dterm(dfx, dfx_l1, dfx_h1, &
     re, rc)

  implicit none

  integer, intent(in) :: dfx_l1, dfx_h1
  double precision, intent(inout) :: dfx(dfx_l1:dfx_h1)
  double precision, intent(in) :: re(dfx_l1:dfx_h1), rc(1)

  integer :: i

  do i=dfx_l1, dfx_h1
     dfx(i) = dfx(i) / (re(i) + 1.d-50)
  end do

end subroutine ca_correct_dterm

subroutine ca_estdt_rad(u,u_l1,u_h1, gpr,gpr_l1,gpr_h1, &
  lo,hi,dx,dt)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UEINT, UTEMP, UFS, UFX, &
       allow_negative_energy
  implicit none
  
  integer u_l1,u_h1
  integer gpr_l1,gpr_h1
  integer lo(1), hi(1)
  double precision u(u_l1:u_h1,NVAR)
  double precision gpr(gpr_l1:gpr_h1)
  double precision dx(1), dt

  double precision :: rhoInv,ux,dt1,c
  integer          :: i
  type(eos_t) :: eos_state

  !     Translate to primitive variables, compute sound speed (call eos), get dtmax
  do i = lo(1),hi(1)

     rhoInv = 1.d0 / u(i,URHO)

     eos_state % rho = u(i,URHO)
     eos_state % T   = u(i,UTEMP)
     eos_state % e   = u(i,UEINT)*rhoInv
     eos_state % xn  = u(i,UFS:UFS+nspec-1) * rhoInv
     eos_state % aux = u(i,UFX:UFX+naux-1) * rhoInv

     ! Protect against negative e
     if (eos_state % e .gt. 0.d0 .or. allow_negative_energy .eq. 1) then
        call eos(eos_input_re, eos_state)
        c = eos_state % cs
     else
        c = 0.d0
     end if

     c = sqrt(c**2 + gpr(i)*rhoInv)

     ux = u(i,UMX)*rhoInv
     
     dt1 = dx(1) /( c + abs(ux) )
     dt  = min(dt,dt1)

  enddo
  
end subroutine ca_estdt_rad


! this is tiling safe
subroutine ca_est_gpr0(Er, Er_l1, Er_h1, gPr, gPr_l1, gPr_h1)

  use rad_params_module, only : ngroups
  
  implicit none
  
  integer, intent(in) :: Er_l1, Er_h1, gpr_l1, gpr_h1
  double precision, intent(in) :: Er(Er_l1:Er_h1, 0:ngroups-1)
  double precision, intent(out) :: gPr(gPr_l1:gPr_h1)

  integer :: i, g

  gPr = 0.d0

  do g = 0, ngroups-1
     do i = gPr_l1, gPr_h1
        gPr(i) = gPr(i) + 4.d0/9.d0*Er(i,g)
     end do
  end do

end subroutine ca_est_gpr0


! this is tiling safe
subroutine ca_est_gpr2(kap, kap_l1, kap_h1, Er, Er_l1, Er_h1, &
     gPr, gPr_l1, gPr_h1, vlo, vhi, dx, limiter, comoving)

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDlambda, Edd_factor

  implicit none

  integer, intent(in) :: kap_l1, kap_h1 
  integer, intent(in) ::  Er_l1,  Er_h1
  integer, intent(in) :: gPr_l1, gPr_h1
  integer, intent(in) :: vlo(1), vhi(1)  ! region with valid Er
  integer, intent(in) :: limiter, comoving
  double precision, intent(in) :: dx
  double precision, intent(in) ::  kap(kap_l1:kap_h1, 0:ngroups-1) 
  double precision, intent(in) ::   Er( Er_l1: Er_h1, 0:ngroups-1)
  double precision, intent(out) :: gPr(gPr_l1:gPr_h1)

  integer :: i, g
  double precision :: r, lam, f, gamr
  integer :: im, ip
  double precision :: xm, xp

  if (gPr_l1-1 .ge. vlo(1)) then
     im = 1
     xm = 2.d0
  else
     im = 0
     xm = 1.d0
  end if

  if (gPr_h1+1 .le. vhi(1)) then
     ip = 1
     xp = 2.d0
  else
     ip = 0
     xp = 1.d0
  end if

  gPr = 0.0d0

  do g = 0, ngroups-1
     i = gPr_l1
     r = abs(Er(i+1,g) - Er(i-im,g)) / (xm*dx)
     r = r / (kap(i,g) * max(Er(i,g), 1.d-50))
     lam = FLDlambda(r, limiter)
     if (comoving .eq. 1) then
        f = Edd_factor(lam)
        gamr = (3.d0-f)/2.d0
     else
        gamr = lam + 1.d0
     end if
     gPr(i) = gPr(i) + gamr * lam * Er(i,g)
     
     do i = gPr_l1+1, gPr_h1-1
        r = abs(Er(i+1,g) - Er(i-1,g)) / (2.d0*dx)
        r = r / (kap(i,g) * max(Er(i,g), 1.d-50))
        lam = FLDlambda(r, limiter)
        if (comoving .eq. 1) then
           f = Edd_factor(lam)
           gamr = (3.d0-f)/2.d0
        else
           gamr = lam + 1.d0
        end if
        gPr(i) = gPr(i) + gamr * lam * Er(i,g)
     end do

     i = gPr_h1
     r = abs(Er(i+ip,g) - Er(i-1,g)) / (xp*dx)
     r = r / (kap(i,g) * max(Er(i,g), 1.d-50))
     lam = FLDlambda(r, limiter)
     if (comoving .eq. 1) then
        f = Edd_factor(lam)
        gamr = (3.d0-f)/2.d0
     else
        gamr = lam + 1.d0
     end if
     gPr(i) = gPr(i) + gamr * lam * Er(i,g)
  end do

end subroutine ca_est_gpr2

