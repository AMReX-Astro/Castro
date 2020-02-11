! This subroutine cannot be tiled
subroutine ca_compute_lamborder(Er, Er_l1, Er_h1, &
     kap, kap_l1, kap_h1, &
     lam, lam_l1, lam_h1, &
     dx, ngrow, limiter, filter_T, S)

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDlambda
  use filter_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: Er_l1, Er_h1, kap_l1, kap_h1, lam_l1, lam_h1
  integer, intent(in) :: ngrow, limiter, filter_T, S
  real(rt)        , intent(in) :: dx
  real(rt)        , intent(in) :: kap(kap_l1:kap_h1, 0:ngroups-1)
  real(rt)        , intent(in) :: Er(Er_l1:Er_h1, 0:ngroups-1)
  real(rt)        , intent(out) :: lam(lam_l1:lam_h1, 0:ngroups-1)

  integer :: i, reg_l1, reg_h1, g
  real(rt)         :: r

  real(rt)        , allocatable :: lamfil(:)

  lam = -1.e50_rt

  reg_l1 = lam_l1 + ngrow
  reg_h1 = lam_h1 - ngrow

  if (filter_T .gt. 0) then
     allocate(lamfil(reg_l1:reg_h1))
  end if

  do g = 0, ngroups-1
     do i=lam_l1, lam_h1
        r = abs(Er(i+1,g) - Er(i-1,g)) / (2.*dx)
        r = r / (kap(i,g) * max(Er(i,g), 1.e-50_rt))
        lam(i,g) = FLDlambda(r, limiter)
     end do

     if (Er(reg_l1-1,g) .eq. -1.e0_rt) then
        r = abs(Er(reg_l1+1,g) - Er(reg_l1,g)) / dx
        r = r / (kap(reg_l1,g) * max(Er(reg_l1,g), 1.e-50_rt))
        lam(reg_l1,g) = FLDlambda(r, limiter)
     end if

     if (Er(reg_h1+1,g) .eq. -1.e0_rt) then
        r = abs(Er(reg_h1,g) - Er(reg_h1-1,g)) / dx
        r = r / (kap(reg_h1,g) * max(Er(reg_h1,g), 1.e-50_rt))
        lam(reg_h1,g) = FLDlambda(r, limiter)
     end if

     ! filter
     if (filter_T .eq. 1) then

        do i=reg_l1, reg_h1
           lamfil(i) = ff1(0) * lam(i,g) &
                &    + ff1(1) * (lam(i-1,g)+lam(i+1,g))
           lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))
         end do

         if (Er(reg_l1-1,g) .eq. -1.e0_rt) then
            i = reg_l1
            lamfil(i) = dot_product(ff1b, lam(i:i+1,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))
         end if

         if (Er(reg_h1+1,g) .eq. -1.e0_rt) then
            i = reg_h1
            lamfil(i) = dot_product(ff1b(1:0:-1), lam(i-1:i,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))
         end if

         lam(reg_l1:reg_h1,g) = lamfil(reg_l1:reg_h1)

      else if (filter_T .eq. 2) then

         do i=reg_l1, reg_h1
            lamfil(i) = ff2(0,S) * lam(i,g) &
                 &    + ff2(1,S) * (lam(i-1,g)+lam(i+1,g)) &
                 &    + ff2(2,S) * (lam(i-2,g)+lam(i+2,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))
         end do

         if (Er(reg_l1-1,g) .eq. -1.e0_rt) then
            i = reg_l1
            lamfil(i) = dot_product(ff2b0, lam(i:i+2,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))

            i = reg_l1 + 1
            lamfil(i) = dot_product(ff2b1, lam(i-1:i+2,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))
         end if

         if (Er(reg_h1+1,g) .eq. -1.e0_rt) then
            i = reg_h1-1
            lamfil(i) = dot_product(ff2b1(2:-1:-1), lam(i-2:i+1,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))

            i = reg_h1
            lamfil(i) = dot_product(ff2b0(2:0:-1), lam(i-2:i,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))
         end if

         lam(reg_l1:reg_h1,g) = lamfil(reg_l1:reg_h1)

      else if (filter_T .eq. 3) then

         do i=reg_l1, reg_h1
            lamfil(i) = ff3(0,S) * lam(i,g) &
                 &    + ff3(1,S) * (lam(i-1,g)+lam(i+1,g)) &
                 &    + ff3(2,S) * (lam(i-2,g)+lam(i+2,g)) &
                 &    + ff3(3,S) * (lam(i-3,g)+lam(i+3,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))
         end do

         if (Er(reg_l1-1,g) .eq. -1.e0_rt) then
            i = reg_l1
            lamfil(i) = dot_product(ff3b0, lam(i:i+3,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))

            i = reg_l1 + 1
            lamfil(i) = dot_product(ff3b1, lam(i-1:i+3,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))

            i = reg_l1 + 2
            lamfil(i) = dot_product(ff3b2, lam(i-2:i+3,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))
         end if

         if (Er(reg_h1+1,g) .eq. -1.e0_rt) then
            i = reg_h1 - 2
            lamfil(i) = dot_product(ff3b2(3:-2:-1), lam(i-3:i+2,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))

            i = reg_h1 - 1
            lamfil(i) = dot_product(ff3b1(3:-1:-1), lam(i-3:i+1,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))

            i = reg_h1
            lamfil(i) = dot_product(ff3b0(3:0:-1), lam(i-3:i,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))
         end if

         lam(reg_l1:reg_h1,g) = lamfil(reg_l1:reg_h1)

      else if (filter_T .eq. 4) then

         do i=reg_l1, reg_h1
            lamfil(i) = ff4(0,S) * lam(i,g) &
                 &    + ff4(1,S) * (lam(i-1,g)+lam(i+1,g)) &
                 &    + ff4(2,S) * (lam(i-2,g)+lam(i+2,g)) &
                 &    + ff4(3,S) * (lam(i-3,g)+lam(i+3,g)) &
                 &    + ff4(4,S) * (lam(i-4,g)+lam(i+4,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))
         end do

         if (Er(reg_l1-1,g) .eq. -1.e0_rt) then
            i = reg_l1
            lamfil(i) = dot_product(ff4b0, lam(i:i+4,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))

            i = reg_l1 + 1
            lamfil(i) = dot_product(ff4b1, lam(i-1:i+4,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))

            i = reg_l1 + 2
            lamfil(i) = dot_product(ff4b2, lam(i-2:i+4,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))

            i = reg_l1 + 3
            lamfil(i) = dot_product(ff4b3, lam(i-3:i+4,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))
         end if

         if (Er(reg_h1+1,g) .eq. -1.e0_rt) then
            i = reg_h1 - 3
            lamfil(i) = dot_product(ff4b3(4:-3:-1), lam(i-4:i+3,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))

            i = reg_h1 - 2
            lamfil(i) = dot_product(ff4b2(4:-2:-1), lam(i-4:i+2,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))

            i = reg_h1 - 1
            lamfil(i) = dot_product(ff4b1(4:-1:-1), lam(i-4:i+1,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))

            i = reg_h1
            lamfil(i) = dot_product(ff4b0(4:0:-1), lam(i-4:i,g))
            lamfil(i) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lamfil(i)))
         end if

         lam(reg_l1:reg_h1,g) = lamfil(reg_l1:reg_h1)

     end if

     ! boundary

     if (Er(reg_l1-1,g) .eq. -1.e0_rt) then
        do i=lam_l1, reg_l1-1
           lam(i,g) = lam(reg_l1,g)
        end do
     end if

     if (Er(reg_h1+1,g) .eq. -1.e0_rt) then
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
     dcf, dcf_l1, dcf_h1) bind(C, name="ca_get_v_dcf")

  use meth_params_module, only : NVAR, URHO, UMX

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: er_l1,er_h1,s_l1,s_h1,T_l1,T_h1,c_v_l1,c_v_h1
  integer, intent(in) :: kr_l1,kr_h1,kp_l1,kp_h1,kp2_l1,kp2_h1
  integer, intent(in) :: v_l1,v_h1,dcf_l1,dcf_h1
  real(rt)        , intent(in)  ::  er( er_l1: er_h1)
  real(rt)        , intent(in)  ::   s(  s_l1:  s_h1, NVAR)
  real(rt)        , intent(in)  ::   T(  T_l1:  T_h1)
  real(rt)        , intent(in)  :: c_v(c_v_l1:c_v_h1)
  real(rt)        , intent(in ) ::  kr( kr_l1: kr_h1)
  real(rt)        , intent(in ) ::  kp( kp_l1: kp_h1)
  real(rt)        , intent(in ) :: kp2(kp2_l1:kp2_h1)
  real(rt)        , intent(in) :: dtemp, dtime, sigma, c
  real(rt)         ::   v(  v_l1:  v_h1)
  real(rt)         :: dcf(dcf_l1:dcf_h1)

  integer :: i
  real(rt)         :: etainv, fac0, fac2, alpha, frc

  fac0 = 4.e0_rt * sigma * dtime / dtemp
  fac2 = c * dtime / dtemp

  do i=lo(1),hi(1)
     v(i) = s(i,UMX)/s(i,URHO)

     alpha = fac0 * (kp2(i) * (T(i) + dtemp) ** 4    &
          -          kp (i) * (T(i)        ) ** 4)   &
          -  fac2 * (kp2(i) - kp(i)) * er(i)

     frc = s(i,URHO) * c_v(i) + 1.0e-50_rt
     etainv = frc / (alpha + frc)

     dcf(i) = 2.e0_rt * etainv * (kp(i) / kr(i))
  end do

end subroutine ca_get_v_dcf

subroutine ca_update_dcf(lo, hi, &
     dcf, dcf_l1, dcf_h1, &
     etainv, eti_l1, eti_h1, &
     kp, kp_l1, kp_h1, kr, kr_l1, kr_h1) bind(C, name="ca_update_dcf")

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: dcf_l1, dcf_h1, eti_l1, eti_h1, &
       kp_l1, kp_h1, kr_l1, kr_h1
  real(rt)        , intent(in) :: etainv(eti_l1:eti_h1)
  real(rt)        , intent(in) :: kp(kp_l1:kp_h1)
  real(rt)        , intent(in) :: kr(kr_l1:kr_h1)
  real(rt)                     :: dcf(dcf_l1:dcf_h1)

  integer :: i

  do i=lo(1),hi(1)
     dcf(i) = 2.e0_rt * etainv(i) * (kp(i)/kr(i))
  end do

end subroutine ca_update_dcf

subroutine ca_set_dterm_face( lo, hi, &
     Er, Er_l1, Er_h1, dc, dc_l1, dc_h1, &
     dtf, dtf_l1, dtf_h1, dx, idir) bind(C, name="ca_set_dterm_face")

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: Er_l1, Er_h1, dc_l1, dc_h1, dtf_l1, dtf_h1, idir
  real(rt)        , intent(in) :: dx
  real(rt)        , intent(in) :: Er(Er_l1:Er_h1)
  real(rt)        , intent(in) :: dc(dc_l1:dc_h1)
  real(rt)                     :: dtf(dtf_l1:dtf_h1)
  integer :: i

  do i=lo(1),hi(1)
     dtf(i) = (Er(i) - Er(i-1)) / dx * dc(i)
  end do

end subroutine ca_set_dterm_face

subroutine ca_face2center( lo, hi, &
     scomp, dcomp, ncomp, nx, nc, &
     foox, foox_l1, foox_h1, &
     fooc, fooc_l1, fooc_h1) bind(C, name="ca_face2center")

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1), scomp,dcomp,ncomp,nx,nc
  integer, intent(in) :: foox_l1, foox_h1
  integer, intent(in) :: fooc_l1, fooc_h1
  real(rt)        , intent(in) :: foox(foox_l1:foox_h1,0:nx-1)
  real(rt)                     :: fooc(fooc_l1:fooc_h1,0:nc-1)

  integer :: i, n

  do n = 0, ncomp-1
     do i=lo(1), hi(1)
        fooc(i,dcomp+n) = (foox(i,scomp+n) + foox(i+1,scomp+n)) * 0.5e0_rt;
     end do
  end do

end subroutine ca_face2center

! no tiling
subroutine ca_correct_dterm(dfx, dfx_l1, dfx_h1, &
     re, rc) bind(C, name="ca_correct_dterm")

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: dfx_l1, dfx_h1
  real(rt)        , intent(inout) :: dfx(dfx_l1:dfx_h1)
  real(rt)        , intent(in) :: re(dfx_l1:dfx_h1), rc(1)

  integer :: i

  do i=dfx_l1, dfx_h1
     dfx(i) = dfx(i) / (re(i) + 1.e-50_rt)
  end do

end subroutine ca_correct_dterm

subroutine ca_estdt_rad(lo, hi, u,u_l1,u_h1, gpr,gpr_l1,gpr_h1, &
  dx,dt) bind(C)

  use network, only : nspec, naux
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_re
  use meth_params_module, only : NVAR, URHO, UMX, UEINT, UTEMP, UFS, UFX
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer u_l1,u_h1
  integer gpr_l1,gpr_h1
  integer lo(1), hi(1)
  real(rt)         u(u_l1:u_h1,NVAR)
  real(rt)         gpr(gpr_l1:gpr_h1)
  real(rt)         dx(1), dt

  real(rt)         :: rhoInv,ux,dt1,c
  integer          :: i
  type(eos_t) :: eos_state

  !     Translate to primitive variables, compute sound speed (call eos), get dtmax
  do i = lo(1),hi(1)

     rhoInv = 1.e0_rt / u(i,URHO)

     eos_state % rho = u(i,URHO)
     eos_state % T   = u(i,UTEMP)
     eos_state % e   = u(i,UEINT)*rhoInv
     eos_state % xn  = u(i,UFS:UFS+nspec-1) * rhoInv
     eos_state % aux = u(i,UFX:UFX+naux-1) * rhoInv

     call eos(eos_input_re, eos_state)
     c = eos_state % cs

     c = sqrt(c**2 + gpr(i)*rhoInv)

     ux = u(i,UMX)*rhoInv

     dt1 = dx(1) /( c + abs(ux) )
     dt  = min(dt,dt1)

  enddo

end subroutine ca_estdt_rad


! this is tiling safe
subroutine ca_est_gpr0(Er, Er_l1, Er_h1, gPr, gPr_l1, gPr_h1)

  use rad_params_module, only : ngroups

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: Er_l1, Er_h1, gpr_l1, gpr_h1
  real(rt)        , intent(in) :: Er(Er_l1:Er_h1, 0:ngroups-1)
  real(rt)        , intent(out) :: gPr(gPr_l1:gPr_h1)

  integer :: i, g

  gPr = 0.e0_rt

  do g = 0, ngroups-1
     do i = gPr_l1, gPr_h1
        gPr(i) = gPr(i) + 4.e0_rt/9.e0_rt*Er(i,g)
     end do
  end do

end subroutine ca_est_gpr0


! this is tiling safe
subroutine ca_est_gpr2(kap, kap_l1, kap_h1, Er, Er_l1, Er_h1, &
     gPr, gPr_l1, gPr_h1, vlo, vhi, dx, limiter, comoving)

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDlambda, Edd_factor

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: kap_l1, kap_h1
  integer, intent(in) ::  Er_l1,  Er_h1
  integer, intent(in) :: gPr_l1, gPr_h1
  integer, intent(in) :: vlo(1), vhi(1)  ! region with valid Er
  integer, intent(in) :: limiter, comoving
  real(rt)        , intent(in) :: dx
  real(rt)        , intent(in) ::  kap(kap_l1:kap_h1, 0:ngroups-1)
  real(rt)        , intent(in) ::   Er( Er_l1: Er_h1, 0:ngroups-1)
  real(rt)        , intent(out) :: gPr(gPr_l1:gPr_h1)

  integer :: i, g
  real(rt)         :: r, lam, f, gamr
  integer :: im, ip
  real(rt)         :: xm, xp

  if (gPr_l1-1 .ge. vlo(1)) then
     im = 1
     xm = 2.e0_rt
  else
     im = 0
     xm = 1.e0_rt
  end if

  if (gPr_h1+1 .le. vhi(1)) then
     ip = 1
     xp = 2.e0_rt
  else
     ip = 0
     xp = 1.e0_rt
  end if

  gPr = 0.0e0_rt

  do g = 0, ngroups-1
     i = gPr_l1
     r = abs(Er(i+1,g) - Er(i-im,g)) / (xm*dx)
     r = r / (kap(i,g) * max(Er(i,g), 1.e-50_rt))
     lam = FLDlambda(r, limiter)
     if (comoving .eq. 1) then
        f = Edd_factor(lam)
        gamr = (3.e0_rt-f)/2.e0_rt
     else
        gamr = lam + 1.e0_rt
     end if
     gPr(i) = gPr(i) + gamr * lam * Er(i,g)

     do i = gPr_l1+1, gPr_h1-1
        r = abs(Er(i+1,g) - Er(i-1,g)) / (2.e0_rt*dx)
        r = r / (kap(i,g) * max(Er(i,g), 1.e-50_rt))
        lam = FLDlambda(r, limiter)
        if (comoving .eq. 1) then
           f = Edd_factor(lam)
           gamr = (3.e0_rt-f)/2.e0_rt
        else
           gamr = lam + 1.e0_rt
        end if
        gPr(i) = gPr(i) + gamr * lam * Er(i,g)
     end do

     i = gPr_h1
     r = abs(Er(i+ip,g) - Er(i-1,g)) / (xp*dx)
     r = r / (kap(i,g) * max(Er(i,g), 1.e-50_rt))
     lam = FLDlambda(r, limiter)
     if (comoving .eq. 1) then
        f = Edd_factor(lam)
        gamr = (3.e0_rt-f)/2.e0_rt
     else
        gamr = lam + 1.e0_rt
     end if
     gPr(i) = gPr(i) + gamr * lam * Er(i,g)
  end do

end subroutine ca_est_gpr2
