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
