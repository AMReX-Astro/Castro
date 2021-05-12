! This subroutine cannot be tiled
subroutine ca_compute_lamborder(Er, Er_l1, Er_l2, Er_h1, Er_h2, &
                                kap, kap_l1, kap_l2, kap_h1, kap_h2, &
                                lam, lam_l1, lam_l2, lam_h1, lam_h2, &
                                dx, ngrow, limiter, filter_T, S)

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDlambda
  use filter_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: Er_l1, Er_l2, Er_h1, Er_h2, kap_l1, kap_l2, kap_h1, kap_h2, &
       lam_l1, lam_l2, lam_h1, lam_h2
  integer, intent(in) :: ngrow, limiter, filter_T, S
  real(rt)        , intent(in) :: dx(2)
  real(rt)        , intent(in) :: kap(kap_l1:kap_h1, kap_l2:kap_h2)
  real(rt)        , intent(in) :: Er(Er_l1:Er_h1, Er_l2:Er_h2, 0:ngroups-1)
  real(rt)        , intent(out) :: lam(lam_l1:lam_h1, lam_l2:lam_h2, 0:ngroups-1)

  integer :: i, j, reg_l1, reg_l2, reg_h1, reg_h2, g
  real(rt)         :: r, r1, r2

  real(rt)        , allocatable :: lamfil(:,:)

  lam = -1.e50_rt

  reg_l1 = lam_l1 + ngrow
  reg_l2 = lam_l2 + ngrow
  reg_h1 = lam_h1 - ngrow
  reg_h2 = lam_h2 - ngrow

  if (filter_T .gt. 0) then
     allocate(lamfil(reg_l1:reg_h1,lam_l2:lam_h2))
  end if

  do g = 0, ngroups-1

  do j=lam_l2, lam_h2
     do i=lam_l1, lam_h1
        if (Er(i,j,g) .eq. -1.e0_rt) then
           cycle
        end if

        if (Er(i-1,j,g) .eq. -1.e0_rt) then
           r1 = (Er(i+1,j,g) - Er(i,j,g)) / (dx(1))
        else if (Er(i+1,j,g) .eq. -1.e0_rt) then
           r1 = (Er(i,j,g) - Er(i-1,j,g)) / (dx(1))
        else
           r1 = (Er(i+1,j,g) - Er(i-1,j,g)) / (2.e0_rt*dx(1))
        end if

        if (Er(i,j-1,g) .eq. -1.e0_rt) then
           r2 = (Er(i,j+1,g) - Er(i,j,g)) / (dx(2))
        else if (Er(i,j+1,g) .eq. -1.e0_rt) then
           r2 = (Er(i,j,g) - Er(i,j-1,g)) / (dx(2))
        else
           r2 = (Er(i,j+1,g) - Er(i,j-1,g)) / (2.e0_rt*dx(2))
        end if

        r = sqrt(r1**2 + r2**2)
        r = r / (kap(i,j) * max(Er(i,j,g), 1.e-50_rt))

        lam(i,j,g) = FLDlambda(r, limiter)
     end do
  end do

  ! filter
  if (filter_T .eq. 1) then

     do j=lam_l2, lam_h2
        if (Er(reg_l1,j,g) .eq. -1.e0_rt) then
           lamfil(:,j) = -1.e-50_rt
           cycle
        endif

        do i=reg_l1, reg_h1
           lamfil(i,j) = ff1(0) * lam(i,j,g) &
                &      + ff1(1) * (lam(i-1,j,g)+lam(i+1,j,g))
        end do

        if (Er(reg_l1-1,j,g) .eq. -1.e0_rt) then
            i = reg_l1
            lamfil(i,j) = dot_product(ff1b, lam(i:i+1,j,g))
         end if

         if (Er(reg_h1+1,j,g) .eq. -1.e0_rt) then
            i = reg_h1
            lamfil(i,j) = dot_product(ff1b(1:0:-1), lam(i-1:i,j,g))
         end if
     end do

     do j=reg_l2, reg_h2
        do i=reg_l1, reg_h1
           lam(i,j,g) = ff1(0) * lamfil(i,j) &
                &     + ff1(1) * (lamfil(i,j-1)+lamfil(i,j+1))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do
     end do

     if (Er(reg_l1,reg_l2-1,g) .eq. -1.e0_rt) then
        j = reg_l2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff1b, lamfil(i,j:j+1))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do
     end if

     if (Er(reg_l1,reg_h2+1,g) .eq. -1.e0_rt) then
        j = reg_h2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff1b(1:0:-1), lamfil(i,j-1:j))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do
     end if

  else if (filter_T .eq. 2) then

     do j=lam_l2, lam_h2
        if (Er(reg_l1,j,g) .eq. -1.e0_rt) then
           lamfil(:,j) = -1.e-50_rt
           cycle
        endif

        do i=reg_l1, reg_h1
           lamfil(i,j) = ff2(0,S) * lam(i,j,g) &
                &      + ff2(1,S) * (lam(i-1,j,g)+lam(i+1,j,g)) &
                &      + ff2(2,S) * (lam(i-2,j,g)+lam(i+2,j,g))
        end do

        if (Er(reg_l1-1,j,g) .eq. -1.e0_rt) then
            i = reg_l1
            lamfil(i,j) = dot_product(ff2b0, lam(i:i+2,j,g))

            i = reg_l1 + 1
            lamfil(i,j) = dot_product(ff2b1, lam(i-1:i+2,j,g))
         end if

         if (Er(reg_h1+1,j,g) .eq. -1.e0_rt) then
            i = reg_h1 - 1
            lamfil(i,j) = dot_product(ff2b1(2:-1:-1), lam(i-2:i+1,j,g))

            i = reg_h1
            lamfil(i,j) = dot_product(ff2b0(2:0:-1), lam(i-2:i,j,g))
         end if
     end do

     do j=reg_l2, reg_h2
        do i=reg_l1, reg_h1
           lam(i,j,g) = ff2(0,S) * lamfil(i,j) &
                &     + ff2(1,S) * (lamfil(i,j-1)+lamfil(i,j+1)) &
                &     + ff2(2,S) * (lamfil(i,j-2)+lamfil(i,j+2))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do
     end do

     if (Er(reg_l1,reg_l2-1,g) .eq. -1.e0_rt) then
        j = reg_l2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff2b0, lamfil(i,j:j+2))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do

        j = reg_l2 + 1
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff2b1, lamfil(i,j-1:j+2))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do
     end if

     if (Er(reg_l1,reg_h2+1,g) .eq. -1.e0_rt) then
        j = reg_h2 - 1
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff2b1(2:-1:-1), lamfil(i,j-2:j+1))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do

        j = reg_h2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff2b0(2:0:-1), lamfil(i,j-2:j))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do
     end if

  else if (filter_T .eq. 3) then

     do j=lam_l2, lam_h2
        if (Er(reg_l1,j,g) .eq. -1.e0_rt) then
           lamfil(:,j) = -1.e-50_rt
           cycle
        endif

        do i=reg_l1, reg_h1
           lamfil(i,j) = ff3(0,S) * lam(i,j,g) &
                &      + ff3(1,S) * (lam(i-1,j,g)+lam(i+1,j,g)) &
                &      + ff3(2,S) * (lam(i-2,j,g)+lam(i+2,j,g)) &
                &      + ff3(3,S) * (lam(i-3,j,g)+lam(i+3,j,g))
        end do

        if (Er(reg_l1-1,j,g) .eq. -1.e0_rt) then
            i = reg_l1
            lamfil(i,j) = dot_product(ff3b0, lam(i:i+3,j,g))

            i = reg_l1 + 1
            lamfil(i,j) = dot_product(ff3b1, lam(i-1:i+3,j,g))

            i = reg_l1 + 2
            lamfil(i,j) = dot_product(ff3b2, lam(i-2:i+3,j,g))
         end if

         if (Er(reg_h1+1,j,g) .eq. -1.e0_rt) then
            i = reg_h1 - 2
            lamfil(i,j) = dot_product(ff3b2(3:-2:-1), lam(i-3:i+2,j,g))

            i = reg_h1 - 1
            lamfil(i,j) = dot_product(ff3b1(3:-1:-1), lam(i-3:i+1,j,g))

            i = reg_h1
            lamfil(i,j) = dot_product(ff3b0(3:0:-1), lam(i-3:i,j,g))
         end if
     end do

     do j=reg_l2, reg_h2
        do i=reg_l1, reg_h1
           lam(i,j,g) = ff3(0,S) * lamfil(i,j) &
                &     + ff3(1,S) * (lamfil(i,j-1)+lamfil(i,j+1)) &
                &     + ff3(2,S) * (lamfil(i,j-2)+lamfil(i,j+2)) &
                &     + ff3(3,S) * (lamfil(i,j-3)+lamfil(i,j+3))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do
     end do

     if (Er(reg_l1,reg_l2-1,g) .eq. -1.e0_rt) then
        j = reg_l2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff3b0, lamfil(i,j:j+3))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do

        j = reg_l2 + 1
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff3b1, lamfil(i,j-1:j+3))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do

        j = reg_l2 + 2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff3b2, lamfil(i,j-2:j+3))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do
     end if

     if (Er(reg_l1,reg_h2+1,g) .eq. -1.e0_rt) then
        j = reg_h2 - 2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff3b2(3:-2:-1), lamfil(i,j-3:j+2))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do

        j = reg_h2 - 1
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff3b1(3:-1:-1), lamfil(i,j-3:j+1))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do

        j = reg_h2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff3b0(3:0:-1), lamfil(i,j-3:j))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do
     end if

  else if (filter_T .eq. 4) then

     do j=lam_l2, lam_h2
        if (Er(reg_l1,j,g) .eq. -1.e0_rt) then
           lamfil(:,j) = -1.e-50_rt
           cycle
        endif

        do i=reg_l1, reg_h1
           lamfil(i,j) = ff4(0,S) * lam(i,j,g) &
                &      + ff4(1,S) * (lam(i-1,j,g)+lam(i+1,j,g)) &
                &      + ff4(2,S) * (lam(i-2,j,g)+lam(i+2,j,g)) &
                &      + ff4(3,S) * (lam(i-3,j,g)+lam(i+3,j,g)) &
                &      + ff4(4,S) * (lam(i-4,j,g)+lam(i+4,j,g))
        end do

        if (Er(reg_l1-1,j,g) .eq. -1.e0_rt) then
            i = reg_l1
            lamfil(i,j) = dot_product(ff4b0, lam(i:i+4,j,g))

            i = reg_l1 + 1
            lamfil(i,j) = dot_product(ff4b1, lam(i-1:i+4,j,g))

            i = reg_l1 + 2
            lamfil(i,j) = dot_product(ff4b2, lam(i-2:i+4,j,g))

            i = reg_l1 + 3
            lamfil(i,j) = dot_product(ff4b3, lam(i-3:i+4,j,g))
         end if

         if (Er(reg_h1+1,j,g) .eq. -1.e0_rt) then
            i = reg_h1 - 3
            lamfil(i,j) = dot_product(ff4b3(4:-3:-1), lam(i-4:i+3,j,g))

            i = reg_h1 - 2
            lamfil(i,j) = dot_product(ff4b2(4:-2:-1), lam(i-4:i+2,j,g))

            i = reg_h1 - 1
            lamfil(i,j) = dot_product(ff4b1(4:-1:-1), lam(i-4:i+1,j,g))

            i = reg_h1
            lamfil(i,j) = dot_product(ff4b0(4:0:-1), lam(i-4:i,j,g))
         end if
     end do

     do j=reg_l2, reg_h2
        do i=reg_l1, reg_h1
           lam(i,j,g) = ff4(0,S) * lamfil(i,j) &
                &     + ff4(1,S) * (lamfil(i,j-1)+lamfil(i,j+1)) &
                &     + ff4(2,S) * (lamfil(i,j-2)+lamfil(i,j+2)) &
                &     + ff4(3,S) * (lamfil(i,j-3)+lamfil(i,j+3)) &
                &     + ff4(4,S) * (lamfil(i,j-4)+lamfil(i,j+4))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do
     end do

     if (Er(reg_l1,reg_l2-1,g) .eq. -1.e0_rt) then
        j = reg_l2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b0, lamfil(i,j:j+4))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do

        j = reg_l2 + 1
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b1, lamfil(i,j-1:j+4))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do

        j = reg_l2 + 2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b2, lamfil(i,j-2:j+4))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do

        j = reg_l2 + 3
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b3, lamfil(i,j-3:j+4))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do
     end if

     if (Er(reg_l1,reg_h2+1,g) .eq. -1.e0_rt) then
        j = reg_h2 - 3
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b3(4:-3:-1), lamfil(i,j-4:j+3))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do

        j = reg_h2 - 2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b2(4:-2:-1), lamfil(i,j-4:j+2))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do

        j = reg_h2 - 1
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b1(4:-1:-1), lamfil(i,j-4:j+1))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do

        j = reg_h2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b0(4:0:-1), lamfil(i,j-4:j))
           lam(i,j,g) = min(1.e0_rt/3.e0_rt, max(1.e-25_rt, lam(i,j,g)))
        end do
     end if

  end if

  ! lo-x lo-y
  do j=lam_l2,reg_l2-1
     do i=lam_l1,reg_l1-1
        if (Er(i,j,g).eq.-1.e0_rt) then
           lam(i,j,g) = lam(reg_l1,reg_l2,g)
        end if
     end do
  end do

  ! reg-x lo-y
  do j=lam_l2,reg_l2-1
     do i=reg_l1,reg_h1
        if (Er(i,j,g).eq.-1.e0_rt) then
           lam(i,j,g) = lam(i,reg_l2,g)
        end if
     end do
  end do

  ! hi-x lo-y
  do j=lam_l2,reg_l2-1
     do i=reg_h1+1,lam_h1
        if (Er(i,j,g).eq.-1.e0_rt) then
           lam(i,j,g) = lam(reg_h1,reg_l2,g)
        end if
     end do
  end do

  ! lo-x reg-y
  do j=reg_l2,reg_h2
     do i=lam_l1,reg_l1-1
        if (Er(i,j,g).eq.-1.e0_rt) then
           lam(i,j,g) = lam(reg_l1,j,g)
        end if
     end do
  end do

  ! hi-x reg-y
  do j=reg_l2,reg_h2
     do i=reg_h1+1,lam_h1
        if (Er(i,j,g).eq.-1.e0_rt) then
           lam(i,j,g) = lam(reg_h1,j,g)
        end if
     end do
  end do

  ! lo-x hi-y
  do j=reg_h2+1,lam_h2
     do i=lam_l1,reg_l1-1
        if (Er(i,j,g).eq.-1.e0_rt) then
           lam(i,j,g) = lam(reg_l1,reg_h2,g)
        end if
     end do
  end do

  ! reg-x hi-y
  do j=reg_h2+1,lam_h2
     do i=reg_l1,reg_h1
        if (Er(i,j,g).eq.-1.e0_rt) then
           lam(i,j,g) = lam(i,reg_h2,g)
        end if
     end do
  end do

  ! hi-x hi-y
  do j=reg_h2+1,lam_h2
     do i=reg_h1+1,lam_h1
        if (Er(i,j,g).eq.-1.e0_rt) then
           lam(i,j,g) = lam(reg_h1,reg_h2,g)
        end if
     end do
  end do

  end do

  if (filter_T .gt. 0) then
     deallocate(lamfil)
  end if

  return
end subroutine ca_compute_lamborder

subroutine ca_set_dterm_face(lo, hi, &
                             Er, Er_l1, Er_l2, Er_h1, Er_h2, &
                             dc, dc_l1, dc_l2, dc_h1, dc_h2, &
                             dtf, dtf_l1, dtf_l2, dtf_h1, dtf_h2, dx, idir) bind(C, name="ca_set_dterm_face")

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: Er_l1, Er_l2, Er_h1, Er_h2,  &
       dc_l1, dc_l2, dc_h1, dc_h2, dtf_l1, dtf_l2, dtf_h1, dtf_h2, idir
  real(rt)        , intent(in) :: dx(2)
  real(rt)        , intent(in) :: Er(Er_l1:Er_h1,Er_l2:Er_h2)
  real(rt)        , intent(in) :: dc(dc_l1:dc_h1,dc_l2:dc_h2)
  real(rt)                     :: dtf(dtf_l1:dtf_h1,dtf_l2:dtf_h2)
  integer :: i, j

  if (idir .eq. 0) then
     do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           dtf(i,j) = (Er(i,j) - Er(i-1,j)) / dx(1) * dc(i,j)
        end do
     end do
  else
     do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           dtf(i,j) = (Er(i,j) - Er(i,j-1)) / dx(2) * dc(i,j)
        end do
     end do
  end if

end subroutine ca_set_dterm_face

! no tiling
subroutine ca_correct_dterm(  &
                            dfx, dfx_l1, dfx_l2, dfx_h1, dfx_h2, &
                            dfy, dfy_l1, dfy_l2, dfy_h1, dfy_h2, &
                            re, rc) bind(C, name="ca_correct_dterm")

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: dfx_l1, dfx_l2, dfx_h1, dfx_h2
  integer, intent(in) :: dfy_l1, dfy_l2, dfy_h1, dfy_h2
  real(rt)        , intent(inout) :: dfx(dfx_l1:dfx_h1,dfx_l2:dfx_h2)
  real(rt)        , intent(inout) :: dfy(dfy_l1:dfy_h1,dfy_l2:dfy_h2)
  real(rt)        , intent(in) :: re(dfx_l1:dfx_h1), rc(dfy_l1:dfy_h1)

  integer :: i, j

  do j=dfx_l2, dfx_h2
     do i=dfx_l1, dfx_h1
        dfx(i,j) = dfx(i,j) / (re(i) + 1.e-50_rt)
     end do
  end do

  do j=dfy_l2, dfy_h2
     do i=dfy_l1, dfy_h1
        dfy(i,j) = dfy(i,j) / rc(i)
     end do
  end do

end subroutine ca_correct_dterm
