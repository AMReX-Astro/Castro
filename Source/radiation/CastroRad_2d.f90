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


subroutine ca_get_v_dcf(lo, hi, &
                        er ,  er_l1,  er_l2,  er_h1,  er_h2, &
                        s  ,   s_l1,   s_l2,   s_h1,   s_h2, &
                        T  ,   T_l1,   T_l2,   T_h1,   T_h2, &
                        c_v, c_v_l1, c_v_l2, c_v_h1, c_v_h2, &
                        kr ,  kr_l1,  kr_l2,  kr_h1,  kr_h2, &
                        kp ,  kp_l1,  kp_l2,  kp_h1,  kp_h2, &
                        kp2, kp2_l1, kp2_l2, kp2_h1, kp2_h2, &
                        dtemp, dtime, sigma, c, &
                        v  ,   v_l1,   v_l2,   v_h1,   v_h2, &
                        dcf, dcf_l1, dcf_l2, dcf_h1, dcf_h2) bind(C, name="ca_get_v_dcf")

  use meth_params_module, only : NVAR, URHO, UMX, UMY

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: er_l1,er_l2,er_h1,er_h2,s_l1,s_l2,s_h1,s_h2
  integer, intent(in) :: T_l1,T_l2,T_h1,T_h2,c_v_l1,c_v_l2,c_v_h1,c_v_h2
  integer, intent(in) :: kr_l1,kr_l2,kr_h1,kr_h2
  integer, intent(in) :: kp_l1,kp_l2,kp_h1,kp_h2,kp2_l1,kp2_l2,kp2_h1,kp2_h2
  integer, intent(in) :: v_l1,v_l2,v_h1,v_h2,dcf_l1,dcf_l2,dcf_h1,dcf_h2
  real(rt)        , intent(in)  ::  er( er_l1: er_h1,  er_l2: er_h2)
  real(rt)        , intent(in)  ::   s(  s_l1:  s_h1,   s_l2:  s_h2, NVAR)
  real(rt)        , intent(in)  ::   T(  T_l1:  T_h1,   T_l2:  T_h2)
  real(rt)        , intent(in)  :: c_v(c_v_l1:c_v_h1, c_v_l2:c_v_h2)
  real(rt)        , intent(in ) ::  kr( kr_l1: kr_h1,  kr_l2: kr_h2)
  real(rt)        , intent(in ) ::  kp( kp_l1: kp_h1,  kp_l2: kp_h2)
  real(rt)        , intent(in ) :: kp2(kp2_l1:kp2_h1, kp2_l2:kp2_h2)
  real(rt)        , intent(in) :: dtemp, dtime, sigma, c
  real(rt)                      ::   v(  v_l1:  v_h1,   v_l2:  v_h2, 2)
  real(rt)                      :: dcf(dcf_l1:dcf_h1, dcf_l2:dcf_h2)

  integer :: i, j
  real(rt)         :: etainv, fac0, fac2, alpha, frc

  fac0 = 4.e0_rt * sigma * dtime / dtemp
  fac2 = c * dtime / dtemp

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        v(i,j,1) = s(i,j,UMX)/s(i,j,URHO)
        v(i,j,2) = s(i,j,UMY)/s(i,j,URHO)

        alpha = fac0 * (kp2(i,j) * (T(i,j) + dtemp) ** 4    &
             -          kp (i,j) * (T(i,j)        ) ** 4)   &
             -  fac2 * (kp2(i,j) - kp(i,j)) * er(i,j)

        frc = s(i,j,URHO) * c_v(i,j) + 1.0e-50_rt
        etainv = frc / (alpha + frc)

        dcf(i,j) = 2.e0_rt * etainv * (kp(i,j) / kr(i,j))
     end do
  end do

end subroutine ca_get_v_dcf

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


subroutine ca_face2center(lo, hi, &
                          scomp, dcomp, ncomp, nf, nc, &
                          foox, foox_l1, foox_l2, foox_h1, foox_h2, &
                          fooy, fooy_l1, fooy_l2, fooy_h1, fooy_h2, &
                          fooc, fooc_l1, fooc_l2, fooc_h1, fooc_h2) bind(C, name="ca_face2center")

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(2), hi(2), scomp,dcomp,ncomp,nf,nc
  integer, intent(in) :: foox_l1, foox_l2, foox_h1, foox_h2
  integer, intent(in) :: fooy_l1, fooy_l2, fooy_h1, fooy_h2
  integer, intent(in) :: fooc_l1, fooc_l2, fooc_h1, fooc_h2
  real(rt)        , intent(in)  :: foox(foox_l1:foox_h1,foox_l2:foox_h2,0:nf-1)
  real(rt)        , intent(in)  :: fooy(fooy_l1:fooy_h1,fooy_l2:fooy_h2,0:nf-1)
  real(rt)                      :: fooc(fooc_l1:fooc_h1,fooc_l2:fooc_h2,0:nc-1)

  integer :: i,j,n

  do n = 0, ncomp-1
     do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           fooc(i,j,dcomp+n) = (foox(i,j,scomp+n) + foox(i+1,j,scomp+n) &
                &             + fooy(i,j,scomp+n) + fooy(i,j+1,scomp+n)) * 0.25e0_rt
        end do
     end do
  end do

end subroutine ca_face2center


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


subroutine ca_estdt_rad(lo, hi, u,u_l1,u_l2,u_h1,u_h2, &
                        gpr,gpr_l1,gpr_l2,gpr_h1,gpr_h2, &
                        dx,dt) bind(C)

  use network, only : nspec, naux
  use eos_module, only : eos
  use eos_type_module, only: eos_t, eos_input_re
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: u_l1,u_l2,u_h1,u_h2
  integer          :: gpr_l1,gpr_l2,gpr_h1,gpr_h2
  integer          :: lo(2), hi(2)
  real(rt)         :: u(u_l1:u_h1,u_l2:u_h2,NVAR)
  real(rt)         :: gpr(gpr_l1:gpr_h1,gpr_l2:gpr_h2)
  real(rt)         :: dx(2),dt

  real(rt)         :: rhoInv,ux,uy,dt1,dt2,c
  integer          :: i,j
  type(eos_t) :: eos_state

  !    Translate to primitive variables, compute sound speed (call eos), get dtmax
  do j = lo(2),hi(2)
     do i = lo(1),hi(1)

        rhoInv = 1.e0_rt / u(i,j,URHO)

        eos_state % rho = u(i,j,URHO)
        eos_state % T   = u(i,j,UTEMP)
        eos_state % e   = u(i,j,UEINT)*rhoInv
        eos_state % xn  = u(i,j,UFS:UFS+nspec-1) * rhoInv
        eos_state % aux = u(i,j,UFX:UFX+naux -1) * rhoInv

        call eos(eos_input_re, eos_state)
        c = eos_state % cs

        c = sqrt(c**2 + gpr(i,j)*rhoInv)

        ux = u(i,j,UMX)*rhoInv
        uy = u(i,j,UMY)*rhoInv

        dt1 = dx(1)/(c + abs(ux))
        dt2 = dx(2)/(c + abs(uy))
        dt = min(dt,dt1,dt2)
     enddo
  enddo

end subroutine ca_estdt_rad


! this is tiling safe
subroutine ca_est_gpr0(Er, Er_l1, Er_l2, Er_h1, Er_h2, &
                       gPr, gPr_l1, gPr_l2, gPr_h1, gPr_h2)

  use rad_params_module, only : ngroups

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: Er_l1, Er_l2, Er_h1, Er_h2, gpr_l1, gpr_l2, gpr_h1, gpr_h2
  real(rt)        , intent(in) :: Er(Er_l1:Er_h1,Er_l2:Er_h2, 0:ngroups-1)
  real(rt)        , intent(out) :: gPr(gPr_l1:gPr_h1,gPr_l2:gPr_h2)

  integer :: i, j, g

  gPr = 0.e0_rt

  do g = 0, ngroups-1
     do j = gPr_l2, gPr_h2
        do i = gPr_l1, gPr_h1
           gPr(i,j) = gPr(i,j) + 4.e0_rt/9.e0_rt*Er(i,j,g)
        end do
     end do
  end do

end subroutine ca_est_gpr0


! this is tiling safe
subroutine ca_est_gpr2(kap, kap_l1, kap_l2, kap_h1, kap_h2, &
                       Er, Er_l1, Er_l2, Er_h1, Er_h2, &
                       gPr, gPr_l1, gPr_l2, gPr_h1, gPr_h2, vlo, vhi, dx, limiter, comoving)

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDlambda, Edd_factor

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: kap_l1, kap_l2, kap_h1, kap_h2, &
       Er_l1, Er_l2, Er_h1, Er_h2, gPr_l1, gPr_l2, gPr_h1, gPr_h2
  integer, intent(in) :: vlo(2), vhi(2)  ! the region with valid Er
  integer, intent(in) :: limiter, comoving
  real(rt)        , intent(in) :: dx(2)
  real(rt)        , intent(in) :: kap(kap_l1:kap_h1,kap_l2:kap_h2, 0:ngroups-1), &
       Er(Er_l1:Er_h1,Er_l2:Er_h2, 0:ngroups-1)
  real(rt)        , intent(out) :: gPr(gPr_l1:gPr_h1,gPr_l2:gPr_h2)

  integer :: i, j, g
  real(rt)         :: gE(gPr_l1:gPr_h1,gPr_l2:gPr_h2)
  real(rt)         :: lam, gE1, gE2, r, f, gamr
  integer :: im, ip, jm, jp
  real(rt)         :: xm, xp, ym, yp

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

  if (gPr_l2-1 .ge. vlo(2)) then
     jm = 1
     ym = 2.e0_rt
  else
     jm = 0
     ym = 1.e0_rt
  end if

  if (gPr_h2+1 .le. vhi(2)) then
     jp = 1
     yp = 2.e0_rt
  else
     jp = 0
     yp = 1.e0_rt
  end if

  gPr = 0.0e0_rt

  do g = 0, ngroups-1

     do j = gPr_l2+1, gPr_h2-1
        do i = gPr_l1+1, gPr_h1-1
           gE1 = (Er(i+1,j,g) - Er(i-1,j,g)) / (2.e0_rt*dx(1))
           gE2 = (Er(i,j+1,g) - Er(i,j-1,g)) / (2.e0_rt*dx(2))
           gE(i,j) = sqrt(gE1**2 + gE2**2)
        end do
     end do

     ! lo-x lo-y corner
     i = gPr_l1
     j = gPr_l2
     gE1 = (Er(i+1,j  ,g) - Er(i-im,j   ,g)) / (xm*dx(1))
     gE2 = (Er(i  ,j+1,g) - Er(i   ,j-jm,g)) / (ym*dx(2))
     gE(i,j) = sqrt(gE1**2 + gE2**2)

     ! med-x lo-y side
     j = gPr_l2
     do i = gPr_l1+1, gPr_h1-1
        gE1 = (Er(i+1,j  ,g) - Er(i-1,j   ,g)) / (2.e0_rt*dx(1))
        gE2 = (Er(i  ,j+1,g) - Er(i  ,j-jm,g)) / (  ym*dx(2))
        gE(i,j) = sqrt(gE1**2 + gE2**2)
     end do

     ! hi-x lo-y corner
     i = gPr_h1
     j = gPr_l2
     gE1 = (Er(i+ip,j  ,g) - Er(i-1,j   ,g)) / (xp*dx(1))
     gE2 = (Er(i   ,j+1,g) - Er(i  ,j-jm,g)) / (ym*dx(2))
     gE(i,j) = sqrt(gE1**2 + gE2**2)

     ! lo-x med-y side
     i = gPr_l1
     do j = gPr_l2+1, gPr_h2-1
        gE1 = (Er(i+1,j  ,g) - Er(i-im,j  ,g)) / (  xm*dx(1))
        gE2 = (Er(i  ,j+1,g) - Er(i   ,j-1,g)) / (2.e0_rt*dx(2))
        gE(i,j) = sqrt(gE1**2 + gE2**2)
     end do

     ! hi-x med-y side
     i = gPr_h1
     do j = gPr_l2+1, gPr_h2-1
        gE1 = (Er(i+ip,j  ,g) - Er(i-1,j  ,g)) / (  xp*dx(1))
        gE2 = (Er(i   ,j+1,g) - Er(i  ,j-1,g)) / (2.e0_rt*dx(2))
        gE(i,j) = sqrt(gE1**2 + gE2**2)
     end do

     ! lo-x hi-y corner
     i = gPr_l1
     j = gPr_h2
     gE1 = (Er(i+1,j   ,g) - Er(i-im,j  ,g)) / (xm*dx(1))
     gE2 = (Er(i  ,j+jp,g) - Er(i   ,j-1,g)) / (yp*dx(2))
     gE(i,j) = sqrt(gE1**2 + gE2**2)

     ! med-x hi-y side
     j = gPr_h2
     do i = gPr_l1+1, gPr_h1-1
        gE1 = (Er(i+1,j   ,g) - Er(i-1,j,g)) / (2.e0_rt*dx(1))
        gE2 = (Er(i  ,j+jp,g) - Er(i,j-1,g)) / (  yp*dx(2))
        gE(i,j) = sqrt(gE1**2 + gE2**2)
     end do

     ! hi-x hi-y corner
     i = gPr_h1
     j = gPr_h2
     gE1 = (Er(i+ip,j   ,g) - Er(i-1,j  ,g)) / (xp*dx(1))
     gE2 = (Er(i   ,j+jp,g) - Er(i  ,j-1,g)) / (yp*dx(2))
     gE(i,j) = sqrt(gE1**2 + gE2**2)

     do j = gPr_l2, gPr_h2
        do i = gPr_l1, gPr_h1
           r = gE(i,j) / (kap(i,j,g) * max(Er(i,j,g), 1.e-50_rt))
           lam = FLDlambda(r, limiter)
           if (comoving .eq. 1) then
              f = Edd_factor(lam)
              gamr = (3.e0_rt-f)/2.e0_rt
           else
              gamr = lam + 1.e0_rt
           end if
           gPr(i,j) = gPr(i,j) + lam * gamr * Er(i,j,g)
        end do
     end do

  end do

end subroutine ca_est_gpr2
