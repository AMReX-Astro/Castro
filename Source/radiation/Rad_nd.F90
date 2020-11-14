module rad_nd_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), parameter, private :: tiny = 1.e-50_rt
  real(rt), parameter, private :: BIGKR = 1.e25_rt

contains

  subroutine ca_compute_dcoefs(lo, hi, &
                               d, d_lo, d_hi, &
                               lam, l_lo, l_hi, &
                               v, v_lo, v_hi, &
                               dcf, f_lo, f_hi, &
                               dx, idir) &
                               bind(C, name="ca_compute_dcoefs")

    use habec_nd_module, only: edge_center_metric

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    integer,  intent(in   ) :: l_lo(3), l_hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(inout) :: d(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    real(rt), intent(in   ) :: lam(l_lo(1):l_hi(1),l_lo(2):l_hi(2),l_lo(3):l_hi(3))
    real(rt), intent(in   ) :: v(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),3)
    real(rt), intent(in   ) :: dcf(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: idir

    integer  :: i, j, k
    real(rt) :: r, s

    !$gpu

    if (idir == 0) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                call edge_center_metric(i, j, k, idir + 1, dx, r, s)

                if (v(i-1,j,k,1) + v(i,j,k,1) .gt. 0.e0_rt) then
                   d(i,j,k) = dcf(i-1,j,k) * v(i-1,j,k,1) * lam(i,j,k)
                else if (v(i-1,j,k,1) + v(i,j,k,1) .lt. 0.e0_rt) then
                   d(i,j,k) = dcf(i,j,k) * v(i,j,k,1) * lam(i,j,k)
                else
                   d(i,j,k) = 0.e0_rt
                end if

                d(i,j,k) = d(i,j,k) * r

             end do
          end do
       end do

    else if (idir == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                call edge_center_metric(i, j, k, idir + 1, dx, r, s)

                if (v(i,j-1,k,2) + v(i,j,k,2) .gt. 0.e0_rt) then
                   d(i,j,k) = dcf(i,j-1,k) * v(i,j-1,k,2) * lam(i,j,k)
                else if (v(i,j-1,k,2) + v(i,j,k,2) .lt. 0.e0_rt) then
                   d(i,j,k) = dcf(i,j,k) * v(i,j,k,2) * lam(i,j,k)
                else
                   d(i,j,k) = 0.e0_rt
                end if

                d(i,j,k) = d(i,j,k) * r

             end do
          end do
       end do

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                call edge_center_metric(i, j, k, idir + 1, dx, r, s)

                if (v(i,j,k-1,3) + v(i,j,k,3) .gt. 0.e0_rt) then
                   d(i,j,k) = dcf(i,j,k-1) * v(i,j,k-1,3) * lam(i,j,k)
                else if (v(i,j,k-1,3) + v(i,j,k,3) .lt. 0.e0_rt) then
                   d(i,j,k) = dcf(i,j,k) * v(i,j,k,3) * lam(i,j,k)
                else
                   d(i,j,k) = 0.e0_rt
                end if

                d(i,j,k) = d(i,j,k) * r

             end do
          end do
       end do

    end if

  end subroutine ca_compute_dcoefs



  subroutine lacoefmgfld(lo, hi, &
                         a, a_lo, a_hi, &
                         kappa, k_lo, k_hi, &
                         dx, dt, c) &
                         bind(C, name="lacoefmgfld")

    use habec_nd_module, only: cell_center_metric

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    real(rt), intent(inout) :: a(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
    real(rt), intent(in   ) :: kappa(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, c

    integer  :: i, j, k
    real(rt) :: r, s

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             call cell_center_metric(i, j, k, dx, r, s)

             a(i,j,k) = c * kappa(i,j,k) + 1.e0_rt / dt
             a(i,j,k) = r * s * a(i,j,k)

          end do
       end do
    end do

  end subroutine lacoefmgfld



  subroutine multrs(lo, hi, &
                    d, d_lo, d_hi, &
                    dx) &
                    bind(C, name="multrs")

    use habec_nd_module, only: cell_center_metric
  
    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    real(rt), intent(inout) :: d(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    real(rt), intent(in   ) :: dx(3)

    integer  :: i, j, k
    real(rt) :: r, s

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             call cell_center_metric(i, j, k, dx, r, s)

             d(i,j,k) = d(i,j,k) * r * s

          end do
       end do
    end do

  end subroutine multrs



  subroutine ca_rhstoer(lo, hi, &
                        rhs, r_lo, r_hi, &
                        dx, dt) &
                        bind(C, name="ca_rhstoer")

    use habec_nd_module, only: cell_center_metric

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    real(rt), intent(inout) :: rhs(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt

    integer  :: i, j, k
    real(rt) :: r, s

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             call cell_center_metric(i, j, k, dx, r, s)

             rhs(i,j,k) = rhs(i,j,k) * dt / r

          end do
       end do
    end do

  end subroutine ca_rhstoer



  subroutine ca_compute_rhs(lo, hi, &
                            rhs, rhs_lo, rhs_hi, &
                            jg, jg_lo, jg_hi, &
                            mugT, mugT_lo, mugT_hi, &
                            cpT, cpT_lo, cpT_hi, &
                            etaT, etaT_lo, etaT_hi, &
                            Er2, Er2_lo, Er2_hi, &
                            re2, re2_lo, re2_hi, &
                            Ers, Ers_lo, Ers_hi, &
                            res, res_lo, res_hi, &
                            dx, dt, igroup, tau) &
                            bind(C, name="ca_compute_rhs")

    use rad_params_module, only: ngroups, clight
    use habec_nd_module, only: cell_center_metric

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3) 
    integer,  intent(in   ) :: rhs_lo(3), rhs_hi(3)
    integer,  intent(in   ) :: jg_lo(3), jg_hi(3)
    integer,  intent(in   ) :: mugT_lo(3), mugT_hi(3)
    integer,  intent(in   ) :: cpT_lo(3), cpT_hi(3)
    integer,  intent(in   ) :: etaT_lo(3), etaT_hi(3)
    integer,  intent(in   ) :: Er2_lo(3), Er2_hi(3)
    integer,  intent(in   ) :: re2_lo(3), re2_hi(3)
    integer,  intent(in   ) :: Ers_lo(3), Ers_hi(3)
    integer,  intent(in   ) :: res_lo(3), res_hi(3)
    real(rt), intent(inout) :: rhs(rhs_lo(1):rhs_hi(1),rhs_lo(2):rhs_hi(2),rhs_lo(3):rhs_hi(3))
    real(rt), intent(in   ) :: jg (jg_lo(1):jg_hi(1),jg_lo(2):jg_hi(2),jg_lo(3):jg_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: mugT(mugT_lo(1):mugT_hi(1),mugT_lo(2):mugT_hi(2),mugT_lo(3):mugT_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: cpT(cpT_lo(1):cpT_hi(1),cpT_lo(2):cpT_hi(2),cpT_lo(3):cpT_hi(3))
    real(rt), intent(in   ) :: etaT(etaT_lo(1):etaT_hi(1),etaT_lo(2):etaT_hi(2),etaT_lo(3):etaT_hi(3))
    real(rt), intent(in   ) :: Er2(Er2_lo(1):Er2_hi(1),Er2_lo(2):Er2_hi(2),Er2_lo(3):Er2_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: re2(re2_lo(1):re2_hi(1),re2_lo(2):re2_hi(2),re2_lo(3):re2_hi(3))
    real(rt), intent(in   ) :: Ers(Ers_lo(1):Ers_hi(1),Ers_lo(2):Ers_hi(2),Ers_lo(3):Ers_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: res(res_lo(1):res_hi(1),res_lo(2):res_hi(2),res_lo(3):res_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, tau
    integer,  intent(in   ), value :: igroup

    integer  :: i, j, k
    real(rt) :: Hg, dt1, r, s

    !$gpu

    dt1 = 1.e0_rt / dt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             Hg = mugT(i,j,k,igroup) * etaT(i,j,k)

             rhs(i,j,k) = clight * (jg(i,j,k,igroup) + Hg * cpT(i,j,k))  &
                          + dt1 * (Er2(i,j,k,igroup) - Hg * (res(i,j,k) - re2(i,j,k)) &
                          + tau * Ers(i,j,k,igroup))

             call cell_center_metric(i, j, k, dx, r, s)

             rhs(i,j,k) = r * rhs(i,j,k)

          end do
       end do
    end do

  end subroutine ca_compute_rhs



  subroutine ca_accel_acoe(lo, hi, &
                           eta1, eta1_lo, eta1_hi, &
                           spc, spc_lo, spc_hi, &
                           kap, kap_lo, kap_hi, &
                           aco, aco_lo, aco_hi, &
                           dt, tau) &
                           bind(C, name='ca_accel_acoe')

    use rad_params_module, only: ngroups, clight

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: eta1_lo(3), eta1_hi(3)
    integer,  intent(in   ) :: spc_lo(3), spc_hi(3)
    integer,  intent(in   ) :: kap_lo(3), kap_hi(3)
    integer,  intent(in   ) :: aco_lo(3), aco_hi(3)
    real(rt), intent(in   ) :: eta1(eta1_lo(1):eta1_hi(1),eta1_lo(2):eta1_hi(2),eta1_lo(3):eta1_hi(3))
    real(rt), intent(in   ) :: spc(spc_lo(1):spc_hi(1),spc_lo(2):spc_hi(2),spc_lo(3):spc_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: kap(kap_lo(1):kap_hi(1),kap_lo(2):kap_hi(2),kap_lo(3):kap_hi(3),0:ngroups-1)
    real(rt), intent(inout) :: aco(aco_lo(1):aco_hi(1),aco_lo(2):aco_hi(2),aco_lo(3):aco_hi(3))
    real(rt), intent(in   ), value :: dt, tau

    integer  :: i, j, k
    real(rt) :: kbar, H1, dt1

    !$gpu

    dt1 = (1.e0_rt + tau) / dt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             kbar = sum(spc(i,j,k,:) * kap(i,j,k,:))
             H1 = eta1(i,j,k)
             aco(i,j,k) = H1 * kbar * clight + dt1

          end do
       end do
    end do

  end subroutine ca_accel_acoe



  subroutine ca_accel_rhs(lo, hi, &
                          Ern, Ern_lo, Ern_hi, &
                          Erl, Erl_lo, Erl_hi, &
                          kap, kap_lo, kap_hi, &
                          etaT, etaT_lo, etaT_hi, &
                          rhs, rhs_lo, rhs_hi, &
                          dt) &
                          bind(C, name='ca_accel_rhs')

    use rad_params_module, only: ngroups, clight

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: Ern_lo(3), Ern_hi(3)
    integer,  intent(in   ) :: Erl_lo(3), Erl_hi(3)
    integer,  intent(in   ) :: kap_lo(3), kap_hi(3)
    integer,  intent(in   ) :: etaT_lo(3), etaT_hi(3)
    integer,  intent(in   ) :: rhs_lo(3), rhs_hi(3)
    real(rt), intent(in   ) :: Ern(Ern_lo(1):Ern_hi(1),Ern_lo(2):Ern_hi(2),Ern_lo(3):Ern_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: Erl(Erl_lo(1):Erl_hi(1),Erl_lo(2):Erl_hi(2),Erl_lo(3):Erl_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: kap(kap_lo(1):kap_hi(1),kap_lo(2):kap_hi(2),kap_lo(3):kap_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: etaT(etaT_lo(1):etaT_hi(1),etaT_lo(2):etaT_hi(2),etaT_lo(3):etaT_hi(3))
    real(rt), intent(inout) :: rhs(rhs_lo(1):rhs_hi(1),rhs_lo(2):rhs_hi(2),rhs_lo(3):rhs_hi(3))
    real(rt), intent(in   ), value :: dt

    integer  :: i, j, k
    real(rt) :: rt_term, H

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rt_term = sum(kap(i,j,k,:) * (Ern(i,j,k,:) - Erl(i,j,k,:)))
             H = etaT(i,j,k)
             rhs(i,j,k) = clight * H * rt_term

          end do
       end do
    end do

  end subroutine ca_accel_rhs



  subroutine ca_accel_spec(lo, hi, &
                           kap, kap_lo, kap_hi, &
                           mugT, mugT_lo, mugT_hi, &
                           spec, spec_lo, spec_hi, &
                           dt, tau) &
                           bind(C, name='ca_accel_spec')

    use rad_params_module, only: ngroups, clight

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: kap_lo(3), kap_hi(3)
    integer,  intent(in   ) :: mugT_lo(3), mugT_hi(3)
    integer,  intent(in   ) :: spec_lo(3), spec_hi(3)
    real(rt), intent(in   ) :: kap(kap_lo(1):kap_hi(1),kap_lo(2):kap_hi(2),kap_lo(3):kap_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: mugT(mugT_lo(1):mugT_hi(1),mugT_lo(2):mugT_hi(2),mugT_lo(3):mugT_hi(3),0:ngroups-1)
    real(rt), intent(inout) :: spec(spec_lo(1):spec_hi(1),spec_lo(2):spec_hi(2),spec_lo(3):spec_hi(3),0:ngroups-1)
    real(rt), intent(in   ), value :: dt, tau

    integer  :: i, j, k
    real(rt) :: cdt1, sumeps
    real(rt), dimension(0:ngroups-1) :: epsilon, kapt

    !$gpu

    cdt1 = 1.e0_rt / (clight * dt)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             kapt = kap(i,j,k,:) + (1.e0_rt + tau) * cdt1
             epsilon = mugT(i,j,k,:) / kapt
             sumeps = sum(epsilon)
             if (sumeps .eq. 0.e0_rt) then
                spec(i,j,k,:) = 0.e0_rt
             else
                spec(i,j,k,:) = epsilon / sumeps
             end if

          end do
       end do
    end do

  end subroutine ca_accel_spec



  subroutine ca_accel_ccoe(lo, hi, &
                           bcgr, b_lo, b_hi, &
                           spec, s_lo, s_hi, &
                           ccoe, c_lo, c_hi, &
                           dx, idim, igroup) &
                           bind(C, name='ca_accel_ccoe')

    use rad_params_module, only: ngroups

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: b_lo(3), b_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    real(rt), intent(in   ) :: bcgr(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    real(rt), intent(in   ) :: spec(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:ngroups-1)
    real(rt), intent(inout) :: ccoe(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),0:1)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: idim, igroup

    integer  :: i, j, k, ioff, joff, koff
    real(rt) :: grad_spec, foo, h1

    !$gpu

    if (idim .eq. 0) then
       ioff = 1
       joff = 0
       koff = 0
       h1 = 1.e0_rt / dx(1)
    else if (idim .eq. 1) then
       ioff = 0
       joff = 1
       koff = 0
       h1 = 1.e0_rt / dx(2)
    else
       ioff = 0
       joff = 0
       koff = 1
       h1 = 1.e0_rt / dx(3)
    end if

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             grad_spec = (spec(i,j,k,igroup) - spec(i-ioff,j-joff,k-koff,igroup)) * h1
             foo = - 0.5e0_rt * bcgr(i,j,k) * grad_spec
             ccoe(i,j,k,0) = ccoe(i,j,k,0) + foo
             ccoe(i,j,k,1) = ccoe(i,j,k,1) + foo
          end do
       end do
    end do

  end subroutine ca_accel_ccoe



  subroutine ca_local_accel(lo, hi, &
                            Ern, Ern_lo, Ern_hi, &
                            Erl, Erl_lo, Erl_hi, &
                            kap, kap_lo, kap_hi, &
                            etaT, etaT_lo, etaT_hi, &
                            mugT, mugT_lo, mugT_hi, &
                            dt, tau) &
                            bind(C, name='ca_local_accel')

    use rad_params_module, only: ngroups, clight

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: Ern_lo(3), Ern_hi(3)
    integer,  intent(in   ) :: Erl_lo(3), Erl_hi(3)
    integer,  intent(in   ) :: kap_lo(3), kap_hi(3)
    integer,  intent(in   ) :: etaT_lo(3), etaT_hi(3)
    integer,  intent(in   ) :: mugT_lo(3), mugT_hi(3)
    real(rt), intent(inout) :: Ern(Ern_lo(1):Ern_hi(1),Ern_lo(2):Ern_hi(2),Ern_lo(3):Ern_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: Erl(Erl_lo(1):Erl_hi(1),Erl_lo(2):Erl_hi(2),Erl_lo(3):Erl_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: kap(kap_lo(1):kap_hi(1),kap_lo(2):kap_hi(2),kap_lo(3):kap_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: etaT(etaT_lo(1):etaT_hi(1),etaT_lo(2):etaT_hi(2),etaT_lo(3):etaT_hi(3))
    real(rt), intent(in   ) :: mugT(mugT_lo(1):mugT_hi(1),mugT_lo(2):mugT_hi(2),mugT_lo(3):mugT_hi(3),0:ngroups-1)
    real(rt), intent(in   ), value :: dt, tau

    integer  :: i, j, k
    real(rt) :: cdt1, rt_term, p
    real(rt), dimension(0:ngroups-1) :: Hg, epsilon, kapt, kk

    !$gpu

    cdt1 = 1.e0_rt / (clight * dt)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rt_term = sum(kap(i,j,k,:) * (Ern(i,j,k,:) - Erl(i,j,k,:)))

             Hg = mugT(i,j,k,:) * etaT(i,j,k)

             kapt = kap(i,j,k,:) + (1.e0_rt + tau) * cdt1
             kk = kap(i,j,k,:) / kapt

             p = 1.e0_rt - sum(Hg * kk)
             epsilon = (Hg * rt_term) / (kapt * p + 1.e-50_rt)

             Ern(i,j,k,:) = Ern(i,j,k,:) + epsilon

          end do
       end do
    end do

  end subroutine ca_local_accel




  subroutine lbcoefna(lo, hi, &
                      bcoef, bco_lo, bco_hi, &
                      bcgrp, bcg_lo, bcg_hi, &
                      spec, s_lo, s_hi, &
                      idim) &
                      bind(C, name="lbcoefna")

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: bco_lo(3), bco_hi(3)
    integer,  intent(in   ) :: bcg_lo(3), bcg_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: bcoef(bco_lo(1):bco_hi(1),bco_lo(2):bco_hi(2),bco_lo(3):bco_hi(3))
    real(rt), intent(in   ) :: bcgrp(bcg_lo(1):bcg_hi(1),bcg_lo(2):bcg_hi(2),bcg_lo(3):bcg_hi(3))
    real(rt), intent(in   ) :: spec(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))
    integer,  intent(in   ), value :: idim

    integer :: i, j, k

    !$gpu

    if (idim == 0) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                bcoef(i,j,k) = bcoef(i,j,k) &
                               + 0.5e0_rt * (spec(i-1,j,k) + spec(i,j,k)) * bcgrp(i,j,k)
             end do
          end do
       end do

    else if (idim == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                bcoef(i,j,k) = bcoef(i,j,k) &
                               + 0.5e0_rt * (spec(i,j-1,k) + spec(i,j,k)) * bcgrp(i,j,k)
             end do
          end do
       end do

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                bcoef(i,j,k) = bcoef(i,j,k) &
                               + 0.5e0_rt * (spec(i,j,k-1) + spec(i,j,k)) * bcgrp(i,j,k)
             end do
          end do
       end do

    end if

  end subroutine lbcoefna



  subroutine ljupna(lo, hi, &
                    jnew, j_lo, j_hi, &
                    spec, s_lo, s_hi, &
                    accel, a_lo, a_hi, &
                    nTotal) &
                    bind(C, name="ljupna")

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: j_lo(3), j_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    real(rt), intent(inout) :: jnew(j_lo(1):j_hi(1),j_lo(2):j_hi(2),j_lo(3):j_hi(3),0:nTotal-1)
    real(rt), intent(in   ) :: spec(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),0:nTotal-1)
    real(rt), intent(in   ) :: accel(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
    integer,  intent(in   ), value :: nTotal

    integer :: i, j, k, n

    !$gpu

    do n = 0, nTotal - 1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                jnew(i,j,k,n) = jnew(i,j,k,n) + spec(i,j,k,n) * accel(i,j,k)
             end do
          end do
       end do
    end do

  end subroutine ljupna



  subroutine ca_check_conv_er(lo, hi, &
                              Ern, Ern_lo, Ern_hi, &
                              Erl, Erl_lo, Erl_hi, &
                              kap, kap_lo, kap_hi, &
                              etTz, etTz_lo, etTz_hi, &
                              temp, temp_lo, temp_hi, &
                              rela, abso, errr, dt) &
                              bind(C, name='ca_check_conv_er')

    use rad_params_module, only: ngroups, clight
    use reduction_module, only: reduce_max

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: Ern_lo(3), Ern_hi(3)
    integer,  intent(in   ) :: Erl_lo(3), Erl_hi(3)
    integer,  intent(in   ) :: kap_lo(3), kap_hi(3)
    integer,  intent(in   ) :: temp_lo(3), temp_hi(3)
    integer,  intent(in   ) :: etTz_lo(3), etTz_hi(3)
    real(rt), intent(in   ) :: Ern(Ern_lo(1):Ern_hi(1),Ern_lo(2):Ern_hi(2),Ern_lo(3):Ern_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: Erl(Erl_lo(1):Erl_hi(1),Erl_lo(2):Erl_hi(2),Erl_lo(3):Erl_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: kap(kap_lo(1):kap_hi(1),kap_lo(2):kap_hi(2),kap_lo(3):kap_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: etTz(etTz_lo(1):etTz_hi(1),etTz_lo(2):etTz_hi(2),etTz_lo(3):etTz_hi(3))
    real(rt), intent(in   ) :: temp(temp_lo(1):temp_hi(1),temp_lo(2):temp_hi(2),temp_lo(3):temp_hi(3))
    real(rt), intent(inout) :: rela, abso, errr
    real(rt), intent(in   ), value :: dt

    integer  :: i, j, k, g
    real(rt) :: chg, tot, cdt, der, kde, err_T, err

    !$gpu

    cdt = clight * dt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             chg = 0.e0_rt
             tot = 0.e0_rt
             kde = 0.e0_rt

             do g = 0, ngroups-1
                der = Ern(i,j,k,g) - Erl(i,j,k,g)
                chg = chg + abs(der)
                tot = tot + abs(Ern(i,j,k,g))
                kde = kde + kap(i,j,k,g) * der
             end do

             call reduce_max(abso, chg)
             call reduce_max(rela, chg / (tot + 1.e-50_rt))

             err_T =  etTz(i,j,k) * kde
             err = abs(err_T / (temp(i,j,k) + 1.e-50_rt))

             call reduce_max(errr, err)

          end do
       end do
    end do

  end subroutine ca_check_conv_er



  subroutine ca_update_matter(lo, hi,  &
                              re_n, re_n_lo, re_n_hi, &
                              Er_n, Er_n_lo, Er_n_hi, &
                              Er_l, Er_l_lo, Er_l_hi, &
                              re_s, re_s_lo, re_s_hi, &
                              re_2, re_2_lo, re_2_hi, &
                              eta1, eta1_lo, eta1_hi, &
                              cpt, cpt_lo, cpt_hi, &
                              kpp, kpp_lo, kpp_hi, &
                              dt, tau) &
                              bind(C, name='ca_update_matter')

    use amrex_fort_module, only: rt => amrex_real
    use rad_params_module, only: ngroups, clight

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: re_n_lo(3), re_n_hi(3)
    integer,  intent(in   ) :: Er_n_lo(3), Er_n_hi(3)
    integer,  intent(in   ) :: Er_l_lo(3), Er_l_hi(3)
    integer,  intent(in   ) :: re_s_lo(3), re_s_hi(3)
    integer,  intent(in   ) :: re_2_lo(3), re_2_hi(3)
    integer,  intent(in   ) :: eta1_lo(3), eta1_hi(3)
    integer,  intent(in   ) :: cpt_lo(3), cpt_hi(3)
    integer,  intent(in   ) :: kpp_lo(3), kpp_hi(3)
    real(rt), intent(inout) :: re_n(re_n_lo(1):re_n_hi(1),re_n_lo(2):re_n_hi(2),re_n_lo(3):re_n_hi(3))
    real(rt), intent(in   ) :: Er_n(Er_n_lo(1):Er_n_hi(1),Er_n_lo(2):Er_n_hi(2),Er_n_lo(3):Er_n_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: Er_l(Er_l_lo(1):Er_l_hi(1),Er_l_lo(2):Er_l_hi(2),Er_l_lo(3):Er_l_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: re_s(re_s_lo(1):re_s_hi(1),re_s_lo(2):re_s_hi(2),re_s_lo(3):re_s_hi(3))
    real(rt), intent(in   ) :: re_2(re_2_lo(1):re_2_hi(1),re_2_lo(2):re_2_hi(2),re_2_lo(3):re_2_hi(3))
    real(rt), intent(in   ) :: eta1(eta1_lo(1):eta1_hi(1),eta1_lo(2):eta1_hi(2),eta1_lo(3):eta1_hi(3))
    real(rt), intent(in   ) :: cpt(cpt_lo(1):cpt_hi(1),cpt_lo(2):cpt_hi(2),cpt_lo(3):cpt_hi(3))
    real(rt), intent(in   ) :: kpp(kpp_lo(1):kpp_hi(1),kpp_lo(2):kpp_hi(2),kpp_lo(3):kpp_hi(3),0:ngroups-1)
    real(rt), intent(in   ), value :: dt, tau

    integer  :: i, j, k
    real(rt) :: cdt, H1, dkEE, chg

    !$gpu

    cdt = clight * dt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             H1 = eta1(i,j,k)

             dkEE = sum(kpp(i,j,k,:) * (Er_n(i,j,k,:) - Er_l(i,j,k,:)))

             chg = cdt * dkEE + H1 * ((re_2(i,j,k) - re_s(i,j,k)) + cdt * cpt(i,j,k))

             re_n(i,j,k) = re_s(i,j,k) + chg

             re_n(i,j,k) = (re_n(i,j,k) + tau*re_s(i,j,k)) / (1.e0_rt + tau)

             ! temperature will be updated after exiting this subroutine

          end do
       end do
    end do

  end subroutine ca_update_matter



  subroutine flxlim(lo, hi, &
                    lambda, l_lo, l_hi, &
                    limiter) &
                    bind(C, name="flxlim")

    use amrex_fort_module, only: rt => amrex_real
    use rad_util_module, only: FLDlambda ! function

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: l_lo(3), l_hi(3)
    real(rt), intent(inout) :: lambda(l_lo(1):l_hi(1),l_lo(2):l_hi(2),l_lo(3):l_hi(3))
    integer,  intent(in   ), value :: limiter

    integer :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             lambda(i,j,k) = FLDlambda(lambda(i,j,k), limiter)

          end do
       end do
    end do

  end subroutine flxlim



  subroutine scgrd(lo, hi, &
                   r, r_lo, r_hi, &
                   idir, dx, &
                   kappar, k_lo, k_hi, &
                   er, e_lo, e_hi, &
                   include_cross_terms) &
                   bind(C, name="scgrd")

    use amrex_constants_module, only: FOURTH, HALF
    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    integer,  intent(in   ) :: e_lo(3), e_hi(3)
    real(rt), intent(inout) :: r(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: kappar(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3))
    real(rt), intent(in   ) :: er(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: idir, include_cross_terms

    integer  :: i, j, k, d
    real(rt) :: kap
    real(rt) :: rg, dal, dar, dbl, dbr
    real(rt) :: dxInv(3)

    !$gpu

    do d = 1, dim
       dxInv(d) = 1.e0_rt / dx(d)
    end do
    do d = dim+1, 3
       dxInv(d) = 0.e0_rt
    end do

    if (idir == 0) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (include_cross_terms == 1) then

#if (BL_SPACEDIM >= 2)
                   dal = er(i-1,j+1,k) - er(i-1,j-1,k)
                   dar = er(i  ,j+1,k) - er(i  ,j-1,k)

                   if      (er(i-1,j-1,k) == -1.e0_rt) then
                      dal = 2.e0_rt * (er(i-1,j+1,k) - er(i-1,j  ,k))
                   else if (er(i-1,j+1,k) == -1.e0_rt) then
                      dal = 2.e0_rt * (er(i-1,j  ,k) - er(i-1,j-1,k))
                   end if

                   if      (er(i  ,j-1,k) == -1.e0_rt) then
                      dar = 2.e0_rt * (er(i  ,j+1,k) - er(i  ,j  ,k))
                   else if (er(i  ,j+1,k) == -1.e0_rt) then
                      dar = 2.e0_rt * (er(i  ,j  ,k) - er(i  ,j-1,k))
                   end if
#else
                   dal = 0.e0_rt
                   dar = 0.e0_rt
#endif

#if (BL_SPACEDIM == 3)
                   dbl = er(i-1,j,k+1) - er(i-1,j,k-1)
                   dbr = er(i  ,j,k+1) - er(i  ,j,k-1)

                   if      (er(i-1,j,k-1) == -1.e0_rt) then
                      dbl = 2.e0_rt * (er(i-1,j,k+1) - er(i-1,j,k  ))
                   else if (er(i-1,j,k+1) == -1.e0_rt) then
                      dbl = 2.e0_rt * (er(i-1,j,k  ) - er(i-1,j,k-1))
                   end if

                   if      (er(i  ,j,k-1) == -1.e0_rt) then
                      dbr = 2.e0_rt * (er(i  ,j,k+1) - er(i  ,j,k  ))
                   else if (er(i  ,j,k+1) == -1.e0_rt) then
                      dbr = 2.e0_rt * (er(i  ,j,k  ) - er(i  ,j,k-1))
                   end if
#else
                   dbl = 0.e0_rt
                   dbr = 0.e0_rt
#endif

                else

                   dal = 0.e0_rt
                   dar = 0.e0_rt
                   dbl = 0.e0_rt
                   dbr = 0.e0_rt

                end if
                   

                if (er(i-1,j,k) == -1.e0_rt) then

                   rg = ((er(i+1,j,k) - er(i,j,k)) * dxInv(1))**2 + &
                         (HALF * dar * dxInv(2))**2 + (HALF * dbr * dxInv(3))**2

                else if (er(i,j,k) == -1.e0_rt) then

                   rg = ((er(i-1,j,k) - er(i-2,j,k)) * dxInv(1))**2 + &
                        (HALF * dal * dxInv(2))**2 + (HALF * dbl * dxInv(3))**2

                else

                   rg = ((er(i,j,k) - er(i-1,j,k)) * dxInv(1))**2 + &
                        (FOURTH * (dal + dar) * dxInv(2))**2 + (FOURTH * (dbl + dbr) * dxInv(3))**2

                end if

                kap = kavg(kappar(i-1,j,k), kappar(i,j,k), dx(1), -1)
                r(i,j,k) = sqrt(rg) / (kap * max(er(i-1,j,k), er(i,j,k), tiny))

             end do
          end do
       end do

    else if (idir == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (include_cross_terms == 1) then

                   dal = er(i+1,j-1,k  ) - er(i-1,j-1,k  )
                   dar = er(i+1,j  ,k  ) - er(i-1,j  ,k  )

                   if      (er(i-1,j-1,k  ) == -1.e0_rt) then
                      dal = 2.e0_rt * (er(i+1,j-1,k  ) - er(i  ,j-1,k  ))
                   else if (er(i+1,j-1,k  ) == -1.e0_rt) then
                      dal = 2.e0_rt * (er(i  ,j-1,k  ) - er(i-1,j-1,k  ))
                   end if

                   if      (er(i-1,j  ,k  ) == -1.e0_rt) then
                      dar = 2.e0_rt * (er(i+1,j  ,k  ) - er(i  ,j  ,k  ))
                   else if (er(i+1,j  ,k  ) == -1.e0_rt) then
                      dar = 2.e0_rt * (er(i  ,j  ,k  ) - er(i-1,j  ,k  ))
                   end if

#if (BL_SPACEDIM == 3)
                   dbl = er(i  ,j-1,k+1) - er(i  ,j-1,k-1)
                   dbr = er(i  ,j  ,k+1) - er(i  ,j  ,k-1)

                   if      (er(i  ,j-1,k-1) == -1.e0_rt) then
                      dbl = 2.e0_rt * (er(i  ,j-1,k+1) - er(i  ,j-1,k  ))
                   else if (er(i  ,j-1,k+1) == -1.e0_rt) then
                      dbl = 2.e0_rt * (er(i  ,j-1,k  ) - er(i  ,j-1,k-1))
                   end if

                   if      (er(i  ,j  ,k-1) == -1.e0_rt) then
                      dbr = 2.e0_rt * (er(i  ,j  ,k+1) - er(i  ,j,  k  ))
                   else if (er(i  ,j  ,k+1) == -1.e0_rt) then
                      dbr = 2.e0_rt * (er(i  ,j  ,k  ) - er(i  ,j,  k-1))
                   end if
#else
                   dbl = 0.e0_rt
                   dbr = 0.e0_rt
#endif

                else

                   dal = 0.e0_rt
                   dar = 0.e0_rt
                   dbl = 0.e0_rt
                   dbr = 0.e0_rt

                end if
                   

                if (er(i,j-1,k) == -1.e0_rt) then

                   rg = ((er(i,j+1,k) - er(i,j,k)) * dxInv(2))**2 + &
                        (HALF * dar * dxInv(1))**2 + (HALF * dbr * dxInv(3))**2

                else if (er(i,j,k) == -1.e0_rt) then

                   rg = ((er(i,j-1,k) - er(i,j-2,k)) * dxInv(2))**2 + &
                        (HALF * dal * dxInv(1))**2 + (HALF * dbl * dxInv(3))**2

                else

                   rg = ((er(i,j,k) - er(i,j-1,k)) * dxInv(2))**2 + &
                        (FOURTH * (dal + dar) * dxInv(1))**2 + (FOURTH * (dbl + dbr) * dxInv(3))**2

                end if

                kap = kavg(kappar(i,j-1,k), kappar(i,j,k), dx(2), -1)
                r(i,j,k) = sqrt(rg) / (kap * max(er(i,j-1,k), er(i,j,k), tiny))

             end do
          end do
       end do

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (include_cross_terms == 1) then

                   dal = er(i+1,j  ,k-1) - er(i-1,j  ,k-1)
                   dar = er(i+1,j  ,k  ) - er(i-1,j  ,k  )

                   if      (er(i-1,j  ,k-1) == -1.e0_rt) then
                      dal = 2.e0_rt * (er(i+1,j  ,k-1) - er(i  ,j  ,k-1))
                   else if (er(i+1,j  ,k-1) == -1.e0_rt) then
                      dal = 2.e0_rt * (er(i  ,j  ,k-1) - er(i-1,j  ,k-1))
                   end if

                   if      (er(i-1,j  ,k  ) == -1.e0_rt) then
                      dar = 2.e0_rt * (er(i+1,j  ,k  ) - er(i  ,j  ,k  ))
                   else if (er(i+1,j  ,k  ) == -1.e0_rt) then
                      dar = 2.e0_rt * (er(i  ,j  ,k  ) - er(i-1,j  ,k  ))
                   end if

                   dbl = er(i  ,j+1,k-1) - er(i  ,j-1,k-1)
                   dbr = er(i  ,j+1,k  ) - er(i  ,j-1,k  )

                   if      (er(i  ,j-1,k-1) == -1.e0_rt) then
                      dbl = 2.e0_rt * (er(i  ,j+1,k-1) - er(i  ,j  ,k-1))
                   else if (er(i  ,j+1,k-1) == -1.e0_rt) then
                      dbl = 2.e0_rt * (er(i  ,j  ,k-1) - er(i  ,j-1,k-1))
                   end if

                   if      (er(i  ,j-1,k  ) == -1.e0_rt) then
                      dbr = 2.e0_rt * (er(i  ,j+1,k  ) - er(i  ,j,  k  ))
                   else if (er(i  ,j+1,k  ) == -1.e0_rt) then
                      dbr = 2.e0_rt * (er(i  ,j  ,k  ) - er(i  ,j-1,k  ))
                   end if

                else

                   dal = 0.e0_rt
                   dar = 0.e0_rt
                   dbl = 0.e0_rt
                   dbr = 0.e0_rt

                end if

                if (er(i,j,k-1) == -1.e0_rt) then

                   rg = ((er(i,j,k+1) - er(i,j,k)) * dxInv(3))**2 + &
                        (HALF * dar * dxInv(1))**2 + (HALF * dbr * dxInv(2))**2

                else if (er(i,j,k) == -1.e0_rt) then

                   rg = ((er(i,j,k-1) - er(i,j,k-2)) * dxInv(3))**2 + &
                        (HALF * dal * dxInv(1))**2 + (HALF * dbl * dxInv(2))**2

                else

                   rg = ((er(i,j,k) - er(i,j,k-1)) * dxInv(3))**2 + &
                        (FOURTH * (dal + dar) * dxInv(1))**2 + (FOURTH * (dbl + dbr) * dxInv(2))**2

                end if

                kap = kavg(kappar(i,j,k-1), kappar(i,j,k), dx(3), -1)
                r(i,j,k) = sqrt(rg) / (kap * max(er(i,j,k-1), er(i,j,k), tiny))

             end do
          end do
       end do

    end if

  end subroutine scgrd



  function kavg(a, b, d, opt) result(k)

    use amrex_fort_module, only: rt => amrex_real
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    implicit none

    real(rt), intent(in   ) :: a, b, d
    integer,  intent(in   ) :: opt

    real(rt) :: k

    !$gpu

#ifndef AMREX_USE_GPU
    if (opt > 2) then
       call castro_error("kavg: invalid averaging option")
    end if
#endif

    if (opt == 0) then

       ! arithmetic average, geometrically correct(?) but underestimates surface flux
       k = 0.5e0_rt * (a + b + tiny)

    else if (opt == 1) then

       ! harmonic average, overestimates surface flux
       k = (2.e0_rt * a * b) / (a + b + tiny) + tiny

    else

       ! chooses arithmetic where optically thin, harmonic where optically thick,
       ! surface flux approximation at a thick/thin boundary
       k = min(0.5e0_rt * (a + b + tiny), &
               max((2.e0_rt * a * b) / (a + b + tiny) + tiny, &
                   4.e0_rt / (3.e0_rt * d)))

    end if

  end function kavg



  subroutine lrhs(lo, hi, &
                  rhs, r_lo, r_hi, &
                  temp, t_lo, t_hi, &
                  fkp, fk_lo, fk_hi, &
                  eta, et_lo, et_hi, &
                  etainv, ei_lo, ei_hi, &
                  frhoem, fm_lo, fm_hi, &
                  frhoes, fs_lo, fs_hi, &
                  dfo, d_lo, d_hi, &
                  ero, er_lo, er_hi, &
                  edot, ed_lo, ed_hi, &
                  dt, dx, sigma, c, theta) &
                  bind(C, name="lrhs")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim
    use habec_nd_module, only: cell_center_metric

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: fk_lo(3), fk_hi(3)
    integer,  intent(in   ) :: et_lo(3), et_hi(3)
    integer,  intent(in   ) :: ei_lo(3), ei_hi(3)
    integer,  intent(in   ) :: fm_lo(3), fm_hi(3)
    integer,  intent(in   ) :: fs_lo(3), fs_hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    integer,  intent(in   ) :: er_lo(3), er_hi(3)
    integer,  intent(in   ) :: ed_lo(3), ed_hi(3)
    real(rt), intent(inout) :: rhs(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: fkp(fk_lo(1):fk_hi(1),fk_lo(2):fk_hi(2),fk_lo(3):fk_hi(3))
    real(rt), intent(in   ) :: eta(et_lo(1):et_hi(1),et_lo(2):et_hi(2),et_lo(3):et_hi(3))
    real(rt), intent(in   ) :: etainv(ei_lo(1):ei_hi(1),ei_lo(2):ei_hi(2),ei_lo(3):ei_hi(3))
    real(rt), intent(in   ) :: frhoem(fm_lo(1):fm_hi(1),fm_lo(2):fm_hi(2),fm_lo(3):fm_hi(3))
    real(rt), intent(in   ) :: frhoes(fs_lo(1):fs_hi(1),fs_lo(2):fs_hi(2),fs_lo(3):fs_hi(3))
    real(rt), intent(in   ) :: dfo(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    real(rt), intent(in   ) :: ero(er_lo(1):er_hi(1),er_lo(2):er_hi(2),er_lo(3):er_hi(3))
    real(rt), intent(in   ) :: edot(ed_lo(1):ed_hi(1),ed_lo(2):ed_hi(2),ed_lo(3):ed_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, sigma, c, theta

    integer  :: i, j, k
    real(rt) :: dtm, ek, bs, es, ekt, r, s

    !$gpu

    dtm = 1.e0_rt / dt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ek = fkp(i,j,k) * eta(i,j,k)
             bs = etainv(i,j,k) * 4.e0_rt * sigma * fkp(i,j,k) * temp(i,j,k)**4
             es = eta(i,j,k) * (frhoem(i,j,k) - frhoes(i,j,k))
             ekt = (1.e0_rt - theta) * eta(i,j,k)

             call cell_center_metric(i, j, k, dx, r, s)

             if (dim == 1) then
                s = 1.e0_rt
             end if

             rhs(i,j,k) = (rhs(i,j,k) + r * s * &
                           (bs + dtm * (ero(i,j,k) + es) + &
                            ek * c * edot(i,j,k) - &
                            ekt * dfo(i,j,k))) / (1.e0_rt - ekt)

          end do
       end do
    end do

  end subroutine lrhs



  subroutine lacoef(lo, hi, &
                    a, a_lo, a_hi, &
                    fkp, f_lo, f_hi, &
                    eta, et_lo, et_hi, &
                    etainv, ei_lo, ei_hi, &
                    dx, c, dt, theta) &
                    bind(C, name="lacoef")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim
    use habec_nd_module, only: cell_center_metric

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: et_lo(3), et_hi(3)
    integer,  intent(in   ) :: ei_lo(3), ei_hi(3)
    real(rt), intent(inout) :: a(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
    real(rt), intent(in   ) :: fkp(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(in   ) :: eta(et_lo(1):et_hi(1),et_lo(2):et_hi(2),et_lo(3):et_hi(3))
    real(rt), intent(in   ) :: etainv(ei_lo(1):ei_hi(1),ei_lo(2):ei_hi(2),ei_lo(3):ei_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: c, dt, theta

    integer  :: i, j, k
    real(rt) :: dtm, r, s

    !$gpu

    dtm = 1.e0_rt / dt

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             call cell_center_metric(i, j, k, dx, r, s)

             if (dim == 1) then
                s = 1.e0_rt
             end if

             a(i,j,k) = r * s * &
                  (fkp(i,j,k) * etainv(i,j,k) * c + dtm) / &
                  (1.e0_rt - (1.e0_rt - theta) * eta(i,j,k))

          end do
       end do
    end do

  end subroutine lacoef



  subroutine bclim(lo, hi, &
                   b, b_lo, b_hi, &
                   lambda, l_lo, l_hi, &
                   n, &
                   kappar, k_lo, k_hi, &
                   c, dx) bind(C, name="bclim")

    use amrex_fort_module, only: rt => amrex_real
    use prob_params_module, only: dim
    use habec_nd_module, only: edge_center_metric

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: b_lo(3), b_hi(3)
    integer,  intent(in   ) :: l_lo(3), l_hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    real(rt), intent(inout) :: b(b_lo(1):b_hi(1),b_lo(2):b_hi(2),b_lo(3):b_hi(3))
    real(rt), intent(in   ) :: lambda(l_lo(1):l_hi(1),l_lo(2):l_hi(2),l_lo(3):l_hi(3))
    real(rt), intent(in   ) :: kappar(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: n
    real(rt), intent(in   ), value :: c

    integer  :: i, j, k
    real(rt) :: kap, r, s

    !$gpu

    if (n == 0) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                call edge_center_metric(i, j, k, n+1, dx, r, s)

                if (dim == 1) then
                   s = 1.e0_rt
                end if

                kap = kavg(kappar(i-1,j,k), kappar(i,j,k), dx(1), -1)
                b(i,j,k) = r * s * c * lambda(i,j,k) / kap

             end do
          end do
       end do

    else if (n == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                call edge_center_metric(i, j, k, n+1, dx, r, s)

                if (dim == 1) then
                   s = 1.e0_rt
                end if

                kap = kavg(kappar(i,j-1,k), kappar(i,j,k), dx(2), -1)
                b(i,j,k) = r * s * c * lambda(i,j,k) / kap

             end do
          end do
       end do

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                call edge_center_metric(i, j, k, n+1, dx, r, s)

                if (dim == 1) then
                   s = 1.e0_rt
                end if
                
                kap = kavg(kappar(i,j,k-1), kappar(i,j,k), dx(3), -1)
                b(i,j,k) = r * s * c * lambda(i,j,k) / kap

             end do
          end do
       end do

    end if

  end subroutine bclim



  subroutine ceup(lo, hi, &
                  relres, absres, &
                  frhoes, fs_lo, fs_hi, &
                  frhoem, fm_lo, fm_hi, &
                  eta, et_lo, et_hi, &
                  etainv, ei_lo, ei_hi, &
                  dfo, do_lo, do_hi, &
                  dfn, dn_lo, dn_hi, &
                  exch, ex_lo, ex_hi, &
                  dt, theta) &
                  bind(C, name="ceup")

    use amrex_fort_module, only: rt => amrex_real
    use reduction_module, only: reduce_max

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: fs_lo(3), fs_hi(3)
    integer,  intent(in   ) :: fm_lo(3), fm_hi(3)
    integer,  intent(in   ) :: et_lo(3), et_hi(3)
    integer,  intent(in   ) :: ei_lo(3), ei_hi(3)
    integer,  intent(in   ) :: do_lo(3), do_hi(3)
    integer,  intent(in   ) :: dn_lo(3), dn_hi(3)
    integer,  intent(in   ) :: ex_lo(3), ex_hi(3)
    real(rt), intent(inout) :: frhoes(fs_lo(1):fs_hi(1),fs_lo(2):fs_hi(2),fs_lo(3):fs_hi(3))
    real(rt), intent(in   ) :: frhoem(fm_lo(1):fm_hi(1),fm_lo(2):fm_hi(2),fm_lo(3):fm_hi(3))
    real(rt), intent(in   ) :: eta(et_lo(1):et_hi(1),et_lo(2):et_hi(2),et_lo(3):et_hi(3))
    real(rt), intent(in   ) :: etainv(ei_lo(1):ei_hi(1),ei_lo(2):ei_hi(2),ei_lo(3):ei_hi(3))
    real(rt), intent(in   ) :: dfo(do_lo(1):do_hi(1),do_lo(2):do_hi(2),do_lo(3):do_hi(3))
    real(rt), intent(in   ) :: dfn(dn_lo(1):dn_hi(1),dn_lo(2):dn_hi(2),dn_lo(3):dn_hi(3))
    real(rt), intent(in   ) :: exch(ex_lo(1):ex_hi(1),ex_lo(2):ex_hi(2),ex_lo(3):ex_hi(3))
    real(rt), intent(inout) :: relres, absres
    real(rt), intent(in   ), value :: dt, theta

    integer  :: i, j, k
    real(rt) :: tmp, chg, tot

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             chg = 0.e0_rt
             tot = 0.e0_rt

             tmp = eta(i,j,k) * frhoes(i,j,k) + &
                   etainv(i,j,k) * &
                   (frhoem(i,j,k) - &
                    dt * ((1.e0_rt - theta) * &
                    (dfo(i,j,k) - dfn(i,j,k)) + &
                    exch(i,j,k)))

             chg = abs(tmp - frhoes(i,j,k))
             tot = abs(frhoes(i,j,k))

             frhoes(i,j,k) = tmp

             call reduce_max(absres, chg)
             call reduce_max(relres, chg / (tot + tiny))

          end do
       end do
    end do

  end subroutine ceup



 
  subroutine cetot(lo, hi, &
                   state, s_lo, s_hi, &
                   frhoe, f_lo, f_hi) &
                   bind(C, name="cetot")

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, UEINT, UEDEN

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: frhoe(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))

    integer  :: i, j, k
    real(rt) :: kin

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             kin = state(i,j,k,UEDEN) - state(i,j,k,UEINT)
             state(i,j,k,UEINT) = frhoe(i,j,k)
             state(i,j,k,UEDEN) = frhoe(i,j,k) + kin
          end do
       end do
    end do

  end subroutine cetot



  subroutine ca_compute_rosseland(lo, hi, &
                                  kpr, k_lo, k_hi, &
                                  state, s_lo, s_hi, &
                                  first_group, last_group, num_groups) &
                                  bind(C, name="ca_compute_rosseland")

    use rad_params_module, only: nugroup
    use opacity_table_module, only: get_opacities
    use network, only: naux
    use meth_params_module, only: NVAR, URHO, UTEMP, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: kpr(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3),0:num_groups-1)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer,  intent(in   ), value :: first_group, last_group, num_groups

    integer  :: i, j, k, g
    real(rt) :: kp, kr, nu, rho, temp, Ye
    logical, parameter :: comp_kp = .false. 
    logical, parameter :: comp_kr = .true.

    !$gpu

    ! Note: the group index here will always correspond to the actual frequency
    ! group (so that the access to nugroup makes sense), while the indexing of
    ! the opacity array does not correspond to that. It will always be true when
    ! calling this that num_groups == last_group - first_group + 1.

    do g = first_group, last_group

       nu = nugroup(g)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                rho = state(i,j,k,URHO)
                temp = state(i,j,k,UTEMP)
                if (naux > 0) then
                   Ye = state(i,j,k,UFX)
                else
                   Ye = 0.e0_rt
                end if

                call get_opacities(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr)

                kpr(i,j,k,g-first_group) = kr

             end do
          end do
       end do
    end do

  end subroutine ca_compute_rosseland



  subroutine ca_compute_planck(lo, hi, &
                               kpp, k_lo, k_hi, &
                               state, s_lo, s_hi, &
                               first_group, last_group, num_groups, &
                               temp_offset) &
                               bind(C, name="ca_compute_planck")

    use rad_params_module, only: nugroup
    use opacity_table_module, only: get_opacities
    use network, only: naux
    use meth_params_module, only: NVAR, URHO, UTEMP, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: kpp(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3),0:num_groups-1)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer,  intent(in   ), value :: first_group, last_group, num_groups
    real(rt), intent(in   ), value :: temp_offset

    integer  :: i, j, k, g
    real(rt) :: kp, kr, nu, rho, temp, Ye
    logical, parameter :: comp_kp = .true. 
    logical, parameter :: comp_kr = .false.

    !$gpu

    ! Note: the group index here will always correspond to the actual frequency
    ! group (so that the access to nugroup makes sense), while the indexing of
    ! the opacity array does not correspond to that. It will always be true when
    ! calling this that num_groups == last_group - first_group + 1.

    do g = first_group, last_group

       nu = nugroup(g)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                rho = state(i,j,k,URHO)
                temp = state(i,j,k,UTEMP) + temp_offset
                if (naux > 0) then
                   Ye = state(i,j,k,UFX)
                else
                   Ye = 0.e0_rt
                end if

                call get_opacities(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr)

                kpp(i,j,k,g-first_group) = kp

             end do
          end do
       end do

    end do

  end subroutine ca_compute_planck



  subroutine ca_compute_scattering(lo, hi, &
                                   kps, kps_lo, kps_hi, &
                                   sta, sta_lo, sta_hi) &
                                   bind(C, name='ca_compute_scattering')

    use rad_params_module, only: nugroup
    use opacity_table_module, only: get_opacities
    use network, only: naux
    use meth_params_module, only: NVAR, URHO, UTEMP, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: kps_lo(3), kps_hi(3)
    integer,  intent(in   ) :: sta_lo(3), sta_hi(3)
    real(rt), intent(inout) :: kps(kps_lo(1):kps_hi(1),kps_lo(2):kps_hi(2),kps_lo(3):kps_hi(3))
    real(rt), intent(in   ) :: sta(sta_lo(1):sta_hi(2),sta_lo(2):sta_hi(2),sta_lo(3):sta_hi(3),NVAR)

    integer  :: i, j, k
    real(rt) :: kp, kr, nu, rho, temp, Ye
    logical, parameter :: comp_kp = .true.
    logical, parameter :: comp_kr = .true.

    !$gpu

    ! scattering is assumed to be independent of nu.

    nu = nugroup(0)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rho = sta(i,j,k,URHO)
             temp = sta(i,j,k,UTEMP)
             if (naux > 0) then
                Ye = sta(i,j,k,UFX)
             else
                Ye = 0.e0_rt
             end if

             call get_opacities(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr)

             kps(i,j,k) = max(kr - kp, 0.e0_rt)

          end do
       end do
    end do

  end subroutine ca_compute_scattering



  subroutine ca_opacs(lo, hi, &
                      Snew, s_lo, s_hi, &
                      T,    t_lo, t_hi, &
                      Ts,   ts_lo, ts_hi, &
                      kpp,  kpp_lo, kpp_hi, &
                      kpr,  kpr_lo, kpr_hi, &
                      dkdT, dkdT_lo, dkdT_hi, &
                      use_dkdT, validStar, lag_opac) &
                      bind(C, name='ca_opacs')

    use rad_params_module, only: ngroups, nugroup
    use opacity_table_module, only: get_opacities
    use network, only: naux
    use meth_params_module, only: NVAR, URHO, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: ts_lo(3), ts_hi(3)
    integer,  intent(in   ) :: kpp_lo(3),  kpp_hi(3)
    integer,  intent(in   ) :: kpr_lo(3),  kpr_hi(3)
    integer,  intent(in   ) :: dkdT_lo(3), dkdT_hi(3)
    real(rt), intent(in   ) :: Snew(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: T(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: Ts(ts_lo(1):ts_hi(1),ts_lo(2):ts_hi(2),ts_lo(3):ts_hi(3))
    real(rt), intent(inout) :: kpp ( kpp_lo(1): kpp_hi(1), kpp_lo(2): kpp_hi(2), kpp_lo(3): kpp_hi(3),0:ngroups-1)
    real(rt), intent(inout) :: kpr ( kpr_lo(1): kpr_hi(1), kpr_lo(2): kpr_hi(2), kpr_lo(3): kpr_hi(3),0:ngroups-1)
    real(rt), intent(inout) :: dkdT(dkdT_lo(1):dkdT_hi(1),dkdT_lo(2):dkdT_hi(2),dkdT_lo(3):dkdT_hi(3),0:ngroups-1)
    integer,  intent(in   ), value :: use_dkdT, validStar, lag_opac

    integer  :: i, j, k, g
    real(rt) :: kp, kr, nu, rho, temp, Ye
    real(rt) :: kp1, kr1
    real(rt) :: kp2, kr2
    real(rt) :: dT
    logical  :: comp_kp, comp_kr
    real(rt), parameter :: fac = 0.5e0_rt, minfrac = 1.e-8_rt

    !$gpu

    if (lag_opac .eq. 1) then
       dkdT(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.e0_rt
       return
    end if

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rho = Snew(i,j,k,URHO)
             temp = T(i,j,k)
             if (naux > 0) then
                Ye = Snew(i,j,k,UFX)
             else
                Ye = 0.e0_rt
             end if

             if (validStar > 0) then
                dT = fac * abs(Ts(i,j,k) - T(i,j,k))
                dT = max(dT, minfrac * T(i,j,k))
             else
                dT = T(i,j,k) * 1.e-3_rt + 1.e-50_rt
             end if

             do g = 0, ngroups - 1

                nu = nugroup(g)

                comp_kp = .true.
                comp_kr = .true.

                call get_opacities(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr)
                kpp(i,j,k,g) = kp
                kpr(i,j,k,g) = kr

                if (use_dkdT .eq. 0) then        

                   dkdT(i,j,k,g) = 0.e0_rt

                else

                   comp_kp = .true.
                   comp_kr = .false.

                   call get_opacities(kp1, kr1, rho, temp-dT, Ye, nu, comp_kp, comp_kr)
                   call get_opacities(kp2, kr2, rho, temp+dT, Ye, nu, comp_kp, comp_kr)

                   dkdT(i,j,k,g) = (kp2 - kp1) / (2.e0_rt * dT)

                end if

             end do
          end do
       end do
    end do

  end subroutine ca_opacs



  subroutine ca_compute_emissivity(lo, hi, &
                                   jg, jg_lo, jg_hi, &
                                   djdT, djdT_lo, djdT_hi, &
                                   T, T_lo, T_hi, &
                                   kap, kap_lo, kap_hi, &
                                   dkdT, dkdT_lo, dkdT_hi) &
                                   bind(C, name='ca_compute_emissivity')

    use amrex_fort_module, only: rt => amrex_real
    use rad_params_module, only: ngroups, nugroup, dnugroup, xnu, arad
    use blackbody_module, only: BdBdTIndefInteg
    use emissivity_override_module, only: emissivity_override

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3) 
    integer,  intent(in   ) :: jg_lo(3), jg_hi(3)
    integer,  intent(in   ) :: djdT_lo(3), djdT_hi(3)
    integer,  intent(in   ) :: T_lo(3), T_hi(3)
    integer,  intent(in   ) :: kap_lo(3), kap_hi(3)
    integer,  intent(in   ) :: dkdT_lo(3), dkdT_hi(3)
    real(rt), intent(inout) :: jg(jg_lo(1):jg_hi(1),jg_lo(2):jg_hi(2),jg_lo(3):jg_hi(3),0:ngroups-1)
    real(rt), intent(inout) :: djdT(djdT_lo(1):djdT_hi(1),djdT_lo(2):djdT_hi(2),djdT_lo(3):djdT_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: T(T_lo(1):T_hi(1),T_lo(2):T_hi(2),T_lo(3):T_hi(3))
    real(rt), intent(in   ) :: kap(kap_lo(1):kap_hi(1),kap_lo(2):kap_hi(2),kap_lo(3):kap_hi(3),0:ngroups-1)
    real(rt), intent(in   ) :: dkdT(dkdT_lo(1):dkdT_hi(1),dkdT_lo(2):dkdT_hi(2),dkdT_lo(3):dkdT_hi(3),0:ngroups-1)

    integer  :: i, j, k, g
    real(rt) :: dBdT, Bg
    real(rt) :: Teff
    real(rt) :: B0, B1, dBdT0, dBdT1
    real(rt) :: xnup

    !$gpu

    ! Integrate the Planck distribution upward from zero frequency.
    ! This handles both the single-group and multi-group cases.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             Teff = max(T(i,j,k), 1.e-50_rt)
             call BdBdTIndefInteg(Teff, 0.0_rt, B1, dBdT1)

             do g = 0, ngroups-1

                xnup = xnu(g+1)

                ! For the last group, make sure that we complete
                ! the integral up to "infinity".

                if (g == ngroups - 1) then
                   xnup = max(xnup, 1.e25_rt)
                end if

                B0 = B1
                dBdT0 = dBdT1
                call BdBdTIndefInteg(Teff, xnup, B1, dBdT1)
                Bg = B1 - B0
                dBdT = dBdT1 - dBdT0

                jg(i,j,k,g) = Bg * kap(i,j,k,g)
                djdT(i,j,k,g) = dkdT(i,j,k,g) * Bg + dBdT * kap(i,j,k,g)

                ! Allow a problem to override this emissivity.

                call emissivity_override(i, j, k, g, T(i,j,k), kap(i,j,k,g), dkdT(i,j,k,g), jg(i,j,k,g), djdT(i,j,k,g))

             end do

          end do
       end do
    end do

  end subroutine ca_compute_emissivity




  subroutine nfloor(lo, hi, &
                    dest, d_lo, d_hi, &
                    nvar) &
                    bind(C, name="nfloor")

    use amrex_fort_module, only: rt => amrex_real

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    real(rt), intent(inout) :: dest(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),0:nvar-1)
    integer,  intent(in   ), value :: nvar

    integer :: i, j, k, n

    real(rt), parameter :: temp_floor = 1.e-10_rt

    !$gpu

    do n = 0, nvar-1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (dest(i,j,k,n) < temp_floor) then
                   dest(i,j,k,n) = temp_floor
                end if
             end do
          end do
       end do
    end do

  end subroutine nfloor



  subroutine ceta2(lo, hi, &
                   eta, et_lo, et_hi, &
                   etainv, ei_lo, ei_hi, &
                   frho, fr_lo, fr_hi, &
                   temp, t_lo, t_hi, &
                   cv, c_lo, c_hi, &
                   fkp, fk_lo, fk_hi, &
                   er, er_lo, er_hi, &
                   dtemp, dtime, sigma, &
                   c, underr, lagpla) &
                   bind(C, name="ceta2")

    use amrex_fort_module, only: rt => amrex_real

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: et_lo(3), et_hi(3)
    integer,  intent(in   ) :: ei_lo(3), ei_hi(3)
    integer,  intent(in   ) :: fr_lo(3), fr_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    integer,  intent(in   ) :: fk_lo(3), fk_hi(3)
    integer,  intent(in   ) :: er_lo(3), er_hi(3)
    real(rt), intent(inout) :: eta(et_lo(1):et_hi(1),et_lo(2):et_hi(2),et_lo(3):et_hi(3))
    real(rt), intent(inout) :: etainv(ei_lo(1):ei_hi(1),ei_lo(2):ei_hi(2),ei_lo(3):ei_hi(3))
    real(rt), intent(in   ) :: frho(fr_lo(1):fr_hi(1),fr_lo(2):fr_hi(2),fr_lo(3):fr_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: cv(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(rt), intent(in   ) :: fkp(fk_lo(1):fk_hi(1),fk_lo(2):fk_hi(2),fk_lo(3):fk_hi(3))
    real(rt), intent(in   ) :: er(er_lo(1):er_hi(1),er_lo(2):er_hi(2),er_lo(3):er_hi(3))
    real(rt), intent(in   ), value :: dtemp, dtime, sigma, c, underr
    integer,  intent(in   ), value :: lagpla

    real(rt) :: d, frc, fac0, fac1, fac2
    integer  :: i, j, k

    !$gpu

    fac1 = 16.e0_rt * sigma * dtime
    if (lagpla == 0) then
       fac0 = 0.25e0_rt * fac1 / dtemp
       fac2 = dtime * c / dtemp
    endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (lagpla /= 0) then

                ! assume eta and fkp are the same
                d = fac1 * fkp(i,j,k) * temp(i,j,k) ** 3

             else

                d = fac0 * (eta(i,j,k) * (temp(i,j,k) + dtemp) ** 4 - &
                            fkp(i,j,k) * (temp(i,j,k)        ) ** 4) - &
                            fac2 * (eta(i,j,k) - fkp(i,j,k)) * er(i,j,k)
                ! alternate form, sometimes worse, sometimes better:
                !                  d = fac1 * fkp(i,j,k) * temp(i,j,k) ** 3 +
                !     @                fac0 * (eta(i,j,k) - fkp(i,j,k)) * temp(i,j,k) ** 4 -
                !     @                fac2 * (eta(i,j,k) - fkp(i,j,k)) * er(i,j,k)
                ! another alternate form (much worse):
                !                  d = fac1 * fkp(i,j,k) * (temp(i,j,k) + dtemp) ** 3 +
                !     @                fac0 * (eta(i,j,k) - fkp(i,j,k))
                !     @                     * (temp(i,j,k) + dtemp) ** 4 -
                !     @                fac2 * (eta(i,j,k) - fkp(i,j,k)) * er(i,j,k)
             end if

             frc = frho(i,j,k) * cv(i,j,k) + tiny
             eta(i,j,k) = d / (d + frc)
             etainv(i,j,k) = underr * frc / (d + frc)
             eta(i,j,k) = 1.e0_rt - etainv(i,j,k)
             !               eta(i,j,k) = 1.e0_rt - underr * (1.e0_rt - eta(i,j,k))

          end do
       end do
    end do

  end subroutine ceta2



  subroutine ca_compute_c_v(lo, hi, &
                            cv, c_lo, c_hi, &
                            temp, t_lo, t_hi, &
                            state, s_lo, s_hi) &
                            bind(C, name="ca_compute_c_v")

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, UFS, UFX, &
                                  prop_temp_floor
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: cv(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer     :: i, j, k
    real(rt)    :: rhoInv
    type(eos_t) :: eos_state
    real(rt)    :: alpha, teff

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = 1.e0_rt / state(i,j,k,URHO)
             eos_state % rho = state(i,j,k,URHO)
             eos_state % T   = temp(i,j,k)
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos(eos_input_rt, eos_state)

             cv(i,j,k) = eos_state % cv

          end do
       end do
    end do

  end subroutine ca_compute_c_v



  subroutine ca_get_rhoe(lo, hi, &
                         rhoe, r_lo, r_hi, &
                         temp, t_lo, t_hi, &
                         state, s_lo, s_hi) &
                         bind(C, name="ca_get_rhoe")

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, UFS, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: rhoe(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer     :: i, j, k
    real(rt)    :: rhoInv
    type(eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = 1.e0_rt / state(i,j,k,URHO)
             eos_state % rho = state(i,j,k,URHO)
             eos_state % T   =  temp(i,j,k)
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux = state(i,j,k,UFX:UFX+naux -1) * rhoInv

             call eos(eos_input_rt, eos_state)

             rhoe(i,j,k) = eos_state % rho * eos_state % e

          end do
       end do
    end do

  end subroutine ca_get_rhoe



  subroutine ca_compute_temp_given_rhoe(lo, hi, &
                                        temp, t_lo, t_hi, &
                                        state, s_lo, s_hi, &
                                        update_state) &
                                        bind(C, name="ca_compute_temp_given_rhoe")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only: NVAR, URHO, UTEMP, UFS, UFX, small_temp
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! temp contains rhoe as input
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer,  intent(in   ), value :: update_state

    integer      :: i, j, k
    real(rt)     :: rhoInv
    type (eos_t) :: eos_state
    real(rt)     :: ex, alpha, rhoal, teff

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (temp(i,j,k) .le. 0.e0_rt) then

                temp(i,j,k) = small_temp

             else

                rhoInv = 1.e0_rt / state(i,j,k,URHO)
                eos_state % rho = state(i,j,k,URHO)
                eos_state % T   = state(i,j,k,UTEMP)
                eos_state % e   =  temp(i,j,k)*rhoInv 
                eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
                eos_state % aux = state(i,j,k,UFX:UFX+naux -1) * rhoInv

                call eos(eos_input_re, eos_state)

                temp(i,j,k) = eos_state % T

             end if

             if (update_state == 1) then
                state(i,j,k,UTEMP) = temp(i,j,k)
             end if

          end do
       end do
    end do

  end subroutine ca_compute_temp_given_rhoe



  subroutine cfrhoe(lo, hi, &
                    frhoe, f_lo, f_hi, &
                    state, s_lo, s_hi) &
                    bind(C, name='cfrhoe')

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, UEINT

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: frhoe(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ! kin = 0.5e0_rt * (state(i,j,k,XMOM)   ** 2 +
             !                   state(i,j,k,XMOM+1) ** 2 +
             !                   state(i,j,k,XMOM+2) ** 2) /
             !                   state(i,j,k,DEN)
             ! frhoe(i,j,k) = state(i,j,k,EDEN) - kin
             frhoe(i,j,k) = state(i,j,k,UEINT)
          end do
       end do
    end do

  end subroutine cfrhoe

end module rad_nd_module



!! -----------------------------------------------------------
!> @brief This routine is called at problem setup time and is used
!! to initialize values of physical constants used by the
!! radiation package.
!! -----------------------------------------------------------
subroutine ca_initradconstants(p, c, h, k, s, a, m, J_is_used) bind(C, name="ca_initradconstants")

  use fundamental_constants_module, only: c_fcm=>c_light, h_fcm=>hplanck, &
                                          k_fcm=>k_B, s_fcm=>sigma_SB, a_fcm=>n_A, ev2erg_fcm=>ev2erg
  use rad_params_module, only: pi, clight, hplanck
  use rad_params_module, only: kboltz, stefbol, arad, avogadro
  use rad_params_module, only: Hz2MeV, mev2erg, tiny
  use rad_params_module, only: radtoE  !, radtoJ, Etorad, radfluxtoF
  use rad_params_module, only: etafactor
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt) :: p, c, h, k, s, a, m
  integer :: J_is_used

  c = c_fcm
  h = h_fcm
  k = k_fcm
  s = s_fcm
  a = a_fcm
  m = 1.e6_rt * ev2erg_fcm

  allocate(pi)
  allocate(clight)
  allocate(hplanck)
  allocate(kboltz)
  allocate(stefbol)
  allocate(arad)
  allocate(avogadro)
  allocate(Hz2MeV)
  allocate(mev2erg)
  allocate(tiny)

  pi       = p
  clight   = c
  hplanck  = h
  kboltz   = k
  stefbol  = s
  arad     = 4.*stefbol/clight
  avogadro = a
  mev2erg  = m
  Hz2MeV   = h / m
  tiny     = 1.e-50_rt

  allocate(radtoE)
  allocate(etafactor)

  if (J_is_used > 0) then
     radtoE = 4.e0_rt*pi/clight
     !           radtoJ = 1.0e0_rt
     !           Etorad = 1.e0_rt/radtoE
     !           radfluxtoF = 4.e0_rt*pi
     etafactor = 1.e0_rt
  else
     radtoE = 1.0e0_rt
     !           radtoJ = clight/(4.e0_rt*pi)
     !           Etorad = 1.0e0_rt
     !           radfluxtoF = 1.e0_rt
     etafactor = 4.e0_rt*pi/clight
  end if

end subroutine ca_initradconstants

! For single group, let set ngroups to 1.
subroutine ca_initsinglegroup(ngr) bind(C, name="ca_initsinglegroup")

  use rad_params_module, only: ngroups, nugroup, dnugroup, ng0, ng1
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer :: ngr

  ! Local variables
  integer :: i

  ng0 = 0
  ng1 = 0

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))

  do i = 0, ngroups-1
     nugroup(i)  = 1.e0_rt  ! dummy
     dnugroup(i) = 1.e0_rt
  enddo

end subroutine ca_initsinglegroup

!! -----------------------------------------------------------
!> @brief This routine is called at problem setup time and is used
!! to initialize the arrays nugroup and dnugroup in
!! probdata with the neutrino group energies and widths.
!!
!! The widths are used to derive neutrino spectrum for plot files
!! -----------------------------------------------------------
subroutine ca_initgroups(nugr, dnugr, ngr, ngr0, ngr1)

  use rad_params_module, only: ngroups, ng0, ng1, nugroup, dnugroup, &
                               current_group
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt) :: nugr(0:ngr-1), dnugr(0:ngr-1)
  integer :: ngr, ngr0, ngr1

  ! Local variables
  integer :: i

  allocate(current_group, ng0, ng1)

  ng0     = ngr0
  ng1     = ngr1

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))

  do i = 0, ngroups-1
     nugroup(i)  = nugr(i)
     dnugroup(i) = dnugr(i)
  enddo

end subroutine ca_initgroups

subroutine ca_initgroups2(nugr, dnugr, xnugr, ngr)

  use rad_params_module, only: ngroups, nugroup, dnugroup, xnu, dlognu, lognugroup, &
                               current_group, ng0, ng1
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), intent(in) :: nugr(0:ngr-1), dnugr(0:ngr-1), xnugr(0:ngr)
  integer :: ngr

  ! Local variables
  integer   :: i

  allocate(current_group, ng0, ng1)

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))
  allocate(xnu(0:ngroups))
  allocate(dlognu(0:ngroups-1))
  allocate(lognugroup(0:ngroups-1))

  nugroup(:) = nugr(:)
  dnugroup(:) = dnugr(:)
  xnu(:) = xnugr(:)
  lognugroup(:) = log(nugroup)

  dlognu(0:ngroups-1) = log(xnu(1:ngroups)) - log(xnu(0:ngroups-1))

end subroutine ca_initgroups2

subroutine ca_initgroups3(nugr, dnugr, dlognugr, xnugr, ngr, ngr0, ngr1)
  ! used by MGFLDSolver

  use rad_params_module, only: ngroups, ng0, ng1, nugroup, dnugroup, &
                               xnu, dlognu, lognugroup, erg2rhoYe, avogadro, hplanck, &
                               current_group
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), intent(in) :: nugr(0:ngr-1), dnugr(0:ngr-1), dlognugr(0:ngr-1), xnugr(0:ngr)
  integer :: ngr, ngr0, ngr1

  ! Local variables
  integer :: i

  allocate(current_group, ng0, ng1)

  ng0     = ngr0
  ng1     = ngr1

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))
  allocate(xnu(0:ngroups))
  allocate(dlognu(0:ngroups-1))
  allocate(erg2rhoYe(0:ngroups-1))
  allocate(lognugroup( 0:ngroups-1))

  nugroup(:) = nugr(:)
  dnugroup(:) = dnugr(:)
  xnu(:) = xnugr(:)
  dlognu(:) = dlognugr(:)
  lognugroup(:) = log(nugroup)

  erg2rhoYe = 0.e0_rt
  if (ng0 > 0) then
     erg2rhoYe(0:ng0-1) = 1.e0_rt / (avogadro*hplanck*nugroup(0:ng0-1))
     if (ng1 > 0) then
        erg2rhoYe(ng0:ng0+ng1-1) = -1.e0_rt / (avogadro*hplanck*nugroup(ng0:ng0+ng1-1))
     end if
  end if

end subroutine ca_initgroups3

subroutine ca_get_dlognu(dlognu_out) bind(C, name="ca_get_dlognu")

  use amrex_fort_module, only: rt => amrex_real
  use rad_params_module, only: ngroups, dlognu
  implicit none

  real(rt), intent(out) :: dlognu_out(0:ngroups-1)

  dlognu_out(:) = dlognu(:)

end subroutine ca_get_dlognu

subroutine ca_get_nugroup(nugroup_out) bind(C, name="ca_get_nugroup")

  use amrex_fort_module, only: rt => amrex_real
  use rad_params_module, only: ngroups, nugroup
  implicit none

  real(rt), intent(out) :: nugroup_out(0:ngroups-1)

  nugroup_out(:) = nugroup(:)

end subroutine ca_get_nugroup

subroutine ca_get_dnugroup(dnugroup_out) bind(C, name="ca_get_dnugroup")

  use amrex_fort_module, only: rt => amrex_real
  use rad_params_module, only: ngroups, dnugroup
  implicit none

  real(rt), intent(out) :: dnugroup_out(0:ngroups-1)

  dnugroup_out(:) = dnugroup(:)

end subroutine ca_get_dnugroup

!! -----------------------------------------------------------

subroutine ca_inelastic_sct(lo, hi, &
                            uu,uu_lo,uu_hi, &
                            Er,Er_lo,Er_hi, &
                            ks,ks_lo,ks_hi, &
                            dt) bind(C)

  use meth_params_module, only: NVAR, UEDEN, UEINT, UTEMP
  use rad_params_module, only: ngroups, nugroup, dlognu
  use radhydro_nd_module, only: inelastic_scatter
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: uu_lo(3),uu_hi(3)
  integer, intent(in) :: Er_lo(3),Er_hi(3)
  integer, intent(in) :: ks_lo(3),ks_hi(3)
  real(rt), intent(inout) :: uu(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3),NVAR)
  real(rt), intent(inout) :: Er(Er_lo(1):Er_hi(1),Er_lo(2):Er_hi(2),Er_lo(3):Er_hi(3),0:ngroups-1)
  real(rt), intent(in   ) :: ks(ks_lo(1):ks_hi(1),ks_lo(2):ks_hi(2),ks_lo(3):ks_hi(3))
  real(rt), intent(in) :: dt

  integer :: i, j, k
  real(rt) :: Ertotold, Ertmp(0:ngroups-1), dEr
  real(rt) :: Erscale(0:ngroups-1)

  Erscale = nugroup*dlognu

  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           Ertmp = Er(i,j,k,:)
           Ertotold = sum(Ertmp)
           Ertmp = Ertmp / Erscale

           call inelastic_scatter(uu(i,j,k,UTEMP), Ertmp, ks(i,j,k), dt, (/i,j,k/))

           Ertmp = Ertmp * Erscale
           dEr = sum(Ertmp) - Ertotold
           Er(i,j,k,:) = Ertmp
           uu(i,j,k,UEINT) = uu(i,j,k,UEINT) - dEr
           uu(i,j,k,UEDEN) = uu(i,j,k,UEDEN) - dEr
        end do
     end do
  end do

end subroutine ca_inelastic_sct
