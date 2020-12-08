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
