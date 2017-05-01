
! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize values of physical constants used by the
! ::: radiation package.
! ::: -----------------------------------------------------------

subroutine ca_initradconstants(p, c, h, k, s, a, m, J_is_used) bind(C, name="ca_initradconstants")

  use fundamental_constants_module, only : c_fcm=>c_light, h_fcm=>hplanck, &
       k_fcm=>k_B, s_fcm=>sigma_SB, a_fcm=>n_A, ev2erg_fcm=>ev2erg

  use rad_params_module, only: pi, clight, hplanck
  use rad_params_module, only: kboltz, stefbol, arad, avogadro
  use rad_params_module, only: Hz2MeV, mev2erg, tiny
  use rad_params_module, only: radtoE  !, radtoJ, Etorad, radfluxtoF
  use rad_params_module, only: etafactor

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)         p, c, h, k, s, a, m
  integer J_is_used

  c = c_fcm
  h = h_fcm
  k = k_fcm
  s = s_fcm
  a = a_fcm
  m = 1.e6_rt * ev2erg_fcm

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

  use rad_params_module, only : ngroups, nugroup, dnugroup, ng0, ng1
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer ngr

  ! Local variables
  integer   :: i

  ngroups = ngr
  ng0 = 0
  ng1 = 0

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))

  do i = 0, ngroups-1
     nugroup(i)  = 1.e0_rt  ! dummy
     dnugroup(i) = 1.e0_rt
  enddo
end subroutine ca_initsinglegroup

! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize the arrays nugroup and dnugroup in
! ::: probdata with the neutrino group energies and widths.
! :::
! ::: The widths are used to derive neutrino spectrum for plot files
! ::: -----------------------------------------------------------
subroutine ca_initgroups(nugr, dnugr, ngr, ngr0, ngr1)

  use rad_params_module, only: ngroups, ng0, ng1, nugroup, dnugroup

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)         nugr(0:ngr-1), dnugr(0:ngr-1)
  integer ngr, ngr0, ngr1

  ! Local variables
  integer   :: i

  ngroups = ngr
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

  use rad_params_module, only: ngroups, nugroup, dnugroup, xnu, dlognu, lognugroup

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , intent(in) :: nugr(0:ngr-1), dnugr(0:ngr-1), xnugr(0:ngr)
  integer ngr

  ! Local variables
  integer   :: i

  ngroups = ngr

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

  use rad_params_module, only: ngroups, ng0, ng1, nnuspec, nradspec, nugroup, dnugroup, &
       xnu, dlognu, lognugroup, erg2rhoYe, avogadro, hplanck

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , intent(in) :: nugr(0:ngr-1), dnugr(0:ngr-1), dlognugr(0:ngr-1), xnugr(0:ngr+2)
  integer ngr, ngr0, ngr1

  ! Local variables
  integer   :: i

  ngroups = ngr
  ng0     = ngr0
  ng1     = ngr1

  if (ng0 > 0) then
     if (ng1 .eq. 0) then
        nnuspec = 1  ! one neutrino species
     else if (ngroups .eq. ng0+ng1) then
        nnuspec = 2  ! two neutrino species
     else
        nnuspec = 3  ! three neutrino species
     end if
  else
     nnuspec = 0
  end if

  nradspec = max(nnuspec, 1)

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))
  allocate(xnu(0:ngroups+2))
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

! ::: -----------------------------------------------------------

subroutine ca_setgroup(igroup)

  use rad_params_module, only: current_group

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer igroup

  current_group = igroup

end subroutine ca_setgroup

! ::: -----------------------------------------------------------

subroutine ca_inelastic_sct (lo, hi, &
     uu,uu_l1,uu_l2,uu_l3,uu_h1,uu_h2,uu_h3, &
     Er,Er_l1,Er_l2,Er_l3,Er_h1,Er_h2,Er_h3, &
     ks,ks_l1,ks_l2,ks_l3,ks_h1,ks_h2,ks_h3, &
     dt) bind(C)

  use meth_params_module, only : NVAR, UEDEN, UEINT, UTEMP
  use rad_params_module, only : ngroups, nugroup, dlognu
  use radhydro_nd_module, only: inelastic_scatter

  use amrex_fort_module, only : rt => amrex_real
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: uu_l1,uu_l2,uu_l3,uu_h1,uu_h2,uu_h3
  integer, intent(in) :: Er_l1,Er_l2,Er_l3,Er_h1,Er_h2,Er_h3
  integer, intent(in) :: ks_l1,ks_l2,ks_l3,ks_h1,ks_h2,ks_h3
  real(rt)        , intent(inout) :: uu(uu_l1:uu_h1,uu_l2:uu_h2,uu_l3:uu_h3,NVAR)
  real(rt)        , intent(inout) :: Er(Er_l1:Er_h1,Er_l2:Er_h2,Er_l3:Er_h3,0:ngroups-1)
  real(rt)        , intent(in   ) :: ks(ks_l1:ks_h1,ks_l2:ks_h2,ks_l3:ks_h3)
  real(rt)        , intent(in) :: dt

  integer :: i, j, k
  real(rt)         :: Ertotold, Ertmp(0:ngroups-1), dEr
  real(rt)         :: Erscale(0:ngroups-1)

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

! ::: -----------------------------------------------------------

subroutine ca_compute_scattering(lo, hi, &
     kps,kps_l1,kps_l2,kps_l3,kps_h1,kps_h2,kps_h3, &
     sta,sta_l1,sta_l2,sta_l3,sta_h1,sta_h2,sta_h3)

  use rad_params_module, only : ngroups, nugroup
  use opacity_table_module, only : get_opacities
  use network, only : naux
  use meth_params_module, only : NVAR, URHO, UTEMP, UFX

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: kps_l1,kps_l2,kps_l3,kps_h1,kps_h2,kps_h3
  integer, intent(in) :: sta_l1,sta_l2,sta_l3,sta_h1,sta_h2,sta_h3
  real(rt)        , intent(inout) :: kps(kps_l1:kps_h1,kps_l2:kps_h2,kps_l3:kps_h3)
  real(rt)        , intent(in   ) :: sta(sta_l1:sta_h1,sta_l2:sta_h2,sta_l3:sta_h3,NVAR)

  integer :: i, j, k
  real(rt)         :: kp, kr, nu, rho, temp, Ye
  logical, parameter :: comp_kp = .true.
  logical, parameter :: comp_kr = .true.

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

subroutine ca_compute_scattering_2(lo, hi, &
     kps,kps_l1,kps_l2,kps_l3,kps_h1,kps_h2,kps_h3, &
     sta,sta_l1,sta_l2,sta_l3,sta_h1,sta_h2,sta_h3, &
     k0_p, m_p, n_p, &
     k0_r, m_r, n_r, &
     Tfloor, kfloor)

  use rad_params_module, only : ngroups, nugroup
  use meth_params_module, only : NVAR, URHO, UTEMP

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: kps_l1,kps_l2,kps_l3,kps_h1,kps_h2,kps_h3
  integer, intent(in) :: sta_l1,sta_l2,sta_l3,sta_h1,sta_h2,sta_h3
  real(rt)        , intent(inout) :: kps(kps_l1:kps_h1,kps_l2:kps_h2,kps_l3:kps_h3)
  real(rt)        , intent(in   ) :: sta(sta_l1:sta_h1,sta_l2:sta_h2,sta_l3:sta_h3,NVAR)
  real(rt)        , intent(in) :: k0_p, m_p, n_p
  real(rt)        , intent(in) :: k0_r, m_r, n_r
  real(rt)        , intent(in) :: Tfloor, kfloor

  integer :: i, j, k
  real(rt)        , parameter :: tiny = 1.0e-50_rt
  real(rt)         :: Teff, k_p, k_r

  ! scattering is assumed to be independent of nu.

  if ( m_p.eq.0.e0_rt .and. n_p.eq.0.e0_rt .and. &
       m_r.eq.0.e0_rt .and. n_r.eq.0.e0_rt ) then
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              kps(i,j,k) = k0_r - k0_p
           end do
        end do
     end do
  else 
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              Teff = max(sta(i,j,k,UTEMP), tiny)
              Teff = Teff + Tfloor * exp(-Teff / (Tfloor + tiny))
              k_p = k0_p * (sta(i,j,k,URHO) ** m_p) * (Teff ** (-n_p))
              k_r = k0_r * (sta(i,j,k,URHO) ** m_r) * (Teff ** (-n_r))
              kps(i,j,k) = max(k_r-k_p, kfloor)
           end do
        end do
     end do
  end if

end subroutine ca_compute_scattering_2

subroutine init_godunov_indices_rad() bind(C)

  use meth_params_module, only : GDRHO, GDU, GDV, GDW, GDPRES, GDGAME, ngdnv, &
       GDLAMS, GDERADS, &
       QU, QV, QW
  use rad_params_module, only: ngroups

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  ngdnv = 6 + 2*ngroups
  GDRHO = 1
  GDU = 2
  GDV = 3
  GDW = 4
  GDPRES = 5
  GDGAME = 6
  GDLAMS = GDGAME+1            ! starting index for rad lambda
  GDERADS = GDLAMS + ngroups   ! starting index for rad energy

  ! sanity check
  if ((QU /= GDU) .or. (QV /= GDV) .or. (QW /= GDW)) then
     call bl_error("ERROR: velocity components for godunov and primitive state are not aligned")
  endif

end subroutine init_godunov_indices_rad
