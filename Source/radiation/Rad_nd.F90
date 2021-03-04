

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
