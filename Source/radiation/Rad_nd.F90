subroutine ca_init_fort_constants(hplanck_in, avogadro_in) bind(C, name="ca_init_fort_constants")

  use rad_params_module, only: hplanck, avogadro
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), intent(in), value :: hplanck_in, avogadro_in

  hplanck = hplanck_in
  avogadro = avogadro_in

end subroutine ca_init_fort_constants



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


