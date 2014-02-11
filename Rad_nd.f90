
! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize values of physical constants used by the
! ::: radiation package.
! ::: -----------------------------------------------------------

      subroutine ca_initradconstants(p, c, h, k, s, a, m, J_is_used)

        use fundamental_constants_module, only : c_fcm=>c_light, h_fcm=>hplanck, &
             k_fcm=>k_B, s_fcm=>sigma_SB, a_fcm=>n_A, ev2erg_fcm=>ev2erg

        use rad_params_module, only: pi, clight, hplanck
        use rad_params_module, only: kboltz, stefbol, arad, avogadro
        use rad_params_module, only: Hz2MeV, mev2erg, tiny
        use rad_params_module, only: radtoE  !, radtoJ, Etorad, radfluxtoF
        use rad_params_module, only: etafactor

        implicit none

        double precision p, c, h, k, s, a, m
        integer J_is_used

        c = c_fcm
        h = h_fcm
        k = k_fcm
        s = s_fcm
        a = a_fcm
        m = 1.d6 * ev2erg_fcm

        pi       = p
        clight   = c
        hplanck  = h
        kboltz   = k
        stefbol  = s
        arad     = 4.*stefbol/clight
        avogadro = a
        mev2erg  = m
        Hz2MeV   = h / m
        tiny     = 1.d-50

        if (J_is_used > 0) then
           radtoE = 4.d0*pi/clight
!           radtoJ = 1.0d0
!           Etorad = 1.d0/radtoE
!           radfluxtoF = 4.d0*pi
           etafactor = 1.d0
        else
           radtoE = 1.0d0
!           radtoJ = clight/(4.d0*pi)
!           Etorad = 1.0d0
!           radfluxtoF = 1.d0
           etafactor = 4.d0*pi/clight
        end if

      end subroutine ca_initradconstants

! For single group, let set ngroups to 1.
      subroutine ca_initsinglegroup(ngr)
        use rad_params_module, only : ngroups, nugroup, dnugroup, ng0, ng1
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
           nugroup(i)  = 1.d0  ! dummy
           dnugroup(i) = 1.d0
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

        implicit none

        double precision nugr(0:ngr-1), dnugr(0:ngr-1)
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

        implicit none

        double precision, intent(in) :: nugr(0:ngr-1), dnugr(0:ngr-1), xnugr(0:ngr)
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

        implicit none

        double precision, intent(in) :: nugr(0:ngr-1), dnugr(0:ngr-1), dlognugr(0:ngr-1), xnugr(0:ngr+2)
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

        erg2rhoYe = 0.d0
        if (ng0 > 0) then
           erg2rhoYe(0:ng0-1) = 1.d0 / (avogadro*hplanck*nugroup(0:ng0-1))
           if (ng1 > 0) then
              erg2rhoYe(ng0:ng0+ng1-1) = -1.d0 / (avogadro*hplanck*nugroup(ng0:ng0+ng1-1))
           end if
        end if

      end subroutine ca_initgroups3

! ::: -----------------------------------------------------------

      subroutine ca_setgroup(igroup)

        use rad_params_module, only: current_group

        implicit none
        integer igroup

        current_group = igroup

      end subroutine ca_setgroup

