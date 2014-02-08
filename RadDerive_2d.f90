
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

!        print *,'executing ca_initgroups'
!
!        print *, "groups are: "
!        do i = 0, ngroups - 1
!           print *, i, nugroup(i), dnugroup(i)
!        enddo

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

        use rad_params_module, only: ngroups, ng0, ng1, nnuspec, nugroup, nradspec, dnugroup, &
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

!-----------------------------------------------------------------------
!       Compute neutrino group energy density per MeV
!-----------------------------------------------------------------------

      subroutine ca_derneut(neut,n_l1,n_l2,n_h1,n_h2,nv, &
           dat,dat_l1,dat_l2,dat_h1,dat_h2,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: dnugroup, current_group
        use rad_params_module, only: hplanck, mev2erg
        use rad_params_module, only: radtoE

        implicit none

        integer n_l1,n_l2,n_h1,n_h2,nv
        integer dat_l1,dat_l2,dat_h1,dat_h2,ncomp
        integer lo(2), hi(2), domlo(2), domhi(2)
        double precision neut(n_l1:n_h1,n_l2:n_h2,nv)
        double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,ncomp)
        double precision dx(2), xlo(2), time, dt
        integer bc(2,2,ncomp), level, grid_no

        double precision :: spectrum_energy_factor
        integer          :: i,j,n

        spectrum_energy_factor = radtoE / (hplanck / mev2erg)

        n = current_group

        do j = lo(2),hi(2)
           do i = lo(1),hi(1)
              neut(i,j,1) = spectrum_energy_factor * dat(i,j,1) / dnugroup(n)
           enddo
        enddo

      end subroutine ca_derneut

!-----------------------------------------------------------------------

      subroutine ca_derrhoyl(rhoyl,y_l1,y_l2,y_h1,y_h2,nv, &
           dat,dat_l1,dat_l2,dat_h1,dat_h2,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        implicit none

        integer y_l1,y_l2,y_h1,y_h2,nv
        integer dat_l1,dat_l2,dat_h1,dat_h2,ncomp
        integer lo(2), hi(2), domlo(2), domhi(2)
        double precision rhoyl(y_l1:y_h1,y_l2:y_h2,nv)
        double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,ncomp)
        double precision dx(2), xlo(2), time, dt
        integer bc(2,2,ncomp), level, grid_no

        integer          :: i,j,n

        double precision :: fac

        fac = radtoE / (hplanck * avogadro)

        do j = lo(2),hi(2)
           do i = lo(1),hi(1)
              rhoyl(i,j,1) = dat(i,j,2)
              do n = 0, ng0-1
                 rhoyl(i,j,1) = rhoyl(i,j,1) + fac * dat(i,j,3+n) / nugroup(n)
              enddo
              do n = ng0, ng0 + ng1 - 1
                 rhoyl(i,j,1) = rhoyl(i,j,1) - fac * dat(i,j,3+n) / nugroup(n)
              enddo
           enddo
        enddo

      end subroutine ca_derrhoyl

!-----------------------------------------------------------------------

      subroutine ca_deryl(yl,y_l1,y_l2,y_h1,y_h2,nv, &
           dat,dat_l1,dat_l2,dat_h1,dat_h2,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        implicit none

        integer y_l1,y_l2,y_h1,y_h2,nv
        integer dat_l1,dat_l2,dat_h1,dat_h2,ncomp
        integer lo(2), hi(2), domlo(2), domhi(2)
        double precision yl(y_l1:y_h1,y_l2:y_h2,nv)
        double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,ncomp)
        double precision dx(2), xlo(2), time, dt
        integer bc(2,2,ncomp), level, grid_no

        integer          :: i,j,n

        double precision :: fac

        fac = radtoE / (hplanck * avogadro)

        do j = lo(2),hi(2)
           do i = lo(1),hi(1)
              yl(i,j,1) = dat(i,j,2)
              do n = 0, ng0-1
                 yl(i,j,1) = yl(i,j,1) + fac * dat(i,j,3+n) / nugroup(n)
              enddo
              do n = ng0, ng0 + ng1 - 1
                 yl(i,j,1) = yl(i,j,1) - fac * dat(i,j,3+n) / nugroup(n)
              enddo
              yl(i,j,1) = yl(i,j,1) / dat(i,j,1)
           enddo
        enddo

      end subroutine ca_deryl

!-----------------------------------------------------------------------

      subroutine ca_derynue(y,y_l1,y_l2,y_h1,y_h2,nv, &
           dat,dat_l1,dat_l2,dat_h1,dat_h2,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        implicit none

        integer y_l1,y_l2,y_h1,y_h2,nv
        integer dat_l1,dat_l2,dat_h1,dat_h2,ncomp
        integer lo(2), hi(2), domlo(2), domhi(2)
        double precision y(y_l1:y_h1,y_l2:y_h2,nv)
        double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,ncomp)
        double precision dx(2), xlo(2), time, dt
        integer bc(2,2,ncomp), level, grid_no

        integer          :: i,j,n

        double precision :: fac

        fac = radtoE / (hplanck * avogadro)

        do j = lo(2),hi(2)
           do i = lo(1),hi(1)
              y(i,j,1) = 0.d0
              do n = 0, ng0-1
                 y(i,j,1) = y(i,j,1) + fac * dat(i,j,3+n) / nugroup(n)
              enddo
              y(i,j,1) = y(i,j,1) / dat(i,j,1)
           enddo
        enddo

      end subroutine ca_derynue

!-----------------------------------------------------------------------

      subroutine ca_derynuae(y,y_l1,y_l2,y_h1,y_h2,nv, &
           dat,dat_l1,dat_l2,dat_h1,dat_h2,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        implicit none

        integer y_l1,y_l2,y_h1,y_h2,nv
        integer dat_l1,dat_l2,dat_h1,dat_h2,ncomp
        integer lo(2), hi(2), domlo(2), domhi(2)
        double precision y(y_l1:y_h1,y_l2:y_h2,nv)
        double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,ncomp)
        double precision dx(2), xlo(2), time, dt
        integer bc(2,2,ncomp), level, grid_no

        integer          :: i,j,n

        double precision :: fac

        fac = radtoE / (hplanck * avogadro)

        do j = lo(2),hi(2)
           do i = lo(1),hi(1)
              y(i,j,1) = 0.d0
              do n = ng0, ng0 + ng1 - 1
                 y(i,j,1) = y(i,j,1) + fac * dat(i,j,3+n) / nugroup(n)
              enddo
              y(i,j,1) = y(i,j,1) / dat(i,j,1)
           enddo
        enddo

      end subroutine ca_derynuae

!-----------------------------------------------------------------------

      subroutine ca_derertot(Et,Et_l1,Et_l2,Et_h1,Et_h2,ncomp_Et, &
           Er,Er_l1,Er_l2,Er_h1,Er_h2,ncomp_Er,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: radtoE

        implicit none

        integer Et_l1,Et_l2,Et_h1,Et_h2,ncomp_Et
        integer Er_l1,Er_l2,Er_h1,Er_h2,ncomp_Er
        integer lo(2), hi(2), domlo(2), domhi(2)
        double precision Et(Et_l1:Et_h1,Et_l2:Et_h2,ncomp_Et)
        double precision Er(Er_l1:Er_h1,Er_l2:Er_h2,ncomp_Er)
        double precision dx(2), xlo(2), time, dt
        integer bc(2,2,ncomp_Er), level, grid_no
      
        integer :: i, j, g

        Et = 0.d0

        do g = 1, ncomp_Er
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 Et(i,j,1) = Et(i,j,1) + Er(i,j,g)*radtoE
              end do
           end do
        end do

      end subroutine ca_derertot

!-----------------------------------------------------------------------

      subroutine ca_derenue(Enue,Enue_l1,Enue_l2,Enue_h1,Enue_h2,ncomp_Enue, &
           Er,Er_l1,Er_l2,Er_h1,Er_h2,ncomp_Er,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: radtoE, ng0

        implicit none

        integer Enue_l1,Enue_l2,Enue_h1,Enue_h2,ncomp_Enue
        integer Er_l1,Er_l2,Er_h1,Er_h2,ncomp_Er
        integer lo(2), hi(2), domlo(2), domhi(2)
        double precision Enue(Enue_l1:Enue_h1,Enue_l2:Enue_h2,ncomp_Enue)
        double precision Er(Er_l1:Er_h1,Er_l2:Er_h2,ncomp_Er)
        double precision dx(2), xlo(2), time, dt
        integer bc(2,2,ncomp_Er), level, grid_no
      
        integer :: i, j, g

        Enue = 0.d0

        do g = 1, ng0
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 Enue(i,j,1) = Enue(i,j,1) + Er(i,j,g)*radtoE
              end do
           end do
        end do

      end subroutine ca_derenue

!-----------------------------------------------------------------------

      subroutine ca_derenuae(Enuae,Enuae_l1,Enuae_l2,Enuae_h1,Enuae_h2,ncomp_Enuae, &
           Er,Er_l1,Er_l2,Er_h1,Er_h2,ncomp_Er,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: radtoE, ng0, ng1

        implicit none

        integer Enuae_l1,Enuae_l2,Enuae_h1,Enuae_h2,ncomp_Enuae
        integer Er_l1,Er_l2,Er_h1,Er_h2,ncomp_Er
        integer lo(2), hi(2), domlo(2), domhi(2)
        double precision Enuae(Enuae_l1:Enuae_h1,Enuae_l2:Enuae_h2,ncomp_Enuae)
        double precision Er(Er_l1:Er_h1,Er_l2:Er_h2,ncomp_Er)
        double precision dx(2), xlo(2), time, dt
        integer bc(2,2,ncomp_Er), level, grid_no
      
        integer :: i, j, g

        Enuae = 0.d0

        do g = ng0+1, ng0+ng1
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 Enuae(i,j,1) = Enuae(i,j,1) + Er(i,j,g)*radtoE
              end do
           end do
        end do

      end subroutine ca_derenuae

!-----------------------------------------------------------------------

