
!-----------------------------------------------------------------------
!       Compute neutrino group energy density per MeV
!-----------------------------------------------------------------------

      subroutine ca_derneut(neut, n_lo, n_hi, nv, &
                            dat, d_lo, d_hi, ncomp, &
                            lo, hi, domlo, domhi, &
                            dx, xlo, time, dt, bc, level, grid_no) bind(C)

        use rad_params_module, only: dnugroup, current_group
        use rad_params_module, only: hplanck, mev2erg
        use rad_params_module, only: radtoE
        use amrex_fort_module, only: rt => amrex_real

        implicit none

        integer,  intent(in   ) :: n_lo(3), n_hi(3), nv
        integer,  intent(in   ) :: d_lo(3), d_hi(3), ncomp
        integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
        real(rt), intent(inout) :: neut(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3),nv)
        real(rt), intent(in   ) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp)
        real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
        integer,  intent(in   ) :: bc(3,2,ncomp), level, grid_no

        real(rt) :: spectrum_energy_factor
        integer  :: i,j,k,n

        spectrum_energy_factor = radtoE / (hplanck / mev2erg)

        n = current_group

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)

                 neut(i,j,k,1) = spectrum_energy_factor * dat(i,j,k,1) / dnugroup(n)

              end do
           end do
        end do

      end subroutine ca_derneut

!-----------------------------------------------------------------------

      subroutine ca_derrhoyl(rhoyl, y_lo, y_hi, nv, &
                             dat, d_lo, d_hi, ncomp, &
                             lo, hi, domlo, domhi, &
                             dx, xlo, time, dt, bc, level, grid_no) bind(C)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE
        use amrex_fort_module, only: rt => amrex_real

        implicit none

        integer,  intent(in   ) :: y_lo(3), y_hi(3), nv
        integer,  intent(in   ) :: d_lo(3), d_hi(3), ncomp
        integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
        real(rt), intent(inout) :: rhoyl(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nv)
        real(rt), intent(in   ) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp)
        real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
        integer,  intent(in   ) :: bc(3,2,ncomp), level, grid_no

        integer :: i, j, k, n

        real(rt) :: fac

        fac = radtoE / (hplanck * avogadro)

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)

                 rhoyl(i,j,k,1) = dat(i,j,k,2)
                 do n = 0, ng0-1
                    rhoyl(i,j,k,1) = rhoyl(i,j,k,1) + fac*dat(i,j,k,3+n) / nugroup(n)
                 end do
                 do n = ng0, ng0 + ng1 - 1
                    rhoyl(i,j,k,1) = rhoyl(i,j,k,1) - fac*dat(i,j,k,3+n) / nugroup(n)
                 end do

              end do
           end do
        enddo

      end subroutine ca_derrhoyl

!-----------------------------------------------------------------------

      subroutine ca_deryl(yl, y_lo, y_hi, nv, &
                          dat, d_lo, d_hi, ncomp, &
                          lo, hi, domlo, domhi, &
                          dx, xlo, time, dt, bc, level, grid_no) bind(C)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE
        use amrex_fort_module, only: rt => amrex_real

        implicit none

        integer,  intent(in   ) :: y_lo(3), y_hi(3), nv
        integer,  intent(in   ) :: d_lo(3), d_hi(3), ncomp
        integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
        real(rt), intent(inout) :: yl(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nv)
        real(rt), intent(in   ) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp)
        real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
        integer,  intent(in   ) :: bc(3,2,ncomp), level, grid_no

        integer :: i,j,k,n

        real(rt) :: fac

        fac = radtoE / (hplanck * avogadro)

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)

                 yl(i,j,k,1) = dat(i,j,k,2)
                 do n = 0, ng0-1
                    yl(i,j,k,1) = yl(i,j,k,1) + fac * dat(i,j,k,3+n) / nugroup(n)
                 end do
                 do n = ng0, ng0 + ng1 - 1
                    yl(i,j,k,1) = yl(i,j,k,1) - fac * dat(i,j,k,3+n) / nugroup(n)
                 end do
                 yl(i,j,k,1) = yl(i,j,k,1) / dat(i,j,k,1)

              end do
           end do
        end do

      end subroutine ca_deryl

!-----------------------------------------------------------------------

      subroutine ca_derynue(y, y_lo, y_hi, nv, &
                            dat, d_lo, d_hi, ncomp, &
                            lo, hi, domlo, domhi, &
                            dx, xlo, time, dt, bc, level, grid_no) bind(C)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE
        use amrex_fort_module, only: rt => amrex_real

        implicit none

        integer,  intent(in   ) :: y_lo(3), y_hi(3), nv
        integer,  intent(in   ) :: d_lo(3), d_hi(3), ncomp
        integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
        real(rt), intent(inout) :: y(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nv)
        real(rt), intent(in   ) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp)
        real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
        integer,  intent(in   ) :: bc(3,2,ncomp), level, grid_no

        integer :: i,j,k,n

        real(rt) :: fac

        fac = radtoE / (hplanck * avogadro)

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)

                 y(i,j,k,1) = 0.e0_rt
                 do n = 0, ng0-1
                    y(i,j,k,1) = y(i,j,k,1) + fac * dat(i,j,k,3+n) / nugroup(n)
                 end do
                 y(i,j,k,1) = y(i,j,k,1) / dat(i,j,k,1)

              end do
           end do
        end do

      end subroutine ca_derynue

!-----------------------------------------------------------------------

      subroutine ca_derynuae(y, y_lo, y_hi, nv, &
                             dat, d_lo, d_hi, ncomp, &
                             lo, hi, domlo, domhi, &
                             dx, xlo, time, dt, bc, level, grid_no) bind(C)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE
        use amrex_fort_module, only: rt => amrex_real

        implicit none

        integer,  intent(in   ) :: y_lo(3), y_hi(3), nv
        integer,  intent(in   ) :: d_lo(3), d_hi(3), ncomp
        integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
        real(rt), intent(inout) :: y(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nv)
        real(rt), intent(in   ) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp)
        real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
        integer,  intent(in   ) :: bc(3,2,ncomp), level, grid_no

        integer :: i,j,k,n

        real(rt) :: fac

        fac = radtoE / (hplanck * avogadro)

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)

                 y(i,j,k,1) = 0.e0_rt
                 do n = ng0, ng0 + ng1 - 1
                    y(i,j,k,1) = y(i,j,k,1) + fac * dat(i,j,k,3+n) / nugroup(n)
                 end do
                 y(i,j,k,1) = y(i,j,k,1) / dat(i,j,k,1)

              end do
           end do
        end do

      end subroutine ca_derynuae

!-----------------------------------------------------------------------

      subroutine ca_derertot(Et, Et_lo, Et_hi, ncomp_Et, &
                             Er, Er_lo, Er_hi, ncomp_Er, &
                             lo,hi, domlo, domhi, &
                             dx, xlo, time, dt, bc, level, grid_no) bind(C)

        use rad_params_module, only: radtoE
        use amrex_fort_module, only: rt => amrex_real

        implicit none

        integer,  intent(in   ) :: Et_lo(3), Et_hi(3), ncomp_Et
        integer,  intent(in   ) :: Er_lo(3), Er_hi(3), ncomp_Er
        integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
        real(rt), intent(inout) :: Et(Et_lo(1):Et_hi(1),Et_lo(2):Et_hi(2),Et_lo(3):Et_hi(3),ncomp_Et)
        real(rt), intent(in   ) :: Er(Er_lo(1):Er_hi(1),Er_lo(2):Er_hi(2),Er_lo(3):Er_hi(3),ncomp_Er)
        real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
        integer,  intent(in   ) :: bc(3,2,ncomp_Er), level, grid_no
      
        integer :: i, j, k, g

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 Et(i,j,k,1) = 0.e0_rt
              end do
           end do
        end do

        do g = 1, ncomp_Er
           do k = lo(3), hi(3)
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)
                    Et(i,j,k,1) = Et(i,j,k,1) + Er(i,j,k,g)*radtoE
                 end do
              end do
           end do
        end do

      end subroutine ca_derertot

!-----------------------------------------------------------------------

      subroutine ca_derenue(Enue, En_lo, En_hi, ncomp_En, &
                            Er, Er_lo, Er_hi, ncomp_Er, &
                            lo, hi, domlo, domhi, &
                            dx, xlo, time, dt, bc, level, grid_no) bind(C)

        use rad_params_module, only: radtoE, ng0
        use amrex_fort_module, only: rt => amrex_real

        implicit none

        integer,  intent(in   ) :: En_lo(3), En_hi(3), ncomp_En
        integer,  intent(in   ) :: Er_lo(3), Er_hi(3), ncomp_Er
        integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
        real(rt), intent(inout) :: Enue(En_lo(1):En_hi(1),En_lo(2):En_hi(2),En_lo(3):En_hi(3),ncomp_En)
        real(rt), intent(in   ) :: Er(Er_lo(1):Er_hi(1),Er_lo(2):Er_hi(2),Er_lo(3):Er_hi(3),ncomp_Er)
        real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
        integer,  intent(in   ) :: bc(3,2,ncomp_Er), level, grid_no
      
        integer :: i, j, k, g

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 Enue(i,j,k,1) = 0.e0_rt
              end do
           end do
        end do

        do g = 1, ng0
           do k = lo(3), hi(3)
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)
                    Enue(i,j,k,1) = Enue(i,j,k,1) + Er(i,j,k,g)*radtoE
                 end do
              end do
           end do
        end do

      end subroutine ca_derenue

!-----------------------------------------------------------------------

      subroutine ca_derenuae(Enuae, En_lo, En_hi, ncomp_En, &
                             Er, Er_lo, Er_hi, ncomp_Er, &
                             lo, hi, domlo, domhi, &
                             dx, xlo, time, dt, bc, level, grid_no) bind(C)

        use rad_params_module, only: radtoE, ng0, ng1
        use amrex_fort_module, only: rt => amrex_real

        implicit none

        integer,  intent(in   ) :: En_lo(3), En_hi(3), ncomp_En
        integer,  intent(in   ) :: Er_lo(3), Er_hi(3), ncomp_Er
        integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
        real(rt), intent(inout) :: Enuae(En_lo(1):En_hi(1),En_lo(2):En_hi(2),En_lo(3):En_hi(3),ncomp_En)
        real(rt), intent(in   ) :: Er(Er_lo(1):Er_hi(1),Er_lo(2):Er_hi(2),Er_lo(3):Er_hi(3),ncomp_Er)
        real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
        integer,  intent(in   ) :: bc(3,2,ncomp_Er), level, grid_no
      
        integer :: i, j, k, g

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 Enuae(i,j,k,1) = 0.e0_rt
              end do
           end do
        end do
        
        do g = ng0+1, ng0+ng1
           do k = lo(3), hi(3)
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)
                    Enuae(i,j,k,1) = Enuae(i,j,k,1) + Er(i,j,k,g)*radtoE
                 end do
              end do
           end do
        end do

      end subroutine ca_derenuae

!-----------------------------------------------------------------------
