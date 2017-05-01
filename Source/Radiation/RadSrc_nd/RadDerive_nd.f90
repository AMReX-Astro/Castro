
!-----------------------------------------------------------------------
!       Compute neutrino group energy density per MeV
!-----------------------------------------------------------------------

      subroutine ca_derneut(neut,n_lo,n_hi,nv, &
                            dat,d_lo,d_hi,ncomp,lo,hi,domlo,domhi, &
                            dx,xlo,time,dt,bc,level,grid_no) bind(C)

        use rad_params_module, only: dnugroup, current_group
        use rad_params_module, only: hplanck, mev2erg
        use rad_params_module, only: radtoE

        use amrex_fort_module, only : rt => amrex_real
        implicit none

        integer n_lo(3),n_hi(3),nv
        integer d_lo(3),d_hi(3),ncomp
        integer lo(3), hi(3), domlo(3), domhi(3)
        real(rt)         neut(n_lo(1):n_hi(1),n_lo(2):n_hi(2),n_lo(3):n_hi(3),nv)
        real(rt)         dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp)
        real(rt)         dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp), level, grid_no

        real(rt)         :: spectrum_energy_factor
        integer          :: i,j,k,n

        spectrum_energy_factor = radtoE / (hplanck / mev2erg)

        n = current_group

        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           neut(i,j,k,1) = spectrum_energy_factor * dat(i,j,k,1) / dnugroup(n)

        enddo
        enddo
        enddo

      end subroutine ca_derneut

!-----------------------------------------------------------------------

      subroutine ca_derrhoyl(rhoyl,y_lo,y_hi,nv, &
                             dat,d_lo,d_hi,ncomp,lo,hi,domlo,domhi, &
                             dx,xlo,time,dt,bc,level,grid_no) bind(C)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        use amrex_fort_module, only : rt => amrex_real
        implicit none

        integer y_lo(3),y_hi(3),nv
        integer d_lo(3),d_hi(3),ncomp
        integer lo(3), hi(3), domlo(3), domhi(3)
        real(rt)         rhoyl(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nv)
        real(rt)         dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp)
        real(rt)         dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp), level, grid_no

        integer          :: i,j,k,n

        real(rt)         :: fac

        fac = radtoE / (hplanck * avogadro)

        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           rhoyl(i,j,k,1) = dat(i,j,k,2)
           do n = 0, ng0-1
              rhoyl(i,j,k,1) = rhoyl(i,j,k,1) + fac*dat(i,j,k,3+n) / nugroup(n)
           enddo
           do n = ng0, ng0 + ng1 - 1
              rhoyl(i,j,k,1) = rhoyl(i,j,k,1) - fac*dat(i,j,k,3+n) / nugroup(n)
           enddo

        enddo
        enddo
        enddo

      end subroutine ca_derrhoyl

!-----------------------------------------------------------------------

      subroutine ca_deryl(yl,y_lo,y_hi,nv, &
                          dat,d_lo,d_hi,ncomp,lo,hi,domlo,domhi, &
                          dx,xlo,time,dt,bc,level,grid_no) bind(C)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        use amrex_fort_module, only : rt => amrex_real
        implicit none

        integer y_lo(3),y_hi(3),nv
        integer d_lo(3),d_hi(3),ncomp
        integer lo(3), hi(3), domlo(3), domhi(3)
        real(rt)         yl(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nv)
        real(rt)         dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp)
        real(rt)         dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp), level, grid_no

        integer          :: i,j,k,n

        real(rt)         :: fac

        fac = radtoE / (hplanck * avogadro)

        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           yl(i,j,k,1) = dat(i,j,k,2)
           do n = 0, ng0-1
              yl(i,j,k,1) = yl(i,j,k,1) + fac * dat(i,j,k,3+n) / nugroup(n)
           enddo
           do n = ng0, ng0 + ng1 - 1
              yl(i,j,k,1) = yl(i,j,k,1) - fac * dat(i,j,k,3+n) / nugroup(n)
           enddo
           yl(i,j,k,1) = yl(i,j,k,1) / dat(i,j,k,1)

        enddo
        enddo
        enddo

      end subroutine ca_deryl

!-----------------------------------------------------------------------

      subroutine ca_derynue(y,y_lo,y_hi,nv, &
                            dat,d_lo,d_hi,ncomp,lo,hi,domlo,domhi, &
                            dx,xlo,time,dt,bc,level,grid_no) bind(C)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        use amrex_fort_module, only : rt => amrex_real
        implicit none

        integer y_lo(3),y_hi(3),nv
        integer d_lo(3),d_hi(3),ncomp
        integer lo(3), hi(3), domlo(3), domhi(3)
        real(rt)         y(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nv)
        real(rt)         dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp)
        real(rt)         dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp), level, grid_no

        integer          :: i,j,k,n

        real(rt)         :: fac

        fac = radtoE / (hplanck * avogadro)

        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           y(i,j,k,1) = 0.e0_rt
           do n = 0, ng0-1
              y(i,j,k,1) = y(i,j,k,1) + fac * dat(i,j,k,3+n) / nugroup(n)
           enddo
           y(i,j,k,1) = y(i,j,k,1) / dat(i,j,k,1)

        enddo
        enddo
        enddo

      end subroutine ca_derynue

!-----------------------------------------------------------------------

      subroutine ca_derynuae(y,y_lo,y_hi,nv, &
                             dat,d_lo,d_hi,ncomp,lo,hi,domlo,domhi, &
                             dx,xlo,time,dt,bc,level,grid_no) bind(C)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        use amrex_fort_module, only : rt => amrex_real
        implicit none

        integer y_lo(3),y_hi(3),nv
        integer d_lo(3),d_hi(3),ncomp
        integer lo(3), hi(3), domlo(3), domhi(3)
        real(rt)         y(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3),nv)
        real(rt)         dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),ncomp)
        real(rt)         dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp), level, grid_no

        integer          :: i,j,k,n

        real(rt)         :: fac

        fac = radtoE / (hplanck * avogadro)

        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           y(i,j,k,1) = 0.e0_rt
           do n = ng0, ng0 + ng1 - 1
              y(i,j,k,1) = y(i,j,k,1) + fac * dat(i,j,k,3+n) / nugroup(n)
           enddo
           y(i,j,k,1) = y(i,j,k,1) / dat(i,j,k,1)

        enddo
        enddo
        enddo

      end subroutine ca_derynuae

!-----------------------------------------------------------------------

      subroutine ca_derertot(Et,Et_lo,Et_hi,ncomp_Et, &
                             Er,Er_lo,Er_hi,ncomp_Er,lo,hi,domlo, &
                             domhi,dx,xlo,time,dt,bc,level,grid_no) bind(C)

        use rad_params_module, only: radtoE

        use amrex_fort_module, only : rt => amrex_real
        implicit none

        integer Et_lo(3),Et_hi(3),ncomp_Et
        integer Er_lo(3),Er_hi(3),ncomp_Er
        integer lo(3), hi(3), domlo(3), domhi(3)
        real(rt)         Et(Et_lo(1):Et_hi(1),Et_lo(2):Et_hi(2),Et_lo(3):Et_hi(3),ncomp_Et)
        real(rt)         Er(Er_lo(1):Er_hi(1),Er_lo(2):Er_hi(2),Er_lo(3):Er_hi(3),ncomp_Er)
        real(rt)         dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp_Er), level, grid_no
      
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

      subroutine ca_derenue(Enue,Enue_lo,Enue_hi,ncomp_Enue, &
                            Er,Er_lo,Er_hi,ncomp_Er,lo,hi,domlo,domhi, &
                            dx,xlo,time,dt,bc,level,grid_no) bind(C)

        use rad_params_module, only: radtoE, ng0

        use amrex_fort_module, only : rt => amrex_real
        implicit none

        integer Enue_lo(3),Enue_hi(3),ncomp_Enue
        integer Er_lo(3),Er_hi(3),ncomp_Er
        integer lo(3), hi(3), domlo(3), domhi(3)
        real(rt)         Enue(Enue_lo(1):Enue_hi(1),Enue_lo(2):Enue_hi(2),Enue_lo(3):Enue_hi(3),&
             ncomp_Enue)
        real(rt)         Er(Er_lo(1):Er_hi(1),Er_lo(2):Er_hi(2),Er_lo(3):Er_hi(3),ncomp_Er)
        real(rt)         dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp_Er), level, grid_no
      
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

      subroutine ca_derenuae(Enuae,Enuae_lo,Enuae_hi,ncomp_Enuae, &
                             Er,Er_lo,Er_hi,ncomp_Er,lo,hi,domlo,domhi, &
                             dx,xlo,time,dt,bc,level,grid_no) bind(C)

        use rad_params_module, only: radtoE, ng0, ng1

        use amrex_fort_module, only : rt => amrex_real
        implicit none

        integer Enuae_lo(3),Enuae_hi(3),ncomp_Enuae
        integer Er_lo(3),Er_hi(3),ncomp_Er
        integer lo(3), hi(3), domlo(3), domhi(3)
        real(rt)         Enuae(Enuae_lo(1):Enuae_hi(1),Enuae_lo(2):Enuae_hi(2),Enuae_lo(3):Enuae_hi(3),&
             ncomp_Enuae)
        real(rt)         Er(Er_lo(1):Er_hi(1),Er_lo(2):Er_hi(2),Er_lo(3):Er_hi(3),ncomp_Er)
        real(rt)         dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp_Er), level, grid_no
      
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

