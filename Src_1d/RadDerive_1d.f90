

!-----------------------------------------------------------------------
!       Compute neutrino group energy density per MeV
!-----------------------------------------------------------------------

      subroutine ca_derneut(neut,n_l1,n_h1,nv, &
           dat,dat_l1,dat_h1,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: dnugroup, current_group
        use rad_params_module, only: hplanck, mev2erg
        use rad_params_module, only: radtoE

        implicit none

        integer n_l1,n_h1,nv
        integer dat_l1,dat_h1,ncomp
        integer lo(1), hi(1), domlo(1), domhi(1)
        double precision neut(n_l1:n_h1,nv)
        double precision dat(dat_l1:dat_h1,ncomp)
        double precision dx(1), xlo(1), time, dt
        integer bc(1,2,ncomp), level, grid_no

        double precision :: spectrum_energy_factor
        integer          :: i,n

        spectrum_energy_factor = radtoE / (hplanck / mev2erg)

        n = current_group

        do i = lo(1),hi(1)
           neut(i,1) = spectrum_energy_factor * dat(i,1) / dnugroup(n)
        enddo

      end subroutine ca_derneut

!-----------------------------------------------------------------------

      subroutine ca_derrhoyl(rhoyl,y_l1,y_h1,nv, &
           dat,dat_l1,dat_h1,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        implicit none

        integer y_l1,y_h1,nv
        integer dat_l1,dat_h1,ncomp
        integer lo(1), hi(1), domlo(1), domhi(1)
        double precision rhoyl(y_l1:y_h1,nv)
        double precision dat(dat_l1:dat_h1,ncomp)
        double precision dx(1), xlo(1), time, dt
        integer bc(1,2,ncomp), level, grid_no

        integer          :: i,n

        double precision :: fac

        fac = radtoE / (hplanck * avogadro)

        do i = lo(1),hi(1)
           rhoyl(i,1) = dat(i,2)
           do n = 0, ng0-1
              rhoyl(i,1) = rhoyl(i,1) + fac * dat(i,3+n) / nugroup(n)
           enddo
           do n = ng0, ng0 + ng1 - 1
              rhoyl(i,1) = rhoyl(i,1) - fac * dat(i,3+n) / nugroup(n)
           enddo
        enddo

      end subroutine ca_derrhoyl

!-----------------------------------------------------------------------

      subroutine ca_deryl(yl,y_l1,y_h1,nv, &
           dat,dat_l1,dat_h1,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        implicit none

        integer y_l1,y_h1,nv
        integer dat_l1,dat_h1,ncomp
        integer lo(1), hi(1), domlo(1), domhi(1)
        double precision yl(y_l1:y_h1,nv)
        double precision dat(dat_l1:dat_h1,ncomp)
        double precision dx(1), xlo(1), time, dt
        integer bc(1,2,ncomp), level, grid_no

        integer          :: i,n

        double precision :: fac

        fac = radtoE / (hplanck * avogadro)

        do i = lo(1),hi(1)
           yl(i,1) = dat(i,2)
           do n = 0, ng0-1
              yl(i,1) = yl(i,1) + fac * dat(i,3+n) / nugroup(n)
           enddo
           do n = ng0, ng0 + ng1 - 1
              yl(i,1) = yl(i,1) - fac * dat(i,3+n) / nugroup(n)
           enddo
           yl(i,1) = yl(i,1) / dat(i,1)
        enddo

      end subroutine ca_deryl

!-----------------------------------------------------------------------

      subroutine ca_derynue(y,y_l1,y_h1,nv, &
           dat,dat_l1,dat_h1,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        implicit none

        integer y_l1,y_h1,nv
        integer dat_l1,dat_h1,ncomp
        integer lo(1), hi(1), domlo(1), domhi(1)
        double precision y(y_l1:y_h1,nv)
        double precision dat(dat_l1:dat_h1,ncomp)
        double precision dx(1), xlo(1), time, dt
        integer bc(1,2,ncomp), level, grid_no

        integer          :: i,n

        double precision :: fac

        fac = radtoE / (hplanck * avogadro)

        do i = lo(1),hi(1)
           y(i,1) = 0.d0
           do n = 0, ng0-1
              y(i,1) = y(i,1) + fac * dat(i,3+n) / nugroup(n)
           enddo
           y(i,1) = y(i,1) / dat(i,1)
        enddo

      end subroutine ca_derynue

!-----------------------------------------------------------------------

      subroutine ca_derynuae(y,y_l1,y_h1,nv, &
           dat,dat_l1,dat_h1,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        implicit none

        integer y_l1,y_h1,nv
        integer dat_l1,dat_h1,ncomp
        integer lo(1), hi(1), domlo(1), domhi(1)
        double precision y(y_l1:y_h1,nv)
        double precision dat(dat_l1:dat_h1,ncomp)
        double precision dx(1), xlo(1), time, dt
        integer bc(1,2,ncomp), level, grid_no

        integer          :: i,n

        double precision :: fac

        fac = radtoE / (hplanck * avogadro)

        do i = lo(1),hi(1)
           y(i,1) = 0.d0
           do n = ng0, ng0 + ng1 - 1
              y(i,1) = y(i,1) + fac * dat(i,3+n) / nugroup(n)
           enddo
           y(i,1) = y(i,1) / dat(i,1)
        enddo

      end subroutine ca_derynuae

!-----------------------------------------------------------------------

      subroutine ca_derertot(Et,Et_l1,Et_h1,ncomp_Et, &
           Er,Er_l1,Er_h1,ncomp_Er,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: radtoE

        implicit none

        integer Et_l1,Et_h1,ncomp_Et
        integer Er_l1,Er_h1,ncomp_Er
        integer lo(1), hi(1), domlo(1), domhi(1)
        double precision Et(Et_l1:Et_h1,ncomp_Et)
        double precision Er(Er_l1:Er_h1,ncomp_Er)
        double precision dx(1), xlo(1), time, dt
        integer bc(1,2,ncomp_Er), level, grid_no
      
        integer :: i, g

        Et(lo(1):hi(1),:) = 0.d0

        do g = 1, ncomp_Er
           do i = lo(1), hi(1)
              Et(i,1) = Et(i,1) + Er(i,g)*radtoE
           end do
        end do

      end subroutine ca_derertot

!-----------------------------------------------------------------------

      subroutine ca_derenue(Enue,Enue_l1,Enue_h1,ncomp_Enue, &
           Er,Er_l1,Er_h1,ncomp_Er,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: radtoE, ng0

        implicit none

        integer Enue_l1,Enue_h1,ncomp_Enue
        integer Er_l1,Er_h1,ncomp_Er
        integer lo(1), hi(1), domlo(1), domhi(1)
        double precision Enue(Enue_l1:Enue_h1,ncomp_Enue)
        double precision Er(Er_l1:Er_h1,ncomp_Er)
        double precision dx(1), xlo(1), time, dt
        integer bc(1,2,ncomp_Er), level, grid_no
      
        integer :: i, g

        Enue(lo(1):hi(1),:) = 0.d0

        do g = 1, ng0
           do i = lo(1), hi(1)
              Enue(i,1) = Enue(i,1) + Er(i,g)*radtoE
           end do
        end do

      end subroutine ca_derenue

!-----------------------------------------------------------------------

      subroutine ca_derenuae(Enuae,Enuae_l1,Enuae_h1,ncomp_Enuae, &
           Er,Er_l1,Er_h1,ncomp_Er,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: radtoE, ng0, ng1

        implicit none

        integer Enuae_l1,Enuae_h1,ncomp_Enuae
        integer Er_l1,Er_h1,ncomp_Er
        integer lo(1), hi(1), domlo(1), domhi(1)
        double precision Enuae(Enuae_l1:Enuae_h1,ncomp_Enuae)
        double precision Er(Er_l1:Er_h1,ncomp_Er)
        double precision dx(1), xlo(1), time, dt
        integer bc(1,2,ncomp_Er), level, grid_no
      
        integer :: i, g

        Enuae(lo(1):hi(1),:) = 0.d0

        do g = ng0+1, ng0+ng1
           do i = lo(1), hi(1)
              Enuae(i,1) = Enuae(i,1) + Er(i,g)*radtoE
           end do
        end do

      end subroutine ca_derenuae

!-----------------------------------------------------------------------

