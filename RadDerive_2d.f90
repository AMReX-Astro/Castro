
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

