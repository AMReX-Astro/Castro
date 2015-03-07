
!-----------------------------------------------------------------------
!       Compute neutrino group energy density per MeV
!-----------------------------------------------------------------------

      subroutine ca_derneut(neut,n_l1,n_l2,n_l3,n_h1,n_h2,n_h3,nv, &
           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: dnugroup, current_group
        use rad_params_module, only: hplanck, mev2erg
        use rad_params_module, only: radtoE

        implicit none

        integer n_l1,n_l2,n_l3,n_h1,n_h2,n_h3,nv
        integer dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,ncomp
        integer lo(3), hi(3), domlo(3), domhi(3)
        double precision neut(n_l1:n_h1,n_l2:n_h2,n_l3:n_h3,nv)
        double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,ncomp)
        double precision dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp), level, grid_no

        double precision :: spectrum_energy_factor
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

      subroutine ca_derrhoyl(rhoyl,y_l1,y_l2,y_l3,y_h1,y_h2,y_h3,nv, &
           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        implicit none

        integer y_l1,y_l2,y_l3,y_h1,y_h2,y_h3,nv
        integer dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,ncomp
        integer lo(3), hi(3), domlo(3), domhi(3)
        double precision rhoyl(y_l1:y_h1,y_l2:y_h2,y_l3:y_h3,nv)
        double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,ncomp)
        double precision dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp), level, grid_no

        integer          :: i,j,k,n

        double precision :: fac

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

      subroutine ca_deryl(yl,y_l1,y_l2,y_l3,y_h1,y_h2,y_h3,nv, &
           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        implicit none

        integer y_l1,y_l2,y_l3,y_h1,y_h2,y_h3,nv
        integer dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,ncomp
        integer lo(3), hi(3), domlo(3), domhi(3)
        double precision yl(y_l1:y_h1,y_l2:y_h2,y_l3:y_h3,nv)
        double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,ncomp)
        double precision dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp), level, grid_no

        integer          :: i,j,k,n

        double precision :: fac

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

      subroutine ca_derynue(y,y_l1,y_l2,y_l3,y_h1,y_h2,y_h3,nv, &
           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        implicit none

        integer y_l1,y_l2,y_l3,y_h1,y_h2,y_h3,nv
        integer dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,ncomp
        integer lo(3), hi(3), domlo(3), domhi(3)
        double precision y(y_l1:y_h1,y_l2:y_h2,y_l3:y_h3,nv)
        double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,ncomp)
        double precision dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp), level, grid_no

        integer          :: i,j,k,n

        double precision :: fac

        fac = radtoE / (hplanck * avogadro)

        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           y(i,j,k,1) = 0.d0
           do n = 0, ng0-1
              y(i,j,k,1) = y(i,j,k,1) + fac * dat(i,j,k,3+n) / nugroup(n)
           enddo
           y(i,j,k,1) = y(i,j,k,1) / dat(i,j,k,1)

        enddo
        enddo
        enddo

      end subroutine ca_derynue

!-----------------------------------------------------------------------

      subroutine ca_derynuae(y,y_l1,y_l2,y_l3,y_h1,y_h2,y_h3,nv, &
           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,ncomp,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: nugroup, ng0, ng1
        use rad_params_module, only: hplanck, avogadro
        use rad_params_module, only: radtoE

        implicit none

        integer y_l1,y_l2,y_l3,y_h1,y_h2,y_h3,nv
        integer dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,ncomp
        integer lo(3), hi(3), domlo(3), domhi(3)
        double precision y(y_l1:y_h1,y_l2:y_h2,y_l3:y_h3,nv)
        double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,ncomp)
        double precision dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp), level, grid_no

        integer          :: i,j,k,n

        double precision :: fac

        fac = radtoE / (hplanck * avogadro)

        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           y(i,j,k,1) = 0.d0
           do n = ng0, ng0 + ng1 - 1
              y(i,j,k,1) = y(i,j,k,1) + fac * dat(i,j,k,3+n) / nugroup(n)
           enddo
           y(i,j,k,1) = y(i,j,k,1) / dat(i,j,k,1)

        enddo
        enddo
        enddo

      end subroutine ca_derynuae

!-----------------------------------------------------------------------

      subroutine ca_derertot(Et,Et_l1,Et_l2,Et_l3,Et_h1,Et_h2,Et_h3,ncomp_Et, &
           Er,Er_l1,Er_l2,Er_l3,Er_h1,Er_h2,Er_h3,ncomp_Er,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: radtoE

        implicit none

        integer Et_l1,Et_l2,Et_l3,Et_h1,Et_h2,Et_h3,ncomp_Et
        integer Er_l1,Er_l2,Er_l3,Er_h1,Er_h2,Er_h3,ncomp_Er
        integer lo(3), hi(3), domlo(3), domhi(3)
        double precision Et(Et_l1:Et_h1,Et_l2:Et_h2,Et_l3:Et_h3,ncomp_Et)
        double precision Er(Er_l1:Er_h1,Er_l2:Er_h2,Er_l3:Er_h3,ncomp_Er)
        double precision dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp_Er), level, grid_no
      
        integer :: i, j, k, g

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 Et(i,j,k,1) = 0.d0
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

      subroutine ca_derenue(Enue,Enue_l1,Enue_l2,Enue_l3, &
           Enue_h1,Enue_h2,Enue_h3,ncomp_Enue, &
           Er,Er_l1,Er_l2,Er_l3,Er_h1,Er_h2,Er_h3,ncomp_Er,&
           lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: radtoE, ng0

        implicit none

        integer Enue_l1,Enue_l2,Enue_l3,Enue_h1,Enue_h2,Enue_h3,ncomp_Enue
        integer Er_l1,Er_l2,Er_l3,Er_h1,Er_h2,Er_h3,ncomp_Er
        integer lo(3), hi(3), domlo(3), domhi(3)
        double precision Enue(Enue_l1:Enue_h1,Enue_l2:Enue_h2,Enue_l3:Enue_h3,&
             ncomp_Enue)
        double precision Er(Er_l1:Er_h1,Er_l2:Er_h2,Er_l3:Er_h3,ncomp_Er)
        double precision dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp_Er), level, grid_no
      
        integer :: i, j, k, g

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 Enue(i,j,k,1) = 0.d0
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

      subroutine ca_derenuae(Enuae,Enuae_l1,Enuae_l2,Enuae_l3, &
           Enuae_h1,Enuae_h2,Enuae_h3,ncomp_Enuae, &
           Er,Er_l1,Er_l2,Er_l3,Er_h1,Er_h2,Er_h3,ncomp_Er,&
           lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

        use rad_params_module, only: radtoE, ng0, ng1

        implicit none

        integer Enuae_l1,Enuae_l2,Enuae_l3,Enuae_h1,Enuae_h2,Enuae_h3,ncomp_Enuae
        integer Er_l1,Er_l2,Er_l3,Er_h1,Er_h2,Er_h3,ncomp_Er
        integer lo(3), hi(3), domlo(3), domhi(3)
        double precision Enuae(Enuae_l1:Enuae_h1,Enuae_l2:Enuae_h2,Enuae_l3:Enuae_h3,&
             ncomp_Enuae)
        double precision Er(Er_l1:Er_h1,Er_l2:Er_h2,Er_l3:Er_h3,ncomp_Er)
        double precision dx(3), xlo(3), time, dt
        integer bc(3,2,ncomp_Er), level, grid_no
      
        integer :: i, j, k, g

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 Enuae(i,j,k,1) = 0.d0
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

