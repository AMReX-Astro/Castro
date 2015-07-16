!-----------------------------------------------------------------------

      subroutine ca_derstate(state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3,nv, &
                             dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                             domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! The incoming   "dat" vector contains (rho,T,(rho X)_1)
      ! The outgoing "state" vector contains (rho,T,X_1)
      !
      implicit none 

      integer          lo(3), hi(3)
      integer          state_l1,state_l2,state_l3,state_h1,state_h2,state_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no
 
      integer i,j,k

      if (nv .ne. 3) then
          print *,'... confusion in derstate ... nv should be 3 but is ',nv
          call bl_error('Error:: Derive_3d.f90 :: ca_derstate')
      end if

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               !
               ! Density
               !
               state(i,j,k,1) = dat(i,j,k,1)
               !
               ! Temperature
               !
               state(i,j,k,2) = dat(i,j,k,2)
               !
               ! (rho X)_1 --> X_1
               !
               state(i,j,k,3) = dat(i,j,k,3) / dat(i,j,k,1)
            end do
         end do
      end do
 
      end subroutine ca_derstate

!-----------------------------------------------------------------------

      subroutine ca_dervel(vel,vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv, &
                           dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive the velocity from the momentum.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision vel(vel_l1:vel_h1,vel_l2:vel_h2,vel_l3:vel_h3,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no
 
      integer i,j,k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               vel(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
            end do
         end do
      end do
 
      end subroutine ca_dervel

!-----------------------------------------------------------------------

      subroutine ca_dermagvel(magvel,vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive magnitude of velocity.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision magvel(vel_l1:vel_h1,vel_l2:vel_h2,vel_l3:vel_h3,nv)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k
      double precision :: dat1inv

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               dat1inv = 1.d0/dat(i,j,k,1)
               magvel(i,j,k,1) = sqrt( (dat(i,j,k,2) * dat1inv)**2 + &
                                       (dat(i,j,k,3) * dat1inv)**2 + &
                                       (dat(i,j,k,4) * dat1inv)**2 )
            end do
         end do
      end do

      end subroutine ca_dermagvel

!-----------------------------------------------------------------------

      subroutine ca_dermaggrav(maggrav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,ng, &
                               dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                               domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive magnitude of the gravity vector.
      !
      implicit none 

      integer          lo(3), hi(3)
      integer          grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,ng
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision maggrav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3,ng)
      double precision     dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               maggrav(i,j,k,1) = sqrt( dat(i,j,k,1)**2  + &
                                        dat(i,j,k,2)**2  + &
                                        dat(i,j,k,3)**2 )
            end do
         end do
      end do

      end subroutine ca_dermaggrav

!-----------------------------------------------------------------------

      subroutine ca_derradialvel(radvel,vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv, &
                                 dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                                 domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive the radial velocity.
      !
      use bl_constants_module
      use prob_params_module, only : center

      implicit none

      integer          lo(3), hi(3)
      integer          vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision radvel(vel_l1:vel_h1,vel_l2:vel_h2,vel_l3:vel_h3,nv)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer          :: i,j,k
      double precision :: x,y,z,r

      do k = lo(3), hi(3)
         z = xlo(3) + (dble(k-lo(3))+HALF) * delta(3) - center(3)
         do j = lo(2), hi(2)
            y = xlo(2) + (dble(j-lo(2))+HALF) * delta(2) - center(2)
            do i = lo(1), hi(1)
               x = xlo(1) + (dble(i-lo(1))+HALF) * delta(1) - center(1)
               r = sqrt(x*x+y*y+z*z)
               radvel(i,j,k,1) = ( dat(i,j,k,2)*x + &
                                   dat(i,j,k,3)*y + &
                                   dat(i,j,k,4)*z ) / ( dat(i,j,k,1)*r )
            end do
         end do
      end do

      end subroutine ca_derradialvel

!-----------------------------------------------------------------------

      subroutine ca_dermagmom(magmom,mom_l1,mom_l2,mom_l3,mom_h1,mom_h2,mom_h3,nv, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive magnitude of momentum.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          mom_l1,mom_l2,mom_l3,mom_h1,mom_h2,mom_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision magmom(mom_l1:mom_h1,mom_l2:mom_h2,mom_l3:mom_h3,nv)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               magmom(i,j,k,1) = sqrt( dat(i,j,k,1)**2 + dat(i,j,k,2)**2 + dat(i,j,k,3)**2 )
            end do
         end do
      end do

      end subroutine ca_dermagmom

!-----------------------------------------------------------------------

      subroutine ca_derpres(p,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3,ncomp_p, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use eos_type_module
      use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer p_l1,p_l2,p_l3,p_h1,p_h2,p_h3,ncomp_p
      integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer lo(3), hi(3), domlo(3), domhi(3)
      double precision p(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3,ncomp_p)
      double precision u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      double precision dx(3), xlo(3), time, dt
      integer bc(3,2,ncomp_u), level, grid_no

      double precision :: rhoInv
      integer          :: i,j,k

      type (eos_t_3D) :: eos_state

      eos_state = eos_t_3D(lo,hi)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhoInv = ONE / u(i,j,k,URHO)
               
               eos_state % rho(i,j,k) = u(i,j,k,URHO)
               eos_state % T(i,j,k)   = u(i,j,k,UTEMP)
               eos_state % e(i,j,k)   = u(i,j,k,UEINT) * rhoInv

               eos_state % xn(i,j,k,:)  = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux(i,j,k,:) = u(i,j,k,UFX:UFX+naux-1) * rhoInv
            enddo
         enddo
      enddo

      call eos(eos_input_re, eos_state)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               p(i,j,k,1) = eos_state % p(i,j,k)
            enddo
         enddo
      enddo

      end subroutine ca_derpres

!-----------------------------------------------------------------------

      subroutine ca_dereint1(e,e_l1,e_l2,e_l3,e_h1,e_h2,e_h3,ncomp_e, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use bl_constants_module
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN 

      implicit none

      integer e_l1,e_l2,e_l3,e_h1,e_h2,e_h3,ncomp_e
      integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer lo(3), hi(3), domlo(3), domhi(3)
      double precision e(e_l1:e_h1,e_l2:e_h2,e_l3:e_h3,ncomp_e)
      double precision u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      double precision dx(3), xlo(3), time, dt
      integer bc(3,2,ncomp_u), level, grid_no

      double precision :: rhoInv,ux,uy,uz
      integer          :: i,j,k
      !
      ! Compute internal energy from (rho E).
      !
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               rhoInv = ONE/u(i,j,k,URHO)
               ux = u(i,j,k,UMX)*rhoInv
               uy = u(i,j,k,UMY)*rhoInv
               uz = u(i,j,k,UMZ)*rhoInv
               e(i,j,k,1) = u(i,j,k,UEDEN)*rhoInv-HALF*(ux**2+uy**2+uz**2)
            enddo
         enddo
      enddo

      end subroutine ca_dereint1

!-----------------------------------------------------------------------

      subroutine ca_dereint2(e,e_l1,e_l2,e_l3,e_h1,e_h2,e_h3,ncomp_e, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use meth_params_module, only : URHO, UEINT

      implicit none

      integer e_l1,e_l2,e_l3,e_h1,e_h2,e_h3,ncomp_e
      integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer lo(3), hi(3), domlo(3), domhi(3)
      double precision e(e_l1:e_h1,e_l2:e_h2,e_l3:e_h3,ncomp_e)
      double precision u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      double precision dx(3), xlo(3), time, dt
      integer bc(3,2,ncomp_u), level, grid_no

      integer :: i,j,k
      !
      ! Compute internal energy from (rho e).
      !
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               e(i,j,k,1) = u(i,j,k,UEINT) / u(i,j,k,URHO)
            enddo
         enddo
      enddo

      end subroutine ca_dereint2

!-----------------------------------------------------------------------

      subroutine ca_dersoundspeed(c,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,ncomp_c, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use eos_type_module
      use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,ncomp_c
      integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer lo(3), hi(3), domlo(3), domhi(3)
      double precision c(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,ncomp_c)
      double precision u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      double precision dx(3), xlo(3), time, dt
      integer bc(3,2,ncomp_u), level, grid_no

      double precision :: rhoInv
      integer          :: i,j,k

      type (eos_t_3D) :: eos_state

      eos_state = eos_t_3D(lo,hi)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhoInv = ONE / u(i,j,k,URHO)
               
               eos_state % rho(i,j,k) = u(i,j,k,URHO)
               eos_state % T(i,j,k)   = u(i,j,k,UTEMP)
               eos_state % e(i,j,k)   = u(i,j,k,UEINT) * rhoInv

               eos_state % xn(i,j,k,:) = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux(i,j,k,:) = u(i,j,k,UFX:UFX+naux-1) * rhoInv
            enddo
         enddo
      enddo

      call eos(eos_input_re, eos_state)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j,k,1) = eos_state % cs(i,j,k)
            enddo
         enddo
      enddo

      end subroutine ca_dersoundspeed

!-----------------------------------------------------------------------

      subroutine ca_dermachnumber(mach,mach_l1,mach_l2,mach_l3,mach_h1,mach_h2,mach_h3,ncomp_mach, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use eos_type_module
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer          :: mach_l1,mach_l2,mach_l3,mach_h1,mach_h2,mach_h3,ncomp_mach
      integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      double precision :: mach(mach_l1:mach_h1,mach_l2:mach_h2,mach_l3:mach_h3,ncomp_mach)
      double precision :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      double precision :: dx(3), xlo(3), time, dt
      integer          :: bc(3,2,ncomp_u), level, grid_no

      double precision :: rhoInv,ux,uy,uz
      integer          :: i,j,k

      type (eos_t_3D) :: eos_state

      eos_state = eos_t_3D(lo,hi)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhoInv = ONE / u(i,j,k,URHO)
               
               eos_state % rho(i,j,k) = u(i,j,k,URHO)
               eos_state % T(i,j,k)   = u(i,j,k,UTEMP)
               eos_state % e(i,j,k)   = u(i,j,k,UEINT) * rhoInv

               eos_state % xn(i,j,k,:)  = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux(i,j,k,:) = u(i,j,k,UFX:UFX+naux-1) * rhoInv
            enddo
         enddo
      enddo

      call eos(eos_input_re, eos_state)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               mach(i,j,k,1) = (u(i,j,k,UMX)**2 + u(i,j,k,UMY)**2 + u(i,j,k,UMZ)**2)**0.5 / u(i,j,k,URHO) &
                             / eos_state % cs(i,j,k)
            enddo
         enddo
      enddo

      end subroutine ca_dermachnumber

!-----------------------------------------------------------------------

      subroutine ca_derentropy(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ncomp_s, &
           u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use eos_type_module
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ncomp_s
      integer u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
      integer lo(3), hi(3), domlo(3), domhi(3)
      double precision s(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ncomp_s)
      double precision u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
      double precision dx(3), xlo(3), time, dt
      integer bc(3,2,ncomp_u), level, grid_no

      double precision :: rhoInv
      integer i,j,k

      type (eos_t_3D) :: eos_state

      eos_state = eos_t_3D(lo,hi)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhoInv = ONE / u(i,j,k,URHO)
               
               eos_state % rho(i,j,k) = u(i,j,k,URHO)
               eos_state % T(i,j,k)   = u(i,j,k,UTEMP)
               eos_state % e(i,j,k)   = u(i,j,k,UEINT) * rhoInv

               eos_state % xn(i,j,k,:)  = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux(i,j,k,:) = u(i,j,k,UFX:UFX+naux-1) * rhoInv
            enddo
         enddo
      enddo               

      call eos(eos_input_re, eos_state)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               s(i,j,k,1) = eos_state % s(i,j,k)
            enddo
         enddo
      enddo

      end subroutine ca_derentropy

!-----------------------------------------------------------------------

      subroutine ca_derspec(spec,spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3,nv, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                            domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive the velocity from the momentum.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision spec(spec_l1:spec_h1,spec_l2:spec_h2,spec_l3:spec_h3,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no
 
      integer i,j,k
 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               spec(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
            end do
         end do
      end do
 
      end subroutine ca_derspec

!-----------------------------------------------------------------------

      subroutine ca_derlogden(logden,ld_l1,ld_l2,ld_l3,ld_h1,ld_h2,ld_h3,nd, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      implicit none

      integer          lo(3), hi(3)
      integer           ld_l1, ld_l2, ld_l3, ld_h1, ld_h2, ld_h3,nd
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3), level, grid_no
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision logden( ld_l1: ld_h1, ld_l2: ld_h2, ld_l3: ld_h3,nd)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
 
      integer    i,j,k
 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               logden(i,j,k,1) = dlog10(dat(i,j,k,1))
            end do
         end do
      end do
 
      end subroutine ca_derlogden

!-----------------------------------------------------------------------

      subroutine ca_dermagvort(vort,v_l1,v_l2,v_l3,v_h1,v_h2,v_h3,nv, & 
                               dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc,lo,hi,domlo, &
                               domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will calculate vorticity
      !     

      use bl_constants_module

      implicit none

      integer          lo(3), hi(3)
      integer            v_l1,  v_l2,  v_l3,  v_h1,  v_h2,  v_h3,nv
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3), level, grid_no
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision vort(  v_l1:  v_h1,  v_l2:  v_h2,  v_l3:  v_h3,nv)
      double precision, intent(in) :: dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)

      integer          :: i,j,k
      double precision :: uy,uz,vx,vz,wx,wy,v1,v2,v3
      double precision :: ldat(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,2:4)
      !
      ! Convert momentum to velocity.
      !
      do k = lo(3)-1, hi(3)+1
         do j = lo(2)-1, hi(2)+1
            do i = lo(1)-1, hi(1)+1
               ldat(i,j,k,2) = dat(i,j,k,2) / dat(i,j,k,1)
               ldat(i,j,k,3) = dat(i,j,k,3) / dat(i,j,k,1)
               ldat(i,j,k,4) = dat(i,j,k,4) / dat(i,j,k,1)
            end do
         end do
      end do
      !
      ! Calculate vorticity.
      !
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uy = HALF * (ldat(i,j+1,k,2) - ldat(i,j-1,k,2)) / delta(2)
               uz = HALF * (ldat(i,j,k+1,2) - ldat(i,j,k-1,2)) / delta(3)
               vx = HALF * (ldat(i+1,j,k,3) - ldat(i-1,j,k,3)) / delta(1)
               vz = HALF * (ldat(i,j,k+1,3) - ldat(i,j,k-1,3)) / delta(3)
               wx = HALF * (ldat(i+1,j,k,4) - ldat(i-1,j,k,4)) / delta(1)
               wy = HALF * (ldat(i,j+1,k,4) - ldat(i,j-1,k,4)) / delta(2)
               v1 = wy - vz
               v2 = uz - wx
               v3 = vx - uy
               vort(i,j,k,1) = sqrt(v1*v1 + v2*v2 + v3*v3)
            end do
         end do
      end do

      end subroutine ca_dermagvort

!-----------------------------------------------------------------------

      subroutine ca_derdivu(divu,div_l1,div_l2,div_l3,div_h1,div_h2,div_h3,nd, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                            lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will divergence of velocity.
      !

      use bl_constants_module

      implicit none

      integer          lo(3), hi(3)
      integer          div_l1,div_l2,div_l3,div_h1,div_h2,div_h3,nd
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision divu(div_l1:div_h1,div_l2:div_h2,div_l3:div_h3,nd)
      double precision  dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer          :: i,j,k
      double precision :: ulo,uhi,vlo,vhi,wlo,whi

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uhi = dat(i+1,j,k,2) / dat(i+1,j,k,1)
               ulo = dat(i-1,j,k,2) / dat(i-1,j,k,1)
               vhi = dat(i,j+1,k,3) / dat(i,j+1,k,1)
               vlo = dat(i,j-1,k,3) / dat(i,j-1,k,1)
               whi = dat(i,j,k+1,4) / dat(i,j,k+1,1)
               wlo = dat(i,j,k-1,4) / dat(i,j,k-1,1)
               divu(i,j,k,1) = HALF * ( (uhi-ulo) / delta(1) + &
                                         (vhi-vlo) / delta(2) + &
                                         (whi-wlo) / delta(3) )
            end do
         end do
      end do

      end subroutine ca_derdivu

!-----------------------------------------------------------------------

      subroutine ca_derkineng(kineng,ken_l1,ken_l2,ken_l3,ken_h1,ken_h2,ken_h3,nk, &
                              dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive kinetic energy = 1/2 rho (u^2 + v^2 + w^2)
      !

      use bl_constants_module

      implicit none

      integer          lo(3), hi(3)
      integer          ken_l1,ken_l2,ken_l3,ken_h1,ken_h2,ken_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision kineng(ken_l1:ken_h1,ken_l2:ken_h2,ken_l3:ken_h3,nk)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      integer i,j,k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               kineng(i,j,k,1) = HALF / dat(i,j,k,1) * ( dat(i,j,k,2)**2 + &
                                                          dat(i,j,k,3)**2 + &
                                                          dat(i,j,k,4)**2 )
            end do
         end do
      end do

      end subroutine ca_derkineng

      subroutine ca_dernull(kineng,ken_l1,ken_l2,ken_l3,ken_h1,ken_h2,ken_h3,nk, &
                            dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
                             lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine is used by particle_count.  Yes it does nothing.
      !
      implicit none

      integer          lo(3), hi(3)
      integer          ken_l1,ken_l2,ken_l3,ken_h1,ken_h2,ken_h3,nk
      integer          dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
      integer          domlo(3), domhi(3)
      integer          bc(3,2,nc)
      double precision delta(3), xlo(3), time, dt
      double precision kineng(ken_l1:ken_h1,ken_l2:ken_h2,ken_l3:ken_h3,nk)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
      integer    level, grid_no

      end subroutine ca_dernull

