!-----------------------------------------------------------------------

      subroutine ca_derstate(state,state_l1,state_l2,state_h1,state_h2,nv, &
                             dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo, &
                             domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     The incoming   "dat" vector contains (rho,T,(rho X)_1)
!     The outgoing "state" vector contains (rho,T,X_1)
!
      implicit none 

      integer          lo(2), hi(2)
      integer          state_l1,state_l2,state_h1,state_h2,nv
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision state(state_l1:state_h1,state_l2:state_h2,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no
 
      integer    i,j

      if (nv .ne. 3) then
          print *,'... confusion in derstate ... nv should be 3 but is ',nv
          call bl_error('Error:: Derive_2d.f90 :: ca_derstate')
      end if

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         ! Density
         state(i,j,1) = dat(i,j,1)
         ! Temperature
         state(i,j,2) = dat(i,j,2)
         ! (rho X)_1 --> X_1
         state(i,j,3) = dat(i,j,3) / dat(i,j,1)
      end do
      end do
 
      end subroutine ca_derstate

!-----------------------------------------------------------------------

      subroutine ca_dervel(vel,vel_l1,vel_l2,vel_h1,vel_h2,nv,&
                           dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo,&
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive the velocity from the momentum.
!
      implicit none 

      integer          lo(2), hi(2)
      integer          vel_l1,vel_l2,vel_h1,vel_h2,nv
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision vel(vel_l1:vel_h1,vel_l2:vel_h2,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no
 
      integer    i,j
 
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            vel(i,j,1) = dat(i,j,2) / dat(i,j,1)
         end do
      end do
 
      end subroutine ca_dervel

!-----------------------------------------------------------------------

      subroutine ca_dermagvel(magvel,vel_l1,vel_l2,vel_h1,vel_h2,nv, &
                              dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo, &
                              domhi,dx,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive magnitude of velocity.
!
      implicit none 

      integer          lo(2), hi(2)
      integer          vel_l1,vel_l2,vel_h1,vel_h2,nv
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision dx(2), xlo(2), time, dt
      double precision magvel(vel_l1:vel_h1,vel_l2:vel_h2,nv)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no

      integer    i,j

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            magvel(i,j,1) = sqrt( (dat(i,j,2) / dat(i,j,1))**2 + (dat(i,j,3) / dat(i,j,1))**2 )
         end do
      end do

      end subroutine ca_dermagvel

!-----------------------------------------------------------------------

      subroutine ca_derradialvel(radvel,vel_l1,vel_l2,vel_h1,vel_h2,nv, &
                                 dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo, &
                                 domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive the radial velocity.
!
      use probdata_module, only : center
      use bl_constants_module

      implicit none 

      integer          lo(2), hi(2)
      integer          vel_l1,vel_l2,vel_h1,vel_h2,nv
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision radvel(vel_l1:vel_h1,vel_l2:vel_h2,nv)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no

      integer          :: i,j
      double precision :: x,y,r

      do j = lo(2), hi(2)
         y = xlo(2) + (dble(j-lo(2))+HALF) * delta(2) - center(2)
         do i = lo(1), hi(1)
            x = xlo(1) + (dble(i-lo(1))+HALF) * delta(1) - center(1)
            r = sqrt(x*x+y*y)
            radvel(i,j,1) = dat(i,j,2)/dat(i,j,1) * (x/r) + dat(i,j,3)/dat(i,j,1) * (y/r)
         end do
      end do

      end subroutine ca_derradialvel

!-----------------------------------------------------------------------

      subroutine ca_dermagmom(magmom,mom_l1,mom_l2,mom_h1,mom_h2,nv, &
                              dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive magnitude of momentum
!
      implicit none 

      integer          lo(2), hi(2)
      integer          mom_l1,mom_l2,mom_h1,mom_h2,nv
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision magmom(mom_l1:mom_h1,mom_l2:mom_h2,nv)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no

      integer    i,j

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            magmom(i,j,1) = sqrt( dat(i,j,1)**2 + dat(i,j,2)**2 )
         end do
      end do

      end subroutine ca_dermagmom

!-----------------------------------------------------------------------

      subroutine ca_derpres(p,p_l1,p_l2,p_h1,p_h2,ncomp_p, &
           u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy, small_temp
      use bl_constants_module

      implicit none

      integer          :: p_l1,p_l2,p_h1,p_h2,ncomp_p
      integer          :: u_l1,u_l2,u_h1,u_h2,ncomp_u
      integer          :: lo(2), hi(2), domlo(2), domhi(2)
      double precision :: p(p_l1:p_h1,p_l2:p_h2,ncomp_p)
      double precision :: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
      double precision :: dx(2), xlo(2), time, dt
      integer          :: bc(2,2,ncomp_u), level, grid_no

      double precision :: rhoInv
      integer          :: i,j

      type (eos_t) :: eos_state

!     Compute pressure from the EOS
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rhoInv = ONE / u(i,j,URHO)

            eos_state % xn  = u(i,j,UFS:UFS+nspec-1) * rhoInv
            eos_state % aux = u(i,j,UFX:UFX+naux-1) * rhoInv
            eos_state % rho = u(i,j,URHO)
            eos_state % T   = u(i,j,UTEMP)
            eos_state % e   = u(i,j,UEINT) * rhoInv

            ! Protect against negative internal energy
            if (allow_negative_energy .eq. 0 .and. eos_state % e .le. ZERO) then
               call eos(eos_input_rt, eos_state)
            else
               call eos(eos_input_re, eos_state)
            end if
 
            p(i,j,1) = eos_state % p
         enddo
      enddo

      end subroutine ca_derpres

!-----------------------------------------------------------------------

      subroutine ca_dereint1(e,e_l1,e_l2,e_h1,e_h2,ncomp_e, &
           u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use meth_params_module, only : URHO, UMX, UMY, UEDEN 
      use bl_constants_module

      implicit none

      integer e_l1,e_l2,e_h1,e_h2,ncomp_e
      integer u_l1,u_l2,u_h1,u_h2,ncomp_u
      integer lo(2), hi(2), domlo(2), domhi(2)
      double precision e(e_l1:e_h1,e_l2:e_h2,ncomp_e)
      double precision u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
      double precision dx(2), xlo(2), time, dt
      integer bc(2,2,ncomp_u), level, grid_no

      double precision :: rhoInv,ux,uy
      integer          :: i,j

!     Compute internal energy from (rho E)
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rhoInv = ONE/u(i,j,URHO)
            ux = u(i,j,UMX)*rhoInv
            uy = u(i,j,UMY)*rhoInv
            e(i,j,1) = u(i,j,UEDEN)*rhoInv-HALF*(ux**2+uy**2)
         enddo
      enddo

      end subroutine ca_dereint1

!-----------------------------------------------------------------------

      subroutine ca_dereint2(e,e_l1,e_l2,e_h1,e_h2,ncomp_e, &
           u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use meth_params_module, only : URHO, UEINT

      implicit none

      integer e_l1,e_l2,e_h1,e_h2,ncomp_e
      integer u_l1,u_l2,u_h1,u_h2,ncomp_u
      integer lo(2), hi(2), domlo(2), domhi(2)
      double precision e(e_l1:e_h1,e_l2:e_h2,ncomp_e)
      double precision u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
      double precision dx(2), xlo(2), time, dt
      integer bc(2,2,ncomp_u), level, grid_no

      integer          :: i,j

!     Compute internal energy from (rho e)
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            e(i,j,1) = u(i,j,UEINT) / u(i,j,URHO)
         enddo
      enddo

      end subroutine ca_dereint2

!-----------------------------------------------------------------------

      subroutine ca_dersoundspeed(c,c_l1,c_l2,c_h1,c_h2,ncomp_c, &
           u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer c_l1,c_l2,c_h1,c_h2,ncomp_c
      integer u_l1,u_l2,u_h1,u_h2,ncomp_u
      integer lo(2), hi(2), domlo(2), domhi(2)
      double precision c(c_l1:c_h1,c_l2:c_h2,ncomp_c)
      double precision u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
      double precision dx(2), xlo(2), time, dt
      integer bc(2,2,ncomp_u), level, grid_no

      double precision :: rhoInv
      integer          :: i,j

      type (eos_t) :: eos_state

!     Compute soundspeed from the EOS
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rhoInv = ONE / u(i,j,URHO)

            eos_state % e  = u(i,j,UEINT) * rhoInv

            ! Protect against negative internal energy
            if (allow_negative_energy .eq. 1 .or. eos_state % e .gt. ZERO) then

               eos_state % xn  = u(i,j,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux = u(i,j,UFX:UFX+naux-1) * rhoInv
               eos_state % rho = u(i,j,URHO)
               eos_state % T   = u(i,j,UTEMP)

               call eos(eos_input_re, eos_state)

               c(i,j,1) = eos_state % cs

            else
               
               c(i,j,1) = zero

            endif
         enddo
      enddo

      end subroutine ca_dersoundspeed

!-----------------------------------------------------------------------

      subroutine ca_dermachnumber(mach,mach_l1,mach_l2,mach_h1,mach_h2,ncomp_mach, &
           u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer          :: mach_l1,mach_l2,mach_h1,mach_h2,ncomp_mach
      integer          :: u_l1,u_l2,u_h1,u_h2,ncomp_u
      integer          :: lo(2), hi(2), domlo(2), domhi(2)
      double precision :: mach(mach_l1:mach_h1,mach_l2:mach_h2,ncomp_mach)
      double precision :: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
      double precision :: dx(2), xlo(2), time, dt
      integer          :: bc(2,2,ncomp_u), level, grid_no

      double precision :: rhoInv,ux,uy
      integer          :: i,j

      type (eos_t) :: eos_state

!     Compute Mach number of the flow
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rhoInv = ONE / u(i,j,URHO)

            eos_state % e = u(i,j,UEINT) * rhoInv

            if (allow_negative_energy .eq. 1 .or. eos_state % e .gt. ZERO) then
               ux = u(i,j,UMX  ) * rhoInv
               uy = u(i,j,UMY  ) * rhoInv
               eos_state % rho = u(i,j,URHO)
               eos_state % T   = u(i,j,UTEMP)
               eos_state % xn  = u(i,j,UFS:UFS+nspec-1) * rhoInv 
               eos_state % aux = u(i,j,UFX:UFX+naux-1) * rhoInv

               call eos(eos_input_re, eos_state)
               mach(i,j,1) = sqrt(ux**2 + uy**2) / eos_state % cs
            else
               mach(i,j,1) = zero
            end if

         enddo
      enddo

      end subroutine ca_dermachnumber


!-----------------------------------------------------------------------

      subroutine ca_derentropy(s,s_l1,s_l2,s_h1,s_h2,ncomp_s, &
           u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer s_l1,s_l2,s_h1,s_h2,ncomp_s
      integer u_l1,u_l2,u_h1,u_h2,ncomp_u
      integer lo(2), hi(2), domlo(2), domhi(2)
      double precision s(s_l1:s_h1,s_l2:s_h2,ncomp_s)
      double precision u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
      double precision dx(2), xlo(2), time, dt
      integer bc(2,2,ncomp_u), level, grid_no

      double precision :: rhoInv
      integer          :: i,j

      type (eos_t) :: eos_state

!     Compute entropy from the EOS
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rhoInv = ONE / u(i,j,URHO)

            eos_state % e = u(i,j,UEINT) * rhoInv

            if (allow_negative_energy .eq. 1 .or. eos_state % e .gt. ZERO) then
               eos_state % rho = u(i,j,URHO)
               eos_state % T   = u(i,j,UTEMP)
               eos_state % xn  = u(i,j,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux = u(i,j,UFX:UFX+naux-1) * rhoInv

               call eos(eos_input_re, eos_state)
               s(i,j,1) = eos_state % s
            else
               s(i,j,1) = zero
            endif
         enddo
      enddo

      end subroutine ca_derentropy

!-----------------------------------------------------------------------

      subroutine ca_derspec(spec,spec_l1,spec_l2,spec_h1,spec_h2,nv,&
                            dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo,&
                            domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive the X from the (rho X)
!
      implicit none 

      integer          lo(2), hi(2)
      integer          spec_l1,spec_l2,spec_h1,spec_h2,nv
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision spec(spec_l1:spec_h1,spec_l2:spec_h2,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no
 
      integer    i,j
 
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            spec(i,j,1) = dat(i,j,2) / dat(i,j,1)
         end do
      end do
 
      end subroutine ca_derspec

!-----------------------------------------------------------------------

      subroutine ca_derlogden(logden,ld_l1,ld_l2,ld_h1,ld_h2,nd, &
                              dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,&
                              domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive magnitude of velocity.
!
      implicit none 

      integer          lo(2), hi(2)
      integer           ld_l1, ld_l2, ld_h1, ld_h2,nd
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision logden( ld_l1: ld_h1, ld_l2: ld_h2,nd)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no

      integer    i,j

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            logden(i,j,1) = dlog10(dat(i,j,1))
         end do
      end do

      end subroutine ca_derlogden

!-----------------------------------------------------------------------

      subroutine ca_dermaggrav(maggrav,grav_l1,grav_l2,grav_h1,grav_h2,ng, &
                               dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo, &
                               domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive magnitude of the gravity vector.
!
      implicit none 

      integer          lo(2), hi(2)
      integer          grav_l1,grav_l2,grav_h1,grav_h2,ng
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision maggrav(grav_l1:grav_h1,grav_l2:grav_h2,ng)
      double precision     dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no

      integer    i,j

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            maggrav(i,j,1) = sqrt( dat(i,j,1)**2  + dat(i,j,2)**2 )
         end do
      end do

      end subroutine ca_dermaggrav

!-----------------------------------------------------------------------

      subroutine ca_dermagvort(vort,v_l1,v_l2,v_h1,v_h2,nv, &
                               dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo, &
                               domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will calculate vorticity
!

      use bl_constants_module

      implicit none

      integer          lo(2), hi(2)
      integer            v_l1,  v_l2,  v_h1,  v_h2,nv
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision vort(  v_l1:  v_h1,  v_l2:  v_h2,nv)
      double precision, intent(in) :: dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no

      integer          :: i,j
      double precision :: vx,uy
      double precision :: ldat(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2:3)

      ! Convert momentum to velocity
      do j = lo(2)-1, hi(2)+1
      do i = lo(1)-1, hi(1)+1
        ldat(i,j,2) = dat(i,j,2) / dat(i,j,1)
        ldat(i,j,3) = dat(i,j,3) / dat(i,j,1)
      end do
      end do

      ! Calculate vorticity
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         vx = HALF * (ldat(i+1,j,3) - ldat(i-1,j,3)) / delta(1) 
         uy = HALF * (ldat(i,j+1,2) - ldat(i,j-1,2)) / delta(2) 
         vort(i,j,1) = abs(vx - uy)
      end do
      end do

      end subroutine ca_dermagvort

!-----------------------------------------------------------------------

      subroutine ca_derdivu(divu,div_l1,div_l2,div_h1,div_h2,nd, &
                            dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo, &
                            domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will divergence of velocity.
!

      use bl_constants_module

      implicit none

      integer          lo(2), hi(2)
      integer          div_l1,div_l2,div_h1,div_h2,nd
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision divu(div_l1:div_h1,div_l2:div_h2,nd)
      double precision  dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no

      integer          :: i,j
      double precision :: ulo,uhi,vlo,vhi

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         uhi = dat(i+1,j,2) / dat(i+1,j,1)
         ulo = dat(i-1,j,2) / dat(i-1,j,1)
         vhi = dat(i,j+1,3) / dat(i,j+1,1)
         vlo = dat(i,j-1,3) / dat(i,j-1,1)
         divu(i,j,1) = HALF * ((uhi-ulo)/delta(1) + (vhi-vlo)/delta(2))
      end do
      end do

      end subroutine ca_derdivu

!-----------------------------------------------------------------------

      subroutine ca_derkineng(kineng,ken_l1,ken_l2,ken_h1,ken_h2,nk, &
                              dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo, &
                              domhi,dx,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive kinetic energy = 1/2 rho (u^2 + v^2)
!

      use bl_constants_module

      implicit none

      integer          lo(2), hi(2)
      integer          ken_l1,ken_l2,ken_h1,ken_h2,nk
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision dx(2), xlo(2), time, dt
      double precision kineng(ken_l1:ken_h1,ken_l2:ken_h2,nk)
      double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no

      integer    i,j

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            kineng(i,j,1) = HALF / dat(i,j,1) * (dat(i,j,2)**2 + dat(i,j,3)**2)
         end do
      end do

      end subroutine ca_derkineng

      subroutine ca_dernull(kineng,ken_l1,ken_l2,ken_h1,ken_h2,nk, &
                            dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo, &
                            domhi,dx,xlo,time,dt,bc,level,grid_no)
        !
        ! This is a derived routine that does nothing.
        !
        implicit none

        integer          lo(2), hi(2)
        integer          ken_l1,ken_l2,ken_h1,ken_h2,nk
        integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
        integer          domlo(2), domhi(2)
        integer          bc(2,2,nc)
        double precision dx(2), xlo(2), time, dt
        double precision kineng(ken_l1:ken_h1,ken_l2:ken_h2,nk)
        double precision    dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
        integer    level, grid_no

      end subroutine ca_dernull

