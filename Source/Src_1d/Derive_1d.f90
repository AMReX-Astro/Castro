
!-----------------------------------------------------------------------

      subroutine ca_derlapvar(var,var_l1,var_h1,nv, &
                              dat,dat_l1,dat_h1,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive the weighted-Laplacian of the variable for
!       the purposes of error estimation
!
      implicit none

      integer          lo(1), hi(1)
      integer          var_l1,var_h1,nv
      integer          dat_l1,dat_h1,nc
      integer          domlo(1), domhi(1)
      integer          bc(1,2,nc)
      double precision delta(1), xlo(1), time, dt
      double precision var(var_l1:var_h1,nv)
      double precision dat(dat_l1:dat_h1,nc)
      integer    level, grid_no
 
      ! added for my routine
      double precision ::  delu(var_l1:var_h1)
      double precision :: delua(var_l1:var_h1)
      double precision :: delu2, delu3, delu4
      double precision :: num, denom
      integer          :: i

      ! This value is taken from FLASH
      double precision, parameter:: epsil=0.02

      ! adapted from ref_marking.f90 in FLASH2.5

      ! d/dx
      do i=lo(1)-1,hi(1)+1
         delu(i)  =     dat(i+1,1)  -     dat(i-1,1) 
         delua(i) = abs(dat(i+1,1)) + abs(dat(i-1,1))
      end do

      ! d/dxdx
      do i = lo(1),hi(1)

         delu2 =     delu(i+1)  -     delu(i-1)
         delu3 = abs(delu(i+1)) + abs(delu(i-1))
         delu4 =    delua(i+1)  +    delua(i-1)

         ! compute the error
         num   = abs(delu2)
         denom = abs(delu3 + (epsil*delu4+1.d-99))

         var(i,1) = num/denom

      end do
 
      end subroutine ca_derlapvar

!-----------------------------------------------------------------------

      subroutine ca_derstate(state,state_l1,state_h1,nv, &
                             dat,dat_l1,dat_h1,nc,lo,hi,domlo, &
                             domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     The incoming   "dat" vector contains (rho,T,(rho X)_1)
!     The outgoing "state" vector contains (rho,T,X_1)
!
      implicit none

      integer          lo(1), hi(1)
      integer          state_l1,state_h1,nv
      integer          dat_l1,dat_h1,nc
      integer          domlo(1), domhi(1)
      integer          bc(1,2,nc)
      double precision delta(1), xlo(1), time, dt
      double precision state(state_l1:state_h1,nv)
      double precision dat(dat_l1:dat_h1,nc)
      integer    level, grid_no

      integer    i

      if (nv .ne. 3) then
          print *,'... confusion in derstate ... nv should be 3 but is ',nv
          call bl_error('Error:: Derive_1d.f90 :: ca_derstate')
      end if

      do i = lo(1), hi(1)
         ! Density
         state(i,1) = dat(i,1)
         ! Temperature
         state(i,2) = dat(i,2)
         ! (rho X)_1 --> X_1
         state(i,3) = dat(i,3) / dat(i,1)
      end do

      end subroutine ca_derstate

!-----------------------------------------------------------------------

      subroutine ca_dervel(vel,vel_l1,vel_h1,nv, &
                           dat,dat_l1,dat_h1,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive the velocity from the momentum.
!
      implicit none

      integer          lo(1), hi(1)
      integer          vel_l1,vel_h1,nv
      integer          dat_l1,dat_h1,nc
      integer          domlo(1), domhi(1)
      integer          bc(1,2,nc)
      double precision delta(1), xlo(1), time, dt
      double precision vel(vel_l1:vel_h1,nv)
      double precision dat(dat_l1:dat_h1,nc)
      integer    level, grid_no
 
      integer    i
 
      do i = lo(1), hi(1)
         vel(i,1) = dat(i,2) / dat(i,1)
      end do
 
      end subroutine ca_dervel

!-----------------------------------------------------------------------

      subroutine ca_dermagvel(magvel,vel_l1,vel_h1,nv, &
                              dat,dat_l1,dat_h1,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will magnitude of velocity.
!
      implicit none

      integer          lo(1), hi(1)
      integer          vel_l1,vel_h1,nv
      integer          dat_l1,dat_h1,nc
      integer          domlo(1), domhi(1)
      integer          bc(1,2,nc)
      double precision delta(1), xlo(1), time, dt
      double precision magvel(vel_l1:vel_h1,nv)
      double precision    dat(dat_l1:dat_h1,nc)
      integer    level, grid_no

      integer    i

      do i = lo(1), hi(1)
         magvel(i,1) = abs(dat(i,2) / dat(i,1))
      end do

      end subroutine ca_dermagvel

!-----------------------------------------------------------------------

      subroutine ca_dermagmom(magmom,mom_l1,mom_h1,nv, &
                              dat,dat_l1,dat_h1,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will magnitude of momentum.
!
      implicit none

      integer          lo(1), hi(1)
      integer          mom_l1,mom_h1,nv
      integer          dat_l1,dat_h1,nc
      integer          domlo(1), domhi(1)
      integer          bc(1,2,nc)
      double precision delta(1), xlo(1), time, dt
      double precision magmom(mom_l1:mom_h1,nv)
      double precision    dat(dat_l1:dat_h1,nc)
      integer    level, grid_no

      integer    i

      do i = lo(1), hi(1)
         magmom(i,1) = abs(dat(i,1))
      end do

      end subroutine ca_dermagmom

!-----------------------------------------------------------------------

      subroutine ca_derpres(p,p_l1,p_h1,ncomp_p, &
           u,u_l1,u_h1,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : URHO, UMX, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer p_l1,p_h1,ncomp_p
      integer u_l1,u_h1,ncomp_u
      integer lo(1), hi(1), domlo(1), domhi(1)
      double precision p(p_l1:p_h1,ncomp_p)
      double precision u(u_l1:u_h1,ncomp_u)
      double precision dx(1), xlo(1), time, dt
      integer bc(1,2,ncomp_u), level, grid_no

      double precision :: rhoInv
      integer          :: i

      type (eos_t) :: eos_state

!     Compute pressure from the EOS
      do i = lo(1),hi(1)

         rhoInv = ONE / u(i,URHO)
         eos_state % rho = u(i,URHO)
         eos_state % e   = u(i,UEINT) * rhoInv
         eos_state % T   = u(i,UTEMP)
         eos_state % xn  = u(i,UFS:UFS+nspec-1) * rhoInv
         eos_state % aux = u(i,UFX:UFX+naux-1) * rhoInv

         ! Protect against negative internal energy
         if (allow_negative_energy .eq. 0 .and. eos_state % e .le. ZERO) then
            call eos(eos_input_rt, eos_state)
         else
            call eos(eos_input_re, eos_state)
         end if

         p(i,1) = eos_state % p

      enddo

      end subroutine ca_derpres

!-----------------------------------------------------------------------

      subroutine ca_dereint1(e,e_l1,e_h1,ncomp_e, &
           u,u_l1,u_h1,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use meth_params_module, only : URHO, UMX, UEDEN 

      implicit none

      integer e_l1,e_h1,ncomp_e
      integer u_l1,u_h1,ncomp_u
      integer lo(1), hi(1), domlo(1), domhi(1)
      double precision e(e_l1:e_h1,ncomp_e)
      double precision u(u_l1:u_h1,ncomp_u)
      double precision dx(1), xlo(1), time, dt
      integer bc(1,2,ncomp_u), level, grid_no

      double precision :: rhoInv,ux
      integer          :: i

!     Compute internal energy from (rho E)
      do i = lo(1),hi(1)
         rhoInv = 1.d0/u(i,URHO)
         ux = u(i,UMX)*rhoInv
         e(i,1) = u(i,UEDEN)*rhoInv-0.5d0*(ux**2)
      enddo

      end subroutine ca_dereint1

!-----------------------------------------------------------------------

      subroutine ca_dereint2(e,e_l1,e_h1,ncomp_e, &
           u,u_l1,u_h1,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use meth_params_module, only : URHO, UEINT

      implicit none

      integer e_l1,e_h1,ncomp_e
      integer u_l1,u_h1,ncomp_u
      integer lo(1), hi(1), domlo(1), domhi(1)
      double precision e(e_l1:e_h1,ncomp_e)
      double precision u(u_l1:u_h1,ncomp_u)
      double precision dx(1), xlo(1), time, dt
      integer bc(1,2,ncomp_u), level, grid_no

      integer          :: i

!     Compute internal energy from (rho e)
      do i = lo(1),hi(1)
         e(i,1) = u(i,UEINT) / u(i,URHO)
      enddo

      end subroutine ca_dereint2

!-----------------------------------------------------------------------

      subroutine ca_dersoundspeed(c,c_l1,c_h1,ncomp_c, &
           u,u_l1,u_h1,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : URHO, UMX, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer c_l1,c_h1,ncomp_c
      integer u_l1,u_h1,ncomp_u
      integer lo(1), hi(1), domlo(1), domhi(1)
      double precision c(c_l1:c_h1,ncomp_c)
      double precision u(u_l1:u_h1,ncomp_u)
      double precision dx(1), xlo(1), time, dt
      integer bc(1,2,ncomp_u), level, grid_no

      double precision :: rhoInv
      integer          :: i

      type (eos_t) :: eos_state

      c = ZERO

!     Compute soundspeed from the EOS
      do i = lo(1),hi(1)

         rhoInv = ONE / u(i,URHO)
         eos_state % e = u(i,UEINT) * rhoInv

         ! Protect against negative internal energy

         if (allow_negative_energy .eq. 1 .or. eos_state % e .gt. ZERO) then

            eos_state % xn  = u(i,UFS:UFS+nspec-1) * rhoInv
            eos_state % aux = u(i,UFX:UFX+naux-1) * rhoInv
            eos_state % rho = u(i,URHO)
            eos_state % T   = u(i,UTEMP)

            call eos(eos_input_re, eos_state)
 
            c(i,1) = eos_state % cs 

         end if

      enddo

      end subroutine ca_dersoundspeed

!-----------------------------------------------------------------------

      subroutine ca_dermachnumber(mach,mach_l1,mach_h1,ncomp_mach, &
           u,u_l1,u_h1,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : URHO, UMX, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer          :: mach_l1,mach_h1,ncomp_mach
      integer          :: u_l1,u_h1,ncomp_u
      integer          :: lo(1), hi(1), domlo(1), domhi(1)
      double precision :: mach(mach_l1:mach_h1,ncomp_mach)
      double precision :: u(u_l1:u_h1,ncomp_u)
      double precision :: dx(1), xlo(1), time, dt
      integer          :: bc(1,2,ncomp_u), level, grid_no

      double precision :: rhoInv,ux
      integer          :: i

      type (eos_t) :: eos_state

      mach = ZERO

!     Compute Mach number of the flow
      do i = lo(1),hi(1)
         rhoInv = ONE / u(i,URHO)
         eos_state % e = u(i,UEINT) * rhoInv

         if (allow_negative_energy .eq. 1 .or. eos_state % e .gt. ZERO) then
 
            ux = u(i,UMX) * rhoInv
            eos_state % rho = u(i,URHO)
            eos_state % T   = u(i,UTEMP)
            eos_state % xn  = u(i,UFS:UFS+nspec-1) * rhoInv
            eos_state % aux = u(i,UFX:UFX+naux-1) * rhoInv

            call eos(eos_input_re, eos_state)

            mach(i,1) = abs(ux) / eos_state % cs
         end if 

      enddo

      end subroutine ca_dermachnumber


!-----------------------------------------------------------------------

      subroutine ca_deruplusc(c,c_l1,c_h1,ncomp_c, &
           u,u_l1,u_h1,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : URHO, UMX, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer c_l1,c_h1,ncomp_c
      integer u_l1,u_h1,ncomp_u
      integer lo(1), hi(1), domlo(1), domhi(1)
      double precision c(c_l1:c_h1,ncomp_c)
      double precision u(u_l1:u_h1,ncomp_u)
      double precision dx(1), xlo(1), time, dt
      integer bc(1,2,ncomp_u), level, grid_no

      double precision :: rhoInv,ux
      integer          :: i

      type (eos_t) :: eos_state

      c = ZERO

!     Compute soundspeed from the EOS
      do i = lo(1),hi(1)

         rhoInv = ONE / u(i,URHO)
         ux = u(i,UMX) * rhoInv
         eos_state % e  = u(i,UEINT) * rhoInv

         c(i,1) = ux

         if (allow_negative_energy .eq. 1 .or. eos_state % e .gt. ZERO) then
            eos_state % T   = u(i,UTEMP)
            eos_state % rho = u(i,URHO)
            eos_state % xn  = u(i,UFS:UFS+nspec-1) * rhoInv
            eos_state % aux = u(i,UFX:UFX+naux-1) * rhoInv

            call eos(eos_input_re, eos_state)
            c(i,1) = c(i,1) + eos_state % cs
         end if

      enddo

      end subroutine ca_deruplusc

!-----------------------------------------------------------------------

      subroutine ca_deruminusc(c,c_l1,c_h1,ncomp_c, &
           u,u_l1,u_h1,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : URHO, UMX, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer c_l1,c_h1,ncomp_c
      integer u_l1,u_h1,ncomp_u
      integer lo(1), hi(1), domlo(1), domhi(1)
      double precision c(c_l1:c_h1,ncomp_c)
      double precision u(u_l1:u_h1,ncomp_u)
      double precision dx(1), xlo(1), time, dt
      integer bc(1,2,ncomp_u), level, grid_no

      double precision :: rhoInv,ux
      integer          :: i

      type (eos_t) :: eos_state

      c = ZERO

!     Compute soundspeed from the EOS
      do i = lo(1),hi(1)

         rhoInv = ONE / u(i,URHO)
         ux = u(i,UMX) * rhoInv
         eos_state % e  = u(i,UEINT) * rhoInv

         c(i,1) = ux

         if (allow_negative_energy .eq. 1 .or. eos_state % e .gt. ZERO) then
            eos_state % T   = u(i,UTEMP)
            eos_state % rho = u(i,URHO)
            eos_state % xn  = u(i,UFS:UFS+nspec-1) * rhoInv
            eos_state % aux = u(i,UFX:UFX+naux-1) * rhoInv

            call eos(eos_input_re, eos_state)
            c(i,1) = c(i,1) - eos_state % cs
         end if

      enddo

      end subroutine ca_deruminusc

!-----------------------------------------------------------------------

      subroutine ca_derentropy(s,s_l1,s_h1,ncomp_s, &
           u,u_l1,u_h1,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : URHO, UMX, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer s_l1,s_h1,ncomp_s
      integer u_l1,u_h1,ncomp_u
      integer lo(1), hi(1), domlo(1), domhi(1)
      double precision s(s_l1:s_h1,ncomp_s)
      double precision u(u_l1:u_h1,ncomp_u)
      double precision dx(1), xlo(1), time, dt
      integer bc(1,2,ncomp_u), level, grid_no

      double precision :: rhoInv
      integer          :: i

      type (eos_t) :: eos_state

      s = ZERO

!     Compute entropy from the EOS
      do i = lo(1),hi(1)

         rhoInv = ONE / u(i,URHO)
         eos_state % e  = u(i,UEINT) * rhoInv

         if (allow_negative_energy .eq. 1 .or. eos_state % e .gt. ZERO) then
            eos_state % rho = u(i,URHO)
            eos_state % T   = u(i,UTEMP)
            eos_state % xn  = u(i,UFS:UFS+nspec-1) * rhoInv
            eos_state % aux = u(i,UFX:UFX+naux-1) * rhoInv

            call eos(eos_input_re, eos_state)
            s(i,1) = eos_state % s
         end if

      enddo

      end subroutine ca_derentropy

!-----------------------------------------------------------------------

      subroutine ca_derspec(spec,spec_l1,spec_h1,nv, &
                            dat,dat_l1,dat_h1,nc,lo,hi,domlo, &
                            domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive the velocity from the momentum.
!
      implicit none

      integer          lo(1), hi(1)
      integer          spec_l1,spec_h1,nv
      integer          dat_l1,dat_h1,nc
      integer          domlo(1), domhi(1)
      integer          bc(1,2,nc)
      double precision delta(1), xlo(1), time, dt
      double precision spec(spec_l1:spec_h1,nv)
      double precision dat(dat_l1:dat_h1,nc)
      integer    level, grid_no

      integer    i

      do i = lo(1), hi(1)
         spec(i,1) = dat(i,2) / dat(i,1)
      end do

      end subroutine ca_derspec

!-----------------------------------------------------------------------

      subroutine ca_derlogden(logden,ld_l1,ld_h1,nd, &
                              dat,dat_l1,dat_h1,nc,lo,hi,&
                              domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive magnitude of velocity.
!
      implicit none

      integer          lo(1), hi(1)
      integer           ld_l1, ld_h1, nd
      integer          dat_l1,dat_h1,nc
      integer          domlo(1), domhi(1)
      integer          bc(1,2,nc)
      double precision delta(1), xlo(1), time, dt
      double precision logden( ld_l1: ld_h1,nd)
      double precision    dat(dat_l1:dat_h1,nc)
      integer    level, grid_no

      integer    i

      do i = lo(1), hi(1)
         logden(i,1) = dlog10(dat(i,1))
      end do

      end subroutine ca_derlogden

!-----------------------------------------------------------------------

      subroutine ca_derdivu(divu,div_l1,div_h1,nv, &
                           dat,dat_l1,dat_h1,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will divergence of velocity.
!
      implicit none

      integer          lo(1), hi(1)
      integer          div_l1,div_h1,nv
      integer          dat_l1,dat_h1,nc
      integer          domlo(1), domhi(1)
      integer          bc(1,2,nc)
      double precision delta(1), xlo(1), time, dt
      double precision divu(div_l1:div_h1,nv)
      double precision  dat(dat_l1:dat_h1,nc)
      integer    level, grid_no

      integer          :: i
      double precision :: vlo,vhi

      do i = lo(1), hi(1)
         vhi = dat(i+1,2) / dat(i+1,1)
         vlo = dat(i-1,2) / dat(i-1,1)
         divu(i,1) = 0.5d0 * (vhi - vlo) / delta(1)
      end do

      end subroutine ca_derdivu

!-----------------------------------------------------------------------

      subroutine ca_derkineng(kineng,ken_l1,ken_h1,nk, &
                              dat,dat_l1,dat_h1,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive kinetic energy = 1/2 rho (u^2 + v^2)
!
      implicit none

      integer          lo(1), hi(1)
      integer          ken_l1,ken_h1,nk
      integer          dat_l1,dat_h1,nc
      integer          domlo(1), domhi(1)
      integer          bc(1,2,nc)
      double precision delta(1), xlo(1), time, dt
      double precision kineng(ken_l1:ken_h1,nk)
      double precision    dat(dat_l1:dat_h1,nc)
      integer    level, grid_no

      integer    i

      do i = lo(1), hi(1)
         kineng(i,1) = 0.5d0 * (dat(i,2)**2) / dat(i,1)
      end do

      end subroutine ca_derkineng

      subroutine ca_dernull(kineng,ken_l1,ken_h1,nk, &
                            dat,dat_l1,dat_h1,nc,lo,hi,domlo, &
                            domhi,delta,xlo,time,dt,bc,level,grid_no)
        !
        ! This is a derived routine that does nothing.
        !
      implicit none

      integer          lo(1), hi(1)
      integer          ken_l1,ken_h1,nk
      integer          dat_l1,dat_h1,nc
      integer          domlo(1), domhi(1)
      integer          bc(1,2,nc)
      double precision delta(1), xlo(1), time, dt
      double precision kineng(ken_l1:ken_h1,nk)
      double precision    dat(dat_l1:dat_h1,nc)
      integer    level, grid_no

      end subroutine ca_dernull

