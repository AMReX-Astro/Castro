! All subroutines in this file must be threadsafe because they are called
! inside OpenMP paralle regions.

!-----------------------------------------------------------------------

      subroutine ca_derstate(state,s_lo,s_hi,nv, &
                             dat,d_lo,d_hi,nc,lo,hi,domlo, &
                             domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! The incoming   "dat" vector contains (rho,T,(rho X)_1)
      ! The outgoing "state" vector contains (rho,T,X_1)
      !
      implicit none 

      integer          :: lo(3), hi(3)
      integer          :: s_lo(3), s_hi(3), nv
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3)
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nv)
      double precision :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
      integer          :: level, grid_no
 
      integer          :: i, j, k

      if (nv .ne. 3) then
          print *,'... confusion in derstate ... nv should be 3 but is ',nv
          call bl_error('Error:: Derive_nd.f90 :: ca_derstate')
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

      subroutine ca_dervel(vel,v_lo,v_hi,nv, &
                           dat,d_lo,d_hi,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive the velocity from the momentum.
      !
      implicit none

      integer          :: lo(3), hi(3)
      integer          :: v_lo(3), v_hi(3), nv
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3)
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
      double precision :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
      integer          :: level, grid_no
 
      integer          :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               vel(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
            end do
         end do
      end do
 
      end subroutine ca_dervel

!-----------------------------------------------------------------------

      subroutine ca_deruplusc(vel,v_lo,v_hi,nv, &
                              dat,d_lo,d_hi,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : URHO, UMX, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: v_lo(3), v_hi(3), nv
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3)
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
      double precision :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
      integer          :: level, grid_no
 
      integer          :: i, j, k
      double precision :: rhoInv

      type (eos_t) :: eos_state

      if (allow_negative_energy .eq. 0) eos_state % reset = .true.

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhoInv = ONE / dat(i,j,k,URHO)
               
               eos_state % e     = dat(i,j,k,UEINT) * rhoInv
               eos_state % T     = dat(i,j,k,UTEMP)
               eos_state % rho   = dat(i,j,k,URHO)
               eos_state % xn  = dat(i,j,k,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux = dat(i,j,k,UFX:UFX+naux-1) * rhoInv
               
               call eos(eos_input_re, eos_state)

               vel(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1) + eos_state % cs
            enddo
         enddo
      enddo

      end subroutine ca_deruplusc

!-----------------------------------------------------------------------

      subroutine ca_deruminusc(vel,v_lo,v_hi,nv, &
                               dat,d_lo,d_hi,nc,lo,hi,domlo, &
                               domhi,delta,xlo,time,dt,bc,level,grid_no)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : URHO, UMX, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy
      use bl_constants_module

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: v_lo(3), v_hi(3), nv
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3)
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
      double precision :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
      integer          :: level, grid_no
 
      integer          :: i, j, k
      double precision :: rhoInv

      type (eos_t) :: eos_state
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhoInv = ONE / dat(i,j,k,URHO)
               
               eos_state % e   = dat(i,j,k,UEINT) * rhoInv
               eos_state % T   = dat(i,j,k,UTEMP)
               eos_state % rho = dat(i,j,k,URHO)
               eos_state % xn  = dat(i,j,k,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux = dat(i,j,k,UFX:UFX+naux-1) * rhoInv

               if (allow_negative_energy .eq. 0) eos_state % reset = .true.      
      
               call eos(eos_input_re, eos_state)

               vel(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1) - eos_state % cs
            end do
         end do
      end do

      end subroutine ca_deruminusc

!-----------------------------------------------------------------------

      subroutine ca_dermagvel(magvel,v_lo,v_hi,nv, &
                              dat,d_lo,d_hi,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive magnitude of velocity.
      !
      implicit none

      integer          :: lo(3), hi(3)
      integer          :: v_lo(3), v_hi(3), nv
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3)
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: magvel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
      double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
      integer          :: level, grid_no

      integer          :: i, j, k
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

      subroutine ca_dermaggrav(maggrav,g_lo,g_hi,ng, &
                               dat,d_lo,d_hi,nc,lo,hi,domlo, &
                               domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive magnitude of the gravity vector.
      !
      implicit none 

      integer          :: lo(3), hi(3)
      integer          :: g_lo(3), g_hi(3), ng
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3)
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: maggrav(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),ng)
      double precision ::     dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
      integer          :: level, grid_no

      integer          :: i, j, k

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

      subroutine ca_derradialvel(radvel,v_lo,v_hi,nv, &
                                 dat,d_lo,d_hi,nc,lo,hi,domlo, &
                                 domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive the radial velocity.
      !
      use bl_constants_module
      use prob_params_module, only: center

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: v_lo(3), v_hi(3), nv
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3)
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: radvel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
      double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
      integer          :: level, grid_no

      integer          :: i, j, k
      double precision :: x, y, z, r

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

      subroutine ca_dermagmom(magmom,m_lo,m_hi,nv, &
                              dat,d_lo,d_hi,nc,lo,hi,domlo, &
                              domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive magnitude of momentum.
      !
      implicit none

      integer          :: lo(3), hi(3)
      integer          :: m_lo(3), m_hi(3), nv
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3)
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: magmom(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),nv)
      double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
      integer          :: level, grid_no

      integer          :: i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               magmom(i,j,k,1) = sqrt( dat(i,j,k,1)**2 + dat(i,j,k,2)**2 + dat(i,j,k,3)**2 )
            end do
         end do
      end do

      end subroutine ca_dermagmom

!-----------------------------------------------------------------------

      subroutine ca_derpres(p,p_lo,p_hi,ncomp_p, &
           u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only: nspec, naux
      use eos_module
      use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX, &
                                    allow_negative_energy
      use bl_constants_module

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: p_lo(3), p_hi(3), ncomp_p
      integer          :: u_lo(3), u_hi(3), ncomp_u
      integer          :: domlo(3), domhi(3)
      double precision :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
      double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
      double precision :: dx(3), xlo(3), time, dt
      integer          :: bc(3,2,ncomp_u), level, grid_no

      double precision :: rhoInv
      integer          :: i, j, k

      type (eos_t) :: eos_state

      if (allow_negative_energy .eq. 0) eos_state % reset = .true.      
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhoInv = ONE / u(i,j,k,URHO)
               
               eos_state % rho  = u(i,j,k,URHO)
               eos_state % T    = u(i,j,k,UTEMP)
               eos_state % e    = u(i,j,k,UEINT) * rhoInv
               eos_state % xn = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux = u(i,j,k,UFX:UFX+naux-1) * rhoInv

               call eos(eos_input_re, eos_state)

               p(i,j,k,1) = eos_state % p
            enddo
         enddo
      enddo

      end subroutine ca_derpres

!-----------------------------------------------------------------------

      subroutine ca_dereint1(e,e_lo,e_hi,ncomp_e, &
           u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use bl_constants_module
      use meth_params_module, only: URHO, UMX, UMY, UMZ, UEDEN 

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: e_lo(3), e_hi(3), ncomp_e
      integer          :: u_lo(3), u_hi(3), ncomp_u
      integer          :: domlo(3), domhi(3)
      double precision :: e(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),ncomp_e)
      double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
      double precision :: dx(3), xlo(3), time, dt
      integer          :: bc(3,2,ncomp_u), level, grid_no

      double precision :: rhoInv, ux, uy, uz
      integer          :: i, j, k
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

      subroutine ca_dereint2(e,e_lo,e_hi,ncomp_e, &
           u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use meth_params_module, only: URHO, UEINT

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: e_lo(3), e_hi(3), ncomp_e
      integer          :: u_lo(3), u_hi(3), ncomp_u
      integer          :: domlo(3), domhi(3)
      double precision :: e(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),ncomp_e)
      double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
      double precision :: dx(3), xlo(3), time, dt
      integer          :: bc(3,2,ncomp_u), level, grid_no

      integer          :: i, j, k
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

      subroutine ca_dersoundspeed(c,c_lo,c_hi,ncomp_c, &
           u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only: nspec, naux
      use eos_module
      use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX, &
                                    allow_negative_energy
      use bl_constants_module

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: c_lo(3), c_hi(3), ncomp_c
      integer          :: u_lo(3), u_hi(3), ncomp_u
      integer          :: domlo(3), domhi(3)
      double precision :: c(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),ncomp_c)
      double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
      double precision :: dx(3), xlo(3), time, dt
      integer          :: bc(3,2,ncomp_u), level, grid_no

      double precision :: rhoInv
      integer          :: i, j, k

      type (eos_t) :: eos_state

      if (allow_negative_energy .eq. 0) eos_state % reset = .true.
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhoInv = ONE / u(i,j,k,URHO)
               
               eos_state % rho  = u(i,j,k,URHO)
               eos_state % T    = u(i,j,k,UTEMP)
               eos_state % e    = u(i,j,k,UEINT) * rhoInv
               eos_state % xn = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux = u(i,j,k,UFX:UFX+naux-1) * rhoInv

               call eos(eos_input_re, eos_state)

               c(i,j,k,1) = eos_state % cs
            enddo
         enddo
      enddo

      end subroutine ca_dersoundspeed

!-----------------------------------------------------------------------

      subroutine ca_dermachnumber(mach,m_lo,m_hi,ncomp_mach, &
           u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only: nspec, naux
      use eos_module
      use meth_params_module, only: URHO, UMX, UMZ, UEINT, UTEMP, UFS, UFX, &
                                    allow_negative_energy
      use bl_constants_module

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: m_lo(3), m_hi(3), ncomp_mach
      integer          :: u_lo(3), u_hi(3), ncomp_u
      integer          :: domlo(3), domhi(3)
      double precision :: mach(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),ncomp_mach)
      double precision ::    u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u   )
      double precision :: dx(3), xlo(3), time, dt
      integer          :: bc(3,2,ncomp_u), level, grid_no

      double precision :: rhoInv, ux, uy, uz
      integer          :: i, j, k

      type (eos_t) :: eos_state

      if (allow_negative_energy .eq. 0) eos_state % reset = .true.

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhoInv = ONE / u(i,j,k,URHO)
               
               eos_state % rho  = u(i,j,k,URHO)
               eos_state % T    = u(i,j,k,UTEMP)
               eos_state % e    = u(i,j,k,UEINT) * rhoInv
               eos_state % xn = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux = u(i,j,k,UFX:UFX+naux-1) * rhoInv

               call eos(eos_input_re, eos_state)

               mach(i,j,k,1) = sum(u(i,j,k,UMX:UMZ)**2)**0.5 / u(i,j,k,URHO) / eos_state % cs
            enddo
         enddo
      enddo

      end subroutine ca_dermachnumber

!-----------------------------------------------------------------------

      subroutine ca_derentropy(s,s_lo,s_hi,ncomp_s, &
           u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
           domhi,dx,xlo,time,dt,bc,level,grid_no)

      use network, only: nspec, naux
      use eos_module
      use meth_params_module, only: URHO, UMX, UMY, UMZ, UEINT, UTEMP, UFS, UFX, &
                                    allow_negative_energy
      use bl_constants_module

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: s_lo(3), s_hi(3), ncomp_s
      integer          :: u_lo(3), u_hi(3), ncomp_u
      integer          :: domlo(3), domhi(3)
      double precision :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ncomp_s)
      double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
      double precision :: dx(3), xlo(3), time, dt
      integer          :: bc(3,2,ncomp_u), level, grid_no

      double precision :: rhoInv
      integer          :: i, j, k

      type (eos_t) :: eos_state

      if (allow_negative_energy .eq. 0) eos_state % reset = .true.      
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhoInv = ONE / u(i,j,k,URHO)
               
               eos_state % rho = u(i,j,k,URHO)
               eos_state % T   = u(i,j,k,UTEMP)
               eos_state % e   = u(i,j,k,UEINT) * rhoInv
               eos_state % xn  = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux = u(i,j,k,UFX:UFX+naux-1) * rhoInv
      
               call eos(eos_input_re, eos_state)

               s(i,j,k,1) = eos_state % s
            enddo
         enddo
      enddo

      end subroutine ca_derentropy

!-----------------------------------------------------------------------

      subroutine ca_derspec(spec,s_lo,s_hi,nv, &
                            dat,d_lo,d_hi,nc,lo,hi,domlo, &
                            domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine derives the mass fractions of the species.
      !
      implicit none

      integer          :: lo(3), hi(3)
      integer          :: s_lo(3), s_hi(3), nv
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3)
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: spec(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nv)
      double precision ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
      integer          :: level, grid_no
 
      integer          :: i, j, k
 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               spec(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
            end do
         end do
      end do
 
      end subroutine ca_derspec

!-----------------------------------------------------------------------

      subroutine ca_derlogden(logden,l_lo,l_hi,nd, &
                              dat,d_lo,d_hi,nc, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      implicit none

      integer          :: lo(3), hi(3)
      integer          :: l_lo(3), l_hi(3), nd
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3), level, grid_no
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: logden(l_lo(1):l_hi(1),l_lo(2):l_hi(2),l_lo(3):l_hi(3),nd)
      double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
 
      integer          :: i, j, k
 
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               logden(i,j,k,1) = dlog10(dat(i,j,k,1))
            end do
         end do
      end do
 
      end subroutine ca_derlogden

!-----------------------------------------------------------------------

      subroutine ca_dermagvort(vort,v_lo,v_hi,nv, & 
                               dat,d_lo,d_hi,nc,lo,hi,domlo, &
                               domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will calculate vorticity
      !     

      use bl_constants_module
      use prob_params_module, only: dg
      
      implicit none

      integer          :: lo(3), hi(3)
      integer          :: v_lo(3), v_hi(3), nv
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3), level, grid_no
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: vort(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
      double precision ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

      integer          :: i, j, k
      double precision :: uy, uz, vx, vz, wx, wy, v1, v2, v3
      double precision :: ldat(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,2:4)

      ldat = ZERO

      uy = ZERO
      uz = ZERO
      vx = ZERO
      vz = ZERO
      wx = ZERO
      wy = ZERO
      
      !
      ! Convert momentum to velocity.
      !
      do k = lo(3)-1*dg(3), hi(3)+1*dg(3)
         do j = lo(2)-1*dg(2), hi(2)+1*dg(2)
            do i = lo(1)-1*dg(1), hi(1)+1*dg(1)
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

               vx = HALF * (ldat(i+1,j,k,3) - ldat(i-1,j,k,3)) / delta(1)
               wx = HALF * (ldat(i+1,j,k,4) - ldat(i-1,j,k,4)) / delta(1)
                  
               if (delta(2) > ZERO) then
                  uy = HALF * (ldat(i,j+1,k,2) - ldat(i,j-1,k,2)) / delta(2)
                  wy = HALF * (ldat(i,j+1,k,4) - ldat(i,j-1,k,4)) / delta(2)
               endif

               if (delta(3) > ZERO) then
                  uz = HALF * (ldat(i,j,k+1,2) - ldat(i,j,k-1,2)) / delta(3)
                  vz = HALF * (ldat(i,j,k+1,3) - ldat(i,j,k-1,3)) / delta(3)
               endif

               v1 = wy - vz
               v2 = uz - wx
               v3 = vx - uy
               vort(i,j,k,1) = sqrt(v1*v1 + v2*v2 + v3*v3)

            end do
         end do
      end do

      end subroutine ca_dermagvort

!-----------------------------------------------------------------------

      subroutine ca_derdivu(divu,u_lo,u_hi,nd, &
                            dat,d_lo,d_hi,nc, &
                            lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will calculate the divergence of velocity.
      !

      use bl_constants_module
      use prob_params_module, only: dg
      
      implicit none

      integer          :: lo(3), hi(3)
      integer          :: u_lo(3), u_hi(3), nd
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3)
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: divu(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
      double precision ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
      integer          :: level, grid_no

      integer          :: i, j, k
      double precision :: ulo, uhi, vlo, vhi, wlo, whi

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               uhi = dat(i+1*dg(1),j,k,2) / dat(i+1*dg(1),j,k,1)
               ulo = dat(i-1*dg(1),j,k,2) / dat(i-1*dg(1),j,k,1)
               vhi = dat(i,j+1*dg(2),k,3) / dat(i,j+1*dg(2),k,1)
               vlo = dat(i,j-1*dg(2),k,3) / dat(i,j-1*dg(2),k,1)
               whi = dat(i,j,k+1*dg(3),4) / dat(i,j,k+1*dg(3),1)
               wlo = dat(i,j,k-1*dg(3),4) / dat(i,j,k-1*dg(3),1)
               divu(i,j,k,1) = HALF * (uhi-ulo) / delta(1)
               if (delta(2) > ZERO) then
                  divu(i,j,k,1) = divu(i,j,k,1) + HALF * (vhi-vlo) / delta(2)
               endif
               if (delta(3) > ZERO) then
                  divu(i,j,k,1) = divu(i,j,k,1) + HALF * (whi-wlo) / delta(3)
               endif
            end do
         end do
      end do

      end subroutine ca_derdivu

!-----------------------------------------------------------------------

      subroutine ca_derkineng(kineng,k_lo,k_hi,nk, &
                              dat,d_lo,d_hi,nc, &
                              lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine will derive kinetic energy = 1/2 rho (u^2 + v^2 + w^2)
      !

      use bl_constants_module
      use prob_params_module, only: problo, center
      use meth_params_module, only: hybrid_hydro

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: k_lo(3), k_hi(3), nk
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3)
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: kineng(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3),nk)
      double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
      integer          :: level, grid_no

      integer          :: i, j, k

      double precision :: u, v, w, rho, rhoInv
      double precision :: x, y, R
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            y = problo(2) + (dble(j) + HALF) * delta(2) - center(2)
            do i = lo(1), hi(1)
               x = problo(1) + (dble(i) + HALF) * delta(1) - center(1)

               R = sqrt( x**2 + y**2 )

               rho = dat(i,j,k,1)
               rhoInv = ONE / rho

               if (hybrid_hydro .eq. 1) then
                  u = (dat(i,j,k,2) * x / R    - dat(i,j,k,3) * y / R**2) * rhoInv
                  v = (dat(i,j,k,3) * x / R**2 + dat(i,j,k,2) * y / R   ) * rhoInv
               else                  
                  u = dat(i,j,k,2) * rhoInv
                  v = dat(i,j,k,3) * rhoInv
               endif
               w = dat(i,j,k,4) * rhoInv

               kineng(i,j,k,1) = HALF * rho * (u**2 + v**2 + w**2)
               
            end do
         end do
      end do

      end subroutine ca_derkineng

!-----------------------------------------------------------------------

      subroutine ca_dernull(kineng,k_lo,k_hi,nk, &
                            dat,d_lo,d_hi,nc, &
                            lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
      !
      ! This routine is used by particle_count.  Yes it does nothing.
      !
      implicit none

      integer          :: lo(3), hi(3)
      integer          :: k_lo(3), k_hi(3), nk
      integer          :: d_lo(3), d_hi(3), nc
      integer          :: domlo(3), domhi(3)
      integer          :: bc(3,2,nc)
      double precision :: delta(3), xlo(3), time, dt
      double precision :: kineng(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3),nk)
      double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
      integer          :: level, grid_no

      end subroutine ca_dernull

