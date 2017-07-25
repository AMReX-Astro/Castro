module derive_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  public

contains
  
! All subroutines in this file must be threadsafe because they are called
! inside OpenMP parallel regions.

  subroutine ca_derstate(state,s_lo,s_hi,nv, &
                         dat,d_lo,d_hi,nc,lo,hi,domlo, &
                         domhi,delta,xlo,time,dt,bc,level,grid_no) &
                         bind(C, name="ca_derstate")
    !
    ! The incoming   "dat" vector contains (rho,T,(rho X)_1)
    ! The outgoing "state" vector contains (rho,T,X_1)
    !
    use amrex_fort_module, only : rt => amrex_real
    implicit none 

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nv)
    real(rt)         :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
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



  subroutine ca_dervel(vel,v_lo,v_hi,nv, &
                       dat,d_lo,d_hi,nc,lo,hi,domlo, &
                       domhi,delta,xlo,time,dt,bc,level,grid_no) &
                       bind(C, name="ca_dervel")
    !
    ! This routine will derive the velocity from the momentum.
    !
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt)         :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
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



  subroutine ca_deruplusc(vel,v_lo,v_hi,nv, &
                          dat,d_lo,d_hi,nc,lo,hi,domlo, &
                          domhi,delta,xlo,time,dt,bc,level,grid_no) &
                          bind(C, name="ca_deruplusc")

    use network, only : nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX
    use bl_constants_module
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt)         :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    real(rt)         :: rhoInv

    type (eos_t) :: eos_state

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



  subroutine ca_deruminusc(vel,v_lo,v_hi,nv, &
                           dat,d_lo,d_hi,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no) &
                           bind(C, name="ca_deruminusc")

    use network, only : nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX
    use bl_constants_module
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt)         :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    real(rt)         :: rhoInv

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

             call eos(eos_input_re, eos_state)

             vel(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1) - eos_state % cs
          end do
       end do
    end do

  end subroutine ca_deruminusc



  subroutine ca_dermagvel(magvel,v_lo,v_hi,nv, &
                          dat,d_lo,d_hi,nc,lo,hi,domlo, &
                          domhi,delta,xlo,time,dt,bc,level,grid_no) &
                          bind(C, name="ca_dermagvel")
    !
    ! This routine will derive magnitude of velocity.
    !
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: magvel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt)         ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    real(rt)         :: dat1inv

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dat1inv = 1.e0_rt/dat(i,j,k,1)
             magvel(i,j,k,1) = sqrt( (dat(i,j,k,2) * dat1inv)**2 + &
                  (dat(i,j,k,3) * dat1inv)**2 + &
                  (dat(i,j,k,4) * dat1inv)**2 )
          end do
       end do
    end do

  end subroutine ca_dermagvel



  subroutine ca_dermaggrav(maggrav,g_lo,g_hi,ng, &
                           dat,d_lo,d_hi,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no) &
                           bind(C, name="ca_dermaggrav")
    !
    ! This routine will derive magnitude of the gravity vector.
    !
    use amrex_fort_module, only : rt => amrex_real
    implicit none 

    integer          :: lo(3), hi(3)
    integer          :: g_lo(3), g_hi(3), ng
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: maggrav(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),ng)
    real(rt)         ::     dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
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



  subroutine ca_derradialvel(radvel,v_lo,v_hi,nv, &
                             dat,d_lo,d_hi,nc,lo,hi,domlo, &
                             domhi,delta,xlo,time,dt,bc,level,grid_no) &
                             bind(C, name="ca_derradialvel")
    !
    ! This routine will derive the radial velocity.
    !
    use bl_constants_module
    use prob_params_module, only: center

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: radvel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt)         ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    real(rt)         :: x, y, z, r

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



  subroutine ca_dermagmom(magmom,m_lo,m_hi,nv, &
                          dat,d_lo,d_hi,nc,lo,hi,domlo, &
                          domhi,delta,xlo,time,dt,bc,level,grid_no) &
                          bind(C, name="ca_dermagmom")
    !
    ! This routine will derive magnitude of momentum.
    !
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: m_lo(3), m_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: magmom(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),nv)
    real(rt)         ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
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



  subroutine ca_derangmomx(L,L_lo,L_hi,ncomp_L, &
                           u,u_lo,u_hi,ncomp_u, &
                           lo,hi,domlo,domhi, &
                           dx,xlo,time,dt,bc,level,grid_no) &
                           bind(C, name="ca_derangmomx")

    use bl_constants_module, only: HALF
    use math_module, only: cross_product

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: L_lo(3), L_hi(3), ncomp_L ! == 1
    integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt)         :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
    real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt)         :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    integer          :: i, j, k
    real(rt)         :: loc(3), mom(3), ang_mom(3), rho

    do k = lo(3), hi(3)
       loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
       do j = lo(2), hi(2)
          loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
          do i = lo(1), hi(1)
             loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

             rho = u(i,j,k,1)
             mom = u(i,j,k,2:4)
             ang_mom = cross_product(loc, mom)

             L(i,j,k,1) = ang_mom(1)

          enddo
       enddo
    enddo

  end subroutine ca_derangmomx



  subroutine ca_derangmomy(L,L_lo,L_hi,ncomp_L, &
                           u,u_lo,u_hi,ncomp_u, &
                           lo,hi,domlo,domhi, &
                           dx,xlo,time,dt,bc,level,grid_no) &
                           bind(C, name="ca_derangmomy")

    use bl_constants_module, only: HALF
    use math_module, only: cross_product

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: L_lo(3), L_hi(3), ncomp_L ! == 1
    integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt)         :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
    real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt)         :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    integer          :: i, j, k
    real(rt)         :: loc(3), mom(3), ang_mom(3), rho

    do k = lo(3), hi(3)
       loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
       do j = lo(2), hi(2)
          loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
          do i = lo(1), hi(1)
             loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

             rho = u(i,j,k,1)
             mom = u(i,j,k,2:4)
             ang_mom = cross_product(loc, mom)

             L(i,j,k,1) = ang_mom(2)

          enddo
       enddo
    enddo

  end subroutine ca_derangmomy



  subroutine ca_derangmomz(L,L_lo,L_hi,ncomp_L, &
                           u,u_lo,u_hi,ncomp_u, &
                           lo,hi,domlo,domhi, &
                           dx,xlo,time,dt,bc,level,grid_no) &
                           bind(C, name="ca_derangmomz")

    use bl_constants_module, only: HALF
    use math_module, only: cross_product

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: L_lo(3), L_hi(3), ncomp_L ! == 1
    integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt)         :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
    real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt)         :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    integer          :: i, j, k
    real(rt)         :: loc(3), mom(3), ang_mom(3), rho

    do k = lo(3), hi(3)
       loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
       do j = lo(2), hi(2)
          loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
          do i = lo(1), hi(1)
             loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

             rho = u(i,j,k,1)
             mom = u(i,j,k,2:4)
             ang_mom = cross_product(loc, mom)

             L(i,j,k,1) = ang_mom(3)

          enddo
       enddo
    enddo

  end subroutine ca_derangmomz



  subroutine ca_derpres(p,p_lo,p_hi,ncomp_p, &
                        u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                        domhi,dx,xlo,time,dt,bc,level,grid_no) &
                        bind(C, name="ca_derpres")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use bl_constants_module
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: p_lo(3), p_hi(3), ncomp_p
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    real(rt)         :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
    real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt)         :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    real(rt)         :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoInv = ONE / u(i,j,k,URHO)

             eos_state % rho  = u(i,j,k,URHO)
             eos_state % T    = u(i,j,k,UTEMP)
             eos_state % e    = u(i,j,k,UEINT) * rhoInv
             eos_state % xn   = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux  = u(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos(eos_input_re, eos_state)

             p(i,j,k,1) = eos_state % p
          enddo
       enddo
    enddo

  end subroutine ca_derpres



  subroutine ca_dereint1(e,e_lo,e_hi,ncomp_e, &
                         u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                         domhi,dx,xlo,time,dt,bc,level,grid_no) &
                         bind(C, name="ca_dereint1")

    use bl_constants_module
    use meth_params_module, only: URHO, UMX, UMY, UMZ, UEDEN 
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: e_lo(3), e_hi(3), ncomp_e
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    real(rt)         :: e(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),ncomp_e)
    real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt)         :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    real(rt)         :: rhoInv, ux, uy, uz
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



  subroutine ca_dereint2(e,e_lo,e_hi,ncomp_e, &
                         u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                         domhi,dx,xlo,time,dt,bc,level,grid_no) &
                         bind(C, name="ca_dereint2")

    use meth_params_module, only: URHO, UEINT

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: e_lo(3), e_hi(3), ncomp_e
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    real(rt)         :: e(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),ncomp_e)
    real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt)         :: dx(3), xlo(3), time, dt
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



  subroutine ca_dersoundspeed(c,c_lo,c_hi,ncomp_c, &
                              u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                              domhi,dx,xlo,time,dt,bc,level,grid_no) &
                              bind(C, name="ca_dersoundspeed")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use bl_constants_module
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: c_lo(3), c_hi(3), ncomp_c
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    real(rt)         :: c(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),ncomp_c)
    real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt)         :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    real(rt)         :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

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



  subroutine ca_dermachnumber(mach,m_lo,m_hi,ncomp_mach, &
                              u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                              domhi,dx,xlo,time,dt,bc,level,grid_no) &
                              bind(C, name="ca_dermachnumber")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use meth_params_module, only: URHO, UMX, UMZ, UEINT, UTEMP, UFS, UFX
    use bl_constants_module
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: m_lo(3), m_hi(3), ncomp_mach
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    real(rt)         :: mach(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),ncomp_mach)
    real(rt)         ::    u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u   )
    real(rt)         :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    real(rt)         :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

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



  subroutine ca_derentropy(s,s_lo,s_hi,ncomp_s, &
                           u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                           domhi,dx,xlo,time,dt,bc,level,grid_no) &
                           bind(C, name="ca_derentropy")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use bl_constants_module
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3), ncomp_s
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    real(rt)         :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ncomp_s)
    real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt)         :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    real(rt)         :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

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



  subroutine ca_derenuctimescale(t,t_lo,t_hi,ncomp_t, &
                                 u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                                 domhi,dx,xlo,time,dt,bc,level,grid_no) &
                                 bind(C, name="ca_derenuctimescale")

    use bl_constants_module, only: ZERO, ONE
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use network, only: nspec, naux
    use prob_params_module, only: dim
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use amrex_fort_module, only : rt => amrex_real

    implicit none
    

    integer          :: lo(3), hi(3)
    integer          :: t_lo(3), t_hi(3), ncomp_t
    integer          :: u_lo(3), u_hi(3), ncomp_u
    integer          :: domlo(3), domhi(3)
    real(rt)         :: t(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3),ncomp_t)
    real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u) ! NVAR, enuc
    real(rt)         :: dx(3), xlo(3), time, dt
    integer          :: bc(3,2,ncomp_u), level, grid_no

    integer          :: i, j, k
    real(rt)         :: rhoInv, eint, enuc, t_s, t_e

    type (eos_t)     :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             enuc = abs(u(i,j,k,ncomp_u))

             if (enuc > 1.e-100_rt) then

                rhoInv = ONE / u(i,j,k,URHO)

                eint = u(i,j,k,UEINT) * rhoInv

                t_e = eint / enuc

                ! Calculate sound-speed

                eos_state % rho  = u(i,j,k,URHO)
                eos_state % T    = u(i,j,k,UTEMP)
                eos_state % e    = eint
                eos_state % xn   = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
                eos_state % aux  = u(i,j,k,UFX:UFX+naux-1) * rhoInv

                call eos(eos_input_re, eos_state)

                t_s = minval(dx(1:dim)) / eos_state % cs

                t(i,j,k,1) = t_s / t_e

             else

                t(i,j,k,1) = ZERO

             endif

          enddo
       enddo
    enddo

  end subroutine ca_derenuctimescale



  subroutine ca_derspec(spec,s_lo,s_hi,nv, &
                        dat,d_lo,d_hi,nc,lo,hi,domlo, &
                        domhi,delta,xlo,time,dt,bc,level,grid_no) &
                        bind(C, name="ca_derspec")
    !
    ! This routine derives the mass fractions of the species.
    !
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: spec(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nv)
    real(rt)         ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
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



  subroutine ca_derlogden(logden,l_lo,l_hi,nd, &
                          dat,d_lo,d_hi,nc, &
                          lo,hi,domlo,domhi,delta, &
                          xlo,time,dt,bc,level,grid_no) &
                          bind(C, name="ca_derlogden")
    
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: l_lo(3), l_hi(3), nd
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3), level, grid_no
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: logden(l_lo(1):l_hi(1),l_lo(2):l_hi(2),l_lo(3):l_hi(3),nd)
    real(rt)         ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             logden(i,j,k,1) = dlog10(dat(i,j,k,1))
          end do
       end do
    end do

  end subroutine ca_derlogden



  subroutine ca_dermagvort(vort,v_lo,v_hi,nv, & 
                           dat,d_lo,d_hi,nc,lo,hi,domlo, &
                           domhi,delta,xlo,time,dt,bc,level,grid_no) &
                           bind(C, name="ca_dermagvort")
    
    !
    ! This routine will calculate vorticity
    !     

    use bl_constants_module
    use prob_params_module, only: dg

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: v_lo(3), v_hi(3), nv
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3), level, grid_no
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: vort(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt)         ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k
    real(rt)         :: uy, uz, vx, vz, wx, wy, v1, v2, v3
    real(rt)         :: ldat(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,2:4)

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



  subroutine ca_derdivu(divu,u_lo,u_hi,nd, &
                        dat,d_lo,d_hi,nc, &
                        lo,hi,domlo,domhi,delta, &
                        xlo,time,dt,bc,level,grid_no) &
                        bind(C, name="ca_derdivu")
    !
    ! This routine will calculate the divergence of velocity.
    !

    use bl_constants_module
    use prob_params_module, only: dg

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_lo(3), u_hi(3), nd
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: divu(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt)         ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k
    real(rt)         :: ulo, uhi, vlo, vhi, wlo, whi

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



  subroutine ca_derkineng(kineng,k_lo,k_hi,nk, &
                          dat,d_lo,d_hi,nc, &
                          lo,hi,domlo,domhi,delta, &
                          xlo,time,dt,bc,level,grid_no) &
                          bind(C, name="ca_derkineng")
    !
    ! This routine will derive kinetic energy = 1/2 rho (u^2 + v^2 + w^2)
    !

    use bl_constants_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: k_lo(3), k_hi(3), nk
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: kineng(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3),nk)
    real(rt)         ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

    integer          :: i, j, k

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



  subroutine ca_dernull(dnull,k_lo,k_hi,nk, &
                        dat,d_lo,d_hi,nc, &
                        lo,hi,domlo,domhi,delta, &
                        xlo,time,dt,bc,level,grid_no) &
                        bind(C, name="ca_dernull")
    !
    ! This routine is used by particle_count.  Yes it does nothing.
    !
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: k_lo(3), k_hi(3), nk
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: dnull(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3),nk)
    real(rt)         ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no

  end subroutine ca_dernull


#ifdef DIFFUSION
  subroutine ca_dercond(cond,u_lo,u_hi,nd, &
                        state,d_lo,d_hi,nc, &
                        lo,hi,domlo,domhi,delta, &
                        xlo,time,dt,bc,level,grid_no) &
                        bind(C, name="ca_dercond")
    !
    ! This routine will calculate the thermal conductivity
    !

    use bl_constants_module
    use meth_params_module, only: diffuse_cutoff_density, &
                                  URHO, UEINT, UTEMP, UFS, UFX
    use eos_type_module, only: eos_input_re, eos_t
    use eos_module, only: eos
    use network
    use conductivity_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_lo(3), u_hi(3), nd
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: cond(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt)         :: state(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no
    real(rt)         :: coeff
    integer          :: i, j, k

    type(eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho    = state(i,j,k,URHO)
             eos_state%T      = state(i,j,k,UTEMP)
             eos_state%e      = state(i,j,k,UEINT)/state(i,j,k,URHO)
             eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)/ state(i,j,k,URHO)
             eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)/ state(i,j,k,URHO)
             call eos(eos_input_re,eos_state)

             if (eos_state%rho > diffuse_cutoff_density) then
                call conductivity(eos_state, coeff)
             else
                coeff = ZERO
             endif

             cond(i,j,k,1) = coeff

          enddo
       enddo
    enddo

  end subroutine ca_dercond


  subroutine ca_derdiffcoeff(diff,u_lo,u_hi,nd, &
                             state,d_lo,d_hi,nc, &
                             lo,hi,domlo,domhi,delta, &
                             xlo,time,dt,bc,level,grid_no) &
                             bind(C, name="ca_derdiffcoeff")
    !
    ! This routine will calculate the thermal conductivity
    !

    use bl_constants_module
    use meth_params_module, only: diffuse_cutoff_density, &
                                  URHO, UEINT, UTEMP, UFS, UFX
    use eos_type_module, only: eos_input_re, eos_t
    use eos_module, only: eos
    use network
    use conductivity_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_lo(3), u_hi(3), nd
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: diff(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt)         :: state(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no
    real(rt)         :: coeff
    integer          :: i, j, k

    type(eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state%rho    = state(i,j,k,URHO)
             eos_state%T      = state(i,j,k,UTEMP)
             eos_state%e      = state(i,j,k,UEINT)/state(i,j,k,URHO)
             eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)/ state(i,j,k,URHO)
             eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)/ state(i,j,k,URHO)
             call eos(eos_input_re,eos_state)

             if (eos_state%rho > diffuse_cutoff_density) then
                call conductivity(eos_state, coeff)
             else
                coeff = ZERO
             endif

             diff(i,j,k,1) = coeff/(eos_state%rho * eos_state%cv)

          enddo
       enddo
    enddo

  end subroutine ca_derdiffcoeff


  subroutine ca_derdiffterm(diff,u_lo,u_hi,nd, &
                            state,d_lo,d_hi,nc, &
                            lo,hi,domlo,domhi,delta, &
                            xlo,time,dt,bc,level,grid_no) &
                            bind(C, name="ca_derdiffterm")
    !
    ! This routine will calculate the thermal conductivity
    !

    use bl_constants_module
    use meth_params_module, only: UTEMP
    use prob_params_module, only: dim
    use diffusion_module
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_lo(3), u_hi(3), nd
    integer          :: d_lo(3), d_hi(3), nc
    integer          :: domlo(3), domhi(3)
    integer          :: bc(3,2,nc)
    real(rt)         :: delta(3), xlo(3), time, dt
    real(rt)         :: diff(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt)         :: state(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer          :: level, grid_no
    real(rt), allocatable  :: coeff_x(:,:,:), coeff_y(:,:,:), coeff_z(:,:,:)
    real(rt) :: diff_term
    integer          :: i, j, k


    ! allocate space for edge-centered conductivities
    allocate(coeff_x(d_lo(1):d_hi(1), d_lo(2):d_hi(2), d_lo(3):d_hi(3)))
    allocate(coeff_y(d_lo(1):d_hi(1), d_lo(2):d_hi(2), d_lo(3):d_hi(3)))
    allocate(coeff_z(d_lo(1):d_hi(1), d_lo(2):d_hi(2), d_lo(3):d_hi(3)))

    call ca_fill_temp_cond(lo, hi, &
                           state, d_lo, d_hi, &
                           coeff_x, d_lo, d_hi, &
                           coeff_y, d_lo, d_hi, &
                           coeff_z, d_lo, d_hi)

    ! create the diff term
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! x
             diff_term = &
                  coeff_x(i+1,j,k)*(state(i+1,j,k,UTEMP) - state(i,  j,k,UTEMP))/delta(1) - &
                  coeff_x(i  ,j,k)*(state(i  ,j,k,UTEMP) - state(i-1,j,k,UTEMP))/delta(1)

             ! y
             if (dim >= 2) then
                diff_term = diff_term + &
                     coeff_y(i,j+1,k)*(state(i,j+1,k,UTEMP) - state(i,j  ,k,UTEMP))/delta(2) - &
                     coeff_y(i,j  ,k)*(state(i,j,  k,UTEMP) - state(i,j-1,k,UTEMP))/delta(2)
             endif
             
             ! z
             if (dim == 3) then
                diff_term = diff_term + &
                     coeff_z(i,j,k+1)*(state(i,j,k+1,UTEMP) - state(i,j,k  ,UTEMP))/delta(3) - &
                     coeff_z(i,j,k  )*(state(i,j,k  ,UTEMP) - state(i,j,k-1,UTEMP))/delta(3)
             endif

             diff(i,j,k,1) = diff_term

          enddo
       enddo
    enddo

    deallocate(coeff_x, coeff_y, coeff_z)
    
  end subroutine ca_derdiffterm
  
#endif


end module derive_module
