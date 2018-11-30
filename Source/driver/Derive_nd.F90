module derive_module

  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  public

contains
  
! All subroutines in this file must be threadsafe because they are called
! inside OpenMP parallel regions.

  subroutine derstate(state,s_lo,s_hi,nv, &
                         dat,d_lo,d_hi,nc,lo,hi,domlo, &
                         domhi,delta,xlo) bind(C, name="derstate")
    !
    ! The incoming   "dat" vector contains (rho,T,(rho X)_1)
    ! The outgoing "state" vector contains (rho,T,X_1)
    !
    use amrex_fort_module, only : rt => amrex_real

    implicit none 

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in), value :: nv, nc
    integer, intent(in)     :: s_lo(3), s_hi(3)
    integer, intent(in)     :: d_lo(3), d_hi(3)
    integer, intent(in)     :: domlo(3), domhi(3)
    real(rt), intent(in)    :: delta(3), xlo(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nv)
    real(rt), intent(in)    :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k

    !$gpu

#ifndef AMREX_USE_CUDA    
    if (nv .ne. 3) then
       print *,'... confusion in derstate ... nv should be 3 but is ',nv
       call amrex_error('Error:: Derive_nd.f90 :: derstate')
    end if
#endif

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

  end subroutine derstate



  subroutine dervel(vel,v_lo,v_hi,nv, &
                       dat,d_lo,d_hi,nc,lo,hi,domlo, &
                       domhi,delta,xlo) bind(C, name="dervel")
    !
    ! This routine will derive the velocity from the momentum.
    !
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in), value :: nv, nc
    integer, intent(in)     :: v_lo(3), v_hi(3)
    integer, intent(in)     :: d_lo(3), d_hi(3)
    integer, intent(in)     :: domlo(3), domhi(3)
    real(rt), intent(in)    :: delta(3), xlo(3)
    real(rt), intent(inout) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt), intent(in)    :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             vel(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
          end do
       end do
    end do

  end subroutine dervel



  subroutine deruplusc(vel,v_lo,v_hi,nv, &
                          dat,d_lo,d_hi,nc,lo,hi,domlo, &
                          domhi,delta,xlo) bind(C, name="deruplusc")

    use network, only : nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX
    use amrex_constants_module, only : ONE
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in)  :: lo(3), hi(3)
    integer, intent(in), value :: nv, nc
    integer, intent(in)  :: v_lo(3), v_hi(3)
    integer, intent(in)  :: d_lo(3), d_hi(3)
    integer, intent(in)  :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3)
    real(rt), intent(inout):: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt), intent(in) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k
    real(rt)         :: rhoInv

    type (eos_t) :: eos_state

    !$gpu

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

  end subroutine deruplusc



  subroutine deruminusc(vel,v_lo,v_hi,nv, &
                           dat,d_lo,d_hi,nc,lo,hi,domlo, &
                           domhi,delta,xlo) bind(C, name="deruminusc")

    use network, only : nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX
    use amrex_constants_module, only : ONE
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in), value :: nv, nc
    integer, intent(in)     :: v_lo(3), v_hi(3)
    integer, intent(in)     :: d_lo(3), d_hi(3)
    integer, intent(in)     :: domlo(3), domhi(3)
    real(rt), intent(in)    :: delta(3), xlo(3)
    real(rt), intent(inout) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt), intent(in)    :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k
    real(rt)         :: rhoInv

    type (eos_t) :: eos_state

    !$gpu

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

  end subroutine deruminusc



  subroutine dermagvel(magvel,v_lo,v_hi,nv, &
                          dat,d_lo,d_hi,nc,lo,hi,domlo, &
                          domhi,delta,xlo) bind(C, name="dermagvel")
    !
    ! This routine will derive magnitude of velocity.
    !
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in)  :: lo(3), hi(3)
    integer, intent(in), value :: nv, nc
    integer, intent(in)  :: v_lo(3), v_hi(3)
    integer, intent(in)  :: d_lo(3), d_hi(3)
    integer, intent(in)  :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3)
    real(rt), intent(inout) :: magvel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt), intent(in) ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k
    real(rt)         :: dat1inv

    !$gpu

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

  end subroutine dermagvel



  subroutine dermaggrav(maggrav,g_lo,g_hi,ng, &
                           dat,d_lo,d_hi,nc,lo,hi,domlo, &
                           domhi,delta,xlo) bind(C, name="dermaggrav")
    !
    ! This routine will derive magnitude of the gravity vector.
    !
    use amrex_fort_module, only : rt => amrex_real

    implicit none 

    integer, intent(in)    :: lo(3), hi(3)
    integer, intent(in), value :: ng, nc
    integer, intent(in)    :: g_lo(3), g_hi(3)
    integer, intent(in)    :: d_lo(3), d_hi(3)
    integer, intent(in)    :: domlo(3), domhi(3)
    real(rt), intent(in)   :: delta(3), xlo(3)
    real(rt), intent(inout):: maggrav(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3),ng)
    real(rt), intent(in)   ::     dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             maggrav(i,j,k,1) = sqrt( dat(i,j,k,1)**2  + &
                  dat(i,j,k,2)**2  + &
                  dat(i,j,k,3)**2 )
          end do
       end do
    end do

  end subroutine dermaggrav



  subroutine derradialvel(radvel,v_lo,v_hi,nv, &
                             dat,d_lo,d_hi,nc,lo,hi,domlo, &
                             domhi,delta,xlo) bind(C, name="derradialvel")
    !
    ! This routine will derive the radial velocity.
    !
    use amrex_constants_module, only : HALF
    use prob_params_module, only: center
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in), value :: nv, nc
    integer, intent(in) :: v_lo(3), v_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3)
    real(rt), intent(inout) :: radvel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt), intent(in) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k
    real(rt)         :: x, y, z, r

    !$gpu

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

  end subroutine derradialvel



  subroutine dermagmom(magmom,m_lo,m_hi,nv, &
                          dat,d_lo,d_hi,nc,lo,hi,domlo, &
                          domhi,delta,xlo) bind(C, name="dermagmom")
    !
    ! This routine will derive magnitude of momentum.
    !
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in), value :: nv, nc
    integer, intent(in) :: m_lo(3), m_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3)
    real(rt), intent(inout) :: magmom(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),nv)
    real(rt), intent(in) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             magmom(i,j,k,1) = sqrt( dat(i,j,k,1)**2 + dat(i,j,k,2)**2 + dat(i,j,k,3)**2 )
          end do
       end do
    end do

  end subroutine dermagmom



  subroutine derangmomx(L,L_lo,L_hi,ncomp_L, &
                           u,u_lo,u_hi,ncomp_u, &
                           lo,hi,domlo,domhi, &
                           dx,xlo) bind(C, name="derangmomx")

    use amrex_constants_module, only: HALF
    use math_module, only: cross_product ! function
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: ncomp_L, ncomp_u ! == 1, 4
    integer, intent(in) :: L_lo(3), L_hi(3) ! == 1
    integer, intent(in) :: u_lo(3), u_hi(3) ! == 4
    integer, intent(in) :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt), intent(inout) :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt), intent(in) :: dx(3), xlo(3)

    integer          :: i, j, k
    real(rt)         :: loc(3), mom(3), ang_mom(3), rho

    !$gpu

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

  end subroutine derangmomx



  subroutine derangmomy(L,L_lo,L_hi,ncomp_L, &
                           u,u_lo,u_hi,ncomp_u, &
                           lo,hi,domlo,domhi, &
                           dx,xlo) bind(C, name="derangmomy")

    use amrex_constants_module, only: HALF
    use math_module, only: cross_product ! function
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: ncomp_L, ncomp_u ! == 1, 4
    integer, intent(in) :: L_lo(3), L_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt), intent(inout) :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt), intent(in) :: dx(3), xlo(3)

    integer          :: i, j, k
    real(rt)         :: loc(3), mom(3), ang_mom(3), rho

    !$gpu

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

  end subroutine derangmomy



  subroutine derangmomz(L,L_lo,L_hi,ncomp_L, &
                           u,u_lo,u_hi,ncomp_u, &
                           lo,hi,domlo,domhi, &
                           dx,xlo) bind(C, name="derangmomz")

    use amrex_constants_module, only: HALF
    use math_module, only: cross_product ! function
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: ncomp_L, ncomp_u ! == 1, 4
    integer, intent(in) :: L_lo(3), L_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt), intent(inout) :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt), intent(in) :: dx(3), xlo(3)

    integer          :: i, j, k
    real(rt)         :: loc(3), mom(3), ang_mom(3), rho

    !$gpu

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

  end subroutine derangmomz



  subroutine derpres(p,p_lo,p_hi,ncomp_p, &
                     u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                     domhi,dx,xlo) bind(C, name="derpres")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use amrex_constants_module, only : ONE
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in), value :: ncomp_p, ncomp_u
    integer, intent(in) :: p_lo(3), p_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt), intent(in) :: dx(3), xlo(3)

    real(rt) :: rhoInv
    integer  :: i, j, k

    type (eos_t) :: eos_state

    !$gpu

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

  end subroutine derpres



  subroutine dereint1(e,e_lo,e_hi,ncomp_e, &
                         u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                         domhi,dx,xlo) bind(C, name="dereint1")

    use amrex_constants_module, only : ONE, HALF
    use meth_params_module, only: URHO, UMX, UMY, UMZ, UEDEN 
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in), value :: ncomp_e, ncomp_u
    integer, intent(in) :: e_lo(3), e_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: e(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),ncomp_e)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt), intent(in) :: dx(3), xlo(3)

    real(rt)         :: rhoInv, ux, uy, uz
    integer          :: i, j, k

    !$gpu


    ! Compute internal energy from (rho E).

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

  end subroutine dereint1



  subroutine dereint2(e,e_lo,e_hi,ncomp_e, &
                         u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                         domhi,dx,xlo) bind(C, name="dereint2")

    use meth_params_module, only: URHO, UEINT
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: ncomp_e, ncomp_u
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: e_lo(3), e_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: e(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),ncomp_e)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt), intent(in) :: dx(3), xlo(3)

    integer          :: i, j, k

    !$gpu

    ! Compute internal energy from (rho e).

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             e(i,j,k,1) = u(i,j,k,UEINT) / u(i,j,k,URHO)
          enddo
       enddo
    enddo

  end subroutine dereint2



  subroutine dersoundspeed(c,c_lo,c_hi,ncomp_c, &
                              u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                              domhi,dx,xlo) bind(C, name="dersoundspeed")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use amrex_constants_module, only : ONE
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: ncomp_c, ncomp_u
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: c_lo(3), c_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: c(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3),ncomp_c)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt), intent(in) :: dx(3), xlo(3)

    real(rt)         :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

    !$gpu

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

  end subroutine dersoundspeed

  subroutine dergamma1(g1,g1_lo,g1_hi,ncomp_g1, &
                          u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                          domhi,dx,xlo) bind(C, name="dergamma1")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use amrex_constants_module, only : ONE
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: ncomp_g1, ncomp_u
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: g1_lo(3), g1_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: g1(g1_lo(1):g1_hi(1),g1_lo(2):g1_hi(2),g1_lo(3):g1_hi(3),ncomp_g1)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt), intent(in) :: dx(3), xlo(3)

    real(rt)         :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

    !$gpu

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

             g1(i,j,k,1) = eos_state % gam1
          enddo
       enddo
    enddo

  end subroutine dergamma1



  subroutine dermachnumber(mach,m_lo,m_hi,ncomp_mach, &
                              u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                              domhi,dx,xlo) bind(C, name="dermachnumber")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use meth_params_module, only: URHO, UMX, UMZ, UEINT, UTEMP, UFS, UFX
    use amrex_constants_module, only : ONE
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: ncomp_mach, ncomp_u
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: m_lo(3), m_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: mach(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),ncomp_mach)
    real(rt), intent(in) ::    u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u   )
    real(rt), intent(in) :: dx(3), xlo(3)

    real(rt)         :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

    !$gpu

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

  end subroutine dermachnumber



  subroutine derentropy(s,s_lo,s_hi,ncomp_s, &
                           u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                           domhi,dx,xlo) bind(C, name="derentropy")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use amrex_constants_module, only : ONE
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: ncomp_s, ncomp_u
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: s(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),ncomp_s)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
    real(rt), intent(in) :: dx(3), xlo(3)

    real(rt)         :: rhoInv
    integer          :: i, j, k

    type (eos_t) :: eos_state

    !$gpu

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

  end subroutine derentropy



  subroutine derenuctimescale(t,t_lo,t_hi,ncomp_t, &
                                 u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                                 domhi,dx,xlo) bind(C, name="derenuctimescale")

    use amrex_constants_module, only: ZERO, ONE
    use meth_params_module, only: URHO, UEINT, UTEMP, UFS, UFX
    use network, only: nspec, naux
    use prob_params_module, only: dim
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use amrex_fort_module, only : rt => amrex_real

    implicit none
    
    integer, intent(in), value :: ncomp_t, ncomp_u
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: t_lo(3), t_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: t(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3),ncomp_t)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u) ! NVAR, enuc
    real(rt), intent(in) :: dx(3), xlo(3)

    integer          :: i, j, k
    real(rt)         :: rhoInv, eint, enuc, t_s, t_e

    type (eos_t)     :: eos_state

    !$gpu

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

  end subroutine derenuctimescale



  subroutine derspec(spec,s_lo,s_hi,nv, &
                        dat,d_lo,d_hi,nc,lo,hi,domlo, &
                        domhi,delta,xlo) bind(C, name="derspec")
    !
    ! This routine derives the mass fractions of the species.
    !
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: nv, nc
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3)
    real(rt), intent(inout) :: spec(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nv)
    real(rt), intent(in) ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             spec(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
          end do
       end do
    end do

  end subroutine derspec



  subroutine derlogden(logden,l_lo,l_hi,nd, &
                          dat,d_lo,d_hi,nc, &
                          lo,hi,domlo,domhi,delta, &
                          xlo) bind(C, name="derlogden")

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: nd, nc
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: l_lo(3), l_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3)
    real(rt), intent(inout) :: logden(l_lo(1):l_hi(1),l_lo(2):l_hi(2),l_lo(3):l_hi(3),nd)
    real(rt), intent(in) ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             logden(i,j,k,1) = dlog10(dat(i,j,k,1))
          end do
       end do
    end do

  end subroutine derlogden



  subroutine dermagvort(vort,v_lo,v_hi,nv, & 
                           dat,d_lo,d_hi,nc,lo,hi,domlo, &
                           domhi,delta,xlo) bind(C, name="dermagvort")
    
    !
    ! This routine will calculate vorticity
    !     

    use amrex_constants_module, only : ZERO, HALF
    use prob_params_module, only: dg
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: nv, nc
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: v_lo(3), v_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3)
    real(rt), intent(inout) :: vort(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt), intent(in) ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k
    real(rt)         :: uy, uz, vx, vz, wx, wy, v1, v2, v3

    !$gpu

    uy = ZERO
    uz = ZERO
    vx = ZERO
    vz = ZERO
    wx = ZERO
    wy = ZERO

    !
    ! Calculate vorticity.
    !
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             vx = HALF * (dat(i+1*dg(1),j,k,3) / dat(i+1*dg(1),j,k,1) - dat(i-1*dg(1),j,k,3) / dat(i-1*dg(1),j,k,1)) / delta(1)
             wx = HALF * (dat(i+1*dg(1),j,k,4) / dat(i+1*dg(1),j,k,1) - dat(i-1*dg(1),j,k,4) / dat(i-1*dg(1),j,k,1)) / delta(1)

             if (delta(2) > ZERO) then
                uy = HALF * (dat(i,j+1*dg(2),k,2) / dat(i,j+1*dg(2),k,1) - dat(i,j-1*dg(2),k,2) / dat(i,j-1*dg(2),k,1)) / delta(2)
                wy = HALF * (dat(i,j+1*dg(2),k,4) / dat(i,j+1*dg(2),k,1) - dat(i,j-1*dg(2),k,4) / dat(i,j-1*dg(2),k,1)) / delta(2)
             endif

             if (delta(3) > ZERO) then
                uz = HALF * (dat(i,j,k+1*dg(3),2) / dat(i,j,k+1*dg(3),1) - dat(i,j,k-1*dg(3),2) / dat(i,j,k-1*dg(3),1)) / delta(3)
                vz = HALF * (dat(i,j,k+1*dg(3),3) / dat(i,j,k+1*dg(3),1) - dat(i,j,k-1*dg(3),3) / dat(i,j,k-1*dg(3),1)) / delta(3)
             endif

             v1 = wy - vz
             v2 = uz - wx
             v3 = vx - uy
             vort(i,j,k,1) = sqrt(v1*v1 + v2*v2 + v3*v3)

          end do
       end do
    end do

  end subroutine dermagvort



  subroutine derdivu(divu,u_lo,u_hi,nd, &
                        dat,d_lo,d_hi,nc, &
                        lo,hi,domlo,domhi,delta, &
                        xlo) bind(C, name="derdivu")
    !
    ! This routine will calculate the divergence of velocity.
    !

    use amrex_constants_module, only : ZERO, HALF
    use prob_params_module, only: dg
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: nd, nc
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3)
    real(rt), intent(inout) :: divu(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt), intent(in) ::  dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k
    real(rt)         :: ulo, uhi, vlo, vhi, wlo, whi

    !$gpu

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

  end subroutine derdivu



  subroutine derkineng(kineng,k_lo,k_hi,nk, &
                          dat,d_lo,d_hi,nc, &
                          lo,hi,domlo,domhi,delta, &
                          xlo) bind(C, name="derkineng")
    !
    ! This routine will derive kinetic energy = 1/2 rho (u^2 + v^2 + w^2)
    !

    use amrex_constants_module, only : HALF
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: nk, nc
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: k_lo(3), k_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3)
    real(rt), intent(inout) :: kineng(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3),nk)
    real(rt), intent(in) ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

    integer          :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             kineng(i,j,k,1) = HALF / dat(i,j,k,1) * ( dat(i,j,k,2)**2 + &
                  dat(i,j,k,3)**2 + &
                  dat(i,j,k,4)**2 )
          end do
       end do
    end do

  end subroutine derkineng



#ifdef DIFFUSION
  subroutine dercond(cond,u_lo,u_hi,nd, &
                        state,d_lo,d_hi,nc, &
                        lo,hi,domlo,domhi,delta, &
                        xlo) bind(C, name="dercond")
    !
    ! This routine will calculate the thermal conductivity
    !

    use amrex_constants_module, only : ZERO
    use meth_params_module, only: diffuse_cutoff_density, &
                                  URHO, UEINT, UTEMP, UFS, UFX
    use eos_type_module, only: eos_input_re, eos_t
    use eos_module, only: eos
    use network, only : nspec, naux
    use conductivity_module, only : conductivity
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: nd, nc
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3)
    real(rt), intent(inout) :: cond(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt), intent(in) :: state(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

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
                call conductivity(eos_state)
             else
                eos_state % conductivity = ZERO
             endif

             cond(i,j,k,1) = eos_state % conductivity

          enddo
       enddo
    enddo

  end subroutine dercond


  subroutine derdiffcoeff(diff,u_lo,u_hi,nd, &
                             state,d_lo,d_hi,nc, &
                             lo,hi,domlo,domhi,delta, &
                             xlo) bind(C, name="derdiffcoeff")
    !
    ! This routine will calculate the thermal conductivity
    !

    use amrex_constants_module, only : ZERO
    use meth_params_module, only: diffuse_cutoff_density, &
                                  URHO, UEINT, UTEMP, UFS, UFX
    use eos_type_module, only: eos_input_re, eos_t
    use eos_module, only: eos
    use network, only : nspec, naux
    use conductivity_module, only : conductivity
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: nd, nc
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3)
    real(rt), intent(inout) :: diff(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt), intent(in) :: state(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

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
                call conductivity(eos_state)
             else
                eos_state % conductivity = ZERO
             endif

             diff(i,j,k,1) = eos_state % conductivity/(eos_state%rho * eos_state%cv)

          enddo
       enddo
    enddo

  end subroutine derdiffcoeff


  subroutine derdiffterm(diff,u_lo,u_hi,nd, &
                            state,d_lo,d_hi,nc, &
                            lo,hi,domlo,domhi,delta, &
                            xlo) bind(C, name="derdiffterm")
    !
    ! This routine will calculate the thermal conductivity
    !

    use meth_params_module, only: UTEMP
    use prob_params_module, only: dim
    use diffusion_module, only : ca_fill_temp_cond
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: nd, nc
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: d_lo(3), d_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: delta(3), xlo(3)
    real(rt), intent(inout) :: diff(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt), intent(in) :: state(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)

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
    
  end subroutine derdiffterm
  
#endif

end module derive_module
