module derive_module
    ! All subroutines in this file must be threadsafe because they are called
    ! inside OpenMP parallel regions.

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  public

contains


  subroutine deruplusc(vel,v_lo,v_hi,nv, &
                          dat,d_lo,d_hi,nc,lo,hi,domlo, &
                          domhi,delta) bind(C, name="deruplusc")

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
    real(rt), intent(in) :: delta(3)
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
                           domhi,delta) bind(C, name="deruminusc")

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
    real(rt), intent(in)    :: delta(3)
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


  subroutine derpres(p,p_lo,p_hi,ncomp_p, &
                     u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                     domhi,dx) bind(C, name="derpres")

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
    real(rt), intent(in) :: dx(3)

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


  subroutine dersoundspeed(c,c_lo,c_hi,ncomp_c, &
                              u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                              domhi,dx) bind(C, name="dersoundspeed")

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
    real(rt), intent(in) :: dx(3)

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
                          domhi,dx) bind(C, name="dergamma1")

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
    real(rt), intent(in) :: dx(3)

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
                              domhi,dx) bind(C, name="dermachnumber")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_input_re, eos_t
    use meth_params_module, only: URHO, UMX, UMZ, UEINT, UTEMP, UFS, UFX
    use amrex_constants_module, only : ONE, HALF
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in), value :: ncomp_mach, ncomp_u
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: m_lo(3), m_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: mach(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3),ncomp_mach)
    real(rt), intent(in) ::    u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u   )
    real(rt), intent(in) :: dx(3)

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

             mach(i,j,k,1) = sum(u(i,j,k,UMX:UMZ)**2)**HALF / u(i,j,k,URHO) / eos_state % cs
          enddo
       enddo
    enddo

  end subroutine dermachnumber



  subroutine derentropy(s,s_lo,s_hi,ncomp_s, &
                           u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                           domhi,dx) bind(C, name="derentropy")

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
    real(rt), intent(in) :: dx(3)

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
                                 domhi,dx) bind(C, name="derenuctimescale")

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
    real(rt), intent(in) :: dx(3)

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


#ifdef DIFFUSION
  subroutine dercond(lo, hi, &
                     cond, u_lo, u_hi, nd, &
                     state, d_lo, d_hi, nc, &
                     domlo, domhi, delta) bind(C, name="dercond")
    !
    ! This routine will calculate the thermal conductivity
    !

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO
    use meth_params_module, only: diffuse_cutoff_density, diffuse_cutoff_density_hi, &
                                  URHO, UEINT, UTEMP, UFS, UFX
    use eos_type_module, only: eos_input_re, eos_t
    use eos_module, only: eos
    use network, only: nspec, naux
    use conductivity_module, only: conductivity

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3)
    real(rt), intent(inout) :: cond(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt), intent(in   ) :: state(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer,  intent(in   ), value :: nd, nc

    integer :: i, j, k

    type(eos_t) :: eos_state
    real(rt) :: multiplier

    !$gpu

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

                if (eos_state%rho < diffuse_cutoff_density_hi) then
                    multiplier = (eos_state%rho - diffuse_cutoff_density) / &
                            (diffuse_cutoff_density_hi - diffuse_cutoff_density)
                    eos_state % conductivity = eos_state % conductivity * multiplier
                endif
             else
                eos_state % conductivity = ZERO
             endif

             cond(i,j,k,1) = eos_state % conductivity

          enddo
       enddo
    enddo

  end subroutine dercond


  subroutine derdiffcoeff(lo, hi, &
                          diff, u_lo, u_hi, nd, &
                          state, d_lo, d_hi, nc, &
                          domlo, domhi, delta) bind(C, name="derdiffcoeff")
    !
    ! This routine will calculate the thermal conductivity
    !

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO
    use meth_params_module, only: diffuse_cutoff_density, diffuse_cutoff_density_hi, &
                                  URHO, UEINT, UTEMP, UFS, UFX
    use eos_type_module, only: eos_input_re, eos_t
    use eos_module, only: eos
    use network, only: nspec, naux
    use conductivity_module, only: conductivity

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3)
    real(rt), intent(inout) :: diff(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt), intent(in   ) :: state(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer,  intent(in   ), value :: nd, nc

    integer :: i, j, k

    type(eos_t) :: eos_state
    real(rt) :: multiplier

    !$gpu

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

                if (eos_state%rho < diffuse_cutoff_density_hi) then
                    multiplier = (eos_state%rho - diffuse_cutoff_density) / &
                            (diffuse_cutoff_density_hi - diffuse_cutoff_density)
                    eos_state % conductivity = eos_state % conductivity * multiplier
                endif
             else
                eos_state % conductivity = ZERO
             endif

             diff(i,j,k,1) = eos_state % conductivity/(eos_state%rho * eos_state%cv)

          enddo
       enddo
    enddo

  end subroutine derdiffcoeff


  subroutine derdiffterm(lo, hi, &
                         diff, u_lo, u_hi, nd, &
                         state, d_lo, d_hi, nc, &
                         coeff_x, x_lo, x_hi, &
#if AMREX_SPACEDIM >= 2
                         coeff_y, y_lo, y_hi, &
#endif
#if AMREX_SPACEDIM == 3
                         coeff_z, z_lo, z_hi, &
#endif
                         domlo, domhi, delta) bind(C, name="derdiffterm")
    !
    ! This routine will calculate the thermal conductivity
    !

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, ONE
    use meth_params_module, only: UTEMP
    use prob_params_module, only: problo, coord_type

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    integer,  intent(in   ) :: x_lo(3), x_hi(3)
#if AMREX_SPACEDIM >= 2
    integer,  intent(in   ) :: y_lo(3), y_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer,  intent(in   ) :: z_lo(3), z_hi(3)
#endif
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3)
    real(rt), intent(inout) :: diff(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nd)
    real(rt), intent(in   ) :: state(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    real(rt), intent(in   ) :: coeff_x(x_lo(1):x_hi(1),x_lo(2):x_hi(2),x_lo(3):x_hi(3))
#if AMREX_SPACEDIM >= 2
    real(rt), intent(in   ) :: coeff_y(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(in   ) :: coeff_z(z_lo(1):z_hi(1),z_lo(2):z_hi(2),z_lo(3):z_hi(3))
#endif
    integer,  intent(in   ), value :: nd, nc

    real(rt) :: diff_term
    real(rt) :: kgradT_xhi, kgradT_xlo, kgradT_yhi, kgradT_ylo, kgradT_zhi, kgradT_zlo
    integer  :: i, j, k

    real(rt) :: r, rp1, rm1

    !$gpu

    ! create the diff term
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             kgradT_xhi = coeff_x(i+1,j,k)*(state(i+1,j,k,UTEMP) - state(i,  j,k,UTEMP))/delta(1)
             kgradT_xlo = coeff_x(i  ,j,k)*(state(i  ,j,k,UTEMP) - state(i-1,j,k,UTEMP))/delta(1)
#if AMREX_SPACEDIM >= 2
             kgradT_yhi = coeff_y(i,j+1,k)*(state(i,j+1,k,UTEMP) - state(i,j  ,k,UTEMP))/delta(2)
             kgradT_ylo = coeff_y(i,j  ,k)*(state(i,j,  k,UTEMP) - state(i,j-1,k,UTEMP))/delta(2)
#endif
#if AMREX_SPACEDIM == 3
             kgradT_zhi = coeff_z(i,j,k+1)*(state(i,j,k+1,UTEMP) - state(i,j,k  ,UTEMP))/delta(3)
             kgradT_zlo = coeff_z(i,j,k  )*(state(i,j,k  ,UTEMP) - state(i,j,k-1,UTEMP))/delta(3)
#endif

             if (coord_type == 0) then
                diff_term = (kgradT_xhi - kgradT_xlo)/delta(1)
#if AMREX_SPACEDIM >= 2
                diff_term = diff_term + (kgradT_yhi - kgradT_ylo)/delta(2)
#endif
#if AMREX_SPACEDIM == 3
                diff_term = diff_term + (kgradT_zhi - kgradT_zlo)/delta(3)
#endif

             else if (coord_type == 1) then
                ! axisymmetric coords (2-d)
                r = dble(i + HALF)*delta(1) + problo(1)
                rm1 = dble(i - ONE + HALF)*delta(1) + problo(1)
                rp1 = dble(i + ONE + HALF)*delta(1) + problo(1)

                diff_term = (rp1*kgradT_xhi - rm1*kgradT_xlo)/(r*delta(1)) + &
                            (kgradT_yhi - kgradT_ylo)/delta(2)

             else if (coord_type == 2) then
                ! spherical coords (1-d)
                r = dble(i + HALF)*delta(1) + problo(1)
                rm1 = dble(i - ONE + HALF)*delta(1) + problo(1)
                rp1 = dble(i + ONE + HALF)*delta(1) + problo(1)

                diff_term = (rp1**2*kgradT_xhi - rm1**2*kgradT_xlo)/(r**2*delta(1))

             endif

             diff(i,j,k,1) = diff_term

          enddo
       enddo
    enddo

  end subroutine derdiffterm

#endif

end module derive_module
