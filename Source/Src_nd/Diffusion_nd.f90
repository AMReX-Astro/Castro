module diffusion_module

  implicit none

  public

contains

  subroutine ca_tempdiffextrap(lo, hi, tdif, t_lo, t_hi) &
       bind(C, name="ca_tempdiffextrap")

    use prob_params_module, only: dg

    implicit none

    integer :: lo(3), hi(3)
    integer :: t_lo(3), t_hi(3)
    double precision :: tdif(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))

    ! Local variables

    integer :: i, j, k

    ! left side
    i = lo(1)-1*dg(1)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          tdif(i,j,k) = tdif(i+1*dg(1),j,k)
       end do
    end do

    ! right side
    i = hi(1)+1*dg(1)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          tdif(i,j,k) = tdif(i-1*dg(1),j,k)
       end do
    end do

    ! bottom side
    j = lo(2)-1*dg(2)
    do k = lo(3), hi(3)
       do i = lo(1), hi(1)
          tdif(i,j,k) = tdif(i,j+1*dg(2),k)
       end do
    end do

    ! top side
    j = hi(2)+1*dg(2)
    do k = lo(3), hi(3)
       do i = lo(1), hi(1)
          tdif(i,j,k) = tdif(i,j-1*dg(2),k)
       end do
    end do

    ! down side
    k = lo(3)-1*dg(3)
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          tdif(i,j,k) = tdif(i,j,k+1*dg(3))
       end do
    end do

    ! up side
    k = hi(3)+1*dg(3)
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          tdif(i,j,k) = tdif(i,j,k-1*dg(3))
       end do
    end do

    ! k-edges
    i = lo(1)-1*dg(1)
    j = lo(2)-1*dg(2)
    do k = lo(3), hi(3)
       tdif(i,j,k) = tdif(i+1*dg(1),j+1*dg(2),k)
    end do

    i = lo(1)-1*dg(1)
    j = hi(2)+1*dg(2)
    do k = lo(3), hi(3)
       tdif(i,j,k) = tdif(i+1*dg(1),j-1*dg(2),k)
    end do

    i = hi(1)+1*dg(1)
    j = lo(2)-1*dg(2)
    do k = lo(3), hi(3)
       tdif(i,j,k) = tdif(i-1*dg(1),j+1*dg(2),k)
    end do

    i = hi(1)+1*dg(1)
    j = hi(2)+1*dg(2)
    do k = lo(3), hi(3)
       tdif(i,j,k) = tdif(i-1*dg(1),j-1*dg(2),k)
    end do

    ! j-edges
    i = lo(1)-1*dg(1)
    k = lo(3)-1*dg(3)
    do j = lo(2), hi(2)
       tdif(i,j,k) = tdif(i+1*dg(1),j,k+1*dg(3))
    end do

    i = lo(1)-1*dg(1)
    k = hi(3)+1*dg(3)
    do j = lo(2), hi(2)
       tdif(i,j,k) = tdif(i+1*dg(1),j,k-1*dg(3))
    end do

    i = hi(1)+1*dg(1)
    k = lo(3)-1*dg(3)
    do j = lo(2), hi(2)
       tdif(i,j,k) = tdif(i-1*dg(1),j,k+1*dg(3))
    end do

    i = hi(1)+1*dg(1)
    k = hi(3)+1*dg(3)
    do j = lo(2), hi(2)
       tdif(i,j,k) = tdif(i-1*dg(1),j,k-1*dg(3))
    end do

    ! i-edges
    j = lo(2)-1*dg(2)
    k = lo(3)-1*dg(3)
    do i = lo(1), hi(1)
       tdif(i,j,k) = tdif(i,j+1*dg(2),k+1*dg(3))
    end do

    j = lo(2)-1*dg(2)
    k = hi(3)+1*dg(3)
    do i = lo(1), hi(1)
       tdif(i,j,k) = tdif(i,j+1*dg(2),k-1*dg(3))
    end do

    j = hi(2)+1*dg(2)
    k = lo(3)-1*dg(3)
    do i = lo(1), hi(1)
       tdif(i,j,k) = tdif(i,j-1*dg(2),k+1*dg(3))
    end do

    j = hi(2)+1*dg(2)
    k = hi(3)+1*dg(3)
    do i = lo(1), hi(1)
       tdif(i,j,k) = tdif(i,j-1*dg(2),k-1*dg(3))
    end do

    ! corners
    i = lo(1)-1*dg(1)
    j = lo(2)-1*dg(2)
    k = lo(3)-1*dg(3)
    tdif(i,j,k) = tdif(i+1*dg(1),j+1*dg(2),k+1*dg(3))

    i = lo(1)-1*dg(1)
    j = hi(2)+1*dg(2)
    k = lo(3)-1*dg(3)
    tdif(i,j,k) = tdif(i+1*dg(1),j-1*dg(2),k+1*dg(3))

    i = hi(1)+1*dg(1)
    j = hi(2)+1*dg(2)
    k = lo(3)-1*dg(3)
    tdif(i,j,k) = tdif(i-1*dg(1),j-1*dg(2),k+1*dg(3))

    i = hi(1)+1*dg(1)
    j = lo(2)-1*dg(2)
    k = lo(3)-1*dg(3)
    tdif(i,j,k) = tdif(i-1*dg(1),j+1*dg(2),k+1*dg(3))

    i = lo(1)-1*dg(1)
    j = lo(2)-1*dg(2)
    k = hi(3)+1*dg(3)
    tdif(i,j,k) = tdif(i+1*dg(1),j+1*dg(2),k-1*dg(3))

    i = lo(1)-1*dg(1)
    j = hi(2)+1*dg(2)
    k = hi(3)+1*dg(3)
    tdif(i,j,k) = tdif(i+1*dg(1),j-1*dg(2),k-1*dg(3))

    i = hi(1)+1*dg(1)
    j = lo(2)-1*dg(2)
    k = hi(3)+1*dg(3)
    tdif(i,j,k) = tdif(i-1*dg(1),j+1*dg(2),k-1*dg(3))

    i = hi(1)+1*dg(1)
    j = hi(2)+1*dg(2)
    k = hi(3)+1*dg(3)
    tdif(i,j,k) = tdif(i-1*dg(1),j-1*dg(2),k-1*dg(3))

  end subroutine ca_tempdiffextrap


  
  ! This routine fills the species diffusion coefficients on the edges of a zone
  ! by calling the cell-centered conductivity routine and averaging to
  ! the interfaces

  subroutine ca_fill_spec_coeff(lo,hi, &
       state,s_lo,s_hi, &
       coefx,cx_lo,cx_hi, &
       coefy,cy_lo,cy_hi, &
       coefz,cz_lo,cz_hi, dx) &
       bind(C, name="ca_fill_spec_coeff")

    use bl_constants_module
    use network, only: nspec, naux
    use meth_params_module, only : NVAR, URHO, UTEMP, UFS, UFX, diffuse_cutoff_density
    use prob_params_module, only : dg
    use conductivity_module
    use eos_type_module

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    integer         , intent(in   ) :: cx_lo(3), cx_hi(3), cy_lo(3), cy_hi(3), cz_lo(3), cz_hi(3)
    real (kind=dp_t), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real (kind=dp_t), intent(inout) :: coefx(cx_lo(1):cx_hi(1),cx_lo(2):cx_hi(2),cx_lo(3):cx_hi(3))
    real (kind=dp_t), intent(inout) :: coefy(cy_lo(1):cy_hi(1),cy_lo(2):cy_hi(2),cy_lo(3):cy_hi(3))
    real (kind=dp_t), intent(inout) :: coefz(cz_lo(1):cz_hi(1),cz_lo(2):cz_hi(2),cz_lo(3):cz_hi(3))
    real (kind=dp_t), intent(in   ) :: dx(3)

    ! local variables
    integer          :: i, j, k
    double precision :: coef_cc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

    type (eos_t) :: eos_state
    double precision :: coeff

    ! fill the cell-centered diffusion coefficient

    do k = lo(3)-1*dg(3),hi(3)+1*dg(3)
       do j = lo(2)-1*dg(2),hi(2)+1*dg(2)
          do i = lo(1)-1*dg(1),hi(1)+1*dg(1)
             eos_state%rho    = state(i,j,k,URHO)
             eos_state%T      = state(i,j,k,UTEMP)
             eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)
             eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)

             if (eos_state%rho > diffuse_cutoff_density) then
                call thermal_conductivity(eos_state, coeff)
                coeff = coeff / eos_state%cp
             else
                coeff = ZERO
             endif

             coef_cc(i,j,k) = coeff
          enddo
       enddo
    enddo

    ! average to the interfaces
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1*dg(1)
             coefx(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i-1*dg(1),j,k))
          end do
       end do
    enddo

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1*dg(2)
          do i = lo(1),hi(1)
             coefy(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j-1*dg(2),k))
          end do
       end do
    enddo

    do k = lo(3),hi(3)+1*dg(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             coefz(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j,k-1*dg(3)))
          end do
       end do
    enddo

  end subroutine ca_fill_spec_coeff



  ! This routine fills the thermal conductivity on the edges of a zone
  ! by calling the cell-centered conductivity routine and averaging to
  ! the interfaces
  
  subroutine ca_fill_temp_cond(lo,hi, &
       state,s_lo,s_hi, &
       coefx,cx_lo,cx_hi, &
       coefy,cy_lo,cy_hi, &
       coefz,cz_lo,cz_hi, dx) &
       bind(C, name="ca_fill_temp_cond")

    use bl_constants_module
    use network, only: nspec, naux
    use meth_params_module, only : NVAR, URHO, UTEMP, UFS, UFX, diffuse_cutoff_density
    use prob_params_module, only : dg
    use conductivity_module
    use eos_type_module

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    integer         , intent(in   ) :: cx_lo(3), cx_hi(3), cy_lo(3), cy_hi(3), cz_lo(3), cz_hi(3)
    real (kind=dp_t), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real (kind=dp_t), intent(inout) :: coefx(cx_lo(1):cx_hi(1),cx_lo(2):cx_hi(2),cx_lo(3):cx_hi(3))
    real (kind=dp_t), intent(inout) :: coefy(cy_lo(1):cy_hi(1),cy_lo(2):cy_hi(2),cy_lo(3):cy_hi(3))
    real (kind=dp_t), intent(inout) :: coefz(cz_lo(1):cz_hi(1),cz_lo(2):cz_hi(2),cz_lo(3):cz_hi(3))
    real (kind=dp_t), intent(in   ) :: dx(3)

    ! local variables
    integer          :: i, j, k
    double precision :: coef_cc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

    type (eos_t) :: eos_state
    double precision :: cond

    ! fill the cell-centered conductivity

    do k = lo(3)-1*dg(3),hi(3)+1*dg(3)
       do j = lo(2)-1*dg(2),hi(2)+1*dg(2)
          do i = lo(1)-1*dg(1),hi(1)+1*dg(1)
             eos_state%rho    = state(i,j,k,URHO)
             eos_state%T      = state(i,j,k,UTEMP)
             eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)
             eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)

             if (eos_state%rho > diffuse_cutoff_density) then
                call thermal_conductivity(eos_state, cond)
             else
                cond = ZERO
             endif

             coef_cc(i,j,k) = cond
          enddo
       enddo
    enddo

    ! average to the interfaces
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1*dg(1)
             coefx(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i-1*dg(1),j,k))
          end do
       end do
    enddo

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1*dg(2)
          do i = lo(1),hi(1)
             coefy(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j-1*dg(2),k))
          end do
       end do
    enddo

    do k = lo(3),hi(3)+1*dg(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             coefz(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j,k-1*dg(3)))
          end do
       end do
    enddo

  end subroutine ca_fill_temp_cond



  ! This routine fills the viscous coefficients "mu" on the edges of a zone
  ! by calling the cell-centered coefficient routine and averaging to
  ! the interfaces
  
  subroutine ca_fill_first_visc_coeff(lo,hi, &
       state,s_lo,s_hi, &
       coefx,cx_lo,cx_hi, &
       coefy,cy_lo,cy_hi, &
       coefz,cz_lo,cz_hi, dx) bind(C, name="ca_fill_first_visc_coeff")

    use bl_constants_module
    use network, only: nspec, naux
    use meth_params_module, only : NVAR, URHO, UTEMP, UFS, UFX, diffuse_cutoff_density
    use prob_params_module, only : dg
    use viscosity_module
    use eos_type_module

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    integer         , intent(in   ) :: cx_lo(3), cx_hi(3), cy_lo(3), cy_hi(3), cz_lo(3), cz_hi(3)
    real (kind=dp_t), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real (kind=dp_t), intent(inout) :: coefx(cx_lo(1):cx_hi(1),cx_lo(2):cx_hi(2),cx_lo(3):cx_hi(3))
    real (kind=dp_t), intent(inout) :: coefy(cy_lo(1):cy_hi(1),cy_lo(2):cy_hi(2),cy_lo(3):cy_hi(3))
    real (kind=dp_t), intent(inout) :: coefz(cz_lo(1):cz_hi(1),cz_lo(2):cz_hi(2),cz_lo(3):cz_hi(3))
    real (kind=dp_t), intent(in   ) :: dx(3)

    ! local variables
    integer          :: i, j, k
    double precision :: coef_cc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

    type (eos_t) :: eos_state
    double precision :: coeff

    ! fill the cell-centered viscous coefficient

    do k = lo(3)-1*dg(3),hi(3)+1*dg(3)
       do j = lo(2)-1*dg(2),hi(2)+1*dg(2)
          do i = lo(1)-1*dg(1),hi(1)+1*dg(1)
             eos_state%rho    = state(i,j,k,URHO)
             eos_state%T      = state(i,j,k,UTEMP)
             eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)
             eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)

             if (eos_state%rho > diffuse_cutoff_density) then
                call viscous_coeff(eos_state, coeff)
             else
                coeff = ZERO
             endif
             coef_cc(i,j,k) = 2.d0 * coeff
          enddo
       enddo
    enddo

    ! average to the interfaces
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1*dg(1)
             coefx(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i-1*dg(1),j,k))
          end do
       end do
    enddo

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1*dg(2)
          do i = lo(1),hi(1)
             coefy(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j-1*dg(2),k))
          end do
       end do
    enddo

    do k = lo(3),hi(3)+1*dg(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             coefz(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j,k-1*dg(3)))
          end do
       end do
    enddo

  end subroutine ca_fill_first_visc_coeff


  
  subroutine ca_fill_secnd_visc_coeff(lo,hi, &
       state,s_lo,s_hi, &
       coefx,cx_lo,cx_hi, &
       coefy,cy_lo,cy_hi, &
       coefz,cz_lo,cz_hi, dx) bind(C, name="ca_fill_secnd_visc_coeff")

    use bl_constants_module
    use network, only: nspec, naux
    use meth_params_module, only : NVAR, URHO, UTEMP, UFS, UFX, diffuse_cutoff_density
    use prob_params_module, only : dg
    use viscosity_module
    use eos_type_module

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    integer         , intent(in   ) :: cx_lo(3), cx_hi(3), cy_lo(3), cy_hi(3), cz_lo(3), cz_hi(3)
    real (kind=dp_t), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real (kind=dp_t), intent(inout) :: coefx(cx_lo(1):cx_hi(1),cx_lo(2):cx_hi(2),cx_lo(3):cx_hi(3))
    real (kind=dp_t), intent(inout) :: coefy(cy_lo(1):cy_hi(1),cy_lo(2):cy_hi(2),cy_lo(3):cy_hi(3))
    real (kind=dp_t), intent(inout) :: coefz(cz_lo(1):cz_hi(1),cz_lo(2):cz_hi(2),cz_lo(3):cz_hi(3))
    real (kind=dp_t), intent(in   ) :: dx(3)

    ! local variables
    integer          :: i, j, k
    double precision :: coef_cc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

    type (eos_t) :: eos_state
    double precision :: bulk_visc, mu, twothirds

    bulk_visc = 0.d0
    twothirds = 2.d0 / 3.d0

    ! fill the cell-centered viscous coefficient

    do k = lo(3)-1*dg(3),hi(3)+1*dg(3)
       do j = lo(2)-1*dg(2),hi(2)+1*dg(2)
          do i = lo(1)-1*dg(1),hi(1)+1*dg(1)
             eos_state%rho    = state(i,j,k,URHO)
             eos_state%T      = state(i,j,k,UTEMP)
             eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)
             eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)

             if (eos_state%rho > diffuse_cutoff_density) then
                call viscous_coeff(eos_state, mu)
             else
                mu = ZERO
             endif
             !          coef_cc(i,j,k) = coeff

             coef_cc(i,j,k) = (bulk_visc - twothirds*mu)

          enddo
       enddo
    enddo

    ! average to the interfaces
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1*dg(1)
             coefx(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i-1*dg(1),j,k))
          end do
       end do
    enddo

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1*dg(2)
          do i = lo(1),hi(1)
             coefy(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j-1*dg(2),k))
          end do
       end do
    enddo

    do k = lo(3),hi(3)+1*dg(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             coefz(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j,k-1*dg(3)))
          end do
       end do
    enddo

  end subroutine ca_fill_secnd_visc_coeff


  
  subroutine ca_compute_div_tau_u(lo,hi,&
       div_tau_u,d_lo,d_hi, &
       state,s_lo,s_hi,dx,coord_type) bind(C, name="ca_compute_div_tau_u")

    use bl_constants_module
    use network, only: nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UTEMP, UFS, UFX
    use prob_params_module, only : dg
    use viscosity_module
    use eos_type_module

    implicit none

    integer         , intent(in   ) ::   lo(3),   hi(3)
    integer         , intent(in   ) :: d_lo(3), d_hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    real (kind=dp_t), intent(  out) :: div_tau_u(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    real (kind=dp_t), intent(in   ) :: state    (s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real (kind=dp_t), intent(in   ) :: dx(3)
    integer         , intent(in   ) :: coord_type

    ! local variables
    real (kind=dp_t) :: twothirds
    type (eos_t)     :: eos_state

    ! These are cell-centered
    real (kind=dp_t) :: rm,rc,rp,vm,vc,vp,mu_temp
    real (kind=dp_t), allocatable :: mu(:), kappa(:)

    ! These are edge-centered
    real (kind=dp_t) :: mu_em, kap_em, tau1m, tau2m, x_em, v_em
    real (kind=dp_t) :: mu_ep, kap_ep, tau1p, tau2p, x_ep, v_ep

    integer          :: i, j, k, imin

    twothirds = 2.d0 / 3.d0

    allocate(   mu(lo(1)-2:hi(1)+2))
    allocate(kappa(lo(1)-2:hi(1)+2))

    ! Fill cell-centered values of div(tau u) as a viscous source term for the energy equation

    ! NOTE THIS IS ONLY CORRECT IN 1-D SPHERICAL !!!!!!

    k = lo(3)
    j = lo(2)

    do i = lo(1)-2,hi(1)+2
       eos_state%rho    = state(i,j,k,URHO)
       eos_state%T      = state(i,j,k,UTEMP)
       eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)
       eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)

       call viscous_coeff(eos_state, mu_temp)
       mu(i) = mu_temp

       kappa(i) = 0.d0
    end do

    imin = max(lo(1)-1,0)

    if (coord_type .eq. 2) then

       do i = imin,hi(1)+1

          ! These are cell-centered
          rm = (dble(i)-0.5d0) * dx(1)
          rc = (dble(i)+0.5d0) * dx(1)
          rp = (dble(i)+1.5d0) * dx(1)

          ! These are cell-centered
          vp = state(i+1,j,k,UMX)/state(i+1,j,k,URHO)
          vc = state(i  ,j,k,UMX)/state(i  ,j,k,URHO)
          vm = state(i-1,j,k,UMX)/state(i-1,j,k,URHO)

          ! These are edge-centered
          mu_em = (mu(i-1)+mu(i)) * 0.5d0
          mu_ep = (mu(i+1)+mu(i)) * 0.5d0

          ! These are edge-centered
          kap_em = (kappa(i-1)+kappa(i)) * 0.5d0
          kap_ep = (kappa(i+1)+kappa(i)) * 0.5d0

          ! These are edge-centered
          tau1p = 2.d0 * mu_ep * (vp - vc) /dx(1)
          tau1m = 2.d0 * mu_em * (vc - vm) /dx(1)

          ! These are edge-centered
          tau2p = (kap_ep - twothirds*mu_ep) * (rp**2*vp - rc**2* vc) /dx(1)
          tau2m = (kap_em - twothirds*mu_em) * (rc**2*vc - rm**2* vm) /dx(1)

          ! These are edge-centered
          x_ep = dble(i+1)*dx(1)
          x_em = dble(i  )*dx(1)

          ! These are edge-centered
          v_ep = (vc + vp) * 0.5d0
          v_em = (vc + vm) * 0.5d0

          div_tau_u(i,j,k) = ( (x_ep**2*tau1p*v_ep) - (x_em**2*tau1m*v_em) &
               +(tau2p*v_ep)-(tau2m*v_em) ) / (rc**2 * dx(1))

       enddo

       if (lo(1) .eq. 0) &
            div_tau_u(lo(1)-1,j,k) = div_tau_u(lo(1),j,k)

    else

       call bl_abort("compute_div_tau_u currently hard-wired for 1-d spherical")

    end if

    deallocate(mu,kappa)

  end subroutine ca_compute_div_tau_u
  
end module diffusion_module
