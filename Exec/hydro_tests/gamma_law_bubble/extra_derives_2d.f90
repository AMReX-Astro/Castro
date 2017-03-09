!-----------------------------------------------------------------------

subroutine ca_derpi(p,p_l1,p_l2,p_h1,p_h2,ncomp_p, &
     u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
     domhi,dx,xlo,time,dt,bc,level,grid_no) bind(C, name="ca_derpi")

  use network, only : nspec, naux
  use eos_module
  use eos_type_module
  use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX, &
       allow_negative_energy, const_grav
  use probdata_module, only: pres_base, dens_base, do_isentropic
  use prob_params_module, only: center
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: p_l1,p_l2,p_h1,p_h2,ncomp_p
  integer          :: u_l1,u_l2,u_h1,u_h2,ncomp_u
  integer          :: lo(2), hi(2), domlo(2), domhi(2)
  real(rt)         :: p(p_l1:p_h1,p_l2:p_h2,ncomp_p)
  real(rt)         :: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
  real(rt)         :: dx(2), xlo(2), time, dt
  integer          :: bc(2,2,ncomp_u), level, grid_no

  real(rt)         :: e, T
  real(rt)         :: rhoInv
  integer          :: i,j,npts_1d
  real(rt)         H,z,xn(1), const
  real(rt)        , allocatable :: pressure(:), density(:), temp(:), eint(:)

  type (eos_t) :: eos_state      

  ! first make a 1D initial model for the entire domain
  npts_1d = (2.e0_rt*center(2)+1.e-8_rt) / dx(2)

  allocate(pressure(0:npts_1d-1))
  allocate(density (0:npts_1d-1))
  allocate(temp    (0:npts_1d-1))
  allocate(eint    (0:npts_1d-1))

  const = pres_base/dens_base**gamma_const

  pressure(0) = pres_base
  density(0)  = dens_base

  ! only initialize the first species
  xn(1) = 1.e0_rt

  ! compute the pressure scale height (for an isothermal, ideal-gas
  ! atmosphere)
  H = pres_base / dens_base / abs(const_grav)

  do j=0,npts_1d-1

     ! initial guess
     temp(j) = 1000.e0_rt

     if (do_isentropic) then
        z = dble(j) * dx(2)
        density(j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)*z/ &
             (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
     else
        z = (dble(j)+0.5e0_rt) * dx(2)
        density(j) = dens_base * exp(-z/H)
     end if

     if (j .gt. 0) then
        pressure(j) = pressure(j-1) - &
             dx(2) * 0.5e0_rt * (density(j)+density(j-1)) * abs(const_grav)
     end if

     eos_state%rho = density(j)
     eos_state%T = temp(j)
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)

     eint(j) = eos_state%e
     temp(j) = eos_state%T

  end do

  ! Compute pressure from the EOS
  do j = lo(2),hi(2)
     do i = lo(1),hi(1)
        rhoInv = 1.e0_rt/u(i,j,URHO)
        T = u(i,j,UTEMP)

        eos_state%rho = u(i,j,URHO)
        eos_state%T = u(i,j,UTEMP)
        eos_state%e = u(i,j,UEINT)*rhoInv
        eos_state%xn(:) = u(i,j,UFS:UFS-1+nspec)/u(i,j,URHO)
        eos_state%aux(:) = u(i,j,UFX:UFX-1+naux)/u(i,j,URHO)            

        ! Protect against negative internal energy
        if (allow_negative_energy .eq. 0 .and. e .le. 0.e0_rt) then
           call eos(eos_input_rt, eos_state)
           p(i,j,1) = eos_state%p

        else
           call eos(eos_input_re, eos_state)
           p(i,j,1) = eos_state%p

        end if

        p(i,j,1) = p(i,j,1) - pressure(j)

     enddo
  enddo

end subroutine ca_derpi

!-----------------------------------------------------------------------

subroutine ca_derpioverp0(p,p_l1,p_l2,p_h1,p_h2,ncomp_p, &
     u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
     domhi,dx,xlo,time,dt,bc,level,grid_no) bind(C, name="ca_derpioverp0")

  use network, only : nspec, naux
  use eos_module
  use eos_type_module
  use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX, &
       allow_negative_energy, const_grav
  use probdata_module, only: pres_base, dens_base, do_isentropic
  use prob_params_module, only: center
  
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: p_l1,p_l2,p_h1,p_h2,ncomp_p
  integer          :: u_l1,u_l2,u_h1,u_h2,ncomp_u
  integer          :: lo(2), hi(2), domlo(2), domhi(2)
  real(rt)         :: p(p_l1:p_h1,p_l2:p_h2,ncomp_p)
  real(rt)         :: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
  real(rt)         :: dx(2), xlo(2), time, dt
  integer          :: bc(2,2,ncomp_u), level, grid_no

  real(rt)         :: e, T
  real(rt)         :: rhoInv
  integer          :: i,j,npts_1d
  real(rt)         H,z,xn(1), const
  real(rt)        , allocatable :: pressure(:), density(:), temp(:), eint(:)

  type (eos_t) :: eos_state      

  ! first make a 1D initial model for the entire domain
  npts_1d = (2.e0_rt*center(2)+1.e-8_rt) / dx(2)

  allocate(pressure(0:npts_1d-1))
  allocate(density (0:npts_1d-1))
  allocate(temp    (0:npts_1d-1))
  allocate(eint    (0:npts_1d-1))

  const = pres_base/dens_base**gamma_const

  pressure(0) = pres_base
  density(0)  = dens_base

  ! only initialize the first species
  xn(1) = 1.e0_rt

  ! compute the pressure scale height (for an isothermal, ideal-gas
  ! atmosphere)
  H = pres_base / dens_base / abs(const_grav)

  do j=0,npts_1d-1

     ! initial guess
     temp(j) = 1000.e0_rt

     if (do_isentropic) then
        z = dble(j) * dx(2)
        density(j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)*z/ &
             (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
     else
        z = (dble(j)+0.5e0_rt) * dx(2)
        density(j) = dens_base * exp(-z/H)
     end if

     if (j .gt. 0) then
        pressure(j) = pressure(j-1) - &
             dx(2) * 0.5e0_rt * (density(j)+density(j-1)) * abs(const_grav)
     end if

     eos_state%rho = density(j)
     eos_state%T = temp(j)
     eos_state%xn(:) = xn(:)
     eos_state%p = pressure(j)

     call eos(eos_input_rp, eos_state)

     eint(j) = eos_state%e

  end do

  ! Compute pressure from the EOS
  do j = lo(2),hi(2)
     do i = lo(1),hi(1)
        rhoInv = 1.e0_rt/u(i,j,URHO)
        e = u(i,j,UEINT)*rhoInv
        T = u(i,j,UTEMP)

        eos_state%rho = u(i,j,URHO)
        eos_state%T = T
        eos_state%xn(:) = u(i,j,UFS:UFS-1+nspec)*rhoInv
        eos_state%aux(:) = u(i,j,UFX:UFX-1+naux)*rhoInv    
        eos_state%e = e

        ! Protect against negative internal energy
        if (allow_negative_energy .eq. 0 .and. e .le. 0.e0_rt) then
           call eos(eos_input_rt, eos_state)
           p(i,j,1) = eos_state%p

        else
           call eos(eos_input_re, eos_state)
           p(i,j,1) = eos_state%p

        end if

        p(i,j,1) = (p(i,j,1) - pressure(j)) / pressure(j)

     enddo
  enddo

end subroutine ca_derpioverp0

!-----------------------------------------------------------------------

subroutine ca_derrhopert(p,p_l1,p_l2,p_h1,p_h2,ncomp_p, &
     u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
     domhi,dx,xlo,time,dt,bc,level,grid_no) bind(C, name="ca_derrhopert")

  use network, only : nspec, naux
  use meth_params_module, only : URHO, const_grav
  use eos_module, only: gamma_const
  use probdata_module
  use interpolate_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: p_l1,p_l2,p_h1,p_h2,ncomp_p
  integer          :: u_l1,u_l2,u_h1,u_h2,ncomp_u
  integer          :: lo(2), hi(2), domlo(2), domhi(2)
  real(rt)         :: p(p_l1:p_h1,p_l2:p_h2,ncomp_p)
  real(rt)         :: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
  real(rt)         :: dx(2), xlo(2), time, dt
  integer          :: bc(2,2,ncomp_u), level, grid_no

  ! local
  integer :: i,j

  real(rt)         :: y,dens,H

  do j=lo(2),hi(2)

     if (do_isentropic) then
        y = xlo(2) + dx(2)*float(j-lo(2))
        dens = dens_base*(const_grav*dens_base*(gamma_const - 1.0)*y/ &
             (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
     else
        y = xlo(2) + dx(2)*(float(j-lo(2)) + 0.5e0_rt)
        dens = dens_base * exp(-y/H)
     end if

     do i=lo(1),hi(1)
        p(i,j,1) = u(i,j,URHO) - dens
     end do

  end do

end subroutine ca_derrhopert

!-----------------------------------------------------------------------

subroutine ca_dertpert(p,p_l1,p_l2,p_h1,p_h2,ncomp_p, &
     u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
     domhi,dx,xlo,time,dt,bc,level,grid_no) bind(C, name="ca_dertpert")

  use network, only : nspec, naux
  use eos_module
  use eos_type_module
  use meth_params_module, only : UTEMP, const_grav
  use probdata_module, only: pres_base, dens_base, do_isentropic
  use prob_params_module, only: center
  
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: p_l1,p_l2,p_h1,p_h2,ncomp_p
  integer          :: u_l1,u_l2,u_h1,u_h2,ncomp_u
  integer          :: lo(2), hi(2), domlo(2), domhi(2)
  real(rt)         :: p(p_l1:p_h1,p_l2:p_h2,ncomp_p)
  real(rt)         :: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
  real(rt)         :: dx(2), xlo(2), time, dt
  integer          :: bc(2,2,ncomp_u), level, grid_no

  integer          :: i,j,npts_1d
  real(rt)         H,z,xn(1), const
  real(rt)        , allocatable :: pressure(:), density(:), temp(:), eint(:)

  type (eos_t) :: eos_state      

  ! first make a 1D initial model for the entire domain
  npts_1d = (2.e0_rt*center(2)+1.e-8_rt) / dx(2)

  allocate(pressure(0:npts_1d-1))
  allocate(density (0:npts_1d-1))
  allocate(temp    (0:npts_1d-1))
  allocate(eint    (0:npts_1d-1))

  const = pres_base/dens_base**gamma_const

  pressure(0) = pres_base
  density(0)  = dens_base

  ! only initialize the first species
  xn(1) = 1.e0_rt

  ! compute the pressure scale height (for an isothermal, ideal-gas
  ! atmosphere)
  H = pres_base / dens_base / abs(const_grav)

  do j=0,npts_1d-1

     ! initial guess
     temp(j) = 1000.e0_rt

     if (do_isentropic) then
        z = dble(j) * dx(2)
        density(j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)*z/ &
             (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
     else
        z = (dble(j)+0.5e0_rt) * dx(2)
        density(j) = dens_base * exp(-z/H)
     end if

     if (j .gt. 0) then
        pressure(j) = pressure(j-1) - &
             dx(2) * 0.5e0_rt * (density(j)+density(j-1)) * abs(const_grav)
     end if

     eos_state%rho = density(j)
     eos_state%T = temp(j)
     eos_state%xn(:) = xn(:)
     eos_state%p = pressure(j)

     call eos(eos_input_rp, eos_state)

     eint(j) = eos_state%e
     temp(j) = eos_state%T

  end do

  do j = lo(2),hi(2)
     do i = lo(1),hi(1)

        p(i,j,1) = u(i,j,UTEMP) - temp(j)

     enddo
  enddo

end subroutine ca_dertpert

