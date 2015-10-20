! construct the timestep limit for diffusion -- this is
! dt < 0.5 dx**2 / D
! where D = k/(rho c_v), and k is the conductivity

subroutine ca_estdt_diffusion(lo,hi,state,s_lo,s_hi,dx,dt)

  use network, only: nspec, naux
  use eos_module
  use eos_type_module
  use meth_params_module, only: NVAR, URHO, UEINT, UTEMP, UFS, UFX, &
       allow_negative_energy
  use prob_params_module, only: dim
  use bl_constants_module
  use conductivity_module

  implicit none

  integer          :: lo(3), hi(3)
  integer          :: s_lo(3), s_hi(3)
  double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
  double precision :: dx(3), dt

  double precision :: dt1, dt2, dt3, rho_inv
  integer          :: i, j, k, n
  double precision :: cond, D

  type (eos_t) :: eos_state

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           rho_inv = ONE/state(i,j,k,URHO)

           ! we need cv
           eos_state % rho = state(i,j,k,URHO )
           eos_state % T   = state(i,j,k,UTEMP)
           eos_state % e   = state(i,j,k,UEINT) * rho_inv

           eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rho_inv
           eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rho_inv

           call eos(eos_input_re, eos_state)

           ! we also need the conductivity
           call thermal_conductivity(eos_state, cond)

           ! maybe we should check (and take action) on negative cv here?
           D = cond*rho_inv/eos_state%cv
           
           dt1 = HALF*dx(1)**2/D

           if (dim >= 2) then
              dt2 = HALF*dx(2)**2/D
           else
              dt2 = dt1
           endif

           if (dim == 3) then
              dt3 = HALF*dx(3)**2/D
           else
              dt3 = dt1
           endif

           dt  = min(dt,dt1,dt2,dt3)

        enddo
     enddo
  enddo

end subroutine ca_estdt_diffusion
