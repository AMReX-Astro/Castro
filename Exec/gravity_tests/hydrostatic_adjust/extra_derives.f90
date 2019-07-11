subroutine ca_derpi(p,p_l1,p_h1,ncomp_p, &
                    u,u_l1,u_h1,ncomp_u,lo,hi,domlo, &
                    domhi,dx,xlo,time,dt,bc,level,grid_no) bind(C, name="ca_derpi")

  use network, only : nspec, naux
  use eos_module
  use eos_type_module
  use meth_params_module, only : URHO, UMX, UEINT, UTEMP, UFS, UFX
  use model_parser_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer p_l1,p_h1,ncomp_p
  integer u_l1,u_h1,ncomp_u
  integer lo(1), hi(1), domlo(1), domhi(1)
  real(rt)         p(p_l1:p_h1,ncomp_p)
  real(rt)         u(u_l1:u_h1,ncomp_u)
  real(rt)         dx(1), xlo(1), time, dt
  integer bc(1,2,ncomp_u), level, grid_no

  real(rt)         :: e, temp, X(nspec+naux)
  real(rt)         :: rhoInv
  integer          :: i,n

  type (eos_t) :: eos_state

  !     Compute pressure from the EOS
  do i = lo(1),hi(1)

     rhoInv = 1.e0_rt/u(i,URHO)
     e  = u(i,UEINT)*rhoInv
     temp  = u(i,UTEMP)
     do n=1,nspec
        X(n) = u(i,UFS+n-1)*rhoInv
     enddo
     do n=1,naux
        X(nspec+n) = u(i,UFX+n-1)*rhoInv
     enddo

     ! Protect against negative internal energy
     if (e .le. 0.e0_rt) then
        eos_state%rho = u(i,URHO)
        eos_state%T = temp
        eos_state%xn(:) = X

        call eos(eos_input_rt, eos_state)

        p(i,1) = eos_state%p
     else
        eos_state%rho = u(i,URHO)
        eos_state%T = temp
        eos_state%xn(:) = X
        eos_state%e = e

        call eos(eos_input_re, eos_state)

        p(i,1) = eos_state%p

     end if

     p(i,1) = p(i,1) - model_state(i, ipres_model)

  enddo

end subroutine ca_derpi

!-----------------------------------------------------------------------

subroutine ca_derpioverp0(p,p_l1,p_h1,ncomp_p, &
                          u,u_l1,u_h1,ncomp_u,lo,hi,domlo, &
                          domhi,dx,xlo,time,dt,bc,level,grid_no) &
                          bind(C, name="ca_derpioverp0")

  use network, only : nspec, naux
  use eos_module
  use eos_type_module
  use meth_params_module, only : URHO, UMX, UEINT, UTEMP, UFS, UFX
  use model_parser_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer p_l1,p_h1,ncomp_p
  integer u_l1,u_h1,ncomp_u
  integer lo(1), hi(1), domlo(1), domhi(1)
  real(rt)         p(p_l1:p_h1,ncomp_p)
  real(rt)         u(u_l1:u_h1,ncomp_u)
  real(rt)         dx(1), xlo(1), time, dt
  integer bc(1,2,ncomp_u), level, grid_no

  real(rt)         :: e, temp, X(nspec+naux)
  real(rt)         :: rhoInv
  integer          :: i,n

  type (eos_t) :: eos_state

  !     Compute pressure from the EOS
  do i = lo(1),hi(1)

     rhoInv = 1.e0_rt/u(i,URHO)
     e  = u(i,UEINT)*rhoInv
     temp  = u(i,UTEMP)
     do n=1,nspec
        X(n) = u(i,UFS+n-1)*rhoInv
     enddo
     do n=1,naux
        X(nspec+n) = u(i,UFX+n-1)*rhoInv
     enddo

     ! Protect against negative internal energy
     if (e .le. 0.e0_rt) then
        eos_state%rho = u(i,URHO)
        eos_state%T = temp
        eos_state%xn(:) = X

        call eos(eos_input_rt, eos_state)

        p(i,1) = eos_state%p
     else
        eos_state%rho = u(i,URHO)
        eos_state%T = temp
        eos_state%xn(:) = X
        eos_state%e = e

        call eos(eos_input_re, eos_state)

        p(i,1) = eos_state%p

     end if

     p(i,1) = abs(p(i,1) - model_state(i, ipres_model)) / model_state(i, ipres_model)

  enddo

end subroutine ca_derpioverp0
