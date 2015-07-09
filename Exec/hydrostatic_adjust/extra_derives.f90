subroutine ca_derpi(p,p_l1,p_h1,ncomp_p, &
                    u,u_l1,u_h1,ncomp_u,lo,hi,domlo, &
                    domhi,dx,xlo,time,dt,bc,level,grid_no)

  use network, only : nspec, naux
  use eos_module
  use eos_type_module
  use eos_data_module
  use meth_params_module, only : URHO, UMX, UEINT, UTEMP, UFS, UFX, &
                                 allow_negative_energy
  use probdata_module, only: hse_p

  implicit none
  
  integer p_l1,p_h1,ncomp_p
  integer u_l1,u_h1,ncomp_u
  integer lo(1), hi(1), domlo(1), domhi(1)
  double precision p(p_l1:p_h1,ncomp_p)
  double precision u(u_l1:u_h1,ncomp_u)
  double precision dx(1), xlo(1), time, dt
  integer bc(1,2,ncomp_u), level, grid_no
  
  double precision :: e, T, X(nspec+naux)
  double precision :: rhoInv
  integer          :: i,n

  type (eos_t) :: eos_state

  !     Compute pressure from the EOS
  do i = lo(1),hi(1)
     
     rhoInv = 1.d0/u(i,URHO)
     e  = u(i,UEINT)*rhoInv
     T  = u(i,UTEMP)
     do n=1,nspec
        X(n) = u(i,UFS+n-1)*rhoInv
     enddo
     do n=1,naux
        X(nspec+n) = u(i,UFX+n-1)*rhoInv
     enddo
     
     ! Protect against negative internal energy
     if (allow_negative_energy .eq. 0 .and. e .le. 0.d0) then
        eos_state%rho = u(i,URHO)
        eos_state%T = T
        eos_state%xn(:) = X

        call eos(eos_input_rt, eos_state)

        p(i,1) = eos_state%p
     else
        eos_state%rho = u(i,URHO)
        eos_state%T = T
        eos_state%xn(:) = X
        eos_state%e = e

        call eos(eos_input_re, eos_state)

        p(i,1) = eos_state%p

     end if
     
     p(i,1) = p(i,1) - hse_p(i)
     
  enddo
  
end subroutine ca_derpi

!-----------------------------------------------------------------------

subroutine ca_derpioverp0(p,p_l1,p_h1,ncomp_p, &
                          u,u_l1,u_h1,ncomp_u,lo,hi,domlo, &
                          domhi,dx,xlo,time,dt,bc,level,grid_no)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : URHO, UMX, UEINT, UTEMP, UFS, UFX, &
                                 allow_negative_energy
  use probdata_module, only: hse_p

  implicit none

  integer p_l1,p_h1,ncomp_p
  integer u_l1,u_h1,ncomp_u
  integer lo(1), hi(1), domlo(1), domhi(1)
  double precision p(p_l1:p_h1,ncomp_p)
  double precision u(u_l1:u_h1,ncomp_u)
  double precision dx(1), xlo(1), time, dt
  integer bc(1,2,ncomp_u), level, grid_no
  
  double precision :: e, T, X(nspec+naux)
  double precision :: rhoInv
  integer          :: i,n

  type (eos_t) :: eos_state
  
  !     Compute pressure from the EOS
  do i = lo(1),hi(1)

     rhoInv = 1.d0/u(i,URHO)
     e  = u(i,UEINT)*rhoInv
     T  = u(i,UTEMP)
     do n=1,nspec
        X(n) = u(i,UFS+n-1)*rhoInv
     enddo
     do n=1,naux
        X(nspec+n) = u(i,UFX+n-1)*rhoInv
     enddo
     
     ! Protect against negative internal energy
     if (allow_negative_energy .eq. 0 .and. e .le. 0.d0) then
        eos_state%rho = u(i,URHO)
        eos_state%T = T
        eos_state%xn(:) = X

        call eos(eos_input_rt, eos_state)

        p(i,1) = eos_state%p
     else
        eos_state%rho = u(i,URHO)
        eos_state%T = T
        eos_state%xn(:) = X
        eos_state%e = e

        call eos(eos_input_re, eos_state)

        p(i,1) = eos_state%p

     end if
     
     p(i,1) = abs(p(i,1) - hse_p(i)) / hse_p(i)
     
  enddo
  
end subroutine ca_derpioverp0

