subroutine ca_derpi(p,p_lo,p_hi,ncomp_p, &
                    u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                    domhi,dx,xlo,time,dt,bc,level,grid_no) &
                    bind(C, name="ca_derpi")

  use network, only : nspec, naux
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX, &
                                 allow_negative_energy
  use probdata_module
  use interpolate_module
  use model_parser_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: p_lo(3),p_hi(3),ncomp_p
  integer          :: u_lo(3),u_hi(3),ncomp_u
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt)         :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt)         :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  ! local
  integer :: i,j,k

  real(rt)         :: y,pres,pres_local

  type (eos_t) :: eos_state

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        y = xlo(2) + dx(2)*(float(j-lo(2)) + 0.5e0_rt)
        pres = interpolate(y,npts_model,model_r, &
                           model_state(:,ipres_model))

        do i=lo(1),hi(1)

           eos_state%rho = u(i,j,k,URHO)
           eos_state%T = u(i,j,k,UTEMP)
           eos_state%xn(:) = u(i,j,k,UFS:UFS-1+nspec)/u(i,j,k,URHO)

           call eos(eos_input_rt, eos_state)

           u(i,j,k,UEINT) = eos_state%e
           pres_local = eos_state%p

           p(i,j,k,1) = pres_local - pres

        end do
     end do
  end do

end subroutine ca_derpi

!-----------------------------------------------------------------------

subroutine ca_derpioverp0(p,p_lo,p_hi,ncomp_p, &
                          u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                          domhi,dx,xlo,time,dt,bc,level,grid_no) &
                          bind(C, name="ca_derpioverp0")

  use network, only : nspec, naux
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX, &
                                 allow_negative_energy
  use probdata_module
  use interpolate_module
  use model_parser_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: p_lo(3),p_hi(3),ncomp_p
  integer          :: u_lo(3),u_hi(3),ncomp_u
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt)         :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt)         :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  ! local
  integer :: i,j,k

  real(rt)         :: y,pres,pres_local

  type (eos_t) :: eos_state

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        y = xlo(2) + dx(2)*(float(j-lo(2)) + 0.5e0_rt)
        pres = interpolate(y,npts_model,model_r, &
                           model_state(:,ipres_model))

        do i=lo(1),hi(1)

           eos_state%rho = u(i,j,k,URHO)
           eos_state%T = u(i,j,k,UTEMP)
           eos_state%xn(:) = u(i,j,k,UFS:UFS-1+nspec)/u(i,j,k,URHO)

           call eos(eos_input_rt, eos_state)

           u(i,j,k,UEINT) = eos_state%e
           pres_local = eos_state%p

           p(i,j,k,1) = (pres_local - pres) / pres

        end do
     end do
  end do

end subroutine ca_derpioverp0

!-----------------------------------------------------------------------

subroutine ca_derrhopert(p,p_lo,p_hi,ncomp_p, &
                         u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                         domhi,dx,xlo,time,dt,bc,level,grid_no) &
                         bind(C, name="ca_derrhopert")

  use network, only : nspec, naux
  use meth_params_module, only : URHO, UEINT, UTEMP, &
                                 allow_negative_energy
  use probdata_module
  use interpolate_module
  use model_parser_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: p_lo(3),p_hi(3),ncomp_p
  integer          :: u_lo(3),u_hi(3),ncomp_u
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt)         :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt)         :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  ! local
  integer :: i,j,k

  real(rt)         :: y,dens

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        y = xlo(2) + dx(2)*(float(j-lo(2)) + 0.5e0_rt)

        dens = interpolate(y,npts_model,model_r, &
                           model_state(:,idens_model))

        do i=lo(1),hi(1)
           p(i,j,k,1) = u(i,j,k,URHO) - dens
        end do

     end do
  end do

end subroutine ca_derrhopert

!-----------------------------------------------------------------------

subroutine ca_dertpert(p,p_lo,p_hi,ncomp_p, &
                       u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                       domhi,dx,xlo,time,dt,bc,level,grid_no) &
                       bind(C, name="ca_dertpert")

  use network, only : nspec, naux
  use meth_params_module, only : URHO, UEINT, UTEMP, &
                                 allow_negative_energy
  use probdata_module
  use interpolate_module
  use model_parser_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: p_lo(3),p_hi(3),ncomp_p
  integer          :: u_lo(3),u_hi(3),ncomp_u
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt)         :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt)         :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt)         :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  ! local
  integer :: i,j,k

  real(rt)         :: y,temp

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        y = xlo(2) + dx(2)*(float(j-lo(2)) + 0.5e0_rt)

        temp = interpolate(y,npts_model,model_r, &
                           model_state(:,itemp_model))

        do i=lo(1),hi(1)
           p(i,j,k,1) = u(i,j,k,UTEMP) - temp
        end do

     end do
  end do

end subroutine ca_dertpert
