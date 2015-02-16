subroutine ca_derpi(p,p_l1,p_l2,p_h1,p_h2,ncomp_p, &
                    u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
                    domhi,dx,xlo,time,dt,bc,level,grid_no)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                 allow_negative_energy
  use probdata_module
  use interpolate_module
  use model_parser_module
  
  implicit none

  integer          :: p_l1,p_l2,p_h1,p_h2,ncomp_p
  integer          :: u_l1,u_l2,u_h1,u_h2,ncomp_u
  integer          :: lo(2), hi(2), domlo(2), domhi(2)
  double precision :: p(p_l1:p_h1,p_l2:p_h2,ncomp_p)
  double precision :: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
  double precision :: dx(2), xlo(2), time, dt
  integer          :: bc(2,2,ncomp_u), level, grid_no

  ! local
  integer :: i,j

  double precision :: y,pres,pres_local

  type (eos_t) :: eos_state

  do j=lo(2),hi(2)
     y = xlo(2) + dx(2)*(float(j-lo(2)) + 0.5d0)
     pres = interpolate(y,npts_model,model_r, &
                        model_state(:,ipres_model))

     do i=lo(1),hi(1)

        eos_state%rho = u(i,j,URHO)
        eos_state%T = u(i,j,UTEMP)
        eos_state%xn(:) = u(i,j,UFS:)/u(i,j,URHO)

        call eos(eos_input_rt, eos_state)

        u(i,j,UEINT) = eos_state%e
        pres_local = eos_state%p

        p(i,j,1) = pres_local - pres

     end do
  end do

end subroutine ca_derpi

!-----------------------------------------------------------------------

subroutine ca_derpioverp0(p,p_l1,p_l2,p_h1,p_h2,ncomp_p, &
                          u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
                          domhi,dx,xlo,time,dt,bc,level,grid_no)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                 allow_negative_energy
  use probdata_module
  use interpolate_module
  use model_parser_module
  
  implicit none

  integer          :: p_l1,p_l2,p_h1,p_h2,ncomp_p
  integer          :: u_l1,u_l2,u_h1,u_h2,ncomp_u
  integer          :: lo(2), hi(2), domlo(2), domhi(2)
  double precision :: p(p_l1:p_h1,p_l2:p_h2,ncomp_p)
  double precision :: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
  double precision :: dx(2), xlo(2), time, dt
  integer          :: bc(2,2,ncomp_u), level, grid_no

  ! local
  integer :: i,j

  double precision :: y,pres,pres_local

  type (eos_t) :: eos_state

  do j=lo(2),hi(2)
     y = xlo(2) + dx(2)*(float(j-lo(2)) + 0.5d0)
     pres = interpolate(y,npts_model,model_r, &
                        model_state(:,ipres_model))

     do i=lo(1),hi(1)

        eos_state%rho = u(i,j,URHO)
        eos_state%T = u(i,j,UTEMP)
        eos_state%xn(:) = u(i,j,UFS:)/u(i,j,URHO)

        call eos(eos_input_rt, eos_state)

        u(i,j,UEINT) = eos_state%e
        pres_local = eos_state%p

        p(i,j,1) = (pres_local - pres) / pres

     end do
  end do

end subroutine ca_derpioverp0

!-----------------------------------------------------------------------

subroutine ca_derrhopert(p,p_l1,p_l2,p_h1,p_h2,ncomp_p, &
                         u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
                         domhi,dx,xlo,time,dt,bc,level,grid_no)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                 allow_negative_energy
  use probdata_module
  use interpolate_module
  use model_parser_module
  
  implicit none

  integer          :: p_l1,p_l2,p_h1,p_h2,ncomp_p
  integer          :: u_l1,u_l2,u_h1,u_h2,ncomp_u
  integer          :: lo(2), hi(2), domlo(2), domhi(2)
  double precision :: p(p_l1:p_h1,p_l2:p_h2,ncomp_p)
  double precision :: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
  double precision :: dx(2), xlo(2), time, dt
  integer          :: bc(2,2,ncomp_u), level, grid_no

  ! local
  integer :: i,j

  double precision :: y,dens

  do j=lo(2),hi(2)
     y = xlo(2) + dx(2)*(float(j-lo(2)) + 0.5d0)

     dens = interpolate(y,npts_model,model_r, &
                        model_state(:,idens_model))

     do i=lo(1),hi(1)
        p(i,j,1) = u(i,j,URHO) - dens
     end do

  end do

end subroutine ca_derrhopert

!-----------------------------------------------------------------------

subroutine ca_dertpert(p,p_l1,p_l2,p_h1,p_h2,ncomp_p, &
                       u,u_l1,u_l2,u_h1,u_h2,ncomp_u,lo,hi,domlo, &
                       domhi,dx,xlo,time,dt,bc,level,grid_no)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                 allow_negative_energy
  use probdata_module
  use interpolate_module
  use model_parser_module
  
  implicit none

  integer          :: p_l1,p_l2,p_h1,p_h2,ncomp_p
  integer          :: u_l1,u_l2,u_h1,u_h2,ncomp_u
  integer          :: lo(2), hi(2), domlo(2), domhi(2)
  double precision :: p(p_l1:p_h1,p_l2:p_h2,ncomp_p)
  double precision :: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
  double precision :: dx(2), xlo(2), time, dt
  integer          :: bc(2,2,ncomp_u), level, grid_no

  ! local
  integer :: i,j

  double precision :: y,temp

  do j=lo(2),hi(2)
     y = xlo(2) + dx(2)*(float(j-lo(2)) + 0.5d0)

     temp = interpolate(y,npts_model,model_r, &
                        model_state(:,itemp_model))

     do i=lo(1),hi(1)
        p(i,j,1) = u(i,j,UTEMP) - temp
     end do

  end do

end subroutine ca_dertpert
