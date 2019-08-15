subroutine derpi(p, p_lo, p_hi, ncomp_p, &
                 u, u_lo, u_hi, ncomp_u, &
                 lo, hi, domlo, domhi, &
                 dx, time) &
                 bind(C, name="derpi")

  use network, only : nspec, naux
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX

  use probdata_module
  use model_parser_module
  use prob_params_module, only: problo
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ), value :: ncomp_p
  integer,  intent(in   ), value :: ncomp_u
  real(rt), intent(in   ), value :: time

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)

  ! local
  integer :: i,j,k

  real(rt) :: y, z, pres, pres_local, height

  type (eos_t) :: eos_state

  do k = lo(3), hi(3)
     z = problo(3) + dx(3)*(dble(k) + 0.5e0_rt)

     do j = lo(2), hi(2)
        y = problo(2) + dx(2)*(dble(j) + 0.5e0_rt)

#if AMREX_SPACEDIM == 2
        height = y
#else
        height = z
#endif

        call interpolate_sub(pres, height, ipres_model)

        do i = lo(1), hi(1)

           eos_state%rho = u(i,j,k,URHO)
           eos_state%T = u(i,j,k,UTEMP)
           eos_state%xn(:) = u(i,j,k,UFS:UFS-1+nspec)/u(i,j,k,URHO)

           call eos(eos_input_rt, eos_state)

           pres_local = eos_state%p

           p(i,j,k,1) = pres_local - pres

        end do
     end do
  end do

end subroutine derpi

!-----------------------------------------------------------------------

subroutine derpioverp0(p, p_lo, p_hi, ncomp_p, &
                       u, u_lo, u_hi, ncomp_u, &
                       lo, hi, domlo, domhi, &
                       dx, time) &
                       bind(C, name="derpioverp0")

  use network, only : nspec, naux
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX
  use probdata_module
  use interpolate_module
  use model_parser_module
  use prob_params_module, only: problo
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ), value :: ncomp_p
  integer,  intent(in   ), value :: ncomp_u
  real(rt), intent(in   ), value :: time

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)


  ! local
  integer :: i, j, k

  real(rt) :: y, z, pres, pres_local, height

  type (eos_t) :: eos_state

  do k = lo(3), hi(3)
     z = problo(3) + dx(3)*(dble(k) + 0.5e0_rt)

     do j = lo(2), hi(2)
        y = problo(2) + dx(2)*(dble(j) + 0.5e0_rt)

#if AMREX_SPACEDIM == 2
        height = y
#else
        height = z
#endif

        call interpolate_sub(pres, height, ipres_model)

        do i = lo(1), hi(1)

           eos_state%rho = u(i,j,k,URHO)
           eos_state%T = u(i,j,k,UTEMP)
           eos_state%xn(:) = u(i,j,k,UFS:UFS-1+nspec)/u(i,j,k,URHO)

           call eos(eos_input_rt, eos_state)

           pres_local = eos_state%p

           p(i,j,k,1) = (pres_local - pres) / pres

        end do
     end do
  end do

end subroutine derpioverp0

!-----------------------------------------------------------------------

subroutine derrhopert(p, p_lo, p_hi, ncomp_p, &
                      u, u_lo, u_hi, ncomp_u, &
                      lo, hi, domlo, domhi, &
                      dx, time) &
                      bind(C, name="derrhopert")

  use network, only : nspec, naux
  use meth_params_module, only : URHO, UEINT, UTEMP
  use probdata_module
  use interpolate_module
  use model_parser_module
  use prob_params_module, only: problo
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ), value :: ncomp_p
  integer,  intent(in   ), value :: ncomp_u
  real(rt), intent(in   ), value :: time

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)

  ! local
  integer :: i, j, k

  real(rt) :: y, z, dens, height

  do k = lo(3), hi(3)
     z = problo(3) + dx(3)*(dble(k) + 0.5e0_rt)

     do j = lo(2), hi(2)
        y = problo(2) + dx(2)*(dble(j) + 0.5e0_rt)

#if AMREX_SPACEDIM == 2
        height = y
#else
        height = z
#endif

        call interpolate_sub(dens, height, idens_model)

        do i = lo(1), hi(1)
           p(i,j,k,1) = u(i,j,k,URHO) - dens
        end do

     end do
  end do

end subroutine derrhopert

!-----------------------------------------------------------------------

subroutine dertpert(p, p_lo, p_hi, ncomp_p, &
                    u, u_lo, u_hi, ncomp_u, &
                    lo, hi, domlo, domhi, &
                    dx, time) &
                    bind(C, name="dertpert")

  use network, only : nspec, naux
  use meth_params_module, only : URHO, UEINT, UTEMP
  use probdata_module
  use interpolate_module
  use model_parser_module
  use prob_params_module, only: problo
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ), value :: ncomp_p
  integer,  intent(in   ), value :: ncomp_u
  real(rt), intent(in   ), value :: time

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)

  ! local
  integer :: i, j, k

  real(rt) :: y, z, temp, height

  do k = lo(3), hi(3)
     z = problo(3) + dx(3)*(dble(k) + 0.5e0_rt)

     do j = lo(2), hi(2)
        y = problo(2) + dx(2)*(dble(j) + 0.5e0_rt)

#if AMREX_SPACEDIM == 2
        height = y
#else
        height = z
#endif

        call interpolate_sub(temp, height, itemp_model)

        do i = lo(1), hi(1)
           p(i,j,k,1) = u(i,j,k,UTEMP) - temp
        end do

     end do
  end do

end subroutine dertpert
