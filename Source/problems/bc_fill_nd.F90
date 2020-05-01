module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine hypfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="hypfill")

    use amrex_filcc_module, only: amrex_filccn
    use meth_params_module, only: fill_ambient_bc, ambient_fill_dir, ambient_outflow_vel, UMX, UMY, UMZ, UEDEN, URHO, UEINT
    use ambient_module, only: ambient_state
    use amrex_bc_types_module, only: amrex_bc_foextrap, amrex_bc_hoextrap
    use amrex_constants_module, only: ZERO, HALF

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2, NVAR)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer :: i, j, k
    logical :: ambient_x_lo, ambient_y_lo, ambient_z_lo
    logical :: ambient_x_hi, ambient_y_hi, ambient_z_hi

    !$gpu

    ambient_x_lo = (ambient_fill_dir == 0 .or. ambient_fill_dir == -1) .and. &
         (bc(1,1,1) == amrex_bc_foextrap .or. bc(1,1,1) == amrex_bc_hoextrap)
    ambient_x_hi = (ambient_fill_dir == 0 .or. ambient_fill_dir == -1) .and. &
         (bc(1,2,1) == amrex_bc_foextrap .or. bc(1,2,1) == amrex_bc_hoextrap)

#if AMREX_SPACEDIM >= 2
    ambient_y_lo = (ambient_fill_dir == 1 .or. ambient_fill_dir == -1) .and. &
         (bc(2,1,1) == amrex_bc_foextrap .or. bc(2,1,1) == amrex_bc_hoextrap)
    ambient_y_hi = (ambient_fill_dir == 1 .or. ambient_fill_dir == -1) .and. &
         (bc(2,2,1) == amrex_bc_foextrap .or. bc(2,2,1) == amrex_bc_hoextrap)
#endif

#if AMREX_SPACEDIM == 3
    ambient_z_lo = (ambient_fill_dir == 2 .or. ambient_fill_dir == -1) .and. &
         (bc(3,1,1) == amrex_bc_foextrap .or. bc(3,1,1) == amrex_bc_hoextrap)
    ambient_z_hi = (ambient_fill_dir == 2 .or. ambient_fill_dir == -1) .and. &
         (bc(3,2,1) == amrex_bc_foextrap .or. bc(3,2,1) == amrex_bc_hoextrap)
#endif

    if (fill_ambient_bc == 1) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (ambient_x_lo .and. i < domlo(1) .or. &
                    ambient_x_hi .and. i > domhi(1) &
#if AMREX_SPACEDIM >= 2
                    .or. &
                    ambient_y_lo .and. j < domlo(2) .or. &
                    ambient_y_hi .and. j > domhi(2) &
#endif
#if AMREX_SPACEDIM == 3
                    .or. &
                    ambient_z_lo .and. k < domlo(3) .or. &
                    ambient_z_hi .and. k > domhi(3) &
#endif
                    ) then
                   adv(i,j,k,:) = ambient_state(:)

                   if (ambient_outflow_vel == 1) then

                      ! extrapolate the normal velocity only if it is outgoing
                      if (i < domlo(1)) then
                         adv(i,j,k,UMX) = min(ZERO, adv(domlo(1),j,k,UMX))
                         adv(i,j,k,UMY) = ZERO
                         adv(i,j,k,UMZ) = ZERO

                      else if (i > domhi(1)) then
                         adv(i,j,k,UMX) = max(ZERO, adv(domhi(1),j,k,UMX))
                         adv(i,j,k,UMY) = ZERO
                         adv(i,j,k,UMZ) = ZERO

                      else if (j < domlo(2)) then
                         adv(i,j,k,UMX) = ZERO
                         adv(i,j,k,UMY) = min(ZERO, adv(i,domlo(2),k,UMY))
                         adv(i,j,k,UMZ) = ZERO

                      else if (j > domhi(2)) then
                         adv(i,j,k,UMX) = ZERO
                         adv(i,j,k,UMY) = max(ZERO, adv(i,domhi(2),k,UMY))
                         adv(i,j,k,UMZ) = ZERO

                      else if (k < domlo(3)) then
                         adv(i,j,k,UMX) = ZERO
                         adv(i,j,k,UMY) = ZERO
                         adv(i,j,k,UMZ) = min(ZERO, adv(i,j,domlo(3),UMZ))

                      else
                         adv(i,j,k,UMX) = ZERO
                         adv(i,j,k,UMY) = ZERO
                         adv(i,j,k,UMZ) = max(ZERO, adv(i,j,domhi(3),UMZ))

                      end if

                      ! now make the energy consistent
                      adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) + &
                           HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)

                   end if

                end if
             end do
          end do
       end do
    end if

  end subroutine hypfill



  subroutine denfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="denfill")

    use amrex_filcc_module, only: amrex_filccn
    use meth_params_module, only: fill_ambient_bc, URHO, ambient_fill_dir
    use ambient_module, only: ambient_state
    use amrex_bc_types_module, only: amrex_bc_foextrap, amrex_bc_hoextrap

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: i, j, k

    logical :: ambient_x_lo, ambient_y_lo, ambient_z_lo
    logical :: ambient_x_hi, ambient_y_hi, ambient_z_hi

    !$gpu

    ambient_x_lo = (ambient_fill_dir == 0 .or. ambient_fill_dir == -1) .and. &
         (bc(1,1) == amrex_bc_foextrap .or. bc(1,1) == amrex_bc_hoextrap)
    ambient_x_hi = (ambient_fill_dir == 0 .or. ambient_fill_dir == -1) .and. &
         (bc(1,2) == amrex_bc_foextrap .or. bc(1,2) == amrex_bc_hoextrap)

#if AMREX_SPACEDIM >= 2
    ambient_y_lo = (ambient_fill_dir == 1 .or. ambient_fill_dir == -1) .and. &
         (bc(2,1) == amrex_bc_foextrap .or. bc(2,1) == amrex_bc_hoextrap)
    ambient_y_hi = (ambient_fill_dir == 1 .or. ambient_fill_dir == -1) .and. &
         (bc(2,2) == amrex_bc_foextrap .or. bc(2,2) == amrex_bc_hoextrap)
#endif

#if AMREX_SPACEDIM == 3
    ambient_z_lo = (ambient_fill_dir == 2 .or. ambient_fill_dir == -1) .and. &
         (bc(3,1) == amrex_bc_foextrap .or. bc(3,1) == amrex_bc_hoextrap)
    ambient_z_hi = (ambient_fill_dir == 2 .or. ambient_fill_dir == -1) .and. &
         (bc(3,2) == amrex_bc_foextrap .or. bc(3,2) == amrex_bc_hoextrap)
#endif

    if (fill_ambient_bc == 1) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (ambient_x_lo .and. i < domlo(1) .or. &
                    ambient_x_hi .and. i > domhi(1) &
#if AMREX_SPACEDIM >= 2
                    .or. &
                    ambient_y_lo .and. j < domlo(2) .or. &
                    ambient_y_hi .and. j > domhi(2) &
#endif
#if AMREX_SPACEDIM == 3
                    .or. &
                    ambient_z_lo .and. k < domlo(3) .or. &
                    ambient_z_hi .and. k > domhi(3) &
#endif
                    ) then
                   adv(i,j,k) = ambient_state(URHO)
                end if
             end do
          end do
       end do
    end if

  end subroutine denfill

#ifdef MHD
  subroutine face_fillx(lo, hi, var, var_lo, var_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="face_fillx")
          
    use amrex_fort_module, only : rt => amrex_real
    use fc_fill_module

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: var_lo(3), var_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: var(var_lo(1):var_hi(1), var_lo(2):var_hi(2), var_lo(3):var_hi(3))
    real(rt), intent(in   ), value :: time
    integer dir

    dir = 1

    call filfc(var,var_lo(1),var_lo(2),var_lo(3),var_hi(1),var_hi(2),var_hi(3),domlo,domhi,delta,xlo,bc,dir)

  end subroutine face_fillx


  subroutine face_filly(lo, hi, var, var_lo, var_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="face_filly")
                       
    use amrex_fort_module, only : rt => amrex_real
    use fc_fill_module

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: var_lo(3), var_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: var(var_lo(1):var_hi(1), var_lo(2):var_hi(2), var_lo(3):var_hi(3))
    real(rt), intent(in   ), value :: time 
    integer dir

    dir = 2

    call filfc(var,var_lo(1),var_lo(2),var_lo(3),var_hi(1),var_hi(2),var_hi(3),domlo,domhi,delta,xlo,bc,dir)

  end subroutine face_filly

  subroutine face_fillz(lo, hi, var, var_lo, var_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="face_fillz")
                     
    use amrex_fort_module, only : rt => amrex_real
    use fc_fill_module

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: var_lo(3), var_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM,2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: var(var_lo(1):var_hi(1), var_lo(2):var_hi(2), var_lo(3):var_hi(3))
    real(rt), intent(in   ), value :: time
    integer dir

    dir = 3

    call filfc(var,var_lo(1),var_lo(2),var_lo(3),var_hi(1),var_hi(2),var_hi(3),domlo,domhi,delta,xlo,bc,dir)

  end subroutine face_fillz
#endif

end module bc_fill_module
