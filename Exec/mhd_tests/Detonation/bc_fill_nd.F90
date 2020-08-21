module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

contains



  subroutine hypfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="hypfill")

    use amrex_constants_module, only: HALF
    use meth_params_module, only: NVAR, fill_ambient_bc
    use prob_params_module, only : problo

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2, NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer  :: i, j, k, n
    real(rt) :: x

    !$gpu

  end subroutine hypfill

  subroutine denfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="denfill")

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

    !$gpu

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
