module generic_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR

  implicit none

  integer :: nvar_bc(AMREX_SPACEDIM, 2, NVAR)
  integer :: bc(AMREX_SPACEDIM, 2)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: nvar_bc, bc
#endif

contains

  subroutine prepare_nvar_bc(bc_in) bind(C, name="prepare_nvar_bc")

    implicit none

    integer, intent(in) :: bc_in(AMREX_SPACEDIM, 2, NVAR)

    nvar_bc(:,:,:) = bc_in(:,:,:)

  end subroutine prepare_nvar_bc

  subroutine prepare_bc(bc_in) bind(C, name="prepare_bc")

    implicit none

    integer, intent(in) :: bc_in(AMREX_SPACEDIM, 2)

    bc(:,:) = bc_in(:,:)

  end subroutine prepare_bc



  ! Used for a generic fill of any StateData.
  subroutine generic_single_fill(lo, hi, state, s_lo, s_hi, domlo, domhi, delta, xlo) bind(C, name="generic_single_fill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    !$gpu

    call amrex_filccn(lo, hi, state, s_lo, s_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine generic_single_fill



  subroutine generic_multi_fill(lo, hi, state, s_lo, s_hi, domlo, domhi, delta, xlo) bind(C, name="generic_multi_fill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    !$gpu

    call amrex_filccn(lo, hi, state, s_lo, s_hi, NVAR, domlo, domhi, delta, xlo, nvar_bc)

  end subroutine generic_multi_fill

end module generic_fill_module
