module generic_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use prob_params_module, only: dim
  use meth_params_module, only: NVAR

  implicit none

contains

  ! Used for a generic fill of any StateData.
  subroutine generic_single_fill(s_lo, s_hi, state, domlo, domhi, delta, xlo, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(dim), xlo(dim)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    call amrex_filccn(s_lo, s_hi, state, s_lo, s_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine generic_single_fill


  subroutine ca_generic_single_fill(state, s_lo, s_hi, &
                                    domlo, domhi, delta, xlo, time, bc) &
                                    bind(C, name="ca_generic_single_fill")

    implicit none

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(dim), xlo(dim), time
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    call generic_single_fill(s_lo, s_hi, state, domlo, domhi, delta, xlo, bc)

  end subroutine ca_generic_single_fill


  subroutine generic_multi_fill(s_lo, s_hi, state, domlo, domhi, delta, xlo, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: bc(dim,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(dim), xlo(dim)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    call amrex_filccn(s_lo, s_hi, state, s_lo, s_hi, NVAR, domlo, domhi, delta, xlo, bc)

  end subroutine generic_multi_fill


  subroutine ca_generic_multi_fill(state, s_lo, s_hi, &
                                   domlo, domhi, delta, xlo, time, bc) &
                                   bind(C, name="ca_generic_multi_fill")

    implicit none

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: bc(dim,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(dim), xlo(dim), time
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    call generic_multi_fill(s_lo, s_hi, state, domlo, domhi, delta, xlo, bc)

  end subroutine ca_generic_multi_fill

end module generic_fill_module
