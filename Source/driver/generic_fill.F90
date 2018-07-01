module generic_fill_module

  ! this module implements the boundary conditions for the Sources
  ! (either hydro sources or the SDC reaction sources)

  use amrex_fort_module, only: rt => amrex_real
  use prob_params_module, only: dim
#if (defined(AMREX_USE_CUDA) && !defined(AMREX_NO_DEVICE_LAUNCH))
  use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

  implicit none

contains

  ! Used for a generic fill of any StateData.
  AMREX_LAUNCH subroutine generic_single_fill(s_lo, s_hi, state, domlo, domhi, delta, xlo, bc)
    use prob_params_module, only: dim
    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(dim), xlo(dim)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, s_lo, s_hi)

    call amrex_filccn(blo, bhi, state, s_lo, s_hi, 1, domlo, domhi, delta, xlo, bc)

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

#if (defined(AMREX_USE_CUDA) && !defined(AMREX_NO_DEVICE_LAUNCH))
    attributes(device) :: state, s_lo, s_hi, bc, domlo, domhi, delta, xlo
#endif

    call generic_single_fill &
#if (defined(AMREX_USE_CUDA) && !defined(AMREX_NO_DEVICE_LAUNCH))
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (s_lo, s_hi, state, domlo, domhi, delta, xlo, bc)

  end subroutine ca_generic_single_fill


  AMREX_LAUNCH subroutine generic_multi_fill(s_lo, s_hi, state, domlo, domhi, delta, xlo, bc)
    use meth_params_module, only: NSRC
    use prob_params_module, only: dim
    use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: bc(dim,2,NSRC)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(dim), xlo(dim)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NSRC)

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, s_lo, s_hi)

    call amrex_filccn(blo, bhi, state, s_lo, s_hi, NSRC, domlo, domhi, delta, xlo, bc)

  end subroutine generic_multi_fill


  subroutine ca_generic_multi_fill(state, s_lo, s_hi, &
                                   domlo, domhi, delta, xlo, time, bc) &
                                   bind(C, name="ca_generic_multi_fill")

    use meth_params_module, only: NSRC

    implicit none

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: bc(dim,2,NSRC)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(dim), xlo(dim), time
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NSRC)

#if (defined(AMREX_USE_CUDA) && !defined(AMREX_NO_DEVICE_LAUNCH))
    attributes(device) :: state, s_lo, s_hi, bc, domlo, domhi, delta, xlo
#endif

    call generic_multi_fill &
#if (defined(AMREX_USE_CUDA) && !defined(AMREX_NO_DEVICE_LAUNCH))
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (s_lo, s_hi, state, domlo, domhi, delta, xlo, bc)

  end subroutine ca_generic_multi_fill

end module generic_fill_module
