  attributes(global) &
  subroutine cuda_enforce_consistent_e(lo,hi,state,s_lo,s_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR => NVAR_d
    use castro_util_module, only: enforce_consistent_e

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    attributes(device) :: lo, hi, s_lo, s_hi, state

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call enforce_consistent_e(idx, idx, state, s_lo, s_hi)

  end subroutine cuda_enforce_consistent_e




