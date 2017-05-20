module c_interface_modules

  use meth_params_module, only: NVAR
  use amrex_fort_module, only: rt => amrex_real  

#ifdef CUDA
    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, &
                       cudaDeviceSynchronize, dim3, cuda_stream_kind
    use cuda_module, only: threads_and_blocks, cuda_streams, max_cuda_streams
#endif

contains

  subroutine ca_enforce_consistent_e(lo,hi,state,s_lo,s_hi,idx) &
       bind(c, name='ca_enforce_consistent_e')

    use castro_util_module, only: enforce_consistent_e

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer, intent(in)     :: idx

#ifdef CUDA

    attributes(device) :: state

    integer, device :: lo_d(3), hi_d(3)
    integer, device :: s_lo_d(3), s_hi_d(3)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice, stream)

    cuda_result = cudaMemcpyAsync(s_lo_d, s_lo, 3, cudaMemcpyHostToDevice, stream)
    cuda_result = cudaMemcpyAsync(s_hi_d, s_hi, 3, cudaMemcpyHostToDevice, stream)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_enforce_consistent_e<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, state, s_lo_d, s_hi_d)

#else

    call enforce_consistent_e(lo, hi, state, s_lo, s_hi)

#endif

  end subroutine ca_enforce_consistent_e


  subroutine ca_compute_temp(lo, hi, state, s_lo, s_hi, idx) &
       bind(C, name="ca_compute_temp")

    use castro_util_module, only: compute_temp

    implicit none

    integer, intent(in   ) :: lo(3),hi(3)
    integer, intent(in   ) :: s_lo(3),s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer, intent(in)     :: idx

    call compute_temp(lo, hi, state, s_lo, s_hi)

  end subroutine ca_compute_temp


  subroutine ca_reset_internal_e(lo, hi, u, u_lo, u_hi, verbose, idx) &
       bind(C, name="ca_reset_internal_e")

    use castro_util_module, only: reset_internal_e

    implicit none

    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    integer, intent(in)     :: idx

    call reset_internal_e(lo, hi, u, u_lo, u_hi, verbose)

  end subroutine ca_reset_internal_e

end module c_interface_modules
