module c_interface_modules

  use meth_params_module, only: NVAR, NQAUX, NQ, QVAR
  use amrex_fort_module, only: rt => amrex_real  

#ifdef CUDA
    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, &
                       cudaDeviceSynchronize, dim3, cuda_stream_kind
    use cuda_module, only: threads_and_blocks, cuda_streams, max_cuda_streams
#endif

contains

  subroutine ca_enforce_consistent_e(lo,hi, &
#ifdef MHD
                                     bx, bx_lo, bx_hi,&
                                     by, by_lo, by_hi,&
                                     bz, bz_lo, bz_hi,&
#endif
                                     state,s_lo,s_hi,idx) &
       bind(c, name='ca_enforce_consistent_e')

    use castro_util_module, only: enforce_consistent_e

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer, intent(in)     :: idx
#ifdef MHD
     integer, intent(in)     :: bx_lo(3), bx_hi(3)
    integer, intent(in)     :: by_lo(3), by_hi(3)
    integer, intent(in)     :: bz_lo(3), bz_hi(3)
    real(rt), intent(in)    :: bx(bx_lo(1):bx_hi(1), bx_lo(2):bx_hi(2), bx_lo(3):bx_hi(3))
    real(rt), intent(in)    :: by(by_lo(1):by_hi(1), by_lo(2):by_hi(2), by_lo(3):by_hi(3))
    real(rt), intent(in)    :: bz(bz_lo(1):bz_hi(1), bz_lo(2):bz_hi(2), bz_lo(3):bz_hi(3))
#endif

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

    call enforce_consistent_e(lo, hi,&
#ifdef MHD
                              bx, bx_lo, bx_hi, &
                              by, by_lo, by_hi, &
                              bz, bz_lo, bz_hi, &
#endif
                              state, s_lo, s_hi)

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


  subroutine ca_reset_internal_e(lo, hi, &
#ifdef MHD
                                 bx, bx_lo, bx_hi, &
                                 by, by_lo, by_hi, &
                                 bz, bz_lo, bz_hi, &
#endif
                                 u, u_lo, u_hi, verbose, idx) &
       bind(C, name="ca_reset_internal_e")

    use castro_util_module, only: reset_internal_e

    implicit none

    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
#ifdef MHD
    integer, intent(in) :: bx_lo(3), bx_hi(3)
    integer, intent(in) :: by_lo(3), by_hi(3)
    integer, intent(in) :: bz_lo(3), bz_hi(3)
    real(rt), intent(in):: bx(bx_lo(1):bx_hi(1),bx_lo(2):bx_hi(2),bx_lo(3):bx_hi(3))
    real(rt), intent(in):: by(by_lo(1):by_hi(1),by_lo(2):by_hi(2),by_lo(3):by_hi(3))
    real(rt), intent(in):: bz(bz_lo(1):bz_hi(1),bz_lo(2):bz_hi(2),bz_lo(3):bz_hi(3))
#endif
    integer, intent(in)     :: idx

    call reset_internal_e(lo, hi, &
#ifdef MHD
                          bx, bx_lo, bx_hi, &
                          by, by_lo, by_hi, &
                          bz, bz_lo, bz_hi, &
#endif
                          u, u_lo, u_hi, verbose)

  end subroutine ca_reset_internal_e


  subroutine ca_normalize_species(lo,hi,u, u_lo, u_hi, idx) &
       bind(C, name="ca_normalize_species")

    use castro_util_module, only: normalize_species

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    integer, intent(in)     :: idx

    call normalize_species(u, u_lo, u_hi, lo, hi)

  end subroutine ca_normalize_species


  subroutine ca_enforce_minimum_density(lo, hi, uin, uin_lo, uin_hi, &
                                        uout, uout_lo, uout_hi, &
                                        vol, vol_lo, vol_hi, &
                                        frac_change, verbose, idx) &
                                        bind(C, name="ca_enforce_minimum_density")
#ifndef MHD 
    use advection_util_module, only: enforce_minimum_density
#else
    use mhd_util_module, only: enforce_minimum_density
#endif
    implicit none

    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) ::  uin_lo(3),  uin_hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) ::  vol_lo(3),  vol_hi(3)

    real(rt), intent(in) ::  uin( uin_lo(1): uin_hi(1), uin_lo(2): uin_hi(2), uin_lo(3): uin_hi(3),NVAR)
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt), intent(in) ::  vol( vol_lo(1): vol_hi(1), vol_lo(2): vol_hi(2), vol_lo(3): vol_hi(3))
    real(rt), intent(inout) :: frac_change
    integer, intent(in)     :: idx

    call enforce_minimum_density(uin, uin_lo, uin_hi, &
                                 uout, uout_lo, uout_hi, &
                                 vol, vol_lo, vol_hi, &
                                 lo, hi, frac_change, verbose)                

  end subroutine ca_enforce_minimum_density


  subroutine ca_check_initial_species(lo, hi, state, state_lo, state_hi, idx) &
                                      bind(C, name="ca_check_initial_species")

    use castro_util_module, only: check_initial_species

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: state_lo(3), state_hi(3)
    real(rt), intent(in) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    integer, intent(in)     :: idx

    call check_initial_species(lo, hi, state, state_lo, state_hi)

  end subroutine ca_check_initial_species

#ifndef MHD
  subroutine ca_ctoprim(lo, hi, &
                        uin, uin_lo, uin_hi, &
#ifdef RADIATION
                        Erin, Erin_lo, Erin_hi, &
                        lam, lam_lo, lam_hi, &
#endif
                        q,     q_lo,   q_hi, &
                        qaux, qa_lo,  qa_hi, idx) bind(C, name = "ca_ctoprim")

    use advection_util_module, only: ctoprim

#ifdef RADIATION
    use rad_params_module, only : ngroups
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
#ifdef RADIATION
    integer, intent(in) :: Erin_lo(3), Erin_hi(3)
    integer, intent(in) :: lam_lo(3), lam_hi(3)
#endif
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    real(rt)        , intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
#ifdef RADIATION
    real(rt)        , intent(in   ) :: Erin(Erin_lo(1):Erin_hi(1),Erin_lo(2):Erin_hi(2),Erin_lo(3):Erin_hi(3),0:ngroups-1)
    real(rt)        , intent(in   ) :: lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
#endif

    real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    integer, intent(in)     :: idx

    call ctoprim(lo, hi, &
                 uin, uin_lo, uin_hi, &
#ifdef RADIATION
                 Erin, Erin_lo, Erin_hi, &
                 lam, lam_lo, lam_hi, &
#endif
                 q,     q_lo,   q_hi, &
                 qaux, qa_lo,  qa_hi)

  end subroutine ca_ctoprim


  subroutine ca_srctoprim(lo, hi, &
                          q,     q_lo,   q_hi, &
                          qaux, qa_lo,  qa_hi, &
                          src, src_lo, src_hi, &
                          srcQ,srQ_lo, srQ_hi, idx) bind(C, name = "ca_srctoprim")

    use advection_util_module, only: srctoprim

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3),   qa_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: srQ_lo(3), srQ_hi(3)

    real(rt)        , intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)        , intent(in   ) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt)        , intent(in   ) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
    real(rt)        , intent(inout) :: srcQ(srQ_lo(1):srQ_hi(1),srQ_lo(2):srQ_hi(2),srQ_lo(3):srQ_hi(3),QVAR)
    integer, intent(in)     :: idx

    call srctoprim(lo, hi, &
                   q,     q_lo,   q_hi, &
                   qaux, qa_lo,  qa_hi, &
                   src, src_lo, src_hi, &
                   srcQ,srQ_lo, srQ_hi)

  end subroutine ca_srctoprim
#endif

end module c_interface_modules
