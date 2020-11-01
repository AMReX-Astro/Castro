module reduction_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

#if !(defined(AMREX_USE_CUDA) && defined(AMREX_GPU_PRAGMA_NO_HOST))

  subroutine reduce_max(x, y)

    implicit none

    ! Set in x the maximum of x and y.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    x = max(x, y)

  end subroutine reduce_max

#endif

#ifdef AMREX_USE_CUDA
  ! Select function name based on GPU pragma choice.

#ifdef AMREX_GPU_PRAGMA_NO_HOST
  attributes(device) subroutine reduce_max(x, y)
#else
  attributes(device) subroutine reduce_max_device(x, y)
#endif

    implicit none

    ! Set in x the maximum of x and y atomically on the GPU.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    real(rt) :: t

    t = atomicMax(x, y)

#ifdef AMREX_GPU_PRAGMA_NO_HOST
  end subroutine reduce_max
#else
  end subroutine reduce_max_device
#endif
#endif



#if !(defined(AMREX_USE_CUDA) && defined(AMREX_GPU_PRAGMA_NO_HOST))

  subroutine reduce_min(x, y)

    implicit none

    ! Set in x the minimum of x and y.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

#ifdef _OPENMP
    !$omp atomic update
#endif
    x = min(x, y)

  end subroutine reduce_min

#endif

#ifdef AMREX_USE_CUDA

#ifdef AMREX_GPU_PRAGMA_NO_HOST
  attributes(device) subroutine reduce_min(x, y)
#else
  attributes(device) subroutine reduce_min_device(x, y)
#endif

    implicit none

    ! Set in x the minimum of x and y atomically on the GPU.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    real(rt) :: t

    t = atomicMin(x, y)

#ifdef AMREX_GPU_PRAGMA_NO_HOST
  end subroutine reduce_min
#else
  end subroutine reduce_min_device
#endif
#endif



#if !(defined(AMREX_USE_CUDA) && defined(AMREX_GPU_PRAGMA_NO_HOST))

  subroutine reduce_add(x, y, blockReduce)

    implicit none

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x
    logical,  intent(in   ), optional :: blockReduce ! Only used in the CUDA version

    x = x + y

  end subroutine reduce_add

#endif

#ifdef AMREX_USE_CUDA

#ifdef AMREX_GPU_PRAGMA_NO_HOST
  attributes(device) subroutine reduce_add(x, y)
#else
  attributes(device) subroutine reduce_add_device(x, y, blockReduce)
#endif
    ! Do a shared memory reduction within a threadblock,
    ! then do an atomic add with a single thread in the block.

    implicit none

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x
    logical,  intent(in   ), optional :: blockReduce

    real(rt) :: t

    t = atomicAdd(x, y)

#ifdef AMREX_GPU_PRAGMA_NO_HOST
  end subroutine reduce_add
#else
  end subroutine reduce_add_device
#endif

#endif

end module reduction_module
