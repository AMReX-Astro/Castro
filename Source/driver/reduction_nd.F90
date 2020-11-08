module reduction_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine reduce_max(x, y)

    implicit none

    ! Set in x the maximum of x and y.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    x = max(x, y)

  end subroutine reduce_max

#if (defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA))

  attributes(device) subroutine reduce_max_device(x, y)

    implicit none

    ! Set in x the maximum of x and y atomically on the GPU.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    real(rt) :: t

    t = atomicMax(x, y)

  end subroutine reduce_max_device

#endif



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

#if (defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA))

  attributes(device) subroutine reduce_min_device(x, y)

    implicit none

    ! Set in x the minimum of x and y atomically on the GPU.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    real(rt) :: t

    t = atomicMin(x, y)

  end subroutine reduce_min_device

#endif



  subroutine reduce_add(x, y, blockReduce)

    implicit none

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x
    logical,  intent(in   ), optional :: blockReduce ! Only used in the CUDA version

    x = x + y

  end subroutine reduce_add

#if (defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA))

  attributes(device) subroutine reduce_add_device(x, y, blockReduce)

    ! Do a shared memory reduction within a threadblock,
    ! then do an atomic add with a single thread in the block.

    implicit none

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x
    logical,  intent(in   ), optional :: blockReduce

    real(rt) :: t

    t = atomicAdd(x, y)

  end subroutine reduce_add_device

#endif

end module reduction_module
