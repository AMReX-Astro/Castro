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

#ifdef AMREX_USE_CUDA
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

    x = min(x, y)

  end subroutine reduce_min

#ifdef AMREX_USE_CUDA
  attributes(device) subroutine reduce_min_device(x, y)

    implicit none

    ! Set in x the minimum of x and y atomically on the GPU.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    real(rt) :: t

    t = atomicMin(x, y)

  end subroutine reduce_min_device
#endif



  subroutine reduce_add(x, y)

    implicit none

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    x = x + y

  end subroutine reduce_add

#ifdef AMREX_USE_CUDA
  attributes(device) function warpReduceSum(x) result(y)
    ! Reduce within a warp.
    ! https://devblogs.nvidia.com/faster-parallel-reductions-kepler/

    implicit none

    real(rt), intent(in) :: x

    real(rt) :: y

    integer :: offset

    offset = warpsize / 2

    y = x

    do while (offset > 0)

       y = y + __shfl_down(y, offset)

       offset = offset / 2

    end do

  end function warpReduceSum

  attributes(device) function blockReduceSum(x) result(y)
    ! Reduce within a threadblock.
    ! https://devblogs.nvidia.com/faster-parallel-reductions-kepler/

    implicit none

    real(rt), intent(in) :: x

    real(rt) :: y

    real(rt), shared :: s(0:(AMREX_GPU_MAX_THREADS/warpsize) - 1)

    integer :: lane, wid

    lane = mod(threadIdx%x-1, warpsize)
    wid = (threadIdx%x-1) / warpsize

    y = x
    y = warpReduceSum(y)

    ! syncthreads() prior to writing to shared memory is necessary
    ! if this reduction call is occurring multiple times in a kernel,
    ! and since we don't know how many times the user is calling it,
    ! we do it always to be safe.

    call syncthreads()

    if (lane == 0) then
       s(wid) = y
    end if

    call syncthreads()

    if ((threadIdx%x-1) < blockDim%x / warpsize) then
       y = s(lane)
    else
       y = 0
    end if

    if (wid == 0) then
       y = warpReduceSum(y)
    end if

  end function blockReduceSum

  attributes(device) subroutine reduce_add_device(x, y)
    ! Do a shared memory reduction within a threadblock,
    ! then do an atomic add with a single thread in the block.

    implicit none

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    real(rt) :: t

    t = y

    t = blockReduceSum(t)

    if (threadIdx%x == 1) then

       t = atomicAdd(x, t)

    end if

  end subroutine reduce_add_device
#endif

end module reduction_module
