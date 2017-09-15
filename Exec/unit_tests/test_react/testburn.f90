subroutine do_burn() bind(C)

  use network
  use eos_type_module, only : eos_t, eos_input_rt
  use eos_module
  use burner_module
  use actual_burner_module
  use actual_rhs_module, only: actual_rhs_init
  use meth_params_module
  use reactions_module, only: ca_react_state
  use extern_probin_module, only: ncell, dt, dens_min, dens_max, temp_min, temp_max

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, parameter :: nv = 7 + nspec

  real(rt)        , parameter :: time = 0.0e0_rt

  integer :: lo(3), hi(3), w(3)

  real(rt)         :: dlogrho, dlogT

  real(rt)        , allocatable :: state(:,:,:,:), reactions(:,:,:,:)
  integer, allocatable :: mask(:,:,:)
  real(rt), allocatable :: weights(:,:,:)

  integer :: i, j, k

  type (eos_t) :: eos_state

  real(rt)         :: start, finish

  character (len=32) :: probin_file
  integer :: probin_pass(32)
  integer :: n

  probin_file = "probin"
  do n = 1, len(trim(probin_file))
     probin_pass(n) = ichar(probin_file(n:n))
  enddo

  call runtime_init(probin_pass(1:len(trim(probin_file))), len(trim(probin_file)))
  call ca_set_castro_method_params()

  call network_init()
  call actual_rhs_init()
  call burner_init()
  call eos_init()

  ! Normally the state indices get loaded in through the main C++ interface,
  ! so we have to hack around that by explicitly setting them here.

  NVAR = nv
  URHO = 1
  UTEMP = 2
  UEINT = 3
  UEDEN = 4
  UMX = 5
  UMY = 6
  UMZ = 7
  UFS = 8
  UFX = -1

  !$acc update device(URHO, UTEMP, UEINT, UEDEN, UMX, UMY, UMZ, UFS, UFX)

  lo = [0, 0, 0]
  hi = [ncell, ncell, ncell]
  w = hi - lo + 1

  allocate(state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),NVAR))
  allocate(reactions(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nspec+2))
  allocate(mask(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
  allocate(weights(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

  dlogrho = (log10(dens_max) - log10(dens_min)) / w(1)
  dlogT   = (log10(temp_max) - log10(temp_min)) / w(2)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           state(i,j,k,:) = ZERO

           eos_state % rho = 10.0e0_rt**(log10(dens_min) + dble(i)*dlogrho)
           eos_state % T   = 10.0e0_rt**(log10(temp_min) + dble(j)*dlogT  )
           eos_state % xn  = 1.e-12_rt
           eos_state % xn(1 + INT( (dble(k) / w(3)) * nspec)) = ONE  - (nspec - 1) * 1.e-12_rt

           call eos(eos_input_rt, eos_state)

           state(i,j,k,URHO)            = eos_state % rho
           state(i,j,k,UTEMP)           = eos_state % T
           state(i,j,k,UFS:UFS+nspec-1) = eos_state % rho * eos_state % xn
           state(i,j,k,UEINT)           = eos_state % rho * eos_state % e
           state(i,j,k,UEDEN)           = eos_state % rho * eos_state % e
           state(i,j,k,UMX:UMZ)         = ZERO

           mask(i,j,k) = 1
           weights(i,j,k) = 0.0

        enddo
     enddo
  enddo

  call cpu_time(start)

  call ca_react_state(lo, hi, state, lo, hi, reactions, lo, hi, &
                      weights, lo, hi, &
                      mask, lo, hi, time, dt, 0)

  call cpu_time(finish)

  print *, 'done!'
  print *, 'execution time = ', finish - start
  print *, 'sum of reactions energy = ', sum(reactions(:,:,:,nspec+1))

end subroutine do_burn
