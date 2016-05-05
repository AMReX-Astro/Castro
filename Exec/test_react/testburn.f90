subroutine do_burn() bind(C)

  use network
  use eos_module
  use burner_module
  use actual_burner_module
  use meth_params_module
  use reactions_module, only: ca_react_state

  implicit none

  integer, parameter :: nv = 7 + nspec

  double precision, parameter :: time = 0.0d0, dt = 1.0d-3

  integer, parameter :: lo(3) = [0, 0, 0], hi(3) = [63, 63, 63], w(3) = hi - lo

  double precision, parameter :: dens_min = 1.0d6, dens_max = 1.0d9
  double precision, parameter :: temp_min = 1.0d6, temp_max = 1.0d12

  double precision :: dlogrho, dlogT

  double precision :: state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nv)
  double precision :: reactions(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nspec+2)

  integer :: i, j, k

  type (eos_t) :: eos_state

  double precision :: start, finish

  character (len=32) :: probin_file
  integer :: probin_pass(32)
  integer :: n
  double precision :: start, finish

  probin_file = "probin"
  do n = 1, len(trim(probin_file))
     probin_pass(n) = ichar(probin_file(n:n))
  enddo

  call runtime_init(probin_pass(1:len(trim(probin_file))), len(trim(probin_file)))

  call network_init()
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

  ! Need to play the same hack for some of the other meth_params variables.

  react_T_min = temp_min / 2.0
  react_T_max = temp_max * 2.0

  react_rho_min = dens_min / 2.0
  react_rho_max = dens_max * 2.0

  disable_shock_burning = 0

  smallt = react_T_min
  smalld = react_rho_min

  !$acc update device(URHO, UTEMP, UEINT, UEDEN, UMX, UMY, UMZ, UFS, UFX)
  !$acc update device(react_T_min, react_T_max, react_rho_min, react_rho_max, disable_shock_burning)

  dlogrho = (log10(dens_max) - log10(dens_min)) / w(1)
  dlogT   = (log10(temp_max) - log10(temp_min)) / w(2)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           state(i,j,k,:) = ZERO

           eos_state % rho = 10.0d0**(log10(dens_min) + dble(i)*dlogrho)
           eos_state % T   = 10.0d0**(log10(temp_min) + dble(j)*dlogT  )
           eos_state % xn  = ZERO
           eos_state % xn(1 + INT( (dble(k) / (w(3) + 1)) * nspec)) = ONE

           call eos(eos_input_rt, eos_state)

           state(i,j,k,URHO)            = eos_state % rho
           state(i,j,k,UTEMP)           = eos_state % T
           state(i,j,k,UFS:UFS+nspec-1) = eos_state % rho * eos_state % xn
           state(i,j,k,UEINT)           = eos_state % rho * eos_state % e
           state(i,j,k,UEDEN)           = eos_state % rho * eos_state % e

        enddo
     enddo
  enddo

  call cpu_time(start)

  call ca_react_state(lo, hi, state, lo, hi, reactions, lo, hi, time, dt)

  call cpu_time(finish)

  print *, 'done!'
  print *, 'execution time = ', finish - start
  print *, 'sum of reactions energy = ', sum(reactions(:,:,:,nspec+1))

end subroutine do_burn
