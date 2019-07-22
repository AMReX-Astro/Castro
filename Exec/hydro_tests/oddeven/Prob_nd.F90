subroutine amrex_probinit (init, name, namlen, problo, probhi) bind(c)

  use prob_params_module, only: center
  use probdata_module, only: p_ambient, dens_ambient, dens_pert_factor, vel_pert
  use castro_error_module, only: castro_error
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in   ) :: init, namlen
  integer,  intent(in   ) :: name(namlen)
  real(rt), intent(in   ) :: problo(3), probhi(3)

  integer :: untin, i

  namelist /fortin/ p_ambient, dens_ambient, dens_pert_factor, vel_pert

  ! Build "probin" filename -- the name of file containing fortin namelist.

  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) call castro_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set center, domain extrema
  center(:) = (problo(:) + probhi(:)) / 2.e0_rt

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin,fortin)
  close(unit=untin)

end subroutine amrex_probinit

subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       dx, xlo, xhi)

  use amrex_constants_module, only: ZERO, HALF, ONE
  use probdata_module, only: p_ambient, dens_ambient, dens_pert_factor, vel_pert
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_rp
  use network, only: nspec
  use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module, only: center
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, dx(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: xcen, ycen, zcen
  real(rt) :: dens, eint, xvel, X(nspec), temp
  integer  :: i, j, k, icen, jcen, kcen

  type (eos_t) :: eos_state

  ! compute the integer location of the center of the domain
  icen = center(1) / dx(1)

#if AMREX_SPACEDIM >= 2
  jcen = center(2) / dx(2)
#else
  jcen = 0
#endif

#if AMREX_SPACEDIM == 3
  kcen = center(3) / dx(3)
#else
  kcen = 0
#endif

  do k = lo(3), hi(3)
     zcen = xlo(3) + dx(3)*(dble(k-lo(3)) + HALF)

     do j = lo(2), hi(2)
        ycen = xlo(2) + dx(2)*(dble(j-lo(2)) + HALF)

        do i = lo(1), hi(1)
           xcen = xlo(1) + dx(1)*(dble(i-lo(1)) + HALF)

           if (i == icen .and. j == jcen .and. k == kcen) then
              dens = dens_ambient*dens_pert_factor
           else
              dens = dens_ambient
           endif

           state(i,j,k,URHO) = dens

           ! velocity perturbation
           if (xcen < center(1)) then 
              xvel = vel_pert
           else if (xcen > center(1)) then
              xvel = -vel_pert
           else
              xvel = ZERO
           endif

           state(i,j,k,UMX) = dens*xvel
           state(i,j,k,UMY) = ZERO
           state(i,j,k,UMZ) = ZERO

           ! set the composition
           X(:) = ZERO
           X(1) = ONE

           ! compute the internal energy and temperature
           eos_state%T = ONE ! initial guess
           eos_state%rho = dens
           eos_state%p = p_ambient
           eos_state%xn(:) = X

           call eos(eos_input_rp, eos_state)

           temp = eos_state%T
           eint = eos_state%e

           state(i,j,k,UEDEN) = dens * eint +  &
                                HALF * (state(i,j,k,UMX)**2 / state(i,j,k,URHO) + &
                                        state(i,j,k,UMY)**2 / state(i,j,k,URHO) + &
                                        state(i,j,k,UMZ)**2 / state(i,j,k,URHO))

           state(i,j,k,UEINT) = dens * eint
           state(i,j,k,UTEMP) = temp

           state(i,j,k,UFS:UFS-1+nspec) = dens * X(:)

        end do
     end do
  end do

end subroutine ca_initdata

