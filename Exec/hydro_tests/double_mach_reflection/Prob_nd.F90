subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use eos_module
  use eos_type_module
  use castro_error_module
  use network
  use probdata_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer :: untin, i

  real(rt) :: xn(nspec)

  type (eos_t) :: eos_state

  namelist /fortin/ p_l, u_l, v_l, rho_l, rhoe_l, &
                    p_r, u_r, v_r, rho_r, rhoe_r, &
                    T_l, T_r, use_Tinit

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen=256
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call castro_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

#if AMREX_SPACEDIM == 1 || AMREX_SPACEDIM == 3
  call castro_error("ERROR: this problem only works for 2-d")
#endif

  ! set namelist defaults

  allocate(p_l)
  allocate(u_l)
  allocate(v_l)
  allocate(rho_l)
  allocate(T_l)

  p_l = 116.5             ! left pressure (erg/cc)
  u_l = 7.1447096          ! left u (cm/s)
  v_l = -4.125          ! left v (cm/s)
  rho_l = 8.0             ! left density (g/cc)
  T_l = 1.0

  allocate(p_r)
  allocate(u_r)
  allocate(v_r)
  allocate(rho_r)
  allocate(T_r)

  p_r = 1.0               ! right pressure (erg/cc)
  u_r = 0.0               ! right u (cm/s)
  v_r = 0.0               ! right v (cm/s)
  rho_r = 1.4             ! right density (g/cc)
  T_r = 1.0

  allocate(use_Tinit)

  use_Tinit = .false.     ! optionally use T_l/r instead of p_l/r for initialization

  ! Read namelists
  open(newunit=untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin, fortin)
  close(unit=untin)

  ! set local variable defaults -- the 'center' variables are the location of the

  ! compute the internal energy (erg/cc) for the left and right state
  xn(:) = 0.0e0_rt
  xn(1) = 1.0e0_rt

  allocate(rhoe_l)
  allocate(rhoe_r)
  
  if (use_Tinit) then

     eos_state%rho = rho_l
     eos_state%T = T_l
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state)

     rhoe_l = rho_l*eos_state%e
     p_l = eos_state%p

     eos_state%rho = rho_r
     eos_state%T = T_r
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state)

     rhoe_r = rho_r*eos_state%e
     p_r = eos_state%p

  else

     eos_state%rho = rho_l
     eos_state%p = p_l
     eos_state%T = 10.e0_rt  ! initial guess
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)

     rhoe_l = rho_l*eos_state%e
     T_l = eos_state%T

     eos_state%rho = rho_r
     eos_state%p = p_r
     eos_state%T = 10.e0_rt  ! initial guess
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)

     rhoe_r = rho_r*eos_state%e
     T_r = eos_state%T

  endif

end subroutine amrex_probinit


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.
! :::
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: level     => amr level of grid
! ::: time      => time at which to init data
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       delta, xlo, xhi)

  use network, only: nspec
  use probdata_module
  use prob_params_module, only : problo
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS
  use amrex_constants_module, only : M_PI, SIXTH, HALF, THREE
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt), parameter :: pi_over_3 = M_pi / THREE
  real(rt), parameter :: ff = 0.25d0

  real(rt) :: x,y,xcen,ycen,shockfront
  integer  :: i, j, k, ii,jj

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        ycen = problo(2) + delta(2) * (dble(j) + HALF)

        do i = lo(1), hi(1)
           xcen = problo(1) + delta(1) * (dble(i) + HALF)

           state(i,j,k,URHO   ) = 0.d0
           state(i,j,k,UMX:UMZ) = 0.d0
           state(i,j,k,UEDEN  ) = 0.d0
           state(i,j,k,UEINT  ) = 0.d0

           do jj = -1, 1
              if (jj == 0) cycle
              y = ycen + HALF * delta(2) * (dble(jj) / sqrt(THREE))

              do ii = -1, 1
                 if (ii == 0) cycle
                 x = xcen + HALF * delta(1) * (dble(ii) / sqrt(THREE))

                 shockfront = tan(pi_over_3) * (x - SIXTH) ! initial shock front

                 if (y .ge. shockfront ) then
                    state(i,j,k,URHO) = state(i,j,k,URHO) + ff*rho_l
                    state(i,j,k,UMX ) = state(i,j,k,UMX ) + ff*rho_l*u_l
                    state(i,j,k,UMY ) = state(i,j,k,UMY ) + ff*rho_l*v_l

                    state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + ff*(rhoe_l + 0.5*rho_l*(u_l*u_l + v_l*v_l))
                    state(i,j,k,UEINT) = state(i,j,k,UEINT) + ff*rhoe_l
                    state(i,j,k,UTEMP) = state(i,j,k,UTEMP) + ff*T_l
                 else
                    state(i,j,k,URHO) = state(i,j,k,URHO) + ff*rho_r
                    state(i,j,k,UMX ) = state(i,j,k,UMX ) + ff*rho_r*u_r
                    state(i,j,k,UMY ) = state(i,j,k,UMY ) + ff*rho_r*v_r

                    state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + ff*(rhoe_r + 0.5*rho_r*(u_r*u_r + v_r*v_r))
                    state(i,j,k,UEINT) = state(i,j,k,UEINT) + ff*rhoe_r
                    state(i,j,k,UTEMP) = state(i,j,k,UTEMP) + ff*T_r
                 end if

              end do
           end do

           state(i,j,k,UFS:UFS-1+nspec) = 0.0e0_rt
           state(i,j,k,UFS  ) = state(i,j,k,URHO)

        end do
     end do
  end do

end subroutine ca_initdata
