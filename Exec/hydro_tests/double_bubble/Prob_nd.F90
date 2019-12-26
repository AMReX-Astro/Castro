subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use prob_params_module, only: center
  use probdata_module
  use castro_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  namelist /fortin/ dens_base, pres_base, &
                    pert_factor, r_pert_center, pert_width, &
                    do_isentropic, &
                    boundary_type, &
                    single


  call probdata_init(name, namlen)

  ! model composition
  xn_model(:) = 0.0e0_rt
  xn_model(1) = 1.0e0_rt

  ! set local variable defaults
  center(1) = 0.5e0_rt*(problo(1)+probhi(1))
  center(2) = 0.5e0_rt*(problo(2)+probhi(2))
  center(3) = 0.5e0_rt*(problo(3)+probhi(3))

#if AMREX_SPACEDIM == 1
  rmin = problo(1)
  rmax = probhi(1)
#elif AMREX_SPACEDIM == 2
  rmin = problo(2)
  rmax = probhi(2)
#else
  rmin = problo(3)
  rmax = probhi(3)
#endif

  if (single) then
     left_bubble_x_center = problo(1) + 0.5e0_rt*(probhi(1)-problo(1))
  else
     left_bubble_x_center = problo(1) + (probhi(1)-problo(1))/3.e0_rt
     right_bubble_x_center = problo(1) + 2.e0_rt*(probhi(1)-problo(1))/3.e0_rt
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

  use amrex_constants_module, only : ZERO, HALF, ONE, TWO
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP, &
       const_grav
  use eos_module
  use eos_type_module
  use network
  use prob_params_module, only : problo, center
  use model_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  integer i, j, k, npts_1d
  real(rt) :: xn(nspec), x, y, z, x1, y1, z1, x2, y2, z2, r1, r2, const

  real(rt), allocatable :: r_model(:), rho_model(:), T_model(:), &
                           e_model(:), p_model(:)

  integer :: lo_model, hi_model, idx_height

  type (eos_t) :: eos_state

  ! we'll generate the initial model at the needed resolution
  call get_model_size(rmin, rmax, delta(AMREX_SPACEDIM), lo_model, hi_model)

  allocate(  r_model(lo_model:hi_model))
  allocate(rho_model(lo_model:hi_model))
  allocate(  T_model(lo_model:hi_model))
  allocate(  e_model(lo_model:hi_model))
  allocate(  p_model(lo_model:hi_model))

  call get_model(rmin, rmax, delta(AMREX_SPACEDIM), &
                 pres_base, dens_base, do_isentropic, &
                 xn_model, &
                 r_model, rho_model, T_model, e_model, p_model, &
                 lo_model, hi_model)

  if (.not. single) then

#if AMREX_SPACEDIM == 1
     call castro_error("Error: 1-d not supported")

#elif AMREX_SPACEDIM == 2
     x1 = left_bubble_x_center
     y1 = r_pert_center
     z1 = ZERO

     x2 = right_bubble_x_center
     y2 = r_pert_center
     z2 = ZERO

#else
     x1 = left_bubble_x_center
     y1 = center(2)
     z1 = r_pert_center

     x2 = right_bubble_x_center
     y2 = center(2)
     z2 = r_pert_center
#endif

     do k = lo(3), hi(3)
        z = (dble(k)+HALF)*delta(3) + problo(3)

        do j=lo(2),hi(2)
           y = (dble(j)+HALF)*delta(2) + problo(2)

           do i=lo(1),hi(1)
              x = (dble(i)+HALF)*delta(1) + problo(1)

              r1 = sqrt( (x-x1)**2 + (y-y1)**2 + (z-z1)**2) / pert_width
              r2 = sqrt( (x-x2)**2 + (y-y2)**2 + (z-z2)**2) / pert_width

#if AMREX_SPACEDIM == 2
              idx_height = j
#else
              idx_height = k
#endif

              eos_state % xn(:) = xn_model(:)
              state(i,j,k,UTEMP) = T_model(idx_height)
              state(i,j,k,URHO) = rho_model(idx_height)

              ! which bubble are we in? -- we want their rho perturbations to be the
              ! same so they have the same buoyancy

              if (r1 < 2.0e0_rt) then
                 state(i,j,k,URHO) = rho_model(idx_height) * (ONE - (pert_factor * (ONE + tanh(TWO-r1))))
                 eos_state % xn(:) = ZERO
                 eos_state % xn(2) = ONE
              endif

              if (r2 < 2.0e0_rt) then
                 state(i,j,k,URHO) = rho_model(idx_height) * (ONE - (pert_factor * (ONE + tanh(TWO-r2))))
                 eos_state % xn(:) = ZERO
                 eos_state % xn(3) = ONE
              endif

              eos_state % p = p_model(idx_height)
              eos_state % rho = state(i,j,k,URHO)

              call eos(eos_input_rp, eos_state)

              state(i,j,k,UEINT) = eos_state % e
              state(i,j,k,UTEMP) = eos_state % T

              ! make state conservative
              state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,URHO)*eos_state % xn(:)
              state(i,j,k,UEINT) = state(i,j,k,URHO)*state(i,j,k,UEINT)

              ! assumes ke=0
              state(i,j,k,UEDEN) = state(i,j,k,UEINT)

              state(i,j,k,UMX:UMZ) = ZERO

           end do
        end do
     end do
  else

#if AMREX_SPACEDIM == 1
     call castro_error("Error: 1-d not supported")
#elif AMREX_SPACEDIM == 2
     x1 = left_bubble_x_center
     y1 = r_pert_center
     z1 = center(3)
#else
     x1 = left_bubble_x_center
     y1 = center(2)
     z1 = r_pert_center
#endif

     do k = lo(3), hi(3)
        z = (dble(j)+HALF)*delta(3) + problo(3)

        do j = lo(2), hi(2)
           y = (dble(j)+HALF)*delta(2) + problo(2)

           do i = lo(1), hi(1)
              x = (dble(i)+HALF)*delta(1) + problo(1)

              r1 = sqrt( (x-x1)**2 + (y-y1)**2 + (z-z1)**2 ) / pert_width

#if AMREX_SPACEDIM == 2
              idx_height = j
#else
              idx_height = k
#endif

              eos_state % xn(:) = xn_model(:)
              state(i,j,k,UTEMP) = T_model(idx_height)
              state(i,j,k,URHO) = rho_model(idx_height)

              ! which bubble are we in? -- we want their rho perturbations to be the
              ! same so they have the same buoyancy
              if (r1 < 2.0e0_rt) then
                 state(i,j,k,URHO) = rho_model(idx_height) * (ONE - (pert_factor * (ONE + tanh(TWO-r1))))
                 eos_state % xn(:) = ZERO
                 eos_state % xn(2) =ONE
              endif

              eos_state % p = p_model(idx_height)
              eos_state % rho = state(i,j,k,URHO)

              call eos(eos_input_rp, eos_state)

              state(i,j,k,UEINT) = eos_state % e
              state(i,j,k,UTEMP) = eos_state % T


              ! make state conservative
              state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,URHO)*eos_state % xn(:)
              state(i,j,k,UEINT) = state(i,j,k,URHO)*state(i,j,k,UEINT)

              ! assumes ke=0
              state(i,j,k,UEDEN) = state(i,j,k,UEINT)

              state(i,j,k,UMX:UMZ) = ZERO

           end do
        end do
     end do

  endif

  deallocate(r_model, rho_model, T_model, p_model, e_model)

end subroutine ca_initdata


