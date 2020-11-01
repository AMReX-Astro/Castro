subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use probdata_module
  use model_parser_module
  use castro_error_module
  use ambient_module, only: ambient_state
  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module, only: URHO, UTEMP, UFS, UMX, UMZ, UEINT, UEDEN
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_rt

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  type(eos_t) :: eos_state

  ! Read initial model
  call read_model_file(model_name)

  ! set the ambient state for the upper boundary condition
  ambient_state(URHO) = model_state(npts_model, idens_model)
  ambient_state(UTEMP) = model_state(npts_model, itemp_model)
  ambient_state(UFS:UFS-1+nspec) = &
       ambient_state(URHO) * model_state(npts_model, ispec_model:ispec_model-1+nspec)

  ambient_state(UMX:UMZ) = 0.0_rt

  ! make the ambient state thermodynamically consistent
  eos_state % rho = ambient_state(URHO)
  eos_state % T = ambient_state(UTEMP)
  eos_state % xn(:) = ambient_state(UFS:UFS-1+nspec) / eos_state % rho

  call eos(eos_input_rt, eos_state)

  ambient_state(UEINT) = eos_state % rho * eos_state % e
  ambient_state(UEDEN) = eos_state % rho * eos_state % e

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

  use probdata_module
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt, eos_input_tp
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use network, only: nspec
  use model_parser_module
  use prob_params_module, only : problo
  use amrex_constants_module, only : ZERO, ONE, HALF, TWO
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: dist, x, y, z, height
  integer :: i, j, k, n

  real(rt) :: t0,x1,y1,z1,r1,x2,y2,z2,r2,x3,y3,z3,r3,x4,y4,r4,temp

  real(rt) :: temppres(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

  type (eos_t) :: eos_state

  do k = lo(3), hi(3)
     z = problo(3) + delta(3)*(dble(k) + HALF)

     do j = lo(2), hi(2)
        y = problo(2) + delta(2)*(dble(j) + HALF)

#if AMREX_SPACEDIM == 2
        height = y
#else
        height = z
#endif

        do i = lo(1), hi(1)

           call interpolate_sub(state(i,j,k,URHO), height, idens_model)
           call interpolate_sub(state(i,j,k,UTEMP), height, itemp_model)

           do n = 1, nspec
              call interpolate_sub(state(i,j,k,UFS-1+n), height, ispec_model-1+n)
           end do

           eos_state%rho = state(i,j,k,URHO)
           eos_state%T = state(i,j,k,UTEMP)
           eos_state%xn(:) = state(i,j,k,UFS:UFS-1+nspec)

           call eos(eos_input_rt, eos_state)

           temppres(i,j,k) = eos_state%p

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,UEINT)

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           enddo

        end do
     end do
  end do

  ! Initial velocities = 0
  state(:,:,:,UMX:UMZ) = ZERO

  ! Now add the perturbation
  do k = lo(3), hi(3)
     z = problo(3) + delta(3)*(dble(k) + HALF)

     do j = lo(2), hi(2)
        y = problo(2) + delta(2)*(dble(j) + HALF)

        do i = lo(1), hi(1)
           x = problo(1) + delta(1)*(dble(i) + HALF)

           t0 = state(i,j,k,UTEMP)

           x1 = 5.0e7_rt

#if AMREX_SPACEDIM == 2
           y1 = 6.5e7_rt
           z1 = ZERO
#else
           y1 = 5.0e7_rt
           z1 = 6.5e7_rt
#endif
           r1 = sqrt( (x-x1)**2 + (y-y1)**2 + (z-z1)**2 ) / (2.5e6_rt*pert_rad_factor)

           x2 = 1.2e8_rt

#if AMREX_SPACEDIM == 2
           y2 = 8.5e7_rt
           z2 = ZERO
#else
           y2 = 1.2e8_rt
           z2 = 8.5e7_rt
#endif

           r2 = sqrt( (x-x2)**2 + (y-y2)**2 + (z-z2)**2 ) / (2.5e6_rt*pert_rad_factor)

           x3 = 2.0e8_rt

#if AMREX_SPACEDIM == 2
           y3 = 7.5e7_rt
           z3 = ZERO
#else
           y3 = 2.0e8_rt
           z3 = 7.5e7_rt
#endif

           r3 = sqrt( (x-x3)**2 + (y-y3)**2 + (z-z3)**2 ) / (2.5e6_rt*pert_rad_factor)

           state(i,j,k,UTEMP) = t0 * (ONE + pert_temp_factor* &
                (0.150e0_rt * (ONE + tanh(TWO-r1)) + &
                 0.300e0_rt * (ONE + tanh(TWO-r2)) + &
                 0.225e0_rt * (ONE + tanh(TWO-r3))))

           do n = 1,nspec
              state(i,j,k,UFS+n-1) =  state(i,j,k,UFS+n-1) / state(i,j,k,URHO)
           enddo

           eos_state % T = state(i,j,k,UTEMP)
           eos_state % p = temppres(i,j,k)
           eos_state % xn(:) = state(i,j,k,UFS:UFS-1+nspec)
           
           call eos(eos_input_tp, eos_state)

           state(i,j,k,URHO) = eos_state%rho

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,UEINT)

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           enddo

        enddo
     enddo
  enddo

end subroutine ca_initdata
