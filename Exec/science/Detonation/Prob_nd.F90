subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use probdata_module, only: T_l, T_r, dens, cfrac, ofrac, idir, w_T, center_T, &
                             xn, ihe4, ic12, io16, smallx, vel, grav_acceleration, fill_ambient_bc, &
                             ambient_dens, ambient_temp, ambient_comp, ambient_e_l, ambient_e_r
  use network, only: network_species_index, nspec
  use castro_error_module, only: castro_error
  use amrex_fort_module, only: rt => amrex_real
  use eos_type_module, only: eos_t, eos_input_rt
  use eos_module, only: eos

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  type(eos_t) :: eos_state

  call probdata_init(name, namlen)

  ! get the species indices
  ihe4 = network_species_index("helium-4")
  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")

  if (ihe4 < 0 .or. ic12 < 0 .or. io16 < 0) then
     call castro_error("ERROR: species indices not found")
  endif

  ! make sure that the carbon fraction falls between 0 and 1
  if (cfrac > 1.e0_rt .or. cfrac < 0.e0_rt) then
     call castro_error("ERROR: cfrac must fall between 0 and 1")
  endif

  ! make sure that the oxygen fraction falls between 0 and 1
  if (ofrac > 1.e0_rt .or. cfrac < 0.e0_rt) then
     call castro_error("ERROR: ofrac must fall between 0 and 1")
  endif

  ! make sure that the C/O fraction sums to no more than 1
  if (cfrac + ofrac > 1.e0_rt) then
     call castro_error("ERROR: cfrac + ofrac cannot exceed 1.")
  end if

  ! set the default mass fractions

  xn(:) = smallx
  xn(ic12) = max(cfrac, smallx)
  xn(io16) = max(ofrac, smallx)
  xn(ihe4) = 1.e0_rt - cfrac - ofrac - (nspec - 2) * smallx

  ! Set the ambient material

  ambient_dens = dens
  ambient_comp = xn

  eos_state % rho = ambient_dens
  eos_state % xn  = ambient_comp

  eos_state % T   = T_l

  call eos(eos_input_rt, eos_state)

  ambient_e_l = eos_state % e

  eos_state % T   = T_r

  call eos(eos_input_rt, eos_state)

  ambient_e_r = eos_state % e

end subroutine amrex_probinit


subroutine ca_initdata(lo, hi, &
                       state, state_lo, state_hi, &
                       dx, problo) bind(c, name='ca_initdata')

  use network, only: nspec
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_rt
  use probdata_module, only: T_l, T_r, center_T, w_T, dens, vel, xn
  use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
  use amrex_fort_module, only: rt => amrex_real
  use prob_params_module, only: probhi

  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: state_lo(3), state_hi(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
  real(rt), intent(in   ) :: dx(3), problo(3)

  real(rt) :: sigma, width, c_T
  real(rt) :: xcen
  integer  :: i, j, k, n

  type (eos_t) :: eos_state

  !$gpu

  width = w_T * (probhi(1) - problo(1))
  c_T = problo(1) + center_T * (probhi(1) - problo(1))

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           xcen = problo(1) + dx(1)*(dble(i) + 0.5e0_rt)

           state(i,j,k,URHO) = dens

           sigma = 1.0 / (1.0 + exp(-(c_T - xcen)/ width))

           state(i,j,k,UTEMP) = T_l + (T_r - T_l) * (1 - sigma)

           do n = 1, nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * xn(n)
           end do

           eos_state%rho = state(i,j,k,URHO)
           eos_state%T = state(i,j,k,UTEMP)
           eos_state%xn(:) = xn

           call eos(eos_input_rt, eos_state)

           state(i,j,k,UMX  ) = state(i,j,k,URHO) * (vel - 2 * vel * (1.0e0_rt - sigma))
           state(i,j,k,UMY  ) = 0.e0_rt
           state(i,j,k,UMZ  ) = 0.e0_rt
           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state%e
           state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5e0_rt * sum(state(i,j,k,UMX:UMZ)**2) / state(i,j,k,URHO)
        enddo
     enddo
  enddo

end subroutine ca_initdata
