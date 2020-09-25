subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use amrex_constants_module
  use probdata_module
  use prob_params_module, only : center
  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  use eos_type_module, only : eos_t, eos_input_rt
  use eos_module, only : eos
  use extern_probin_module, only : small_x
  use network

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  real(rt) :: Ye, abar, dq, dyedt
  real (rt) :: xn(nspec)
  type(eos_t) :: eos_state

  ! set explosion center
  center(:) = ZERO
  center(1) = HALF*(problo(1) + probhi(1))
#if BL_SPACEDIM >= 2
  center(2) = HALF*(problo(2) + probhi(2))
#endif
#if BL_SPACEDIM == 3
  center(3) = HALF*(problo(3) + probhi(3))
#endif

  ! we only work in NSE mode, so the mass fractions are not used by the EOS
#ifndef NSE_THERMO
  call castro_error("Error: this problem requires USE_NSE=TRUE")
#endif

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
subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_lo,state_hi, &
                       delta,xlo,xhi)

  use probdata_module
  use amrex_constants_module, only: M_PI, FOUR3RD, ZERO, HALF, ONE
  use meth_params_module , only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UFX, UTEMP
  use prob_params_module, only : center, coord_type, problo, probhi
  use amrex_fort_module, only : rt => amrex_real
  use network
  use eos_type_module, only : eos_t, eos_input_rt
  use eos_module, only : eos
  use nse_module

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  real(rt) :: xx, yy, zz
  real(rt) :: dist, T, eint
  real(rt) :: ye0, dye, ye, abar, dq, dyedt
  real(rt) :: xn(nspec)
  integer :: i, j, k

  type(eos_t) :: eos_state

  ye0 = 0.5
  dye = -0.05

  do k = lo(3), hi(3)
     zz = problo(3) + delta(3)*(dble(k) + HALF)

     do j = lo(2), hi(2)
        yy = problo(2) + delta(2)*(dble(j) + HALF)

        do i = lo(1), hi(1)
           xx = problo(1) + delta(1)*(dble(i) + HALF)

           dist = sqrt((center(1)-xx)**2 + (center(2)-yy)**2 + (center(3)-zz)**2)

           if (dist <= center(1)) then
              T = T0*(ONE + dT_fact*exp(-(dist/L_pert)**2) * &
                   cos(M_PI*(dist/(probhi(1)-problo(1))))**6)
              ye = ye0*(ONE + dye*exp(-(dist/L_pert)**2) * &
                   cos(M_PI*(dist/(probhi(1)-problo(1))))**6)
           else
              ye = ye0
           endif

           state(i,j,k,UMX:UMZ) = 0.e0_rt

           ! we are isentropic, so find rho
           eos_state % T = T
           eos_state % rho = rho0  ! initial guess

           call nse_interp(T, rho0, Ye, abar, dq, dyedt, xn)

           ! since the species are interpolated, normalize them
           xn(:) = xn(:) / sum(xn(:))

           eos_state % xn(:) = xn(:)
           eos_state % aux(iye) = Ye
           eos_state % aux(iabar) = abar
           eos_state % aux(ibea) = dq

           call eos(eos_input_rt, eos_state)

           state(i,j,k,URHO) = eos_state % rho

           state(i,j,k,UEDEN) = eos_state % rho * eos_state % e
           state(i,j,k,UEINT) = eos_state % rho * eos_state % e

           state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,URHO) * xn(:)

           state(i,j,k,UFX:UFX-1+naux) = state(i,j,k,URHO) * eos_state % aux(:)

           state(i,j,k,UTEMP) = eos_state % T
        enddo
     enddo
  enddo

end subroutine ca_initdata
