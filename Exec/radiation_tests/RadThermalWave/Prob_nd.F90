subroutine amrex_probinit(init,name,namlen,problo,probhi) bind(c)

  use amrex_constants_module
  use probdata_module
  use amrex_fort_module, only : rt => amrex_real
  use castro_error_module

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer untin, i

  namelist /fortin/ rhocv, T0, Eexp, rexp

  ! Build "probin" filename -- the name of file containing fortin namelist.

  integer, parameter :: maxlen = 127
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call castro_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! Set defaults
  rhocv = -1.0e50_rt
  T0 = -1.0e50_rt
  Eexp = -1.e-50_rt
  rexp = -1.e50_rt

  ! domain extrema and center
  xmin = problo(1)
  xmax = probhi(1)
#if AMREX_SPACEDIM >= 2
  ymin = problo(2)
  ymax = probhi(2)
#endif
#if AMREX_SPACEDIM == 3
  zmin = problo(3)
  zmax = probhi(3)
#endif

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin,fortin)
  close(unit=untin)

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
                       state, state_lo, state_hi, &
                       delta, xlo, xhi)

  use amrex_constants_module
  use probdata_module
  use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS

  use network, only : nspec
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use prob_params_module, only : dg, coord_type
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1), &
                                   state_lo(2):state_hi(2), &
                                   state_lo(3):state_hi(3), nscal)

  integer i, j, k, ii, jj, kk
  type(eos_t) :: eos_state
  real(rt) :: rho, cv, T, p, eint
  real(rt) :: rhoeexp, Vexp, rhoe0
  real(rt) :: xx, xcl, xcr, dx_sub
  real(rt) :: yy, ycl, ycr, dy_sub
  real(rt) :: zz, zcl, zcr, dz_sub
  real(rt) :: rr2, xcmin, xcmax, ycmin, ycmax, zcmin, zcmax , rcmin, rcmax
  real(rt) :: vol_pert, vol_ambient, T_zone, rhoe_zone

  integer, parameter :: nsub = 64

  Vexp = 4.e0_rt/3.e0_rt*Pi * rexp**3

  rhoeexp = Eexp / Vexp

  rho = 1.0e0_rt
  T = 1.0e0_rt

  eos_state % rho = rho
  eos_state % T   = T
  eos_state % xn  = ZERO
  eos_state % xn(1)  = 1.e0_rt

  call eos(eos_input_rt, eos_state)

  cv = eos_state % cv

  rho = rhocv / cv

  eint = rhoeexp / rho
  T = eint / cv

  rhoe0 = rho * cv * T0

  dx_sub = delta(1)/dble(nsub)
  dy_sub = delta(2)/dble(nsub)
  dz_sub = delta(3)/dble(nsub)

  do k = lo(3), hi(3)

     zcl = xlo(3) + delta(3) * dble(k-lo(3))
     zcr = zcl + delta(3)
     zcmin = min(abs(zcl), abs(zcr))
     zcmax = max(abs(zcl), abs(zcr))

     do j = lo(2), hi(2)

        ycl = xlo(2) + delta(2)*dble(j-lo(2))
        ycr = ycl + delta(2)
        ycmin = min(abs(ycl), abs(ycr))
        ycmax = max(abs(ycl), abs(ycr))

        do i = lo(1), hi(1)

           xcl = xlo(1) + delta(1)*dble(i-lo(1))
           xcr = xcl + delta(1)
           xcmin = min(abs(xcl), abs(xcr))
           xcmax = max(abs(xcl), abs(xcr))

           rcmin = sqrt(xcmin**2+ycmin**2+zcmin**2)
           rcmax = sqrt(xcmax**2+ycmax**2+zcmax**2)

           state(i,j,k,URHO) = rho
           state(i,j,k,UMX)  = ZERO
           state(i,j,k,UMY)  = ZERO
           state(i,j,k,UMZ)  = ZERO

           if (rcmin .ge. rexp) then
              T_zone = T0
              rhoe_zone = rhoe0
           else if (rcmax .le. rexp) then
              T_zone = T
              rhoe_zone = rhoeexp
           else
              vol_pert    = ZERO
              vol_ambient = ZERO

              do kk = 0, dg(3)*(nsub-1)
                 zz = zcl + (dble(kk) + HALF) * dz_sub

                 do jj = 0, dg(2)*(nsub-1)
                    yy = ycl + (dble(jj) + HALF) * dy_sub

                    do ii = 0, nsub-1
                       xx = xcl + (dble(ii) + HALF) * dx_sub

                       rr2 = xx**2 + yy**2 + zz**2

                       if (coord_type == 1) then
                          ! we are axisymmetric so weight by the radial coordinate
                          if (rr2 <= rexp**2) then
                             vol_pert = vol_pert + xx
                          else
                             vol_ambient = vol_ambient + xx
                          endif
                       else
                          ! Cartesian
                          if (rr2 <= rexp**2) then
                             vol_pert = vol_pert + ONE
                          else
                             vol_ambient = vol_ambient + ONE
                          endif
                       endif

                    end do
                 end do
              end do

              T_zone = (vol_pert*T + vol_ambient*T0)/(vol_pert+vol_ambient)
              rhoe_zone = (vol_pert*rhoeexp + vol_ambient*rhoe0)/(vol_pert+vol_ambient)

           end if

           state(i,j,k,UTEMP) = T_zone
           state(i,j,k,UEDEN) = rhoe_zone
           state(i,j,k,UEINT) = rhoe_zone

           state(i,j,k,UFS) = state(i,j,k,URHO)
        enddo
     enddo
  end do

end subroutine ca_initdata


subroutine ca_initrad(level, time, lo, hi, nrad, &
                      rad_state, rad_state_lo, rad_state_hi, &
                      delta, xlo, xhi)

  use probdata_module
  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: level, nrad
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: rad_state_lo(3), rad_state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: rad_state(rad_state_lo(1):rad_state_hi(1), &
                                       rad_state_lo(2):rad_state_hi(2), &
                                       rad_state_lo(3):rad_state_hi(3), 0:nrad-1)

  rad_state(:,:,:,:) = ZERO

end subroutine ca_initrad
