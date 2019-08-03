subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(C, name="amrex_probinit")

  use probdata_module
  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer untin, i
  character(1) dummy

  namelist /fortin/ &
       rwind0, rwind1, rhowind1, Twind1, rbasefac, filter_rhomax, filter_timemax, model_file

  ! Build "probin" filename -- the name of file containing fortin namelist.

  integer, parameter :: maxlen=127
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call castro_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  rbasefac = 0.99e0_rt
  rwind0 = 0.7e14_rt
  rwind1 = 1.e14_rt
  rhowind1 = 1.e-14_rt
  Twind1 = 1.1e3_rt

  filter_rhomax = -1.e20_rt
  filter_timemax = -1.e20_rt

  xmin = problo(1)
  xmax = probhi(1)

  model_file = "model.input"

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin,fortin)
  close(unit=untin)

  open(newunit=untin, file=trim(model_file))
  print*,'reading model inputs'
  read(untin,*) dummy
  read(untin,*) npts_model
  read(untin,*) dummy
  if (npts_model > npts_max) then
     call castro_error('npts_max in probdata.f90 is too small')
  end if

  do i = 1, npts_model
     read(untin,*)model_r(i), model_rho(i), model_v(i), &
          model_T(i), model_Ye(i), model_Abar(i)
  enddo

  print *,'done reading model inputs'

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

  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UFX, UTEMP
  use network, only : nspec, naux
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use interpolate_module
  use fundamental_constants_module, only: k_B, n_A
  use prob_params_module, only : problo

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1), &
                                   state_lo(2):state_hi(2), &
                                   state_lo(3):state_hi(3), nscal)

  integer :: i, j, k, ii
  real(rt) :: rhoInv, xcl, xx, vtot, vsub, rho, T, u, rhosub, Tsub, usub, dx_sub
  real(rt) :: rhowind0, rlast, rholast, Twind0, Tlast
  integer, parameter :: nsub = 16
  real(rt) :: rho_tmp, rbase, T_tmp
  real(rt) :: Ye, Abar, invmu
  real(rt), parameter :: Tindex=0.5
  type(eos_t) :: eos_state

  if (naux .ne. 2) then
     call castro_error("naux in network is not equal to 2")
  end if

  if (nspec .ne. 1) then
     call castro_error("nspec in network is not equal to 1")
  end if

  dx_sub = delta(1) / dble(nsub)

  rhowind0 = rhowind1 * (rwind1/rwind0)**2
  rlast = model_r(npts_model)
  rholast = model_rho(npts_model)
  rbase = rlast * rbasefac

  Twind0 = Twind1 * (rwind1/rwind0)**Tindex
  Tlast = model_T(npts_model)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           xcl = problo(1) + delta(1) * dble(i)

           vtot = 0.e0_rt
           rho = 0.e0_rt
           T = 0.e0_rt
           u = 0.e0_rt
           Ye = 0.e0_rt
           Abar = 0.e0_rt

           do ii=0,nsub-1
              xx = xcl + (dble(ii)+0.5e0_rt) * dx_sub
              vsub = xx**2
              vtot = vtot + vsub
              if (xx .ge. model_r(npts_model)) then
                 if (xx .ge. rwind0 ) then
                    rho_tmp = rhowind1 * (rwind1/xx)**2
                 else
                    rho_tmp = rholast * (rhowind0/rholast)** &
                         ((log(xx-rbase)-log(rlast-rbase))/(log(rwind0-rbase)-log(rlast-rbase)))
                 end if
                 rho = rho + rho_tmp * vsub

                 if (xx .ge. rwind0 ) then
                    T_tmp = Twind1 * (rwind1/xx)**Tindex
                 else
                    T_tmp = Tlast * (Twind0/Tlast)** &
                         ((log(xx-rbase)-log(rlast-rbase))/(log(rwind0-rbase)-log(rlast-rbase)))
                 end if
                 T = T + vsub * T_tmp

                 Ye = Ye + model_Ye(npts_model) * vsub
                 Abar = Abar + model_Abar(npts_model) * vsub

                 u = u + 0.e0_rt

              else if (xx .le. model_r(1)) then
                 rho = rho + model_rho(1) * vsub
                 T = T + model_T(1) * vsub
                 u = u + 0.e0_rt
                 Ye = Ye + model_Ye(1) * vsub
                 Abar = Abar + model_Abar(1) * vsub
              else
                 rho = rho + interpolate(xx,npts_model,model_r,model_rho) * vsub
                 T   = T   + interpolate(xx,npts_model,model_r,model_T  ) * vsub
                 u   = u   + interpolate(xx,npts_model,model_r,model_v  ) * vsub
                 Ye  = Ye  + interpolate(xx,npts_model,model_r,model_Ye ) * vsub
                 Abar=Abar + interpolate(xx,npts_model,model_r,model_Abar) * vsub
              end if
           end do

           rho = rho / vtot
           T = T / vtot
           u = u / vtot
           Ye = Ye / vtot
           Abar = Abar / vtot

           invmu = (1.e0_rt+Abar*Ye)/Abar

           state(i,j,k,URHO)  = rho
           state(i,j,k,UTEMP)  = T
           state(i,j,k,UMX)   = rho * u
           state(i,j,k,UMY:UMZ) = 0.e0_rt

           ! set the composition to be all in the first species
           state(i,j,k,UFS:UFS-1+nspec) = 0.e0_rt
           state(i,j,k,UFS  ) = state(i,j,k,URHO)
           state(i,j,k,UFX) = Ye*rho
           state(i,j,k,UFX+1) = invmu*rho

           ! set the internal energy via the EOS
           rhoInv = 1.e0_rt / state(i,j,k,URHO)
           eos_state % rho = state(i,j,k,URHO)
           eos_state % T   = state(i,j,k,UTEMP)
           eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
           eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rhoInv

           call eos(eos_input_rt, eos_state)

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                0.5e0_rt*sum(state(i,j,k,UMX:UMZ)**2)/state(i,j,k,URHO)
        end do
     end do
  end do

end subroutine ca_initdata


! :::
! ::: -----------------------------------------------------------
! :::

subroutine ca_initrad(level, time, lo, hi, nrad, &
                      rad_state, rad_state_lo, rad_state_hi, &
                      delta, xlo, xhi)

  use probdata_module
  use fundamental_constants_module, only: a_rad
  use interpolate_module
  use rad_params_module, only : xnu
  use blackbody_module, only : BGroup

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: level, nrad
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: rad_state_lo(3), rad_state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: rad_state(rad_state_lo(1):rad_state_hi(1), &
                                       rad_state_lo(2):rad_state_hi(2), &
                                       rad_state_lo(3):rad_state_hi(3), 0:nrad-1)

  ! local variables
  integer :: i, j, k, ii, igroup
  real(rt) :: xcl, T, Tsub, xx, vtot, vsub, dx_sub
  integer, parameter :: nsub = 16

  real(rt) :: rlast, rbase, T_tmp, Twind0, Tlast
  real(rt), parameter :: Tindex=0.5

  dx_sub = delta(1) / dble(nsub)

  rlast = model_r(npts_model)
  rbase = rlast * rbasefac

  Twind0 = Twind1 * (rwind1/rwind0)**Tindex
  Tlast = model_T(npts_model)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           xcl = problo(1) + delta(1) * dble(i)

           vtot = 0.e0_rt
           T = 0.e0_rt

           do ii = 0, nsub-1
              xx = xcl + (dble(ii)+0.5e0_rt) * dx_sub
              vsub = xx**2
              vtot = vtot + vsub
              if (xx .ge. model_r(npts_model)) then
                 if (xx .ge. rwind0 ) then
                    T_tmp = Twind1 * (rwind1/xx)**Tindex
                 else
                    T_tmp = Tlast * (Twind0/Tlast)** &
                         ((log(xx-rbase)-log(rlast-rbase))/(log(rwind0-rbase)-log(rlast-rbase)))
                 end if
                 T = T + vsub * T_tmp

              else if (xx .le. model_r(1)) then
                 T = T + model_T(1) * vsub
              else
                 T = T + interpolate(xx,npts_model,model_r,model_T  ) * vsub
              end if
           end do
           T = T / vtot
           ! set radiation energy density to a T**4

           if (nrad .eq. 1) then
              rad_state(i,j,k,0) = a_rad*T**4
           else
              do igroup=0,nrad-1
                 rad_state(i,j,k,igroup) = BGroup(T, xnu(igroup), xnu(igroup+1))
              end do
           end if

        end do
     end do
  end do

end subroutine ca_initrad
