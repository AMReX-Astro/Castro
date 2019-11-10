subroutine ca_initdata_maestro(lo,hi,MAESTRO_init_type, &
                               state,state_l1,state_l2,state_h1,state_h2, &
                               dx,dr,xlo,xhi,p0,MAESTRO_npts_model,level)

  use probdata_module
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP
  use network, only: nspec

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer lo(2), hi(2), MAESTRO_init_type, level
  integer MAESTRO_npts_model
  integer state_l1,state_l2,state_h1,state_h2
  real(rt)         xlo(2), xhi(2), dx(2), dr
  real(rt)         p0(0:MAESTRO_npts_model-1)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,NVAR)

  ! local variables
  real(rt)         ekin
  real(rt)         pressure,entropy,minpres

  integer i,j,n

  type (eos_t) eos_state

  minpres = p0(MAESTRO_npts_model-1)

  ! compute p0 and add pi if necessary
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        if (MAESTRO_init_type .eq. 1) then
           ! set pressure = p0
           state(i,j,UEDEN) = p0(j)
        else
           ! set pressure = p0+pi
           state(i,j,UEDEN) = state(i,j,UEDEN) + p0(j)
        end if

     end do
  end do

  if (MAESTRO_init_type .eq. 1 .or. MAESTRO_init_type .eq. 2) then

     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           ! load pressure from our temporary storage field
           pressure = max(state(i,j,UEDEN),minpres)

           ! compute e and T
           eos_state % p   = pressure
           eos_state % rho = state(i,j,URHO)
           eos_state % T   = state(i,j,UTEMP)
           eos_state % xn  = state(i,j,UFS:UFS+nspec-1)

           call eos_on_host(eos_input_rp, eos_state)

           state(i,j,UTEMP) = state(i,j,UTEMP)
           state(i,j,UEINT) = eos_state % e

           ! compute kinetic energy
           ekin = 0.5*state(i,j,URHO)*(state(i,j,UMX)**2+state(i,j,UMY)**2)

           ! convert velocity to momentum
           state(i,j,UMX:UMY) = state(i,j,UMX:UMY)*state(i,j,URHO)

           ! compute rho*e
           state(i,j,UEINT) = state(i,j,URHO) * state(i,j,UEINT)

           ! compute rho*E = rho*e + ke
           state(i,j,UEDEN) = state(i,j,UEINT) + ekin

           ! convert X to rhoX
           do n = 1,nspec
              state(i,j,UFS+n-1) = state(i,j,URHO) * state(i,j,UFS+n-1)
           end do

        end do
     end do

  else if (MAESTRO_init_type .eq. 3) then

     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           ! load pressure from our temporary storage field
           pressure = max(state(i,j,UEDEN),minpres)

           ! compute rho and e
           eos_state % T  = state(i,j,UTEMP)
           eos_state % p  = pressure
           eos_state % xn = state(i,j,UFS:UFS+nspec-1)

           call eos_on_host(eos_input_pt, eos_state)

           state(i,j,URHO)  = eos_state % rho
           state(i,j,UEINT) = eos_state % e

           ! compute kinetic energy
           ekin = 0.5*state(i,j,URHO)*(state(i,j,UMX)**2+state(i,j,UMY)**2)

           ! convert velocity to momentum
           state(i,j,UMX:UMY) = state(i,j,UMX:UMY)*state(i,j,URHO)

           ! compute rho*e
           state(i,j,UEINT) = state(i,j,URHO) * state(i,j,UEINT)

           ! compute rho*E = rho*e + ke
           state(i,j,UEDEN) = state(i,j,UEINT) + ekin

           ! convert X to rhoX
           do n = 1,nspec
              state(i,j,UFS+n-1) = state(i,j,URHO) * state(i,j,UFS+n-1)
           end do

        end do
     end do

  else if (MAESTRO_init_type .eq. 4) then

     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           ! load pressure from our temporary storage field
           pressure = max(state(i,j,UEDEN),minpres)

           ! load entropy from our temporary storage field
           entropy = state(i,j,UEINT)

           ! compute kinetic energy
           ekin = 0.5*state(i,j,URHO)*(state(i,j,UMX)**2+state(i,j,UMY)**2)

           ! compute rho, T, and e
           eos_state % p  = pressure
           eos_state % s  = entropy
           eos_state % xn = state(i,j,UFS:UFS+nspec-1)

           call eos_on_host(eos_input_ps, eos_state)

           state(i,j,URHO)  = eos_state % rho
           state(i,j,UEINT) = eos_state % e
           state(i,j,UTEMP) = eos_state % T

           ! convert velocity to momentum
           state(i,j,UMX:UMY) = state(i,j,UMX:UMY)*state(i,j,URHO)

           ! compute rho*e
           state(i,j,UEINT) = state(i,j,URHO) * state(i,j,UEINT)

           ! compute rho*E = rho*e + ke
           state(i,j,UEDEN) = state(i,j,UEINT) + ekin

           ! convert X to rhoX
           do n = 1,nspec
              state(i,j,UFS+n-1) = state(i,j,URHO) * state(i,j,UFS+n-1)
           end do

        end do
     end do

  end if

end subroutine ca_initdata_maestro

subroutine ca_initdata_makemodel(model,model_size,MAESTRO_npts_model, &
                                 rho0,tempbar,dx,dr,r_model_start)

  use network, only: nspec
  use eos_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer model_size,MAESTRO_npts_model
  real(rt)         model(model_size,0:MAESTRO_npts_model-1)
  real(rt)         rho0   (0:MAESTRO_npts_model-1)
  real(rt)         tempbar(0:MAESTRO_npts_model-1)
  real(rt)         dx(2), dr
  integer r_model_start

  ! local
  integer i,n,iter,comp
  integer MAX_ITER

  real(rt)         TOL

  real(rt)         pres(0:MAESTRO_npts_model-1)
  real(rt)         entropy_want,entropy
  real(rt)         dens_zone,temp_zone,xn(nspec)

  logical converged_hse,isentropic,fluff
  real(rt)         g_zone,low_density_cutoff,temp_fluff

  real(rt)         p_want,pres_zone,dpt,dpd,dst,dsd,A,B,drho,dtemp

  type (eos_t) eos_state

  MAX_ITER = 250
  TOL = 1.e-10

  g_zone = -1.5e10_rt

  low_density_cutoff = 1.e-4_rt
  temp_fluff = 1.e7_rt

  !-----------------------------------------------------------------------------
  ! put the model onto our new uniform grid
  !-----------------------------------------------------------------------------

  fluff = .false.

  ! hard code species
  xn(1) = 0.3e0_rt
  xn(2) = 0.7e0_rt
  xn(3) = 0.e0_rt

  eos_state % rho  = rho0(r_model_start)
  eos_state % T    = tempbar(r_model_start)
  eos_state % xn   = xn(:)

  call eos_on_host(eos_input_rt, eos_state)

  model(1,r_model_start) = rho0(r_model_start)
  model(2,r_model_start) = tempbar(r_model_start)

  pres(r_model_start) = eos_state % p
  entropy_want = eos_state % s

  !-----------------------------------------------------------------------------
  ! HSE + entropy solve
  !-----------------------------------------------------------------------------

  ! the HSE state will be done putting creating an isentropic state until
  ! the temperature goes below temp_fluff -- then we will do isothermal.
  ! also, once the density goes below low_density_cutoff, we stop HSE

  isentropic = .true.

  !---------------------------------------------------------------------------
  ! integrate up
  !---------------------------------------------------------------------------
  do i=r_model_start+1,MAESTRO_npts_model-1

     ! as the initial guess for the temperature and density, use the previous
     ! zone
     dens_zone = model(1,i-1)
     temp_zone = model(2,i-1)

     !-----------------------------------------------------------------------
     ! iteration loop
     !-----------------------------------------------------------------------

     ! start off the Newton loop by saying that the zone has not converged
     converged_hse = .FALSE.

     if (.not. fluff) then

        do iter = 1, MAX_ITER

           if (isentropic) then

              ! get the pressure we want from the HSE equation, just the
              ! zone below the current.  Note, we are using an average of
              ! the density of the two zones as an approximation of the
              ! interface value -- this means that we need to iterate for
              ! find the density and pressure that are consistent

              ! furthermore, we need to get the entropy that we need,
              ! which will come from adjusting the temperature in
              ! addition to the density.

              ! HSE differencing
              p_want = pres(i-1) + &
                   dr*0.5*(dens_zone + model(1,i-1))*g_zone

              ! now we have two functions to zero:
              !   A = p_want - p(rho,T)
              !   B = entropy_want - s(rho,T)
              ! We use a two dimensional Taylor expansion and find the deltas
              ! for both density and temperature


              ! now we know the pressure and the entropy that we want, so we
              ! need to find the temperature and density through a two
              ! dimensional root find

              ! (t, rho) -> (p, s)
              eos_state % T   = temp_zone
              eos_state % rho = dens_zone
              eos_state % xn  = xn(:)

              call eos_on_host(eos_input_rt, eos_state)

              entropy = eos_state % s
              pres_zone = eos_state % p

              dpt = eos_state % dpdt
              dpd = eos_state % dpdr
              dst = eos_state % dsdt
              dsd = eos_state % dsdt

              A = p_want - pres_zone
              B = entropy_want - entropy

              dtemp = ((dsd/(dpd-0.5*dr*g_zone))*A - B)/ &
                       (dsd*dpt/(dpd -0.5*dr*g_zone) - dst)

              drho = (A - dpt*dtemp)/(dpd - 0.5*dr*g_zone)

              dens_zone = max(0.9e0_rt*dens_zone, &
                   min(dens_zone + drho, 1.1e0_rt*dens_zone))

              temp_zone = max(0.9e0_rt*temp_zone, &
                   min(temp_zone + dtemp, 1.1e0_rt*temp_zone))


              ! check if the density falls below our minimum cut-off --
              ! if so, floor it
              if (dens_zone < low_density_cutoff) then

                 dens_zone = low_density_cutoff
                 temp_zone = temp_fluff
                 converged_hse = .TRUE.
                 fluff = .TRUE.
                 exit

              endif

              ! if (A < TOL .and. B < ETOL) then
              if (abs(drho) < TOL*dens_zone .and. abs(dtemp) < TOL*temp_zone) then
                 converged_hse = .TRUE.
                 exit
              endif

           else

              ! do isothermal
              p_want = pres(i-1) + &
                   dr*0.5*(dens_zone + model(1,i-1))*g_zone

              temp_zone = temp_fluff

              ! (t, rho) -> (p)
              eos_state % T   = temp_zone
              eos_state % rho = dens_zone
              eos_state % xn  = xn(:)

              call eos_on_host(eos_input_rt, eos_state)

              entropy = eos_state % s
              pres_zone = eos_state % p

              dpd = eos_state % dpdr

              drho = (p_want - pres_zone)/(dpd - 0.5*dr*g_zone)

              dens_zone = max(0.9*dens_zone, &
                   min(dens_zone + drho, 1.1*dens_zone))

              if (abs(drho) < TOL*dens_zone) then
                 converged_hse = .TRUE.
                 exit
              endif

              if (dens_zone < low_density_cutoff) then

                 dens_zone = low_density_cutoff
                 temp_zone = temp_fluff
                 converged_hse = .TRUE.
                 fluff = .TRUE.
                 exit

              endif

           endif

        enddo

        if (.NOT. converged_hse) then

           print *, 'Error zone', i, ' did not converge in init_1d'
           print *, 'integrate up'
           print *, dens_zone, temp_zone
           print *, p_want, entropy_want, entropy
           print *, drho, dtemp
           call castro_error('Error: HSE non-convergence')

        endif

        if (temp_zone < temp_fluff) then
           temp_zone = temp_fluff
           isentropic = .false.
        endif

     else
        dens_zone = low_density_cutoff
        temp_zone = temp_fluff
     endif

     ! call the EOS one more time for this zone and then go on to the next
     ! (t, rho) -> (p)
     eos_state % T   = temp_zone
     eos_state % rho = dens_zone
     eos_state % xn  = xn(:)

     call eos_on_host(eos_input_rt, eos_state)

     pres_zone = eos_state % p

     ! update the thermodynamics in this zone
     model(1,i) = dens_zone
     model(2,i) = temp_zone

     pres(i) = pres_zone

     ! to make this process converge faster, set the density in the next zone to
     ! the density in this zone
     ! model_hse(i+1,idens) = dens_zone

  end do

end subroutine ca_initdata_makemodel

subroutine ca_initdata_overwrite(lo,hi, &
                                 state,state_l1,state_l2,state_h1,state_h2, &
                                 model,model_size,MAESTRO_npts_model,dx,dr, &
                                 xlo,xhi,r_model_start)

  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP
  use network, only: nspec
  use eos_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer lo(2), hi(2), r_model_start
  integer model_size,MAESTRO_npts_model
  integer state_l1,state_l2,state_h1,state_h2
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,NVAR)
  real(rt)         model(model_size,0:MAESTRO_npts_model-1)
  real(rt)         dx(2), xlo(2), xhi(2), dr

  ! local
  integer i,j,n
  real(rt)         ekin,radius,temppres

  type (eos_t) eos_state

  do j=lo(2),hi(2)

     if (j .ge. r_model_start) then

        do i=lo(1),hi(1)
           ! need to zero momentum since we're going to be
           ! lowering the density by many orders of magnitude
           state(i,j,UMX:UMY) = 0.e0_rt

           state(i,j,URHO) = model(1,j)
           state(i,j,UTEMP) = model(2,j)

           state(i,j,UFS  ) = 0.3e0_rt
           state(i,j,UFS+1) = 0.7e0_rt
           state(i,j,UFS+2) = 0.0e0_rt

           ! compute e
           eos_state % rho = state(i,j,URHO)
           eos_state % T   = state(i,j,UTEMP)
           eos_state % xn  = state(i,j,UFS:UFS+nspec-1)

           call eos_on_host(eos_input_rt, eos_state)

           state(i,j,UEINT) = eos_state % e

           ! compute rho*e
           state(i,j,UEINT) = state(i,j,URHO) * state(i,j,UEINT)

           ! compute rho*E = rho*e
           state(i,j,UEDEN) = state(i,j,UEINT)

           ! convert X to rhoX
           do n = 1,nspec
              state(i,j,UFS+n-1) = state(i,j,URHO) * state(i,j,UFS+n-1)
           end do
        end do

     end if

  end do

end subroutine ca_initdata_overwrite
