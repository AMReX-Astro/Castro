subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use probdata_module
  use model_parser_module
  use amrex_error_module
  use amrex_constants_module
  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module, only : const_grav

  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt) :: problo(3), probhi(3)

  real(rt) :: offset
  integer untin,i

  namelist /fortin/ model_name, apply_vel_field,shear_vel_field, &
       velpert_scale, velpert_amplitude, velpert_height_loc, num_vortices, &
       shear_height_loc, shear_amplitude, shear_height, shear_width_x,&
       shear_width_y, cutoff_density

  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  real(rt) :: max_hse_err, dpdr, rhog, hse_err, dr_model

  ! Build "probin" filename from C++ land --
  ! the name of file containing fortin namelist.


  if (namlen .gt. maxlen) call amrex_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do


  ! Namelist defaults
  apply_vel_field = .false.
  shear_vel_field = .false.
  velpert_scale = 1.0e2_rt
  velpert_amplitude = 1.0e2_rt
  velpert_height_loc = 6.5e3_rt
  num_vortices = 1
  cutoff_density = 50.e0_rt
  shear_height_loc = 0.0_rt
  shear_height     = 0.0_rt
  shear_width_x    = 0.0_rt
  shear_width_y    = 0.0_rt
  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)


  ! Read initial model
  call read_model_file(model_name)

  ! HSE check
  max_hse_err = -1.e30
  dr_model = model_r(2) - model_r(1)
  do i = 2, npts_model-2
     dpdr = (model_state(i,ipres_model) - model_state(i-1,ipres_model))/dr_model
     rhog = HALF*(model_state(i,idens_model) + model_state(i-1,idens_model))*const_grav
     hse_err = abs(dpdr - rhog)

     ! only count the error if we have a non-zero gradient in the next cell
     if (model_state(i+1,ipres_model) /= model_state(i,ipres_model)) then
        max_hse_err = max(max_hse_err, hse_err)
     endif
  enddo

  print *, 'maximum hse error = ', max_hse_err

  if (parallel_IOProcessor()) then
     do i = 1, npts_model
        print *, i, model_r(i), model_state(i,idens_model)
     enddo
  endif

  ! set local variable defaults
  center(1) = HALF*(problo(1)+probhi(1))
  center(2) = HALF*(problo(2)+probhi(2))
  center(3) = HALF*(problo(3)+probhi(3))

  ! velocity perturbation stuff
  offset = (probhi(1) - problo(1)) / (num_vortices)

  allocate(xloc_vortices(num_vortices))

  do i = 1, num_vortices
     xloc_vortices(i) = (dble(i-1) + HALF) * offset + problo(1)
  enddo

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
     state,state_l1,state_l2,state_l3,state_h1,state_h2, &
     state_h3,delta,xlo,xhi)

  use amrex_constants_module
  use probdata_module
  use interpolate_module
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY,UMZ, UEDEN, UEINT, UFS, UTEMP
  use network, only: nspec
  use model_parser_module
  use eos_type_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer level, nscal
  integer lo(3), hi(3)
  integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  real(rt) :: xlo(3), xhi(3), time, delta(3)
  real(rt) :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

  real(rt) :: xdist,ydist,zdist,x,y,z,r,upert(2),shear_height_r,velocity_gradient
  integer i,j,k,n,vortex,int_shear_width_x,int_shear_width_y,shear_bottom_index

  type (eos_t) :: eos_state

  shear_height_r = shear_height * delta(3)

  int_shear_width_x=int(shear_width_x)
  int_shear_width_y=int(shear_width_y)
  velocity_gradient=shear_amplitude/shear_height

  shear_bottom_index=int((shear_height_loc)/delta(3))-shear_height

  do k = lo(3), hi(3)
     z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           state(i,j,k,URHO)  = interpolate(z,npts_model,model_r, &
                model_state(:,idens_model))
           state(i,j,k,UTEMP) = interpolate(z,npts_model,model_r, &
                model_state(:,itemp_model))


           do n = 1, nspec
              state(i,j,k,UFS-1+n) = interpolate(z,npts_model,model_r, &
                   model_state(:,ispec_model-1+n))
           enddo

           eos_state%rho = state(i,j,k,URHO)
           eos_state%T = state(i,j,k,UTEMP)
           eos_state%xn(:) = state(i,j,k,UFS:)

           call eos(eos_input_rt, eos_state)
           state(i,j,k,UEINT) = eos_state%e

        end do
     end do
  end do
  ! switch to conserved quantities
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           state(i,j,k,UEDEN) = state(i,j,k,URHO) * state(i,j,k,UEDEN)
           state(i,j,k,UEINT) = state(i,j,k,URHO) * state(i,j,k,UEINT)

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

        enddo
     end do
  enddo

  ! Initial velocities
  state(:,:,:,UMX:UMZ) = ZERO
  if (shear_vel_field)then
  ! First give shear velocity at h>=shear_height_loc - shear_height/2


    do k= lo(3),hi(3)
      z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)

      if(k > shear_bottom_index .and. k <= shear_bottom_index + shear_height)then
        state(:,:,k,UMX) = state(:,:,k,URHO)*velocity_gradient*abs(k-shear_bottom_index)
      elseif (k> shear_bottom_index + shear_height)then
        state(:,:,k,UMX) = state(:,:,k,URHO) * shear_amplitude
      end if
    end do


    do k= lo(3),hi(3)
      z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)
      if(int_shear_width_y>2)then
        do j = lo(2),hi(2)
          if (k > shear_bottom_index .and. k <= shear_bottom_index + shear_height .and. mod(j,int_shear_width_y)<=4) then
            state(:,j,k,UMX:UMZ) = ZERO
          end if
        end do
      end if
      if(int_shear_width_x>2)then
          do i = lo(1),hi(1)
            if (k > shear_bottom_index .and. k <= shear_bottom_index + shear_height  .and. mod(i,int_shear_width_x)<=4) then
              state(i,:,k,UMX:UMZ) = ZERO
            endif
          end do
        end if
    end do
  end if

  ! Now add the velocity perturbation (update the kinetic energy too)
  if (apply_vel_field) then
     do k = lo(3), hi(3)
        z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)
        zdist = z - velpert_height_loc
        do i = lo(1), hi(1)
           x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

           upert = ZERO

           ! loop over each vortex
           if(velpert_height_loc>ZERO .and. velpert_scale>ZERO .and. velpert_amplitude>ZERO &
             .and. num_vortices>0) then
           do vortex = 1, num_vortices

              xdist = x - xloc_vortices(vortex)
              r = sqrt(xdist**2 + zdist**2)

              upert(1) = upert(1) - (zdist/velpert_scale) * &
                   velpert_amplitude * exp(-r**2/(TWO*velpert_scale**2)) &
                   * (-ONE)**vortex

              upert(2) = upert(2) + (xdist/velpert_scale) * &
                   velpert_amplitude * exp(-r**2/(TWO*velpert_scale**2)) &
                   * (-ONE)**vortex
           enddo
           else
            upert=ZERO
           end if

           do j = lo(2), hi(2)
             state(i,j,k,UMX) = state(i,j,k,UMX) + state(i,j,k,URHO) * upert(1)
             state(i,j,k,UMZ) = state(i,j,k,UMZ) + state(i,j,k,URHO) * upert(2)
             state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + HALF*sum(state(i,j,k,UMX:UMZ)**2)/state(i,j,k,URHO)
           end do
        end do
     end do
  endif



end subroutine ca_initdata


subroutine ca_initrad(level,time,lo,hi,nrad, &
     rad_state,rad_state_l1,rad_state_l2, rad_state_l3, &
     rad_state_h1,rad_state_h2, rad_state_h3, &
     delta,xlo,xhi)

  use amrex_constants_module
  use probdata_module
  use amrex_fort_module, only : rt => amrex_real

  integer level, nrad
  integer lo(3), hi(3)
  integer rad_state_l1,rad_state_l2,rad_state_l3
  integer rad_state_h1,rad_state_h2,rad_state_h3
  real(rt) :: xlo(3), xhi(3), time, delta(3)
  real(rt) :: rad_state(rad_state_l1:rad_state_h1, &
       rad_state_l2:rad_state_h2,rad_state_l3:rad_state_h3, nrad)

  integer i,j,k
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           rad_state(i,j,k,:) = ZERO
        end do
     end do
  end do
end subroutine ca_initrad
