
      subroutine PROBINIT (init,name,namlen,problo,probhi)
      use probdata_module
      use network, only : network_init
      use meth_params_module, only : outflow_data_old_time, outflow_data_new_time
      use interpolate_module
      use fundamental_constants_module, only : k_B, ev2erg

      implicit none
      integer, intent(in)          :: init, namlen
      integer, intent(in)          :: name(namlen)
      double precision, intent(in) :: problo(2), probhi(2)

      integer untin, i, j, k, kdummy
      double precision dummy_mass
      double precision :: x_distsq, y_distsq, max_dist_for_interp

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
      integer maxlen
      parameter (maxlen=256)
      character probin*(maxlen)

!     Filename for "modelInput" file
      character model_filename*(maxlen)

      namelist /fortin/  &
           model_filename, npts_model, model_has_neut_data, model_temp_in_K, &
           denerr,dengrad,max_denerr_lev,max_dengrad_lev, &
           presserr,pressgrad,max_presserr_lev,max_pressgrad_lev, &
           raderr,radgrad,max_raderr_lev,max_radgrad_lev,&
           enterr, entgrad, max_enterr_lev, max_entgrad_lev, &
           yeerr, yegrad, max_yeerr_lev, max_yegrad_lev, &
           masserr, max_masserr_lev

      if (init .eq. 2) then

!        This option should be called after initialization is complete.
!        Currently this is called from Radiation::post_init

         deallocate(model_rad)
         deallocate(model_rho)
         deallocate(model_vel)
         deallocate(model_pres)
         deallocate(model_temp)
         deallocate(model_energy)
         deallocate(model_entropy)
         deallocate(model_ye)
         if (model_has_neut_data .eq. 1) then
            deallocate(model_neut)
         endif

         return

      endif

      call network_init()

      if (namlen .gt. maxlen) then
         write(6,*) 'probin file name too long'
         stop
      end if

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do

!     Set defaults

      model_filename = "modelInput"
      npts_model = 500
      model_has_neut_data = 0
      model_temp_in_K = 0

      denerr = 1.d20
      dengrad = 1.d20
      max_denerr_lev = -1
      max_dengrad_lev = -1

      presserr = 1.d20
      pressgrad = 1.d20
      max_presserr_lev = -1
      max_pressgrad_lev = -1

      raderr = 1.d20
      radgrad = 1.d20
      max_raderr_lev = -1
      max_radgrad_lev = -1

      enterr = 1.d20
      entgrad = 1.d20
      max_enterr_lev = -1
      max_entgrad_lev = -1

      yeerr = 1.d20
      yegrad = 1.d20
      max_yeerr_lev = -1
      max_yegrad_lev = -1

      masserr = 0.d0
      max_masserr_lev = -1

!     Read namelist
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

!      if (init .eq. 1) then
      if (.true.) then  ! otherwise restart will fail

!        We only need this on first initialization, not on restart.

!        Order of model input file is:
!        radius density velocity pressure temp energy entropy mass ye zone

         allocate(model_rad(npts_model))
         allocate(model_rho(npts_model))
         allocate(model_vel(npts_model))
         allocate(model_pres(npts_model))
         allocate(model_temp(npts_model))
         allocate(model_energy(npts_model))
         allocate(model_entropy(npts_model))
         allocate(model_ye(npts_model))

         i = index(model_filename, ' ')
         open(unit=untin,file=model_filename(1:i-1))

         do k = 1, npts_model
            read(untin,*) model_rad(k),     model_rho(k),  model_vel(k),    &
                          model_pres(k),    model_temp(k), model_energy(k), &
                          model_entropy(k), dummy_mass,    model_ye(k),     &
                          kdummy
         enddo

         if (model_temp_in_K .eq. 1) then
            model_temp = model_temp * k_B / (ev2erg * 1.d6)
         end if

         if (model_has_neut_data .eq. 1) then

            read(untin,*) ngr_model
            allocate(model_neut(npts_model,ngr_model))

            do k = 1, npts_model
               do i = 1, ngr_model - 4, 4
                  read(untin,*) (model_neut(k,j), j = i,i+3)
               enddo
               read(untin,*) (model_neut(k,j), j = i,ngr_model), kdummy
            enddo
         endif

         close(unit=untin)

      endif

      ! Define the center of the star
      center(1) = 0.d0
      center(2) = 0.5d0 * (problo(2) + probhi(2))

      ! lower corner
      corner(1) = problo(1)
      corner(2) = problo(2)

      ! Define the maximum distance from the center of the star to a corner
      x_distsq = max( (problo(1)-center(1))**2 , (probhi(1)-center(1))**2 )
      y_distsq = max( (problo(2)-center(2))**2 , (probhi(2)-center(2))**2 )
      max_dist_for_interp = sqrt( x_distsq + y_distsq)

      ! This fills rho_bndry, etc, with the values at the distance of the far corner
      rho_bndry = interpolate(max_dist_for_interp,npts_model,model_rad,model_rho)
        T_bndry = interpolate(max_dist_for_interp,npts_model,model_rad,model_temp)
       Ye_bndry = interpolate(max_dist_for_interp,npts_model,model_rad,model_ye)
        v_bndry = interpolate(max_dist_for_interp,npts_model,model_rad,model_vel)

      outflow_data_old_time = -1.d0
      outflow_data_new_time = -1.d0

      end subroutine PROBINIT


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
           state,state_l1,state_l2,state_h1,state_h2, &
           delta,xlo,xhi)

      use probdata_module
      use eos_module
      use network, only : nspec, naux
      use meth_params_module, only : URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFS, UFX
      use interpolate_module

      implicit none

      integer level, nscal
      integer lo(2), hi(2)
      integer state_l1, state_l2, state_h1, state_h2
      double precision xlo(2), xhi(2), time, delta(2)
      double precision state(state_l1:state_h1,state_l2:state_h2,nscal)

      integer i,j,n
      double precision :: u, v, dist
      double precision :: x1, x2, y1, y2, x, y
      double precision :: a, c, vel, a_vel, b_vel
      type(eos_t) :: eos_state

      ! Compute formula for interpolating velocity near the origin
      x1 = model_rad(1)
      x2 = model_rad(2)
      y1 = model_vel(1)
      y2 = model_vel(2)
      a_vel = (y1/x1 - y2/x2) / (x1-x2)
      b_vel = (y1 - a_vel*x1**2) / x1

      do j = lo(2), hi(2)
         y = xlo(2) + (dble(j - lo(2)) + 0.5d0) * delta(2) - center(2)

         do i = lo(1), hi(1)   
            x = xlo(1) + (dble(i - lo(1)) + 0.5d0) * delta(1) - center(1)

            dist = sqrt(x**2 + y**2)

            if (dist .lt. x1) then

               ! All variables but velocity -- zero slope at origin
               a = (model_rho(2)-model_rho(1)) / (x2**2 - x1**2)
               c =  model_rho(1) - a * x1**2
               state(i,j,URHO) = a * dist**2 + c

               a = (model_temp(2)-model_temp(1)) / (x2**2 - x1**2)
               c =  model_temp(1) - a * x1**2
               state(i,j,UTEMP) = a * dist**2 + c

               a = (model_ye(2)-model_ye(1)) / (x2**2 - x1**2)
               c =  model_ye(1) - a * x1**2
               state(i,j,UFX) = a * dist**2 + c

               ! Special case for velocity -- value must go to zero at the origin
               vel = a_vel*dist**2 + b_vel*dist
               state(i,j,UMX) = state(i,j,URHO) * vel * x / dist
               state(i,j,UMY) = state(i,j,URHO) * vel * y / dist

            else

               state(i,j,URHO)  = interpolate(dist,npts_model,model_rad,model_rho)
               state(i,j,UTEMP) = interpolate(dist,npts_model,model_rad,model_temp)

               state(i,j,UFX)   = interpolate(dist,npts_model,model_rad,model_ye)

               ! Conert momenutm rather than velocity
               vel = interpolate(dist,npts_model,model_rad,model_rho*model_vel)
               state(i,j,UMX) = vel * (x/dist)
               state(i,j,UMY) = vel * (y/dist)
               ! ignore hydro for now
               ! state(i,j,UMX:UMY)   = 0.d0

            endif

            ! Convert Kelvin to MeV
            !            state(i,UTEMP) = state(i,UTEMP) / 1.160445d10
            state(i,j,UTEMP) = state(i,j,UTEMP)

            ! Set default species concentration to 1.
            state(i,j,UFS) = 1.d0

         enddo
      enddo

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            eos_state % rho = state(i,j,URHO)
            eos_state % T   = state(i,j,UTEMP)
            eos_state % xn  = state(i,j,UFS:UFS+nspec-1)
            eos_state % aux = state(i,j,UFX:UFX+naux-1)

            call eos(eos_input_rt, eos_state)

            u = state(i,j,UMX) / state(i,j,URHO)
            v = state(i,j,UMY) / state(i,j,URHO)

            ! Convert e to rho*e
            state(i,j,UEINT) = state(i,j,URHO) * eos_state % e

            ! Convert e to rho*E = rho*(e + 1/2 (u^2 + v^2) )
            state(i,j,UEDEN) = state(i,j,URHO) * (eos_state % e  &
                 + 0.5d0*(u**2 + v**2)  )

            ! Convert Y to (rho*Y)
            do n = 1,nspec
               state(i,j,UFS+n-1) = state(i,j,URHO) * state(i,j,UFS+n-1)
            end do

            ! Convert Ye to (rho*Ye)
            do n = 1,naux
               state(i,j,UFX+n-1) = state(i,j,URHO) * state(i,j,UFX+n-1)
            end do
               
         enddo
      enddo

      end subroutine ca_initdata

! ::: 
! ::: -----------------------------------------------------------
! :::
      subroutine ca_initrad(level,time,lo,hi,nrad, &
           rad_state,rad_state_l1,rad_state_l2,rad_state_h1,rad_state_h2, &
           delta,xlo,xhi)

        use probdata_module
        use rad_params_module, only: ngroups
        use interpolate_module

        implicit none
        integer level, nrad
        integer lo(2), hi(2)
        integer rad_state_l1,rad_state_h1, rad_state_l2, rad_state_h2

        real(kind=8) xlo(2), xhi(2), time, delta(2)
        real(kind=8) rad_state(rad_state_l1:rad_state_h1,&
                               rad_state_l2:rad_state_h2,nrad)

! begin local variables
        integer i, j, n
        real(kind=8) x1, x2, x, y
        real(kind=8) dist, a, c
!  end  local variables

        if (model_has_neut_data .eq. 0) then

           rad_state(:,:,:) = 0.d0

        else

           if (nrad .ne. ngroups .or. ngroups .ne. ngr_model) then
              print *, "ca_initrad: numbers of groups do not match"
              stop
           endif

           ! For interpolating near the origin
           x1 = model_rad(1)
           x2 = model_rad(2)

           do n = 1, nrad

              do j = lo(2), hi(2)
                 y = xlo(2) + (dble(j-lo(2))+0.5d0) * delta(2) - center(2)

                 do i = lo(1), hi(1)
                    x = xlo(1) + (dble(i-lo(1))+0.5d0) * delta(1) - center(1)

                    dist = sqrt(x**2 + y**2)

                    if (dist .lt. x1) then

                       ! All variables but velocity -- zero slope at origin
                       a = (model_neut(2,n)-model_neut(1,n)) / (x2**2 - x1**2)
                       c =  model_neut(1,n) - a * x1**2
                       rad_state(i,j,n) = a * dist**2 + c

                    else

                       rad_state(i,j,n) = &
                         interpolate(dist,npts_model,model_rad,model_neut(:,n))

                    endif

                 enddo
              enddo
           enddo

        endif

      end subroutine ca_initrad

! ::: 
! ::: -----------------------------------------------------------
! :::

