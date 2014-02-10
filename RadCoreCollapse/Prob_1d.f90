
      subroutine PROBINIT (init,name,namlen,problo,probhi)
      use probdata_module
      use network, only : network_init
      use fundamental_constants_module, only : k_B, ev2erg

      implicit none
      integer, intent(in)          :: init, namlen
      integer, intent(in)          :: name(namlen)
      double precision, intent(in) :: problo(1), probhi(1)

      integer untin, i, j, k, kdummy
      double precision dummy_mass

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
      integer maxlen
      parameter (maxlen=256)
      character probin*(maxlen)

!     Filename for "modelInput" file
      character model_filename*(maxlen)

      namelist /fortin/ do_neutrino_test, &
           model_filename, npts_model, model_has_neut_data, model_temp_in_K, &
           denerr,dengrad,max_denerr_lev,max_dengrad_lev, &
           presserr,pressgrad,max_presserr_lev,max_pressgrad_lev, &
           direction,raderr,radgrad,max_raderr_lev,max_radgrad_lev,&
           enterr, entgrad, max_enterr_lev, max_entgrad_lev, &
           yeerr, yegrad, max_yeerr_lev, max_yegrad_lev, &
           masserr, max_masserr_lev, init_smoothing

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

      init_smoothing = 0

!     Read namelist
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

      if (init .eq. 1) then

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
           state,state_l1,state_h1,delta,xlo,xhi)

      use probdata_module
      use eos_module
      use network, only : nspec, naux
      use meth_params_module, only : URHO, UMX, UEDEN, UEINT, UTEMP, UFS, UFX
      use interpolate_module
      use filter_module

      implicit none

      integer level, nscal
      integer lo(1), hi(1)
      integer state_l1,state_h1
      double precision xlo(1), xhi(1), time, delta(1)
      double precision state(state_l1:state_h1,nscal)

      integer i,n
      double precision :: u, pres, dist
      double precision :: x1, x2, y1, y2
      double precision :: a, c, a_vel, b_vel
      double precision :: x_in(nspec+naux)

      double precision, dimension(lo(1)-4:hi(1)+4) :: t_rho, t_temp, t_ye, t_mx
      integer, parameter :: S=3

      ! Compute formula for interpolating velocity near the origin
      x1 = model_rad(1)
      x2 = model_rad(2)
      y1 = model_vel(1)
      y2 = model_vel(2)
      a_vel = (y1/x1 - y2/x2) / (x1-x2)
      b_vel = (y1 - a_vel*x1**2) / x1

      do i = lo(1)-4, hi(1)+4

         dist = (dble(i) + 0.5d0) * delta(1)

         if (dist .lt. x1) then

            ! All variables but velocity -- zero slope at origin
            a = (model_rho(2)-model_rho(1)) / (x2**2 - x1**2)
            c =  model_rho(1) - a * x1**2
            t_rho(i) = a * dist**2 + c

            a = (model_temp(2)-model_temp(1)) / (x2**2 - x1**2)
            c =  model_temp(1) - a * x1**2
            t_temp(i) = a * dist**2 + c

            a = (model_ye(2)-model_ye(1)) / (x2**2 - x1**2)
            c =  model_ye(1) - a * x1**2
            t_ye(i) = a * dist**2 + c

            ! Special case for velocity -- value must go to zero at the origin
            if (dist > 0.d0) then
               t_mx(i) = a_vel*dist**2 + b_vel*dist
            else
               t_mx(i) = -a_vel*dist**2 + b_vel*dist
            end if
            t_mx(i) = t_rho(i) * t_mx(i)

         else

            t_rho(i)  = interpolate(dist,npts_model,model_rad,model_rho)
            t_temp(i) = interpolate(dist,npts_model,model_rad,model_temp)

            t_ye(i) = interpolate(dist,npts_model,model_rad,model_ye)

            ! Conert momenutm rather than velocity
            t_mx(i)   = interpolate(dist,npts_model,model_rad,model_rho*model_vel)

         endif

      enddo

      do i = lo(1), hi(1)

         if (init_smoothing > 0) then
            state(i,URHO)  = ff4(0,S) *  t_rho(i) &
                 &         + ff4(1,S) * (t_rho(i-1)+t_rho(i+1)) &
                 &         + ff4(2,S) * (t_rho(i-2)+t_rho(i+2)) &
                 &         + ff4(3,S) * (t_rho(i-3)+t_rho(i+3)) &
                 &         + ff4(4,S) * (t_rho(i-4)+t_rho(i+4))
            state(i,UTEMP) = ff4(0,S) *  t_temp(i) &
                 &         + ff4(1,S) * (t_temp(i-1)+t_temp(i+1)) &
                 &         + ff4(2,S) * (t_temp(i-2)+t_temp(i+2)) &
                 &         + ff4(3,S) * (t_temp(i-3)+t_temp(i+3)) &
                 &         + ff4(4,S) * (t_temp(i-4)+t_temp(i+4))
            state(i,UFX)   = ff4(0,S) *  t_ye(i) &
                 &         + ff4(1,S) * (t_ye(i-1)+t_ye(i+1)) &
                 &         + ff4(2,S) * (t_ye(i-2)+t_ye(i+2)) &
                 &         + ff4(3,S) * (t_ye(i-3)+t_ye(i+3)) &
                 &         + ff4(4,S) * (t_ye(i-4)+t_ye(i+4))
            state(i,UMX)   = ff4(0,S) *  t_mx(i) &
                 &         + ff4(1,S) * (t_mx(i-1)+t_mx(i+1)) &
                 &         + ff4(2,S) * (t_mx(i-2)+t_mx(i+2)) &
                 &         + ff4(3,S) * (t_mx(i-3)+t_mx(i+3)) &
                 &         + ff4(4,S) * (t_mx(i-4)+t_mx(i+4))
         else
            state(i,URHO) = t_rho(i) 
            state(i,UTEMP) = t_temp(i)
            state(i,UFX) = t_ye(i)
            state(i,UMX) = t_mx(i)
         end if

         ! Set default species concentration to 1.
         state(i,UFS) = 1.d0

         call eos_given_RTX(state(i,UEDEN),pres,state(i,URHO),state(i,UTEMP),state(i,UFS:))

         u = state(i,UMX) / state(i,URHO)

         ! Convert e to rho*e
         state(i,UEINT) = state(i,URHO) * state(i,UEDEN)

         ! Convert e to rho*E = rho*(e + 1/2 u^2)
         state(i,UEDEN) = state(i,URHO) * (state(i,UEDEN) + 0.5d0*(u**2))

         ! Convert Y to (rho*Y)
         do n = 1,nspec
            state(i,UFS+n-1) = state(i,URHO) * state(i,UFS+n-1)
         end do

         ! Convert Ye to (rho*Ye)
         do n = 1,naux
            state(i,UFX+n-1) = state(i,URHO) * state(i,UFX+n-1)
         end do

      enddo

    end subroutine ca_initdata

! ::: 
! ::: -----------------------------------------------------------
! :::
      subroutine ca_initrad(level,time,lo,hi,nrad, &
           rad_state,rad_state_l1,rad_state_h1, &
           delta,xlo,xhi)

        use probdata_module
        use rad_params_module, only: ngroups
        use interpolate_module

        implicit none
        integer level, nrad
        integer lo(1), hi(1)
        integer rad_state_l1,rad_state_h1

        real(kind=8) xlo(1), xhi(1), time, delta(1)
        real(kind=8) rad_state(rad_state_l1:rad_state_h1,nrad)

! begin local variables
        integer i, n
        real(kind=8) x1, x2
        real(kind=8) dist, a, c
!  end  local variables

        if (model_has_neut_data .eq. 0) then

           do i = lo(1), hi(1)
              rad_state(i,:) = 0.d0
           enddo

        else

           if (nrad .ne. ngroups .or. ngroups .ne. ngr_model) then
              print *, "ca_initrad: numbers of groups do not match"
              stop
           endif

           ! For interpolating near the origin
           x1 = model_rad(1)
           x2 = model_rad(2)

           do n = 1, nrad

              do i = lo(1), hi(1)

                 dist = (dble(i) + 0.5d0) * delta(1)

                 if (dist .lt. x1) then

                    ! All variables but velocity -- zero slope at origin
                    a = (model_neut(2,n)-model_neut(1,n)) / (x2**2 - x1**2)
                    c =  model_neut(1,n) - a * x1**2
                    rad_state(i,n) = a * dist**2 + c

                 else

                    rad_state(i,n) = &
                       interpolate(dist,npts_model,model_rad,model_neut(:,n))

                 endif

              enddo

           enddo

        endif

      end subroutine ca_initrad

! ::: 
! ::: -----------------------------------------------------------
! :::
