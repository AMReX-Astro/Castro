
      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
      use probdata_module
      use network, only : network_init
      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer init, namlen
      integer name(namlen)
      real(rt)         problo(2), probhi(2)

      integer untin,i,j,k,dir
      real(rt)         rn

      namelist /fortin/ rhocv, T0, Eexp, rexp

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
      integer maxlen
      parameter (maxlen=256)
      character probin*(maxlen)

      call network_init()

      if (namlen .gt. maxlen) then
         write(6,*) 'probin file name too long'
         stop
      end if

      do i = 1, namlen
         probin(i:i) = char(name(i))
      end do

!     Set defaults
      rhocv = -1.0e50_rt
      T0 = -1.0e50_rt
      Eexp = -1.e-50_rt
      rexp = -1.e50_rt

      ! domain extrema
      xmin = problo(1)
      xmax = probhi(1)

      ymin = problo(2)
      ymax = probhi(2)

!     Read namelists
      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      read(untin,fortin)
      close(unit=untin)

      end

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
      use eos_module, only : eos
      use eos_type_module, only : eos_t, eos_input_rt
      use meth_params_module, only : URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFS, UFX

      use network, only : nspec

      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer level, nscal
      integer lo(2), hi(2)
      integer state_l1,state_l2,state_h1,state_h2
      real(rt)         xlo(2), xhi(2), time, delta(2)
      real(rt)         state(state_l1:state_h1,state_l2:state_h2, &
          nscal)

      integer i,j, ii,jj
      type(eos_t) :: eos_state
      real(rt)         :: rho, cv, T, p, eint
      real(rt)         :: rhoeexp, Vexp, rhoe0
      real(rt)         :: xx, xcl, xcr, dx_sub
      real(rt)         :: yy, ycl, ycr, dy_sub 
      real(rt)         :: rr2, xcmin, xcmax, ycmin, ycmax, rcmin, rcmax
      real(rt)         :: vol_pert, vol_ambient, T_zone, rhoe_zone
      integer, parameter :: nsub = 64

      Vexp = 4.e0_rt/3.e0_rt*Pi * rexp**3

      rhoeexp = Eexp / Vexp

      rho = 1.0e0_rt
      T = 1.0e0_rt

      eos_state % rho = rho
      eos_state % T   = T
      eos_state % xn  = 0.e0_rt
      eos_state % xn(1)  = 1.e0_rt

      call eos(eos_input_rt, eos_state)

      cv = eos_state % cv

      rho = rhocv / cv

      eint = rhoeexp / rho
      T = eint / cv

      rhoe0 = rho * cv * T0

      dx_sub = delta(1)/dble(nsub)
      dy_sub = delta(2)/dble(nsub)

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

            rcmin = sqrt(xcmin**2+ycmin**2)
            rcmax = sqrt(xcmax**2+ycmax**2)

            state(i,j,URHO) = rho
            state(i,j,UMX)  = 0.e0_rt
            state(i,j,UMY)  = 0.e0_rt

            if (rcmin .ge. rexp) then
               T_zone = T0
               rhoe_zone = rhoe0
            else if (rcmax .le. rexp) then
               T_zone = T
               rhoe_zone = rhoeexp
            else 

               vol_pert    = 0.e0_rt
               vol_ambient = 0.e0_rt

               do jj = 0, nsub-1
                  yy = ycl + (dble(jj) + 0.5e0_rt) * dy_sub
                  do ii = 0, nsub-1
                     xx = xcl + (dble(ii) + 0.5e0_rt) * dx_sub
                           
                     rr2 = xx**2 + yy**2
                           
                     if (rr2 <= rexp**2) then
                        vol_pert = vol_pert + xx
                     else
                        vol_ambient = vol_ambient + xx
                     endif
                  end do
               end do
               
               T_zone = (vol_pert*T + vol_ambient*T0)/(vol_pert+vol_ambient)
               rhoe_zone = (vol_pert*rhoeexp + vol_ambient*rhoe0)/(vol_pert+vol_ambient)

            end if

            state(i,j,UTEMP) = T_zone
            state(i,j,UEDEN) = rhoe_zone
            state(i,j,UEINT) = rhoe_zone

            state(i,j,UFS) = state(i,j,URHO)

         enddo
      enddo

      end subroutine ca_initdata


! ::: 
! ::: -----------------------------------------------------------
! :::
      subroutine ca_initrad(level,time,lo,hi,nrad, &
           rad_state,rad_state_l1,rad_state_l2, &
           rad_state_h1,rad_state_h2, &
           delta,xlo,xhi)

        use probdata_module

        use amrex_fort_module, only : rt => amrex_real
        implicit none
        integer level, nrad
        integer lo(2), hi(2)
        integer rad_state_l1,rad_state_l2
        integer rad_state_h1,rad_state_h2
        real(kind=8) xlo(2), xhi(2), time, delta(2)
        real(kind=8) rad_state(rad_state_l1:rad_state_h1, &
             rad_state_l2:rad_state_h2, nrad)

        integer i,j

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              rad_state(i,j,:) = 0.e0_rt
           end do
        end do
      end subroutine ca_initrad

