      
      subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
      use probdata_module
      use network, only : network_init
      use bl_fort_module, only : rt => c_real
      implicit none
      integer init, namlen
      integer name(namlen)
      real(rt)         problo(3), probhi(3)

      integer untin,i,j,k,dir
      real(rt)         rn

      namelist /fortin/ rhocv, T0, Eexp, rexp

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
      integer maxlen
      parameter (maxlen=127)
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

      ! domain extrema and center
      xmin = problo(1)
      xmax = probhi(1)

      ymin = problo(2)
      ymax = probhi(2)

      zmin = problo(3)
      zmax = probhi(3)

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
       state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
       delta,xlo,xhi)

      use probdata_module
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS

      use network, only : nspec
      use eos_module

      use bl_fort_module, only : rt => c_real
      implicit none
      integer level, nscal
      integer lo(3), hi(3)
      integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
      real(rt)         xlo(3), xhi(3), time, delta(3)
      real(rt)         state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3, &
          nscal)

      integer i,j,k, ii, jj, kk
      type(eos_t) :: eos_state
      real(rt)         :: rho, cv, T, p, eint
      real(rt)         :: rhoeexp, Vexp, rhoe0
      real(rt)         :: xx, xcl, xcr, dx_sub
      real(rt)         :: yy, ycl, ycr, dy_sub 
      real(rt)         :: zz, zcl, zcr, dz_sub 
      real(rt)         :: rr2, xcmin, xcmax, ycmin, ycmax, zcmin, zcmax , rcmin, rcmax
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
               state(i,j,k,UMX)  = 0.e0_rt
               state(i,j,k,UMY)  = 0.e0_rt
               state(i,j,k,UMZ)  = 0.e0_rt

               if (rcmin .ge. rexp) then
                  T_zone = T0
                  rhoe_zone = rhoe0
               else if (rcmax .le. rexp) then
                  T_zone = T
                  rhoe_zone = rhoeexp
               else 
                  vol_pert    = 0.e0_rt
                  vol_ambient = 0.e0_rt
                  
                  do kk = 0, nsub-1
                     zz = zcl + (dble(kk) + 0.5e0_rt) * dz_sub
                     do jj = 0, nsub-1
                        yy = ycl + (dble(jj) + 0.5e0_rt) * dy_sub
                        do ii = 0, nsub-1
                           xx = xcl + (dble(ii) + 0.5e0_rt) * dx_sub
                           
                           rr2 = xx**2 + yy**2 + zz**2
                           
                           if (rr2 <= rexp**2) then
                              vol_pert = vol_pert + 1.e0_rt
                           else
                              vol_ambient = vol_ambient + 1.e0_rt
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

!     -----------------------------------------------------------


      subroutine ca_initrad(level,time,lo,hi,nrad, &
           rad_state,rad_state_l1,rad_state_l2,rad_state_l3, &
           rad_state_h1,rad_state_h2,rad_state_h3, &
           delta,xlo,xhi)

        use probdata_module

        use bl_fort_module, only : rt => c_real
        implicit none
        integer level, nrad
        integer lo(1), hi(1)
        integer rad_state_l1,rad_state_l2,rad_state_l3
        integer rad_state_h1,rad_state_h2,rad_state_h3
        real(kind=8) xlo(3), xhi(3), time, delta(3)
        real(kind=8) rad_state(rad_state_l1:rad_state_h1, &
             rad_state_l2:rad_state_h2, & 
             rad_state_l3:rad_state_h3,nrad)

        rad_state = 0.e0_rt

      end subroutine ca_initrad


