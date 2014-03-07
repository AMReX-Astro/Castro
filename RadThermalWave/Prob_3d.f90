      
      subroutine PROBINIT (init,name,namlen,problo,probhi)
      use probdata_module
      use network, only : network_init
      implicit none
      integer init, namlen
      integer name(namlen)
      double precision problo(3), probhi(3)

      integer untin,i,j,k,dir
      double precision rn

      namelist /fortin/ rhocv, T0, Eexp, rexp, &
           denerr,   dengrad,  max_denerr_lev,  max_dengrad_lev, &
           velerr,   velgrad,  max_velerr_lev,  max_pressgrad_lev, &
           presserr, pressgrad,max_presserr_lev,max_pressgrad_lev, &
           temperr,  tempgrad, max_temperr_lev, max_tempgrad_lev, &
           raderr,   radgrad,  max_raderr_lev,  max_radgrad_lev
      
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
      rhocv = -1.0d50
      T0 = -1.0d50
      Eexp = -1.d-50
      rexp = -1.d50

      denerr = 1.d20
      dengrad = 1.d20
      max_denerr_lev = -1
      max_dengrad_lev = -1
      
      velerr = 1.d20
      velgrad = 1.d20
      max_velerr_lev = -1
      max_velgrad_lev = -1
      
      presserr = 1.d20
      pressgrad = 1.d20
      max_presserr_lev = -1
      max_pressgrad_lev = -1
      
      temperr = -1.d20
      tempgrad = 1.d20
      max_temperr_lev = -1
      max_tempgrad_lev = -1
      
      raderr = 1.d20
      radgrad = 1.d20
      max_raderr_lev = -1
      max_radgrad_lev = -1
      
      ! domain extrema and center
      xmin = problo(1)
      xmax = probhi(1)
      center(1) = 0.5d0*(xmin + xmax)

      ymin = problo(2)
      ymax = probhi(2)
      center(2) = 0.5d0*(ymin + ymax)

      zmin = problo(3)
      zmax = probhi(3)
      center(3) = 0.5d0*(zmin + zmax)

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

      implicit none
      integer level, nscal
      integer lo(3), hi(3)
      integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
      double precision xlo(3), xhi(3), time, delta(3)
      double precision state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3, &
          nscal)

      integer i,j,k, ii, jj, kk
      type(eos_t) :: eos_state
      double precision :: rho, cv, T, p, eint
      double precision :: rhoeexp, Vexp, rhoe0
      double precision :: xx, xcl, xcr, dx_sub
      double precision :: yy, ycl, ycr, dy_sub 
      double precision :: zz, zcl, zcr, dz_sub 
      double precision :: rr2, xcmin, xcmax, ycmin, ycmax, zcmin, zcmax , rcmin, rcmax
      double precision :: vol_pert, vol_ambient, T_zone, rhoe_zone
      integer, parameter :: nsub = 64

      Vexp = 4.d0/3.d0*Pi * rexp**3

      rhoeexp = Eexp / Vexp

      rho = 1.0d0
      T = 1.0d0

      eos_state % rho = rho
      eos_state % T   = T
      eos_state % xn  = 0.d0
      eos_state % xn(1)  = 1.d0
      
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
               state(i,j,k,UMX)  = 0.d0
               state(i,j,k,UMY)  = 0.d0
               state(i,j,k,UMZ)  = 0.d0

               if (rcmin .ge. rexp) then
                  T_zone = T0
                  rhoe_zone = rhoe0
               else if (rcmax .le. rexp) then
                  T_zone = T
                  rhoe_zone = rhoeexp
               else 
                  vol_pert    = 0.d0
                  vol_ambient = 0.d0
                  
                  do kk = 0, nsub-1
                     zz = zcl + (dble(kk) + 0.5d0) * dz_sub
                     do jj = 0, nsub-1
                        yy = ycl + (dble(jj) + 0.5d0) * dy_sub
                        do ii = 0, nsub-1
                           xx = xcl + (dble(ii) + 0.5d0) * dx_sub
                           
                           rr2 = xx**2 + yy**2 + zz**2
                           
                           if (rr2 <= rexp**2) then
                              vol_pert = vol_pert + 1.d0
                           else
                              vol_ambient = vol_ambient + 1.d0
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

        implicit none
        integer level, nrad
        integer lo(1), hi(1)
        integer rad_state_l1,rad_state_l2,rad_state_l3
        integer rad_state_h1,rad_state_h2,rad_state_h3
        real(kind=8) xlo(3), xhi(3), time, delta(3)
        real(kind=8) rad_state(rad_state_l1:rad_state_h1, &
             rad_state_l2:rad_state_h2, & 
             rad_state_l3:rad_state_h3,nrad)

        rad_state = 0.d0

      end subroutine ca_initrad


!     -----------------------------------------------------------


      subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
          adv_h3,domlo,domhi,delta,xlo,time,bc)

      use meth_params_module, only : NVAR

      implicit none

      include 'bc_types.fi'
      integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      double precision delta(3), xlo(3), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3, &
                           NVAR)

      integer n

      do n = 1,NVAR
         call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
             adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
             domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

      end

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
          adv_h3,domlo,domhi,delta,xlo,time,bc)
      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      double precision delta(3), xlo(3), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

      call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
          domlo,domhi,delta,xlo,bc)

      end

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_radfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
          adv_h3,domlo,domhi,delta,xlo,time,bc)
      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      double precision delta(3), xlo(3), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

      call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
          domlo,domhi,delta,xlo,bc)

      end

! ::: -----------------------------------------------------------
