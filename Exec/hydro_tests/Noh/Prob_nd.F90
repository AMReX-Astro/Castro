   subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

     use amrex_fort_module, only : rt => amrex_real
     implicit none

     integer :: init, namlen
     integer :: name(namlen)
     real(rt)         :: problo(3), probhi(3)

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

     use eos_module, only: eos
     use eos_type_module, only: eos_t, eos_input_rp
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
                                    UEDEN, UEINT, UFS
     use network, only : nspec
     use amrex_constants_module, only: ZERO, HALF, ONE, TWO, M_PI
     use fundamental_constants_module, only: Gconst, M_solar
     use prob_params_module, only: center, dim, problo

     use amrex_fort_module, only : rt => amrex_real
     implicit none

     integer :: level, nscal
     integer :: lo(3), hi(3)
     integer :: state_lo(3), state_hi(3)
     real(rt)         :: xlo(3), xhi(3), time, delta(3)
     real(rt)         :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

     real(rt)         :: loc(3), r, vel(3)

     type (eos_t) :: zone_state

     integer :: i,j,k

     !$OMP PARALLEL DO PRIVATE(i, j, k, loc, r, zone_state)
     do k = lo(3), hi(3)
        loc(3) = problo(3) + delta(3)*(dble(k)+HALF) - center(3)

        do j = lo(2), hi(2)
           loc(2) = problo(2) + delta(2)*(dble(j)+HALF) - center(2)

           do i = lo(1), hi(1)
              loc(1) = problo(1) + delta(1)*(dble(i)+HALF) - center(1)

              ! Uniform density, negligible pressure.

              zone_state % rho = 1.0e0_rt
              zone_state % P   = 1.0e-6_rt
              zone_state % xn(:) = ONE / nspec

              call eos(eos_input_rp, zone_state)

              ! Radial inflow with |v| = 1.

              r = sqrt( sum(loc**2) )
              vel(:) = -loc(:) / r

              state(i,j,k,URHO)  = zone_state % rho
              state(i,j,k,UTEMP) = zone_state % T
              state(i,j,k,UEINT) = zone_state % e * zone_state % rho
              state(i,j,k,UFS:UFS+nspec-1) = zone_state % xn(:) * zone_state % rho

              state(i,j,k,UMX:UMZ) = state(i,j,k,URHO) * vel(:)

              state(i,j,k,UEDEN) = state(i,j,k,UEINT) + sum( state(i,j,k,UMX:UMZ)**2 ) / ( TWO * state(i,j,k,URHO) )

           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO

   end subroutine ca_initdata
