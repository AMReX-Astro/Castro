   subroutine PROBINIT (init,name,namlen,problo,probhi)

     use probdata_module
     use bl_constants_module
     use fundamental_constants_module
     use eos_module

     implicit none

     integer :: init, namlen
     integer :: name(namlen)
     double precision :: problo(3), probhi(3)

     integer :: untin
     integer :: i

     namelist /fortin/ density, diameter, ambient_dens, problem

     integer, parameter :: maxlen=127
     character :: probin*(maxlen)
     character :: model*(maxlen)
     integer :: ipp, ierr, ipp1

     ! Build "probin" filename -- the name of file containing fortin namelist.
     if (namlen .gt. maxlen) then
        call bl_error("ERROR: probin file name too long")
     end if

     do i = 1, namlen
        probin(i:i) = char(name(i))
     end do

     density  = 1.0d0
     diameter = 1.0d0

     ambient_dens = 1.0d-8

     problem = 1

     ! Read namelists -- override the defaults

     untin = 9
     open(untin,file=probin(1:namlen),form='formatted',status='old')
     read(untin,fortin)
     close(unit=untin)

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
                          state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                          delta,xlo,xhi)

     use probdata_module
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
          UEDEN, UEINT, UFS, UFA
     use network, only : nspec
     use bl_constants_module
     use prob_params_module, only: problo, probhi, center

     implicit none

     integer :: level, nscal
     integer :: lo(3), hi(3)
     integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     double precision :: xlo(3), xhi(3), time, delta(3)
     double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

     double precision :: xx, yy, zz

     integer :: i, j, k

     !$OMP PARALLEL DO PRIVATE(i, j, k, xx, yy, zz)
     do k = lo(3), hi(3)
        zz = xlo(3) + delta(3) * (dble(k-lo(3))+HALF) - center(3)

        do j = lo(2), hi(2)
           yy = xlo(2) + delta(2) * (dble(j-lo(2))+HALF) - center(2)

           do i = lo(1), hi(1)
              xx = xlo(1) + delta(1) * (dble(i-lo(1))+HALF) - center(1)

              ! Establish the cube or sphere

              if (problem .eq. 1 .or. problem .eq. 2) then

                 if ((xx**2 + yy**2 + zz**2)**0.5 < diameter / 2) then
                    state(i,j,k,URHO) = density
                 else
                    state(i,j,k,URHO) = ambient_dens
                 endif

              else if (problem .eq. 3) then

                 if (abs(xx) < diameter/2 .and. abs(yy) < diameter/2 .and. abs(zz) < diameter/2) then
                    state(i,j,k,URHO) = density
                 else
                    state(i,j,k,URHO) = ambient_dens
                 endif

              else

                 call bl_error("Problem not defined.")

              endif

              ! Establish the thermodynamic quantities. They don't have to be
              ! valid because this test will never do a hydro step.

              state(i,j,k,UTEMP) = 1.0d0
              state(i,j,k,UEINT) = 1.0d0
              state(i,j,k,UEDEN) = 1.0d0

              state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,URHO) / nspec

           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO

   end subroutine ca_initdata
