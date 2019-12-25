subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  call probdata_init(name, namlen)

end subroutine amrex_probinit



module initdata_module

  implicit none

contains

  subroutine ca_initdata(lo, hi, &
                         state, s_lo, s_hi, &
                         dx, problo) bind(C, name='ca_initdata')

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: HALF
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif
    use probdata_module, only: density, ambient_dens, diameter, problem
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
                                  UEDEN, UEINT, UFS
    use network, only: nspec
    use prob_params_module, only: center

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: dx(3), problo(3)

    integer  :: i, j, k, n
    real(rt) :: xx, yy, zz

    !$gpu

    !$OMP PARALLEL DO PRIVATE(i, j, k, n, xx, yy, zz)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             xx = problo(1) + dx(1) * (dble(i) + HALF) - center(1)
             yy = problo(2) + dx(2) * (dble(j) + HALF) - center(2)
             zz = problo(3) + dx(3) * (dble(k) + HALF) - center(3)

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

#ifndef AMREX_USE_GPU
             else

                call castro_error("Problem not defined.")
#endif

             endif

             ! Establish the thermodynamic quantities. They don't have to be
             ! valid because this test will never do a hydro step.

             state(i,j,k,UTEMP) = 1.0e0_rt
             state(i,j,k,UEINT) = 1.0e0_rt
             state(i,j,k,UEDEN) = 1.0e0_rt

             do n = 1, nspec
                state(i,j,k,UFS+n-1) = state(i,j,k,URHO) / nspec
             end do

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine ca_initdata

end module initdata_module



module problem_module

  implicit none

contains

  ! Return the problem type.

  subroutine get_problem_number(problem_out) bind(C,name='get_problem_number')

    use probdata_module, only: problem
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer :: problem_out

    problem_out = problem

  end subroutine get_problem_number



  ! Return the diameter.

  subroutine get_diameter(diameter_out) bind(C,name='get_diameter')

    use probdata_module, only: diameter
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt) :: diameter_out

    diameter_out = diameter

  end subroutine get_diameter



  ! Return the density.

  subroutine get_density(density_out) bind(C,name='get_density')

    use probdata_module, only: density
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt) :: density_out

    density_out = density

  end subroutine get_density



  ! Update the density field. This ensures
  ! that the sum of the mass on the domain
  ! is what we intend it to be.

  subroutine update_density(lo, hi, &
                            state, s_lo, s_hi, &
                            dx, update_factor) bind(C, name='update_density')

    use amrex_constants_module, only: HALF
    use network, only: nspec
    use meth_params_module, only: NVAR, URHO, UFS
    use prob_params_module, only: problo, center
    use probdata_module, only: problem, diameter
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: update_factor

    integer  :: i, j, k
    real(rt) :: xx, yy, zz

    !$gpu

    if (problem .eq. 2) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                xx = problo(1) + dx(1) * (dble(i) + HALF) - center(1)
                yy = problo(2) + dx(2) * (dble(j) + HALF) - center(2)
                zz = problo(3) + dx(3) * (dble(k) + HALF) - center(3)

                if ((xx**2 + yy**2 + zz**2)**0.5 < diameter / 2) then

                   state(i,j,k,URHO) = state(i,j,k,URHO) * update_factor
                   state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,UFS:UFS-1+nspec) * update_factor

                endif

             enddo
          enddo
       enddo

    endif

  end subroutine update_density

end module problem_module
