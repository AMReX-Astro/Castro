subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

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
