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

    real(rt), intent(inout) :: density_out

    density_out = density

  end subroutine get_density



  ! Set the density.

  subroutine set_density(density_in) bind(C,name='set_density')

    use probdata_module, only: density
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in), value :: density_in

    density = density_in

  end subroutine set_density



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
                end if

             end do
          end do
       end do

    end if

  end subroutine update_density




  ! Compute the norm of the difference between the calculate potential
  ! and the analytical solution.

  subroutine compute_norm(lo, hi, &
                          phi, p_lo, p_hi, &
                          vol, v_lo, v_hi, &
                          dx, norm_power, &
                          norm_diff, norm_exact) bind(C, name='compute_norm')

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO, HALF, FOUR3RD, M_PI
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif
    use fundamental_constants_module, only: Gconst
    use prob_params_module, only: problo, center
    use probdata_module, only: problem, diameter, density
    use reduction_module, only: reduce_add

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: p_lo(3), p_hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3)
    real(rt), intent(in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    real(rt), intent(in   ) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: norm_power
    real(rt), intent(inout) :: norm_diff, norm_exact

    integer  :: i, j, k
    integer  :: ii, jj, kk
    real(rt) :: xx, yy, zz, rr
    real(rt) :: x(0:1), y(0:1), z(0:1), r

    real(rt) :: mass, radius, phiExact

    radius = HALF * diameter
    mass = FOUR3RD * M_PI * radius**3 * density

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             xx = problo(1) + dx(1) * (dble(i) + HALF) - center(1)
             yy = problo(2) + dx(2) * (dble(j) + HALF) - center(2)
             zz = problo(3) + dx(3) * (dble(k) + HALF) - center(3)

             rr = (xx**2 + yy**2 + zz**2)**HALF

             if (problem .eq. 1 .or. problem .eq. 2) then

                if (rr <= radius) then
                   phiExact = -Gconst * mass * (3 * radius**2 - rr**2) / (2 * radius**3)
                else
                   phiExact = -Gconst * mass / rr
                end if

             else if (problem .eq. 3) then

                x(0) = radius + xx
                x(1) = radius - xx
                y(0) = radius + yy
                y(1) = radius - yy
                z(0) = radius + zz
                z(1) = radius - zz

                phiExact = ZERO

                do ii = 0, 1
                   do jj = 0, 1
                      do kk = 0, 1

                         r = (x(ii)**2 + y(jj)**2 + z(kk)**2)**HALF

                         ! This is Equation 20 in Katz et al. (2016). There is a special case
                         ! where, e.g., x(ii) = y(jj) = 0, in which case z(kk) = r and the
                         ! atanh will evaluate to infinity, making the product ill-defined.
                         ! Handle this case by only doing the atanh evaluation away from those
                         ! points, since the product should be zero anyway. We also want to
                         ! avoid a divide by zero, so we'll skip all the cases where r = 0.

                         if (r / radius .gt. 1.e-6_rt) then

                            if (abs(x(ii)) / radius .gt. 1.e-6_rt .and. abs(y(jj)) / radius .gt. 1.e-6_rt) then
                               phiExact = phiExact - Gconst * density * (x(ii) * y(jj) * atanh(z(kk) / r))
                            end if

                            if (abs(y(jj)) / radius .gt. 1.e-6_rt .and. abs(z(kk)) / radius .gt. 1.e-6_rt) then
                               phiExact = phiExact - Gconst * density * (y(jj) * z(kk) * atanh(x(ii) / r))
                            end if

                            if (abs(z(kk)) / radius .gt. 1.e-6_rt .and. abs(x(ii)) / radius .gt. 1.e-6_rt) then
                               phiExact = phiExact - Gconst * density * (z(kk) * x(ii) * atanh(y(jj) / r))
                            end if

                            ! Also, avoid a divide-by-zero for the atan terms.

                            if (abs(x(ii)) / radius .gt. 1.e-6_rt) then
                               phiExact = phiExact + Gconst * density * (x(ii)**2 / 2.0_rt * atan(y(jj) * z(kk) / (x(ii) * r)))
                            end if

                            if (abs(y(jj)) / radius .gt. 1.e-6_rt) then
                               phiExact = phiExact + Gconst * density * (y(jj)**2 / 2.0_rt * atan(z(kk) * x(ii) / (y(jj) * r)))
                            end if

                            if (abs(z(kk)) / radius .gt. 1.e-6_rt) then
                               phiExact = phiExact + Gconst * density * (z(kk)**2 / 2.0_rt * atan(x(ii) * y(jj) / (z(kk) * r)))
                            end if

                         end if

                      end do
                   end do
                end do


#ifndef AMREX_USE_GPU
             else

                call castro_error("Problem not defined.")
#endif

             end if

             call reduce_add(norm_diff, vol(i,j,k) * (phi(i,j,k) - phiExact)**norm_power)
             call reduce_add(norm_exact, vol(i,j,k) * phiExact**norm_power)

          end do
       end do
    end do

  end subroutine compute_norm

end module problem_module
