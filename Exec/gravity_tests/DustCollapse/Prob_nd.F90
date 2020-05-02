subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: ZERO, ONE
  use castro_error_module, only: castro_error
  use probdata_module, only: rho_0, T_0, X_0, p_0, rho_ambient, T_ambient, &
                             r_old, r_old_s, r_0, smooth_delta, r_offset, offset_smooth_delta, &
                             center_x, center_y, center_z, nsub
  use eos_type_module, only: eos_t, eos_input_rp
  use eos_module, only: eos
  use network, only: nspec
  use meth_params_module, only: small_temp
  use prob_params_module, only: center

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  type (eos_t) :: eos_state

  call probdata_init(name, namlen)

  r_old_s = r_old

  ! in 3-d we center the sphere at (center_x, center_y, center_z)

  ! in 2-d we are going to enforce that the lower left corner of the
  ! domain is 0.0 (i.e., we only model a quadrant)

  ! in 1-d, we enforce that the center is the origin (since we are
  ! spherical)

  center(:) = ZERO

#if AMREX_SPACEDIM == 1
  center(1) = ZERO
#else
  center(1) = center_x
  center(2) = center_y
#if AMREX_SPACEDIM == 3
  center(3) = center_z
#endif
#endif

  if (problo(1) /= ZERO) call castro_error("ERROR: xmin should be 0!")
  if (problo(2) /= ZERO) call castro_error("ERROR: ymin should be 0!")
  if (problo(3) /= ZERO) call castro_error("ERROR: zmin should be 0!")

  ! set the composition to be uniform
  X_0(:) = ZERO
  X_0(1) = ONE

  ! get the ambient temperature and sphere temperature, T_0

  eos_state % rho = rho_0
  eos_state % p   = p_0
  eos_state % xn  = x_0
  eos_state % T   = small_temp ! Initial guess for the EOS

  call eos(eos_input_rp, eos_state)

  T_0 = eos_state % T

  eos_state % rho = rho_ambient

  call eos(eos_input_rp, eos_state)

  T_ambient = eos_state % T

end subroutine amrex_probinit



module initdata_module

  implicit none

contains

  subroutine ca_initdata(lo, hi, &
                         state, state_lo, state_hi, &
                         dx, problo) bind(C, name='ca_initdata')

    use amrex_constants_module, only: ZERO, HALF, ONE
    use probdata_module, only: rho_0, X_0, p_0, rho_ambient, r_0, smooth_delta, r_offset, offset_smooth_delta, nsub
    use eos_type_module, only: eos_t, eos_input_rp
    use eos_module, only: eos
    use network, only: nspec
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS, small_temp
    use prob_params_module, only: center, dg
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: state_lo(3), state_hi(3)
    real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    real(rt), intent(in   ) :: dx(3), problo(3)

    real(rt) :: xl, yl, zl, xx, yy, zz
    real(rt) :: dist
    real(rt) :: pres, eint, temp, avg_rho, rho_n
    real(rt) :: volinv
    real(rt) :: dx_sub, dy_sub, dz_sub
    integer  :: i, j, k, ii, jj, kk, n

    type (eos_t) :: eos_state

    !$gpu

#if AMREX_SPACEDIM == 1
    volinv = ONE / dble(nsub)
#elif AMREX_SPACEDIM == 2
    volinv = ONE / dble(nsub*nsub)
#else
    volinv = ONE / dble(nsub*nsub*nsub)
#endif

    dx_sub = dx(1) / dble(nsub)
    dy_sub = dx(2) / dble(nsub)
    dz_sub = dx(3) / dble(nsub)

    do k = lo(3), hi(3)
       zl = problo(1) + dble(k) * dx(3)

       do j = lo(2), hi(2)
          yl = problo(2) + dble(j) * dx(2)

          do i = lo(1), hi(1)
             xl = problo(3) + dble(i) * dx(1)

             avg_rho = ZERO

             do kk = 0, dg(3) * (nsub-1)
                zz = zl + (dble(kk) + HALF) * dz_sub

                do jj = 0, dg(2) * (nsub-1)
                   yy = yl + (dble(jj) + HALF) * dy_sub

                   do ii = 0, nsub-1
                      xx = xl + (dble(ii) + HALF) * dx_sub

                      dist = sqrt((xx-center(1))**2 + (yy-center(2))**2 + (zz-center(3))**2)

                      ! use a tanh profile to smooth the transition between rho_0
                      ! and rho_ambient
                      rho_n = rho_0 - HALF * (rho_0 - rho_ambient) * (ONE + tanh((dist - r_0) / smooth_delta))

                      ! allow for the center to be empty
                      if (r_offset > ZERO) then
                         rho_n = rho_n - HALF * (rho_n - rho_ambient) * (ONE + tanh((r_offset - dist) / offset_smooth_delta))
                      end if

                      avg_rho = avg_rho + rho_n

                   end do
                end do
             end do

             state(i,j,k,URHO) = avg_rho * volinv

             eos_state % rho = state(i,j,k,URHO)
             eos_state % p   = p_0
             eos_state % T   = small_temp ! Initial guess for the EOS
             eos_state % xn  = X_0

             call eos(eos_input_rp, eos_state)

             temp = eos_state % T
             eint = eos_state % e

             state(i,j,k,UTEMP) = temp
             state(i,j,k,UMX) = ZERO
             state(i,j,k,UMY) = ZERO
             state(i,j,k,UMZ) = ZERO
             state(i,j,k,UEDEN) = state(i,j,k,URHO) * eint
             state(i,j,k,UEINT) = state(i,j,k,URHO) * eint
             do n = 1, nspec
                state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * X_0(n)
             end do

          end do
       end do
    end do

  end subroutine ca_initdata

end module initdata_module

#ifdef MHD
subroutine ca_initmag(level, time, lo, hi, &
                      nbx, mag_x, bx_lo, bx_hi, &
                      nby, mag_y, by_lo, by_hi, &
                      nbz, mag_z, bz_lo, bz_hi, &
                      delta, xlo, xhi)

  use probdata_module, only : B_x, B_y, B_z
  use prob_params_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer :: level, nbx, nby, nbz
  integer :: lo(3), hi(3)
  integer :: bx_lo(3), bx_hi(3)
  integer :: by_lo(3), by_hi(3)
  integer :: bz_lo(3), bz_hi(3)
  real(rt) :: xlo(3), xhi(3), time, delta(3)

  real(rt) :: mag_x(bx_lo(1):bx_hi(1), bx_lo(2):bx_hi(2), bx_lo(3):bx_hi(3), nbx)
  real(rt) :: mag_y(by_lo(1):by_hi(1), by_lo(2):by_hi(2), by_lo(3):by_hi(3), nby)
  real(rt) :: mag_z(bz_lo(1):bz_hi(1), bz_lo(2):bz_hi(2), bz_lo(3):bz_hi(3), nbz)

  real(rt) :: xcen, ycen, zcen
  integer  :: i, j, k


  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
              mag_x(i,j,k,1) = B_x
        enddo
     enddo
  enddo
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
           mag_y(i,j,k,1) = B_y
         enddo
     enddo
  enddo

  do k = lo(3), hi(3)+1
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           mag_z(i,j,k,1) = B_z
        enddo
     enddo
  enddo

end subroutine ca_initmag
#endif
