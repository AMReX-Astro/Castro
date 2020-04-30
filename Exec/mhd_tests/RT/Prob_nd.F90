subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use probdata_module, only: frac, rho_1, rho_2, p0_base, split, L_x
  use castro_error_module, only: castro_error
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  call probdata_init(name, namlen)

  ! set local variable defaults
  split(:) = frac * (problo(:) + probhi(:))

  L_x = probhi(1) - problo(1)

end subroutine amrex_probinit


subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       dx, xlo, xhi)

  use probdata_module, only: p0_base, split, L_x, rho_1, rho_2
  use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, &
                                UEDEN, UEINT, UFS, UTEMP, small_temp
  use amrex_constants_module, only: ZERO, HALF, M_PI
  use actual_eos_module, only: gamma_const
  use amrex_fort_module, only: rt => amrex_real
  use prob_params_module, only: dim, problo

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, dx(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  integer  :: i, j, k
  real(rt) :: r(3), r2d, pres, presmid, pertheight

  presmid  = p0_base - rho_1 * split(dim)

  do k = lo(3), hi(3)
     r(3) = problo(3) + (dble(k) + HALF) * dx(3)

     do j = lo(2), hi(2)
        r(2) = problo(2) + (dble(j) + HALF) * dx(2)

        do i = lo(1), hi(1)
           r(1) = problo(1) + (dble(i) + HALF) * dx(1)

           if (r(dim) .lt. split(dim)) then
              pres = p0_base - rho_1 * r(dim)
              state(i,j,k,UEDEN) = pres / (gamma_const - 1.0e0_rt)
              state(i,j,k,UEINT) = pres / (gamma_const - 1.0e0_rt)
           else
              pres = presmid - rho_2 * (r(dim) - split(dim))
              state(i,j,k,UEDEN) = pres / (gamma_const - 1.0e0_rt)
              state(i,j,k,UEINT) = pres / (gamma_const - 1.0e0_rt)
           end if
        
           !doing it similar to 2d, will be something in x-z though 

           ! We explicitly make the perturbation symmetric here.
           ! This prevents the RT from bending.
           pertheight = 0.01e0_rt * HALF * (cos(2.0e0_rt * M_PI * r(1) / L_x) + &
                                            cos(2.0e0_rt * M_PI * (L_x - r(1)) / L_x)) + 0.5e0_rt
 
                         
           state(i,j,k,URHO) = rho_1 + ((rho_2 - rho_1) / 2.0e0_rt) * &
                                        (1 + tanh((r(dim) - pertheight) / 0.005e0_rt))
           state(i,j,k,UFS) = state(i,j,k,URHO)

           state(i,j,k,UMX)   = ZERO
           state(i,j,k,UMY)   = ZERO
           state(i,j,k,UMZ)   = ZERO
           state(i,j,k,UTEMP) = small_temp

        end do
     end do
  end do

end subroutine ca_initdata

! :::
! ::: --------------------------------------------------------------------
! :::
subroutine ca_initmag(level, time, lo, hi, &
                      nbx, mag_x, bx_lo, bx_hi, &
                      nby, mag_y, by_lo, by_hi, &
                      nbz, mag_z, bz_lo, bz_hi, &
                      delta, xlo, xhi)

  use probdata_module
  use prob_params_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module, only : M_PI, TWO, FOUR


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

  real(rt) :: x, y, z
  integer  :: i, j, k

  print *, "Initializing magnetic field!!"

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
                 mag_x(i,j,k,1) = B_0

           enddo
        enddo
     enddo

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)
                 mag_y(i,j,k,1) = 0.0e0_rt
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)+1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
                 mag_z(i,j,k,1) = 0.0e0_rt
           enddo
        enddo
     enddo

  end subroutine ca_initmag


