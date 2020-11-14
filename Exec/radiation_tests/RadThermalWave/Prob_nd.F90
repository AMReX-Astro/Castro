subroutine amrex_probinit(init,name,namlen,problo,probhi) bind(c)

  use amrex_fort_module, only: rt => amrex_real
  use probdata_module

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
                         dx, problo) &
                         bind(C, name='ca_initdata')

    use amrex_constants_module, only: ZERO, HALF, ONE, M_PI
    use probdata_module, only: rhocv, T0, Eexp, rexp
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS
    use network, only: nspec
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use prob_params_module, only: dg, coord_type
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(in   ) :: dx(3), problo(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer  :: i, j, k, ii, jj, kk
    real(rt) :: rho, cv, T, p, eint
    real(rt) :: rhoeexp, Vexp, rhoe0
    real(rt) :: xx, xcl, xcr, dx_sub
    real(rt) :: yy, ycl, ycr, dy_sub
    real(rt) :: zz, zcl, zcr, dz_sub
    real(rt) :: rr2, xcmin, xcmax, ycmin, ycmax, zcmin, zcmax , rcmin, rcmax
    real(rt) :: vol_pert, vol_ambient, T_zone, rhoe_zone

    type(eos_t) :: eos_state

    integer, parameter :: nsub = 64

    !$gpu

    Vexp = 4.e0_rt / 3.e0_rt * M_PI * rexp**3

    rhoeexp = Eexp / Vexp

    rho = 1.0e0_rt
    T = 1.0e0_rt

    eos_state % rho = rho
    eos_state % T   = T
    eos_state % xn  = ZERO
    eos_state % xn(1)  = 1.e0_rt

    call eos(eos_input_rt, eos_state)

    cv = eos_state % cv

    rho = rhocv / cv

    eint = rhoeexp / rho
    T = eint / cv

    rhoe0 = rho * cv * T0

    dx_sub = dx(1) / dble(nsub)
    dy_sub = dx(2) / dble(nsub)
    dz_sub = dx(3) / dble(nsub)

    do k = lo(3), hi(3)

       zcl = problo(3) + dx(3) * dble(k)
       zcr = zcl + dx(3)
       zcmin = min(abs(zcl), abs(zcr))
       zcmax = max(abs(zcl), abs(zcr))

       do j = lo(2), hi(2)

          ycl = problo(2) + dx(2) * dble(j)
          ycr = ycl + dx(2)
          ycmin = min(abs(ycl), abs(ycr))
          ycmax = max(abs(ycl), abs(ycr))

          do i = lo(1), hi(1)

             xcl = problo(1) + dx(1) * dble(i)
             xcr = xcl + dx(1)
             xcmin = min(abs(xcl), abs(xcr))
             xcmax = max(abs(xcl), abs(xcr))

             rcmin = sqrt(xcmin**2+ycmin**2+zcmin**2)
             rcmax = sqrt(xcmax**2+ycmax**2+zcmax**2)

             state(i,j,k,URHO) = rho
             state(i,j,k,UMX)  = ZERO
             state(i,j,k,UMY)  = ZERO
             state(i,j,k,UMZ)  = ZERO

             if (rcmin .ge. rexp) then
                T_zone = T0
                rhoe_zone = rhoe0
             else if (rcmax .le. rexp) then
                T_zone = T
                rhoe_zone = rhoeexp
             else
                vol_pert    = ZERO
                vol_ambient = ZERO

                do kk = 0, dg(3)*(nsub-1)
                   zz = zcl + (dble(kk) + HALF) * dz_sub

                   do jj = 0, dg(2)*(nsub-1)
                      yy = ycl + (dble(jj) + HALF) * dy_sub

                      do ii = 0, nsub-1
                         xx = xcl + (dble(ii) + HALF) * dx_sub

                         rr2 = xx**2 + yy**2 + zz**2

                         if (coord_type == 1) then
                            ! we are axisymmetric so weight by the radial coordinate
                            if (rr2 <= rexp**2) then
                               vol_pert = vol_pert + xx
                            else
                               vol_ambient = vol_ambient + xx
                            endif
                         else
                            ! Cartesian
                            if (rr2 <= rexp**2) then
                               vol_pert = vol_pert + ONE
                            else
                               vol_ambient = vol_ambient + ONE
                            endif
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



  subroutine ca_initrad(lo, hi, &
                        rad_state, r_lo, r_hi, &
                        dx, problo) &
                        bind(C, name='ca_initrad')

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: ZERO

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    real(rt), intent(in   ) :: dx(3), problo(3)
    real(rt), intent(inout) :: rad_state(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),0:NGROUPS-1)

    integer :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rad_state(i,j,k,:) = ZERO
          end do
       end do
    end do

  end subroutine ca_initrad

end module initdata_module
