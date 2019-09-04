module ppm_module
  !
  ! this does the parabolic reconstruction on a variable and the (optional)
  ! integration under the characteristic domain of the parabola

  use amrex_constants_module, only: ZERO, HALF, ONE, TWO, SIXTH, &
       TWO3RD, THREE, SIX, SEVEN12TH, TWELFTH
  use prob_params_module, only : dg
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains

  subroutine ca_ppm_reconstruct(i, j, k, &
                                idir, &
                                s, s_lo, s_hi, nc, nstart, nend, &
                                flatn, f_lo, f_hi, &
                                qm, qp, &
                                ncq, nqstart, nqend) bind(c, name='ca_ppm_reconstruct')
    ! this routine does the reconstruction of the zone data into
    ! parabola.  The loops are over zone centers and for zone center,
    ! it will compute the left and right values of the parabola and
    ! store these in the qm and qp arrays.  
    !
    ! here s has nc components and we reconstruct component n
    !
    ! we store the result in the arrays qm and qp, in component nq out
    ! of ncq components



#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif
    use prob_params_module, only : dim

    implicit none

    integer, intent(in), value :: idir
    integer, intent(in   ) :: s_lo(3), s_hi(3)
    integer, intent(in   ) :: f_lo(3), f_hi(3)
    integer, intent(in   ), value :: nstart, nend, nc
    integer, intent(in   ), value :: nqstart, nqend, ncq

    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), nc)
    real(rt), intent(in   ) :: flatn(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3))
    real(rt), intent(inout) :: qm(ncq, AMREX_SPACEDIM)
    real(rt), intent(inout) :: qp(ncq, AMREX_SPACEDIM)

    integer, intent(in) :: i, j, k

    ! local
    integer :: n, nq, ic

    real(rt) :: dsl, dsr, dsc
    real(rt) :: dsvl_l, dsvl_r

    ! s_{\ib,+}, s_{\ib,-}
    real(rt) :: sm, sp

    !$gpu

#ifndef AMREX_USE_GPU
    if (nend - nstart /= nqend - nqstart) then
       call castro_error("Number of components do not match in ca_ppm_reconstruct")
    end if
#endif

    if (idir == 1) then

       do ic = 0, (nend - nstart)

          n = nstart + ic
          nq = nqstart + ic

          !!!!!!!!!!!!!!!!!!!!
          ! x-direction
          !!!!!!!!!!!!!!!!!!!!
          ! Compute van Leer slopes

          dsl = TWO  * (s(i-1,j,k,n) - s(i-2,j,k,n))
          dsr = TWO  * (s(i  ,j,k,n) - s(i-1,j,k,n))
          if (dsl*dsr .gt. ZERO) then
             dsc = HALF * (s(i  ,j,k,n) - s(i-2,j,k,n))
             dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl_l = ZERO
          end if

          dsl = TWO  * (s(i  ,j,k,n) - s(i-1,j,k,n))
          dsr = TWO  * (s(i+1,j,k,n) - s(i  ,j,k,n))
          if (dsl*dsr .gt. ZERO) then
             dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
             dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl_r = ZERO
          end if

          ! Interpolate s to x-edges

          sm = HALF*(s(i,j,k,n)+s(i-1,j,k,n)) - SIXTH*(dsvl_r - dsvl_l)

          ! Make sure sedge lies in between adjacent cell-centered values

          sm = max(sm, min(s(i,j,k,n),s(i-1,j,k,n)))
          sm = min(sm, max(s(i,j,k,n),s(i-1,j,k,n)))

          ! Compute van Leer slopes

          dsl = TWO  * (s(i  ,j,k,n) - s(i-1,j,k,n))
          dsr = TWO  * (s(i+1,j,k,n) - s(i  ,j,k,n))
          if (dsl*dsr .gt. ZERO) then
             dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
             dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl_l = ZERO
          end if

          dsl = TWO  * (s(i+1,j,k,n) - s(i  ,j,k,n))
          dsr = TWO  * (s(i+2,j,k,n) - s(i+1,j,k,n))
          if (dsl*dsr .gt. ZERO) then
             dsc = HALF * (s(i+2,j,k,n) - s(i  ,j,k,n))
             dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl_r = ZERO
          end if

          ! Interpolate s to x-edges

          sp = HALF*(s(i+1,j,k,n)+s(i,j,k,n)) - SIXTH*(dsvl_r - dsvl_l)

          ! Make sure sedge lies in between adjacent cell-centered values

          sp = max(sp, min(s(i+1,j,k,n),s(i,j,k,n)))
          sp = min(sp, max(s(i+1,j,k,n),s(i,j,k,n)))

          ! Flatten the parabola
          sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
          sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

          ! Modify using quadratic limiters -- note this version of the limiting comes
          ! from Colella and Sekora (2008), not the original PPM paper.
          if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then

             sp = s(i,j,k,n)
             sm = s(i,j,k,n)

          else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then

             sp = THREE*s(i,j,k,n) - TWO*sm

          else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then

             sm = THREE*s(i,j,k,n) - TWO*sp

          end if

          qp(nq,1) = sp
          qm(nq,1) = sm

       end do

    else if (idir == 2) then

       do ic = 0, (nend - nstart)

          n = nstart + ic
          nq = nqstart + ic

          !!!!!!!!!!!!!!!!!!!!
          ! y-direction
          !!!!!!!!!!!!!!!!!!!!

          ! Compute van Leer slopes

          dsl = TWO  * (s(i,j-1,k,n) - s(i,j-2,k,n))
          dsr = TWO  * (s(i,j  ,k,n) - s(i,j-1,k,n))
          if (dsl*dsr .gt. ZERO) then
             dsc = HALF * (s(i,j  ,k,n) - s(i,j-2,k,n))
             dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl_l = ZERO
          end if

          dsl = TWO  * (s(i,j  ,k,n) - s(i,j-1,k,n))
          dsr = TWO  * (s(i,j+1,k,n) - s(i,j  ,k,n))
          if (dsl*dsr .gt. ZERO) then
             dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
             dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl_r = ZERO
          end if

          ! Interpolate s to y-edges

          sm = HALF*(s(i,j,k,n)+s(i,j-1,k,n)) - SIXTH*(dsvl_r - dsvl_l)

          ! Make sure sedge lies in between adjacent cell-centered values

          sm = max(sm, min(s(i,j,k,n),s(i,j-1,k,n)))
          sm = min(sm, max(s(i,j,k,n),s(i,j-1,k,n)))

          ! Compute van Leer slopes

          dsl = TWO  * (s(i,j  ,k,n) - s(i,j-1,k,n))
          dsr = TWO  * (s(i,j+1,k,n) - s(i,j  ,k,n))
          if (dsl*dsr .gt. ZERO) then
             dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
             dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl_l = ZERO
          end if

          dsl = TWO  * (s(i,j+1,k,n) - s(i,j  ,k,n))
          dsr = TWO  * (s(i,j+2,k,n) - s(i,j+1,k,n))
          if (dsl*dsr .gt. ZERO) then
             dsc = HALF * (s(i,j+2,k,n) - s(i,j  ,k,n))
             dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl_r = ZERO
          end if

          ! Interpolate s to y-edges

          sp = HALF*(s(i,j+1,k,n)+s(i,j,k,n)) - SIXTH*(dsvl_r - dsvl_l)

          ! Make sure sedge lies in between adjacent cell-centered values

          sp = max(sp, min(s(i,j+1,k,n),s(i,j,k,n)))
          sp = min(sp, max(s(i,j+1,k,n),s(i,j,k,n)))

          ! Flatten the parabola

          sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
          sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

          ! Modify using quadratic limiters

          if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then

             sp = s(i,j,k,n)
             sm = s(i,j,k,n)

          else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then

             sp = THREE*s(i,j,k,n) - TWO*sm

          else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then

             sm = THREE*s(i,j,k,n) - TWO*sp

          end if

          qp(nq,2) = sp
          qm(nq,2) = sm

       end do

    else

       do ic = 0, (nend - nstart)

          n = nstart + ic
          nq = nqstart + ic

          !!!!!!!!!!!!!!!!!!!!
          ! z-direction
          !!!!!!!!!!!!!!!!!!!!

          ! Compute van Leer slopes

          dsl = TWO  * (s(i,j,k-1,n) - s(i,j,k-2,n))
          dsr = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
          if (dsl*dsr .gt. ZERO) then
             dsc = HALF * (s(i,j,k  ,n) - s(i,j,k-2,n))
             dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl_l = ZERO
          end if

          dsl = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
          dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
          if (dsl*dsr .gt. ZERO) then
             dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
             dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl_r = ZERO
          end if

          ! Interpolate s to z-edges

          sm = HALF*(s(i,j,k,n)+s(i,j,k-1,n)) - SIXTH*(dsvl_r - dsvl_l)

          ! Make sure sedge lies in between adjacent cell-centered values

          sm = max(sm, min(s(i,j,k,n),s(i,j,k-1,n)))
          sm = min(sm, max(s(i,j,k,n),s(i,j,k-1,n)))

          ! Compute van Leer slopes

          dsl = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
          dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
          if (dsl*dsr .gt. ZERO) then
             dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
             dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl_l = ZERO
          end if

          dsl = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
          dsr = TWO  * (s(i,j,k+2,n) - s(i,j,k+1,n))
          if (dsl*dsr .gt. ZERO) then
             dsc = HALF * (s(i,j,k+2,n) - s(i,j,k  ,n))
             dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          else
             dsvl_r = ZERO
          end if

          ! Interpolate s to z-edges

          sp = HALF*(s(i,j,k+1,n)+s(i,j,k,n)) - SIXTH*(dsvl_r - dsvl_l)

          ! Make sure sedge lies in between adjacent cell-centered values

          sp = max(sp, min(s(i,j,k+1,n),s(i,j,k,n)))
          sp = min(sp, max(s(i,j,k+1,n),s(i,j,k,n)))

          ! Flatten the parabola

          sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
          sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

          ! Modify using quadratic limiters

          if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then

             sp = s(i,j,k,n)
             sm = s(i,j,k,n)

          else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then

             sp = THREE*s(i,j,k,n) - TWO*sm

          else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then

             sm = THREE*s(i,j,k,n) - TWO*sp

          end if

          qp(nq,3) = sp
          qm(nq,3) = sm

       end do
    end if

  end subroutine ca_ppm_reconstruct



  subroutine ppm_int_profile(i, j, k, &
                             idir, &
                             s, s_lo, s_hi, ncomp, n, &
                             q, qd_lo, qd_hi, &
                             qaux, qa_lo, qa_hi, &
                             sm_in, sp_in, &
                             Ip, Im, icomp, ic, &
                             dx, dt)

    use meth_params_module, only : NQAUX, QC, NQ, QU, QV, QW

    implicit none

    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) :: ncomp, n
    integer, intent(in) :: icomp, ic
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in), value :: idir

    real(rt), intent(in) ::     s( s_lo(1): s_hi(1), s_lo(2): s_hi(2), s_lo(3): s_hi(3), ncomp)
    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3), NQ)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(in) :: sm_in(AMREX_SPACEDIM)
    real(rt), intent(in) :: sp_in(AMREX_SPACEDIM)
    real(rt), intent(inout) :: Ip(1:3, icomp)
    real(rt), intent(inout) :: Im(1:3, icomp)

    real(rt), intent(in) :: dx(3), dt
    integer, intent(in) :: i,j,k

    ! local
    real(rt) :: speed

    real(rt)         dtdx, dtdy, dtdz
    real(rt)         sigma, s6
    real(rt)         :: sm, sp

    !$gpu

    dtdx = dt/dx(1)
#if (AMREX_SPACEDIM >= 2)
    dtdy = dt/dx(2)
#endif
#if (AMREX_SPACEDIM == 3)
    dtdz = dt/dx(3)
#endif

    if (idir == 1) then

       !!!!!!!!!!!!!!!!!!!!
       ! x-direction
       !!!!!!!!!!!!!!!!!!!!

       ! copy sedge into sp and sm
       sp = sp_in(1)
       sm = sm_in(1)

       ! compute x-component of Ip and Im
       s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)

       ! Ip/m is the integral under the parabola for the extent
       ! that a wave can travel over a timestep
       !
       ! Ip integrates to the right edge of a cell
       ! Im integrates to the left edge of a cell

       ! u-c wave
       speed = q(i,j,k,QU)-qaux(i,j,k,QC)
       sigma = abs(speed)*dtdx

       ! if speed == ZERO, then either branch is the same
       if (speed <= ZERO) then
          Ip(1,ic) = sp
          Im(1,ic) = sm + &
               HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
       else
          Ip(1,ic) = sp - &
               HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
          Im(1,ic) = sm
       endif

       ! u wave
       speed = q(i,j,k,QU)
       sigma = abs(speed)*dtdx

       if (speed <= ZERO) then
          Ip(2,ic) = sp
          Im(2,ic) = sm + &
               HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
       else
          Ip(2,ic) = sp - &
               HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
          Im(2,ic) = sm
       endif

       ! u+c wave
       speed = q(i,j,k,QU)+qaux(i,j,k,QC)
       sigma = abs(speed)*dtdx

       if (speed <= ZERO) then
          Ip(3,ic) = sp
          Im(3,ic) = sm + &
               HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
       else
          Ip(3,ic) = sp - &
               HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
          Im(3,ic) = sm
       endif

    else if (idir == 2) then

       !!!!!!!!!!!!!!!!!!!!
       ! y-direction
       !!!!!!!!!!!!!!!!!!!!

       ! copy sedge into sp and sm
       sp = sp_in(2)
       sm = sm_in(2)

       ! compute y-component of Ip and Im
       s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)

       ! v-c wave
       speed = q(i,j,k,QV)-qaux(i,j,k,QC)
       sigma = abs(speed)*dtdy

       if (speed <= ZERO) then
          Ip(1,ic) = sp
          Im(1,ic) = sm + &
               HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
       else
          Ip(1,ic) = sp - &
               HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
          Im(1,ic) = sm
       endif

       ! v wave
       speed = q(i,j,k,QV)
       sigma = abs(speed)*dtdy

       if (speed <= ZERO) then
          Ip(2,ic) = sp
          Im(2,ic) = sm + &
               HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
       else
          Ip(2,ic) = sp - &
               HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
          Im(2,ic) = sm
       endif

       ! v+c wave
       speed = q(i,j,k,QV)+qaux(i,j,k,QC)
       sigma = abs(speed)*dtdy

       if (speed <= ZERO) then
          Ip(3,ic) = sp
          Im(3,ic) = sm + &
               HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
       else
          Ip(3,ic) = sp - &
               HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
          Im(3,ic) = sm
       endif

    else
       !!!!!!!!!!!!!!!!!!!!
       ! z-direction
       !!!!!!!!!!!!!!!!!!!!

       sp = sp_in(3)
       sm = sm_in(3)

       ! compute z-component of Ip and Im
       s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)

       ! w-c wave
       speed = q(i,j,k,QW)-qaux(i,j,k,QC)
       sigma = abs(speed)*dtdz

       if (speed <= ZERO) then
          Ip(1,ic) = sp
          Im(1,ic) = sm + &
               HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
       else
          Ip(1,ic) = sp - &
               HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
          Im(1,ic) = sm
       endif

       ! w wave
       speed = q(i,j,k,QW)
       sigma = abs(speed)*dtdz

       if (speed <= ZERO) then
          Ip(2,ic) = sp
          Im(2,ic) = sm + &
               HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
       else
          Ip(2,ic) = sp - &
               HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
          Im(2,ic) = sm
       endif

       ! w+c wave
       speed = q(i,j,k,QW)+qaux(i,j,k,QC)
       sigma = abs(speed)*dtdz

       if (speed <= ZERO) then
          Ip(3,ic) = sp
          Im(3,ic) = sm + &
               HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
       else
          Ip(3,ic) = sp - &
               HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
          Im(3,ic) = sm
       endif

    endif

  end subroutine ppm_int_profile




  subroutine ppm_reconstruct_with_eos(Ip, Im, Ip_gc, Im_gc)
    ! temperature-based PPM -- if desired, take the Ip(T)/Im(T)
    ! constructed above and use the EOS to overwrite Ip(p)/Im(p)
    ! get an edge-based gam1 here if we didn't get it from the EOS
    ! call above (for ppm_temp_fix = 1)

    use meth_params_module, only : NQ, QRHO, QTEMP, QPRES, QREINT, QFS, QFX, small_temp
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use network, only : nspec, naux

    real(rt), intent(inout) :: Ip(1:3, NQ)
    real(rt), intent(inout) :: Im(1:3, NQ)
    real(rt), intent(inout) :: Ip_gc(1:3, 1)
    real(rt), intent(inout) :: Im_gc(1:3, 1)

    integer :: iwave

    type(eos_t) :: eos_state

    !$gpu

    do iwave = 1, 3

       eos_state % rho = Ip(iwave,QRHO)
       eos_state % T   = max(Ip(iwave,QTEMP), small_temp)

       eos_state % xn  = Ip(iwave,QFS:QFS+nspec-1)
       eos_state % aux = Ip(iwave,QFX:QFX+naux-1)

       call eos(eos_input_rt, eos_state)

       Ip(iwave,QPRES)  = eos_state % p
       Ip(iwave,QREINT) = Ip(iwave,QRHO) * eos_state % e
       Ip_gc(iwave,1)   = eos_state % gam1


       eos_state % rho = Im(iwave,QRHO)
       eos_state % T   = max(Im(iwave,QTEMP), small_temp)

       eos_state % xn  = Im(iwave,QFS:QFS+nspec-1)
       eos_state % aux = Im(iwave,QFX:QFX+naux-1)

       call eos(eos_input_rt, eos_state)

       Im(iwave,QPRES)  = eos_state % p
       Im(iwave,QREINT) = Im(iwave,QRHO) * eos_state % e
       Im_gc(iwave,1)   = eos_state % gam1

    end do

  end subroutine ppm_reconstruct_with_eos

end module ppm_module
