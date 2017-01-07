module ppm_module

  use bl_fort_module, only : rt => c_real
  implicit none

contains

  ! characteristics based on u
  subroutine ppm(s,qd_l1,qd_h1, &
                 u,cspd, &
                 flatn, &
                 Ip,Im,ilo,ihi,dx,dt)
       
    use meth_params_module, only : ppm_type
    use bl_constants_module

    use bl_fort_module, only : rt => c_real
    implicit none
       
    integer          qd_l1,qd_h1
    integer          ilo,ihi
    real(rt)         s(qd_l1:qd_h1)
    real(rt)         u(qd_l1:qd_h1)
    real(rt)         cspd(qd_l1:qd_h1)
    real(rt)         flatn(qd_l1:qd_h1)
    real(rt)         Ip(ilo-1:ihi+1,1:3)
    real(rt)         Im(ilo-1:ihi+1,1:3)
    real(rt)         dx,dt

    ! local
    integer i
    logical extremum, bigp, bigm

    real(rt)         dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    real(rt)         sgn, sigma, s6, amax, delam, delap
    real(rt)         dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp

    ! s_{\ib,+}, s_{\ib,-}
    real(rt)        , allocatable :: sp(:)
    real(rt)        , allocatable :: sm(:)
    
    ! \delta s_{\ib}^{vL}
    real(rt)        , allocatable :: dsvl(:)
    
    ! s_{i+\half}^{H.O.}
    real(rt)        , allocatable :: sedge(:)

    ! constant used in Colella 2008
    real(rt), parameter :: C = 1.25e0_rt

    ! a constant used for testing extrema
    real(rt), parameter :: SMALL = 1.e-10_rt    


    ! cell-centered indexing
    allocate(sp(ilo-1:ihi+1))
    allocate(sm(ilo-1:ihi+1))
    

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(ilo-2:ihi+2))

    ! edge-centered indexing for x-faces
    if (ppm_type .eq. 1) then
       allocate(sedge(ilo-1:ihi+2))
    else
       allocate(sedge(ilo-2:ihi+3))
    end if
    
    ! compute s at x-edges
    if (ppm_type .eq. 1) then

       ! compute van Leer slopes in x-direction (CW Eq. 1.7, 1.8 
       ! w/ zone widths (dxi) all equal)
       dsvl = ZERO
       do i=ilo-2,ihi+2
          dsc = HALF * (s(i+1) - s(i-1))
          dsl = TWO  * (s(i  ) - s(i-1))
          dsr = TWO  * (s(i+1) - s(i  ))
          if (dsl*dsr .gt. ZERO) dsvl(i) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
       end do
       
       ! interpolate s to x-edges (CW 1.6)
       do i=ilo-1,ihi+2
          sedge(i) = HALF*(s(i)+s(i-1)) - SIXTH*(dsvl(i)-dsvl(i-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i) = max(sedge(i),min(s(i),s(i-1)))
          sedge(i) = min(sedge(i),max(s(i),s(i-1)))
       end do

       ! copy sedge into sp and sm
       do i=ilo-1,ihi+1
          sp(i) = sedge(i+1)
          sm(i) = sedge(i  )
       end do

       ! flatten the parabola BEFORE doing the other
       ! monotonozation -- this is the method that Flash does
       do i=ilo-1,ihi+1
          sm(i) = flatn(i)*sm(i) + (ONE-flatn(i))*s(i)
          sp(i) = flatn(i)*sp(i) + (ONE-flatn(i))*s(i)
       enddo

       ! modify using quadratic limiters (CW 1.10) -- with a slightly
       ! different form from Colella & Sekora (Eqs. 14, 15)
       do i=ilo-1,ihi+1
          if ((sp(i)-s(i))*(s(i)-sm(i)) .le. ZERO) then
             sp(i) = s(i)
             sm(i) = s(i)
          else if (abs(sp(i)-s(i)) .ge. TWO*abs(sm(i)-s(i))) then
             sp(i) = THREE*s(i) - TWO*sm(i)
          else if (abs(sm(i)-s(i)) .ge. TWO*abs(sp(i)-s(i))) then
             sm(i) = THREE*s(i) - TWO*sp(i)
          end if
       end do

       ! flatten the parabola AFTER doing the monotonization --
       ! this is the method that Miller & Colella do
       !do i=ilo-1,ihi+1
       !      sm(i) = flatn(i)*sm(i) + (ONE-flatn(i))*s(i)
       !      sp(i) = flatn(i)*sp(i) + (ONE-flatn(i))*s(i)
       !enddo


    else if (ppm_type .eq. 2) then
       
       ! interpolate s to x-edges
       do i=ilo-2,ihi+3
          sedge(i) = SEVEN12TH*(s(i-1)+s(i)) - TWELFTH*(s(i-2)+s(i+1))
          ! limit sedge
          if ((sedge(i)-s(i-1))*(s(i)-sedge(i)) .lt. ZERO) then
             D2  = THREE*(s(i-1)-TWO*sedge(i)+s(i))
             D2L = s(i-2)-TWO*s(i-1)+s(i)
             D2R = s(i-1)-TWO*s(i)+s(i+1)
             sgn = sign(ONE,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),ZERO)
             sedge(i) = HALF*(s(i-1)+s(i)) - SIXTH*D2LIM
          end if
       end do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
       do i=ilo-1,ihi+1

          alphap = sedge(i+1)-s(i)
          alpham = sedge(i  )-s(i)
          bigp = abs(alphap).gt.TWO*abs(alpham)
          bigm = abs(alpham).gt.TWO*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. ZERO) then
             extremum = .true.
          else if (bigp .or. bigm) then
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             dafacem = sedge(i) - sedge(i-1)
             dafacep = sedge(i+2) - sedge(i+1)
             dabarm = s(i) - s(i-1)
             dabarp = s(i+1) - s(i)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin= min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. ZERO)
          end if

          if (extremum) then
             D2  = SIX*(alpham + alphap)
             D2L = s(i-2)-TWO*s(i-1)+s(i)
             D2R = s(i)-TWO*s(i+1)+s(i+2)
             D2C = s(i-1)-TWO*s(i)+s(i+1)
             sgn = sign(ONE,D2)
             D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),ZERO)
             alpham = alpham*D2LIM/max(abs(D2),SMALL)
             alphap = alphap*D2LIM/max(abs(D2),SMALL)
          else
             if (bigp) then
                sgn = sign(ONE,alpham)
                amax = -alphap**2 / (4*(alpham + alphap))
                delam = s(i-1) - s(i)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge. SMALL) then
                      alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -TWO*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn = sign(ONE,alphap)
                amax = -alpham**2 / (4*(alpham + alphap))
                delap = s(i+1) - s(i)
               if (sgn*amax .ge. sgn*delap) then
                  if (sgn*(delap - alphap).ge. SMALL) then
                     alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                  else
                     alpham = -TWO*alphap
                  endif
               endif
             end if
          end if

          sm(i) = s(i) + alpham
          sp(i) = s(i) + alphap

       end do

       ! flatten the parabola AFTER doing the monotonization --
       ! this is the method that Miller & Colella do
       do i=ilo-1,ihi+1
             sm(i) = flatn(i)*sm(i) + (ONE-flatn(i))*s(i)
             sp(i) = flatn(i)*sp(i) + (ONE-flatn(i))*s(i)
       enddo

    end if

        ! compute x-component of Ip and Im
       do i=ilo-1,ihi+1

          ! Ip/m is the integral under the parabola for the extent
          ! that a wave can travel over a timestep
          !                                              
          ! Ip integrates to the right edge of a cell   
          ! Im integrates to the left edge of a cell
        
          s6 = SIX*s(i) - THREE*(sm(i)+sp(i))

          ! u-c wave
          sigma = abs(u(i)-cspd(i))*dt/dx

          if (u(i)-cspd(i) <= ZERO) then
             Ip(i,1) = sp(i)
          else
             Ip(i,1) = sp(i) - &
                  HALF*sigma*(sp(i)-sm(i)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i)-cspd(i) >= ZERO) then
             Im(i,1) = sm(i)
          else
             Im(i,1) = sm(i) + &
                  HALF*sigma*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! u wave
          sigma = abs(u(i))*dt/dx

          if (u(i) <= ZERO) then
             Ip(i,2) = sp(i)
          else
             Ip(i,2) = sp(i) - &
                  HALF*sigma*(sp(i)-sm(i)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i) >= ZERO) then
             Im(i,2) = sm(i)
          else
             Im(i,2) = sm(i) + &
                  HALF*sigma*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          endif

          ! u+c wave
          sigma = abs(u(i)+cspd(i))*dt/dx

          if (u(i) + cspd(i) <= ZERO) then
             Ip(i,3) = sp(i)
          else
             Ip(i,3) = sp(i) - &
                  HALF*sigma*(sp(i)-sm(i)-(ONE-TWO3RD*sigma)*s6)
          endif

          if (u(i) + cspd(i) >= ZERO) then
             Im(i,3) = sm(i)
          else
             Im(i,3) = sm(i) + &
                  HALF*sigma*(sp(i)-sm(i)+(ONE-TWO3RD*sigma)*s6)
          endif
       end do
       
       deallocate(sedge,dsvl,sp,sm)

     end subroutine ppm

end module ppm_module
