module ppm_module

  implicit none

contains

     ! characteristics based on u
     subroutine ppm(s,qd_l1,qd_h1,u,cspd,Ip,Im,ilo,ihi,dx,dt)
       
       use meth_params_module, only : ppm_type

       implicit none
       
       integer          qd_l1,qd_h1
       integer          ilo,ihi
       double precision s(qd_l1:qd_h1)
       double precision u(qd_l1:qd_h1)
       double precision cspd(qd_l1:qd_h1)
       double precision Ip(ilo-1:ihi+1,1:3)
       double precision Im(ilo-1:ihi+1,1:3)
       double precision dx,dt

       ! local
       integer i
       logical extremum, bigp, bigm

       double precision dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham
       double precision sgn, sigma, s6, w0cc, amax, delam, delap
       double precision dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp

       ! s_{\ib,+}, s_{\ib,-}
       double precision, allocatable :: sp(:)
       double precision, allocatable :: sm(:)

       ! \delta s_{\ib}^{vL}
       double precision, allocatable :: dsvl(:)

       ! s_{i+\half}^{H.O.}
       double precision, allocatable :: sedge(:)

       ! cell-centered indexing
       allocate(sp(ilo-1:ihi+1))
       allocate(sm(ilo-1:ihi+1))

       ! constant used in Colella 2008
       C = 1.25d0

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

       ! compute van Leer slopes in x-direction
       dsvl = 0.d0
       do i=ilo-2,ihi+2
          dsc = 0.5d0 * (s(i+1) - s(i-1))
          dsl = 2.d0  * (s(i  ) - s(i-1))
          dsr = 2.d0  * (s(i+1) - s(i  ))
          if (dsl*dsr .gt. 0.d0) dsvl(i) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
       end do
       
       ! interpolate s to x-edges
       do i=ilo-1,ihi+2
          sedge(i) = 0.5d0*(s(i)+s(i-1)) - (1.d0/6.d0)*(dsvl(i)-dsvl(i-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i) = max(sedge(i),min(s(i),s(i-1)))
          sedge(i) = min(sedge(i),max(s(i),s(i-1)))
       end do

       ! copy sedge into sp and sm
       do i=ilo-1,ihi+1
          sp(i) = sedge(i+1)
          sm(i) = sedge(i  )
       end do

       ! modify using quadratic limiters
       do i=ilo-1,ihi+1
          if ((sp(i)-s(i))*(s(i)-sm(i)) .le. 0.d0) then
             sp(i) = s(i)
             sm(i) = s(i)
          else if (abs(sp(i)-s(i)) .ge. 2.d0*abs(sm(i)-s(i))) then
             sp(i) = 3.d0*s(i) - 2.d0*sm(i)
          else if (abs(sm(i)-s(i)) .ge. 2.d0*abs(sp(i)-s(i))) then
             sm(i) = 3.d0*s(i) - 2.d0*sp(i)
          end if
       end do

    else if (ppm_type .eq. 2) then
       
       ! interpolate s to x-edges
       do i=ilo-2,ihi+3
          sedge(i) = (7.d0/12.d0)*(s(i-1)+s(i)) - (1.d0/12.d0)*(s(i-2)+s(i+1))
          ! limit sedge
          if ((sedge(i)-s(i-1))*(s(i)-sedge(i)) .lt. 0.d0) then
             D2  = 3.d0*(s(i-1)-2.d0*sedge(i)+s(i))
             D2L = s(i-2)-2.d0*s(i-1)+s(i)
             D2R = s(i-1)-2.d0*s(i)+s(i+1)
             sgn = sign(1.d0,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
             sedge(i) = 0.5d0*(s(i-1)+s(i)) - (1.d0/6.d0)*D2LIM
          end if
       end do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
       do i=ilo-1,ihi+1

          alphap = sedge(i+1)-s(i)
          alpham = sedge(i  )-s(i)
          bigp = abs(alphap).gt.2.d0*abs(alpham)
          bigm = abs(alpham).gt.2.d0*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. 0.d0) then
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
             extremum = (dachkm*dachkp .le. 0.d0)
          end if

          if (extremum) then
             D2  = 6.d0*(alpham + alphap)
             D2L = s(i-2)-2.d0*s(i-1)+s(i)
             D2R = s(i)-2.d0*s(i+1)+s(i+2)
             D2C = s(i-1)-2.d0*s(i)+s(i+1)
             sgn = sign(1.d0,D2)
             D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn = sign(1.d0,alpham)
                amax = -alphap**2 / (4*(alpham + alphap))
                delam = s(i-1) - s(i)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -2.d0*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn = sign(1.d0,alphap)
                amax = -alpham**2 / (4*(alpham + alphap))
                delap = s(i+1) - s(i)
               if (sgn*amax .ge. sgn*delap) then
                  if (sgn*(delap - alphap).ge.1.d-10) then
                     alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                  else
                     alpham = -2.d0*alphap
                  endif
               endif
             end if
          end if

          sm(i) = s(i) + alpham
          sp(i) = s(i) + alphap

       end do

    end if

        ! compute x-component of Ip and Im
       do i=ilo-1,ihi+1
          s6 = 6.0d0*s(i) - 3.0d0*(sm(i)+sp(i))
          sigma = abs(u(i)-cspd(i))*dt/dx
          Ip(i,1) = sp(i) - &
               (sigma/2.0d0)*(sp(i)-sm(i)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,1) = sm(i) + &
               (sigma/2.0d0)*(sp(i)-sm(i)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i))*dt/dx
          Ip(i,2) = sp(i) - &
               (sigma/2.0d0)*(sp(i)-sm(i)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,2) = sm(i) + &
               (sigma/2.0d0)*(sp(i)-sm(i)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i)+cspd(i))*dt/dx
          Ip(i,3) = sp(i) - &
               (sigma/2.0d0)*(sp(i)-sm(i)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,3) = sm(i) + &
               (sigma/2.0d0)*(sp(i)-sm(i)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
       end do
       
       deallocate(sedge,dsvl,sp,sm)

     end subroutine ppm

end module ppm_module
